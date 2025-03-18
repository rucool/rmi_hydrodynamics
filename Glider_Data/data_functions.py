import pandas as pd
import matplotlib.pyplot as plt
from gsw import SA_from_SP, CT_from_t, rho, p_from_z, pot_rho_t_exact
import time
import numpy as np
import scipy
import datetime as dt
from collections import namedtuple
from erddapy import ERDDAP
from pprint import pprint
import multiprocessing
import matplotlib.dates as mdates
import geopandas



def density(temperature, depth, salinity, latitude, longitude):
    """
    Calculates density given practical salinity, depth, latitude,
    and longitude using Gibbs gsw SA_from_SP and rho functions.

    Args:
        temperature (_type_): in-situ temperature (C)
        depth (_type_): depth, positive up (m)
        salinity (array): salinity
        latitude (array): latitude (decimal degrees)
        longitude (array): longitude (decimal degrees)

    Returns:
        density: Density calculated using the Gibbs GSW
    """

    # Calculates sea pressure from height using computationally-efficient 
    # 75-term expression for density, in terms of SA, CT and p 
    # (Roquet et al., 2015). 
    pressure = p_from_z(
        depth,
        latitude,
    )

    # Calculates Absolute Salinity from Practical Salinity. 
    # Since SP is non-negative by definition,
    # this function changes any negative input values of SP to be zero.
    absolute_salinity = SA_from_SP(
        salinity,
        pressure,
        longitude,
        latitude
    )

    # Calculates Conservative Temperature of seawater from in-situ temperature.
    conservative_temperature = CT_from_t(
        absolute_salinity,
        temperature,
        pressure
    )

    # Calculates in-situ density from Absolute Salinity and
    # Conservative Temperature, using the computationally-efficient expression 
    # for specific volume in terms of SA, CT and p (Roquet et al., 2015).
    density = rho(
        absolute_salinity,
        conservative_temperature,
        pressure
    )

    return density

def prof_dens(pdf):
    pdf['density'] = density (pdf.temperature.values,-pdf.depth.values,pdf.salinity.values,pdf.lat.values,pdf.lon.values)
    return pdf

def potential_density(temperature, depth, salinity, latitude, longitude):
    """
    Calculates density given practical salinity, depth, latitude,
    and longitude using Gibbs gsw SA_from_SP and rho functions.

    Args:
        temperature (_type_): in-situ temperature (C)
        depth (_type_): depth, positive up (m)
        salinity (array): salinity
        latitude (array): latitude (decimal degrees)
        longitude (array): longitude (decimal degrees)

    Returns:
        potential density: Density calculated using the Gibbs GSW
    """

    # Calculates sea pressure from height using computationally-efficient 
    # 75-term expression for density, in terms of SA, CT and p 
    # (Roquet et al., 2015). 
    pressure = p_from_z(
        depth,
        latitude,
    )

    # Calculates Absolute Salinity from Practical Salinity. 
    # Since SP is non-negative by definition,
    # this function changes any negative input values of SP to be zero.
    absolute_salinity = SA_from_SP(
        salinity,
        pressure,
        longitude,
        latitude
    )

   

    # Calculates potential density from Absolute Salinity and pressure
    # using the computationally-efficient expression 
    # for specific volume in terms of SA, t and p (Roquet et al., 2015).
    potential_density = pot_rho_t_exact(
        absolute_salinity,
        temperature,
        pressure,
        0)

    return potential_density

def prof_potdens(pdf):
    pdf['potential_density'] = potential_density (pdf.temperature.values,-pdf.depth.values,pdf.salinity.values,pdf.lat.values,pdf.lon.values)
    return pdf

def potential_energy_anomaly100(depth, dens,z_thresh):
    """
    This function calculates the potential energy anomaly
    (Simpson J, Brown J, Matthews J, Allen G (1990) Tidal straining, density
    currents and stirring in the control of estuarine stratification.
    Estuaries 13(2):125-132), in the top 100 meters

    Args:
        depth (numpy.ndarray or pandas.Series or xarray.Series): depth
        dens (numpy.ndarray or pandas.Series or xarray.Series): density

    Returns:
        np.ndarray: potential energy anomaly in J/m^3
    """

    g = 9.8 #m/s
    # dindex = np.fliplr(np.where(np.asarray(np.abs(depth)) <= 100))[0]
    dindex = np.fliplr(np.where(np.asarray(np.abs(depth)) <= z_thresh))[0]
    if len(dindex) == 0:
        PEA = np.nan
    else:
        zz = np.asarray(np.abs(depth.iloc[dindex]))
        denss = np.asarray(dens.iloc[dindex])
        ok = np.isfinite(denss)
        z = zz[ok]
        densi = denss[ok]
        depth_diff=abs(np.nanmax(zz)-np.nanmin(zz))
    
        if len(z)==0 or len(densi)==0 or depth_diff < 10 or np.isnan(depth_diff):
            PEA = np.nan
        else:
            if z[-1] - z[0] > 0:
                # So PEA is < 0
                # sign = -1
                # Adding 0 to sigma integral is normalized
                z = np.append(0,z)
                # densit=densi
            else:
                # So PEA is < 0
                # sign = 1
                # Adding 0 to sigma integral is normalized
                z = np.flipud(z)
                z = np.append(0,z)
            densit = np.flipud(densi)

            # adding density at depth = 0
            densitt = np.interp(z,z[1:],densit)
            density = np.flipud(densitt)

            # defining sigma
            max_depth = np.nanmax(zz[ok])
            sigma = -1*z/max_depth
            sigma = np.flipud(sigma)

            rhomean = np.trapz(density,sigma,axis=0)
            drho = rhomean - density
            torque = drho * sigma
            PEA = g * max_depth * np.trapz(torque,sigma,axis=0)

    return PEA
def prof_pea(pdf):
    pdf=pdf.sort_values('depth',ascending=True)
    tdf=depth_interpolate(pdf,depth_var='depth',depth_min=0,depth_max=pdf.depth.max(),stride=0.25,method='linear')
    z_thresh=np.nanmax(tdf.depth.values)
    pea = potential_energy_anomaly100 (tdf.depth,tdf.density,z_thresh)
    return pea



def phi_carpenter(rho,z,dz,max_depth):
    """
    Based on Carpenter et al 2016
    rho - rho profile
    z - depth profile (z oriented upwards +)
    dz - step for interpolation (meters)

    returns pea kj/m^2
    """

    g = 9.8 #m/s

    # convert z to height above bottom
    z_new=-1*(z-max_depth)

    # get index of sorted array for monotonically increasing z
    sorted_i = np.argsort(z_new)

    # sort z and rho by increasing z
    z_new = z_new[sorted_i]
    rho_new = rho[sorted_i]  
  
    #calculate depth averaged rho
    rho_mix=np.trapz(rho_new, z_new, axis=0)*(1/max_depth)

    drho = (rho_mix-rho_new)*(z_new)*g
    phi=np.trapz(drho, z_new, axis=0)/1000

    return phi




def prof_phi(pdf):
    """
    pdf - single glider profile
    """
    
    pdf=pdf.sort_values('depth',ascending=True)
    tdf=depth_interpolate(pdf,depth_var='depth',depth_min=0,depth_max=pdf.depth.max(),stride=0.25,method='linear')
    
    # rho_avg=np.nanmean(tdf.density.values)
    rho =tdf.potential_density.values
    z=tdf.depth.values
    max_depth=z.max()
    dz=0.25

    pea=phi_carpenter(rho,-z,dz,max_depth)
    return pea


def depth_interpolate(df, depth_var='depth', depth_min=0, depth_max=1000, stride=10, method='linear', index=None):
    """

    Args:
        df (pd.DataFrame): Depth profile in the form of a pandas dataframe
        depth_var (str, optional): Name of the depth variable in the dataframe. Defaults to 'depth'.
        depth_min (float or string, optional): Shallowest bin depth. Pass 'round' to round to nearest minimumdepth. Defaults to None.
        depth_max (float or string, optional): Deepest bin depth. Pass 'round' to round to nearest maximum depth. Defaults to None.
        stride (int, optional): Amount of space between bins. Defaults to 10.
        method (str, optional): Interpolation type. Defaults to 'linear'.

    Returns:
        pd.DataFrame: dataframe with depth interpolated
    """
    if df.empty:
        print("Dataframe empty. Returning to original function")
        return

    if isinstance(depth_min, str):
        if depth_min == 'round':
            depth_min = round(df[depth_var].min())
        else:
            depth_min = int(depth_min)

    if isinstance(depth_min, str):
        if depth_min == 'round':
            depth_max = round(df[depth_var].max())
        else:
            depth_max = int(depth_max)

    bins = np.arange(depth_min, depth_max+stride, stride)

    # Create temporary dataframe to interpolate to dz m depths
    temp = df.set_index(depth_var)  #set index to depth
    temp = temp[~temp.index.duplicated()]  #Remove duplicated indexs (not sure why there would be duplicates)
    temp = temp.reindex(temp.index.union(bins))  # reindex to depths in bins
    # temp = temp.drop('time', axis=1).interpolate(method=method, limit_direction='both')  # drop time and interpolate new depth indexes
    temp = temp.interpolate(method=method, limit_direction='both')
    temp = temp.reindex(index=bins)  # only want to see new_index data
    temp = temp.reset_index()  # reset index so you can access the depth variable

    if index:
        temp = temp.set_index(index)

    return temp





def get_vars_bydepth(pdf):
    pdf=pdf.sort_values('depth',ascending=True)

    surf_depth = pdf.iloc[0].depth
    bot_depth = pdf.iloc[-1].depth

    surf_temp = pdf.iloc[0].temperature
    bot_temp = pdf.iloc[-1].temperature

    surf_sal = pdf.iloc[0].salinity
    bot_sal = pdf.iloc[-1].salinity

    lat=pdf.iloc[0].lat
    lon=pdf.iloc[0].lon

    return surf_depth, bot_depth , surf_temp, bot_temp, surf_sal, bot_sal,lat,lon

def calc_mld_upper(df,ref_depth,deltaT,_z,_temp,_sal,lat,lon):
    #find the depth index closest to reference depth
    zdiff = _z - ref_depth
    dref=np.where((zdiff>0) & (zdiff<=1) & (~np.isnan(zdiff)))[0]
    
    # dref=np.argmin(np.abs(np.abs(_z)-ref_depth))
    if (dref.size!=0) and (~np.isnan(_z[dref[0]])):
        dref=dref[0]
    # if (np.abs(_z[dref]-ref_depth)<=1) and (~np.isnan(_z[dref])):
        
        pres=p_from_z(-_z,lat)
        SA = SA_from_SP(_sal,pres,lon,lat)
        #calculate intial rho
        ini_rho=pot_rho_t_exact(SA,_temp,pres,0)

        if (~np.isnan(SA[dref])) and (~np.isnan(_temp[dref])) and (~np.isnan(pres[dref])):
            # ini_rho=pot_rho_t_exact(SA[dref:],_temp[dref:],pres[dref:],0)
            #calculate rho profile using reference depth
            refdepth_rho=pot_rho_t_exact(SA[dref],_temp[dref],pres[dref],0)
            #calculate rho -deltaT reference depth
            refdepth_rho_deltaT =pot_rho_t_exact(SA[dref],_temp[dref]-deltaT,pres[dref],0)
            #calculate delta sigma
            delta_sig=np.abs(refdepth_rho_deltaT - refdepth_rho)
            #calculate difference between rho profile and rho at reference depth
            dens_diff=np.abs(ini_rho[dref:] - refdepth_rho)

            #grab first index where this is true
            if np.where(dens_diff>=delta_sig)[0].size>0:
                
                mld_idx=np.where(dens_diff>=delta_sig)[0][0]
                mldU=_z[dref:][mld_idx]
                qcU=1
                
            else:
                mldU=df.depth.iloc[-1]
                qcU=2 #profile doesn't surpass delta_sigma threshold

        else:
            mldU=np.nan
            qcU= 5 #one of T S P being nan cuased related density vars to be nan       
    else:
        mldU=np.nan
        qcU= 3 #profile doesn't have depth close enough to ref_depth or depth[ref_depth] is nan
    
    return mldU,qcU

def calc_mld_lower(df,deltaT,_z,_temp,_sal,lat,lon):
    bot_depth = df.depth.iloc[-1]
    dref=np.where(_z==bot_depth)[0]
    if (dref.size!=0) and (~np.isnan(dref[-1])):

        dref=dref[-1]
        pres=p_from_z(-_z,lat)
        SA = SA_from_SP(_sal,pres,lon,lat)
        #calculate intial rho
        ini_rho=pot_rho_t_exact(SA,_temp,pres,0)

        if (~np.isnan(SA[dref])) and (~np.isnan(_temp[dref])) and (~np.isnan(pres[dref])):
            
            #calculate rho profile using reference depth
            refdepth_rho=pot_rho_t_exact(SA[dref],_temp[dref],pres[dref],0)
            #calculate rho -deltaT reference depth
            refdepth_rho_deltaT =pot_rho_t_exact(SA[dref],_temp[dref]-deltaT,pres[dref],0)
            #calculate delta sigma
            delta_sig=np.abs(refdepth_rho_deltaT - refdepth_rho)
            #calculate difference between rho profile and rho at reference depth
            dens_diff=np.abs(ini_rho - refdepth_rho)
            # dens_diff=np.abs(ini_rho - ini_rho[-1])

            #grab first index where this is true
            if np.where(dens_diff>delta_sig)[0].size>0:
                
                mld_idx=np.where(dens_diff>delta_sig)[0][-1]
                mldL=_z[mld_idx]
                qcL=1
                
            else:
                mldL=df.depth.iloc[-1]
                qcL=2 #profile doesn't surpass delta_sigma threshold
        else:
            mldL=np.nan
            qcL=5 #one of T S P being nan cuased related density vars to be nan    
    else:
        mldL=np.nan
        qcL=3 #ref depth is nan
    return mldL,qcL




def calc_mld(df,ref_depth,deltaT):
    """
    based off of de Boyer Montegut (2007) and Rudzin (2017)
    df - single glider profile
    ref_depth - reference depth to calculate potential desnity at 
    deltaT - change in temperature threshold

    qc: 1 - profile is good and has mld
        2 - profile doesn't surpass delta_signma threshold
        3 - profile doesn't have a depth close enough to ref_depth
        4 - profile is less than 10m long
        5 - T, S, P being nan cuased related density vars to be nan
        6 - mldL < mldU, relatively unstratified
    """
    ttdf=df.sort_values('depth',ascending=True)
    tdf=depth_interpolate(ttdf,depth_var='depth',depth_min=ttdf.depth.min(),depth_max=ttdf.depth.max(),stride=0.25,method='linear')
    # print(tdf.gliders.iloc[0],tdf.time.iloc[0])
    depth_diff=tdf.depth.iloc[-1]-tdf.depth.iloc[0]
    if (depth_diff > 10) or (np.isnan(depth_diff)):
        _z=tdf.depth.values
        _temp=tdf.temperature.values
        _sal=tdf.salinity.values
        lat=tdf.lat.values
        lon=tdf.lon.values
    
        if (~np.isnan(_sal).all()) and (~np.isnan(_temp).all()) and (~np.isnan(_z).all()):

            mldU,qcU = calc_mld_upper(tdf,ref_depth,deltaT,_z,_temp,_sal,lat,lon)

            mldL,qcL = calc_mld_lower(tdf,deltaT,_z,_temp,_sal,lat,lon)

            if mldL < mldU:
                qcU=6
                qcL=6
                mldU=df.depth.iloc[-1]
                mldL=df.depth.iloc[-1]

        else:
            mldU=np.nan
            mldL=np.nan
            qcU=5 # one of T S P being nan cuased related density vars to be nan
            qcL=5

    else:
        mldU=np.nan
        mldL=np.nan
        qcU=4 #profile less than 10m
        qcL=4
    return mldU, mldL, qcU, qcL



## Get gliders from erddap
def get_active_gliders(bbox=None, t0=None, t1=dt.date.today(), variables=None, 
                       timeout=5, parallel=False):
    variables = variables or ['time', 'latitude', 'longitude','pressure','depth','temperature','conductivity']
    bbox = bbox or [-100, -40, 18, 60]
    t0 = t0 or (t1 - dt.timedelta(days=1))

    # Convert dates to strings
    t0 = t0.strftime('%Y-%m-%dT%H:%M:%SZ')
    t1 = t1.strftime('%Y-%m-%dT%H:%M:%SZ')

    # Initialize GliderDAC Object
    # e = ERDDAP(server='NGDAC')
    # e = ERDDAP(server="https://gliders.ioos.us/erddap")
    e = ERDDAP(server='http://slocum-data.marine.rutgers.edu//erddap')
    

    # Set timeout (seconds)
    e.requests_kwargs['timeout'] = timeout

    # Grab every dataset available
    # Search constraints
    kw = dict()
    kw['min_time'] = t0
    kw['max_time'] = t1

    if bbox:
        kw['min_lon'] = bbox[0]
        kw['max_lon'] = bbox[1]
        kw['min_lat'] = bbox[2]
        kw['max_lat'] = bbox[3]

    search_url = e.get_search_url(search_for=None, response='csv', **kw)

    try:
        # Grab the results
        search = pd.read_csv(search_url)
    except uHTTPError as error:
        print(f"{inspect.currentframe().f_code.co_name} - Error: {error}")
        # return empty dataframe if there are no results
        return pd.DataFrame()
    except URLError as e:
        print(f"{inspect.currentframe().f_code.co_name} - Error: {error}")
        # return empty dataframe if there are no results
        return pd.DataFrame()

    # Extract the IDs
    gliders = search['Dataset ID'].values

    msg = f"Found {len(gliders)} Glider Datasets: "
    pprint(msg + ', '.join(gliders.tolist()))

    # Setting constraints
    constraints = {
            'time>=': t0,
            'time<=': t1,
            # 'longitude>=': bbox[0],
            # 'longitude<=': bbox[1],
            # 'latitude>=': bbox[2],
            # 'latitude<=': bbox[3],
            }
    def find_duplicates_keep_delayed(df,server='rutgers'):
        """
        Find where deployment names are duplicated meaning real time and delayed datasets are returned.
        If so, keep delayed mode dataset
        server is rutgers or ioos
        """
        
        
        if server == 'ioos':
            #get unique deployment names
            deploy_names=np.unique(df.gliders)
            #strip '-delayed' from name
            strip_names=pd.Series(deploy_names).str.strip('-delayed').values
        elif server == 'rutgers':
            df=df.iloc[np.where(df.gliders.str.contains('-profile-sci'))[0]]
            #get unique deployment names
            deploy_names=np.unique(df.gliders)
            strip_names=pd.Series(df.gliders).str.replace('-profile-sci-rt','').values
            strip_names=pd.Series(strip_names).str.replace('-profile-sci-delayed','').values
        namedf=pd.DataFrame({'orig':deploy_names,'strip':strip_names,'dup':np.zeros(len(deploy_names))})
        #find where original name and stripped name are duplicated, set dup to 1 to indicated true
        namedf.dup[namedf.duplicated(subset=['strip'],keep=False)]=1
        #reset dup value back to zero to keep the delayed dataset
        namedf.dup.iloc[np.where((pd.Series(namedf.orig).str.contains('-delayed'))&(namedf.dup==1))[0]]=0
        # get names of deployments that we want to keep, dup==0
        n=namedf[namedf.dup==0].orig.values
        #subset original dataframe by list of names we want to keep
        df=df[df['gliders'].isin(n)]
        return df
    dfn=pd.DataFrame({'gliders':gliders})
    dfn=find_duplicates_keep_delayed(dfn,server='ioos')
    def request_multi(dataset_id, protocol="tabledap", variables=None):
   
        e.constraints = constraints
        e.protocol = protocol
        # e.variables = ['time', 'latitude', 'longitude']
        e.variables = variables
        e.dataset_id = dataset_id
        print(dataset_id)
        
        # Drop units in the first line and Nans
        
        df = e.to_pandas(
            response="csv", 
            index_col="time",
            parse_dates=True,
            skiprows=(1,)
            ).dropna().tz_localize(None)
        return (dataset_id, df)
    dfs = {glider: df for (glider, df) in [request_multi(id) for id in dfn.gliders.values]}
    df = pd.concat(dfs)
    return df.DataFrame()
    return df

def remove_profiles(pdf,thresh=3):
    return np.nanmin(pdf.sort_values('depth',ascending=True).depth.values) <=thresh



def mld_stats(mld_df,df):
    std_mldu = np.nanstd(mld_df.mld_upper.values)
    std_mldl = np.nanstd(mld_df.mld_lower.values)
    
    avg_mldu = np.nanmean(mld_df.mld_upper.values)
    avg_mldl = np.nanmean(mld_df.mld_lower.values)
    
    med_mldu = np.nanmedian(mld_df.mld_upper.values)
    med_mldl = np.nanmedian(mld_df.mld_lower.values)

    # mld_stat_df = pd.DataFrame(data={'std_u':std_u,'std_l':std_l,'avg_u':avg_u,'avg_l':avg_l,'med_u':med_u,'med_l':med_l},index=[0])
    

    threshu=avg_mldu

    std_densu = np.nanstd(df[df.depth<=threshu].potential_density.values)
    avg_densu = np.nanmean(df[df.depth<=threshu].potential_density.values)
    med_densu = np.nanmedian(df[df.depth<=threshu].potential_density.values)


    threshl=avg_mldl

    std_densl = np.nanstd(df[df.depth>=threshl].potential_density.values)
    avg_densl = np.nanmean(df[df.depth>=threshl].potential_density.values)
    med_densl = np.nanmedian(df[df.depth>=threshl].potential_density.values)
    
    upper = {'avg_mldu':avg_mldu,'med_mldu':med_mldu,'std_mldu':std_mldu,'avg_densu':avg_densu,'med_densu':med_densu,'std_densu':std_densu}
    lower = {'avg_mldl':avg_mldl,'med_mldl':med_mldl,'std_mldl':std_mldl,'avg_densl':avg_densl,'med_densl':med_densl,'std_densl':std_densl}
    print(upper)
    print(lower)

    std_pea=np.nanstd(mld_df.pea.values)
    avg_pea=np.nanmean(mld_df.pea.values)
    med_pea=np.nanmedian(mld_df.pea.values)

    pp = {'avg':avg_pea,'med':med_pea,'std':std_pea}
    print(pp)
    return upper,lower,pp


def create_dens_prof(df,stat_up,stat_low,dmin,dmax,stride,vtype):
    """
    df - one of the "stuff" dataframes - setup_stuff
    stat_up,stat_low stats dict for given period
    stride - depth increment 
    vtype - avg, med, std to get value from stat_up,stat_low
    """
    # dmin = dmin or np.round(df.surf_depth.mean())
    # dmax = dmax or np.round(df.bot_depth.mean())
    
    z=np.arange(dmin,dmax+stride,stride)
    
    
    upper_mld= np.round(stat_up[vtype+'_mldu'],0)
    lower_mld = np.round(stat_low[vtype+'_mldl'],0)
    upper_den = np.round(stat_up[vtype+'_densu'],0)
    lower_den = np.round(stat_low[vtype+'_densl'],0)
    
    idxU = np.abs(z - upper_mld).argmin()
    idxL = np.abs(z - lower_mld).argmin()
    # Get the closest value
    mldU= z[idxU]
    mldL= z[idxL]
    p_thick = mldL-mldU
    
    
    
    dens=np.array([np.nan]*len(z))
    
    #set upper layer to upper dens
    dens[:idxU+1]=upper_den
    #set lower layer to lower dens
    dens[idxL:]=lower_den
    
    
    
    # Indices where dens is not NaN
    not_nan_indices = np.where(~np.isnan(dens))[0]
    
    # Values in dens that are not NaN
    not_nan_values = dens[not_nan_indices]
    
    # Linear interpolation to fill NaNs
    interpolated_values = np.interp(np.arange(len(dens)), not_nan_indices, not_nan_values)
    
    # Replace NaNs with interpolated values
    dens = np.where(np.isnan(dens), interpolated_values, dens)
    
    # pea = phi_carpenter(dens,z,stride,dmax)
    depth=pd.Series(z)
    dens = pd.Series(dens)
    z_thresh=dmax
    # pea = potential_energy_anomaly100(depth, dens,z_thresh)
    pea = phi_carpenter(dens,z,stride,dmax)

    return mldU,mldL,p_thick,dens,z,pea
