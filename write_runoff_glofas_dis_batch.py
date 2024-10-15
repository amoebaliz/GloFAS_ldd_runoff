import numpy as np
import xarray
import xesmf
import sys
import copy
import warnings
warnings.filterwarnings("ignore")

def update_stencil_sum(ocean_mask):
    
    # Make 3x3 stensil 
    sten_nlat = 3
    sten_nlon = 3

    glofas_stensil_sum = np.zeros(ocean_mask.shape)
    for nlat in range(sten_nlat):
        ilat_last = -(2-nlat)
        if ilat_last == 0:
            ilat_last = None
        
        for nlon in range(sten_nlon):
            ilon_last = -(2-nlon)
            if ilon_last == 0:
                ilon_last = None
       
            glofas_stensil_sum[1:-1,1:-1] += ocean_mask[nlat:ilat_last,nlon:ilon_last]
        
    return glofas_stensil_sum 

def get_glofas_pour_points():
    ldd_modified = copy.deepcopy(ldd)
    # GloFAS VERSION 4: editing pour points along map seam in Russia Chukotka Region at date line.
    lat_idx = np.where(np.round(glofas_lat,3) == 68.875)[0][0] 
    lon_idx = np.where(np.round(glofas_lon,3) == 180.025)[0][0]
    ldd_modified[lat_idx,lon_idx]=1 # change from 5 to stop advancing halo point 
    
    glofas_ocean_mask = np.isnan(ldd)

    # initialize the while loop
    n_updates = 1

    # continue 
    while n_updates>0:
        n_updates = 0
        ocean_stencil_sum = update_stencil_sum(glofas_ocean_mask)    
        for jj in range(len(glofas_lat)):
            for ii in range(len(glofas_lon)):
            
                if (ldd_modified[jj,ii]==5) and (glofas_ocean_mask[jj,ii]==0) and (ocean_stencil_sum[jj,ii]>0):
                    # update ocean mask so counted in next calculation of the stencil sum
                    glofas_ocean_mask[jj,ii] = 1
                    n_updates+=1
                
    pour_points = (ldd_modified==5)*(ocean_stencil_sum>0)
    
    return pour_points

def get_mom_mask_for_glofas():
    # eliminate pour points outside of nep domain so they won't be mapped on the boundary in the d2s phase
    mom_to_glofas = xesmf.Regridder(
        {'lat': lat, 'lon': lon, 'lat_b': latb, 'lon_b': lonb},
        {'lon': glofas_lon, 'lat': glofas_lat, 'lon_b': glofas_lonb , 'lat_b': glofas_latb },
        method='conservative',
        periodic=True,
        reuse_weights=False)
    mom_mask_for_glofas = mom_to_glofas(np.ones((ocn_mask.shape)))
    
    mom_mask_for_glofas[mom_mask_for_glofas<.5]=0
    mom_mask_for_glofas[mom_mask_for_glofas>.5]=1
    
    return mom_mask_for_glofas.astype(bool)


def write_runoff(glofas, hgrid, coast_mask, out_file):
    # From Alistair
    area = (hgrid.area[::2, ::2] + hgrid.area[1::2, 1::2]) + (hgrid.area[1::2, ::2] + hgrid.area[::2, 1::2])

    pour_points = get_glofas_pour_points()
    glofas_mom_pour_points = pour_points*get_mom_mask_for_glofas()

    # Identify nearest coastal land point to GloFAS pour points
    # Source pour points indices
    flat_pour_point_mask = glofas_mom_pour_points.ravel().astype('bool')
    pour_lon = glo_lons.ravel()[flat_pour_point_mask]
    pour_lat = glo_lats.ravel()[flat_pour_point_mask]
    #linear index of glofas_elements
    glo_id = np.arange(np.prod(glofas_mom_pour_points.shape))
    pour_id = glo_id[flat_pour_point_mask]

    # Coastal destination indices
    # Flatten mask and coordinates to 1D
    flat_coast_mask = coast_mask.ravel().astype('bool')
    coast_lon = lon.ravel()[flat_coast_mask]
    coast_lat = lat.ravel()[flat_coast_mask]
    #linear index of mom elements
    mom_id = np.arange(np.prod(coast_mask.shape))

    # Use xesmf to find the index of the nearest coastal cell
    # for every grid cell in the MOM domain
    coast_to_mom = xesmf.Regridder(
       {'lat': coast_lat, 'lon': coast_lon},
       {'lat': pour_lat, 'lon': pour_lon},
       method='nearest_s2d',
       locstream_in=True,
       locstream_out=True,
       reuse_weights=False)

    coast_id = mom_id[flat_coast_mask]
    nearest_coast = coast_to_mom(coast_id).ravel()

    # discharge on GloFAS grid, reshaped to 2D (time, grid_id)
    raw = glofas.values.reshape([glofas.shape[0],-1])

    # Zero array that will be filled with runoff at coastal cells
    mom_filled = np.zeros((raw.shape[0],np.prod(ocn_mask.shape)))

    # Loop over each GloFAS pour point and add it to the nearest coastal cell

    for mom_i, glo_i in zip(nearest_coast,pour_id):
        mom_filled[:, mom_i] += raw[:, glo_i] 

    # Reshape back to 3D and convert to kg m-2 s-1
    filled_reshape = 1000*mom_filled.reshape((glofas.shape[0],ocn_mask.shape[0],ocn_mask.shape[1]))/area.values

    # Convert to xarray
    ds = xarray.Dataset({
        'runoff': (['time', 'y', 'x'], filled_reshape),
        'area': (['y', 'x'], area.data),
        'lat': (['y', 'x'], lat.data),
        'lon': (['y', 'x'], lon.data)
        },
        coords={'time': glofas['time'].data, 'y': np.arange(filled_reshape.shape[1]), 'x': np.arange(filled_reshape.shape[2])}
    )

    # Drop '_FillValue' from all variables when writing out
    all_vars = list(ds.data_vars.keys()) + list(ds.coords.keys())
    encodings = {v: {'_FillValue': None} for v in all_vars}

    # Make sure time has the right units and datatype
    # otherwise it will become an int and MOM will fail. 
    encodings['time'].update({
        'units': 'days since 1950-01-01',
        'dtype': np.float, 
        'calendar': 'gregorian'
    })

    ds['time'].attrs = {'cartesian_axis': 'T'}
    ds['x'].attrs = {'cartesian_axis': 'X'}
    ds['y'].attrs = {'cartesian_axis': 'Y'}
    ds['lat'].attrs = {'units': 'degrees_north'}
    ds['lon'].attrs = {'units': 'degrees_east'}
    ds['runoff'].attrs = {'units': 'kg m-2 s-1'}

    # Write out
    ds.to_netcdf(
        out_file,
        unlimited_dims=['time'],
        format='NETCDF3_64BIT',
        encoding=encodings,
        engine='netcdf4'
    )
    ds.close()


if __name__ == '__main__':
    ## Adapted from code originally developped by Andrew C. Ross 
    # Determine coastal points in NEP domain
    ocn_mask = xarray.open_dataarray(<MOM6 OCEAN MASK FILE STRING>).values.astype(bool)
    stencil_sum = 0 * ocn_mask
    # sum ocean mask values to the:     north          south                 east                   west
    stencil_sum[1:-1,1:-1] = ~ocn_mask[2:,1:-1] + ~ocn_mask[:-2,1:-1] + ~ocn_mask[1:-1,:-2] + ~ocn_mask[1:-1,2:]
    coast = (ocn_mask)*(stencil_sum>0)

    # Load regional ocean hgrid
    hgrid = xarray.open_dataset(<MOM6 OCEAN_HGRID FILE STRING>)
    lon = hgrid.x[1::2, 1::2].values
    lonb = hgrid.x[::2, ::2].values
    lat = hgrid.y[1::2, 1::2].values
    latb = hgrid.y[::2, ::2].values 
     
    # Load GloFAS local drain direction map
    ldd = xarray.open_dataarray(<GLOFAS LDD MAP FILE STRING>).values 

    y=int(sys.argv[1])
    # GloFAS 3.1 subset to model region:
    files = [f'<GOFAS FILE STRING CONTAINING YEAR y>' for y in [y-1, y, y+1]]
    glofas = (
         xarray.open_mfdataset(files, combine='by_coords')
         .rename({'latitude': 'lat', 'longitude': 'lon'})
	 .sel(time=slice(f'{y-1}-12-31 12:00:00', f'{y+1}-01-01 12:00:00'))
         .dis24)
    
    glofas_lat = glofas['lat'].values 
    glofas_lon = glofas['lon'].values
    glofas_lon[glofas_lon<0]=glofas_lon[glofas_lon<0]+360
    deg_incr = abs(np.unique(np.diff(glofas_lat))[0])
    glofas_latb = np.arange(glofas_lat[0] + deg_incr/2., glofas_lat[-1]-deg_incr, -deg_incr) 
    glofas_lonb = np.arange(glofas_lon[0] - deg_incr/2., glofas_lon[-1]+deg_incr, deg_incr)
    glo_lons,glo_lats = np.meshgrid(glofas_lon,glofas_lat) 
    
    out_file = f'<OUTPUT FILE STRING>'
    write_runoff(glofas, hgrid, coast, out_file)


