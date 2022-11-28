import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from matplotlib.offsetbox import AnchoredText
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import numpy as np

import xarray as xr

def main(filename='map_cluster_2017080100_2018103100_N2_02_0.80.nc4', outfile=None):
    ds = xr.open_dataset(filename)
    # ds = xr.open_dataset('map_cluster_2017080100_2018103100_S2_02_0.85.nc4')
    midrlon = int(len(ds.rlon)/2)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.AlbersEqualArea(
        ds.lon.data[0,midrlon], ds.lat.data[0,midrlon]
    ))

    # ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mercator())

    # ax.set_extent([np.min(ds.lon.data), np.max(ds.lon.data),
    #                np.min(ds.lat.data), np.max(ds.lat.data)],
    #               crs=ccrs.Mercator())

    # Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
        facecolor='none')

    SOURCE = 'Natural Earth'
    LICENSE = 'public domain'

    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.LAKES)
    ax.add_feature(cfeature.RIVERS)
    ax.add_feature(states_provinces)

    # plt.contourf(ds.lon, ds.lat, ds.clusters, transform=ccrs.PlateCarree())
    tmpcm = ListedColormap(np.tile(cm.get_cmap('tab20').colors, (5,1)))
    pcm = ax.pcolormesh(ds.lon, ds.lat, ds.clusters,
                        transform=ccrs.PlateCarree(),
                        cmap=tmpcm)
    
    # Add a text annotation for the license information to the
    # the bottom right corner.
    # text = AnchoredText('\u00A9 {}; license: {}'
    #                     ''.format(SOURCE, LICENSE),
    #                     loc=4, prop={'size': 12}, frameon=True)
    # ax.add_artist(text)
    
    # mark a known place to help us geo-locate ourselves
    ax.plot(-113.48244, 53.56462, 'bo', markersize=4, transform=ccrs.Geodetic())
    ax.text(-113.1, 53.7, 'Edmonton', transform=ccrs.Geodetic())
    ax.plot(-114.055482, 51.069558, 'bo', markersize=4, transform=ccrs.Geodetic())
    ax.text(-113.8, 51.2, 'Calgary', transform=ccrs.Geodetic())

    reader = shpreader.Reader('/space/hall5/sitestore/eccc/aq/r1/juz001/cnfs_pxarqdr2/Shapefile/OS_SubLayers_w_Org/OS_Cleared_D.shp')
    # There is no way to automatically load the projection file with a shapefile
    # so I am doing this manually and harcoding it >:*|
    tm = ccrs.TransverseMercator(false_easting=500000.0, false_northing=0.0, central_longitude=-115.0, scale_factor=0.9992, central_latitude=0.0)
    
    os_ops = cfeature.ShapelyFeature(reader.geometries(), tm, edgecolor='red')
    ax.add_feature(os_ops)

    cb = fig.colorbar(pcm, ax=ax)
    if outfile is None:
        plt.show()
    else:
        plt.savefig(outfile, bbox_inches='tight')


if __name__ == '__main__':
    import sys
    if len(sys.argv) == 1:
        main()
    elif len(sys.argv) == 2:
        main(filename=sys.argv[1])
    else:
        main(filename=sys.argv[1], outfile=sys.argv[2])
