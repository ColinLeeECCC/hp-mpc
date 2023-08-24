import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from matplotlib.offsetbox import AnchoredText
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import numpy as np
from math import ceil

import xarray as xr

def main(filename='map_cluster_2017080100_2018103100_N2_02_0.80.nc4', outfile=None, rval=None, spec=None):
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

    # ax.set_extent([-111.2, -109.8,
    #               57.25, 58.25], crs=ccrs.Mercator())

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

    clusters = ds.clusters.data.astype('int64') - 1
    numClust=int(np.max(ds.clusters)) + 1
    # Assign a colour index from 0-19 to each cluster
    cl_cm = np.arange(numClust)
    cl_cm = np.remainder( cl_cm, 20 )

    while True:
        colidx = np.take(cl_cm, clusters)
        # we're using a repeated, random, 20-colour colourmap, which means
        # if we have any adjacent (in lat-lon space) clusters with IDs that
        # differ by a multiple of 20, they get the same colour and look like
        # one cluster. So let's try to detect that situation.
        # not_same_idx = np.logical_not(
        #     np.logical_and(
        #         np.equal( clusters[:-1,:-1], clusters[:-1,1:] ),
        #         np.equal( clusters[1:,:-1],  clusters[1:,1:] )
        #     )
        # )
        # same_col = np.logical_and(
        #     np.equal( colidx[:-1,:-1],  colidx[:-1,1:] ),
        #     np.equal( colidx[1:,:-1], colidx[1:,1:] )
        # )


        # tmp = clusters[1:,1:]
        # print('tmp.shape: ', tmp.shape)
        # print('not_same_idx.shape: ', not_same_idx.shape)
        # print('same_col.shape: ', same_col.shape)
        # blenders= np.unique(tmp[np.logical_and(not_same_idx, same_col)])
        mismatch_left = np.logical_and(
            np.logical_not( np.equal( clusters[:-1,:-1], clusters[:-1, 1:] ) ),
            np.equal( colidx[:-1,:-1], colidx[:-1, 1:] )
        )
        mismatch_down = np.logical_and(
            np.logical_not( np.equal( clusters[:-1,:-1], clusters[1:, :-1] ) ),
            np.equal( colidx[:-1,:-1], colidx[1:, :-1] )
        )
        tmp = clusters[:-1,:-1]
        blenders = np.unique(tmp[np.logical_or(mismatch_left, mismatch_down)])
        
        if len(blenders) == 0:
            break
        print('These colours clash: ',blenders)
        for i,b in enumerate(blenders):
            # if i % 2 == 1:
            #     continue
            cnew = np.floor(np.random.rand(1) * 20)
            while cnew == cl_cm[b]:
                cnew = np.floor(np.random.rand(1) * 20)
            print('Swapping {} from {} to {}'.format(b, cl_cm[b], cnew))
            cl_cm[b] = cnew
    
    # plt.contourf(ds.lon, ds.lat, ds.clusters, transform=ccrs.PlateCarree())
    mult = ceil(float(numClust) / 20.0)
    print('numClust = {:}, mult = {:}'.format(numClust, mult))
    tmpcm = ListedColormap(np.tile(cm.get_cmap('tab20').colors, (mult,1))[0:numClust,:])
    pcm = ax.pcolormesh(ds.lon, ds.lat, colidx,
                        transform=ccrs.PlateCarree(),
                        cmap=tmpcm)
    
    # Add a text annotation for the license information to the
    # the bottom right corner.
    # text = AnchoredText('\u00A9 {}; license: {}'
    #                     ''.format(SOURCE, LICENSE),
    #                     loc=4, prop={'size': 12}, frameon=True)
    # ax.add_artist(text)
    
    # # mark a known place to help us geo-locate ourselves
    # ax.plot(-113.48244, 53.56462, 'bo', markersize=4, transform=ccrs.Geodetic())
    # ax.text(-113.1, 53.7, 'Edmonton', transform=ccrs.Geodetic())
    # ax.plot(-114.055482, 51.069558, 'bo', markersize=4, transform=ccrs.Geodetic())
    # ax.text(-113.8, 51.2, 'Calgary', transform=ccrs.Geodetic())

    # # Mark the current location of Mildred lake and the 2 alternatives
    # ax.plot(-111.564147, 57.050006, marker='o', markersize=10,
    #         markerfacecolor='green', color='white',
    #         transform=ccrs.Geodetic())
    # ax.text(-111.557, 57.050006, 'Mildred Lake', transform=ccrs.Geodetic())
    # ax.plot(-111.57798, 57.045828, marker='o', markersize=10,
    #         markerfacecolor='blue', color='white',
    #         transform=ccrs.Geodetic())
    # ax.text(-111.9, 57.040, 'ML alt', transform=ccrs.Geodetic())
    # ax.plot(-111.61209, 57.105280, marker='o', markersize=10,
    #         markerfacecolor='red', color='white',
    #         transform=ccrs.Geodetic())
    # ax.text(-111.62, 57.105, 'Poplar Creek', transform=ccrs.Geodetic())
    # Mark the current location of Poplar Creek and the 2 alternatives
    # ax.plot(-111.61209, 57.105280, marker='o', markersize=10,
    #         markerfacecolor='green', color='white',
    #         transform=ccrs.Geodetic())
    # ax.text(-111.62, 57.105, 'Poplar Creek', transform=ccrs.Geodetic())
    # ax.plot(-111.61132, 57.13059, marker='o', markersize=10,
    #         markerfacecolor='green', color='white',
    #         transform=ccrs.Geodetic())
    # ax.text(-111.61132, 57.13059, 'Option 1', transform=ccrs.Geodetic())
    # ax.plot(-111.61662, 57.12364, marker='o', markersize=10,
    #         markerfacecolor='blue', color='white',
    #         transform=ccrs.Geodetic())
    # ax.text(-111.61662, 57.12364, 'Option 2', transform=ccrs.Geodetic())
    # # Mark the current location of Fort McKay South and the alternative
    # ax.plot(-111.642344, 57.149175, marker='o', markersize=10,
    #         markerfacecolor='green', color='white',
    #         transform=ccrs.Geodetic())
    # ax.text(-111.742344, 57.149175, 'Fort McKay South', transform=ccrs.Geodetic())
    # ax.plot(-111.63349, 57.14138, marker='o', markersize=10,
    #         markerfacecolor='blue', color='white',
    #         transform=ccrs.Geodetic())
    # ax.text(-111.63349, 57.14138, 'FMS Alt', transform=ccrs.Geodetic())
    # Mark the current location of Firebag and Wapasu
    ax.plot(-110.89799, 57.239525, marker='o', markersize=10,
            markerfacecolor='green', color='white',
            transform=ccrs.Geodetic())
    ax.text(-110.89799, 57.239525, 'Firebag', transform=ccrs.Geodetic())
    ax.plot(-111.0385386, 57.25913818, marker='o', markersize=10,
            markerfacecolor='blue', color='white',
            transform=ccrs.Geodetic())
    ax.text(-111.0385386, 57.25913818, 'Wapasu', transform=ccrs.Geodetic())
    # # Mark the current location of all the southern regional sites
    # ax.plot(-110.749758, 55.903286, marker='o', markersize=10,
    #         markerfacecolor='green', color='white',
    #         transform=ccrs.Geodetic())
    # ax.plot(-111.440343, 55.812883, marker='o', markersize=10,
    #         markerfacecolor='green', color='white',
    #         transform=ccrs.Geodetic())
    # ax.plot(-111.172798, 55.621487, marker='o', markersize=10,
    #         markerfacecolor='green', color='white',
    #         transform=ccrs.Geodetic())
    # ax.plot(-111.078877, 55.632330, marker='o', markersize=10,
    #         markerfacecolor='green', color='white',
    #         transform=ccrs.Geodetic())
    # ax.plot(-110.705815, 55.667804, marker='o', markersize=10,
    #         markerfacecolor='green', color='white',
    #         transform=ccrs.Geodetic())
    # ax.plot(-110.876009, 55.579153, marker='o', markersize=10,
    #         markerfacecolor='green', color='white',
    #         transform=ccrs.Geodetic())
    # ax.plot(-110.9759839, 55.51870937, marker='o', markersize=10,
    #         markerfacecolor='green', color='white',
    #         transform=ccrs.Geodetic())
    # ax.plot(-110.865308, 55.523807, marker='o', markersize=10,
    #         markerfacecolor='green', color='white',
    #         transform=ccrs.Geodetic())
    # ax.plot(-111.0139229, 55.343734, marker='o', markersize=10,
    #         markerfacecolor='green', color='white',
    #         transform=ccrs.Geodetic())
    # ax.plot(-111.218823, 55.461647, marker='o', markersize=10,
    #         markerfacecolor='green', color='white',
    #         transform=ccrs.Geodetic())

    reader = shpreader.Reader('/space/hall5/sitestore/eccc/aq/r1/juz001/cnfs_pxarqdr2/Shapefile/OS_SubLayers_w_Org/OS_Cleared_D.shp')
    # There is no way to automatically load the projection file with a shapefile
    # so I am doing this manually and harcoding it >:*|
    tm = ccrs.TransverseMercator(false_easting=500000.0, false_northing=0.0, central_longitude=-115.0, scale_factor=0.9992, central_latitude=0.0)
    
    os_ops = cfeature.ShapelyFeature(reader.geometries(), tm, edgecolor='red', facecolor='none', alpha=0.5)
    ax.add_feature(os_ops)

    ax.set_extent([-111.4, -110.2,
                   57.0, 57.5])

    textstr = ""

    if spec is not None:
        textstr = textstr + f"{spec}\n"
        
    textstr = textstr + f"N = {numClust}"

    if rval is not None:
        textstr = textstr + f"\nR \u2265 {rval}"
        
    text = AnchoredText(textstr,
                        loc=1, prop={'size': 12}, frameon=True)
    ax.add_artist(text)
    
    # cb = fig.colorbar(pcm, ax=ax)
    plt.tight_layout()
    if outfile is None:
        plt.show()
    else:
        plt.savefig(outfile, bbox_inches='tight')


if __name__ == '__main__':
    import sys
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('outfile', nargs='?', default=None)
    parser.add_argument('-r', '--rval', default=None)
    parser.add_argument('-s', '--species', default=None)

    args = parser.parse_args()

    main(filename=args.infile, outfile=args.outfile, rval=args.rval, spec=args.species)
         
    # if len(sys.argv) == 1:
    #     main()
    # elif len(sys.argv) == 2:
    #     main(filename=sys.argv[1])
    # else:
    #     main(filename=sys.argv[1], outfile=sys.argv[2])
