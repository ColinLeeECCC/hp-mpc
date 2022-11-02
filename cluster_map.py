#!/bin/python


if __name__ == "__main__":
    import argparse
    import glob
    import xarray as xr
    from os import path

    desc="Display cluster map"
    usage = """
    %(prog)s cluster.dat xsize ysize [-t|--threshold f] [-n|--num-clusters n]
         Take the cluster.dat input and combine clusters until the 
         dissimilarity threshold is reached or until there are only 
         n clusters remaining
    """
    epilog = "Please specificy exactly one of threshold or num-clusters"

    def file_exists(fn):
        return path.exists(fn)
    
    parser = argparse.ArgumentParser(
        description=desc, usage=usage, epilog=epilog,
        prefix_chars='-+', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("cluster_filename",  help="Cluster output file")
    parser.add_argument("nc_filename", help="One of the netcdf files used to create the clusters")
    parser.add_argument("xsize",     type=int,
                        help="The X side of the display grid")
    parser.add_argument("ysize",     type=int,
                        help="The Y side of the display grid")
    parser.add_argument("-t", "--threshold", type=float,
                        default=9999.999,
                        help="The dissimilarity score cutoff")
    parser.add_argument("-n", "--num-clusters", type=int,
                        default=-1,
                        help="The final number of clusters you would like")
    parser.add_argument("-g", "--debug", action='store_true',
                        help="Turn on debug printing")
    parser.add_argument("--start", help="The starting X and Y for extracting")

                        

    args = parser.parse_args()

    if args.start is not None:
        start_x = int(args.start.split(',')[0])
        start_y = int(args.start.split(',')[1])

    if (args.threshold == 9999.999 and args.num_clusters == -1) or \
       (args.threshold != 9999.999 and args.num_clusters != -1):
        print("You must specificy a dissimilarity threshold or a"
              " final number of clusters, but not both.")
        exit(-1)


    import numpy as np
    from math import floor

    num_clusters = args.xsize * args.ysize
    clusters = []
    cluster_idx = np.arange(num_clusters)
    # Create a list of lists where each cluster
    # just contains its own number initially
    for i in range(num_clusters):
        clusters.append([i])

    with open(args.cluster_filename, newline = '') as cl_file:
        headers = cl_file.readline()

        i = 0
        d = -999
        while d < args.threshold and \
              num_clusters > args.num_clusters:
            i    = i+1
            line = next(cl_file)
            
            pair = line.split()
            d  = float(pair[0])
            ik = int(pair[1])-1
            jk = int(pair[2])-1

            # If cluster ik doesn't contain itself,
            # swap ik and jk becuase we always want
            # to add to the cluster that contains
            # itself
            if (cluster_idx[ik] != ik):
                if args.debug:
                    print("{} belongs to {}... swapping".format(
                        ik, cluster_idx[ik]
                    ))
                ik = cluster_idx[ik]
            if (cluster_idx[jk] != jk):
                if args.debug:
                    print("{} belongs to {}... swapping".format(
                        jk, cluster_idx[jk]
                    ))
                jk = cluster_idx[jk]

            if (cluster_idx[ik] != ik):
                print("Something is wrong. Ik does not contain itself.")
                exit(-2)
            if (cluster_idx[jk] != jk):
                print("Something is wrong. Jk does not contain itself.")
                exit(-2)
            if args.debug:
                print("Adding cluster {}({}) to cluster {}({})".format(
                    ik, len(clusters[jk]), jk, len(clusters[ik])))

            tmp = clusters[jk]
            for lk in tmp:
                cluster_idx[lk] = ik
                
            clusters[ik].extend(tmp)
            clusters[jk] = []
            num_clusters = num_clusters - 1
            if args.debug:
                print("clusters[ik] now has {} points".format(len(clusters[ik])))
        print("Clustering ended with dissimilarity score of {:f} and {:d} clusters remaining".format(d,num_clusters))
        # Okay, now we need to sort the remaining clusters so that the most similar ones are next to one another
        # Read rest of file
        for line in cl_file:
            pair = line.split()
            d  = float(pair[0])
            ik = int(pair[1])-1
            jk = int(pair[2])-1
            # If cluster ik doesn't contain itself,
            # swap ik and jk becuase we always want
            # to add to the cluster that contains
            # itself
            if (cluster_idx[ik] != ik):
                if args.debug:
                    print("{} belongs to {}... swapping".format(
                        ik, cluster_idx[ik]
                    ))
                ik = cluster_idx[ik]
            if (cluster_idx[jk] != jk):
                if args.debug:
                    print("{} belongs to {}... swapping".format(
                        jk, cluster_idx[jk]
                    ))
                jk = cluster_idx[jk]

            if (cluster_idx[ik] != ik):
                print("Something is wrong. Ik does not contain itself.")
                exit(-2)
            if (cluster_idx[jk] != jk):
                print("Something is wrong. Jk does not contain itself.")
                exit(-2)
            if args.debug:
                print("Matching cluster {}({}) with cluster {}({})".format(
                    ik, len(clusters[jk]), jk, len(clusters[ik])))
            

    dims = (args.ysize, args.xsize)
    ## This may be the wrong order, may need to switch x and y sizes
    clmap = np.full(dims, np.nan)
    # Grab the remaining clusters
    for clidx, cl in enumerate(filter(lambda x: len(x) > 0, clusters)):
        for i in cl:
            gridi = floor(float(i) / float(dims[1]))
            gridj = i - gridi * dims[1]
            clmap[gridi,gridj] = clidx+1
    
    geodata = xr.open_dataset(args.nc_filename)
    ## There are different netcdf inputs used so we need to check what ours looks like
    #  by examing the metadata
    if 'rlat1' in geodata.dims.keys():
        coords={
            "rlat": geodata.rlat1.data[start_x:(start_x+dims[0])],
            "lat":  (('rlat', 'rlon'), geodata.lat_1.data[start_x:(start_x+dims[0]), start_y:(start_y+dims[1])]),
            "rlon": geodata.rlon1.data[start_y:(start_y+dims[1])],
            "lon":  (('rlat', 'rlon'), geodata.lon_1.data[start_x:(start_x+dims[0]), start_y:(start_y+dims[1])])
                 }
    elif 'lat' in geodata.dims.keys():
        # Colin's test netcdf files
        coords={
                 "lat": geodata.lat[start_x:(start_x+dims[0])],
                 "lon": geodata.lon[start_y:(start_y+dims[1])]
                 }

    clmap_ds = xr.Dataset(
        {'clusters': (("rlat", "rlon"), clmap, {'grid_mapping':"rotated_pole"} )},
        coords=coords
    )
    clmap_ds.to_netcdf("map_cluster.nc4")
