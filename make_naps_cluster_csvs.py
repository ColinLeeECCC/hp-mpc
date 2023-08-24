#!/usr/bin/python

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
    parser.add_argument("num_points",        type=int,
                        help="The total number of points in the file")
    parser.add_argument("-t", "--threshold", type=float,
                        default=9999.999,
                        help="The dissimilarity score cutoff")
    parser.add_argument("-n", "--num-clusters", type=int,
                        default=-1,
                        help="The final number of clusters you would like")
    parser.add_argument("-g", "--debug", action='store_true',
                        help="Turn on debug printing")

    args = parser.parse_args()

    import numpy as np
    
    num_clusters = args.num_points
    
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

    for i,cluster_id in enumerate(np.unique(cluster_idx)):
        cluster_idx[cluster_idx == cluster_id] = i
        
    for id in cluster_idx:
        print(id)
