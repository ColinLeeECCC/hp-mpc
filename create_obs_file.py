#!/usr/bin/python3

import argparse
import csv
import numpy as np
import sys

parser = argparse.ArgumentParser(description='Take a dendrograma and aircraft geolocation data to create an obs file with clusters.')
parser.add_argument('N', type=int, nargs='+',
                    help='Number of final clusters')
parser.add_argument('dendrogram_file', 
                    help='File containing dendrogram data.')

parser.add_argument('geolocation_file', 
                    help='File containing geolocation data.')


args = parser.parse_args()


time = []
lat  = []
lon  = []
alt  = []
mask = []
# Open the geolocation file
with open(args.geolocation_file, newline='') as geofile :
    geo_reader = csv.reader(geofile)
    for row in geo_reader:
        if (geo_reader.line_num > 1):
            time.append(float(row[0]))
            lat.append( float(row[1]))
            lon.append( float(row[2]))
            alt.append( float(row[3]))
            mask.append( bool(int(row[4])))

nreal_clusters = np.sum(mask)
nbad_clusters = len(mask) - nreal_clusters

# print("There are {} valid clusters ({} bad)".format(nreal_clusters, len(mask) - nreal_clusters))

d = []
lki = []
lkj = []

# Open the dendogram file
with open(args.dendrogram_file, newline='') as dendfile:
    dend_reader = csv.reader(dendfile, delimiter='\t')
    for row in dend_reader:
        if (dend_reader.line_num > 1):
            # print(float(row[0]), int(row[1]), int(row[2]))
            d.append(float(row[0]))
            lki.append(int(row[1]))
            lkj.append(int(row[2]))

# print('Read in {} rows.'.format(len(d)))

num_clusters = len(d)
idx = np.repeat(np.array([[x for x in range(1,num_clusters+2)]]), len(args.N), axis=0)

# For now, since the clustering program isn't handling bad data
# in any special way, we need to watch in case those bad rows
# get clustered. They *should* only be clustered at the *end*
# but just in case...

bg_cluster=[-1]*len(args.N)

for k, M in enumerate(args.N):
    lists = [[x] for x in idx[k,:]]
    for i in range(nreal_clusters - M):
        ii = lki[i] - 1
        jj = lkj[i] - 1

        if not ( mask[ii] and mask[jj] ):
            print('Step {:05d}: Clustering {} and {} but mask = {},{}'.format(
                i, ii+1, jj+1, mask[ii], mask[jj]), file=sys.stderr)

        l = lists[jj]
        for j in l:
            idx[k,j-1] = ii + 1
            lists[ii].append(j)

    final_clusters, cluster_sizes = np.unique(np.array(idx[k,:])[mask], return_counts=True)
    # print('There are {} final clusters and M = {}'.format(len(final_clusters), M))
    bg_cluster[k] = final_clusters[np.argmax(cluster_sizes)]

# # let's reorder the remaining clusters for relatedness. This is a first
# # stab at this but I'm not confident the method is correct
# for i in range(nreal_clusters - args.N, nreal_clusters + 1):
#     ii = lki[i] - 1
#     jj = lkj[i] - 1
    
#     if not ( mask[ii] and mask[jj] ):
#         print('Step {:05d}: Clustering {} and {} but mask = {},{}'.format(
#             i, ii+1, jj+1, mask[ii], mask[jj]), file=sys.stderr)

#     l = lists[jj]
#     for j in l:
#         idx[j-1] = ii + 1
#         lists[ii].append(j)

print('Obs 3.1')
print('Aircraft clusters')
print('ID\tLAT\tLON\tELEV', end='')
for k,M in enumerate(args.N):
    print('\tDATA.CLUSTER_{}'.format(M), end='')
print('')

for i in range(num_clusters):
    if (lat[i] < -90 or not mask[i]):
        continue
    print('CLUSTER\t{:.04f}\t{:.04f}\t{}'.format(lat[i], lon[i], int(alt[i])), end='')
    for k,M in enumerate(args.N):
        if ((idx[k,i] != bg_cluster[k]) or ((i-1) % 60 == 0) ):
            print('\t{}'.format( idx[k,i] ), end='')
        else:
            print('\tNaN', end='')
    print('')

    
