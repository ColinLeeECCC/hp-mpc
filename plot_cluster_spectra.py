#!/usr/bin/python3

import argparse
import csv
import numpy as np
import sys

parser = argparse.ArgumentParser(description='Take a dendrograma and mass spec data to create plots of the spectrum of each cluster.')
parser.add_argument('N', type=int,
                    help='Number of final clusters')
parser.add_argument('dendrogram_file', 
                    help='File containing dendrogram data.')
parser.add_argument('spectra_file', 
                    help='File containing mass spec data.')
parser.add_argument('--output-filename', 
                    help='File name for saving output figure.')
parser.add_argument('--linkage', type=int, 
                    help='Linkage for title (integer 0-6).')

linkages = ['Single', 'Complete', 'Average', 'Weighted', 'Centroid', 'Median', 'Ward''s']

args = parser.parse_args()

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
idx = np.array([x for x in range(1,num_clusters+2)])

spec = []
mask = np.ones( num_clusters+1, dtype=bool)
# Open the geolocation file
with open(args.spectra_file, newline='') as specfile :
    spec_reader = csv.reader(specfile)
    for row in spec_reader:
        n = spec_reader.line_num - 2
        if n == -1:
            continue
        if n == 0:
            spec = np.zeros([num_clusters+1, len(row)])
            
        spec[n,:] = np.array(row[:], dtype=float)
        if np.any(spec[n,:] == -999999.0):
            print('Row {:d} has NA value = {:g}'.format( n, spec[n, 0] ))
            mask[n] = False
            
nreal_clusters = num_clusters - sum( np.logical_not(mask) ) + 1
# For now, since the clustering program isn't handling bad data
# in any special way, we need to watch in case those bad rows
# get clustered. They *should* only be clustered at the *end*
# but just in case...

lists = [[x] for x in idx]
clustered = 0
i = 0
# for i in range(nreal_clusters - args.N):
while clustered < nreal_clusters - args.N:
    ii = lki[i] - 1
    jj = lkj[i] - 1
    
    i = i + 1
    if not ( mask[ii] and mask[jj] ):
        print('Step {:05d}: Clustering {} and {} but mask = {},{}'.format(
            i-1, ii+1, jj+1, mask[ii], mask[jj]), file=sys.stderr)
        continue
    clustered = clustered + 1
    l = lists[jj]
    for j in l:
        idx[j-1] = ii + 1
        lists[ii].append(j)

final_clusters, cluster_sizes = np.unique(idx[mask], return_counts=True)
if (len(final_clusters) != args.N):
    raise Exception('We didn''t get the number of clusters we were expecting.')

cluster_spectra = np.zeros([args.N, spec.shape[1]])
cluster_size    = np.zeros([args.N])
for i in range(num_clusters):
    if (not mask[i]):
        continue
    n = np.argwhere(final_clusters == idx[i])
    cluster_spectra[n, :] = cluster_spectra[n, :] + spec[i,:]
    cluster_size[n]       = cluster_size[n] + 1

for i in range(args.N):
    cluster_spectra[i, :] = cluster_spectra[i, :] / cluster_size[i]

bg_cluster_i = np.argmax(cluster_size)
bg_cluster = cluster_spectra[bg_cluster_i, :]
bg_size    = cluster_size[bg_cluster_i]
cluster_spectra = np.delete(cluster_spectra, bg_cluster_i, 0)
cluster_size    = np.delete(cluster_size,    bg_cluster_i)
dl = np.percentile(bg_cluster, 10)
bg_cluster = np.maximum(bg_cluster, dl)
cluster_spectra = np.maximum(cluster_spectra, dl)


for i in range(cluster_spectra.shape[0]):
    cluster_spectra[i, :] = (cluster_spectra[i, :] - bg_cluster) / bg_cluster

import matplotlib.pyplot as plt

x = np.repeat(range(1,cluster_spectra.shape[1]+1), cluster_spectra.shape[0]).T.reshape(cluster_spectra.shape[1::-1]).T
x = x - np.repeat(np.arange(cluster_spectra.shape[0]) / cluster_spectra.shape[0] - 0.5, x.shape[1]).reshape(x.shape)

# plt.plot(x.T, cluster_spectra.T, linewidth=0.5)
# plt.legend(['{:d} ({:.0f})'.format(x, cluster_size[x]) for x in range(args.N)])
# plt.show()

figsize = (20,5)
h = plt.figure(figsize=figsize, dpi=196, tight_layout=True)

w = 1.0 / (x.shape[0] + 1.0)
for i in range(x.shape[0]):
    ax = plt.bar(x[i,:], cluster_spectra[i,:], w)
    
title = 'Relative Diff from Bg Cluster ({:d})'.format(int(bg_size))
if (args.linkage is not None):
    title = title + ': {:s} Linkage'.format(linkages[args.linkage])

plt.title(title)
plt.legend(['{:d} ({:.0f})'.format(x, cluster_size[x]) for x in range(len(cluster_size))])

if (args.output_filename is not None):
    plt.savefig(args.output_filename)
else:
    plt.show()
