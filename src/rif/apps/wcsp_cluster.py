#!/usr/bin/env python

from rif.cluster.cookie_cutter import *
import pandas as pd
import numpy as np
import sys
from io import StringIO


def update_clusters(thresh, lines, clusters, centers):
    norig = len(clusters)
    xincr = pd.read_csv(StringIO(''.join(lines)), ' ', header=None,
                        dtype='i4', engine='c').as_matrix()
    if len(clusters) > 0:
        xincr = np.concatenate((centers, xincr))
        clusters = cookie_cutter_update_i4(
            xincr, thresh, clusters, 'ndiff')
    else:
        clusters = cookie_cutter_i4(xincr, thresh, 'ndiff')
    for i in clusters:
        if i >= norig:
            # print(i, len(clusters))
            sys.stdout.write(lines[i - norig])
            sys.stdout.flush()
    centers = xincr[clusters, :]
    clusters = list(range(len(clusters)))
    lines.clear()
    return lines, clusters, centers


def main():
    if len(sys.argv) > 1 and sys.argv[1].endswith('.wcsp'):
        for fname in sys.argv[1:]:
            print('reading...')
            if fname.endswith('.wcsp'):
                x = pd.read_csv(fname, ' ', dtype='i4', header=None)
                # print(x.head(1))
            else:
                print('dont know how to handle', fname)
                sys.exit(-1)
            for thresh in (10, 15, 20, 25, 30, 35, 40, 45, 50):
                print('clustering', thresh)
                clust = cookie_cutter_i4(x, thresh, 'ndiff')
                x = x.iloc[clust, :]
                x.to_csv(fname + "_clust_" + str(int(thresh)) + "_" +
                         str(clust.shape[0]) + ".wcsp", ' ', header=None,
                         index=False)
    else:
        try:
            thresh = int(sys.argv[1])
        except IndexError:
            thresh = 10
        lines = list()
        clusters = list()
        centers = None
        for line in sys.stdin:
            # icolon = line.find(':')
            # tag = line[:icolon].split()
            # if len(tag) != 2:
                # continue
            # tag = tag[1]
            # if not tag.startswith('solution'):
                # continue
            # lines.append(line[icolon + 3:])
            if line.startswith('New rotamers: '):
                lines.append(line[14:])
            if len(lines) >= 1000:
                lines, clusters, centers = update_clusters(
                    thresh, lines, clusters, centers)
        if lines:
            lines, clusters, centers = update_clusters(
                thresh, lines, clusters, centers)


if __name__ == '__main__':
    main()
