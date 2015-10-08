#!/usr/bin/python3

# Takes in an edgelist and corresponding communities, and
# outputs a gml file to be read by, for example, Gephi

## For example, when I create plots, I read the gml file into
## Gephi, color by community using the Partition->Nodes tab,
## layout using ForceAtlas2, and then tweaking parameters
## and manually dragging nodes around

import networkx as nx
import numpy as np
import sys
from random import random
from optparse import OptionParser

# gives gml strings for nodes
def node_to_str(name, comm, size, proportional_size=False):
    if proportional_size:
        size = size*.7
    else:
        size = 20

    ret = 'node\n[\nid %d\ncomm %d\ngraphics\n[\nw %f\n'%(name, comm, size)
    ret += ']\n]\n'
    return ret

def combine(edgelist, comms, out_file, proportional_size):
    G = nx.read_edgelist(edgelist, nodetype = int, data=(('weight',float),))

    # get the community for each vertex, taken as the max weight in each line
    comm_list = []
    for l in open(comms):
        split = list(map(float,l.split()))
        comm_list.append((np.argmax(split), max(split)))
        '''
        try:
            print(np.argmax(split), max(split))
        except Exception:
            print(l)
        '''

    #output the graph file
    out = open(out_file, 'w+')
    out.write('graph\n[\ndirected 0\n')

    # first write nodes and their communities
    for n in G.nodes():
        out.write(node_to_str(n, comm_list[n][0], comm_list[n][1], proportional_size))

    # then write the edges
    for e in G.edges(data=True):
        if len(e) == 3 and 'weight' in e[2]:
            out.write('edge\n[\nsource %d\ntarget %d\nweight %f\n]\n'%(e[0],e[1],e[2]['weight']))
        else:
            out.write('edge\n[\nsource %d\ntarget %d\n\n]\n'%(e[0],e[1]))

    out.write(']')

if __name__ == '__main__':
    # parse command line options
    parser = OptionParser()
    parser.add_option('-i', type=str, dest = 'edgelist', help='input file edgelist', default='data/protein.edges')
    parser.add_option('-c', type=str, dest = 'comms', help='input communities (output from bp algorithm)', default='out/protein.out')
    parser.add_option('-o', type=str, dest = 'out_file', help='output file, containing vertices and their communities in gml format', default='out/protein.gml')
    parser.add_option('-s', action='store_true', dest='proportional_size', help='flag to make node size proportional to its community membership', default=False)

    (options, _) = parser.parse_args()

    combine(options.edgelist, options.comms, options.out_file, options.proportional_size)
