#!/usr/bin/python3

from networkx import Graph, write_edgelist, write_weighted_edgelist, non_edges
import networkx as nx
import numpy as np
import os
import sys
from optparse import OptionParser
from multiprocessing import Process
from fast_sbm import fast_sbm
from time import time
from datetime import datetime
from random import random, shuffle

def make_net(c_in, c_out, n):
  return fast_sbm(c_in, c_out, n)

# make a network uncertain
def fuzz_network(G_orig, threshold, b, edge_frac=1., nonedge_mult=5.):
  G = G_orig.copy()
  n = len(G.nodes())
  H = Graph()
  H.add_nodes_from(range(n))
  pairs = n * (n-1) / 2
  actual_edges = len(G.edges())
  edges = int(edge_frac * actual_edges)
  nonedges = int(edges * nonedge_mult)

  a = b / nonedge_mult

  # though these distributions are normalized to one, by selecting the appropriate number of edges
  # and nonedges, we make these 'distributions' correct
  edge_probs = np.random.beta(a+1, b, edges)
  nonedge_probs = np.random.beta(a,b+1, nonedges)

  # picking the right number of edges from the appropriate list
  edge_list = G.edges()
  nonedge_list = list(non_edges(G))
  shuffle(edge_list)
  shuffle(nonedge_list)
  for i in range(len(edge_probs)):
    G[edge_list[i][0]][edge_list[i][1]]['weight'] = edge_probs[i]
    if edge_probs[i] > threshold:
      H.add_edge(edge_list[i][0],edge_list[i][1])
  for i in range(len(nonedge_probs)):
    G.add_edge(nonedge_list[i][0],nonedge_list[i][1],weight = nonedge_probs[i])
    if nonedge_probs[i] > threshold:
      H.add_edge(nonedge_list[i][0], nonedge_list[i][1])

  return G, H


# looks through filename, sees what fraction of the first half of nodes are in the same community
# and likewise for the second half of nodes. n should be odd
def evaluate(filename, n):
  correct = 0
  for i,l in enumerate(open(filename)):
    line_floats = [float(x) for x in l.split()]
    if (i < n/2) == (line_floats[0] < .5):
      correct += 1
  correct = float(correct)/n
  return max(correct, 1-correct)
            
'''
Define a range of numbers by the count and increment
'''
def custom_range(start, count, inc):
    return np.arange(start, start + count*inc, inc)

'''
trials: number of times to run
'''
def main(n_start, n_count=1, n_inc=1, c_in_start=10, c_in_count=1, c_in_inc=1, c_out_start=5, c_out_count=1, c_out_inc=1, comm_count = 2, DC=False, i=0):
    bp_uncertain = 'src/bp'

    edge_frac = 1.
    nonedge_mult = 5.
    b = 2

    trials = 2
    comm_count = 2

    os.makedirs('out', exist_ok=True)
    os.makedirs('data', exist_ok=True)

    for n in custom_range(n_start, n_count, n_inc):
        for c_in in custom_range(c_in_start, c_in_count, c_in_inc):
            for c_out in custom_range(c_out_start, c_out_count, c_out_inc):
                uncertain_net = 'data/wnet-%d-%f-%f-%f-%f-%f-%d.edges'%(n,c_in,c_out,b,edge_frac,nonedge_mult, i)
                uncertain_comms = 'out/uncertain_comms-%d-%f-%f-%f-%f-%f-%d.txt'%(n,c_in,c_out,b,edge_frac,nonedge_mult, i)
 
                G_orig = make_net(c_in, c_out, n)
                G, _ = fuzz_network(G_orig, 1, b, edge_frac, nonedge_mult)
                write_weighted_edgelist(G, uncertain_net)
 
                start1 = time()
                print(uncertain_net)
                os.system('%s -i %s -o %s -c %d -l %d -n %d' % (bp_uncertain, uncertain_net, uncertain_comms, comm_count, 3, trials))
                end1 = time()

                with open('out/results.txt', 'a+') as out_file:
                    out_file.write("%d %f %f\t%f %f %f\t %f %f\t %s %d\n" %(n,
                                    c_in, c_out,
                                    b,edge_frac,nonedge_mult,
                                    evaluate(uncertain_comms, n), end1-start1,
                                    str(datetime.now()), i))


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-n', type=int, dest = 'n', help='number of nodes in network', default=1000)
    parser.add_option('--c_in', type=float, dest = 'c_in', help='average within-community degree', default=10)
    parser.add_option('--c_out', type=float, dest = 'c_out', help='average between-community degree', default=5)
    parser.add_option('-k', type=int, dest = 'comm_count', help='number of communities', default=2)
    parser.add_option('--iters', type=int, help='number of instances to run, multithreaded', default = 1)

    (options, _) = parser.parse_args()
    n = options.n

    ps = []
    for i in range(options.iters):
        p = Process(target=main(n))#c_cc, c_cp, x, n, i))
        ps.append(p)
        p.start()
    for p in ps:
        p.join()
