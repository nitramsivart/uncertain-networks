from fast_sbm import fast_sbm
from networkx import Graph, non_edges
from random import random, shuffle
import numpy as np

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

