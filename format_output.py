#!/user/bin/python3

# Takes in an edgelist and corresponding communities, and
# outputs a gml file to be read by, for example, Gephi

import networkx as nx
import numpy as np
import sys
from random import random
from optparse import OptionParser

comm_count = int(sys.argv[2])
id_num = int(sys.argv[1])

weighted = (sys.argv[3] == 'True')
use_old_ordering = True
only_one_component = False
cluster_comms = False

downsample_edges = True

cutoff = 0.0
gml_location_name = None
gml_location_name = 'out2/' + sys.argv[1] + '.gml'

edgelist = 'edges2/%d.edges'%id_num
out = 'out2/%d_%d.gml'%(id_num, comm_count)
#posterior = 'out2/%d_%d_new.txt_posterior.edges' % (id_num, comm_count)
posterior = edgelist
comms = 'tem place holder'
comms = 'out2/%d_%d_new.txt'%(id_num, comm_count)

if weighted == False:
  edgelist = 'edges2/%d-%f-%r.edges'%(id_num, cutoff, only_one_component)
  #out = 'out2/%d_%d_%f_%r_un.gml'%(id_num, comm_count, cutoff, only_one_component)
  # the following two don't do the cutoff, the former two do
  #edgelist = 'edges2/%d.edges'%(id_num)
  out = 'out2/%d_%d_%f_%r_un.gml'%(id_num, comm_count, cutoff, only_one_component)
  comms = 'out/%d_%d_%f_%r.txt'%(id_num, comm_count, cutoff, only_one_component)

out = open(out, 'w+')

G = nx.read_edgelist(edgelist, nodetype = int, data=(('weight',float),))

# convert a file of strings to a list of the community for each vertex
comm_list = []
for l in open(comms):
  split = [float(x) for x in l.split()]
  comm_list.append((np.argmax(split), max(split)))

def node_to_str(name, comm, size, loc):
  ret = 'node\n[\nid %d\ncomm %d\ngraphics\n[\nw %f\n'%(name, comm, size)
  if loc is not None:
    ret += 'x %f\ny %f\n' % (loc[0], loc[1])
  ret += ']\n]\n'
  return ret


# takes in dict mapping old names to current names.
# and a location for each old name, as a dict of tuples
# and a graph with edge weights
# and whether the graph is weighted
# and a list of community tuples: (index of most likely community, likelihood of most likely comm)
def print_unweighted(old_names, comm_list, G, locs=None, weighted=True, cluster_comms=False):
  out.write('graph\n[\ndirected 0\n')
  if downsample_edges:
    G_post = nx.read_edgelist(posterior, nodetype = int, data = (('weight', float),))
  for n in G.nodes():
    if locs is not None:
      if old_names is not None:
        loc = locs[old_names[n]]
      else:
        loc = locs[n]
    else:
      loc = None

    out.write(node_to_str(n, comm_list[n][0], comm_list[n][1]*.7, loc))
  for e in G.edges(data=True):
    if weighted:
      if downsample_edges:
        # make the edge 1 with probability proportional to its weight
        r = random()
        if r < G_post[e[0]][e[1]]['weight']:
          out.write('edge\n[\nsource %d\ntarget %d\nweight %f\n]\n'%(e[0],e[1],1))
        #else:
        #  print(r, G_post[e[0]][e[1]]['weight'])
      else:
        out.write('edge\n[\nsource %d\ntarget %d\nweight %f\n]\n'%(e[0],e[1],e[2]['weight']))
    else:
      out.write('edge\n[\nsource %d\ntarget %d\n]\n'%(e[0],e[1]))

  if cluster_comms:
    # add extra edges between nodes in the same community
    for i in range(len(comm_list)):
      for j in range(len(comm_list)):
        if i != j and comm_list[i][0] == comm_list[j][0]:
          out.write('edge\n[\nsource %d\ntarget %d\nweight %f\n]\n'%(i,j,.1))
  out.write(']')

old_names = None
if not weighted and use_old_ordering:
  old_names = {}
  for i, l in enumerate(open('edges2/old-ordering-%d-%f-%r.txt'%(id_num, cutoff, only_one_component))):
    old_names[i] = int(l)

def extract_locs(gml_file_name):
  if gml_file_name is None:
    return None
  curr_id = None
  curr_x = None
  locs = {}
  for l in open(gml_file_name):
    if 'label' in l:
      curr_id = int(float(l.split()[1].replace('"', '')))
    elif 'x' in l:
      curr_x = float(l.split()[1])/27
    elif 'y' in l:
      curr_y = float(l.split()[1])/27
      locs[curr_id] = (curr_x, curr_y)
  return locs

print_unweighted(old_names, comm_list, G, extract_locs(gml_location_name), weighted, cluster_comms)

