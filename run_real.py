import sys
import os

trials = 300
comms = int(sys.argv[2]) #3
id_num = int(sys.argv[1]) #314723
unweighted = False#('True' == sys.argv[2])
home = '../../'
KL = home + 'projects/unobserved_BP/exact/'
#command = KL + 'full_KL_DC'
#command = KL + 'full_KL'
#command = '/home/travisbm/.clion10/system/cmake/generated/a0078319/a0078319/Debug/exact_no_error'
command = home + 'projects/unobserved_BP/KullbackLeibler/full_KL_DC'
command = KL + 'full_KL_DC'

data = home + 'data/string.embdl.exact/'
C_global = 0.
m_global = 0.
nodes = set([])
if unweighted:
  cutoff = float(sys.argv[3])
  big_component = False#('True' == sys.argv[4])
  print("true")
  command += '_unweighted'
  edge_param = 2
  edges = data + 'edges2/%d-%f-%r.edges' % (id_num, cutoff, big_component)
  #edges = data + 'edges/%d.edges' % id_num
  out = data + 'out/%d_%d_%f_%r.txt' % (id_num, comms, cutoff, big_component)
else:
  print("false")
  edge_param = 3
  edges = data + 'edges2/%d.edges' % id_num 
  out = data + 'out2/%d_%d_new.txt' % (id_num, comms)
  #edges = home + 'data/compare/edgest/wnet-1000-40.000000-20.000000-2.000000-1.000000-5.000000-0.edges'

  # calculating C_global, m_global
  for l in open(edges):
    m_global += float(l.split()[2])
    nodes.add(int(l.split()[0]))
    nodes.add(int(l.split()[1]))

  n = len(nodes)
  #m_global = 195000
  C_global = 1 * (m_global / ( n * (n-1)/2. - m_global))


print('%s %s %d %d %d %s %f %f' % (command, edges, comms, edge_param, trials, out, C_global, m_global))
os.system('%s %s %d %d %d %s %f %f' % (command, edges, comms, edge_param, trials, out, C_global, m_global)) #this is for degree corrected
#os.system('%s %s %d %d %d %s' % (command, edges, comms, edge_param, trials, out)) # this is for unweighted
