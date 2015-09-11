#!/usr/bin/python3

from networkx import write_edgelist, write_weighted_edgelist
import os
from utils import *
from optparse import OptionParser
from multiprocessing import Process
from time import time
from datetime import datetime

def main(n_start, n_count=1, n_inc=1, c_in_start=10, c_in_count=1, c_in_inc=1, c_out_start=5, c_out_count=1, c_out_inc=1, comm_count = 2, DC=False, i=0):
    bp_uncertain = 'src/bp'

    edge_frac = 1.
    nonedge_mult = 5.
    b = 2
    trials = 2

    os.makedirs('out', exist_ok=True)
    os.makedirs('data', exist_ok=True)

    for n in custom_range(n_start, n_count, n_inc):
        for c_in in custom_range(c_in_start, c_in_count, c_in_inc):
            for c_out in custom_range(c_out_start, c_out_count, c_out_inc):
                original_net = 'data/original_net-%d-%f-%f-%f-%f-%f-%d.edges'%(n,c_in,c_out,b,edge_frac,nonedge_mult, i)
                uncertain_net = 'data/noisy_net-%d-%f-%f-%f-%f-%f-%d.edges'%(n,c_in,c_out,b,edge_frac,nonedge_mult, i)
                uncertain_comms = 'out/uncertain_comms-%d-%f-%f-%f-%f-%f-%d.out'%(n,c_in,c_out,b,edge_frac,nonedge_mult, i)
 
                print("making and fuzzing network")
                G_orig = make_net(c_in, c_out, n)
                write_edgelist(G_orig, original_net)
                G, _ = fuzz_network(G_orig, 1, b, edge_frac, nonedge_mult)
                write_weighted_edgelist(G, uncertain_net)
 
                start1 = time()
                print("running belief propagation")
                os.system('%s -i %s -o %s -c %d -l %d -n %d' % (bp_uncertain, uncertain_net, uncertain_comms, comm_count, 3, trials))
                end1 = time()

                with open('out/results.txt', 'a+') as out_file:
                    out_file.write("%d %f %f\t%f %f %f\t %f %f\t %s %d\n" %(n,
                                    c_in, c_out,
                                    b,edge_frac,nonedge_mult,
                                    evaluate(uncertain_comms, n), end1-start1,
                                    str(datetime.now()), i))


if __name__ == '__main__':
    # parse command line options
    parser = OptionParser()
    parser.add_option('-n', type=int, dest = 'n', help='number of nodes in network', default=2000)
    parser.add_option('--c_in', type=float, dest = 'c_in', help='average within-community degree', default=40)
    parser.add_option('--c_out', type=float, dest = 'c_out', help='average between-community degree', default=10)
    parser.add_option('-k', type=int, dest = 'comm_count', help='number of communities', default=2)
    parser.add_option('--iters', type=int, help='number of instances to run, multithreaded', default = 1)

    (options, _) = parser.parse_args()

    ps = []
    # run multiple threads
    for i in range(options.iters):
        p = Process(target=main, args = (options.n, 1, 1, options.c_in, 1, 1, options.c_out, 1, 1, options.comm_count, False, i))
        ps.append(p)
        p.start()
    for p in ps:
        p.join()
