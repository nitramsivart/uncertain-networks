#!/usr/bin/python3

import sys
import os
from optparse import OptionParser
from multiprocessing import Process

#from utils import *

'''
trials: number of times to run
'''
def main(file_name, comm_count = 2, DC=False, i=0):
    bp_uncertain = 'src/bp'

    trials = 2
    EM_iterations = 20

    os.makedirs('out', exist_ok=True)
    os.makedirs('data', exist_ok=True)

    out_file = file_name.replace('data/', '')
    out_file = 'out/' + out_file.replace('.edges', '') + '.out'

    os.system('%s -i %s -o %s -c %d -l %d -n %d -d %r' % (bp_uncertain, file_name, out_file, comm_count, 3, trials, DC))


if __name__ == '__main__':
    # parse command line options
    parser = OptionParser()
    parser.add_option('-i', type=str, dest = 'file_name', help='input file name', default='data/protein.edges')
    parser.add_option('-k', type=int, dest = 'comm_count', help='number of communities', default=2)
    parser.add_option('-d', action='store_true', dest='degree_corrected', default=False)
    parser.add_option('--iters', type=int, help='number of instances to run, multithreaded', default = 1)

    (options, _) = parser.parse_args()

    ps = []
    for i in range(options.iters):
        p = Process(target=main, args = (options.file_name, options.comm_count, options.degree_corrected, i))#c_cc, c_cp, x, n, i))
        ps.append(p)
        p.start()
    for p in ps:
        p.join()
