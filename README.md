# uncertain-networks

Cpp files for running belief propagation can be found in src. Currently the only accepts edge lists where node labels are continuous integers starting at 0. For an example, see 'data/protein.edges'.


Must have installed:
- git
- A C++ compiler (default is g++).
- Boost program_options
- python3 (I use python3.4. python2 should work, but code might need tweaking)
- networkx (for synthetic network generation)


Steps to run:
- Clone this repository to your computer: "git clone https://github.com/tbmbob/uncertain-networks.git"
- Compile bp: "make -C src"
- Run on real world datasets with "./run_real.py -i data/protein.edges -k 4"
- Try generating synthetic networks and then running bp with "./run_synthetic.py"


Results are stored in the ./out folder. The 7th column of out/results.py contains information on community recovery - what fraction of nodes have the correct community predicted by belief propagation.


Get help on various methods using:
- ./run_real.py -h
- ./run_synthetic.py -h
- src/bp -h

