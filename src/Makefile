CC=g++
CFLAGS=-O3
BOOST_INCLUDE = -I /usr/include
BOOST_LIB = -L/usr/libs
LIBS = -lboost_program_options -lgsl -lblas
STD = -std=c++11
ALL = $(CC) $(BOOST_INCLUDE) $(BOOST_LIB) $(CFLAGS) $(STD)  

all:
	$(ALL) main.cpp SBM.cpp DCSBM.cpp utils.cpp -o bp $(LIBS)

clean:
	rm *o bp
