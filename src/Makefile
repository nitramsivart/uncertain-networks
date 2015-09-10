CC=g++
CFLAGS=-O3
BOOST_INCLUDE = -I /usr/include
BOOST_LIB = -L/usr/libs
LIBS = -lboost_program_options -lgsl -lblas
STD = -std=c++11
ALL = $(CC) $(BOOST_INCLUDE) $(BOOST_LIB) $(CFLAGS) $(STD)  

#g++ -I /usr/include -L/usr/libs -std=c++11 -o main main.cpp SBM.cpp DCSBM.cpp utils.cpp -lboost_program_options -lgsl -lblas

all:
	$(ALL) main.cpp SBM.cpp DCSBM.cpp utils.cpp -o bp $(LIBS)

bp: main.o SBM.o DCSBM.o
	$(ALL) main.o SBM.o DCSBM.o -o bp

main.o: main.cpp
	$(ALL) main.cpp

SBM.o: SBM.cpp
	$(ALL) SBM.cpp
	
DCSBM.o: DCSBM.cpp
	$(ALL) DCSBM.cpp

clean:
	rm *o bp
