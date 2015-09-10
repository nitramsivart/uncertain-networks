#ifndef INC_UTILS_H
#define INC_UTILS_H

#include <execinfo.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <list>
#include <limits>
#include "math.h"
#include <fstream>

using namespace std;

struct Trio
{
  int out;
  int in;
  double known;
};

// It's strange-looking, but this DOES give the right answer for x = 0, 1.
// Returns log(x!) without overflow-style problems
double log_factorial(int);

double entropy(double);

// Computes the number of vertices and edges in the network.
void FindStats(int&, long int&, ifstream&);

// Sets the network.
void GetTheNetworkEdges(string, int, Trio*, double*, double*, const int, const bool);

#endif
