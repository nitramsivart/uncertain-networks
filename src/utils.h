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

// These GLOBAL variables are defined at program initialization, in main.cpp
extern unsigned long vertices, edges;
extern unsigned communities;
extern double edge_ratio, edge_sum, message_converged_diff, zero_thresh;
extern bool unweighted_degree;

void check(double num, string str);
bool check_zero(double num, string str, double omega, double m1, double m2);

double entropy(double);

double SBM_joint_marginal(double omega, double Q, double m1, double m2);

// Computes the number of vertices and edges in the network.
void FindStats(unsigned long&, unsigned long&, ifstream&);

double GetEdgeRatio(double* degrees);
double GetEdgeSum(double* degrees);

// Sets the network.
void GetTheNetworkEdges(string, int, Trio*, double*, double*);

#endif
