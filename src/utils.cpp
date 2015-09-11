
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <list>
#include <limits>
#include "math.h"
#include "utils.h"


void check(double num, string str) {
  if(isnan(num)) {
    printf("found nan!!!!");
    cout << str << "\n";
    exit(0);
  }
  if(isinf(num)) {
    printf("found inf!!!!");
    cout << str << "\n";
    exit(0);
  }
  return;
  if(0 == num) {
    printf("found zero!!!!");
    cout << str << "\n";
    //exit(0);
  }
}
bool check_zero(double num, string str, double omega, double m1, double m2) {
  if(0 == num) {
    printf("found zero!!!!");
    //cout << str << " " << omega << " " << m1 << " " << m2 << "\n";
    //exit(0);
    return false;
  }
  return true;
}


double entropy(double x)
{
  if (x == 0)
  {
    return 0;
  }
  return x * log(x);
}

// calculates q_rs^ij
double SBM_joint_marginal(double omega, double Q, double m1, double m2) {
  return (omega * Q + edge_ratio * (1-omega) * (1-Q)) * m1 * m2;
}


// Computes the number of vertices and edges in the network.
void FindStats(unsigned long& vertices, unsigned long& edges, ifstream& file)
{
  char* buffer;
  long int entry1 = 0;
  long int entry2 = 0;
  string lineread;

  while(std::getline(file, lineread)) // Read line by line
  {
    buffer = new char [lineread.size()+1];
    strcpy(buffer, lineread.c_str());
    sscanf(buffer, "%ld %ld %*s", &entry1, &entry2); // "%*s" ignores the rest of the line.
    if (entry1 != entry2)
    {
      edges++;
    }
    if (entry1 > vertices)
    {
      vertices = entry1;
    }
    if (entry2 > vertices)
    {
      vertices = entry2;
    }
    delete[] buffer;
  }
  vertices++;
  file.clear();
  file.seekg(0, ios::beg);
  return;
}

// Returns the sum of all edge weights
// vertices MUST be initialized
double GetEdgeSum(double* degrees)
{
  double total = 0;
  for(unsigned long i; i < vertices; i++) {
    total += degrees[i];
  }
  return total / 2;
}

// Get edges / nonedges, corresponding to the constant C from the paper
// vertices MUST be initialized
double GetEdgeRatio(double* degrees)
{
  double edge_count = GetEdgeSum(degrees);
  double nonedge_count = vertices * (vertices - 1) / 2. - edge_count;
  return edge_count / nonedge_count;
}

// Sets the network.
void GetTheNetworkEdges(string fileName, int lines, Trio* EdgeList, double* degrees, double* missing_degrees)
{
  ifstream InputFile;
  string lineread;
  char *buffer;
  int count_edge = 0;
  long int entry1 = 0;
  long int entry2 = 0;
  float entry3 = 0;

  InputFile.open(fileName.c_str());
  if (!InputFile)
  {
    cout << "Error in opening file";
    cin.get();
    return;
  }

  // DON'T COUNT THE DEGREES FOR THE UNKNOWN EDGES!!! The estimation of the omega matrix uses the *OBSERVED* degree
  while(count_edge < edges && std::getline(InputFile, lineread)) // Read line by line
  {
    buffer = new char [lineread.size()+1];
    strcpy(buffer, lineread.c_str());
    // weighted case
    // put a ! here if we want this to be the mean field for the unweighted case
    // no ! means that we will be doing the mean field for the weighted case, by counting actual
    // 		weights even though we've told it not too
    if (lines > 2 && !unweighted_degree)
    {
      sscanf(buffer, "%ld %ld %f %*s", &entry1, &entry2, &entry3);
      if (entry1 != entry2)
      {
        EdgeList[count_edge].out = entry1;
        EdgeList[count_edge].in = entry2;
        EdgeList[count_edge].known = entry3;

        if (EdgeList[count_edge].known)
        {
          degrees[entry1] += entry3;
          degrees[entry2] += entry3;
        }
        else
        {
          missing_degrees[entry1]++;
          missing_degrees[entry2]++;
        }
      }
      else
      {
        count_edge--;
      }
    }
    else
    {
      sscanf(buffer, "%ld %ld %*s", &entry1, &entry2);
      if (entry1 != entry2)
      {
        EdgeList[count_edge].out = entry1;
        EdgeList[count_edge].in = entry2;
        EdgeList[count_edge].known = 1.;

        degrees[entry1]++;
        degrees[entry2]++;
      }
      else
      {
        count_edge--;
      }
    }

    delete[] buffer;
    count_edge++;
  }
  InputFile.close();

  return;
}
