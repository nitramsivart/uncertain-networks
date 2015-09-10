
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

// It's strange-looking, but this DOES give the right answer for x = 0, 1.
// Returns log(x!) without overflow-style problems
double log_factorial(int x)
{
  double result = 0;
  for (int i = 1; i < x; i++)
  {
    result += log(double(i));
  }
  return result;
}

double entropy(double x)
{
  if (x == 0)
  {
    return 0;
  }
  return x * log(x);
}

// Computes the number of vertices and edges in the network.
void FindStats(int& vertices, long int& edges, ifstream& file)
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

// Sets the network.
void GetTheNetworkEdges(string fileName, int lines, Trio* EdgeList, double* degrees, double* missing_degrees, const int edges, const bool unweighted_degree)
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
