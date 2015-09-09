// Computes the EM algorithm (using belief-propagation for the E step) for Kullback-Leibler divergence derived stochastic blockmodels.
// The base programming allows you to program in a mean-field approximation if you know how to. Unfortunately, the degree corrected stochastic blockmodel doesn't have an easy one that we can think of.
//  The standard stochastic blockmodel is programmed correctly, however.

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
#include <gsl/gsl_rng.h>  // REPLACERNG
#include "gsl/gsl_sf_gamma.h"
#include <boost/program_options.hpp>

using namespace std;

namespace po = boost::program_options;

struct Trio
{
  int out;
  int in;
  double known;
};


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


/* The M step has to reflect the function pointers as well! */
double SBM(const int adj_value, const double known, double& omega, double& k1, double& k2)
{
  if (adj_value == -1)
  {
    return 1;
  }
  // the below expression isn't defined for omega == 0, so we short circuit
  if (omega == 0) {
    //printf("short circuiting SBM\n");
    return 0;
  }
  return exp(known * log(omega) + (1. - known) * log(1. - omega));
}

// This is the Poisson formulation, which is why it looks different than the one above.
double DCSBM(const int adj_value, const double known, double& omega, double& k1, double& k2)
{
  if (adj_value == -1)
  {
    return 1;
  }
  // the below expression isn't defined for omega == 0, so we short circuit
  if (omega == 0)
    return 0;
  return exp(known * log(k1 * omega * k2) - k1 * omega * k2);
}




/* M step initialization equations (can't assume the replica symmetric cavity equations since we don't have the messages) */
void M_SBM_init(Trio* EdgeList, double*** current_message, double** group_membership, double** omega, double* degrees, const int& vertices, const int& edges, const int& communities)
{
  double* nr = new double[communities];
  double** denom = new double*[communities];
  for (int i = 0; i < communities; i++)
  {
    nr[i] = 0;
    denom[i] = new double[communities];
    for (int j = 0; j < communities; j++)
    {
      denom[i][j] = 0;
      omega[i][j] = 0;
    }
  }

  for (int i = 0; i < vertices; i++)
  {
    for (int j = 0; j < communities; j++)
    {
      nr[j] += group_membership[i][j];
    }
  }

  for (int i = 0; i < edges; i++)
  {
    // Count towards m_rs
    if (EdgeList[i].known)
    {
      for (int j = 0; j < communities; j++)
      {
        for (int k = 0; k < communities; k++)
        {
          // Edges are undirected. I'm only buffering each edge once, so I need to count it twice, once for each direction.
          omega[j][k] += EdgeList[i].known * group_membership[EdgeList[i].out][j] * group_membership[EdgeList[i].in][k];
          omega[j][k] += EdgeList[i].known * group_membership[EdgeList[i].in][j] * group_membership[EdgeList[i].out][k];
        }
      }
    }
  }

  for (int j = 0; j < communities; j++)
  {
    for (int k = 0; k < communities; k++)
    {
      denom[j][k] += nr[j] * nr[k];
      if (denom[j][k] != 0)
        omega[j][k] /= denom[j][k];
      else
        printf("denom is 0, %d, %d\n", j, k);
    }
  }

  delete[] nr;
  for (int j = 0; j < communities; j++)
  {
    delete[] denom[j];
  }
  delete[] denom;
}


void M_DCSBM_init(Trio* EdgeList, double*** current_message, double** group_membership, double** omega, double* degrees, const int& vertices, const int& edges, const int& communities)
{
  double* kappar = new double[communities];
  double** denom = new double*[communities];
  for (int i = 0; i < communities; i++)
  {
    kappar[i] = 0;
    denom[i] = new double[communities];
    for (int j = 0; j < communities; j++)
    {
      denom[i][j] = 0;
      omega[i][j] = 0;
    }
  }

  for (int i = 0; i < vertices; i++)
  {
    for (int j = 0; j < communities; j++)
    {
      kappar[j] += degrees[i] * group_membership[i][j];
    }
  }

  for (int i = 0; i < edges; i++)
  {
    // Count towards m_rs
    if (EdgeList[i].known)
    {
      for (int j = 0; j < communities; j++)
      {
        for (int k = 0; k < communities; k++)
        {
          // Edges are undirected. I'm only buffering each edge once, so I need to count it twice, once for each direction.
          omega[j][k] += EdgeList[i].known * group_membership[EdgeList[i].out][j] * group_membership[EdgeList[i].in][k];
          omega[j][k] += EdgeList[i].known * group_membership[EdgeList[i].in][j] * group_membership[EdgeList[i].out][k];
        }
      }
    }
  }

  for (int j = 0; j < communities; j++)
  {
    for (int k = 0; k < communities; k++)
    {
      denom[j][k] += kappar[j] * kappar[k];
      if (denom[j][k] != 0)
        omega[j][k] /= denom[j][k];
    }
  }

  delete[] kappar;
  for (int j = 0; j < communities; j++)
  {
    delete[] denom[j];
  }
  delete[] denom;
}




/* M step equations */
void M_SBM(Trio* EdgeList, double*** current_message, double** group_membership, double** omega, double* degrees, const int& vertices, const int& edges, const int& communities)
{
  double* nr = new double[communities];
  double** denom = new double*[communities];
  double** omega_old = new double*[communities];
  double denom_color;
  for (int i = 0; i < communities; i++)
  {
    omega_old[i] = new double[communities];
    nr[i] = 0;
    denom[i] = new double[communities];
    for (int j = 0; j < communities; j++)
    {
      denom[i][j] = 0;
      omega_old[i][j] = omega[i][j];
      omega[i][j] = 0;
    }
  }

  for (int i = 0; i < vertices; i++)
  {
    for (int j = 0; j < communities; j++)
    {
      nr[j] += group_membership[i][j];
    }
  }

  for (int i = 0; i < edges; i++)
  {
    // Count towards m_rs
    if (EdgeList[i].known)
    {
      denom_color = 0;
      for (int j = 0; j < communities; j++)
      {
        for (int k = 0; k < communities; k++)
        {
          denom_color += omega_old[j][k] * current_message[0][i][j] * current_message[1][i][k];
        }
      }
      for (int j = 0; j < communities; j++)
      {
        for (int k = 0; k < communities; k++)
        {
          // Edges are undirected. I'm only buffering each edge once, so I need to count it twice, once for each direction.
          omega[j][k] += EdgeList[i].known * omega_old[j][k] * current_message[0][i][j] * current_message[1][i][k] / denom_color;
          omega[j][k] += EdgeList[i].known * omega_old[k][j] * current_message[0][i][k] * current_message[1][i][j] / denom_color;
        }
      }
    }
  }

  for (int j = 0; j < communities; j++)
  {
    for (int k = 0; k < communities; k++)
    {
      denom[j][k] += nr[j] * nr[k];
      if (denom[j][k] != 0){
        omega[j][k] /= denom[j][k];
      }
      else {
        printf("setting omega to 0\n");
        omega[j][k] = 0.;
      }
    }
  }

  delete[] nr;
  for (int j = 0; j < communities; j++)
  {
    delete[] denom[j];
    delete[] omega_old[j];
  }
  delete[] denom;
  delete[] omega_old;
}


void M_DCSBM(Trio* EdgeList, double*** current_message, double** group_membership, double** omega, double* degrees, const int& vertices, const int& edges, const int& communities)
{
  double* kappar = new double[communities];
  double** denom = new double*[communities];
  double** omega_old = new double*[communities];
  double denom_color;
  for (int i = 0; i < communities; i++)
  {
    omega_old[i] = new double[communities];
    kappar[i] = 0;
    denom[i] = new double[communities];
    for (int j = 0; j < communities; j++)
    {
      denom[i][j] = 0;
      omega_old[i][j] = omega[i][j];
      omega[i][j] = 0;
    }
  }

  for (int i = 0; i < vertices; i++)
  {
    for (int j = 0; j < communities; j++)
    {
      kappar[j] += degrees[i] * group_membership[i][j];
    }
  }

  for (int i = 0; i < edges; i++)
  {
    // Count towards m_rs
    if (EdgeList[i].known)
    {
      denom_color = 0;
      for (int j = 0; j < communities; j++)
      {
        for (int k = 0; k < communities; k++)
        {
          denom_color += degrees[EdgeList[i].out] * omega_old[j][k] * degrees[EdgeList[i].in] * exp(-degrees[EdgeList[i].out] * omega_old[j][k] * degrees[EdgeList[i].in]) * current_message[0][i][j] * current_message[1][i][k];
        }
      }


      for (int j = 0; j < communities; j++)
      {
        for (int k = 0; k < communities; k++)
        {
          // Edges are undirected. I'm only buffering each edge once, so I need to count it twice, once for each direction.
          omega[j][k] += EdgeList[i].known * degrees[EdgeList[i].out] * omega_old[j][k] * degrees[EdgeList[i].in] * exp(-degrees[EdgeList[i].out] * omega_old[j][k] * degrees[EdgeList[i].in]) * current_message[0][i][j] * current_message[1][i][k] / denom_color;
          omega[j][k] += EdgeList[i].known * degrees[EdgeList[i].out] * omega_old[k][j] * degrees[EdgeList[i].in] * exp(-degrees[EdgeList[i].out] * omega_old[k][j] * degrees[EdgeList[i].in]) * current_message[0][i][k] * current_message[1][i][j] / denom_color;
        }
      }
    }
  }

  for (int j = 0; j < communities; j++)
  {
    for (int k = 0; k < communities; k++)
    {
      denom[j][k] += kappar[j] * kappar[k];
      if (denom[j][k] != 0)
        omega[j][k] /= denom[j][k];
      else
        omega[j][k] = 0.;
    }
  }

  delete[] kappar;
  for (int j = 0; j < communities; j++)
  {
    delete[] denom[j];
    delete[] omega_old[j];
  }
  delete[] denom;
  delete[] omega_old;
}






/* Computing the mean field approximation terms */
void MF_SBM(double** degree_only, double** group_membership, double* degrees, double* missing_degrees, bool* degrees_present, double** omega, const int& vertices, const int& communities, double (*model)(const int, const double, double&, double&, double&))
{
  double part;
  double zero = 0.0;
  for (int j = 0; j < communities; j++)
  {
    degree_only[0][j] = 0;
  }
  for (int k = 0; k < communities; k++)
  {
    for (int i = 0; i < vertices; i++)
    {
      part = 0;
      for (int j = 0; j < communities; j++)
      {
        // Not degree-specific, so just throw it whatever for the degree.
        part += group_membership[i][j] * model(0, 0., omega[j][k], zero, degrees[i]);
      }
      if (part == 0)
        part = 1e-10;
      degree_only[0][k] += log(part);
    }
  }
  return;
}


// Not actually a useful mean field algorithm. We have to do the full calculation since weighted degrees aren't integers
void MF_DCSBM(double** degree_only, double** group_membership, double* degrees, double* missing_degrees, bool* degrees_present, double** omega, const int& vertices, const int& communities, double (*model)(const int, const double, double&, double&, double&))
{
  double part;
  for (int i = 0; i < vertices; i++)
  {
    for (int k = 0; k < communities; k++)
    {
      degree_only[i][k] = 0;
    }

    for (int k = 0; k < communities; k++)
    {
      for (int j = 0; j < vertices; j++)
      {
        part = 0;
        for (int l = 0; l < communities; l++)
        {
          part += group_membership[j][l] * model(0, 0., omega[k][l], degrees[i], degrees[j]);
        }
        if (part == 0)
          part = 1e-10;
        degree_only[i][k] += log(part);
      }
    }
  }
  return;
}


// If the network is unweighted, you can make use of the mean field approximation
void uMF_DCSBM(double** degree_only, double** group_membership, double* degrees, double* missing_degrees, bool* degrees_present, double** omega, const int& vertices, const int& communities, double (*model)(const int, const double, double&, double&, double&))
{
  double part;
  double thisdegree;
  for (int i = 0; i < vertices; i++)
  {
    thisdegree = i;
    if (degrees_present[i])
    {
      for (int k = 0; k < communities; k++)
      {
        degree_only[i][k] = 0;
      }

      for (int k = 0; k < communities; k++)
      {
        for (int j = 0; j < vertices; j++)
        {
          part = 0;
          for (int l = 0; l < communities; l++)
          {
            part += group_membership[j][l] * model(0, 0., omega[k][l], thisdegree, degrees[j]);
          }
          degree_only[i][k] += log(part);
        }
      }
    }
  }
  return;
}








/* Returning the mean field approximation terms */
double MFR_SBM(double** degree_only, const int& vertex, const int& community)
{
  return degree_only[0][community];
}

// Acts as the degree or the particular vertex, depending on whether the edges are weighted or not.
double MFR_DCSBM(double** degree_only, const int& vertex, const int& community)
{
  return degree_only[vertex][community];
}



// Computes the full log-likelihood
double LL_SBM(Trio* EdgeList, double** group_membership, double* nKcount, double** omega, double* degrees, const int& vertices, const int& edges, const int& communities)
{
  double likelihood = 0;
  // Should match nKcount since the function has converged, but it costs me nothing to calculate just in case.
  double* nx = new double[communities];
  for (int j = 0; j < communities; j++)
  {
    nx[j] = 0;
    for (int i = 0; i < vertices; i++)
    {
      nx[j] += group_membership[i][j];

      // Takes care of the entropy term.
      likelihood -= entropy(group_membership[i][j]);
    }

    // Takes care of the prior term: Q(Z)log(P(Z))
    if (nKcount[j])
    {
      likelihood += nx[j] * log(nKcount[j]);
    }
  }
  // Edges and the missing chunk from the previous term.
  for (int i = 0; i < edges; i++)
  {
    for (int j = 0; j < communities; j++)
    {
      for (int k = 0; k < communities; k++)
      {
        if (!(omega[j][k] == 0 || omega[j][k] == 1))
        {
          likelihood += group_membership[EdgeList[i].out][j] * group_membership[EdgeList[i].in][k] * (EdgeList[i].known * log(omega[j][k]) + (1. - EdgeList[i].known ) * log(1. - omega[j][k]));
        }
      }
    }
  }

  delete[] nx;
  return likelihood;
}


double LL_DCSBM(Trio* EdgeList, double** group_membership, double* nKcount, double** omega, double* degrees, const int& vertices, const int& edges, const int& communities)
{
  double likelihood = 0;
  double* kappa = new double[communities];
  // Should match nKcount since the function has converged, but it costs me nothing to calculate just in case.
  double* nx = new double[communities];
  for (int j = 0; j < communities; j++)
  {
    kappa[j] = 0;
    nx[j] = 0;
    for (int i = 0; i < vertices; i++)
    {
      kappa[j] += group_membership[i][j] * degrees[i];
      nx[j] += group_membership[i][j];

      // Takes care of the entropy term.
      likelihood -= entropy(group_membership[i][j]);
    }

    // Takes care of the prior term: Q(Z)log(P(Z))
    if (nKcount[j])
    {
      likelihood += nx[j] * log(nKcount[j]);
    }
  }

  // Term in the exponent of the model. Simplified to save time.
  for (int j = 0; j < communities; j++)
  {
    for (int k = 0; k < communities; k++)
    {
      likelihood -= kappa[j] * omega[j][k] * kappa[k];
    }
  }

  // Edges and the missing chunk from the previous term.
  for (int i = 0; i < edges; i++)
  {
    for (int j = 0; j < communities; j++)
    {
      for (int k = 0; k < communities; k++)
      {
        likelihood += group_membership[EdgeList[i].out][j] * group_membership[EdgeList[i].in][k] * (EdgeList[i].known * log(degrees[EdgeList[i].out] * degrees[EdgeList[i].in] * omega[j][k]) - degrees[EdgeList[i].out] * degrees[EdgeList[i].in] * omega[j][k]);
      }
    }
  }

  delete[] kappa;
  delete[] nx;
  return likelihood;
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


double Compute_Maxes(Trio* EdgeList, double*** current_message, double** group_membership, double* nKcount, double** omega, double* degrees, const int& vertices, const int& edges, const int& communities, void (*M_model)(Trio*, double***, double**, double**, double*, const int&, const int&, const int&))
{
  double** omega_past = new double*[communities];
  for (int j = 0; j < communities; j++)
  {
    omega_past[j] = new double[communities];
    for (int k = 0; k < communities; k++)
    {
      omega_past[j][k] = omega[j][k];
      if(std::isnan(omega_past[j][k]))
        printf("omega_past[%d][%d] is nan\n", j,k);
    }

    // should only use this for the simulated SBM
    nKcount[j] = .5;
    // the follwing is for real world stuff
    /*
    for (int i = 0; i < vertices; i++)
    {
      nKcount[j] += group_membership[i][j];
    }
    nKcount[j] /= double(vertices);
    */
  }

  M_model(EdgeList, current_message, group_membership, omega, degrees, vertices, edges, communities);

  double converged = 0;
  for (int j = 0; j < communities; j++)
  {
    for (int k = 0; k < communities; k++)
    {
      if (omega_past[j][k] >= 1e-10 && converged < fabs((omega[j][k] - omega_past[j][k]) / omega_past[j][k]))
      {
        converged = fabs((omega[j][k] - omega_past[j][k]) / omega_past[j][k]);
        if(std::isnan(converged) || std::isnan(omega[j][k]) || std::isnan(omega_past[j][k]))
          printf("found a nan\n");
      }
    }
  }
  // The following is only for the simulated SBM
  omega[0][0] = max(omega[0][0], omega[1][1]);
  omega[1][1] = max(omega[0][0], omega[1][1]);
  for (int j = 0; j < communities; j++)
  {
    delete[] omega_past[j];
  }
  delete[] omega_past;

  return converged;

}




// For the messages, the first index (0) is for "out" sending to "in". (1) is for "in" sending to "out" (remember this is NOT symmetric!)
void BP_algorithm(Trio* EdgeList, double** group_membership, double** general_message, double*** current_message, double*** former_message, double* nKcount, double** omega, double** degree_only, double* degrees, double* missing_degrees, bool* degrees_present, double (*model)(const int, const double, double&, double&, double&), void (*MF_Precompute)(double**, double**, double*, double*, bool*, double**, const int&, const int&, double(const int, const double, double&, double&, double&)), double (*MF_Return)(double**, const int&, const int&), const int& vertices, const int& edges, const int& communities, const double& message_converged_diff, const double& zero_thresh, const bool& unweighted_degree)
{
  double max_diff = vertices;
  double partial_message_out, partial_message_in;
  double partial_message_out_denom, partial_message_in_denom;
  double norm1, norm2;
  double* temp_message1 = new double[communities];
  double* temp_message2 = new double[communities];
  bool set_zero;

  for (int i = 0; i < vertices; i++)
  {
    for (int j = 0; j < communities; j++)
    {
      general_message[i][j] = 0;
    }
  }

  int iteration = 0;
  while (max_diff > message_converged_diff && iteration < 20)
  {
    iteration++;
    // Pre-precompute the mean-field approximated term.
    MF_Precompute(degree_only, group_membership, degrees, missing_degrees, degrees_present, omega, vertices, communities, model);


    // Takes care of the mean-field term for the entire iteration!
    for (int i = 0; i < vertices; i++)
    {
      for (int j = 0; j < communities; j++)
      {
        if (unweighted_degree)
        {
          general_message[i][j] = MF_Return(degree_only, ceil(degrees[i]), j);
        }
        else
        {
          general_message[i][j] = MF_Return(degree_only, i, j);
        }
      }
    }


    // precomupte the general messages.
    max_diff = 0;
    for (int i = 0; i < edges; i++)
    {
      if (EdgeList[i].known)
      {
        for (int j = 0; j < communities; j++)
        {
          partial_message_out = 0;
          partial_message_in = 0;
          partial_message_out_denom = 0;
          partial_message_in_denom = 0;
          for (int k = 0; k < communities; k++)
          {
            partial_message_out += former_message[1][i][k] * model(1, EdgeList[i].known, omega[j][k], degrees[EdgeList[i].out], degrees[EdgeList[i].in]);
            partial_message_in += former_message[0][i][k] * model(1, EdgeList[i].known, omega[j][k], degrees[EdgeList[i].in], degrees[EdgeList[i].out]);

            partial_message_out_denom += group_membership[EdgeList[i].in][k] * model(0, 0., omega[j][k], degrees[EdgeList[i].out], degrees[EdgeList[i].in]);
            partial_message_in_denom += group_membership[EdgeList[i].out][k] * model(0, 0., omega[j][k], degrees[EdgeList[i].in], degrees[EdgeList[i].out]);
          }
          // if the denominator is 0, then this means 'out' will never connect to any nodes if it's in community j
          // either because there are no nodes in k or omega[j][k] is 0. we set this to a default value,
          // need to make sure we set it to the same value in the mean field calculation
          if (partial_message_out_denom == 0)
            partial_message_out_denom = 1e-10;
          general_message[EdgeList[i].out][j] += log(partial_message_out) - log(partial_message_out_denom);
          if (partial_message_in_denom == 0)
            partial_message_in_denom = 1e-10;
          general_message[EdgeList[i].in][j] += log(partial_message_in) - log(partial_message_in_denom);
        }
      }
      else
      {
        for (int j = 0; j < communities; j++)
        {
          partial_message_out = 0;
          partial_message_in = 0;
          partial_message_out_denom = 0;
          partial_message_in_denom = 0;
          for (int k = 0; k < communities; k++)
          {
            partial_message_out_denom += group_membership[EdgeList[i].in][k] * model(0, 0., omega[j][k], degrees[EdgeList[i].out], degrees[EdgeList[i].in]);
            partial_message_in_denom += group_membership[EdgeList[i].out][k] * model(0, 0., omega[j][k], degrees[EdgeList[i].in], degrees[EdgeList[i].out]);
          }
          // is it possible for messages to be not known?
          general_message[EdgeList[i].out][j] -= log(partial_message_out_denom);
          general_message[EdgeList[i].in][j] -= log(partial_message_in_denom);
        }
      }
    }

    // Need to take off the self-edge contribution that we're not counting
    for (int i = 0; i < vertices; i++)
    {
      for (int j = 0; j < communities; j++)
      {
        partial_message_out_denom = 0;
        for (int k = 0; k < communities; k++)
        {
          partial_message_out_denom += group_membership[i][k] * model(0, 0., omega[j][k], degrees[i], degrees[i]);
        }
        if (partial_message_out_denom == 0)
          partial_message_out_denom = 1e-10;
        general_message[i][j] -= log(partial_message_out_denom);
      }
    }


    // Then compute the actual messages
    for (int i = 0; i < edges; i++)
    {
      if (EdgeList[i].known)
      {
        for (int j = 0; j < communities; j++)
        {
          partial_message_out = 0;
          partial_message_in = 0;
          for (int k = 0; k < communities; k++)
          {
            partial_message_out += former_message[1][i][k] * model(1, EdgeList[i].known, omega[j][k], degrees[EdgeList[i].out], degrees[EdgeList[i].in]);
            partial_message_in += former_message[0][i][k] * model(1, EdgeList[i].known, omega[j][k], degrees[EdgeList[i].in], degrees[EdgeList[i].out]);
          }
          current_message[0][i][j] = general_message[EdgeList[i].out][j] - log(partial_message_out);
          current_message[1][i][j] = general_message[EdgeList[i].in][j] - log(partial_message_in);
        }

        for (int j = 0; j < communities; j++)
        {
          norm1 = nKcount[j];
          norm2 = nKcount[j];
          for (int k = 0; k < communities; k++)
          {
            if (j != k)
            {
              norm1 += nKcount[k] * exp(current_message[0][i][k] - current_message[0][i][j]);
              norm2 += nKcount[k] * exp(current_message[1][i][k] - current_message[1][i][j]);
            }
          }
          temp_message1[j] = nKcount[j] / norm1;
          temp_message2[j] = nKcount[j] / norm2;
        }

        for (int j = 0; j < communities; j++)
        {
          current_message[0][i][j] = temp_message1[j];
          current_message[1][i][j] = temp_message2[j];

          if (fabs(current_message[0][i][j] - former_message[0][i][j]) > max_diff)
          {
            max_diff = fabs(current_message[0][i][j] - former_message[0][i][j]);
          }

          if (fabs(current_message[1][i][j] - former_message[1][i][j]) > max_diff)
          {
            max_diff = fabs(current_message[1][i][j] - former_message[1][i][j]);
          }
        }
      }
    }


    for (int i = 0; i < edges; i++)
    {
      for (int j = 0; j < communities; j++)
      {
        former_message[0][i][j] = current_message[0][i][j];
        former_message[1][i][j] = current_message[1][i][j];
      }
    }


    for (int i = 0; i < vertices; i++)
    {
      norm1 = 0;

      for (int j = 0; j < communities; j++)
      {
        norm1 = nKcount[j];
        for (int k = 0; k < communities; k++)
        {
          if (k != j)
          {
            norm1 += nKcount[k] * exp(general_message[i][k] - general_message[i][j]);
          }
        }
        temp_message1[j] = nKcount[j] / norm1;
      }

      set_zero = false;
      for (int j = 0; j < communities; j++)
      {
        group_membership[i][j] = temp_message1[j];
        if (group_membership[i][j] < zero_thresh)
        {
          group_membership[i][j] = 0;
          set_zero = true;
        }
      }

      if (set_zero == true)
      {
        norm2 = 0;
        for (int j = 0; j < communities; j++)
        {
          norm2 += group_membership[i][j];
        }

        for (int j = 0; j < communities; j++)
        {
          group_membership[i][j] /= norm2;
        }
      }
    }
  }


  delete[] temp_message1;
  delete[] temp_message2;

  return;
}

bool Sort_list_probability(Trio& a, Trio& b)
{
	if (a.known < b.known)
	{
		return false;
	}
	return true;
}

// current assumption: if original edge probability was 0, we aren't gonna predict an edge.
void predict_edges(Trio* EdgeList, double*** current_message, double** omega, double* degrees, const int& edges, const int& communities, double (*model)(const int, const double, double&, double&, double&), char * out_file_name) {
  list<Trio> MissingEdgeList;
  for (int i = 0; i < edges; i++) {
    Trio edge = EdgeList[i];
    Trio unknown_edge;
    unknown_edge.out = edge.out;
    unknown_edge.in = edge.in;
    unknown_edge.known = 0; //known is the probability of the edge

    // The first part corresponds to P(A_ij | Z_i = r, Z_j = S, Theta), the probability of an edge given the communities
    // only the second part of the term is normalized, it corresponds to P(Z_i, Z_j | A^obs, Theta), the probability of observing r & s as communities.
    double normalization_factor = 0;
    for (int s = 0; s < communities; s++)
    {
      for (int r = 0; r < communities; r++)
      {
        double a_ij = model(1, edge.known, omega[r][s], degrees[edge.out], degrees[edge.in]);
        //maybe a little inconsistent with paper notation for out/in, but shouldn't matter because we're summing over all r,s
        double zi_zj = current_message[0][i][r] * current_message[1][i][s] * a_ij;
        //printf("msgs %f\t%f\t%f\t%f\n", current_message[0][i][r], current_message[1][i][s], zi_zj, a_ij);
        normalization_factor += zi_zj;
        unknown_edge.known += omega[r][s] * zi_zj;
      }
    }

    unknown_edge.known /= normalization_factor;
    double p = unknown_edge.known;
    double q = edge.known;
    // we are interested in the probability that there actually is an edge between i,j, given that
    // we know the probability according to the model, and according to experiment
    // prob = p / (p + (1-q) * p /q)
    if (q == 0)
      unknown_edge.known = 0;
    else
      unknown_edge.known = p / (p + (1-q) * (1-p) / q);
    //printf("final: %f\t%f\n", unknown_edge.known, normalization_factor);
    MissingEdgeList.push_back(unknown_edge);
  }

  MissingEdgeList.sort(Sort_list_probability);
  ofstream out_file;

  out_file.open((string(out_file_name) + ".edges").c_str());

  while (!MissingEdgeList.empty())
  {
    out_file << (MissingEdgeList.front()).out << "  " << (MissingEdgeList.front()).in << "  " << (MissingEdgeList.front()).known << endl;
    MissingEdgeList.pop_front();
  }

  out_file.close();

}



int main(int argc, char* argv[])
{
  // Likelihood and M-step of the EM algorithm.
  // Only ever change this line to change the program's behavior. It consistently picks the model.
  // Right now this is only doing legacy stuff, but it's still needed to compile and I didn't want to change functionality too much.
  const bool degree_correction = false;
  // True if the edges are unweighted. This lets you use the mean field approximation for the DCSBM.
  const bool unweighted_degree = false;
  // True if we should run the edge prediction algorithm after running the normal community detection
  const bool predict_edges_flag = false;

  double (*model)(const int, const double, double&, double&, double&);
  void (*M_model_init)(Trio*, double***, double**, double**, double*, const int&, const int&, const int&);
  void (*M_model)(Trio*, double***, double**, double**, double*, const int&, const int&, const int&);
  void (*MF_Precompute)(double**, double**, double*, double*, bool*, double**, const int&, const int&, double(const int, const double, double&, double&, double&));
  double (*MF_Return)(double**, const int&, const int&);
  double (*Compute_Likelihood)(Trio*, double**, double*, double**, double*, const int&, const int&, const int&);

  // If you want the degree corrected version, add the changes here. I had been using the flag "degree_correction" above and an if-then statement here.
  if (degree_correction == false)
  {
    model = SBM;
    M_model_init = M_SBM_init;
    M_model = M_SBM;
    MF_Precompute = MF_SBM;
    MF_Return = MFR_SBM;
    Compute_Likelihood = LL_SBM;
  }
  else
  {
    model = DCSBM;
    M_model_init = M_DCSBM_init;
    M_model = M_DCSBM;
    if (unweighted_degree)
    {
      MF_Precompute = uMF_DCSBM;
    }
    else
    {
      MF_Precompute = MF_DCSBM;
    }
    MF_Return = MFR_DCSBM;
    Compute_Likelihood = LL_DCSBM;
  }




  int vertices = 0;
  long int edges = 0;
  const gsl_rng_type * T;
  clock_t cpu_time = time(NULL);
  long int seed = (long int)(cpu_time);
  gsl_rng_default_seed = seed;
  T = gsl_rng_default;
  gsl_rng* numgen = gsl_rng_alloc (T);

  int communities = atoi(argv[2]);
  ifstream InputFile;
  InputFile.open(argv[1]);
  if (!InputFile)
  {
    cout << "Error in opening file";
    cin.get();
    return 0;
  }
  FindStats(vertices, edges, InputFile);
  InputFile.close();
  if (argc > 8)
  {
    vertices = atoi(argv[8]);
  }

  Trio* EdgeList = new Trio[edges];
  // This is the OBSERVED degree. The actual expected degree could be something different, depending on the number of unobserved vertex pairs.
  double* degrees = new double[vertices];
  double* missing_degrees = new double[vertices];
  // No vertex can have degree higher than the number of vertices
  bool* degrees_present = new bool[vertices];

  for (int i = 0; i < vertices; i++)
  {
    degrees[i] = 0;
    degrees_present[i] = false;
    missing_degrees[i] = 0;
  }

  GetTheNetworkEdges(string(argv[1]), atoi(argv[3]), EdgeList, degrees, missing_degrees, edges, unweighted_degree);

  // This is only useful if the degrees are already integers (i.e. the network is unweighted).
  for (int i = 0; i < vertices; i++)
  {
    // The addition of this ceil term is only to handle weighted 'histogram style' mean field approximating
    // I've removed it for the non-mean-field
    //degrees_present[int(ceil(degrees[i]))] = true;
    degrees_present[int(degrees[i])] = true;
  }

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "show help")
    ("cache_file,c", po::value<string>(&cache_file)->default_value(""), "cached input file containing edge list")
    ("input_file,i", po::value<string>(&input_file)->default_value(""), "input file containing edge list")
    ("output_file,o", po::value<string>(&output_file)->default_value(""), "output file for cascade statistics")
    ("mean_offspring,m", po::value< std::vector<float> >(&mean_offspring)->multitoken(), "mean offspring used to set beta = mean_offspring / mean_degree")
    ("offspring_alpha,a", po::value< std::vector<float> >(&offspring_alpha)->multitoken(), "exponent for power law distribution for mean offspring")
    ("mean_offspring_uniform,u", po::value< std::vector<float> >(&mean_offspring_uniform)->multitoken(), "mean offspring used to set beta = mean_offspring / mean_degree")
    ("num_cascades,n", po::value<unsigned>(&num_cascades)->default_value(100000), "number of cascades to generate")
    ("size,s", po::value<unsigned>(&min_size)->default_value(1),"minimum cascade size for which stats will be written")
    ("min_mean_offspring", po::value<float>(&min_mean_offspring)->default_value(0.1),"minimum mean offspring")
    ("max_mean_offspring", po::value<float>(&max_mean_offspring)->default_value(1),"maximum mean offspring")
    ("cascades_per_seed", po::value<int>(&cascades_per_seed)->default_value(1),"how many cascades are run from the same cascade seed")
    ("EM_restarts", po::value<int>(&restarts)->default_value(10),"random restarts of the full EM algorithm");
    ("EM_max_iterations", po::value<int>(&max_iterations)->default_value(20),"cap on number of EM iterations");
    ("BP_convergence_thresh", po::value<double>(&message_converged_diff)->default_value(0.001), "BP is considered converged when no message changes by more than this")
    ("EM_convergence_thresh", po::value<double>(&message_converged_diff)->default_value(5e-3), "EM is considered converged when no parameter changes by more than this")
    ("zero_thresh", po::value<double>(&zero_thresh)->default_value(1e-50), "any number below this is considered to be zero")


  int restarts;
  int max_iterations;
  double message_converged_diff;

  // For overall EM convergence.
  double converged_diff;

  double zero_thresh;




  double** group_membership = new double*[vertices];
  double** best_groups = new double*[vertices];
  double** general_message = new double*[vertices];
  double** degree_only;


  if (degree_correction)
  {
    degree_only = new double*[vertices];
  }
  else
  {
    // Seems kind of pointless, but it makes the program more modular.
    degree_only = new double*[1];
  }
  for (int i = 0; i < vertices; i++)
  {
    group_membership[i] = new double[communities];
    best_groups[i] = new double[communities];
    general_message[i] = new double[communities];
    if (degree_correction)
    {
      degree_only[i] = new double[communities];
    }
  }
  if (!degree_correction)
  {
    degree_only[0] = new double[communities];
  }

  double* nKcount = new double[communities];
  double*** current_message = new double**[2];
  double*** former_message = new double**[2];
  double*** best_message = new double**[2];
  double norm1, norm2;

  current_message[0] = new double*[edges];
  current_message[1] = new double*[edges];
  former_message[0] = new double*[edges];
  former_message[1] = new double*[edges];
  best_message[0] = new double*[edges];
  best_message[1] = new double*[edges];
  for (int i = 0; i < edges; i++)
  {
    current_message[0][i] = new double[communities];
    current_message[1][i] = new double[communities];
    former_message[0][i] = new double[communities];
    former_message[1][i] = new double[communities];
    best_message[0][i] = new double[communities];
    best_message[1][i] = new double[communities];
  }

  double** omega = new double*[communities];
  double** best_omega = new double*[communities];
  for (int j = 0; j < communities; j++)
  {
    omega[j] = new double[communities];
    best_omega[j] = new double[communities];
  }


  double LL, best_LL;
  int EMiterations;

  for (int thisiteration = 0; thisiteration < restarts; thisiteration++)
  {	
    // Begin (randomized) initialization steps
    printf("iteration %d\n", thisiteration);
    for (int j = 0; j < communities; j++)
    {
      nKcount[j] = 0;
    }

    for (int i = 0; i < vertices; i++)
    {
      norm1 = 0;
      for (int j = 0; j < communities; j++)
      {
        group_membership[i][j] = gsl_rng_uniform(numgen);
        norm1 += group_membership[i][j];
      }
      for (int j = 0; j < communities; j++)
      {
        group_membership[i][j] /= norm1;
        nKcount[j] += group_membership[i][j];
      }
    }

    for (int j = 0; j < communities; j++)
    {
      nKcount[j] = .5;// /= double(vertices);
    }

    for (int j = 0; j < edges; j++)
    {
      norm1 = 0;
      norm2 = 0;
      for (int k = 0; k < communities; k++)
      {
        former_message[0][j][k] = gsl_rng_uniform(numgen);
        former_message[1][j][k] = gsl_rng_uniform(numgen);

        //current message needs to be reset each time
        //current_message[0][j][k] = 0;
        //current_message[1][j][k] = 0;

        norm1 += former_message[0][j][k];
        norm2 += former_message[1][j][k];
      }
      for (int k = 0; k < communities; k++)
      {
        former_message[0][j][k] /= norm1;
        former_message[1][j][k] /= norm2;
      }
    }

    // why is omega set to 0 each time?? this gives problems
    // in the first run of compute maxes, which here is essentially
    // just skipping things
    for (int j = 0; j < communities; j++)
    {
      for (int k = 0; k < communities; k++)
      {
        omega[j][k] = 0;
      }
    }

    // Start with the M step. Except this should be the non-BP version, whereas the full EM algorithm one should include the messages as expected.
    double converged = Compute_Maxes(EdgeList, current_message, group_membership, nKcount, omega, degrees, vertices, edges, communities, M_model_init);
    converged = converged_diff * 4;
    for (int i = 0; i < communities; i++)
    {
      for (int j = 0; j < communities; j++)
      {
        if (i != j)
        {
          omega[i][j] *= 0.5;
        }
      }
    }
    EMiterations = 0;
    if (edges != 0)
    {
      while (converged > converged_diff && EMiterations < max_iterations)
      {
        printf("converged %f, EMiteration %d\n", converged, EMiterations);
        BP_algorithm(EdgeList, group_membership, general_message, current_message, former_message, nKcount, omega, degree_only, degrees, missing_degrees, degrees_present, model, MF_Precompute, MF_Return, vertices, edges, communities, message_converged_diff, zero_thresh, unweighted_degree);
        converged = Compute_Maxes(EdgeList, current_message, group_membership, nKcount, omega, degrees, vertices, edges, communities, M_model);
        for(int i=0; i < communities; i++) {
          for(int j=0; j < communities; j++) {
            cout << omega[i][j] << "\t";
          }
          cout << "\n";
        }
        EMiterations++;
      }

      LL = Compute_Likelihood(EdgeList, group_membership, nKcount, omega, degrees, vertices, edges, communities);
    }
    else
    {
      for (int i = 0; i < vertices; i++)
      {
        for (int j = 0; j < communities; j++)
        {
          group_membership[i][j] = 1. / communities;
        }
      }

      LL = Compute_Likelihood(EdgeList, group_membership, nKcount, omega, degrees, vertices, edges, communities);
    }

    if (std::isnan(LL))
    {
      LL = -std::numeric_limits<float>::max();
    }

    if (thisiteration == 0)
    {
      best_LL = LL;
      for (int j = 0; j < communities; j++)
      {
        for (int i = 0; i < vertices; i++)
        {
          best_groups[i][j] = group_membership[i][j];
        }

        for (int k = 0; k < communities; k++)
        {
          best_omega[j][k] = omega[j][k];
        }
        for (int i = 0; i < edges; i++)
        {
          best_message[0][i][j] = current_message[0][i][j];
          best_message[1][i][j] = current_message[1][i][j];
        }
      }
    }
    if (thisiteration == 0 || LL > best_LL)
    {
      best_LL = LL;
      for (int j = 0; j < communities; j++)
      {
        for (int i = 0; i < vertices; i++)
        {
          best_groups[i][j] = group_membership[i][j];
        }

        for (int k = 0; k < communities; k++)
        {
          best_omega[j][k] = omega[j][k];
        }
        for (int i = 0; i < edges; i++)
        {
          best_message[0][i][j] = current_message[0][i][j];
          best_message[1][i][j] = current_message[1][i][j];
        }
      }
      ofstream out_file;
      if (argc > 5)
      {
        out_file.open(argv[5]);
      }
      if (out_file.is_open())
      {
        for (int i = 0; i < vertices; i++)
        {
          for (int j = 0; j < communities; j++)
          {
            if (best_groups[i][j] < zero_thresh)
            {
              out_file << "0  ";
            }
            else
            {
              out_file << best_groups[i][j] << "  ";
            }
          }
          out_file << endl;
        }
      }
      out_file.close();
    }
    printf("LL: %f, Best LL: %f\n", LL, best_LL);
  }

  if(predict_edges_flag == true)
    predict_edges(EdgeList, best_message, omega, degrees, edges, communities, model, argv[5]);


  for (int i = 0; i < communities; i++)
  {
    delete[] omega[i];
    delete[] best_omega[i];
  }
  delete[] omega;
  delete[] best_omega;
  for (int i = 0; i < vertices; i++)
  {
    delete[] group_membership[i];
    delete[] best_groups[i];
    delete[] general_message[i];
  }
  for (int i = 0; i < edges; i++)
  {
    delete[] current_message[0][i];
    delete[] current_message[1][i];
    delete[] former_message[0][i];
    delete[] former_message[1][i];
    delete[] best_message[0][i];
    delete[] best_message[1][i];
  }
  delete[] group_membership;
  delete[] best_groups;
  delete[] general_message;
  delete[] current_message;
  delete[] former_message;
  delete[] best_message;
  delete[] degrees;
  delete[] degrees_present;
  delete[] missing_degrees;
  if (degree_correction)
  {
    for (int i = 0; i < vertices; i++)
    {
      delete[] degree_only[i];
    }
  }
  else
  {
    delete[] degree_only[0];
  }
  delete[] degree_only;
  delete[] EdgeList;
  return 0;
}
