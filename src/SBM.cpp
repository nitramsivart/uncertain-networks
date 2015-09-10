#include "SBM.h"

using namespace std;

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


/* Returning the mean field approximation terms */
double MFR_SBM(double** degree_only, const int& vertex, const int& community)
{
  return degree_only[0][community];
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
