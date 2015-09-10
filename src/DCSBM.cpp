#include "DCSBM.h"

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



// Acts as the degree or the particular vertex, depending on whether the edges are weighted or not.
double MFR_DCSBM(double** degree_only, const int& vertex, const int& community)
{
  return degree_only[vertex][community];
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
