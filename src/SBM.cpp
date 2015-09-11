#include "SBM.h"


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
  return known * omega + edge_ratio * (1. - known) *(1. - omega);
  //OLD: return exp(known * log(omega) + (1. - known) * log(1. - omega));
}

/* M step initialization equations (can't assume the replica symmetric cavity equations since we don't have the messages) */
void M_SBM_init(Trio* EdgeList, double*** current_message, double** group_membership, double** omega, double* degrees)
{
  double* nr = new double[communities];
  double** denom = new double*[communities];
  for (unsigned i = 0; i < communities; i++)
  {
    nr[i] = 0;
    denom[i] = new double[communities];
    for (unsigned j = 0; j < communities; j++)
    {
      denom[i][j] = 0;
      omega[i][j] = 0;
    }
  }

  for (unsigned long i = 0; i < vertices; i++)
  {
    for (unsigned j = 0; j < communities; j++)
    {
      nr[j] += group_membership[i][j];
    }
  }

  for (unsigned long i = 0; i < edges; i++)
  {
    // Count towards m_rs
    if (EdgeList[i].known)
    {
      for (unsigned j = 0; j < communities; j++)
      {
        for (unsigned k = 0; k < communities; k++)
        {
          // Edges are undirected. I'm only buffering each edge once, so I need to count it twice, once for each direction.
          omega[j][k] += EdgeList[i].known * group_membership[EdgeList[i].out][j] * group_membership[EdgeList[i].in][k];
          omega[j][k] += EdgeList[i].known * group_membership[EdgeList[i].in][j] * group_membership[EdgeList[i].out][k];
        }
      }
    }
  }

  for (unsigned j = 0; j < communities; j++)
  {
    for (unsigned k = 0; k < communities; k++)
    {
      denom[j][k] += nr[j] * nr[k];
      if (denom[j][k] != 0)
        omega[j][k] /= denom[j][k];
      else
        printf("denom is 0, %d, %d\n", j, k);
    }
  }

  delete[] nr;
  for (unsigned j = 0; j < communities; j++)
  {
    delete[] denom[j];
  }
  delete[] denom;
}


/* M step equations */
void M_SBM(Trio* EdgeList, double*** current_message, double** group_membership, double** omega, double* degrees)
{
  double* nr = new double[communities];
  double** denom = new double*[communities];
  double** omega_old = new double*[communities]; // this stores omega before M step
  double** omega_temp = new double*[communities]; // this stores omega within M step, from previous loop
  double denom_color;
  for (unsigned i = 0; i < communities; i++)
  {
    omega_old[i] = new double[communities];
    omega_temp[i] = new double[communities];
    nr[i] = 0;
    denom[i] = new double[communities];
    for (unsigned j = 0; j < communities; j++)
    {
      omega_old[i][j] = omega[i][j];
    }
  }

  for (unsigned long i = 0; i < vertices; i++)
  {
    for (unsigned j = 0; j < communities; j++)
    {
      nr[j] += group_membership[i][j];
    }
  }

  for (int counter = 0; counter < 3; counter++) {
    for (unsigned i = 0; i < communities; i++)
    {
      for (unsigned j = 0; j < communities; j++)
      {
        denom[i][j] = 0;
        omega_temp[i][j] = omega[i][j];
        omega[i][j] = 0;
      }
    }
    //we set omega by another round of EM.
    for (unsigned long i = 0; i < edges; i++) {
      // non-edges shouldn't contribute to omega, and they give divide-by-zero errors
      if(EdgeList[i].known == 0) {
        continue;
      }
      denom_color = 0;
      for (unsigned j = 0; j < communities; j++)
      {
        for (unsigned k = 0; k < communities; k++)
        {
          denom_color += SBM_joint_marginal(omega_old[j][k], EdgeList[i].known, current_message[0][i][j], current_message[1][i][k]);
        }
      }
      for (unsigned j = 0; j < communities; j++)
      {
        for (unsigned k = 0; k < communities; k++)
        {
          double q0 = SBM_joint_marginal(omega_old[j][k], EdgeList[i].known, current_message[0][i][j], current_message[1][i][k]) / denom_color;
          double q1 = SBM_joint_marginal(omega_old[k][j], EdgeList[i].known, current_message[0][i][k], current_message[1][i][j]) / denom_color;


          // these are the joint marginals for edge / noneddge
          double qe = omega_temp[j][k] * EdgeList[i].known; 
          double qne = edge_ratio * (1-omega_temp[j][k]) * (1-EdgeList[i].known);
          double norm = qe + qne;
          qe /= norm;
          qne /= norm;


          // Edges are undirected. I'm only buffering each edge once, so I need to count it twice, once for each direction.
          omega[j][k] += q0 * qe;
          omega[j][k] += q1 * qe;
          denom[j][k] += q0;
          denom[j][k] += q1;
        }
      }

    }

    for (unsigned j = 0; j < communities; j++)
    {
      for (unsigned k = 0; k < communities; k++)
      {
        denom[j][k] = nr[j] * nr[k];
        if (denom[j][k] != 0){
          omega[j][k] /= denom[j][k];
        }
        else {
          //printf("setting omega to 0\n");
          omega[j][k] = 0.;
        }
      }
    }
  }


  delete[] nr;
  for (unsigned j = 0; j < communities; j++)
  {
    delete[] denom[j];
    delete[] omega_old[j];
    delete[] omega_temp[j];
  }
  delete[] denom;
  delete[] omega_old;
  delete[] omega_temp;
}




/* Computing the mean field approximation terms */
void MF_SBM(double** degree_only, double** group_membership, double* degrees, double* missing_degrees, bool* degrees_present, double** omega, double (*model)(const int, const double, double&, double&, double&))
{
  double part;
  double zero = 0.0;
  for (unsigned j = 0; j < communities; j++)
  {
    degree_only[0][j] = 0;
  }
  for (unsigned k = 0; k < communities; k++)
  {
    for (unsigned long i = 0; i < vertices; i++)
    {
      part = 0;
      for (unsigned j = 0; j < communities; j++)
      {
        // Not degree-specific, so just throw it whatever for the degree.
        part += group_membership[i][j] * model(0, 0., omega[j][k], zero, degrees[i]);
      }
      if (part == 0)
      {
        //printf("ERROR1");
        part = 1;
      }

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
double LL_SBM(Trio* EdgeList, double** group_membership, double* nKcount, double** omega, double* degrees)
{
  double likelihood = 0;
  // Should match nKcount since the function has converged, but it costs me nothing to calculate just in case.
  double* nx = new double[communities];
  for (unsigned j = 0; j < communities; j++)
  {
    nx[j] = 0;
    for (unsigned long i = 0; i < vertices; i++)
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
  for (unsigned long i = 0; i < edges; i++)
  {
    for (unsigned j = 0; j < communities; j++)
    {
      for (unsigned k = 0; k < communities; k++)
      {
        if (!(omega[j][k] == 0 || omega[j][k] == 1))
        {
          likelihood += group_membership[EdgeList[i].out][j] * group_membership[EdgeList[i].in][k] * log(EdgeList[i].known * omega[j][k] + edge_ratio * (1. - EdgeList[i].known ) * (1. - omega[j][k]));
        }
      }
    }
  }

  delete[] nx;
  return likelihood;
}
