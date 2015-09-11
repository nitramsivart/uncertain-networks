#include "DCSBM.h"


// This calculates the equivalent of omega, adjusted by degree
double DCSBM_P(double& omega, double& k1, double& k2) {
  double p = k1 * omega * k2 / (2.*edge_sum);
  if(p > 1) {
    //printf("%f %f\n",k1,k2);
    //printf("too large");
    return 1;
  }
  else
    return p;
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
  //cout << edge_ratio << "  " << edge_sum << endl;
  double p = DCSBM_P(omega, k1, k2);
  return known * p + edge_ratio * (1. - known) * (1. - p);
  //return exp(log(known * k1 * omega * k2 + edge_ratio * (1-known)) - k1 * omega * k2) ;
}



void M_DCSBM_init(Trio* EdgeList, double*** current_message, double** group_membership, double** omega, double* degrees)
{
  double* kappar = new double[communities];
  double** denom = new double*[communities];
  for (unsigned i = 0; i < communities; i++)
  {
    kappar[i] = 0;
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
      kappar[j] += degrees[i] * group_membership[i][j];
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
      denom[j][k] += kappar[j] * kappar[k];
      if (denom[j][k] != 0)
        omega[j][k] /= (denom[j][k] / (2. * edge_sum));
    }
  }

  delete[] kappar;
  for (unsigned j = 0; j < communities; j++)
  {
    delete[] denom[j];
  }
  delete[] denom;
}






void M_DCSBM(Trio* EdgeList, double*** current_message, double** group_membership, double** omega, double* degrees)
{
  double* kappar = new double[communities];
  double** denom = new double*[communities];
  double** omega_old = new double*[communities]; // stores omega before M step
  double** omega_temp = new double*[communities]; // stores omega within M step, from prev loop
  double denom_color;
  for (unsigned i = 0; i < communities; i++)
  {
    omega_old[i] = new double[communities];
    omega_temp[i] = new double[communities];
    kappar[i] = 0;
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
      kappar[j] += degrees[i] * group_membership[i][j];
    }
  }

  for (int counter = 0; counter < 1; counter++) {
    for (unsigned i = 0; i < communities; i++)
    {
      for (unsigned j = 0; j < communities; j++)
      {
        denom[i][j] = 0;
        omega_temp[i][j] = omega[i][j];
        omega[i][j] = 0;
      }
    }

    for (unsigned long i = 0; i < edges; i++)
    {
      // Count towards m_rs
      if (EdgeList[i].known != 0)
      {
        denom_color = 0;
        for (unsigned j = 0; j < communities; j++)
        {
          for (unsigned k = 0; k < communities; k++)
          {
            double p_old = DCSBM_P(omega_old[j][k], degrees[EdgeList[i].out], degrees[EdgeList[i].in]);
            denom_color += SBM_joint_marginal(p_old, EdgeList[i].known, current_message[0][i][j], current_message[1][i][k]);
          }
        }


        for (unsigned j = 0; j < communities; j++)
        {
          for (unsigned k = 0; k < communities; k++)
          {
            // this was omega for regular SBM. assume matrix is symmetric
            double p_old = DCSBM_P(omega_old[j][k], degrees[EdgeList[i].out], degrees[EdgeList[i].in]);

            // Edges are undirected. I'm only buffering each edge once, so I need to count it twice, once for each direction.
            double q0 = SBM_joint_marginal(p_old, EdgeList[i].known, current_message[0][i][j], current_message[1][i][k]) / denom_color;
            double q1 = SBM_joint_marginal(p_old, EdgeList[i].known, current_message[0][i][k], current_message[1][i][j]) / denom_color;


            double p_temp = DCSBM_P(omega_temp[j][k], degrees[EdgeList[i].out], degrees[EdgeList[i].in]);

            double qe = p_temp * EdgeList[i].known;
            double qne = edge_ratio * (1-p_temp) * (1 - EdgeList[i].known);
            double norm = qe + qne;
            qe /= norm;

            omega[j][k] += q0*qe;
            omega[j][k] += q1*qe;
            denom[j][k] += q0;
            denom[j][k] += q1;
          }
        }
      }
    }

    for (unsigned j = 0; j < communities; j++)
    {
      for (unsigned k = 0; k < communities; k++)
      {
        denom[j][k] = kappar[j] * kappar[k];
        omega[j][k] *= 2. * edge_sum;
        if (denom[j][k] != 0)
          omega[j][k] /= denom[j][k];
        else{
          //printf("setting omega to 0\n");
          omega[j][k] = 0.;
        }
      }
    }
  }

  delete[] kappar;
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






// Not actually a useful mean field algorithm. We have to do the full calculation since weighted degrees aren't integers
void MF_DCSBM(double** degree_only, double** group_membership, double* degrees, double* missing_degrees, bool* degrees_present, double** omega, double (*model)(const int, const double, double&, double&, double&))
{
  double part;
  for (unsigned long i = 0; i < vertices; i++)
  {
    for (unsigned k = 0; k < communities; k++)
    {
      degree_only[i][k] = 0;
    }

    for (unsigned k = 0; k < communities; k++)
    {
      for (unsigned long j = 0; j < vertices; j++)
      {
        part = 0;
        for (unsigned l = 0; l < communities; l++)
        {
          //check_zero(group_membership[j][l], "mf_group_membership",0,0,0);
          part += group_membership[j][l] * model(0, 0., omega[k][l], degrees[i], degrees[j]);
          //printf("gmemb, model: %f %f\n", group_membership[j][l], model(0,0,omega[k][l], degrees[i], degrees[j]));
        }
        if (part == 0) {
/*          if(group_membership[j][0] == 0 and group_membership[j][1] == 0 and group_membership[j][2] == 0)
            printf("G");*/
          //if(model(0,0,omega[k][0],degrees[i],degrees[j]) == 0) {
          //  printf("M %d %d %f %f", i, j, degrees[i], degrees[j]);
          //}
          //part = ERROR;
          part = 1;
          //printf("zero part\n");
        }
        degree_only[i][k] += log(part);
        //printf("i,part: %d %f\n", i, part);
        //check(degree_only[i][k], "degree_only");
      }
    }
  }
  return;
}

// If the network is unweighted, you can make use of the mean field approximation
void uMF_DCSBM(double** degree_only, double** group_membership, double* degrees, double* missing_degrees, bool* degrees_present, double** omega, double (*model)(const int, const double, double&, double&, double&))
{
  double part;
  double thisdegree;
  for (unsigned long i = 0; i < vertices; i++)
  {
    thisdegree = i;
    if (degrees_present[i])
    {
      for (unsigned k = 0; k < communities; k++)
      {
        degree_only[i][k] = 0;
      }

      for (unsigned k = 0; k < communities; k++)
      {
        for (unsigned long j = 0; j < vertices; j++)
        {
          part = 0;
          for (unsigned l = 0; l < communities; l++)
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



double LL_DCSBM(Trio* EdgeList, double** group_membership, double* nKcount, double** omega, double* degrees)
{
  double likelihood = 0;
  double* kappa = new double[communities];
  // Should match nKcount since the function has converged, but it costs me nothing to calculate just in case.
  double* nx = new double[communities];
  for (unsigned j = 0; j < communities; j++)
  {
    kappa[j] = 0;
    nx[j] = 0;
    for (unsigned long i = 0; i < vertices; i++)
    {
      kappa[j] += group_membership[i][j] * degrees[i];
      nx[j] += group_membership[i][j];

      // Takes care of the entropy term.
      likelihood -= entropy(group_membership[i][j]);
    }
    printf("likelihood after entropy: %f\n", likelihood);

    // Takes care of the prior term: Q(Z)log(P(Z))
    if (nKcount[j])
    {
      likelihood += nx[j] * log(nKcount[j]);
    }
    printf("likelihood after prior: %f\n", likelihood);
  }

  printf("likelihood before model exp: %f\n", likelihood);

  // Term in the exponent of the model. Simplified to save time.
  for (unsigned j = 0; j < communities; j++)
  {
    for (unsigned k = 0; k < communities; k++)
    {
      likelihood -= kappa[j] * omega[j][k] * kappa[k] / (2*edge_sum);
    }
  }
  printf("likelihood before edges: %f\n", likelihood);

  // Edges and the missing chunk from the previous term.
  int print_once_count = 0;
  for (unsigned long i = 0; i < edges; i++)
  {
    for (unsigned j = 0; j < communities; j++)
    {
      for (unsigned k = 0; k < communities; k++)
      {
        double p = DCSBM_P(omega[j][k], degrees[EdgeList[i].out], degrees[EdgeList[i].in]);
          likelihood += group_membership[EdgeList[i].out][j] * group_membership[EdgeList[i].in][k] *
                                                               log(EdgeList[i].known * p + edge_ratio * (1. - EdgeList[i].known) * (1. - p));
      }
    }
  }

  delete[] kappa;
  delete[] nx;
  return likelihood;
}

