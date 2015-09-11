#ifndef INC_DCSBM_H
#define INC_DCSBM_H

#include "utils.h"

// This is the Poisson formulation, which is why it looks different than the one above.
double DCSBM(const int adj_value, const double known, double& omega, double& k1, double& k2);

void M_DCSBM_init(Trio* EdgeList, double*** current_message, double** group_membership, double** omega, double* degrees);

void M_DCSBM(Trio* EdgeList, double*** current_message, double** group_membership, double** omega, double* degrees);


// Not actually a useful mean field algorithm. We have to do the full calculation since weighted degrees aren't integers
void MF_DCSBM(double** degree_only, double** group_membership, double* degrees, double* missing_degrees, bool* degrees_present, double** omega, double (*model)(const int, const double, double&, double&, double&));


// If the network is unweighted, you can make use of the mean field approximation
void uMF_DCSBM(double** degree_only, double** group_membership, double* degrees, double* missing_degrees, bool* degrees_present, double** omega, double (*model)(const int, const double, double&, double&, double&));



// Acts as the degree or the particular vertex, depending on whether the edges are weighted or not.
double MFR_DCSBM(double** degree_only, const int& vertex, const int& community);

double LL_DCSBM(Trio* EdgeList, double** group_membership, double* nKcount, double** omega, double* degrees);

#endif
