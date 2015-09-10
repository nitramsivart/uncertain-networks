#ifndef INC_SBM_H
#define INC_SBM_H

#include "utils.h"

/* The M step has to reflect the function pointers as well! */
double SBM(const int adj_value, const double known, double& omega, double& k1, double& k2);


/* M step initialization equations (can't assume the replica symmetric cavity equations since we don't have the messages) */
void M_SBM_init(Trio* EdgeList, double*** current_message, double** group_membership, double** omega, double* degrees, const int& vertices, const int& edges, const int& communities);


/* M step equations */
void M_SBM(Trio* EdgeList, double*** current_message, double** group_membership, double** omega, double* degrees, const int& vertices, const int& edges, const int& communities);


/* Computing the mean field approximation terms */
void MF_SBM(double** degree_only, double** group_membership, double* degrees, double* missing_degrees, bool* degrees_present, double** omega, const int& vertices, const int& communities, double (*model)(const int, const double, double&, double&, double&));


/* Returning the mean field approximation terms */
double MFR_SBM(double** degree_only, const int& vertex, const int& community);


// Computes the full log-likelihood
double LL_SBM(Trio* EdgeList, double** group_membership, double* nKcount, double** omega, double* degrees, const int& vertices, const int& edges, const int& communities);

#endif
