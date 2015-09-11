// Computes the EM algorithm (using belief-propagation for the E step) for Kullback-Leibler divergence derived stochastic blockmodels.
// The base programming allows you to program in a mean-field approximation if you know how to. Unfortunately, the degree corrected stochastic blockmodel doesn't have an easy one that we can think of.
//  The standard stochastic blockmodel is programmed correctly, however.

#include <execinfo.h>
#include <stdio.h>
#include <string>
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <list>
#include <limits>
#include "math.h"
#include <boost/program_options.hpp>
#include "SBM.h"
#include "DCSBM.h"
#include <fstream>

using namespace std;

namespace po = boost::program_options;


double ERROR = 1e-10;

bool unweighted_degree; // THIS IS A GLOBAL VARIABLE
unsigned long vertices;
unsigned long edges;
unsigned communities;
double edge_ratio, edge_sum, message_converged_diff, zero_thresh;


double Compute_Maxes(bool not_init, Trio* EdgeList, double*** current_message, double** group_membership, double* nKcount, double** omega, double* degrees, void (*M_model)(Trio*, double***, double**, double**, double*))
{
  //// POTENTIAL FOR MEMORY LEAK
  double** omega_past = new double*[communities];
  double converged = 0;

  for (unsigned j = 0; j < communities; j++)
  {
    omega_past[j] = new double[communities];
    for (unsigned k = 0; k < communities; k++)
    {
      omega_past[j][k] = omega[j][k];
      if(std::isnan(omega_past[j][k]))
        printf("omega_past[%d][%d] is nan\n", j,k);
    }

    // should only use this for the simulated SBM
    //nKcount[j] = .5;
    nKcount[j] = 0;
    // the follwing is for real world stuff
    for (unsigned long i = 0; i < vertices; i++)
    {
      nKcount[j] += group_membership[i][j];
      //if(not_init)
        //printf("group memb: %f\n", group_membership[i][j]);
    }
    //printf("nKcount, vertices: %f, %d\n", nKcount[j], vertices);
    nKcount[j] /= double(vertices);
    //nKcount[j] = .5;
  }

  M_model(EdgeList, current_message, group_membership, omega, degrees);


  for (unsigned j = 0; j < communities; j++)
  {
    for (unsigned k = 0; k < communities; k++)
    {
      if (omega_past[j][k] >= 1e-10 && converged < fabs((omega[j][k] - omega_past[j][k]) / omega_past[j][k]))
      {
        converged = fabs((omega[j][k] - omega_past[j][k]) / omega_past[j][k]);
        if(std::isnan(converged) || std::isnan(omega[j][k]) || std::isnan(omega_past[j][k]))
          printf("found a nan\n");
      }
    }
  }

  //omega[0][0] = max(omega[0][0], omega[1][1]);
  //omega[1][1] = max(omega[0][0], omega[1][1]);

  //omega[0][0] = 0.04;
  //omega[1][1] = 0.04;
  //omega[0][1] = 0.028;
  //omega[1][0] = 0.028;

  for (unsigned j = 0; j < communities; j++)
  {
    delete[] omega_past[j];
  }
  delete[] omega_past;

  return converged;

}







// For the messages, the first index (0) is for "out" sending to "in". (1) is for "in" sending to "out" (remember this is NOT symmetric!)
void BP_algorithm(Trio* EdgeList, double** group_membership, double** general_message, double*** current_message, double*** former_message, double* nKcount, double** omega, double** degree_only, double* degrees, double* missing_degrees, bool* degrees_present, double (*model)(const int, const double, double&, double&, double&), void (*MF_Precompute)(double**, double**, double*, double*, bool*, double**, double(const int, const double, double&, double&, double&)), double (*MF_Return)(double**, const int&, const int&))
{
  double max_diff = vertices;
  double partial_message_out, partial_message_in;
  double partial_message_out_denom, partial_message_in_denom;
  double norm1, norm2;
  double* temp_message1 = new double[communities];
  double* temp_message2 = new double[communities];
  bool set_zero;

  for (unsigned long i = 0; i < vertices; i++)
  {
    for (unsigned j = 0; j < communities; j++)
    {
      general_message[i][j] = 0;
    }
  }

  int iteration = 0;
  while (max_diff > message_converged_diff && iteration < 20)
  {
    iteration++;
    // Pre-precompute the mean-field approximated term.
    MF_Precompute(degree_only, group_membership, degrees, missing_degrees, degrees_present, omega, model);


    // Takes care of the mean-field term for the entire iteration!
    for (unsigned long i = 0; i < vertices; i++)
    {
      for (unsigned j = 0; j < communities; j++)
      {
        if (unweighted_degree)
        {
          general_message[i][j] = MF_Return(degree_only, ceil(degrees[i]), j);
        }
        else
        {
          general_message[i][j] = MF_Return(degree_only, i, j);
          //printf("%d\n",i);
          //check(general_message[i][j], "gen_msg");
        }
      }
    }


    // precomupte the general messages.
    max_diff = 0;
    for (unsigned long i = 0; i < edges; i++)
    {
      if (EdgeList[i].known)
      {
        for (unsigned j = 0; j < communities; j++)
        {
          partial_message_out = 0;
          partial_message_in = 0;
          partial_message_out_denom = 0;
          partial_message_in_denom = 0;
          for (unsigned k = 0; k < communities; k++)
          {
            partial_message_out += former_message[1][i][k] * model(1, EdgeList[i].known, omega[j][k], degrees[EdgeList[i].out], degrees[EdgeList[i].in]);
            partial_message_in += former_message[0][i][k] * model(1, EdgeList[i].known, omega[j][k], degrees[EdgeList[i].in], degrees[EdgeList[i].out]);

            //check(partial_message_out, "partial_message_out");
            //check(partial_message_in, "partial_message_out");

            partial_message_out_denom += group_membership[EdgeList[i].in][k] * model(0, 0., omega[j][k], degrees[EdgeList[i].out], degrees[EdgeList[i].in]);
            partial_message_in_denom += group_membership[EdgeList[i].out][k] * model(0, 0., omega[j][k], degrees[EdgeList[i].in], degrees[EdgeList[i].out]);

            //check(partial_message_out, "partial_message_out_denom");
            //check(partial_message_in, "partial_message_out_denom");
          }
          // if the denominator is 0, then this means 'out' will never connect to any nodes if it's in community j
          // either because there are no nodes in k or omega[j][k] is 0. we set this to a default value,
          // need to make sure we set it to the same value in the mean field calculation
          if(partial_message_out == 0) {
            //printf("ERROR3.num");
            partial_message_out = 1;
          }
          if (partial_message_out_denom == 0)
          {
            //printf("ERROR3");
            partial_message_out_denom = 1;
          }
          general_message[EdgeList[i].out][j] += log(partial_message_out) - log(partial_message_out_denom);

          if (partial_message_in == 0) {
            //printf("ERROR4.num");
            partial_message_in = 1;
          }
          if (partial_message_in_denom == 0)
          {
            //printf("ERROR4");
            partial_message_in_denom = 1;
          }
          general_message[EdgeList[i].in][j] += log(partial_message_in) - log(partial_message_in_denom);
        }
      }
      else
      {
        for (unsigned j = 0; j < communities; j++)
        {
          partial_message_out = 0;
          partial_message_in = 0;
          partial_message_out_denom = 0;
          partial_message_in_denom = 0;
          for (unsigned k = 0; k < communities; k++)
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
    for (unsigned long i = 0; i < vertices; i++)
    {
      for (unsigned j = 0; j < communities; j++)
      {
        partial_message_out_denom = 0;
        for (unsigned k = 0; k < communities; k++)
        {
          partial_message_out_denom += group_membership[i][k] * model(0, 0., omega[j][k], degrees[i], degrees[i]);
        }
        if (partial_message_out_denom == 0)
        {//printf("ERROR5 %d %f\n",i, degrees[i]);
          partial_message_out_denom = 1;
          }

        general_message[i][j] -= log(partial_message_out_denom);
      }
    }


    // Then compute the actual messages
    for (unsigned long i = 0; i < edges; i++)
    {
      if (EdgeList[i].known)
      {
        for (unsigned j = 0; j < communities; j++)
        {
          partial_message_out = 0;
          partial_message_in = 0;
          for (unsigned k = 0; k < communities; k++)
          {
            partial_message_out += former_message[1][i][k] * model(1, EdgeList[i].known, omega[j][k], degrees[EdgeList[i].out], degrees[EdgeList[i].in]);
            partial_message_in += former_message[0][i][k] * model(1, EdgeList[i].known, omega[j][k], degrees[EdgeList[i].in], degrees[EdgeList[i].out]);
          }
          if(partial_message_out == 0)
            partial_message_out = 1;
          if(partial_message_in == 0)
            partial_message_in = 1;

          current_message[0][i][j] = general_message[EdgeList[i].out][j] - log(partial_message_out);
          current_message[1][i][j] = general_message[EdgeList[i].in][j] - log(partial_message_in);
        }

        for (unsigned j = 0; j < communities; j++)
        {
          norm1 = nKcount[j];
          norm2 = nKcount[j];
          for (unsigned k = 0; k < communities; k++)
          {
            if (j != k)
            {
              norm1 += nKcount[k] * exp(current_message[0][i][k] - current_message[0][i][j]);
              norm2 += nKcount[k] * exp(current_message[1][i][k] - current_message[1][i][j]);
            }
          }
          temp_message1[j] = nKcount[j] / norm1;
          temp_message2[j] = nKcount[j] / norm2;
          /*if(current_message[0][i][j] <= ERROR) {
            printf("no way");
            temp_message1[j] = .5;
          }
          else {
            temp_message1[j] = nKcount[j] / norm1;
          }

          if(current_message[0][i][j] <= ERROR) {
            printf("no way");
            temp_message2[j] = .5;
          }
          else {
            temp_message2[j] = nKcount[j] / norm2;
          }*/
        }

        for (unsigned j = 0; j < communities; j++)
        {
          current_message[0][i][j] = temp_message1[j];
          current_message[1][i][j] = temp_message2[j];

          //check(current_message[0][i][j], "current_message[0]");
          //check(current_message[1][i][j], "current_message[1]");

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


    for (unsigned long i = 0; i < edges; i++)
    {
      for (unsigned j = 0; j < communities; j++)
      {
        former_message[0][i][j] = current_message[0][i][j];
        former_message[1][i][j] = current_message[1][i][j];
      }
    }


    for (unsigned long i = 0; i < vertices; i++)
    {
      for (unsigned j = 0; j < communities; j++)
      {
        norm1 = nKcount[j];
        for (unsigned k = 0; k < communities; k++)
        {
          if (k != j)
          {
            // we do the subtraction this way to avoid exponentiating large values
            // should work out to the same normalization factor
            norm1 += nKcount[k] * exp(general_message[i][k] - general_message[i][j]);
          }
        }
        temp_message1[j] = nKcount[j] / norm1;
      }

      set_zero = false;
      for (unsigned j = 0; j < communities; j++)
      {
        group_membership[i][j] = temp_message1[j];
        if (group_membership[i][j] < zero_thresh)
        {
          group_membership[i][j] = 0;
          set_zero = true;
        }
        //check(group_membership[i][j], "group_membership");
      }

      if (set_zero == true)
      {
        norm2 = 0;
        for (unsigned j = 0; j < communities; j++)
        {
          norm2 += group_membership[i][j];
        }

        for (unsigned j = 0; j < communities; j++)
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
void predict_edges(Trio* EdgeList, double*** current_message, double** omega, double* degrees, double (*model)(const int, const double, double&, double&, double&), string out_file_name) {
  ofstream out_file;
  out_file.open((out_file_name + "_posterior.edges").c_str());

  for (unsigned long i = 0; i < edges; i++) {
    double edge_posterior = 0;

    if (EdgeList[i].known != 0) {
      double denom = 0; // a normalizing factor, for computing q_{ij}
      for (unsigned j = 0; j < communities; ++j) {
        for (unsigned k = 0; k < communities; ++k) {
          denom += SBM_joint_marginal(omega[j][k], EdgeList[i].known, current_message[0][i][j], current_message[1][i][k]);
        }
      }
      for (unsigned j = 0; j < communities; ++j) {
        for (unsigned k = 0; k < communities; ++k) {
          double p = omega[j][k];
          double q = SBM_joint_marginal(omega[j][k], EdgeList[i].known, current_message[0][i][j], current_message[1][i][k]) / denom;

          double t = p * EdgeList[i].known / (p * EdgeList[i].known + edge_ratio * (1 - p) * (1 - EdgeList[i].known));
          edge_posterior += t * q;
        }
      }
    }
    out_file << EdgeList[i].out << " " << EdgeList[i].in << " " << edge_posterior << endl;
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
  unweighted_degree = false; // THIS IS A GLOBAL VARIABLE
  // True if we should run the edge prediction algorithm after running the normal community detection
  const bool predict_edges_flag = false;

  double (*model)(const int, const double, double&, double&, double&);
  void (*M_model_init)(Trio*, double***, double**, double**, double*);
  void (*M_model)(Trio*, double***, double**, double**, double*);
  void (*MF_Precompute)(double**, double**, double*, double*, bool*, double**, double(const int, const double, double&, double&, double&));
  double (*MF_Return)(double**, const int&, const int&);
  double (*Compute_Likelihood)(Trio*, double**, double*, double**, double*);

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



  // these are global variables!
  vertices = 0;
  edges = 0;

  std::random_device rd;
  std::mt19937 generator(rd());
  std::uniform_real_distribution<double> numgen(0.0,1.0);

  // local variables used only in this method
  int restarts, max_iterations, lines;
  string input_file_name, out_file_name;
  double converged_diff;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "show help")
    ("input_file,i", po::value<string>(&input_file_name)->default_value(""), "input file containing edge list")
    ("communities,c", po::value<unsigned>(&communities)->default_value(2), "number of communities to detect")
    ("lines,l", po::value<int>(&lines)->default_value(2), "number of columns in the edge file. allows one to specify a certain network")
    ("out_file,o", po::value<string>(&out_file_name)->default_value(""), "file to output communities to")
    ("EM_restarts,n", po::value<int>(&restarts)->default_value(10),"random restarts of the full EM algorithm")
    ("EM_max_iterations", po::value<int>(&max_iterations)->default_value(20),"cap on number of EM iterations")
    ("BP_convergence_thresh", po::value<double>(&message_converged_diff)->default_value(0.001), "BP is considered converged when no message changes by more than this")
    ("EM_convergence_thresh", po::value<double>(&converged_diff)->default_value(5e-3), "EM is considered converged when no parameter changes by more than this")
    ("zero_thresh", po::value<double>(&zero_thresh)->default_value(1e-50), "any number below this is considered to be zero");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
    
  /* show help and exit */
  if (vm.count("help") || input_file_name == "" || out_file_name == "") {
    cout << desc << "\n";
    return 1;
  }

  

  std::ifstream input_file;
  input_file.open(input_file_name);
  if (!input_file)
  {
    cout << "Error in opening file";
    cin.get();
    return 0;
  }
  FindStats(vertices, edges, input_file);
  input_file.close();

  Trio* EdgeList = new Trio[edges];
  // This is the OBSERVED degree. The actual expected degree could be something different, depending on the number of unobserved vertex pairs.
  double* degrees = new double[vertices];
  double* missing_degrees = new double[vertices];
  // No vertex can have degree higher than the number of vertices
  bool* degrees_present = new bool[vertices];

  for (unsigned long i = 0; i < vertices; i++)
  {
    degrees[i] = 0;
    degrees_present[i] = false;
    missing_degrees[i] = 0;
  }

  GetTheNetworkEdges(input_file_name, lines, EdgeList, degrees, missing_degrees);

  // compute C_global and m_global
  edge_ratio = GetEdgeRatio(degrees);
  edge_sum = GetEdgeSum(degrees);



  // This is only useful if the degrees are already integers (i.e. the network is unweighted).
  for (unsigned long i = 0; i < vertices; i++)
  {
    // The addition of this ceil term is only to handle weighted 'histogram style' mean field approximating
    // I've removed it for the non-mean-field
    //degrees_present[int(ceil(degrees[i]))] = true;
    degrees_present[int(degrees[i])] = true;
  }





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
  for (unsigned long i = 0; i < vertices; i++)
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
  for (unsigned long i = 0; i < edges; i++)
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
  for (unsigned j = 0; j < communities; j++)
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
    for (unsigned j = 0; j < communities; j++)
    {
      nKcount[j] = 0;
    }

    
    for (unsigned long i = 0; i < vertices; i++)
    {
      norm1 = 0;
      for (unsigned j = 0; j < communities; j++)
      {
        group_membership[i][j] = 1; // numgen(generator);
        norm1 += group_membership[i][j];
      }
      for (unsigned j = 0; j < communities; j++)
      {
        group_membership[i][j] /= norm1;
        nKcount[j] += group_membership[i][j];
      }
    }

    for (unsigned j = 0; j < communities; j++)
    {
      nKcount[j] /= double(vertices);
    }

    /*
    for (unsigned j = 0; j < communities; j++)
    {
      nKcount[j] = .5;// /= double(vertices);
    }*/

    for (unsigned long j = 0; j < edges; j++)
    {
      norm1 = 0;
      norm2 = 0;
      for (unsigned k = 0; k < communities; k++)
      {
        former_message[0][j][k] = numgen(generator);
        former_message[1][j][k] = numgen(generator);

        //current message needs to be reset each time
        current_message[0][j][k] = 0;
        current_message[1][j][k] = 0;

        norm1 += former_message[0][j][k];
        norm2 += former_message[1][j][k];
      }
      for (unsigned k = 0; k < communities; k++)
      {
        former_message[0][j][k] /= norm1;
        former_message[1][j][k] /= norm2;
      }
    }

    // why is omega set to 0 each time?? this gives problems
    // in the first run of compute maxes, which here is essentially
    // just skipping things
    for (unsigned j = 0; j < communities; j++)
    {
      for (unsigned k = 0; k < communities; k++)
      {
        omega[j][k] = 0;
      }
    }

    // Start with the M step. Except this should be the non-BP version, whereas the full EM algorithm one should include the messages as expected.
    double converged = Compute_Maxes(false, EdgeList, current_message, group_membership, nKcount, omega, degrees, M_model_init);
    printf("Init omega!\n");
    for(unsigned i=0; i < communities; i++) {
      for(unsigned j=0; j < communities; j++) {
        cout << omega[i][j] << "\t";
      }
      cout << "\n";
    }
    cout << "gamma\n";

    for(unsigned i = 0; i < communities; i++) {
      cout << nKcount[i] << "\t";
    }
    cout << endl;

    converged = converged_diff * 4;
    for (unsigned i = 0; i < communities; i++)
    {
      for (unsigned j = 0; j < communities; j++)
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
        BP_algorithm(EdgeList, group_membership, general_message, current_message, former_message, nKcount, omega, degree_only, degrees, missing_degrees, degrees_present, model, MF_Precompute, MF_Return);
        converged = Compute_Maxes(true, EdgeList, current_message, group_membership, nKcount, omega, degrees, M_model);
        for(unsigned i=0; i < communities; i++) {
          for(unsigned j=0; j < communities; j++) {
            cout << omega[i][j] << "\t";
          }
          cout << "\n";
        }
        cout << "gamma\n";

        for(unsigned i = 0; i < communities; i++) {
          cout << nKcount[i] << "\t";
        }
        cout << endl;

        EMiterations++;
      }

      LL = Compute_Likelihood(EdgeList, group_membership, nKcount, omega, degrees);
    }
    else
    {
      for (unsigned long i = 0; i < vertices; i++)
      {
        for (unsigned j = 0; j < communities; j++)
        {
          group_membership[i][j] = 1. / communities;
        }
      }

      LL = Compute_Likelihood(EdgeList, group_membership, nKcount, omega, degrees);
    }

    if (std::isnan(LL))
    {
      LL = -std::numeric_limits<float>::max();
    }

    if (thisiteration == 0)
    {
      best_LL = LL;
      for (unsigned j = 0; j < communities; j++)
      {
        for (unsigned long i = 0; i < vertices; i++)
        {
          best_groups[i][j] = group_membership[i][j];
        }

        for (unsigned k = 0; k < communities; k++)
        {
          best_omega[j][k] = omega[j][k];
        }
        for (unsigned long i = 0; i < edges; i++)
        {
          best_message[0][i][j] = current_message[0][i][j];
          best_message[1][i][j] = current_message[1][i][j];
        }
      }
    }
    if (thisiteration == 0 || LL > best_LL)
    {
      best_LL = LL;
      for (unsigned j = 0; j < communities; j++)
      {
        for (unsigned long i = 0; i < vertices; i++)
        {
          best_groups[i][j] = group_membership[i][j];
        }

        for (unsigned k = 0; k < communities; k++)
        {
          best_omega[j][k] = omega[j][k];
        }
        for (unsigned long i = 0; i < edges; i++)
        {
          best_message[0][i][j] = current_message[0][i][j];
          best_message[1][i][j] = current_message[1][i][j];
        }
      }
      ofstream out_file;
      out_file.open(out_file_name);
      if (out_file.is_open())
      {
        for (unsigned long i = 0; i < vertices; i++)
        {
          for (unsigned j = 0; j < communities; j++)
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

      if(predict_edges_flag == true)
        predict_edges(EdgeList, best_message, best_omega, degrees, model, out_file_name);
    }
    printf("LL: %f, Best LL: %f\n", LL, best_LL);
  }



  for (unsigned i = 0; i < communities; i++)
  {
    delete[] omega[i];
    delete[] best_omega[i];
  }
  delete[] omega;
  delete[] best_omega;
  for (unsigned long i = 0; i < vertices; i++)
  {
    delete[] group_membership[i];
    delete[] best_groups[i];
    delete[] general_message[i];
  }
  for (unsigned long i = 0; i < edges; i++)
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
    for (unsigned long i = 0; i < vertices; i++)
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
