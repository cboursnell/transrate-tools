#include <iostream>
#include <vector>
#include <math.h> /* modf, tgamma */
#include <float.h> /* FLT_MIN */

#ifndef __segmenter_h__
#define __segmenter_h__

using namespace std;

class Segmenter {

  vector<int> seq;
  const int maxk = 1;
  double nullprior;
  vector<int> states_;
  vector<double> pRk_;

public:

  Segmenter (vector<int> &seq_, double nullprior_);
  Segmenter (vector<int> &seq_);

  // sequence states
  void load_states();
  vector<int> states() { return states_; };

  // probability calculations
  double marginal_likelihood_R();
  double prior_k(int k);
  double prob_R_given_k(int k);
  double prob_R_given_zero_k();
  double prob_R_given_unit_k();
  double prob_R_given_k_rhs(vector<int> &states, int length);
  double prob_k_given_R(int k);

};

#endif
