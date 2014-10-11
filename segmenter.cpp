// make sure to compile with flag: -std=c++11

#include "segmenter.h"

using namespace std;

Segmenter::Segmenter (vector<int> &seq_) {

  seq = seq_;
  nullprior = 0.5;
  load_states();

}

Segmenter::Segmenter (vector<int> &seq_, double nullprior_) {

  seq = seq_;
  nullprior = nullprior_;
  load_states();

}

void Segmenter::load_states() {

  states_.resize(24, 0);
  for (auto p = seq.begin(); p != seq.end(); ++p) {
    ++ (states_[(int)*p]);
  }
  pRk_.resize(2, -1.0);

}

double Segmenter::marginal_likelihood_R() {

  double sum = 0.0;

  for (int k = 0; k < 2; ++k) {
    sum += prior_k(k) * prob_R_given_k(k);
  }

  return sum;

}

double Segmenter::prior_k(int k) {

  if (k == 0) {
    return nullprior;
  } else {
    return ((1.0 - nullprior) / maxk);
  }

}

double Segmenter::prob_R_given_k(int k) {

  float result = pRk_[k];

  if (result == -1.0) {
    if (k == 0) {
      result = prob_R_given_zero_k();
    } else if (k == 1) {
      result = prob_R_given_unit_k();
    }
    pRk_[k] = result;
  }

  return result;

}

double Segmenter::prob_R_given_zero_k() {

  if (pRk_[0] != -1.0) {
    return pRk_[0];
  }

  double result = prob_R_given_k_rhs(states_, seq.size());
  pRk_[0] = result;
  return result;

}

double Segmenter::prob_R_given_unit_k() {

  if (pRk_[1] != -1.0) {
    return pRk_[1];
  }

  vector<int> lstates = states_;
  int total = seq.size();
  vector<vector<double> > pmat(total, vector<double>(total, 0.0));

  for (int i = 0; i < total; ++i) {
    vector<int> segstates = lstates;
    for (int j = total - 1; j >= i; --j) {
      // get the probability of this segment given k=0
      double segresult = prob_R_given_k_rhs(segstates, j + 1 - i);
      // save it for lookup
      pmat[i][j] = segresult;
      // decrease the count of the last base in this segment
      --segstates[seq[j]];
    }
    // decrease the count of the first base in this sequence
    --lstates[seq[i]];
  }
  // we have a flat prior on observing any given segmentation
  // to avoid zero priors, we have a lower limit of FLT_MIN
  double pA = max((1.0 / (double)(total - 1)), (double)FLT_MIN);
  // now we calculate the probability for k=1
  double result = 0.0;
  for (int i = 0; i < total - 1; ++i) {
    // probability for the sequence to the left of the change point
    double left = pmat[0][i];
    // probability for the sequence to the right of the change point
    double right = pmat[i + 1][total - 1];
    // probability for this segmentation
    double product = left * right * pA;
    result += product;
  }
  pRk_[1] = result;
  return result;
}

double Segmenter::prob_R_given_k_rhs(vector<int> &states, int length) {
  int nstates = states.size();
  double l = tgamma(nstates);

  double upper = 1;
  for (auto p : states) {
    upper *= tgamma(p + 1);
  }

  double lower = tgamma(length + nstates);

  return l * (upper / lower);

}

double Segmenter::prob_k_given_R(int k) {

  double pRk = prob_R_given_k(k);
  double pk = prior_k(k);
  double mlR = marginal_likelihood_R();

  return pRk * pk / mlR;

}
