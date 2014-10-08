#include <iostream>
#include <vector>
#include <math.h> /* modf, tgamma */
#include <float.h> /* FLT_MIN */

using namespace std;

class Segmenter {

  vector<int> seq;
  const int maxk = 1;
  float nullprior;
  vector<int> states_;

public:

  Segmenter (vector<int> &seq_, float nullprior_);
  Segmenter (vector<int> &seq_);

  // sequence states
  void load_states();
  vector<int> states() { return states_; };

  // probability calculations
  float marginal_likelihood_R();
  float prior_k(int k);
  float prob_R_given_k(int k);
  float prob_R_given_zero_k();
  float prob_R_given_unit_k();
  float prob_R_given_k_rhs(vector<int> &states, int length);
  float prob_k_given_R(int k);

};

Segmenter::Segmenter (vector<int> &seq_) {

  seq = seq_;
  nullprior = 0.5;
  load_states();

}

Segmenter::Segmenter (vector<int> &seq_, float nullprior_) {

  seq = seq_;
  nullprior = nullprior_;
  load_states();

}

void Segmenter::load_states() {

  states_.resize(24, 0);
  for (auto p = seq.begin(); p != seq.end(); ++p) {
    ++ (states_[(int)*p]);
  }

}

float Segmenter::marginal_likelihood_R() {

  float sum = 0.0;

  for (int k = 0; k < 2; ++k) {
    sum += prior_k(k) * prob_R_given_k(k);
  }

  return sum;

}

float Segmenter::prior_k(int k) {

  if (k == 0) {
    return nullprior;
  } else {
    return ((1.0 - nullprior) / maxk);
  }

}

float Segmenter::prob_R_given_k(int k) {

  if (k == 0) {
    return prob_R_given_zero_k();
  } else if (k == 1) {
    return prob_R_given_unit_k();
  }

  return 0.0;

}

float Segmenter::prob_R_given_zero_k() {

  return prob_R_given_k_rhs(states_, seq.size());

}

float Segmenter::prob_R_given_unit_k() {

  vector<int> lstates = states_;
  int total = seq.size();
  vector<vector<float> > pmat(total - 1, vector<float>(total - 1, 0.0));

  for (int i = 0; i < total; ++i) {
    vector<int> segstates = lstates;
    for (int j = total - 1; j > 0; --j) {
      // get the probability of this segment given k=0
      float segresult = prob_R_given_k_rhs(segstates, j + 1 - i);
      // save it for lookup
      pmat[i][j] = segresult;
      // decrease the count of the last base in this segment
      --segstates[seq[j]];
    }
    // decrease the count of the first base in this sequence
    --segstates[seq[i]];
  }
  // we have a flat prior on observing any given segmentation
  // to avoid zero priots, we have a lower limit of FLT_MIN
  float pA = max((float)(1.0 / (float)total - 1), FLT_MIN);
  // now we calculate the probability for k=1
  float result = 0.0;
  for (int i = 0; i < total; ++i) {
    // probability for the sequence to the left of the change point
    float left = pmat[0][i];
    // probability for the sequence to the right of the change point
    float right = pmat[i + 1][total - 1];
    // probability for this segmentation
    result += left * right * pA;
  }
  return result;
}

float Segmenter::prob_R_given_k_rhs(vector<int> &states, int length) {

  int nstates = 0;

  float upper = 1;
  for (int i = 0; i < 24; ++i) {
    if (states[i] > 0) {
      ++nstates;
      upper *= states[i];
    }
  }

  float l = tgamma(nstates);
  float lower = tgamma(length + nstates);

  return l * upper / lower;

}

float Segmenter::prob_k_given_R(int k) {

  float pRk = prob_R_given_k(k);
  float pk = prior_k(k);
  float mlR = marginal_likelihood_R();

  cout << "pRk: " << pRk << endl;
  cout << "pk: " << pk << endl;
  cout << "mlR: " << mlR << endl;
  return pRk * pk / mlR;

}

int main() {

  // create a test sequence with a single changepoint
  //                    changepoint ---->|
  vector<int> seq = { 0,1,1,0,1,0,0,1,1,0,3,4,5,3,4,5,3,5,3,4,3,5 };
  for (auto p : seq) {
    cout << p << ',';
  }
  cout << endl;

  Segmenter seg{seq, 0.5};

  vector<int> allstates = seg.states();
  for (auto p : allstates) {
    cout << p;
  }
  cout << endl;

  cout << "gamma(4): " << tgamma(4) << endl;
  cout << "gamma(5): " << tgamma(5) << endl;
  cout << "gamma(6): " << tgamma(6) << endl;
  cout << "gamma(30): " << tgamma(30) << endl;

  cout << "p(R|k=0): " << seg.prob_R_given_zero_k() << endl;
  float p0 = seg.prob_k_given_R(0);
  cout << "p(k=0): " << p0 << endl;
  float p1 = seg.prob_k_given_R(1);
  cout << "p(k=1): " << p1 << endl;

  return 0;

}
