#ifndef _BEAMFORM_H_
#define _BEAMFORM_H_
#include<vector>

#include "Dtype.h"
#include "PARAM.h"
using namespace std;



class Beamform {
private:
  unsigned n, type;
  double d, f, c, fs, theta ; // distance, frequency, speed, sampling freq, theta estimation
  vector<vector<signal_t>> sgn; //2d n x len signal used for beamforming
  vector<signal_t> sgn_1d_origin , sgn_beamformed; // original signal for simulation, beamformed result
  bool DoA;
public:
  Beamform();
  ~Beamform();
  Beamform(unsigned n, double d, double f, double c, double fs);
  void get_signal();
  void set_DoA(bool is_set);
  void status();
  double estimate_DoA();
  double gccphat(vector<signal_t> &x, vector<signal_t> &x_ref, double fs);
  vector<complex<double>> get_weight(double F, int Type);
  vector<double> beamform_Rx( );
  vector<vector<double>> beamform_Tx();


};




#endif
