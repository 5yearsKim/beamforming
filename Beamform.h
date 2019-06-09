#ifndef _BEAMFORM_H_
#define _BEAMFORM_H_
#include<vector>

#include "Dtype.h"
#include "PARAM.h"
using namespace std;



class Beamform {
private:
  unsigned n, type, noise_type;
  double d, f, c, fs, theta ; // distance, frequency, speed, sampling freq, theta estimation
  vector<vector<signal_t>> sgn; //2d n x len signal used for beamforming
  bool DoA;
  double estimate_DoA();
  vector<complex<double>> get_weight(double F, int Type);
public:
  vector<signal_t> sgn_1d_origin , sgn_beamformed; // original signal for simulation, beamformed result
  Beamform();
  ~Beamform();
  Beamform(unsigned m_n);
  void get_signal(double s_theta, double n_theta, double snr);
  void set_DoA(bool is_set);
  void status();
  vector<double> beamform_Rx( );
  vector<vector<double>> beamform_Tx();


};




#endif
