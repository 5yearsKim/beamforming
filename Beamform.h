#ifndef _BEAMFORM_H_
#define _BEAMFORM_H_
#include<vector>

#include "Dtype.h"
#include "PARAM.h"




class Beamform {
private:
  unsigned n, type;
  double d, f, c, fs, theta ;
  vector<vector<signal_t>> sgn;
  bool DoA;
public:
  Beamform();
  ~Beamform();
  Beamform(unsigned n, double d, double f, double c, double fs);
  void get_signal();
  void set_DoA(bool is_set);
  double estimate_DoA();
  double gccphat(vector<signal_t> &x, vector<signal_t> &x_ref, double fs);
  vector<complex<double>> get_weight(double F, int Type);
  vector<double> beamform_Rx( );

};


#endif
