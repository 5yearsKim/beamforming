#ifndef _BEAMFORM_H_
#define _BEAMFORM_H_
#include<vector>

#include "Dtype.h"
#include "PARAM.h"




class Beamform {
private:
  unsigned n;
  double d, f, c, fs;
  vector<vector<signal_t>> sgn;
public:
  Beamform();
  ~Beamform();
  Beamform(unsigned n, double d, double f, double c, double fs);
  void get_signal();
  double estimate_DoA();
  double gccphat(vector<signal_t>* x, vector<signal_t>* x_ref, double fs);

};


#endif
