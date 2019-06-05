#ifndef _BEAMFORM_H_
#define _BEAMFORM_H_

#include "Dtype.h"
#include "PARAM.h"




class Beamform {
private:
  unsigned n, sgn_len;
  double d, f, c, fs;
  signal_t **sgn;
public:
  Beamform();
  ~Beamform();
  Beamform(unsigned n, double d, double f, double c, double fs);
  void get_signal();
  double estimate_DoA();
  double gccphat(signal_t* x, signal_t* x_ref, size_t N, double fs);

};


#endif
