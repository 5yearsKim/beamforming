#include "Dtype.h"
#include "PARAM.h"


#ifndef _BEAMFORM_H_
#define _BEAMFORM_H_

class Beamform {
private:
  unsigned n;
  double d, f, c, fs;
  signal_t **sgn;
public:
  Beamform();
  ~Beamform();
  Beamform(unsigned n, double d, double f, double c, double fs);
  void get_signal();
  double estimate_DoA();
};


#endif
