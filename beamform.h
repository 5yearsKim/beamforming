#include "Dtype.h"
#include "PARAM.h"


#ifndef _BEAMFORM_H_
#define _BEAMFORM_H_

class beamform {
private:
  int n;
  double d, f, c, fs;
public:
  beamform();
  beamform(int n, double d, double f, double c, double fs);
  ~beamform();
  void initialize();

};


#endif
