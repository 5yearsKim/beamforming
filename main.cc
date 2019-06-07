#include <iostream>
#include <cstring>
#include <cmath>
#include <complex>
#include <vector>

using namespace std;

#include "fft.h"
#include "Beamform.h"
#include "PARAM.h"
#include "Dtype.h"



int main(){



  Beamform bf;
  bf.get_signal();
  vector<double> beam_rec;
  beam_rec = bf.beamform_Rx( );
  bf.status();


  return 0;
}
