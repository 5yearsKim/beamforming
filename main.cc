#include <iostream>
#include <cstring>
#include <cmath>
#include <complex>

using namespace std;
#include "fft.h"
#include "Beamform.h"
#include "PARAM.h"
#include "Dtype.h"

int main(){
/*  complex<double> temp[30];
  temp[0] = 0;
  temp[1] = 10;
  fft(temp, 30 , FORWARD);

  for (unsigned i = 0; i<30; i++){
    printf("%f + %fi \n" , temp[i].real(), temp[i].imag() );
  }

*/
  Beamform bf;
  bf.get_signal();


  return 0;
}
