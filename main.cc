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
/*   signal_t* tmp = new signal_t[30];
  signal_t* tmp2 = new signal_t[30];
  tmp[0] = 3;
  tmp[1] = 10;
  tmp[2] = 20;
  tmp2[7] = 3;
  tmp2[8] = 10;
  tmp2[9] = 20;
*/
/*  fftrealloc(tmp, 30, nextPow2(30) );
  fft(tmp, nextPow2(30) , FORWARD);

  for (unsigned i = 0; i<32; i++){
    printf("%f + %fi \n" , tmp[i].real(), tmp[i].imag() );
  }
*/
  Beamform bf;
  bf.get_signal();
  cout<<bf.estimate_DoA()<<endl;

  return 0;
}
