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


//  fft(&tmp,FORWARD);
//  cout<<tmp[3]<<endl;
/*  vector<complex<double>>::iterator iter;
  for (iter = tmp.begin();iter!=tmp.end();iter++){
    cout<<*iter<<endl;

  }*/
/*  fftrealloc(tmp, 30, nextPow2(30) );
  fft(tmp, nextPow2(30) , FORWARD);

  for (unsigned i = 0; i<32; i++){
    printf("%f + %fi \n" , tmp[i].real(), tmp[i].imag() );
  }
*/


  Beamform bf;
  bf.get_signal();
  vector<double> beam_rec;
  beam_rec = bf.beamform_Rx( );


/*
vector<complex<double>> tst(128);
tst[0] = 3;
tst[1] = 10;
tst[2] = 20;
tst[3] = 7;
tst[4] = 3;
tst[5] = 22;
fft(&tst, FORWARD);

for (unsigned i = 0 ; i < tst.size() ; i++){
  cout<<tst[i]<<endl;

}
cout<<endl;
fft(&tst,INVERSE);
#ifdef DEBUG
  for (unsigned i = 0 ; i < tst.size() ; i++){
    cout<<tst[i]<<endl;

  }
  cout<<endl;
#endif
*/
  return 0;
}
