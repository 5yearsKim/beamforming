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
  vector<complex<double>> tmp(64);
  for (unsigned i = 0; i<64; i++){
    tmp[i] = i;
  }

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
  vector<complex<double>> wnd = hann(16);
//  vector<complex<double>> wnd;
//  wnd.assign(16,1);
  #ifdef DEBUG
  /*  for (unsigned i = 0 ; i < wnd.size() ; i++){
      cout<<wnd[i].real() <<"   ";

    }
    cout<<endl;*/
  #endif

  vector<vector<complex<double>>> res = STFT(&tmp,16,8, &wnd );
#ifdef DEBUG

  for (unsigned j =0; j<res.size(); j++){
    for (unsigned i =0; i<res[0].size(); i++){
      cout<<res[j][i]<<" ";
    }
    cout<<endl;
  }
#endif

  vector<complex<double>> ires = iSTFT(&res , 16, 8);




#ifdef DEBUG
  for (unsigned i = 0 ; i < ires.size() ; i++){
    cout<<ires[i]<<endl;

  }
  cout<<endl;
#endif

//  bf.get_signal();

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
