#include <iostream>
#include <cstring>
#include <cmath>
#include <complex>
#include <vector>

using namespace std;

#include "PARAM.h"
#include "fft.h"
#include "Beamform.h"
#include "Dtype.h"



int main(){

int experiment = 1;
vector<int> NoS; //number of sensor
vector<double> theta; //number of sensor
for (unsigned i = 2; i<=5; i++){
  NoS.push_back(i);
}
for (int i = 1; i<=90; i = i + 10){
  theta.push_back(m_PI/180 * double(i));
}
if (experiment == 1){
  double snr = 10, snri;
  vector<vector<double>> results;
  for (unsigned i = 0 ; i < NoS.size(); i++){
    vector<double> resi; // result vector of ith index
    for (unsigned j = 0; j<theta.size(); j++){
      Beamform bf(NoS[i]); // create beamform that has N[i] sensors
      bf.get_signal(m_PI/2, theta[j], snr);
      vector<signal_t> beam_rec = bf.beamform_Rx( ); // beamforming signal

      //power calculation to get snr
      double pow_out_noise(0), pow_out_sig(0);
      for (unsigned m = 0; m < min(beam_rec.size(), bf.sgn_1d_origin.size()); m++){
        pow_out_sig += beam_rec[m] * beam_rec[m];
        pow_out_noise += (beam_rec[m] - bf.sgn_1d_origin[m])*(beam_rec[m] - bf.sgn_1d_origin[m]);
      }
      snri = 10 * log10(pow_out_sig / pow_out_noise) - snr;
      resi.push_back(snri);
    }
    results.push_back(resi);
  }

  for (unsigned i = 0; i< NoS.size(); i++){
    cout<<"\n\n=============n==========="<<NoS[i]<<"=========="<<endl;
    for (unsigned j = 0; j< theta.size(); j++){
      cout <<"  (theta: "<< j << " >> " <<results[i][j];
    }
  }
}



/*
  Beamform bf;
  bf.get_signal( m_PI/2, m_PI/3 , 10);
  vector<double> beam_rec;
  beam_rec = bf.beamform_Rx( );
  bf.status();
*/

  return 0;
}
