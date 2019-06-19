#include <iostream>
#include <cstring>
#include <cmath>
#include <complex>
#include <vector>
#include <fstream>
using namespace std;

#include "fft.h"
#include "Beamform.h"

int main(){

  int experiment = 1;  //experiment 1: changing theta of noise, experiment 2 : chainging input snr
  vector<int> NoS; //number of sensor
  vector<double> theta, target_snr; //number of sensor, snr as variable
  for (unsigned i = 2; i<=5; i++){
    NoS.push_back(i);
  }
  for (int i = 1; i<=90; i = i + 10){
    theta.push_back(m_PI/180 * double(i));
  }
  for (int i = -10; i<=20; i = i + 2){
    target_snr.push_back(i);
  }
  if (experiment == 0){
    double snr = 10;
    Beamform bf;
    bf.get_signal(THETA_INPUT, m_PI/3, snr);
    vector<signal_t> beam_rec = bf.beamform_Rx( );
    //record the result
    ofstream f("Beamformed_result.txt");
    for (unsigned i = 0; i< beam_rec.size(); i++){
      f<<beam_rec[i]<<endl;
    }

    if(f.is_open()==true){
    	f.close();
    }

  }
  else if (experiment == 1){
    double snr = 10, snri;
    vector<vector<double>> results;
    for (unsigned i = 0 ; i < NoS.size(); i++){
      vector<double> resi; // result vector of ith index
      for (unsigned j = 0; j<theta.size(); j++){
        Beamform bf(NoS[i]); // create beamform that has N[i] sensors
        bf.get_signal(THETA_INPUT, theta[j], snr);
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

//record the result

    ofstream f("Experiment1_result.txt");
  	for (unsigned i = 0; i< NoS.size(); i++){
      f<<"\n\n=============n==========="<<NoS[i]<<"=========="<<endl;
  		for (unsigned j = 0; j< theta.size(); j++){
  			f <<"  (theta: "<< 10*j + 1 << " >> " <<results[i][j]<<endl;
  		}
  		f<<"\n";
  	}

  	if(f.is_open()==true)
  	{
  		f.close();
  	}
  }



  else if (experiment == 2){
    double  snri;
    vector<vector<double>> results;
    for (unsigned i = 0 ; i < NoS.size(); i++){
      vector<double> resi; // result vector of ith index
      for (unsigned j = 0; j<target_snr.size(); j++){
        Beamform bf(NoS[i]); // create beamform that has N[i] sensors
        bf.get_signal(THETA_INPUT, m_PI, target_snr[j]);
        vector<signal_t> beam_rec = bf.beamform_Rx( ); // beamforming signal

        //power calculation to get snr
        double pow_out_noise(0), pow_out_sig(0);
        for (unsigned m = 0; m < min(beam_rec.size(), bf.sgn_1d_origin.size()); m++){
          pow_out_sig += beam_rec[m] * beam_rec[m];
          pow_out_noise += (beam_rec[m] - bf.sgn_1d_origin[m])*(beam_rec[m] - bf.sgn_1d_origin[m]);
        }
        cout<<"N = "<<NoS[i]<<" SNR = "<<target_snr[j] <<":: power = "<<pow_out_sig<<"  "<<pow_out_noise<<"\n\n";
        snri = 10 * log10(pow_out_sig / pow_out_noise) ;
        resi.push_back(snri);
      }
      results.push_back(resi);
    }

//record the result
    ofstream f("Experiment2_result.txt");
    for (unsigned i = 0; i< NoS.size(); i++){
      f<<"\n\n=============n==========="<<NoS[i]<<"=========="<<endl;
    	for (unsigned j = 0; j< target_snr.size(); j++){
    		f <<"  (target_snr: "<< -10 + 2*int(j) << " >> " <<results[i][j]<<endl;
    	}
    	f<<"\n";
    }

    if(f.is_open()==true){
    	f.close();
    }
  }

  return 0;
}
