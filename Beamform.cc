#include <iostream>
#include <cmath>
#include <complex>

using namespace std;

#include "Beamform.h"
#include "fft.h"


Beamform::Beamform(){
	n = 4;
	d = 0.02;
	f = 7200;
	c = 340;
	fs = 44100;
}
Beamform::~Beamform(){

}

Beamform::Beamform(unsigned m_n, double m_d, double m_f, double m_c, double m_fs):
n(m_n), d(m_d), f(m_f), c(m_c), fs(m_fs) {

}

void Beamform::get_signal(){
	unsigned dlength = 10000;


	for (unsigned i = 0 ; i<n; i++){
		vector<signal_t> v(dlength);
		sgn.push_back(v);
	}

	FILE *pFile = fopen("noisy_signal.txt", "r");
	rewind(pFile);
	for (unsigned i =0; i<dlength; i++){
		fscanf(pFile, "%lf %lf %lf %lf", &sgn[0][i], &sgn[1][i],&sgn[2][i],&sgn[3][i] );
	}


	fclose(pFile);


}

double Beamform::estimate_DoA(){
	if (!sgn.size()){
		cout<<"error::get_signal first"<<endl;
		return 0;
	}
	double tau_est(0);
	for (unsigned i = 0; i<n-1 ; i++){
		tau_est += gccphat(&sgn[i+1], &sgn[i], fs);
	}
	tau_est /= (n-1);

	double Theta_est = acos(c * tau_est / d);

	return Theta_est;
}

double Beamform::gccphat(vector<signal_t>* x, vector<signal_t>* x_ref, double fs){
	size_t N = x->size();
	unsigned new_N = nextPow2(N);
	vector<complex<double>> comp_x(new_N) ;
	vector<complex<double>> comp_x_ref(new_N);

	for (unsigned i = 0; i<N; i++){
		comp_x[i] = (*x)[i];
		comp_x_ref[i] = (*x_ref)[i];
	}
	fft(&comp_x,  FORWARD);
	fft(&comp_x_ref,FORWARD);
	for (unsigned i = 0; i<new_N; i++){
		comp_x[i] = comp_x[i]*conj(comp_x_ref[i]);
	}
	fft(&comp_x , INVERSE);

	double max = 0;
	int i_max = 0;
	for (unsigned i = 0; i<new_N; i++){
		if (comp_x[i].real()>max){
			max = comp_x[i].real();
			i_max = i;
		}
	}
	if (i_max>N/2){
		i_max =  i_max - (int)new_N;
	}


	return (double)i_max / fs;
}
