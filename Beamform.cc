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
	sgn_len = 0;
}
Beamform::~Beamform(){
	if (!sgn_len){
		for (unsigned i = 0 ; i<n; i++){
			delete [] sgn[i];
		}
		delete [] sgn;
		
	}
}

Beamform::Beamform(unsigned m_n, double m_d, double m_f, double m_c, double m_fs):
n(m_n), d(m_d), f(m_f), c(m_c), fs(m_fs) {

}

void Beamform::get_signal(){
	unsigned dlength = 10000;
	sgn_len = dlength;
	sgn = new signal_t*[n];
	for (unsigned i = 0 ; i<n; i++){
		sgn[i] = new signal_t[dlength];
	}

	FILE *pFile = fopen("noisy_signal.txt", "r");
	rewind(pFile);
	for (unsigned i =0; i<dlength; i++){
		fscanf(pFile, "%lf %lf %lf %lf", &sgn[0][i], &sgn[1][i],&sgn[2][i],&sgn[3][i] );
	}


	fclose(pFile);


}

double Beamform::estimate_DoA(){
	if (!sgn_len){
		cout<<"error::get_signal first"<<endl;
		return 0;
	}
	double tau_est(0);
	for (unsigned i = 0; i<n-1 ; i++){
		tau_est += gccphat(sgn[i+1], sgn[i], sgn_len, fs);
	}
	tau_est /= (n-1);

	double Theta_est = acos(c * tau_est / d);

	return Theta_est;
}

double Beamform::gccphat(signal_t* x, signal_t* x_ref, size_t N, double fs){
	unsigned new_N = nextPow2(N);
	complex<double>* comp_x = new complex<double>[new_N];
	complex<double>* comp_x_ref = new complex<double>[new_N];

	for (unsigned i = 0; i<N; i++){
		comp_x[i] = x[i];
		comp_x_ref[i] = x_ref[i];
	}
	fft(comp_x, new_N, FORWARD);
	fft(comp_x_ref,new_N,FORWARD);
	for (unsigned i = 0; i<new_N; i++){
		comp_x[i] = comp_x[i]*conj(comp_x_ref[i]);
	}
	fft(comp_x , new_N, INVERSE);

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

	delete[] comp_x;
	delete[] comp_x_ref;

	return (double)i_max / fs;
}
