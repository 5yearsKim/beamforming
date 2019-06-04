#include <iostream>
#include <cmath>
#include <complex>

using namespace std;

#include "Beamform.h"


Beamform::Beamform(){
	n = 4;
	d = 0.02;
	f = 7200;
	c = 340;
	fs = 44100;
}
Beamform::~Beamform(){
	for (unsigned i = 0 ; i<n; i++){
		delete [] sgn[i];
	}
	delete [] sgn;

}

Beamform::Beamform(unsigned m_n, double m_d, double m_f, double m_c, double m_fs):
n(m_n), d(m_d), f(m_f), c(m_c), fs(m_fs) {

}

void Beamform::get_signal(){
	unsigned dlength = 10000;

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

double estimate_DoA(){


}
