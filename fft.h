#ifndef _FFT_T_
#define _FFT_T_

#include <cstring>
#include <complex>

using namespace std;

#include "PARAM.h"
#include "Dtype.h"

void fft(complex<double> *sgn, int n, int inv);
int nextPow2(int n);
void fftrealloc(complex<double>* sgn, size_t N, size_t new_N );


#endif
