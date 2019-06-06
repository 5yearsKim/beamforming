#ifndef _FFT_T_
#define _FFT_T_

//macro function, maximum and minimum
#define MIN(a,b) \
({ __typeof__ (a) _a = (a); \
__typeof__ (b) _b = (b); \
_a < _b ? _a : _b; })

#define MAX(a,b) \
({ __typeof__ (a) _a = (a); \
__typeof__ (b) _b = (b); \
_a > _b ? _a : _b; })
///////////////////////////////////////

#include <cstring>
#include <complex>
#include <vector>
#include <iostream>

using namespace std;

#include "PARAM.h"
#include "Dtype.h"

void fft(vector<complex<double>> &sgn, int inv);
int nextPow2(int n);
vector<vector<complex<double>>> buffer_shape(vector<complex<double>> &sgn, size_t frm_len, size_t overlap);
vector<vector<complex<double>>> STFT(vector<complex<double>> &sgn, size_t frm_len, size_t overlap, vector<complex<double>> &wnd);
vector<complex<double>> hann(size_t n);
vector<complex<double>> iSTFT(vector<vector<complex<double>>> &sgn, size_t frm_len, size_t overlap);

#endif
