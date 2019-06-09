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

//helper function for simulation
void fft(vector<complex<double>> &sgn, int inv);
int nextPow2(int n); // getting next power of 2 for fft
vector<vector<complex<double>>> buffer_shape(vector<complex<double>> &sgn, size_t frm_len, size_t overlap); // forming 1D signal to 2D for STFT
vector<vector<complex<double>>> STFT(vector<complex<double>> &sgn, size_t frm_len, size_t overlap, vector<complex<double>> &wnd);
vector<complex<double>> hann(size_t n); // hanning window
vector<complex<double>> iSTFT(vector<vector<complex<double>>> &sgn, size_t frm_len, size_t overlap);


vector<vector<double>> gen_arr_sig(vector<double> &in_sgn, unsigned N, double D, double Theta, double C, double fs); // copying delayed 1D signal and generate 2D signal
double gccphat(vector<signal_t> &x, vector<signal_t> &x_ref, double fs); // get time delay of two input array

#endif
