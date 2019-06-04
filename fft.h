#ifndef _fft_h_
#define _fft_h_

#include "PARAM.h"

void fft(complex<double> *sgn, int n, int inv);

void fft(complex<double> *sgn, int n, int inv) {
	int i, j, j1, j2, k, iter, irem,  s, m, n1, n2, nxp, mxp, nxp2;
	double arg, wr, wi, tr, ti, wpwr;

	double* fr = new double[n];
	double* fi = new double[n];

	for (unsigned q = 0; q <n; q++){
		fr[q] = sgn[q].real();
		fi[q] = sgn[q].imag();
	}


	iter = 0;
	irem = n;
	while ((irem / 2) >= 1) {
		iter++;
		irem /= 2;
	}
	s = (inv ? 1 : -1);
	nxp2 = n;
	for (i=0; i<iter; i++) {
		nxp = nxp2;
		nxp2 = nxp / 2;
		wpwr = m_PI / (double)nxp2;
		for (m=0; m<nxp2; m++) {
			arg = (double)m * wpwr;
			wr = cos(arg);
			wi = s * sin(arg);
			for (mxp=nxp; mxp<=n; mxp+=nxp) {
				j1 = mxp - nxp + m;
				j2 = j1 + nxp2;
				tr = fr[j1] - fr[j2];
				ti = fi[j1] - fi[j2];
				fr[j1] = fr[j1] + fr[j2];
				fi[j1] = fi[j1] + fi[j2];
				fr[j2] = tr * wr - ti * wi;
				fi[j2] = tr * wi + ti * wr;
			}
		}
	}
	n2 = n / 2;
	n1 = n - 1;
	j = 1;
	for (i=1; i<=n1; i++) {
		if (i < j) {
			tr = fr[j - 1];
			ti = fi[j - 1];
			fr[j - 1] = fr[i - 1];
			fi[j - 1] = fi[i - 1];
			fr[i - 1] = tr;
			fi[i - 1] = ti;
		}
		k = n2;
		while (k < j) {
			j -= k;
			k /= 2;
		} j += k;
	}
	if (inv) {
		for (i=0; i<n; i++) {
			fr[i] /= n;
			fi[i] /= n;
		}
	}

	for (unsigned l = 0; l<n ; l++){
		sgn[l] = fr[l] + fi[l] * 1i;
	}
	delete [] fr;
	delete [] fi;
}


#endif
