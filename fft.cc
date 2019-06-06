#include "fft.h"


void fft(vector<complex<double>> &sgn, int inv) {
	int i, j, j1, j2, k, iter, irem,  s, m, n1, n2, nxp, mxp, nxp2;
	double arg, wr, wi, tr, ti, wpwr;
	int n = nextPow2((sgn).size());
	sgn.resize(n);
	double* fr = new double[n];
	double* fi = new double[n];

	for (int q = 0; q <n; q++){
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
			wi = double(s) * sin(arg);
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

	for (int l = 0; l<n ; l++){
		sgn[l] = fr[l] + fi[l] * 1i;
	}
	delete [] fr;
	delete [] fi;
}

int nextPow2(int n)
{
    if ( n <= 1 ) return n;
    double d = n-1;
    return 1 << ((((int*)&d)[1]>>20)-1022);
}


//sub function of STFT
vector<vector<complex<double>>> buffer_shape(vector<complex<double>> &sgn, size_t frm_len, size_t overlap){
	size_t N = sgn.size();
	if (frm_len<overlap){
		cout<<"underlap is not allowed for this simulation"<<endl;
		exit(1); //we don't consider underlap for this simulation
	}
	size_t row = N/ (frm_len-overlap);

	vector<vector<complex<double>>> buf;
	for (unsigned i = 0; i<row; i++){
			vector<complex<double>> v(frm_len);
			memcpy(&v[0], &sgn[(frm_len-overlap)*i] , min(frm_len, N - i*(frm_len-overlap)) * sizeof(complex<double>));
			buf.push_back(v);
	}
	return buf;

}

//Short Time Fourie Transform(signal, frame_lenght, overlap, window func )
vector<vector<complex<double>>> STFT(vector<complex<double>> &sgn, size_t frm_len, size_t overlap, vector<complex<double>> &wnd){
	if ( wnd.size() != frm_len){
		cout<<"window size not matched!"<<endl;
		exit(1);
	}
	vector<vector<complex<double>>> buf = buffer_shape(sgn, frm_len, overlap);
	vector<vector<complex<double>>>::iterator pos;
	for (pos = buf.begin(); pos != buf.end(); pos++){
		//mul(sgn .* wnd)
		for(unsigned i = 0; i<wnd.size(); i++){
			 (*pos)[i] = (*pos)[i] * wnd[i];
		}
		fft(*pos,FORWARD);
	}
	return buf;
}

//hanning window generation
vector<complex<double>> hann(size_t n){
	vector<complex<double>> wnd(n);
	for (unsigned i = 0; i<n; i++){
		wnd[i] = pow(sin(m_PI*i/n),2);
	}
	return wnd;
}

//inverse short time Fourie tramnsform
vector<complex<double>> iSTFT(vector<vector<complex<double>>> &sgn, size_t frm_len, size_t overlap){
	unsigned row = sgn.size();
	vector<complex<double>> Isgn((frm_len - overlap) * row + overlap);
	//strectchiing 2D vector to 1D vector
	for (unsigned i = 0 ; i < row; i++){
		fft(sgn[i], INVERSE);
		for (unsigned j = 0; j<frm_len; j++){
			#ifdef DEBUG2
			cout<<sgn[i][j]<<"  ";
			#endif
			Isgn[(frm_len - overlap)*i +j] += sgn[i][j];
		}
	}


	return Isgn;
}
