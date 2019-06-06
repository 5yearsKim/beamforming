#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>

#include "Beamform.h"
#include "fft.h"



Beamform::Beamform():
n(N_OF_SENSOR), d(DISTANCE_INPUT), f(FREQ_INPUT), c(SPEED_INPUT), fs(FS_INPUT), theta(THETA_INPUT), type(DSB), DoA(DOA){

}
Beamform::~Beamform(){

}

Beamform::Beamform(unsigned m_n, double m_d, double m_f, double m_c, double m_fs):
n(m_n), d(m_d), f(m_f), c(m_c), fs(m_fs) , DoA(true), type(DSB), theta(THETA_INPUT){

}

void Beamform::get_signal(){
	signal_t s;
	ifstream in("1Dnoisy_signal.txt");
	unsigned len_max = 100000;
	if(in.is_open()){
		while(!in.eof()){
			in >>s;
			sgn_1d_origin.push_back(s);
			if (sgn_1d_origin.size()>len_max)
				break;
		}
	}
	cout<<theta<<endl;
	sgn = gen_arr_sig(sgn_1d_origin, n , d ,theta, c,  fs);

	ofstream f("2Dnoisy_signal.txt");
	for (unsigned  j= 0; j<sgn[0].size(); j++){
		for (unsigned i=0; i<sgn.size(); i++){
			f<< sgn[i][j]<<"  ";
	//		cout<< sgn[i][j]<<"%%";
		}
		f<<"\n";
	}

	if(f.is_open()==true)
	{
		f.close();
	}
}

void Beamform::set_DoA(bool is_set){
	if (is_set == true)
		DoA = true;
	else if( is_set == false)
		DoA = false;
}

void Beamform::status(){
	for (unsigned i = 0 ; i < sgn_beamformed.size() ; i++){
		cout << sgn_beamformed[i] << "  ";
		if (!(i%10))
			cout<<endl;
	}
	cout<<"\n\n\n"<<"theta:"<<theta<<endl;


}

double Beamform::estimate_DoA(){
	if (!sgn.size()){
		cout<<"error::get_signal first"<<endl;
		return 0;
	}
	double tau_est(0);
	for (unsigned i = 0; i<n-1 ; i++){
		tau_est += gccphat(sgn[i+1], sgn[i], fs);
	}
	tau_est /= (n-1);

	double Theta_est = acos(c * tau_est / d);

	return Theta_est;
}

double Beamform::gccphat(vector<signal_t> &x, vector<signal_t> &x_ref, double fs){
	size_t N = x.size();
	unsigned new_N = nextPow2(N);
	vector<complex<double>> comp_x(new_N) ;
	vector<complex<double>> comp_x_ref(new_N);

	for (unsigned i = 0; i<N; i++){
		comp_x[i] = x[i];
		comp_x_ref[i] = x_ref[i];
	}
	fft(comp_x,  FORWARD);
	fft(comp_x_ref,FORWARD);
	for (unsigned i = 0; i<new_N; i++){
		comp_x[i] = comp_x[i]*conj(comp_x_ref[i]);
	}
	fft(comp_x , INVERSE);

	double max = 0;
	int i_max = 0;
	for (unsigned i = 0; i<new_N; i++){
		if (comp_x[i].real()>max){
			max = comp_x[i].real();
			i_max = i;
		}
	}
	if (i_max>(int)N/2){
		i_max =  i_max - (int)new_N;
	}


	return (double)i_max / fs;
}


vector<complex<double>> Beamform::get_weight(double F, int Type){
	vector<complex<double>> n_vec;

	for (unsigned i = 0; i<n; i++){
		n_vec.push_back(double(i));
	}

	vector<complex<double>> W(n);
	if (Type == FIXED){
	/*	vector<complex<double>> spatial_sample(37) ;
		for (unsigned i = 0 ; i<=36 ; i++ ) {
			spatial_sample[i] = double(i)/36 *m_PI;
		}
		vector<complex<double>> r_d(37) ;
		r_d[17] = 1;
		r_d[18] = 1;
		r_d[19] = 1.25;
		r_d[20] = 1;
		r_d[21] = 1;
		//A = exp(-1i * 2 * pi * n_vec * F * d * cos(spatial_sample) / c);
		vector<vector<complex<double>>> A;
		for (unsigned i = 0; i < n; i++){
			vector<complex<double>> v(37), u(37);
			A.push_back(v);
			A_H.push_back(u);
		}
		for (unsigned i = 0 ; i<n ; i++ ) {
			for (unsigned j = 0; j <= 36; j++){
				A[i][j] = exp(-1i * 2 * m_PI * n_vec[i] * F * d * cos(spatial_sample[j]) / c);
				A_H[i][j] = conj(A[i][j]);
			}
		}
		// challanged at matrix inverse.. giving up here

		*/
		}
	else{
		if (Type != DSB){
			cout<<"error! Type "<< Type<< "doesn't exist"<<endl;
		}
		//W = (1 / n) .* exp(1i * 2 * pi * n_vec * F * d * cos(theta) / c);
		complex<double> j(0,1);
		for (unsigned i = 0; i<n; i++){
			W[i] = exp(j * double(2) * m_PI* F * d * cos(theta) / c * n_vec[i])/double(n);
		}
	}
	return W;

}

vector<double> Beamform::beamform_Rx( ){

// converting signal_t to complex
	vector<vector<complex<double>>> m_sgn;
	for (unsigned i = 0; i<sgn.size(); i++){
		vector<complex<double>> v(sgn[0].size());
		m_sgn.push_back(v);
		for (unsigned j = 0; j < sgn[0].size(); j++ ){
			m_sgn[i][j] = sgn[i][j];
		}
	}

//setting theta to estimation value
	if (DoA){
		double theta_est = estimate_DoA();
		// Let theta be in range [80` 100`]
		if (theta_est < 80.0* m_PI / 180.0){
			theta = 80.0 * (m_PI / 180.0);
		}
		else if (theta_est > 100.0* m_PI / 180.0){
			theta = 100.0 * (m_PI / 180.0);
		}
		else {
			theta = theta_est;
		}
	}
	else {
		theta = m_PI/2;
	}

//get beamformer weight
	vector<complex<double>> W = get_weight(f, type);

//Beamforming - compensating the delays of the signals
	unsigned w_len = 512;
	vector<complex<double>> wnd = hann(w_len);
	vector<vector<complex<double>>> sig_stft, tmp;

	for (unsigned i = 0 ; i<n ; i++){
		tmp = STFT(m_sgn[i], w_len, w_len/2, wnd);
		//sig_stft = sig_stft + w(i) .* tmp;
		if ( !i ){
			sig_stft = tmp;
		}
		else {
			for (unsigned x = 0; x < tmp.size(); x++){
				for(unsigned y = 0; y < tmp[0].size(); y++){
					sig_stft[x][y] += tmp[x][y];
				}
			}
		}

	}
	//ISTFT to get beamformed signal
	vector<complex<double>> sig_istft = iSTFT(sig_stft,w_len, w_len/2);

	//Leave only real parts and overlap-and-add
	vector<double> sig_out(sig_istft.size());
	for (unsigned i = 0; i < sig_istft.size(); i++) {
		sig_out[i] = sig_istft[i].real();
	}
	for (unsigned i = 1 ; i < w_len/2; i++){
		sig_out[i] /=  wnd[i].real();
		sig_out[sig_out.size() - i] /= wnd[wnd.size() - i].real();
	}

 // Match the power of input and output signal
	double pow_in(0), pow_out(0);
 //pow_in = sum(Sig_arr(:, 1) .^ 2);
 //pow_out = sum(sig_istft .^ 2);
  for (unsigned i = 0; i < sgn[0].size() ; i++){
	 	pow_in += sgn[0][i] * sgn[0][i];
		pow_out = sig_out[i] * sig_out[i];
  }
 //sig_istft = sig_istft .* sqrt(pow_in / pow_out);
	for (unsigned i = 0; i < sig_out.size() ; i++){
		sig_out[i] *= sqrt(pow_in/ pow_out);
  }
	sgn_beamformed = sig_out;
	return sig_out;
}

//beamforming transimission
vector<vector<double>> Beamform::beamform_Tx(){
	return  gen_arr_sig(sgn_beamformed, n, d,theta, c, fs);
}
