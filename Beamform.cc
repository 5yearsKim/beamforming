#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>

#include "Beamform.h"



Beamform::Beamform():
n(N_OF_SENSOR),type(TEST_TYPE), noise_type(NOISE), d(DISTANCE_INPUT), f(FREQ_INPUT), c(SPEED_INPUT), fs(FS_INPUT), theta(THETA_INPUT), DoA(DOA){

}
Beamform::~Beamform(){

}

Beamform::Beamform(unsigned m_n):
n(m_n),type(TEST_TYPE), noise_type(NOISE), d(DISTANCE_INPUT), f(FREQ_INPUT), c(SPEED_INPUT), fs(FS_INPUT), theta(THETA_INPUT), DoA(DOA){

}
// setting input signal for simulation
void Beamform::get_signal(double s_theta, double n_theta, double snr){
	signal_t s, noise_element;
	ifstream s_in("1d_signal.txt"), n_in;
	if (noise_type == AWGN){
		n_in.open("AWGN.txt");
	}
	else if (noise_type == OFDM){
		n_in.open("OFDM.txt");
	}
	else{
		cout<<"noise file open error"<<endl;
	}

	vector<signal_t> m_noise;
	unsigned len_max = 100000;
	signal_t pow_noise(0), pow_sig(0);
	if(s_in.is_open() && n_in.is_open()){
		while(!s_in.eof() && !n_in.eof()){
			s_in >> s;
			n_in >> noise_element;
			sgn_1d_origin.push_back(s);
			m_noise.push_back(noise_element);
			//power calculation for fitting snr
			pow_sig += s * s;
			pow_noise += noise_element * noise_element;
			if (sgn_1d_origin.size()>len_max){
				cout<<"get_signal:getting signal successful"<<endl;
				break;
			}
		}
	}
	//power setting for noise
	for (unsigned i = 0; i< m_noise.size(); i++){
		m_noise[i] *= sqrt(pow_sig/ pow_noise / pow(10, snr/10));
	}
	//generating signal - noise mixed 2D array for STFT
	cout<<"generating "<<n<<" x len"<<" 2D signal for simulation\n sig theta: = "<<s_theta<<"  noise theta: = "<<n_theta<<endl;
	sgn = gen_arr_sig(sgn_1d_origin, n , d ,s_theta, c,  fs);
	vector<vector<signal_t>> noise_2d =  gen_arr_sig(m_noise, n , d ,n_theta, c,  fs);
	//add noise to original signal for simulation
	for (unsigned i = 0; i < sgn.size() ; i++){
		for (unsigned j = 0; j< sgn[0].size(); j++){
			sgn[i][j] += noise_2d[i][j];
		}
	}

	if (sgn.size() >0 && sgn[0].size()>0){ cout<<"signal received successfully"<<endl;}

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

//show variables of class
void Beamform::status(){
	for (unsigned i = 0 ; i < 1000 ; i++){
		cout << i<<"  "<<sgn_1d_origin[i]<<"  "<<sgn_beamformed[i] << "  "<<endl;;
		if (!(i%10))
			cout<<endl;
	}
	cout<<"\n\n\n"<<"theta:"<<theta<<endl;


}

//estimate direction of source signal
double Beamform::estimate_DoA(){
	if (!sgn.size()){
		cout<<"error::get_signal first"<<endl;
		return 0;
	}
	vector<double> tau_est;
	for (unsigned i = 0; i<n-1 ; i++){
		tau_est.push_back(gccphat(sgn[i+1], sgn[i], fs));
		tau_est[i] /= i+1;
	}
	double max_tau = 0;
	for(unsigned i = 0 ; i < tau_est.size(); i++){
		if (abs(max_tau) < abs(tau_est[i])){
			max_tau  = tau_est[i];
		}
	}

	double Theta_est = acos(c * max_tau / d);

	return Theta_est;
}


//weight value is used for time shift of the signal
vector<vector<complex<double>>> Beamform::get_weight(vector<double> F_vec, int Type){
	vector<complex<double>> n_vec;
	complex<double> Imag(0,1);
	for (unsigned i = 0; i<n; i++){
		n_vec.push_back(complex<double>(i));
	}

	if (Type == FIXED){
		vector<complex<double>> W(n);
		double F(FS_INPUT) ;
		cout<<"TYPE is Fixed"<<endl;
		vector<complex<double>> spatial_sample(37) ;
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
		vector<vector<complex<double>>> A, A_H, temp, inv_temp ;
		vector<complex<double>> temp2;
		for (unsigned i = 0; i < n; i++){
			vector<complex<double>> v(37), u(37), v2(n), v3(n);
			A.push_back(v);
			A_H.push_back(u);
			temp.push_back(v2);
			inv_temp.push_back(v3);
		}
		for (unsigned i = 0 ; i<n ; i++ ) {
			for (unsigned j = 0; j <= 36; j++){
				A[i][j] = exp(-Imag* complex<double>(2) * complex<double>(m_PI) * n_vec[i] *complex<double>(F) * complex<double>(d) * cos(spatial_sample[j]) / complex<double>(c));
				A_H[i][j] = conj(A[i][j]);
			}
		}
		temp = A;
		for (unsigned i = 0 ; i<n ; i++ ) {
			for (unsigned j = 0; j < n; j++){
				for (unsigned k = 0; k <=36; k++){
					temp[i][j] += A[i][k] * A[j][k];
				}
			}
		}
		inverse(temp, inv_temp);
		complex<double> t;
		for (unsigned j = 0; j < n; j++){
			t = 0;
			for (unsigned k = 0; k <=36; k++){
				t += 	A[j][k] * r_d[k];
			}
			temp2.push_back(t);
		}

		for (unsigned j = 0; j < n; j++){
			for (unsigned k = 0; k < n; k++){
				W[j] += inv_temp[j][k]*temp2[k];
			}
		}
		complex<double> sum;
		for (unsigned i = 0; i < n; i++){
			sum += W[i];
		}
		for (unsigned i = 0; i < n; i++){
			W[i] /= sum;
		}
		vector<vector<complex<double>>> W_vec;
		for (unsigned i = 0; i < n; i++){
			vector<complex<double>> v(F_vec.size(), W[i]);
			W_vec.push_back(v);
		}

		return W_vec;


		}
	else{
		if (Type != DSB){
			cout<<"error! Type "<< Type<< "doesn't exist"<<endl;
		}
		//    W = (1 / n) .* exp(1i * 2 * pi * n_vec * F_vec * d * cos(theta) / c);
		vector<vector<complex<double>>> W_vec;
		complex<double> imag(0,1.0);
		for (unsigned i = 0; i<n; i++){
			vector<complex<double>> W;
			W_vec.push_back(W);
			for (unsigned j = 0 ; j < F_vec.size(); j++){
				W_vec[i].push_back(exp(imag * 2.0 * m_PI* F_vec[j] * d * cos(theta) * n_vec[i] / c )/double(n) );
//			cout<<i<<"th weight is "<<W[i]<<endl;
//			cout<<"cos(theta)="<<cos(theta)<<endl;
			}
		}
		return W_vec;
	}
}


//beamforming receive
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
			cout<<"estimated DOA is "<<theta_est<<endl;
		}
	}
	else {
		theta = m_PI/2;
	}
	unsigned w_len = 512;
	vector<double> F_vec;
	for (double i = 0.0; i <= w_len/2; i++){
		F_vec.push_back(i*fs/w_len);
	}
//get beamformer weight
	vector<vector<complex<double>>> W = get_weight(F_vec, type);

//Beamforming - compensating the delays of the signals
	vector<complex<double>> wnd = hann(w_len);
	vector<vector<complex<double>>> sig_stft, tmp;

	for (unsigned i = 0 ; i<n ; i++){
		tmp = STFT(m_sgn[i], w_len, w_len/2, wnd);
		//sig_stft = sig_stft + w(i) .* tmp;
		if ( !i ){ //first initialization
			sig_stft = tmp;
			for (unsigned x = 0; x < tmp.size(); x++){
				for(unsigned y = 0; y < F_vec.size(); y++){
					sig_stft[x][y] *= W[i][y];
				}
			}
		}
		else {
			for (unsigned x = 0; x < tmp.size(); x++){
				for(unsigned y = 0; y < F_vec.size(); y++){
					sig_stft[x][y] += W[i][y]*tmp[x][y];
				}
			}
		}
		for (unsigned x = 0; x < tmp.size(); x++){
			for(unsigned y = 1; y < F_vec.size() - 1; y++){
			 	sig_stft[x][w_len - y] = conj(sig_stft[x][y]);
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
		pow_out += sig_out[i] * sig_out[i];
  }
//	cout<<"power_in = "<<pow_in<< "  pow_out = "<< pow_out<<endl;
 //sig_istft = sig_istft .* sqrt(pow_in / pow_out);
	for (unsigned i = 0; i < sig_out.size() ; i++){
		sig_out[i] *= sqrt(pow_in/ pow_out);
  }
	sgn_beamformed = sig_out;
	return sig_out;
}

//beamforming transimission
vector<vector<double>> Beamform::beamform_Tx(){
	return  gen_arr_sig(sgn_beamformed, n, d,theta	, c, fs);
}
