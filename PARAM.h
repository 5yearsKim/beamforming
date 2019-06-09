#ifndef _PARAM_H_
#define _PARAM_H_

//parameter list
#define m_PI 3.141592
#define FORWARD 0
#define INVERSE 1

#define DSB 0
#define FIXED 1

#define AWGN 0
#define OFDM 1

//simulation default parameter value
#define THETA_INPUT m_PI/2 * 100/90   // used when generating signal for simulation
#define FS_INPUT 44100
#define SPEED_INPUT 340
#define FREQ_INPUT 7200
#define DISTANCE_INPUT 0.02
#define N_OF_SENSOR 4
#define DOA true // whether simulation estimate direction or not

#define NOISE AWGN 


#endif
