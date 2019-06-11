#ifndef _MATRIX_H_
#define _MATRIX_H_
using namespace std;
// C++ program to find adjoint and inverse of a matrix
//#include<bits/stdc++.h>
#include<vector>
#include<complex>
#include<iostream>




void getCofactor(vector<vector<complex<double>>> &A, vector<vector<complex<double>>> &temp, int p, int q, int n) ;
complex<double> determinant(vector<vector<complex<double>>> &A, int n) ;
// Function to get adjoint of A[N][N] in adj[N][N].
void adjoint(vector<vector<complex<double>>> &A, vector<vector<complex<double>>> &adj) ;
// Function to calculate and store inverse, returns false if
// matrix is singular
bool inverse(vector<vector<complex<double>>> &A, vector<vector<complex<double>>> &inverse) ;

// Generic function to display the matrix. We use it to display
// both adjoin and inverse. adjoin is integer matrix and inverse
// is a float.

template<class T>
void display(T A)
{
	int N = A.size();
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
			cout << A[i][j] << " ";
		cout << endl;
	}
}

#endif
