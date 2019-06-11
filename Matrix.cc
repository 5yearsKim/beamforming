#include "Matrix.h"
// Function to get cofactor of A[p][q] in temp[][]. n is current
// dimension of A[][]
void getCofactor(vector<vector<complex<double>>> &A, vector<vector<complex<double>>> &temp, int p, int q, int n)
{
	int i = 0, j = 0;

	// Looping for each element of the matrix
	for (int row = 0; row < n; row++)
	{
		for (int col = 0; col < n; col++)
		{
			// Copying into temporary matrix only those element
			// which are not in given row and column
			if (row != p && col != q)
			{
				temp[i][j++] = A[row][col];

				// Row is filled, so increase row index and
				// reset col index
				if (j == n - 1)
				{
					j = 0;
					i++;
				}
			}
		}
	}
}

/* Recursive function for finding determinant of matrix.
n is current dimension of A[][]. */
complex<double> determinant(vector<vector<complex<double>>> &A, int n)
{
	complex<double> D = 0; // Initialize result

	// Base case : if matrix contains single element
	if (n == 1)
		return A[0][0];

	vector<vector<complex<double>>> temp; // To store cofactors
	for (int i = 0; i<n; i++){
		vector<complex<double>> v(n);
		temp.push_back(v);
	}

	complex<double> sign = 1; // To store sign multiplier

	// Iterate for each element of first row
	for (int f = 0; f < n; f++)
	{
		// Getting Cofactor of A[0][f]
		getCofactor(A, temp, 0, f, n);
		D += sign * A[0][f] * determinant(temp, n - 1);

		// terms are to be added with alternate sign
		sign = -sign;
	}

	return D;
}

// Function to get adjoint of A[N][N] in adj[N][N].
void adjoint(vector<vector<complex<double>>> &A, vector<vector<complex<double>>> &adj)
{
	int N = A.size();
	if (N == 1)
	{
		adj[0][0] = 1;
		return;
	}

	// temp is used to store cofactors of A[][]
	complex<double> sign = 1;

	vector<vector<complex<double>>> temp; // To store cofactors
	for (int i = 0; i<N; i++){
		vector<complex<double>> v(N);
		temp.push_back(v);
	}


	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			// Get cofactor of A[i][j]
			getCofactor(A, temp, i, j, N);

			// sign of adj[j][i] positive if sum of row
			// and column indexes is even.
			sign = ((i+j)%2==0)? 1: -1;

			// Interchanging rows and columns to get the
			// transpose of the cofactor matrix
			adj[j][i] = (sign)*(determinant(temp, N-1));
		}
	}
}

// Function to calculate and store inverse, returns false if
// matrix is singular
bool inverse(vector<vector<complex<double>>> &A, vector<vector<complex<double>>> &inverse)
{
	int N = A.size();
	// Find determinant of A[][]
	complex<double> det = determinant(A, N);
	if (abs(det)<  0.0000000001)
	{
		cout << "Singular matrix, can't find its inverse";
		return false;
	}

	vector<vector<complex<double>>> adj; // To store cofactors
	for (int i = 0; i<N; i++){
		vector<complex<double>> v(N);
		adj.push_back(v);
	}

	// Find adjoint
	adjoint(A, adj);

	// Find Inverse using formula "inverse(A) = adj(A)/det(A)"
	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++)
			inverse[i][j] = adj[i][j]/det;

	return true;
}



// Driver program

/*
int main()
{
	vector<vector<complex<double>>> A;


	vector<complex<double>> v1{11, 2.0, 3.0}, v2{47.0, 55.0, 6.0}, v3{7.0, 8.0, 9.0};
	A.push_back(v1);
	A.push_back(v2);
	A.push_back(v3);

	vector<vector<complex<double>>>  adj; // To store adjoint of A[][]
	for (int i =0; i<3 ; i++){
		vector<complex<double>> v(3);
		adj.push_back(v);
	}
	vector<vector<complex<double>>>  inv; // To store adjoint of A[][]
	for (int i =0; i<3 ; i++){
		vector<complex<double>> v(3);
		inv.push_back(v);
	}


	cout << "Input matrix is :\n";
	display(A);

	cout << "\nThe Adjoint is :\n";
	adjoint(A, adj);
	display(adj);

	cout << "\nThe Inverse is :\n";
	if (inverse(A, inv))
		display(inv);

	return 0;
}
*/
