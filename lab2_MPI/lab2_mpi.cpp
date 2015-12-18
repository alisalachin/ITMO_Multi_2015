// lab_pi_mpi.cpp: ���������� ����� ����� ��� ����������� ����������.
//

#include <mpi.h>
#include "stdafx.h" 
#include <iostream> 
#include <omp.h> 
#include <time.h> 
#include <fstream>
#include <vector>
#include <math.h>
using namespace std;


//const double eps = 0.1; ///< �������� ��������
vector<double> Jacobi(int N, vector<vector<double>> A, vector<double> F, vector<double> X, int argc, char** argv, double eps)
{
	//�������� �� ����������� ������
	if ((F.size()<1)||(A.size()<1)||(X.size()<1)){
		throw invalid_argument("Empty input");
	} else {
		if (A[0].size() < 1) {
			throw invalid_argument("Empty input");
		}
	}
	//�������� ���������� ��������
	if ((F.size() != A.size()) || (F.size() != A[0].size()) || (F.size() != X.size())) {
		throw invalid_argument("Sizes do not match");
	}
	if (F.size() == 1) {
		return F;
	}
	//�������� ������� ����������
	for (int i = 0; i < F.size(); i++) {   
		int sum=0;
		for (int j = 0; j < F.size(); j++) {
			if (i != j) {
				sum += A[i][j];
			}
		}
		if (A[i][i] <= sum) {
			throw invalid_argument("Matrix does not converge");
		}
	}
	int size;
	int id;
	int got_x = 0;
	double time;
	double norm; // �����, ������������ ��� ���������� �������� ��������� ������� ����� �������� ��������.
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	vector<double> TempX(N);
	if (id == 0) {
		time = MPI_Wtime();
	}
	do {
		for (int i = 0; i < N; i++) {
			TempX[i] = F[i];
			for (int g = 0; g < N; g++) {
				if (i != g)
					TempX[i] -= A[i][g] * X[g];
			}
			TempX[i] /= A[i][i];
		}
		norm = abs(X[0] - TempX[0]);
		for (int h = 0; h < N; h++) {
			if (abs(X[h] - TempX[h]) > norm)
			{
				norm = abs(X[h] - TempX[h]);
			}
			X[h] = TempX[h];
		}
		MPI_Bcast(&X[0], N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		if (norm <= eps) {
			got_x = 1;
			MPI_Bcast(&got_x, 1, MPI_INT, 0, MPI_COMM_WORLD);
		}
	//	cout << norm;
	} while (got_x==0);
	if (id == 0) {
		time = MPI_Wtime()-time;
		cout << time<<endl;
	}
	MPI_Finalize();
	return X;
}


int _tmain(int argc, char** argv)
{
	setlocale(LC_ALL, "Russian");
	time_t begin, end;
	int N, M, M2;
	
	if (argc != 5) {
		cout << "Incorrect number of arguments";
		return 0;
	}
	
	ifstream fin1(argv[1]);	//"in300.txt"
	ifstream fin2(argv[2]);//"in2_300.txt"
	double eps = atof(argv[3]);
	if ((eps <= 0)||(eps>=1)){
		cout << "eps must be greater than 0 and less than 1";
		return 0;
	}
	ofstream fout(argv[4]);//"out.txt"

	if (!fin1.is_open() || !fin2.is_open())
	{
		cout << "Files not found"<<endl;
		return 0;
	}
	fin1 >> M >> N;
	if (N != M + 1) {
		cout << "N!=M+1 or an incorrect file format"<<endl;
		return 0;
	}
	vector<vector<double>> matrix(M);						//������� ������� 
	vector<double> x(M);								
	vector<double> b(M);									//������ b (������� ������ ������)

	for (int i = 0; i < M; i++)
	{
		matrix[i].resize(M);
		for (int j = 0; j < M; j++) {
			fin1 >> matrix[i][j];
		}
		fin1 >> b[i];
	}
	fin2 >> M2;
	if (M != M2) {
		cout << "Number of elements in the files does not match" << endl;
		return 0;
	}
	
	for (int i = 0; i < M; i++)
	{
		fin2 >> x[i];
	}
	
	try {
		x = Jacobi(M, matrix, b, x, argc, argv, eps);
	}
	catch (invalid_argument e) {
		cout << e.what();
		return 0;
	}
	
	fout << M << endl;
	for (int i = 0; i < M; i++) {
		fout << x[i] << endl;
	}
	fin1.close();
	fin2.close();
	fout.close();
	return 0;
}