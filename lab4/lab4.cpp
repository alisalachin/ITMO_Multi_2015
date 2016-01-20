

#include "stdafx.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <string>
#include <mpi.h>

using namespace std;

vector<double> Dijkstra(vector<vector<double>> m, int s) {
	
	int n;
	int id;
	vector<double> finRes;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &n);
	double time;

	if (id == 0) {
		time = MPI_Wtime();
	}
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int size = m.size();
	vector<double> Vt1(size);
	fill(Vt1.begin(), Vt1.end(), 0);
	vector<double> Vt(size);
	for (int i = id; i < size; i += n) {
		if (m[i].size() != size) {
			cout << "Размер матрицы по высоте и ширине должен совпадать";
			return finRes;
		}
		for (int j = 0; j < size; j++){
			if (m[i][j] < 0){
				cout << "Веса не могут быть отрицательными";
				return finRes;
			}
			if ((i == j) && (m[i][j] != 0)){
				cout << "Веса на главной диагонали дожны быть равны 0";
				return finRes;
			}
			if ((i < j) && (m[i][j] != m[j][i])){
				cout << "Пути должны быть двунаправленными";
				return finRes;
			}
		}
	}
	if (id == 0) {
		if ((s >= size) || (s < 0)){
			cout << "Выбранной стартовой вершины не существует в графе";
			return finRes;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	vector<double> res1(size);
	fill(res1.begin(), res1.end(), 0);
	vector<double> res(size);
	for (int i = id; i < size; i+=n) { 
		
		if (i == s){
			res1[i]=0.0;
			Vt1[i] = INFINITY;
		}
		else {
			res1[i]=m[s][i];
			Vt1[i] = 1;
		}
	}
	MPI_Allreduce(&Vt1[0], &Vt[0], size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&res1[0], &res[0], size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	struct {
		double min;
		int   imin;
	} in, out;
	//cout << m[s][0] << " " << endl;
	int vtsize = 1;
	double min,min2;
	while (vtsize!=size)
	{
		min = INFINITY;
		in.min = INFINITY;
		in.imin = 0;
		for (int i = id; i < size; i+=n) {
			if (Vt[i] != INFINITY) {
				if (res[i] < in.min) {
					in.min = res[i];
					in.imin = i;
				}
			}
		}
		MPI_Allreduce(&in.min, &out.min, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
		int imin=out.imin;	
		Vt[imin] = INFINITY;
		vtsize++;
		for (int i = id; i < size; i+=n) {
			if (Vt[i] != INFINITY) {
			//	cout << res[i] <<" "<< res[imin] << " " << m[imin][i] << " "<<endl;
				if (res[i]>res[imin] + m[imin][i]) {
					res[i] = res[imin] + m[imin][i];
				}
			}
		}
		vector<double> res2(size);
		MPI_Allreduce(&res[0], &res2[0], size, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		res = res2;
	}
	if (id == 0) {
		time = MPI_Wtime() - time;
		cout << time << endl;
	}
	return res;


}

int main(int argc, char** argv) {
	setlocale(LC_ALL, "Russian");
	/*ofstream fout("matrix100.txt");
	for (int i = 0; i < 100; i++) {
		for (int j = 0; j < 100; j++) {
			if (i==j){
				fout << "0 ";
			}
			else {
				fout << "1 ";
			}
		}
		fout << endl;
	}*/
	if (argc != 4) {
		cout << "Неверное количество агруметов";
		return 0;
	}
	MPI_Init(&argc, &argv);
	ifstream fin1(argv[1]);//"matrix1.txt");
	int s = atoi(argv[2]);
	ofstream fout(argv[3]);
	if (!fin1.is_open())
	{
		cout << "Файл не найден";
		return 0;
	}
	
	
	vector<vector<double>> m;
	double d;
	char str[1000000];
	char* nums;
	while (fin1.getline(str, 1000000))
	{
		nums = strtok(str, " ");
		vector<double> vec1;
		while (nums != NULL) {
			//cout << nums << endl;
			if (!strncmp(nums, "INF",3)){
				vec1.push_back(INFINITY);
			}
			else {
				vec1.push_back(atof(nums));
			}
			nums = strtok(NULL, " ");
		}
		m.push_back(vec1);
	}
	/*for (int i = 0; i < m.size(); i++) {
		for (int j = 0; j < m[0].size(); j++){

			cout << m[i][j] << " ";
		}
		cout << endl;

	}*/
	vector<double> res = Dijkstra(m,s);
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	if (id == 0) {
		for (int i = 0; i < res.size(); i++) {
			if (!(res[i]==INFINITY)){
				fout << res[i] << " ";
			}
			else {
				fout << "INF ";
			}
		}
	}
	fin1.close();
	fout.close();
	MPI_Finalize();
	return 0;
}