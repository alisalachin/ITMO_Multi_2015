// Matrix_OpenMP.cpp: определяет точку входа для консольного приложения.
//

#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <time.h> 
using namespace std;


class Matrix{
private:
	vector<vector<double>> vec;
public:
	/*Matrix(vector<vector<double>> vecIn)
	{
		//проверяем, что одинаковый размер у строк: 
		for (int i = 0; i < vecIn.size(); i++)
		{
			if (vecIn[i].size() != vecIn[0].size())
			{
			//	cout << vecIn[i].size() << "  " << vecIn[0].size() << endl;
				cout << "Неккоректные входные данные для создания матрицы";
				return;
			}
		}

		vec.resize(vecIn.size());
		vector<vector<double>>::iterator iterVecIn = vecIn.begin();
		copy(vecIn.begin(), vecIn.end(), vec.begin());
	}
	*/
	Matrix()
	{

	}

	Matrix operator+(const Matrix& vec_other)
	{
		Matrix vec_res;
		vec_res.vec.resize(this->vec.size());
		if ((vec_other.vec.size() != this->vec.size()) || (vec_other.vec[0].size() != this->vec[0].size()))
		{
			throw exception("Addition impossible. Sizes of the matrices are not suitable.");
			//cout << "Сложение не возможно. Размеры матриц не совпадают.";
			return vec_res;
		}
		for (int i = 0; i < this->vec.size(); i++)
		{
			for (int j = 0; j < this->vec[i].size(); j++)
			{
				vec_res.vec[i].push_back(this->vec[i][j] + vec_other.vec[i][j]);
			}
		}
		return vec_res;
	}

	Matrix operator-(const Matrix& vec_other)
	{
		Matrix vec_res;
		vec_res.vec.resize(this->vec.size());
		if ((vec_other.vec.size() != this->vec.size()) || (vec_other.vec[0].size() != this->vec[0].size()))
		{
			throw exception("Subtraction is impossible. Sizes of the matrices are not suitable.");
		//	cout << "Вычитание не возможно. Размеры матриц не совпадают";
			return vec_res;
		}
		for (int i = 0; i < this->vec.size(); i++)
		{
			for (int j = 0; j < this->vec[i].size(); j++)
			{
				vec_res.vec[i].push_back(this->vec[i][j] - vec_other.vec[i][j]);
			}
		}
		return vec_res;
	}

	Matrix operator*(const Matrix& vec_other)
	{
		Matrix vec_res;
		vec_res.vec.resize(this->vec.size());
		if (vec_other.vec.size() != this->vec[0].size())
		{
			throw exception("Multiplication is impossible.Sizes of the matrices are not suitable.");
		//	cout << "Умножение не возможно. Размеры матриц не совпадают.";
			return vec_res;
		}
		if (this->vec.size() >= vec_other.vec[0].size())
		{
		
#pragma omp parallel 
			{
#pragma omp for ordered schedule(dynamic,1)
				for (int i = 0; i < this->vec.size(); i++)
				{
					vec_res.vec[i].resize(vec_other.vec[0].size());
					for (int j = 0; j < vec_other.vec[0].size(); j++)
					{
						double val = 0;
						for (int z = 0; z < this->vec[0].size(); z++)
						{
							val += this->vec[i][z] * vec_other.vec[z][j];
						}
						vec_res.vec[i][j] = val;
					}
				}
			}
		}
		else
		{
		
#pragma omp parallel 
			{
#pragma omp for ordered schedule(dynamic,1)
				for (int i = 0; i < this->vec.size(); i++)
				{
					vec_res.vec[i].resize(vec_other.vec[0].size());
				}
			}
#pragma omp parallel 
			{
#pragma omp for ordered schedule(dynamic,1)
				for (int j = 0; j < vec_other.vec[0].size(); j++)
				{
					for (int i = 0; i < this->vec.size(); i++)
					{
						double val = 0;
						for (int z = 0; z < this->vec[0].size(); z++)
						{
							val += this->vec[i][z] * vec_other.vec[z][j];
						}
						vec_res.vec[i][j] = val;
					}
				}
			}
			
		}
		
		return vec_res;
	}

	friend ostream& operator<<(ostream& os, const Matrix& m)
	{
		vector<vector<double>>::const_iterator iterVec = m.vec.begin();
		for (iterVec; iterVec != m.vec.end(); iterVec++)
		{
			vector<double>::const_iterator iterVec2 = (*iterVec).begin();
			for (iterVec2; iterVec2 != (*iterVec).end(); iterVec2++)
			{
				os << *iterVec2 << " ";
			}
			os << endl;
		}
		return os;
	}
	friend istream& operator>>(istream& is, Matrix& m)
	{
		char str[1000000];
		char* nums;	
		m = Matrix();
		while (is.getline(str, 1000000))
		{
			nums = strtok(str, " ");
			vector<double> vec1;
			while (nums != NULL) {
				vec1.push_back(atof(nums));
				nums = strtok(NULL, " ");
			}
			m.vec.push_back(vec1);
		}
		//проверяем, что одинаковый размер у строк: 
		for (int i = 0; i < m.vec.size(); i++)
		{
			if (m.vec[i].size() != m.vec[0].size())
			{
				throw exception("lines in matrix hasn't same size");
				return is;
			}
		}

		return is;
	}
};



int main(int argc, char *argv[])
{

	setlocale(LC_ALL, "Russian");
	if (argc != 4) {
		cout << "Неверное количество агруметов";
		return 0;
	}
	
	ifstream fin1(argv[1]);
	ifstream fin2(argv[2]);
	ofstream fout(argv[3]);

	int N, M, N2, M2;
	clock_t begin, end;
	
	//cout << argv[1] << argv[2] << endl;
	
	/*
	
	ofstream fout("out.txt");
	ifstream fin1("1.txt");
	ifstream fin2("2.txt");*/

	if ((!fin1.is_open()) || (!fin2.is_open())) 
	{
		cout << "Файл(ы) входной(ые) не найден(ы)";
		return 0;
	}
	
	Matrix res1;
	try {
		fin1 >> res1;
		//cout << res1 << endl;
		Matrix res2;
		fin2 >> res2;
		//cout << res2<< endl;
	
	begin = clock();
	Matrix matrix3 = res1 * res2;							// перемножаем матрицы 
	end = clock();
	fout << matrix3 << endl;
	cout << (end - begin) / (double)CLOCKS_PER_SEC<<endl;

	} catch (exception & e)
	{
		cout << e.what();
	}
	fin1.close();
	fin2.close();
	fout.close();
	
	return 0;
}