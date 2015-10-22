

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
	Matrix(vector<vector<double>> vecIn)
	{
		//проверяем, что одинаковый размер у строк:
		for (int i = 0; i < vecIn.size(); i++)
		{
			if (vecIn[i].size() != vecIn[0].size())
			{
				cout << "Неккоректные входные данные для создания матрицы";
				return;
			}
		}

		vec.resize(vecIn.size());
		vector<vector<double>>::iterator iterVecIn = vecIn.begin();
		copy(vecIn.begin(), vecIn.end(), vec.begin());
	}

	Matrix()
	{

	}

	Matrix operator+(const Matrix& vec_other)
	{
		Matrix vec_res;
		vec_res.vec.resize(this->vec.size());
		if ((vec_other.vec.size() != this->vec.size()) || (vec_other.vec[0].size() != this->vec[0].size()))
		{
			cout << "Сложение не возможно. Размеры матриц не совпадают.";
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
			cout << "Вычитание не возможно. Размеры матриц не совпадают.";
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
			cout << "Умножение не возможно. Размеры матриц не совпадают.";
			return vec_res;
		}
		if (this->vec.size() >= vec_other.vec[0].size())
		{
		
			#pragma omp parallel 
			{
				#pragma omp for ordered schedule(guided,1)
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
				#pragma omp for ordered schedule(guided,1)
				for (int i = 0; i < this->vec.size(); i++)
				{
					vec_res.vec[i].resize(vec_other.vec[0].size());
				}
			}
			#pragma omp parallel 
			{
				#pragma omp for ordered schedule(guided,1)
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
};



int main(int argc, char *argv[] )
{
	
	setlocale(LC_ALL, "Russian");
	if (argc!=3) {
		cout << "Неверное количество агруметов";
		return 0;
	}
	int N, M, N2, M2;
	clock_t begin, end;
	//cout << argv[1] << argv[2] << endl;
	FILE* f1 = fopen(argv[1], "r");
	FILE* f2 = fopen(argv[2], "r");
	ofstream fout("out.txt");

	if ((f1 == NULL) || (f2 == NULL))
	{
		cout << "Файл(ы) входной(ые) не найден(ы)";
		return 0;
	}


	char str[1000000];
	char* nums;
	vector<vector<double>> myVec1;						//создаем первую матрицу 
	while (!feof(f1))
	{
		fgets(str, 1000000, f1);
		nums = strtok(str, " ");
		vector<double> vec;
		while (nums != NULL) {
			vec.push_back(atof(nums));
			nums = strtok(NULL, " ");
		}
		myVec1.push_back(vec);
	}
	Matrix matrix1(myVec1);
	//cout << matrix1;

	//cout << endl;

	vector<vector<double>> myVec2;						//создаем вторую матрицу 
	while (!feof(f2))
	{
		fgets(str, 1000000, f2);
		nums = strtok(str, " ");
		vector<double> vec;
		while (nums != NULL) {
			vec.push_back(atof(nums));
			nums = strtok(NULL, " ");
		}
		myVec2.push_back(vec);
	}
	Matrix matrix2(myVec2);
	//cout << matrix2 << endl;


	begin = clock();
	Matrix matrix3 = matrix1 * matrix2;							// перемножаем матрицы
	end = clock();
	fout << matrix3 << endl;
	cout << (end - begin)/(double)CLOCKS_PER_SEC;

	return 0;
}

