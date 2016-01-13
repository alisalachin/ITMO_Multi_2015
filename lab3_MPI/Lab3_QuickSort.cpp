#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <iostream> 
using namespace std;

double When();
int sort(const void* a, const void*x);
int getPivot(int* array, int size);


int main(int argc, char * argv[])
{
	double time;
	int nproc, iproc;
	int* vals = (int*)malloc(100000 * sizeof(int));
	
	MPI_Status status;
	MPI_Init(&argc, &argv);
	
	int i = 0;
	double starttime;
	int* pivot = (int*)malloc(sizeof(int));
	int* send = (int*)malloc(sizeof(int));
	int* recv = (int*)malloc(sizeof(int));

	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
	if (argc != 3) {
		if (iproc == 0) {
			cout << "Incorrect number of arguments";
		}
		MPI_Finalize();
		return 0;
	}
	ifstream fin(argv[1]);	// "in.txt"


	char str[1000000];
	fin.getline(str, 1000000);
	char* nums;
	nums = strtok(str, " ");
	int ii = 0;
	while (nums != NULL) {
		
		vals[ii]=atof(nums);
		nums = strtok(NULL, " ");
		ii++;
	}
	int size = ii;
	int mySize = size / nproc;
	int nprocWorld = 0;
	int iprocWorld = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocWorld);
	MPI_Comm_rank(MPI_COMM_WORLD, &iprocWorld);
	if (iprocWorld == 0) {
		time = MPI_Wtime();
	}

	int* curVals = (int*)malloc(size * sizeof(int));
	
	int num = 0;
	for (int i = mySize*iproc; i < mySize*iproc + mySize; i++){
		curVals[num] = vals[i];
		num++;
	}
	int* theirs = (int*)malloc(size * sizeof(int));
	int* tmp = (int*)malloc(size * sizeof(int));
	int* under = (int*)malloc(size * sizeof(int));
	int* over = (int*)malloc(size * sizeof(int));

	qsort(curVals, mySize, sizeof(int), sort);
	MPI_Comm newcomm;
	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Comm_split(MPI_COMM_WORLD, 1, iproc, &newcomm);

	
	int groupSize = nproc;
	starttime = When();
	int numOfIter = 0;
	if (groupSize == 1){
		ofstream fout(argv[2]);//"out2.txt"
		for (int i = 0; i < size; i++){
			fout << *(curVals + i) << " ";
		}
		fout.close();
	}
	
	while (groupSize>1) {
		MPI_Comm_size(newcomm, &nproc);
		MPI_Comm_rank(newcomm, &iproc);

		*pivot = getPivot(curVals, mySize);

			if (iproc == 0){
				for (i = 1; i<nproc; i++) {
					MPI_Recv(recv, 1, MPI_INT, i, 0, newcomm, &status);
					tmp[i] = *recv;
				}
				tmp[0] = *pivot;
				qsort(tmp, nproc, sizeof(int), sort);
				*pivot = tmp[(nproc / 2)];
				for (i = 1; i<nproc; i++) {
					MPI_Send(pivot, 1, MPI_INT, i, 0, newcomm);
				}

			}
			else{
				MPI_Send(pivot, 1, MPI_INT, 0, 0, newcomm);
				MPI_Recv(pivot, 1, MPI_INT, 0, 0, newcomm, &status);
			}

		//вычисление количества отправляемых элементов

		*send = 0;
		for (i = 0; i<mySize; i++) {
			if (iproc < nproc / 2){
				if (curVals[i] >= *pivot) {
					tmp[*send] = curVals[i];
					(*send)++;

				}
			}
			else{
				if (curVals[i] < *pivot) {
					tmp[*send] = curVals[i];
					(*send)++;
				}
			}
		}
		
		//сначала посылают младшие

		if (iproc < nproc / 2){

			MPI_Send(send, 1, MPI_INT, iproc + (nproc / 2), 0, newcomm);
			MPI_Send(tmp, *send, MPI_INT, iproc + (nproc / 2), 0, newcomm);
			
			MPI_Recv(recv, 1, MPI_INT, iproc + (nproc / 2), 0, newcomm, &status);
			MPI_Recv(theirs, *recv, MPI_INT, iproc + (nproc / 2), 0, newcomm, &status);
		}
		else {

			MPI_Recv(recv, 1, MPI_INT, iproc - (nproc / 2), 0, newcomm, &status);
			MPI_Recv(theirs, *recv, MPI_INT, iproc - (nproc / 2), 0, newcomm, &status);

			MPI_Send(send, 1, MPI_INT, iproc - (nproc / 2), 0, newcomm);
			MPI_Send(tmp, *send, MPI_INT, iproc - (nproc / 2), 0, newcomm);
		}

		//соединяем с пришедшими элементами
		
		if (iproc<nproc / 2){
			mySize -= *send;

			for (i = 0; i<*recv; i++) {
				curVals[mySize] = theirs[i];
				mySize++;

			}
		}
		else{
			int off = 0;
			for (i = 0; i<mySize; i++) {
				if (curVals[i] >= *pivot){
					theirs[*recv + off] = curVals[i];
					off++;

				}
			}
			int* temp = curVals;
			curVals = theirs;
			theirs = temp;
			mySize = *recv + (mySize - *send);
		}
		
		qsort(curVals, mySize, sizeof(int), sort);
		
		groupSize /= 2;
		if (groupSize <= 1){
		//	MPI_Barrier(MPI_COMM_WORLD);
			if (iprocWorld == 0){
				int curProcSize = 0;
				int allProcSize = mySize;
				for (int i = 1; i < nprocWorld; i++) {
					MPI_Recv(&curProcSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
					MPI_Recv(theirs, curProcSize, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
					int cur = 0;
					for (int j = allProcSize; j < allProcSize + curProcSize; j++) {
						curVals[j] = theirs[cur];
						cur++;
					}
					allProcSize = allProcSize + curProcSize;
				}
			}
			else {
				MPI_Send(&mySize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
				MPI_Send(curVals, mySize, MPI_INT, 0, 0, MPI_COMM_WORLD);

			}
			//MPI_Barrier(MPI_COMM_WORLD);
			if (iprocWorld == 0){
				ofstream fout(argv[2]);//"out2.txt"
				for (int i = 0; i < size; i++){
					fout << *(curVals + i) << " ";
				}
				fout.close();
			}
		}

		MPI_Comm_split(newcomm, iproc < nproc / 2, iproc, &newcomm);
	}
	
	if (iprocWorld == 0) {
		time = MPI_Wtime() - time;
		cout << time << endl;
	}
	
	free(vals);
	free(theirs);
	free(tmp);
	free(send);
	free(recv);

	MPI_Finalize();

	return 0;
}

int getPivot(int* array, int size){
	if (size == 0) {
		return 0;
	}
	return array[size / 2];
}

int sort(const void *x, const void *y) {
	return (*(int*)x - *(int*)y);
}


double When()
{
	//struct timeval tp;
	//gettimeofday(&tp, NULL);
	return 0;//((double)tp.tv_sec + (double)tp.tv_usec * 1e-6);
}
