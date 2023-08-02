#pragma once
#include "timedef.h"
#include <cublas_v2.h>

typedef unsigned int uint;
template<typename FLOAT_T>
struct DenseMatrix {
	int R, C;
	FLOAT_T* v;
	DenseMatrix(int R, int C) :R(R), C(C) {
		v = (FLOAT_T*)malloc(R * C * sizeof(FLOAT_T));
		for (int i = 0; i < R * C; i++) {
			v[i] = randreal();
			// v[i]=(FLOAT_T)1;
			// v[i]=(i/C);
		}
	}
	void del() {
		free(v);
	}
	void clear() {
		for (int i = 0; i < R * C; i++) v[i] = 0;
	}
	FLOAT_T* operator[](int i) {
		return v + i * C;
	}
	void println() {
		cout << "[";
		for (int i = 0; i < R; i++) {
			cout << '[';
			for (int j = 0; j < C; j++) {
				cout << v[i * C + j];
				if (j == C - 1) {
					if (i == R - 1) cout << "]]\n";
					else cout << "],\n";
				}
				else cout << ",";
			}
		}
		cout << "--------------------------------------\n";
	}
	void print() {
		cout << "[";
		for (int i = 0; i < R; i++) {
			cout << '[';
			for (int j = 0; j < C; j++) {
				cout << v[i * C + j];
				if (j == C - 1) {
					if (i == R - 1) cout << "]]";
					else cout << "],";
				}
				else cout << ",";
			}
		}
		cout << "--------------------------------------\n";
	}
	void printErr() {
		FLOAT_T s = 0;
		for (int i = 0; i < R; i++)for (int j = 0; j < C; j++) {
			s += fabs(v[i * C + j]);
		}
		cout << std::fixed << std::setprecision(5) << "sum=" << s << '\n';
	}
	void printErr(const DenseMatrix& d) {
		FLOAT_T s = 0;
		for (int i = 0; i < R; i++)for (int j = 0; j < C; j++) {
			s += fabs(v[i * C + j] - d.v[i * C + j]);
		}
		cout << std::fixed << std::setprecision(5) << "sum=" << s << '\n';
	}
};

__global__ void DenseMatMul(double* a, double* b, double* c, int n, int k, int m) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int j = blockDim.y * blockIdx.y + threadIdx.y;
	if (i < n && j < m) {
		double* val = c + i * m + j;
		*val = 0;
		for (int l = 0; l < k; l++) {
			*val += a[i * k + l] * b[l * m + j];
		}
	}
}


template<uint block_size>
__global__ void DenseMatMulBlock(double* a, double* b, double* c, int n, int k, int m) {
	int blockX = blockIdx.x;
	int blockY = blockIdx.y;
	int lux = blockX * block_size;
	int luy = blockY * block_size;

	int totala = n * k;
	int totalb = k * m;
	int ti = threadIdx.x, tj = threadIdx.y;
	double val = 0;


	for (int d = 0; d < k; d += block_size) {
		__shared__ double As[block_size][block_size];
		__shared__ double Bs[block_size][block_size];

		int offa = lux * k + d;
		int offb = d * m + luy;

		double* astart = a + offa;
		double* bstart = b + offb;

		if (offa + ti * k + tj < totala) As[ti][tj] = astart[ti * k + tj];
		if (offb + ti * m + tj < totalb) Bs[ti][tj] = bstart[ti * m + tj];

		__syncthreads();
		int eu = block_size < k - d ? block_size : k - d;
		for (int e = 0; e < eu; e++) {
			val += As[ti][e] * Bs[e][tj];
		}
		__syncthreads();
	}
	double* cstart = c + lux * m + luy;
	if (lux + ti < n && luy + tj < m) {
		cstart[ti * m + tj] = val;
	}
}

/*

warmUpTime:0
#
805128192 805128192 805128192
start multiplication
multiplication:0
copy back time:96735
yes
*/



void DenseMatMul(DenseMatrix<double> A, DenseMatrix<double> B, DenseMatrix<double> C) {
	double t0, t1, t2;
	int N = A.R, K = A.C, M = B.C;
	int size_a = N * K * sizeof(double);
	int size_b = K * M * sizeof(double);
	int size_c = N * M * sizeof(double);
	double* a;
	double* b;
	double* c;

	cout << "#\n";
	cout << size_a << ' ' << size_b << ' ' << size_c << endl;
	cudaMalloc(&a, size_a);
	cudaMalloc(&b, size_b);
	cudaMalloc(&c, size_c);


	if (!a) {
		cout << "!a\n";
		return;
	}
	if (!b) {
		cout << "!b\n";
		return;
	}
	if (!c) {
		cout << "!c\n";
		return;
	}
	cudaMemcpy(a, A.v, size_a, cudaMemcpyHostToDevice);
	cudaMemcpy(b, B.v, size_b, cudaMemcpyHostToDevice);

	const int block_size = 32;
	int BX = (N + block_size - 1) / block_size, BY = (M + block_size - 1) / block_size;
	dim3 thread_nums(block_size, block_size);
	dim3 block_nums(BX, BY);

	cout << "start multiplication\n";
	t0 = time();
	DenseMatMulBlock<block_size> << <block_nums, thread_nums >> > (a, b, c, N, K, M);
	t1 = time();
	cout << "multiplication:" << t1 - t0 << endl;
	cudaMemcpy(C.v, c, size_c, cudaMemcpyDeviceToHost);

	t2 = time();
	cout << "copy back time:" << t2 - t1 << endl;
	cudaFree(a);
	cudaFree(b);
	cudaFree(c);

}


void check(DenseMatrix<double> A, DenseMatrix<double> B, DenseMatrix<double> C) {
	cublasHandle_t handle;
	cublasCreate(&handle);
	double alpha = 1.0, beta = -1.0;
	double* a;
	double* b;
	double* c;
	int N = A.R;
	int K = A.C;
	int M = B.C;
	cout << M << " " << K << ' ' << N << endl;
	int size_a = N * K * sizeof(double);
	int size_b = K * M * sizeof(double);
	int size_c = N * M * sizeof(double);
	cudaMalloc(&a, size_a);
	cudaMalloc(&b, size_b);
	cudaMalloc(&c, size_c);

	cudaMemcpy(a, A.v, size_a, cudaMemcpyHostToDevice);
	cudaMemcpy(b, B.v, size_b, cudaMemcpyHostToDevice);
	cudaMemcpy(c, C.v, size_c, cudaMemcpyHostToDevice);

	cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, M, N, K, &alpha, b, M, a, K, &beta, c, M);
	cublasDestroy(handle);

	cudaMemcpy(C.v, c, size_c, cudaMemcpyDeviceToHost);


	cudaFree(a);
	cudaFree(b);
	cudaFree(c);

}