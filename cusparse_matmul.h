#pragma once

#include <cusparse.h>
#define CHECK_ERROR(V) {\
	status=V;\
	if(status!=CUSPARSE_STATUS_SUCCESS){cout<<"error at:"<<__LINE__;}\
}

void check(SparseMatrix A, SparseMatrix B, DCSR* ret) {
	cusparseHandle_t handle = 0;
	cusparseStatus_t status;
	csrgemm2Info_t info = NULL;
	CHECK_ERROR(cusparseCreate(&handle));

	cusparseMatDescr_t descrA, descrB, descrC, descrD;
	cusparseCreateMatDescr(&descrA);
	cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO);

	cusparseCreateMatDescr(&descrB);
	cusparseSetMatType(descrB, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descrB, CUSPARSE_INDEX_BASE_ZERO);

	cusparseCreateMatDescr(&descrC);
	cusparseSetMatType(descrC, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descrC, CUSPARSE_INDEX_BASE_ZERO);

	cusparseCreateMatDescr(&descrD);
	cusparseSetMatType(descrD, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descrD, CUSPARSE_INDEX_BASE_ZERO);


	//////////////////the following described some usages:
	//////////////////note that cublas and cusparse are column prioritized

	/*

		DenseMatrix<double> MAT(8,4);
		MAT[0][1]=MAT[0][2]=0;
		MAT[1][3]=0;
		MAT[4][3]=0;
		MAT[5][2]=0;
		MAT[6][0]=MAT[6][1]=MAT[6][2]=0;

		int*nnz_RowA;
		int nnzA;
		cudaMalloc(&nnz_RowA, sizeof(int)*MAT.R);

		int size_a = sizeof(double)*MAT.R*MAT.C;
		double*a;
		cudaMalloc(&a, size_a);
		cudaMemcpy(a, MAT.v, size_a, cudaMemcpyHostToDevice);


		MAT.println();

		CHECK_ERROR(cusparseDnnz(handle, CUSPARSE_DIRECTION_COLUMN, MAT.C, MAT.R, descrA, a, MAT.C, nnz_RowA, &nnzA));
		cout<<"nnzA:"<<nnzA<<endl;
		int *nnz_row = new int [MAT.R];
		cudaMemcpy(nnz_row, nnz_RowA, sizeof(int)*MAT.R, cudaMemcpyDeviceToHost);
		cout<<"nnz_RowA:"<<endl;
		for(int i=0;i<MAT.R;i++){
			cout<<nnz_row[i]<<" ";
		}cout<<"~"<<endl;

		delete[](nnz_row);

		double*csrValA;
		int*csrIdxA;
		int*csrColPtrA;

		cudaMalloc(&csrValA, sizeof(double)*nnzA);
		cudaMalloc(&csrIdxA, sizeof(int)*MAT.R);
		cudaMalloc(&csrColPtrA, sizeof(int)*nnzA);

		CHECK_ERROR(cusparseDdense2csc(handle, MAT.C, MAT.R, descrA, a, MAT.C, nnz_RowA, csrValA, csrColPtrA, csrIdxA));

		int*csrIdx=new int[MAT.R];
		int*csrColPtr=new int[nnzA];
		double*csrVal=new double[nnzA];
		cudaMemcpy(csrIdx, csrIdxA, sizeof(int)*MAT.R, cudaMemcpyDeviceToHost);
		cudaMemcpy(csrColPtr, csrColPtrA, sizeof(int)*nnzA, cudaMemcpyDeviceToHost);
		cudaMemcpy(csrVal, csrValA, sizeof(double)*nnzA, cudaMemcpyDeviceToHost);

		for(int i=0;i<MAT.R;i++){
			cout<<csrIdx[i]<<' ';
		}
		cout<<endl;
		for(int i=0;i<nnzA;i++){
			cout<<csrColPtr[i]<<' ';
		}
		cout<<endl;
		for(int i=0;i<nnzA;i++){
			cout<<csrVal[i]<<' ';
		}
		cout<<endl;



		delete[]csrIdx;
		delete[]csrColPtr;
		delete[]csrVal;

	*/

	A.print();
	B.print();

	DCSR CSRA(A, 0, 0);
	DCSR CSRB(B, 0, 0);

	// CSRB.print();

	int n, m, k;
	int nnzA = A.nnz;
	int nnzB = B.nnz;
	int nnzC;
	int* nnzTotalDevHostPtr = &nnzC;
	int baseC;

	n = A.n;
	k = A.m;
	m = B.m;

	cusparseSetPointerMode(handle, CUSPARSE_POINTER_MODE_HOST);
	cusparseCreateCsrgemm2Info(&info);

	int* csrRowPtrC;
	int* csrRowPtrA;
	int* csrRowPtrB;
	int* csrRowPtrD;
	int* csrColIdxA;
	int* csrColIdxB;
	int* csrColIdxC;
	int* csrColIdxD;

	cudaMalloc(&csrRowPtrA, sizeof(int) * (n + 1));
	cudaMalloc(&csrRowPtrC, sizeof(int) * (n + 1));
	cudaMalloc(&csrRowPtrB, sizeof(int) * (k + 1));
	cudaMalloc(&csrRowPtrD, sizeof(int) * (n + 1));
	cudaMalloc(&csrColIdxA, sizeof(int) * (nnzA));
	cudaMalloc(&csrColIdxB, sizeof(int) * (nnzB));
	cudaMalloc(&csrColIdxD, sizeof(int));

	cudaMemcpy(csrRowPtrA, CSRA.offset, sizeof(int) * (n + 1), cudaMemcpyHostToDevice);
	cudaMemcpy(csrRowPtrB, CSRB.offset, sizeof(int) * (k + 1), cudaMemcpyHostToDevice);

	cudaMemcpy(csrColIdxA, CSRA.e, sizeof(int) * (nnzA), cudaMemcpyHostToDevice);
	cudaMemcpy(csrColIdxB, CSRB.e, sizeof(int) * (nnzB), cudaMemcpyHostToDevice);

	double alpha = 1.0;
	double beta = 0.0;

	size_t bufferSize;
	void* buffer = NULL;

	cusparseDcsrgemm2_bufferSizeExt(handle, n, m, k, &alpha,
		descrA, nnzA, csrRowPtrA, csrColIdxA,
		descrB, nnzB, csrRowPtrB, csrColIdxB,
		NULL,
		descrD, 0, csrRowPtrD, csrColIdxD,
		info,
		&bufferSize
	);
	cudaMalloc(&buffer, bufferSize);
	cout << "bufferSize=" << bufferSize << "!\n";

	cusparseXcsrgemm2Nnz(handle, n, m, k,
		descrA, nnzA, csrRowPtrA, csrColIdxA,
		descrB, nnzB, csrRowPtrB, csrColIdxB,
		descrD, 0, csrRowPtrD, csrColIdxD,
		descrC, csrRowPtrC, nnzTotalDevHostPtr,
		info, buffer
	);
	if (nnzTotalDevHostPtr != nullptr) {
		nnzC = *nnzTotalDevHostPtr;
	}
	else {
		cudaMemcpy(&nnzC, csrRowPtrC + n, sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(&baseC, csrRowPtrC, sizeof(int), cudaMemcpyDeviceToHost);
		nnzC -= baseC;
	}


	double* csrValA;
	double* csrValB;
	double* csrValC;
	double* csrValD;
	cudaMalloc(&csrColIdxC, sizeof(double) * nnzC);

	cudaMalloc(&csrValA, sizeof(double) * nnzA);
	cudaMalloc(&csrValB, sizeof(double) * nnzB);
	cudaMalloc(&csrValC, sizeof(double) * nnzC);
	cudaMalloc(&csrValD, sizeof(double) * 1);

	cudaMemcpy(csrValA, CSRA.w, sizeof(double) * nnzA, cudaMemcpyHostToDevice);
	cudaMemcpy(csrValB, CSRB.w, sizeof(double) * nnzB, cudaMemcpyHostToDevice);


	cusparseDcsrgemm2(handle, n, m, k, &alpha,
		descrA, nnzA, csrValA, csrRowPtrA, csrColIdxA,
		descrB, nnzB, csrValB, csrRowPtrB, csrColIdxB,
		NULL,
		descrD, 0, csrValD, csrRowPtrD, csrColIdxD,
		descrC, csrValC, csrRowPtrC, csrColIdxC,
		info, buffer
	);
	/*
	int *nzId;
	int *offset;
	int *e;
	double *w;
	int nnz;
	int r,c;

	*/
	ret->r = n;
	ret->c = m;
	ret->nnz = nnzC;
	ret->offset = (int*)malloc(sizeof(int) * (n + 1));
	ret->e = (int*)malloc(sizeof(int) * (nnzC));
	ret->w = (double*)malloc(sizeof(double) * (nnzC));

	cudaMemcpy(ret->offset, csrRowPtrC, sizeof(int) * (n + 1), cudaMemcpyDeviceToHost);
	cudaMemcpy(ret->e, csrColIdxC, sizeof(int) * nnzC, cudaMemcpyDeviceToHost);
	cudaMemcpy(ret->w, csrValC, sizeof(double) * nnzC, cudaMemcpyDeviceToHost);

	cout << nnzC << endl;
	for (int i = 0; i <= n; i++) {
		cout << ret->offset[i] << ' ';
	}cout << endl;
	for (int i = 0; i < nnzC; i++) {
		cout << ret->e[i] << ' ';
	}cout << endl;
	for (int i = 0; i < nnzC; i++) {
		cout << std::fixed << std::setprecision(3) << ret->w[i] << ' ';
	}
	cout << endl;


	cusparseDestroyCsrgemm2Info(info);
	cudaFree(csrRowPtrA);
	cudaFree(csrRowPtrB);
	cudaFree(csrRowPtrC);
	cudaFree(csrRowPtrD);
	cudaFree(csrColIdxA);
	cudaFree(csrColIdxB);
	cudaFree(csrColIdxC);
	cudaFree(csrColIdxD);
	cudaFree(csrValA);
	cudaFree(csrValB);
	cudaFree(csrValC);
	cudaFree(csrValD);
}

#undef CHECK_ERROR
