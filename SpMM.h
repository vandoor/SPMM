#pragma once
#include <iostream>
using std::cout;
using std::endl;

struct tpl {
	int r, c;
	double val;
};
struct SparseMatrix {
	int n, m, nnz;
	tpl* coo;
	SparseMatrix(int n, int m, int nnz) :n(n), m(m), nnz(nnz) {
		coo = (tpl*)malloc(sizeof(tpl) * nnz);
		for (int i = 0; i < nnz; i++) coo[i] = { i / m,i % m,1 };
	}
	SparseMatrix(char* s) {
		int r, c;
		double va;
		std::ifstream fin;
		fin.open(s, std::ios::in);
		fin >> n >> m >> nnz;
		coo = (tpl*)malloc(sizeof(tpl) * nnz);
		for (int i = 0; i < nnz; i++) {
			fin >> r >> c >> va;
			coo[i] = { r,c,va };
		}
		fin.close();
	}
	void print() {
		std::vector<std::vector<double>>a(n, std::vector<double>(m));
		for (int i = 0; i < nnz; i++) {
			a[coo[i].r][coo[i].c] = coo[i].val;
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				cout << std::fixed << std::setprecision(3) << a[i][j] << ' ';
			}cout << endl;
		}
	}
	void del() {
		cout << "!";
		free(coo);
		cout << "ok\n";
	}
};

struct DCSR {
	int* nzId;
	int* offset;
	int* e;
	double* w;
	int nnz;
	int r, c;
	DCSR() {
		nzId=offset=e=nullptr;
		w = nullptr;
		nnz=r=c=0;
	}
	DCSR(SparseMatrix m, bool COL = 0, bool compress = 0) {
		int R = COL ? m.m : m.n;
		int C = COL ? m.n : m.m;
		r = R, c = C;
		int U = m.nnz;
		int r, c;
		double v;
		std::vector<std::vector<std::pair<int, double> > >rec(R);
		if (compress) {
			std::vector<int>ids;
			for (int i = 0; i < U; i++) {
				tpl t = m.coo[i];
				r = t.r;
				c = t.c;
				if (COL) std::swap(r, c);
				v = t.val;
				rec[r].push_back({ c,v });
				if ((int)rec[r].size() == 1) {
					ids.push_back(r);
				}
			}
			int K = ids.size();
			nzId = (int*)malloc(K * sizeof(int));
			offset = (int*)malloc((K + 1) * sizeof(int));
			e = (int*)malloc(U * sizeof(int));
			w = (double*)malloc(U * sizeof(double));
			for (int i = 0; i < K; i++) nzId[i] = ids[i];
			int cur = 0;
			for (int i = 0; i < K; i++) {
				offset[i] = cur;
				sort(rec[nzId[i]].begin(), rec[nzId[i]].end());
				for (auto& [c, v] : rec[nzId[i]]) {
					e[cur] = c;
					w[cur] = v;
					cur++;
				}
			}
			offset[K] = cur;
			assert(cur == U);
		}
		else {
			std::vector<int>ids;
			for (int i = 0; i < U; i++) {
				tpl t = m.coo[i];
				r = t.r;
				c = t.c;
				if (COL) std::swap(r, c);
				v = t.val;
				rec[r].push_back({ c,v });
			}
			offset = (int*)malloc((R + 1) * sizeof(int));
			e = (int*)malloc(U * sizeof(int));
			w = (double*)malloc(U * sizeof(double));
			int cur = 0;
			for (int i = 0; i < R; i++) {
				offset[i] = cur;
				sort(rec[i].begin(), rec[i].end());
				for (auto& [c, v] : rec[i]) {
					e[cur] = c;
					w[cur] = v;
					cur++;
				}
			}
			offset[R] = cur;
			assert(cur == U);
		}
		nnz = U;
	}
	void print() {
		cout << r << ' ' << c << ":" << nnz << endl;
		for (int i = 0; i <= r; i++) {
			cout << offset[i] << ' ';
		}cout << endl;
		for (int i = 0; i < nnz; i++) {
			cout << e[i] << ' ';
		}cout << endl;
		for (int i = 0; i < nnz; i++) {
			cout << w[i] << ' ';
		}cout << endl;
	}
	void del() {
		if (nzId != nullptr) free(nzId);
		free(offset);
		free(e);
		free(w);
	}
	void setspace(int R,int nnz){
		if(offset!=nullptr) this->del();
		offset=new int [R];
		e=new int[nnz];
		w=new double[nnz];
	}
};

