#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<array>
#include<vector>
#include<assert.h>
#include<fstream>
#include<iomanip>

double randreal(){return 0;}

using namespace std;

struct tpl{
	int r,c;
	double val;
};
struct SparseMatrix{
	int n,m,nnz;
	tpl *coo;
	SparseMatrix(int n,int m,int nnz):n(n),m(m),nnz(nnz){
		coo=(tpl*)malloc(sizeof(tpl)*nnz);
		for(int i=0;i<nnz;i++) coo[i]={i/m,i%m,randreal()};
	}
	SparseMatrix(char*s){
		int r,c; 
		double va;
		std::ifstream fin;
		fin.open(s, std::ios::in);
		fin>>n>>m>>nnz;
		coo=(tpl*)malloc(sizeof(tpl)*nnz);
		for(int i=0;i<nnz;i++){
			fin>>r>>c>>va;
			coo[i]={r,c,va};
		}
		fin.close();
	}
	void print(){
		std::vector<std::vector<double>>a(n,std::vector<double>(m));
		for(int i=0;i<nnz;i++){
			a[coo[i].r][coo[i].c]=coo[i].val;
		}
		for(int i=0;i<n;i++){
			for(int j=0;j<m;j++){
				cout<<std::fixed<<std::setprecision(3)<<a[i][j]<<' ';
			}cout<<endl;
		}
	}
	void del(){
		cout<<"!";
		free(coo);
		cout<<"ok\n";
	}
};

struct DCSR{
	int *nzId;
	int *offset;
	int *e;
	double *w;
	int nnz;
	int r,c;
	DCSR(){}
	DCSR(SparseMatrix m, bool COL=0, bool compress=0){
		int R=COL?m.m:m.n;
		int C=COL?m.n:m.m;
		r=R, c=C;
		int U=m.nnz;
		int r,c;
		double v;
		std::vector<std::vector<std::pair<int,double> > >rec(R);
		if(compress){
			std::vector<int>ids;
			for(int i=0;i<U;i++){
				tpl t=m.coo[i];
				r=t.r;
				c=t.c;
				if(COL) std::swap(r,c);
				v=t.val;
				rec[r].push_back({c,v});
				if((int)rec[r].size()==1){
					ids.push_back(r);
				}
			}
			int K=ids.size();
			nzId=(int*)malloc(K*sizeof(int));
			offset=(int*)malloc((K+1)*sizeof(int));
			e=(int*)malloc(U*sizeof(int));
			w=(double*)malloc(U*sizeof(double));
			for(int i=0;i<K;i++) nzId[i]=ids[i];
			int cur=0;
			for(int i=0;i<K;i++){
				offset[i]=cur;
				for(auto&[c,v]:rec[nzId[i]]){
					e[cur]=c;
					w[cur]=v;
					cur++;
				}
			}
			offset[K]=cur;
			assert(cur==U);
		}
		else{
			std::vector<int>ids;
			for(int i=0;i<U;i++){
				tpl t=m.coo[i];
				r=t.r;
				c=t.c;
				if(COL) std::swap(r,c);
				v=t.val;
				rec[r].push_back({c,v});
			}
			offset=(int*)malloc((R+1)*sizeof(int));
			e=(int*)malloc(U*sizeof(int));
			w=(double*)malloc(U*sizeof(double));
			int cur=0;
			for(int i=0;i<R;i++){
				offset[i]=cur;
				for(auto&[c,v]:rec[i]){
					e[cur]=c;
					w[cur]=v;
					cur++;
				}
			}
			offset[R]=cur;
			assert(cur==U);
		}
		nnz=U;
	}
};

template<typename T>
struct KV{
	int k;
	T v;
};


using ull=unsigned long long;
int a[1000];
int find(int x){
	int l=0,r=999,mid;
	while(l<=r){
		mid=(l+r)>>1;
		if(a[mid]==x) {
			return mid;
		}
		else if(a[mid]>x){
			r=mid-1;
		}
		else{ 
			l=mid+1;
		}
	}
	return -1;
}
int main(){
	for(int i=0;i<1000;i++){
		a[i]=2*i;
	}
	int x;cin>>x;
	cout<<find(x);
	return 0;
}