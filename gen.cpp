#include<bits/stdc++.h>
#define rep(i,a,b) for(int i=(a),i##ss=(b);i<=i##ss;i++)
#define dwn(i,a,b) for(int i=(a),i##ss=(b);i>=i##ss;i--)
#define rng(i,a,b) for(int i=(a);i<(b);i++)
#define deb(x) cerr<<(#x)<<":"<<(x)<<'\n'
#define pb push_back
#define mkp make_pair
#define fi first
#define se second
#define hvie '\n'
using namespace std;
typedef pair<int,int> pii;
typedef long long ll;
typedef unsigned long long ull;
typedef double db;
int yh(){
	int ret=0;bool f=0;char c=getchar();
	while(!isdigit(c)){if(c==EOF)return -1;if(c=='-')f=1;c=getchar();}
	while(isdigit(c))ret=(ret<<3)+(ret<<1)+(c^48),c=getchar();
	return f?-ret:ret;
}
const int maxn=3e5+5;
ull randll(){
	static ull seed=1141514;
	seed ^= seed << 13;
	seed ^= seed >> 7;
	seed ^= seed << 17;
	return seed;
}
double randreal(){
	return (long long)randll()*1.0/LLONG_MAX;
}


signed main(){
	int n=yh(),m=yh(),nnz=yh();
	assert(0.9*n*m>=nnz);

	string name;
	cin>>name;

	freopen(name.c_str(),"w", stdout);
	cout<<n<<' '<<m<<' '<<nnz<<'\n';

	set<pair<int,int>>pos;

	for(int i=1;i<=nnz;i++){
		int x,y;
		do{
			x=randll()%n;
			y=randll()%m;
		}while(pos.count({x,y}));
		pos.insert({x,y});
		cout<<x<<' '<<y<<' '<<(x+y)%7+1<<"\n";
		// cout<<x<<' '<<y<<' '<<randreal()<<'\n';
	}

	return 0;
}
