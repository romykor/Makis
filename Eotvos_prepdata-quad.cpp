// Preparation of Eotvos matrix
// based on EGM2008 data
// Version Quad
// 05 August 2018

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <quadmath.h>
//#include <omp.h>

using namespace std;

//#define MATR_DIM	3
#define STR_SIZE	50
typedef __float128 quadfloat;


int main()
{
	int i,j,s,yp,id,n,max;

	quadfloat p[3][3],T[3][3];
	quadfloat phi,lam,dg,tem,geoh,ortho;
	quadfloat pl,pr,zl,zr,dum,ksi,eta;
	quadfloat const zer=0.000000000000000000000000;
	
//  Parameters
	pl=2.00000000;
	pr=2.00000000;
	zl=0.50000000;
	zr=0.50000000;
	max=RAND_MAX;  
	cout<<"\n RAND_MAX = "<<max<<endl;
	
	ifstream dat("mod_data.txt");
	if (!dat.is_open())  return 1;
	
	ofstream res("Eotvos_simuldata-quad.txt");
	if (!res.is_open())  return 2;
	
	res<<fixed;

	char x_s[STR_SIZE], y_s[STR_SIZE], z_s[STR_SIZE], w_s[STR_SIZE];
	memset(x_s, 0, STR_SIZE * sizeof(char));
	memset(y_s, 0, STR_SIZE * sizeof(char));
	memset(z_s, 0, STR_SIZE * sizeof(char));
	memset(w_s, 0, STR_SIZE * sizeof(char));
	
	srand(time(0));
	
  for (n=1;n<=12;n++)  {
		
	for (i=0;i<3;i++) 
		for (j=0;j<3;j++)  {
			T[i][j]=zer;
			
		}
		
	dat>>id>>x_s>>y_s>>z_s>>w_s;
	phi   = strtoflt128(x_s, NULL);
	lam   = strtoflt128(y_s, NULL);
	geoh  = strtoflt128(z_s, NULL);
	ortho = strtoflt128(w_s, NULL);
		
	dat>>x_s>>y_s>>z_s>>w_s;
	dg    = strtoflt128(x_s, NULL);
	dum   = strtoflt128(y_s, NULL);
	ksi   = strtoflt128(z_s, NULL);
	eta   = strtoflt128(w_s, NULL);
		
	for (i=0;i<3;i++) 
		for (j=i;j<3;j++)  {
			dat>>x_s;
			T[i][j] = strtoflt128(x_s, NULL);
		}
		
/*
//  reverse x - y  axes
    tem=T[0][0];
	T[0][0]=T[1][1];
	T[1][1]=tem;
    tem=T[0][2];
	T[0][2]=T[1][2];
	T[1][2]=tem;
//
*/

	for (i=1;i<3;i++) 
		for (j=0;j<i;j++)  
			T[i][j]=T[j][i];

	res<<endl<<setw(3)<<id<<endl;
	
	quadmath_snprintf(x_s, sizeof(x_s), "%.6Qe", phi);
	quadmath_snprintf(y_s, sizeof(y_s), "%.6Qe", lam);
	quadmath_snprintf(z_s, sizeof(z_s), "%.6Qe", geoh);
	quadmath_snprintf(w_s, sizeof(w_s), "%.6Qe", ortho);

	res<<setw(25)<<x_s<<endl;
	res<<setw(25)<<y_s<<endl;
	res<<setw(25)<<z_s<<endl;
	res<<setw(25)<<w_s<<endl;

	quadmath_snprintf(x_s, sizeof(x_s), "%.8Qe", dg);
	res<<endl<<setw(25)<<x_s<<endl;

	s=1;
	yp=rand()%2;
	if (yp==1)  s=-1;
	p[0][0]=s*(pl+pr*rand()/max);
	p[0][2]=zl+zr*rand()/max;

	s=1;
	yp=rand()%2;
	if (yp==1)  s=-1;
	p[1][1]=s*(pl+pr*rand()/max);
	p[1][2]=zl+zr*rand()/max;
	
	p[2][2]=1.+zl+zr*rand()/max;

	for (i=0;i<3;i++) {
		res<<endl;
		for (j=0;j<3;j++) {
			quadmath_snprintf(x_s, sizeof(x_s), "%.6Qe", p[i][j]);
			res<<setw(25)<<x_s;
		}
	}
	res<<endl;

	for (i=0;i<3;i++) {
		res<<endl;
		for (j=0;j<3;j++) {
			quadmath_snprintf(y_s, sizeof(y_s), "%.6Qe", T[i][j]);
			res<<setw(25)<<y_s;
		}
	}
	res<<endl;

	quadmath_snprintf(z_s, sizeof(z_s), "%.6Qe", ksi);
	quadmath_snprintf(w_s, sizeof(w_s), "%.6Qe", eta);

	res<<endl<<setw(19)<<z_s;
	res<<endl<<setw(19)<<w_s;
	res<<endl;
}
	dat.close();
	res.close();
	
	return 0;
}
