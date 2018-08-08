// Determination of deflection of vertical
// based on gravimetry data 
// and the already estimated W matrix
// Final type  with multiple sites
// 07 August 2018

#include <iostream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <quadmath.h>
//#include <omp.h>

using namespace std;

#define M2	2
#define M3	3
#define M6	6
#define STR_SIZE	50
typedef __float128 quadfloat;

quadfloat const zer = 0.000000000000000000000000;
quadfloat const one= 1.000000000000000000000000;


void matinv2(quadfloat a[M2][M2], quadfloat b[M2][M2])
{
//  calculates the inverse of a matrix using Gauss-Jordan elimination
//  the inverse of matrix a is calculated and stored in the matrix b
//  beware:  matrix a  is converted to unit matrix!

    int i,j,k;
    quadfloat dum;

//  Build the unit matrix b
    for (i=0;i<M2;i++)  {
      	for (j=0;j<M2;j++)
      		b[i][j]=zer;
      	b[i][i]=one;
    }

	for (k=0;k<M2;k++)  {

	  dum=a[k][k];
	  for (j=0;j<M2;j++)  {
	  	  a[k][j]=a[k][j]/dum;
	  	  b[k][j]=b[k][j]/dum;
	  }

	  for (i=k+1;i<M2;i++) {
	      dum=a[i][k];
	  	  for (j=0;j<M2;j++) {
	  	  	  a[i][j]=a[i][j]-a[k][j]*dum;
	  	  	  b[i][j]=b[i][j]-b[k][j]*dum;
	  	  }
	  }
	}

	for (k=M2-1;k>0;k--)
	  for (i=k-1;i>=0;i--) {
	      dum=a[i][k];
	  	  for (j=0;j<M2;j++) {
	  	  	  a[i][j]=a[i][j]-a[k][j]*dum;
	  	  	  b[i][j]=b[i][j]-b[k][j]*dum;
		 }
	  }
}

void mulmat22(quadfloat a[M2][M2],quadfloat b[M2][M2], quadfloat c[M2][M2])
{
	int i,j,k;
	quadfloat s;

	for (i=0;i<M2;i++)
	for (j=0;j<M2;j++) {
		s=zer;
		for (k=0;k<M2;k++)
			s=s+a[i][k]*b[k][j];
		c[i][j]=s;
	}
}

void mulmat2(quadfloat a[M2][M2],quadfloat b[M2], quadfloat c[M2])
{
	int i,j;
	quadfloat s;

	for (i=0;i<M2;i++) {
		s=zer;
		for (j=0;j<M2;j++)
			s=s+a[i][j]*b[j];
		c[i]=s;
	}
}

void matinv3(quadfloat a[M3][M3], quadfloat b[M3][M3])
{
//  calculates the inverse of a matrix using Gauss-Jordan elimination
//  the inverse of matrix a is calculated and stored in the matrix b
//  beware:  matrix a  is converted to unit matrix!

    int i,j,k;
    quadfloat dum;

//  Build the unit matrix b
    for (i=0;i<M3;i++)  {
      	for (j=0;j<M3;j++)
      		b[i][j]=zer;
      	b[i][i]=one;
    }

	for (k=0;k<M3;k++)  {

	  dum=a[k][k];
	  for (j=0;j<M3;j++)  {
	  	  a[k][j]=a[k][j]/dum;
	  	  b[k][j]=b[k][j]/dum;
	  }

	  for (i=k+1;i<M3;i++) {
	      dum=a[i][k];
	  	  for (j=0;j<M3;j++) {
	  	  	  a[i][j]=a[i][j]-a[k][j]*dum;
	  	  	  b[i][j]=b[i][j]-b[k][j]*dum;
	  	  }
	  }
	}

	for (k=M3-1;k>0;k--)
	  for (i=k-1;i>=0;i--) {
	      dum=a[i][k];
	  	  for (j=0;j<M3;j++) {
	  	  	  a[i][j]=a[i][j]-a[k][j]*dum;
	  	  	  b[i][j]=b[i][j]-b[k][j]*dum;
		 }
	  }
}

void mulmat33(quadfloat a[M3][M3],quadfloat b[M3][M3], quadfloat c[M3][M3])
{
	int i,j,k;
	quadfloat s;

	for (i=0;i<M3;i++)
	for (j=0;j<M3;j++) {
		s=zer;
		for (k=0;k<M3;k++)
			s=s+a[i][k]*b[k][j];
		c[i][j]=s;
	}
}

void mulmat3(quadfloat a[M3][M3],quadfloat b[M3], quadfloat c[M3])
{
	int i,j;
	quadfloat s;

	for (i=0;i<M3;i++) {
		s=zer;
		for (j=0;j<M3;j++)
			s=s+a[i][j]*b[j];
		c[i]=s;
	}
}


int main()
{
	int i,j,s,sid,sgn;

	quadfloat const omega=0.7292115e-04;
	quadfloat const pi=4.*atanq(1.);
	quadfloat const mf=100000.;        // Conversion from (m/sec^2) to mgal
	quadfloat const eot=1.0e09;  	     // Conversion from (sec^-2) to Eotvos units
	quadfloat const ras=648000./pi; 	 // Conversion fram rad to arcsec

	quadfloat lat,lon,h,gH;
	quadfloat SW[3][3],W[3][3],Win[3][3],Isu[3][3];
	quadfloat U[3][3],T[3][3],TW[3][3],UG[3][3],TG[3][3],WG[3][3];
	quadfloat Q[3][3],dmin[3][3],dmax[3][3],dcW[3][3];
	quadfloat dP[4][3],g[4],gp[4],hN[3],gN[3],wz[3];
	quadfloat gA,gB,gC,gS,gAp,gBp,gC1,dgS,gamms,om2;
	quadfloat ksi,eta,modksi,modeta,difksi,difeta,eps;
	quadfloat Wx,Wy,Wz,Ad,Bd,dum,st,ct,nWx,nWy,nWz;
	quadfloat TM[3],Rh[3],nR[2],nT[2],nW[2][2],nWin[2][2],nTW[2][2],nIsu[2][2];
	quadfloat epsa,epsb,alpha,beta,bet1,bet2;	
	quadfloat Ls,Ns,difg,leps,gP,Uyy,C,V,gApmod,gBpmod;
	quadfloat aq,bq,cq,diq,qr1,qr2,gApp,gBpp;
	quadfloat calp,cbet,disc,zet1,zet2,theta;
	quadfloat det1,det2,rat,gCp,nom,den;
	quadfloat a1,a2,a3,b1,b2,b3,c1,c2,c3;

    char x_s[STR_SIZE], y_s[STR_SIZE], z_s[STR_SIZE], w_s[STR_SIZE];
	memset(x_s, 0, STR_SIZE * sizeof(char));
	memset(y_s, 0, STR_SIZE * sizeof(char));
	memset(z_s, 0, STR_SIZE * sizeof(char));
	memset(w_s, 0, STR_SIZE * sizeof(char));
		
	ofstream res("Def-est-q-detail.txt");
	ofstream sta("Def-est-q-stats.txt");
	ifstream dat("Def_grav-q.txt");
	
	om2=2.0*omega*omega;

	for(i=0;i<3;i++) 
		for (j=0;j<3;j++)  {
			dmin[i][j]=100000.0;
			dmax[i][j]=-100000.0;
		}

for (s=0;s<12;s++) {

    dat>>sid;

	dat>>x_s>>y_s>>z_s>>w_s;
	lat   = strtoflt128(x_s, NULL);
	lon  = strtoflt128(y_s, NULL);
	gH  = strtoflt128(z_s, NULL);
	h    = strtoflt128(w_s, NULL);
    
	for (j=0;j<3;j++)
		dP[0][j]=zer;
		
	dat>>x_s>>y_s>>z_s;
	g[0]   = strtoflt128(x_s, NULL);
	dgS  = strtoflt128(y_s, NULL);
	gamms = strtoflt128(z_s, NULL);

 	dat>>w_s>>z_s;
	C  = strtoflt128(w_s, NULL);
	V  = strtoflt128(z_s, NULL);
   
//  Data in geometric system
	for (i=1;i<4;i++) {
		for (j=0;j<3;j++)  {
            dat>>x_s;
			dP[i][j] = strtoflt128(x_s, NULL);
        }
        dat>>y_s;
		g[i] =  strtoflt128(y_s, NULL);;
	}

 	dat>>w_s>>z_s;
	gApmod  = strtoflt128(w_s, NULL);
	gBpmod  = strtoflt128(z_s, NULL);

	gS=g[0]/mf;
	gA=g[1]/mf;
	gB=g[2]/mf;
	gC=g[3]/mf;
	dgS=dgS/mf;
	gamms=gamms/mf;
	gApmod=gApmod/mf;
	gBpmod=gBpmod/mf;
	leps=dgS/gS;

//  Data in physical system
	for (i=1;i<4;i++) {
		for (j=0;j<3;j++)  {
            dat>>x_s;
			dP[i][j] = strtoflt128(x_s, NULL);
        }
        dat>>y_s;
		gp[i] = strtoflt128(y_s, NULL);
	}

 	dat>>w_s>>z_s;
	gApp  = strtoflt128(w_s, NULL);
	gBpp  = strtoflt128(z_s, NULL);

    gApp=gApp/mf;
	gBpp=gBpp/mf;

	for (i=0;i<3;i++)
		for (j=0;j<3;j++)  {
            dat>>x_s;
			T[i][j] = strtoflt128(x_s, NULL);
        }

    for (i=0;i<3;i++)
		for (j=0;j<3;j++)  {
            dat>>y_s;
			U[i][j] = strtoflt128(y_s, NULL);
        }
	
    for (i=0;i<3;i++)
		for (j=0;j<3;j++)  {
            dat>>z_s;
			TG[i][j] = strtoflt128(z_s, NULL);
        }
	
    for (i=0;i<3;i++)
		for (j=0;j<3;j++)  {
            dat>>w_s;
			UG[i][j] = strtoflt128(w_s, NULL);
        }

	dat>>x_s>>y_s>>z_s;
	dum  = strtoflt128(x_s, NULL);
	modksi  = strtoflt128(y_s, NULL);
	modeta = strtoflt128(z_s, NULL);

	Uyy=dum/eot;

//  end of data

	for (i=0;i<3;i++)
		for (j=0;j<3;j++)  {
			T[i][j]=T[i][j]/eot;
			U[i][j]=U[i][j]/eot;
			UG[i][j]=UG[i][j]/eot;
			TG[i][j]=TG[i][j]/eot;
			SW[i][j]=T[i][j]+U[i][j];
			WG[i][j]=TG[i][j]+UG[i][j];
		}

	res<<fixed<<endl;
	quadmath_snprintf(x_s, sizeof(x_s), "%.6Qe", lat);
	quadmath_snprintf(y_s, sizeof(y_s), "%.6Qe", lon);
	quadmath_snprintf(z_s, sizeof(z_s), "%.6Qe", h);
	quadmath_snprintf(w_s, sizeof(w_s), "%.6Qe", gH);

	res<<"\n Point "<<sid;
    res<<" (S) :  Latitude = "<<setw(14)<<x_s;
	res<<" deg  -  Longitude = "<<setw(14)<<y_s<<" deg \n";
	res<<"  - Orthometric Height = "<<setw(14)<<z_s;
	res<<"  -   Geometric Height = "<<setw(14)<<w_s<<" m \n";

	W[2][2]=-(gC-gS)/dP[3][2];
	W[1][1]=Uyy*gS/gamms;
	W[0][0]=om2-W[1][1]-W[2][2];

	Ls=-W[0][0]/gS;
	Ns=-W[1][1]/gS;

	hN[1]=0.5*Ls*dP[1][0]*dP[1][0];
	epsa=fabsq(0.5*Ls*dP[1][0]);
	gN[1]=gA+W[2][2]*dP[1][2];
	wz[1]=gN[1]*cosq(epsa);
	W[0][2]=-(wz[1]-gS)/dP[1][0];

	hN[2]=0.5*Ns*dP[2][1]*dP[2][1];
	epsb=fabsq(0.5*Ns*dP[2][1]);
	gN[2]=gB+W[2][2]*dP[2][2];
	wz[2]=gN[2]*cosq(epsb);
	W[1][2]=-(wz[2]-gS)/dP[2][1];
	
    calp=W[0][0]-2.0*W[1][1];
    cbet=W[0][0];
    disc=sqrtq(calp*calp+cbet*cbet);

    zet1=(-cbet+disc)/calp;
    zet2=(-cbet-disc)/calp;

	alpha=atanq(zet1);
	bet1=pi/4.-alpha;
	theta=bet1;

	alpha=atanq(zet2);
	bet2=pi/4.-alpha;

	disc=zet1;
	if (fabsq(bet2)<fabsq(bet1))  {  
		disc=zet2;
		theta=bet2;
	}

	beta=zet1*zet2;
	
	W[0][1]=(disc*disc-one)*(W[1][1]-W[0][0])/(4.0*disc);

	nom=W[1][2];
	den=W[0][2];
	if (den!=0.)  
		theta=atanq(nom/den);
	else theta=pi/2.;
	if (den<0.)  theta=theta+pi;
	
	beta=W[1][1];
	alpha=W[0][0];
	dum=2.*theta;
	st=sinq(theta)*sinq(theta);
	ct=cosq(theta)*cosq(theta);
	
	W[0][0]=alpha*ct+beta*st;
	W[1][1]=alpha*st+beta*ct;
	W[0][1]=0.5*(alpha-beta)*sinq(dum);  

//  Iteration

	Ls=-W[0][0]/gS;
	Ns=-W[1][1]/gS;

	hN[1]=0.5*Ls*dP[1][0]*dP[1][0];
	epsa=fabsq(0.5*Ls*dP[1][0]);
	gN[1]=gA+W[2][2]*dP[1][2];
	wz[1]=gN[1]*cosq(epsa);
	W[0][2]=-(wz[1]-gS)/dP[1][0];

	hN[2]=0.5*Ns*dP[2][1]*dP[2][1];
	epsb=fabsq(0.5*Ns*dP[2][1]);
	gN[2]=gB+W[2][2]*dP[2][2];
	wz[2]=gN[2]*cosq(epsb);
	W[1][2]=-(wz[2]-gS)/dP[2][1];
	
    calp=W[0][0]-2.0*W[1][1];
    cbet=W[0][0];
    disc=sqrtq(calp*calp+cbet*cbet);

    zet1=(-cbet+disc)/calp;
    zet2=(-cbet-disc)/calp;

	alpha=atanq(zet1);
	bet1=pi/4.-alpha;
	theta=bet1;

	alpha=atanq(zet2);
	bet2=pi/4.-alpha;

	disc=zet1;
	if (fabsq(bet2)<fabsq(bet1))  {  
		disc=zet2;
		theta=bet2;
	}

	beta=zet1*zet2;
	
	W[0][1]=(disc*disc-one)*(W[1][1]-W[0][0])/(4.0*disc);

/*
	res<<"\n  Replace Wxz & Wyz with model values \n";
	W[0][2]=SW[0][2];
	W[1][2]=SW[1][2];
*/	
		
	for(i=0;i<3;i++) 
		for (j=i;j<3;j++)  
			W[j][i]=W[i][j];

 	res<<"\n  Computed Eotvos matrix W at point P (Eotvos units) \n\n";
	for(i=0;i<3;i++) {
		for (j=0;j<3;j++) {
			TW[i][j]=W[i][j];
            quadmath_snprintf(x_s, sizeof(x_s), "%.12Qe", TW[i][j]*eot);
			res<<setw(24)<<x_s;
		}
		res<<endl;
	}

	difg=(om2-W[0][0]-W[1][1]-W[2][2])*eot;
    quadmath_snprintf(z_s, sizeof(z_s), "%.16Qe", difg);
	res<<"\n 2ù^2 -ÓWii = "<<z_s<<"  Eotvos "<<endl;    
	
//-----------------------------------------

//	gN[1]=(W[0][0]*a1+W[0][1]*a2+W[0][2]*a3)*dP[1][0]/mf; 
//	gAp=gS+gN[1]/gS;
//	gN[2]=(W[1][0]*a1+W[1][1]*a2+W[1][2]*a3)*dP[2][1]/mf; 
//	gBp=gS+gN[2]/gS;

//	gAp=gA+W[2][2]*dP[1][2];
//	gBp=gB+W[2][2]*dP[2][2];

	gAp=gApmod;
	gBp=gBpmod;
	
	res<<"\n Results using computed W and prime geometric gravities\n";
    
	quadmath_snprintf(x_s, sizeof(x_s), "%.16Qe", gAp*mf);
	quadmath_snprintf(y_s, sizeof(y_s), "%.16Qe", gBp*mf);
	res<<"\n  gAp = "<<x_s<<"  mgal ";
	res<<"\n  gBp = "<<y_s<<"  mgal \n";

	Rh[0]=(gAp-gS)*gS/dP[1][0];
	Rh[1]=(gBp-gS)*gS/dP[2][1];
	Rh[2]=(gC-gS)*gS/dP[3][2];

	res<<"\n Right-hand matrix of system at P\n";
    
	quadmath_snprintf(x_s, sizeof(x_s), "%.16Qe", Rh[0]);
	quadmath_snprintf(y_s, sizeof(y_s), "%.16Qe", Rh[1]);
	quadmath_snprintf(z_s, sizeof(z_s), "%.16Qe", Rh[2]);
	res<<"\n G1 = "<<x_s;
	res<<"\n G2 = "<<y_s;
	res<<"\n G3 = "<<z_s;
	
	matinv3(TW,Win);
	
	mulmat3(Win,Rh,TM);

	Wx=TM[0];
	Wy=TM[1];
	Wz=TM[2];

//  check
    Bd=Wx*Wx+Wy*Wy+Wz*Wz-gS*gS;
    quadmath_snprintf(z_s, sizeof(z_s), "%.10Qe", Bd);
    res<<"\n Wx^2 + Wy^2 + Wz^2 - gP^2 = "<<z_s<<"   (zero)\n";
	
	res<<"\n  First-order gradients";
	quadmath_snprintf(x_s, sizeof(x_s), "%.16Qe", Wx);
	quadmath_snprintf(y_s, sizeof(y_s), "%.16Qe", Wy);
	quadmath_snprintf(z_s, sizeof(z_s), "%.16Qe", Wz);    
	res<<"\n  Wx = "<<x_s;
	res<<"\n  Wy = "<<y_s;
	res<<"\n  Wz = "<<z_s<<endl;

	dum=sqrtq(Wx*Wx+Wy*Wy+Wz*Wz);
	gP=dum;
	dum=fabsq(Wz)/dum;
	eps=acosq(dum);
	eps=eps*ras;

	ksi=Wx/Wz;
	ksi=ksi*ras;

	eta=Wy/Wz;
	eta=eta*ras;

	difksi=ksi-modksi;
	difeta=eta-modeta;

	res<<"\n  Estimated Deflection of vertical - geometric \n";
	quadmath_snprintf(x_s, sizeof(x_s), "%.8Qe", eps);
	quadmath_snprintf(y_s, sizeof(y_s), "%.8Qe", ksi);
	quadmath_snprintf(z_s, sizeof(z_s), "%.8Qe", difksi);    
	res<<"\n eps = "<<x_s<<" arcsec \n";
	res<<"\n ksi = "<<y_s<<" arcsec   -  Difference = "<<z_s<<" arcsec";
	quadmath_snprintf(y_s, sizeof(y_s), "%.8Qe", eta);
	quadmath_snprintf(z_s, sizeof(z_s), "%.8Qe", difeta);    
	res<<"\n eta = "<<y_s<<" arcsec   -  Difference = "<<z_s<<" arcsec\n";

	for (j=0;j<3;j++) {
		beta=difksi;
		if (beta<dmin[0][0])  dmin[0][0]=beta;
		if (beta>dmax[0][0])  dmax[0][0]=beta;
		beta=difeta;
		if (beta<dmin[1][1])  dmin[1][1]=beta;
		if (beta>dmax[1][1])  dmax[1][1]=beta;
	}

	//*
	res<<"\n  Computed Disturbance matrix T at point S (in Eotvos units)\n\n";
	for(i=0;i<3;i++) {
		for (j=0;j<3;j++) {
			Q[i][j]=(W[i][j]-U[i][j]);
            quadmath_snprintf(w_s, sizeof(w_s), "%.12Qe", Q[i][j]*eot);
			res<<setw(24)<<w_s;
		}
		res<<endl;
	}
	
	res<<"\n\n Difference between computed-simulated matrix T at point S (in Eotvos units)\n\n";
	for(i=0;i<3;i++) {
		for (j=0;j<3;j++) {
			dcW[i][j]=(Q[i][j]-T[i][j])*eot;
            quadmath_snprintf(x_s, sizeof(x_s), "%.12Qe", dcW[i][j]*eot);
			res<<setw(24)<<x_s;
			beta=(dcW[i][j]);
			if (beta<dmin[i][j])  dmin[i][j]=beta;
			if (beta>dmax[i][j])  dmax[i][j]=beta;
		}
		res<<endl;
	}

//*/

	res<<"\n\n Results using model W \n";

	res<<"\n  Model Eotvos matrix W at point P (Eotvos units) \n";
	for(i=0;i<3;i++) {
		for (j=0;j<3;j++)  {
            quadmath_snprintf(y_s, sizeof(y_s), "%.16Qe", SW[i][j]*eot);
			res<<setw(26)<<y_s;
        }
		res<<endl;
	}
	
	gAp=gApp;
	gBp=gBpp;
	
	res<<"\n\n Results using model W and prime physical gravities\n";
	quadmath_snprintf(x_s, sizeof(x_s), "%.16Qe", gAp*mf);
	quadmath_snprintf(y_s, sizeof(y_s), "%.16Qe", gBp*mf);
	res<<"\n  gAp = "<<x_s<<"  mgal ";
	res<<"\n  gBp = "<<y_s<<"  mgal \n";
    
	Rh[0]=(gAp-gS)*gS/dP[1][0];
	Rh[1]=(gBp-gS)*gS/dP[2][1];
	Rh[2]=(gC-gS)*gS/dP[3][2];

	res<<"\n Model Right-hand matrix of system at P\n";
    
	quadmath_snprintf(x_s, sizeof(x_s), "%.16Qe", Rh[0]);
	quadmath_snprintf(y_s, sizeof(y_s), "%.16Qe", Rh[1]);
	quadmath_snprintf(z_s, sizeof(z_s), "%.16Qe", Rh[2]);
	res<<"\n G1 = "<<x_s;
	res<<"\n G2 = "<<y_s;
	res<<"\n G3 = "<<z_s;

	matinv3(SW,Win);

	mulmat3(Win,Rh,TM);

	Wx=TM[0];
	Wy=TM[1];
	Wz=TM[2];

	res<<endl;
	res<<"\n  First-order gradients";
	quadmath_snprintf(x_s, sizeof(x_s), "%.16Qe", Wx);
	quadmath_snprintf(y_s, sizeof(y_s), "%.16Qe", Wy);
	quadmath_snprintf(z_s, sizeof(z_s), "%.16Qe", Wz);    
	res<<"\n  Wx = "<<x_s;
	res<<"\n  Wy = "<<y_s;
	res<<"\n  Wz = "<<z_s<<endl;
    
//  check
    Bd=Wx*Wx+Wy*Wy+Wz*Wz-gS*gS;
    quadmath_snprintf(z_s, sizeof(z_s), "%.10Qe", Bd);
    res<<"\n Wx^2 + Wy^2 + Wz^2 - gP^2 = "<<z_s<<"   (zero)\n";
    
	dum=sqrtq(Wx*Wx+Wy*Wy+Wz*Wz);
	dum=fabsq(Wz)/dum;
	eps=acosq(dum);
	eps=eps*ras;

	ksi=Wx/Wz;
	ksi=ksi*ras;

	eta=Wy/Wz;
	eta=eta*ras;

	difksi=ksi-modksi;
	difeta=eta-modeta;

	res<<"\n  Deflection of vertical - physical \n";
	quadmath_snprintf(x_s, sizeof(x_s), "%.8Qe", eps);
	quadmath_snprintf(y_s, sizeof(y_s), "%.8Qe", ksi);
	quadmath_snprintf(z_s, sizeof(z_s), "%.8Qe", difksi);    
	res<<"\n eps = "<<x_s<<" arcsec \n";
	res<<"\n ksi = "<<y_s<<" arcsec   -  Difference = "<<z_s<<" arcsec";
	quadmath_snprintf(y_s, sizeof(y_s), "%.8Qe", eta);
	quadmath_snprintf(z_s, sizeof(z_s), "%.8Qe", difeta);    
	res<<"\n eta = "<<y_s<<" arcsec   -  Difference = "<<z_s<<" arcsec\n";
    
}

	sta<<fixed;
	sta<<"\n Range of differences between computed & model vertical deviations - geometric (in arcsec)\n";
	quadmath_snprintf(x_s, sizeof(x_s), "%.8Qe", dmin[0][0]);
	quadmath_snprintf(y_s, sizeof(y_s), "%.8Qe", dmax[0][0]);
	quadmath_snprintf(z_s, sizeof(z_s), "%.8Qe", dmin[1][1]);    
	quadmath_snprintf(w_s, sizeof(w_s), "%.8Qe", dmax[1][1]);    
    
	sta<<"\n ksi : "<<setw(18)<<x_s<<"  to "<<setw(19)<<y_s;
	sta<<"\n eta : "<<setw(18)<<z_s<<"  to "<<setw(19)<<w_s<<endl;
	
	return 0;
}
