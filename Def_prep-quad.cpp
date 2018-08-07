// Simulation of gravimetry
// based on EGM2008 data
// 07 August 2018
// Quad version with multiple sites, 3 net points
// Matrix W in physical system, projection to geoid
// (Rotation of U according to deflection of vertical)

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

#define M3	3
#define M6	6
#define STR_SIZE	50
typedef __float128 quadfloat;

quadfloat const zer = 0.000000000000000000000000;
quadfloat const one= 1.000000000000000000000000;

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

void mulmat6(quadfloat a[M6][M6],quadfloat b[M6], quadfloat c[M6])
{
	int i,j;
	quadfloat s;

	for (i=0;i<M6;i++) {
		s=zer;
		for (j=0;j<M6;j++)
			s=s+a[i][j]*b[j];
		c[i]=s;
	}
}

void matinv6(quadfloat a[M6][M6], quadfloat b[M6][M6])
{
//  calculates the inverse of a matrix using Gauss-Jordan elimination
//  the inverse of matrix a is calculated and stored in the matrix b
//  beware:  matrix a  is converted to unit matrix!

    int i,j,k;
    quadfloat dum;

//  Build the unit matrix b
    for (i=0;i<M6 ;i++)  {
      	for (j=0;j<M6 ;j++)
      		b[i][j]=zer;
      	b[i][i]=one;
    }

	for (k=0;k<M6 ;k++)  {

	  dum=a[k][k];
	  for (j=0;j<M6 ;j++)  {
	  	  a[k][j]=a[k][j]/dum;
	  	  b[k][j]=b[k][j]/dum;
	  }

	  for (i=k+1;i<M6 ;i++) {
	      dum=a[i][k];
	  	  for (j=0;j<M6 ;j++) {
	  	  	  a[i][j]=a[i][j]-a[k][j]*dum;
	  	  	  b[i][j]=b[i][j]-b[k][j]*dum;
	  	  }
	  }
	}

	for (k=M6-1;k>0;k--)
	  for (i=k-1;i>=0;i--) {
	      dum=a[i][k];
	  	  for (j=0;j<M6 ;j++) {
	  	  	  a[i][j]=a[i][j]-a[k][j]*dum;
	  	  	  b[i][j]=b[i][j]-b[k][j]*dum;
		 }
	  }
}


int main()
{
	int i,j,n,s,sid;

	quadfloat const omega=0.7292115e-04;
	quadfloat const pi=4.*atanq(1.);
	quadfloat const mf=100000.;        // Conversion from (m/sec^2) to mgal
	quadfloat const eot=1.0e09;  	     // Conversion from (sec^-2) to Eotvos units
	quadfloat const ras=648000./pi; 	 // Conversion fram rad to arcsec
	quadfloat const a=6378137.0;
	quadfloat const finv=298.257222101;
	quadfloat const b=6356752.3142;
	quadfloat const gamf=9.7803267714e00;
	quadfloat const gamp=9.8321863685e00;
	
	quadfloat lat,lon,h,dg,gamma,gS,gell;
	quadfloat rl,dum,difg,ksi,eta,geksi,geeta,gH;  
	quadfloat R,f,onf,e2,ep2,C,k,k1,k2,J,V;
	quadfloat pgdf,pk2df,pk1df,pJdx;
	quadfloat U[M3][M3],UG[M3][M3],T[M3][M3],TG[M3][M3],W[M3][M3];
	quadfloat dP[M3][M3],WG[M3][M3],w1[M3],g[5],gp[5];
	quadfloat Ux,Uz,kg,om2,Ad,Bd,Uyy;
	quadfloat eps,sk,ck,si,ci,se,ce;
	quadfloat f21,f22,f23,f41,f42,f43,f51,f52;
	quadfloat P[M6][M6],PW[M6][M6],PIN[M6][M6];
	quadfloat gu[M6],ru[M6];

    char x_s[STR_SIZE], y_s[STR_SIZE], z_s[STR_SIZE], w_s[STR_SIZE];
	memset(x_s, 0, STR_SIZE * sizeof(char));
	memset(y_s, 0, STR_SIZE * sizeof(char));
	memset(z_s, 0, STR_SIZE * sizeof(char));
	memset(w_s, 0, STR_SIZE * sizeof(char));
	

//  Parameters
	R=(2*a+b)/3;
	f=one/finv;
	onf=one-f;
	e2=f*(2.-f);
	ep2=e2/(onf*onf);
	C=a/onf;
	om2=2.*omega*omega;
	kg=(b*gamp)/(a*gamf)-one;

	ofstream res("Def_simul-q.txt");
	ifstream dat("Eotvos_simuldata-quad.txt"); 
	ofstream gnor("Def_grav-q.txt");

  for (s=0;s<12;s++)	 {

	for(i=0;i<3;i++) 
		for (j=0;j<3;j++)  {
			T[i][j]=zer;
			U[i][j]=zer;
			W[i][j]=zer;
			TG[i][j]=zer;
			UG[i][j]=zer;
			WG[i][j]=zer;
		}
			
    dat>>sid;
    
	dat>>x_s>>y_s>>z_s>>w_s;
	lat   = strtoflt128(x_s, NULL);
	lon  = strtoflt128(y_s, NULL);
	gH  = strtoflt128(z_s, NULL);
	h    = strtoflt128(w_s, NULL);
    
	dat>>x_s;
    dg = strtoflt128(x_s, NULL);
	dg =dg/mf;
	rl = lat*pi/180.;

	dum=cosq(rl);
	V=sqrtq(one+ep2*dum*dum);
	dum=sinq(rl);
	dum=dum*dum*e2;
	k2=sqrtq(one-dum)/a;
	k1=k2*(one-dum)/(one-e2);
	J=0.5*(k1+k2);
	
	dum=sinq(rl);
	Ad=one-e2*dum*dum;
	Bd=one+kg*dum*dum;
	gell=gamf*Bd/sqrtq(Ad);

	dum=sinq(2.*rl);
	pgdf=gell/Bd*(kg*dum+0.5*e2*dum*Bd/Ad);
	k=pgdf/(R*gell);

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
    res<<"\n Ellipsoid parameters \n";
    
	quadmath_snprintf(x_s, sizeof(x_s), "%.12Qe", k1);
	quadmath_snprintf(y_s, sizeof(y_s), "%.12Qe", k2);
	quadmath_snprintf(z_s, sizeof(z_s), "%.12Qe", J);
	quadmath_snprintf(w_s, sizeof(w_s), "%.12Qe", k);

    res<<"  k1 = "<<x_s<<" m^-1 \n";
	res<<"  k2 = "<<y_s<<" m^-1 \n";
	res<<"   J = "<<z_s<<" m^-1 \n";
	res<<"   k = "<<w_s<<" m^-1 \n";
 	
	pk1df=-1.5*e2*dum*k2/(1.-e2);
	pk2df=-0.5*e2*dum/(k2*a*a);
//	pJdx=(pk1df+pk2df)/R;  //  double value
	pJdx=0.5*(pk1df+pk2df)/R;  
	
	dum=-gell*k*(k2-k1)-gell*pJdx;
	dum=0.5*dum*gH*gH;
	Ux=dum-gell*k*gH; 

	dum=-gell*(k1*k1+k2*k2);  //  p2gdx  omitted
	dum=0.5*dum*gH*gH;
	Uz=(om2+2.*gell*J)*gH;
	Uz=Uz+dum-gell;

    gamma=sqrtq(Ux*Ux+Uz*Uz);
	
	gS=gamma+dg;
	quadmath_snprintf(x_s, sizeof(x_s), "%.16Qe", gell);
	quadmath_snprintf(y_s, sizeof(y_s), "%.16Qe", gamma);
	quadmath_snprintf(z_s, sizeof(z_s), "%.16Qe", dg);
	quadmath_snprintf(w_s, sizeof(w_s), "%.16Qe", gS);

    res<<"\n Foot on ellipsoid\n  gamma (ell) = "<<x_s;
	res<<"\n Point S  (on surface) \n  gamma = "<<y_s;
	res<<"\n dg = "<<z_s;
	res<<"\n  g = "<<w_s<<endl;
    
 	for(i=0;i<3;i++) 
		for (j=0;j<3;j++)  {
            dat>>x_s;
			dP[i][j] = strtoflt128(x_s, NULL);
        }
		
	for(i=0;i<3;i++) 
		for (j=0;j<3;j++)  {      //  on surface
			dat>>y_s;
            T[i][j] = strtoflt128(y_s, NULL);
		}
	
	dat>>z_s>>w_s;         //  on surface
	ksi = strtoflt128(z_s, NULL);
	eta = strtoflt128(w_s, NULL);
    res<<"\n ksi = "<<setw(15)<<z_s<<" arcsec";
    res<<"\n eta = "<<setw(15)<<w_s<<" arcsec\n";
    
	geksi=ksi/ras;
	geeta=eta/ras;
   
	res<<fixed;
	res<<"\n  Disturbance matrix  T at point S - physical system (in Eotvos units)\n\n";
	for(i=0;i<3;i++) {
		for (j=0;j<3;j++) {
            quadmath_snprintf(y_s, sizeof(y_s), "%.12Qe", T[i][j]);
			res<<setw(23)<<y_s;
			T[i][j]=T[i][j]/eot;
		}
		res<<endl;
	}
	
	U[0][0]=-gell*k1;
	U[0][2]=-gell*k;
	U[1][1]=-gell*k2;
	U[2][2]=2.*(omega*omega+gell*J);
	U[2][0]=U[0][2];

	res<<"\n  Constant part of Normal Eotvos matrix U at ellipsoid - geometric system (in  sec^-2)\n\n";
	for(i=0;i<3;i++) {
		for (j=0;j<3;j++) {
            quadmath_snprintf(z_s, sizeof(z_s), "%.16Qe", U[i][j]);
			res<<setw(28)<<z_s;
        }
		res<<endl;
	}

	U[0][0]=U[0][0]+gell*k1*k1*gH; 
	U[0][2]=U[0][2]-(gell*k*(k2-k1)+2*gell*pJdx)*gH;
	U[1][1]=U[1][1]+gell*k2*k2*gH;  //  p2gdx  omitted
	U[2][2]=U[2][2]-gell*(k1*k1+k2*k2)*gH;  //  p2gdx  omitted
//	U[2][0]=U[0][2];
		
	for(i=0;i<3;i++) 
		for (j=i;j<3;j++)  
			U[j][i]=U[i][j];

	res<<"\n  Complete Normal Eotvos matrix U at point S - geometric system (in sec^-2)\n\n";
	for(i=0;i<3;i++) {
		for (j=0;j<3;j++) {
            quadmath_snprintf(x_s, sizeof(x_s), "%.16Qe", U[i][j]);
			res<<setw(28)<<x_s;
			UG[i][j]=U[i][j];
		}
		res<<endl;
	}
	
	Uyy=U[1][1];
   
    eps=geeta*tanq(rl);
    sk=sinq(geksi);
    ck=cosq(geksi);
    si=sinq(geeta);
    ci=cosq(geeta);
    se=sinq(eps);
    ce=cosq(eps);
    
    f21=sk*si*ce-ci*se;
    f22=sk*ci*ce+si*se;
    f23=ck*ce;

    f41=sk*si*se+ci*ce;
    f42=sk*ci*se+si*ce;
    f43=ck*se;

    f51=ck*si;
    f52=ck*ci;
    
    P[0][0]=f23*f23;
    P[0][1]=-2.*f23*f21;
    P[0][2]=2.*f23*f22;
    P[0][3]=f21*f21;
    P[0][4]=-2.*f21*f22;
    P[0][5]=f22*f22;
    
    P[1][0]=-f23*f43;
    P[1][1]=f23*f41+f43*f21;
    P[1][2]=-f23*f42-f43*f22;
    P[1][3]=-f21*f41;
    P[1][4]=f21*f42-f41*f22;
    P[1][5]=-f22*f42;
    
    P[2][0]=-f23*sk;
    P[2][1]=-f23*f51+sk*f21;
    P[2][2]=f23*f52-sk*f22;
    P[2][3]=f21*f51;
    P[2][4]=-f21*f52-f51*f22;
    P[2][5]=f22*f52;
    
    P[3][0]=f43*f43;
    P[3][1]=-2.*f43*f41;
    P[3][2]=2.*f43*f42;
    P[3][3]=f41*f41;
    P[3][4]=-2.*f41*f42;
    P[3][5]=f42*f42;
    
    P[4][0]=f43*sk;
    P[4][1]=f43*f51-sk*f41;
    P[4][2]=f43*f51+sk*f42;
    P[4][3]=-f41*f51;
    P[4][4]=f41*f52+f51*f42;
    P[4][5]=-f42*f52;
    
    P[5][0]=sk*sk;
    P[5][1]=2.*sk*f51;
    P[5][2]=-2.*sk*f52;
    P[5][3]=f51*f51;
    P[5][4]=-2.*f51*f52;
    P[5][5]=f52*f52;
    
//	res<<"\n  Computed Rotation matrix P \n";

    for(i=0;i<6;i++) 
		for (j=0;j<6;j++) 
			PW[i][j]=P[i][j];

    matinv6(PW,PIN);
    
    gu[0]=U[0][0];
    gu[1]=U[0][1];
    gu[2]=U[0][2];
    gu[3]=U[1][1];
    gu[4]=U[1][2];
    gu[5]=U[2][2];

	mulmat6(P,gu,ru);

    U[0][0]=ru[0];
    U[0][1]=ru[1];
    U[0][2]=ru[2];
    U[1][1]=ru[3];
    U[1][2]=ru[4];
    U[2][2]=ru[5];

	for(i=0;i<3;i++) 
		for (j=i;j<3;j++)  
			U[j][i]=U[i][j];
		
	res<<"\n  Rotated Normal Eotvos matrix U at point S - physical system (in sec^-2)\n\n";
	for(i=0;i<3;i++) {
		for (j=0;j<3;j++) {
            quadmath_snprintf(x_s, sizeof(x_s), "%.16Qe", U[i][j]);
			res<<setw(28)<<x_s;
		}
		res<<endl;
	}
	
	res<<"\n  Full Eotvos matrix W at point S - physical system (in sec^-2)\n\n";
	for(i=0;i<3;i++) {
		for (j=0;j<3;j++) {
			W[i][j]=U[i][j]+T[i][j];
            quadmath_snprintf(y_s, sizeof(y_s), "%.16Qe", W[i][j]);
			res<<setw(28)<<y_s;
		}	
		res<<endl;
	}
	res<<endl;
    
	difg=(om2-U[0][0]-U[1][1]-U[2][2])*eot;
    quadmath_snprintf(z_s, sizeof(z_s), "%.16Qe", difg);
	res<<"\n 2ù^2 -ÓUii = "<<z_s<<"  Eotvos ";
    
	difg=(om2-W[0][0]-W[1][1]-W[2][2])*eot;
    quadmath_snprintf(z_s, sizeof(z_s), "%.16Qe", difg);
	res<<"\n 2ù^2 -ÓWii = "<<z_s<<"  Eotvos "<<endl;    
    
    gnor<<fixed<<endl<<setw(3)<<sid<<endl;
	
	quadmath_snprintf(x_s, sizeof(x_s), "%.6Qe", lat);
	quadmath_snprintf(y_s, sizeof(y_s), "%.6Qe", lon);
	quadmath_snprintf(z_s, sizeof(z_s), "%.6Qe", h);
	quadmath_snprintf(w_s, sizeof(w_s), "%.6Qe", gH);
    gnor<<setw(16)<<x_s<<setw(16)<<y_s<<setw(16)<<w_s<<setw(16)<<z_s<<endl;
    
	quadmath_snprintf(y_s, sizeof(y_s), "%.16Qe", gS*mf);
	quadmath_snprintf(z_s, sizeof(z_s), "%.16Qe", dg*mf);
	quadmath_snprintf(w_s, sizeof(w_s), "%.16Qe", gamma*mf);
	gnor<<endl<<setw(26)<<y_s<<setw(26)<<z_s<<setw(26)<<w_s<<endl;            // on surface
	
	quadmath_snprintf(z_s, sizeof(z_s), "%.15Qe", C);
	quadmath_snprintf(w_s, sizeof(w_s), "%.15Qe", V);
	gnor<<setw(26)<<z_s<<setw(26)<<w_s<<endl<<endl;            // on surface
	
    gu[0]=T[0][0];
    gu[1]=T[0][1];
    gu[2]=T[0][2];
    gu[3]=T[1][1];
    gu[4]=T[1][2];
    gu[5]=T[2][2];

	mulmat6(PIN,gu,ru);

    TG[0][0]=ru[0];
    TG[0][1]=ru[1];
    TG[0][2]=ru[2];
    TG[1][1]=ru[3];
    TG[1][2]=ru[4];
    TG[2][2]=ru[5];

	for(i=0;i<3;i++) 
		for (j=i;j<3;j++)  
			TG[j][i]=TG[i][j];

	res<<"\n  Rotated Disturbance matrix T at point S - geometric system (in Eotvos units)\n\n";
	for(i=0;i<3;i++) {
		for (j=0;j<3;j++)  {
            quadmath_snprintf(x_s, sizeof(x_s), "%.16Qe", TG[i][j]*eot);
			res<<setw(28)<<x_s;
 		}
		res<<endl;
	}
	
	res<<"\n  Full Eotvos matrix W at point S - geometric system (in Eotvos units)\n\n";
	for(i=0;i<3;i++) {
		for (j=0;j<3;j++) {
			WG[i][j]=UG[i][j]+TG[i][j];
            quadmath_snprintf(y_s, sizeof(y_s), "%.16Qe", WG[i][j]*eot);
			res<<setw(28)<<y_s;
		}	
		res<<endl;
	}
	
	difg=(om2-UG[0][0]-UG[1][1]-UG[2][2])*eot;
    quadmath_snprintf(z_s, sizeof(z_s), "%.16Qe", difg);
	res<<"\n 2ù^2 -ÓUGii = "<<z_s<<"  Eotvos ";
    
	difg=(om2-WG[0][0]-WG[1][1]-WG[2][2])*eot;
    quadmath_snprintf(z_s, sizeof(z_s), "%.16Qe", difg);
	res<<"\n 2ù^2 -ÓWGii = "<<z_s<<"  Eotvos \n";    
    
	eps=sqrtq(geksi*geksi+geeta*geeta);
	w1[2]=-gS*cosq(eps);   // geometric
	w1[0]=w1[2]*geksi;
	w1[1]=w1[2]*geeta;

	difg=sqrtq(w1[0]*w1[0]+w1[1]*w1[1]+w1[2]*w1[2])-gS;
    quadmath_snprintf(x_s, sizeof(x_s), "%.10Qe", difg*mf);
	res<<"\n\n  ÓWGi^2 - gS = "<<x_s<<"  mgal  (geometric) \n";

//  Points 3 & 4 correspond to A' and B'

	g[3]=gS+(WG[0][0]*w1[0]+WG[0][1]*w1[1]+WG[0][2]*w1[2])*dP[0][0]/gS;
	g[0]=g[3]+(WG[0][2]*w1[0]+WG[1][2]*w1[1]+WG[2][2]*w1[2])*dP[0][2]/gS;

	g[4]=gS+(WG[1][0]*w1[0]+WG[1][1]*w1[1]+WG[1][2]*w1[2])*dP[1][1]/gS;
	g[1]=g[4]+(WG[0][2]*w1[0]+WG[1][2]*w1[1]+WG[2][2]*w1[2])*dP[1][2]/gS;

	g[2]=gS+(WG[0][2]*w1[0]+WG[1][2]*w1[1]+WG[2][2]*w1[2])*dP[2][2]/gS;
    
	res<<"\n Geometric system on surface:\n";

	for (n=0;n<3;n++)  {
//		difg=(g[n]-gS)*mf;
        res<<" Point "<<setw(3)<<n; 
        quadmath_snprintf(x_s, sizeof(x_s), "%.6Qe", dP[n][0]);
        quadmath_snprintf(y_s, sizeof(y_s), "%.6Qe", dP[n][1]);
        quadmath_snprintf(z_s, sizeof(z_s), "%.6Qe", dP[n][2]);
        quadmath_snprintf(w_s, sizeof(w_s), "%.16Qe", g[n]*mf);
        
  		res<<" :  x ="<<setw(16)<<x_s<<" m ,  y ="<<setw(16)<<y_s<<" m ,  z ="<<setw(16)<<z_s<<" m ";
 		res<<"\tg = "<<setw(26)<<w_s<<" mgal\n";

		gnor<<setw(16)<<x_s<<setw(16)<<y_s<<setw(16)<<z_s;
 		gnor<<setw(26)<<w_s<<endl;
	}
	
    quadmath_snprintf(z_s, sizeof(z_s), "%.16Qe", g[3]*mf);
    quadmath_snprintf(w_s, sizeof(w_s), "%.16Qe", g[4]*mf);
 	
	gnor<<endl<<setw(26)<<z_s<<endl;
	gnor<<setw(26)<<w_s<<endl;
	gnor<<endl;

//==================================

	w1[2]=-gS;		// physical
	w1[0]=zer;
	w1[1]=zer;

	difg=sqrtq(w1[0]*w1[0]+w1[1]*w1[1]+w1[2]*w1[2])-gS;
    quadmath_snprintf(x_s, sizeof(x_s), "%.10Qe", difg*mf);
	res<<"\n\n  ÓWi^2 - gS = "<<x_s<<"  mgal  (physical) \n";
    
//  Points 3 & 4 correspond to A' and B'

	gp[3]=gS-W[0][2]*dP[0][0];
	gp[0]=gp[3]-W[2][2]*dP[0][2]; 


	gp[4]=gS-W[1][2]*dP[1][1];
	gp[1]=gp[4]-W[2][2]*dP[1][2]; 

	gp[2]=gS-W[2][2]*dP[2][2];
    
	res<<"\n Physical system on surface:\n";     

	for (n=0;n<3;n++)  {
//		difg=(gp[n]-gS)*mf;
        res<<" Point "<<setw(3)<<n; 
        quadmath_snprintf(x_s, sizeof(x_s), "%.6Qe", dP[n][0]);
        quadmath_snprintf(y_s, sizeof(y_s), "%.6Qe", dP[n][1]);
        quadmath_snprintf(z_s, sizeof(z_s), "%.6Qe", dP[n][2]);
        quadmath_snprintf(w_s, sizeof(w_s), "%.16Qe", gp[n]*mf);
        
  		res<<" :  x ="<<setw(16)<<x_s<<" m ,  y ="<<setw(16)<<y_s<<" m ,  z ="<<setw(16)<<z_s<<" m ";
 		res<<"\tg = "<<setw(26)<<w_s<<" mgal\n";

		gnor<<setw(16)<<x_s<<setw(16)<<y_s<<setw(16)<<z_s;
 		gnor<<setw(26)<<w_s<<endl;
 	}
 	
    quadmath_snprintf(z_s, sizeof(z_s), "%.16Qe", gp[3]*mf);
    quadmath_snprintf(w_s, sizeof(w_s), "%.16Qe", gp[4]*mf);
 	
	gnor<<endl<<setw(26)<<z_s<<endl;
	gnor<<setw(26)<<w_s<<endl;
	gnor<<endl;

	for(i=0;i<3;i++) {             // on surface - physical
		for (j=0;j<3;j++)  {
            quadmath_snprintf(x_s, sizeof(x_s), "%.16Qe", T[i][j]*eot);
            gnor<<setw(28)<<x_s;
        }
		gnor<<endl;
	}
	gnor<<endl;
    
	for(i=0;i<3;i++) {             // on surface - physical
		for (j=0;j<3;j++)  {
            quadmath_snprintf(x_s, sizeof(x_s), "%.16Qe", U[i][j]*eot);
            gnor<<setw(28)<<x_s;
        }
		gnor<<endl;
	}
	gnor<<endl;
    
	for(i=0;i<3;i++) {             // on surface - geometric
		for (j=0;j<3;j++)  {
            quadmath_snprintf(x_s, sizeof(x_s), "%.16Qe", TG[i][j]*eot);
            gnor<<setw(28)<<x_s;
        }
		gnor<<endl;
	}
	gnor<<endl;
    
	for(i=0;i<3;i++) {             // on surface - geometric
		for (j=0;j<3;j++)  {
            quadmath_snprintf(x_s, sizeof(x_s), "%.16Qe", UG[i][j]*eot);
            gnor<<setw(28)<<x_s;
        }
		gnor<<endl;
	}

    quadmath_snprintf(x_s, sizeof(x_s), "%.16Qe", Uyy*eot);
    quadmath_snprintf(y_s, sizeof(y_s), "%.6Qe", ksi);
    quadmath_snprintf(z_s, sizeof(z_s), "%.6Qe", eta);
	
	gnor<<endl<<setw(28)<<x_s<<endl;          // on surface - geometric
	gnor<<setw(15)<<y_s<<setw(15)<<z_s<<endl<<endl;           //on surface
	
  }
  
    dat.close();
    res.close();
    gnor.close();
    
	return 0;
}
