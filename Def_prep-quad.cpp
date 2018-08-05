// Simulation of gravimetry
// based on EGM2008 data
// 05 August 2018
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

#define M3_DIM	3
#define M6_DIM	6
#define STR_SIZE	50
typedef __float128 quadfloat;

quadfloat const zer = 0.000000000000000000000000;
quadfloat const one= 1.000000000000000000000000;

void mulmat3(quadfloat a[M3_DIM][M3_DIM],quadfloat b[M3_DIM], quadfloat c[M3_DIM])
{
	int i,j;
	quadfloat s;

	for (i=0;i<M3_DIM;i++) {
		s=zer;
		for (j=0;j<M3_DIM;j++)
			s=s+a[i][j]*b[j];
		c[i]=s;
	}
}

void mulmat6(quadfloat a[M6_DIM][M6_DIM],quadfloat b[M6_DIM], quadfloat c[M6_DIM])
{
	int i,j;
	quadfloat s;

	for (i=0;i<M6_DIM;i++) {
		s=zer;
		for (j=0;j<M6_DIM;j++)
			s=s+a[i][j]*b[j];
		c[i]=s;
	}
}

void matinv6(quadfloat a[M6_DIM][M6_DIM], quadfloat b[M6_DIM][M6_DIM])
{
//  calculates the inverse of a matrix using Gauss-Jordan elimination
//  the inverse of matrix a is calculated and stored in the matrix b
//  beware:  matrix a  is converted to unit matrix!

    int i,j,k;
    quadfloat dum;

//  Build the unit matrix b
    for (i=0;i<M6_DIM ;i++)  {
      	for (j=0;j<M6_DIM ;j++)
      		b[i][j]=zer;
      	b[i][i]=one;
    }

	for (k=0;k<M6_DIM ;k++)  {

	  dum=a[k][k];
	  for (j=0;j<M6_DIM ;j++)  {
	  	  a[k][j]=a[k][j]/dum;
	  	  b[k][j]=b[k][j]/dum;
	  }

	  for (i=k+1;i<M6_DIM ;i++) {
	      dum=a[i][k];
	  	  for (j=0;j<M6_DIM ;j++) {
	  	  	  a[i][j]=a[i][j]-a[k][j]*dum;
	  	  	  b[i][j]=b[i][j]-b[k][j]*dum;
	  	  }
	  }
	}

	for (k=M6_DIM-1;k>0;k--)
	  for (i=k-1;i>=0;i--) {
	      dum=a[i][k];
	  	  for (j=0;j<M6_DIM ;j++) {
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
	quadfloat U[M3_DIM][M3_DIM],UG[M3_DIM][M3_DIM],T[M3_DIM][M3_DIM],TG[M3_DIM][M3_DIM],W[M3_DIM][M3_DIM];
	quadfloat dP[M3_DIM][M3_DIM],WG[M3_DIM][M3_DIM],w1[M3_DIM],g[5],gp[5];
	quadfloat Ux,Uz,kg,om2,Ad,Bd,Uyy;
	quadfloat eps,sk,ck,si,ci,se,ce;
	quadfloat f21,f22,f23,f41,f42,f43,f51,f52;
	quadfloat P[M6_DIM][M6_DIM],PW[M6_DIM][M6_DIM],PIN[M6_DIM][M6_DIM];
	quadfloat gu[M6_DIM],ru[M6_DIM];

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

	res<<endl;
    
	quadmath_snprintf(x_s, sizeof(x_s), "%.6Qe", lat);
	quadmath_snprintf(y_s, sizeof(y_s), "%.6Qe", lon);
	quadmath_snprintf(z_s, sizeof(z_s), "%.6Qe", h);
	quadmath_snprintf(w_s, sizeof(w_s), "%.6Qe", gH);

	res<<"\n Point "<<sid;
    res<<" (P) :  Latitude = "<<setw(14)<<x_s;
	res<<" deg  -  Longitude = "<<setw(14)<<y_s<<" deg \n";
	res<<"  - Orthometric Height = "<<setw(14)<<z_s;
	res<<"  -   Geometric Height = "<<setw(14)<<w_s<<" m \n"<<endl;
    
    
    
    
//  ==========================
    
    
    
	res<<"\n Ellipsoid parameters \n"<<scientific<<setprecision(12);
	res<<"  k1 = "<<k1<<" m^-1 \n";
	res<<"  k2 = "<<k2<<" m^-1 \n";
	res<<"   J = "<<J<<" m^-1 \n";
	res<<"   k = "<<k<<" m^-1 \n";
	
	res<<"\n Foot on ellipsoid\n  gamma (ell) = "<<gell*mf<<" mgal";
	
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

    gamma=sqrt(Ux*Ux+Uz*Uz);
	
	gS=gamma+dg;

	res<<"\n Point S  (on surface) \n  gamma = "<<gamma*mf<<" mgal";
	res<<"\n dg = "<<dg*mf<<" mgal";
	res<<"\n  g = "<<gS*mf<<" mgal\n";

	for(i=0;i<3;i++) 
		for (j=0;j<3;j++) 
			dat>>dP[i][j];
		
	for(i=0;i<3;i++) 
		for (j=i;j<3;j++)  {  //  on surface
			dat>>T[i][j];
			T[j][i]=T[i][j];
		}
	
	dat>>ksi>>eta;  //  on surface
	geksi=ksi/ras;
	geeta=eta/ras;
	
	res<<fixed<<setprecision(4);
	res<<"\n  Disturbance matrix  T at point P - physical system (in Eotvos units)\n\n";
	for(i=0;i<3;i++) {
		for (j=0;j<3;j++) {
			res<<setw(23)<<T[i][j];
			T[i][j]=T[i][j]/eot;
		}
		res<<endl;
	}
	
	U[0][0]=-gell*k1;
	U[0][2]=-gell*k;
	U[1][1]=-gell*k2;
	U[2][2]=2.*(omega*omega+gell*J);
	U[2][0]=U[0][2];

	res<<setprecision(6);
	res<<"\n  Constant part of Normal Eotvos matrix U at ellipsoid - geometric system (in Eotvos units)\n\n";
	for(i=0;i<3;i++) {
		for (j=0;j<3;j++) 
			res<<setw(25)<<U[i][j]*eot;
		res<<endl;
	}

	U[0][0]=U[0][0]+gell*k1*k1*gH; 
	U[0][2]=U[0][2]-(gell*k*(k2-k1)+2*gell*pJdx)*gH;
	U[1][1]=U[1][1]+gell*k2*k2*gH;  //  p2gdx  omitted
	U[2][2]=U[2][2]-gell*(k1*k1+k2*k2)*gH;  //  p2gdx  omitted
	U[2][0]=U[0][2];
		
	for(i=0;i<3;i++) 
		for (j=i;j<3;j++)  
			U[j][i]=U[i][j];

	res<<"\n  Complete Normal Eotvos matrix U at point P - geometric system (in Eotvos units)\n\n";
	for(i=0;i<3;i++) {
		for (j=0;j<3;j++) {
			res<<setw(25)<<U[i][j]*eot;
			UG[i][j]=U[i][j];
		}
		res<<endl;
	}
	Uyy=U[1][1];

//================================
    
    eps=geeta*tan(rl);
    sk=sin(geksi);
    ck=cos(geksi);
    si=sin(geeta);
    ci=cos(geeta);
    se=sin(eps);
    ce=cos(eps);
    
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
//	res<<scientific<<setprecision(7)<<endl;
	for(i=0;i<6;i++) {
		for (j=0;j<6;j++) {
			PW[i][j]=P[i][j];
//			res<<setw(16)<<PW[i][j];
		}
//		res<<endl;
	}
	matinv6(PW,PIN);
    
//	res<<"\n\n  Computed Inverse Rotation matrix PIN \n";
//	for(i=0;i<6;i++) {
//		for (j=0;j<6;j++) {
//			res<<setw(16)<<PIN[i][j];
//		}
//		res<<endl;
//	}

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
		
	res<<fixed<<setprecision(6);
	res<<"\n  Rotated Normal Eotvos matrix U at point P - physical system (in Eotvos units)\n\n";
	for(i=0;i<3;i++) {
		for (j=0;j<3;j++)  {
			res<<setw(25)<<U[i][j]*eot;
		}
		res<<endl;
	}
	
	res<<"\n  Full Eotvos matrix W at point P - physical system (in Eotvos units)\n\n";
	for(i=0;i<3;i++) {
		for (j=0;j<3;j++) {
			W[i][j]=U[i][j]+T[i][j];
			res<<setw(25)<<W[i][j]*eot;
		}	
		res<<endl;
	}
	res<<endl;
	difg=(om2-U[0][0]-U[1][1]-U[2][2])*eot;
	res<<"\n 2ù^2 -ÓUii = "<<difg<<"  Eotvos ";
	difg=(om2-W[0][0]-W[1][1]-W[2][2])*eot;
	res<<"\n 2ù^2 -ÓWii = "<<difg<<"  Eotvos \n";
	
	res<<endl;

	gnor<<fixed<<setprecision(7);
	gnor<<setw(3)<<sid<<setw(13)<<lat<<setw(13)<<lon<<setw(15)<<gH<<setw(15)<<h<<endl;
	gnor<<setprecision(13)<<setw(22)<<gS*mf<<setw(20)<<dg*mf<<setw(22)<<gamma*mf<<endl;  // on surface
	gnor<<setprecision(8)<<setw(20)<<C<<setprecision(13)<<setw(18)<<V<<endl;
//
//-------------------------------

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
		
	res<<fixed<<setprecision(6);
	res<<"\n  Rotated Disturbance matrix T at point P - geometric system (in Eotvos units)\n\n";
	for(i=0;i<3;i++) {
		for (j=0;j<3;j++)  {
			res<<setw(25)<<TG[i][j]*eot;
		}
		res<<endl;
	}
	
	res<<"\n  Full Eotvos matrix W at point P - geometric system (in Eotvos units)\n\n";
	for(i=0;i<3;i++) {
		for (j=0;j<3;j++) {
			WG[i][j]=UG[i][j]+TG[i][j];
			res<<setw(25)<<WG[i][j]*eot;
		}	
		res<<endl;
	}
	res<<endl;
	difg=(om2-UG[0][0]-UG[1][1]-UG[2][2])*eot;
	res<<"\n 2ù^2 -ÓUGii = "<<difg<<"  Eotvos ";
	difg=(om2-WG[0][0]-WG[1][1]-WG[2][2])*eot;
	res<<"\n 2ù^2 -ÓWGii = "<<difg<<"  Eotvos \n";
	
	res<<endl;

//================================

	eps=sqrt(geksi*geksi+geeta*geeta);
	w1[2]=-gS*cos(eps);   // geometric
	w1[0]=w1[2]*geksi;
	w1[1]=w1[2]*geeta;

	difg=sqrt(w1[0]*w1[0]+w1[1]*w1[1]+w1[2]*w1[2])-gS;
	res<<"\n  ÓWGi^2 - gS = "<<difg*mf<<"  mgal  (geometric) \n";

//  Points 3 & 4 correspond to A' and B'

	g[3]=gS+(WG[0][0]*w1[0]+WG[0][1]*w1[1]+WG[0][2]*w1[2])*dP[0][0]/gS;
	g[0]=g[3]+(WG[0][2]*w1[0]+WG[1][2]*w1[1]+WG[2][2]*w1[2])*dP[0][2]/gS;

	g[4]=gS+(WG[1][0]*w1[0]+WG[1][1]*w1[1]+WG[1][2]*w1[2])*dP[1][1]/gS;
	g[1]=g[4]+(WG[0][2]*w1[0]+WG[1][2]*w1[1]+WG[2][2]*w1[2])*dP[1][2]/gS;

	g[2]=gS+(WG[0][2]*w1[0]+WG[1][2]*w1[1]+WG[2][2]*w1[2])*dP[2][2]/gS;

	for (n=0;n<3;n++)  {
		difg=(g[n]-gS)*mf;
		res<<"\n Geometric system on surface:";
		res<<"\n Point "<<n<<" :  x ="<<setw(10)<<dP[n][0]<<" m ,  y ="<<setw(10)<<dP[n][1]<<" m ,  z ="<<setw(10)<<dP[n][2]<<" m ";
		res<<"\tg = "<<g[n]*mf<<" mgal \t g("<<n<<") - gS = "<<difg<<" mgal \n";
		gnor<<setprecision(3);
		for(i=0;i<3;i++)
			gnor<<setw(10)<<dP[n][i];
		gnor<<setprecision(13)<<setw(23)<<g[n]*mf<<endl;
	}
	gnor<<setw(23)<<g[3]*mf<<endl;
	gnor<<setw(23)<<g[4]*mf<<endl;
	gnor<<setprecision(6)<<endl;

//==================================

	w1[2]=-gS;		// physical
	w1[0]=0.0;
	w1[1]=0.0;

	difg=sqrt(w1[0]*w1[0]+w1[1]*w1[1]+w1[2]*w1[2])-gS;
	res<<"\n  ÓWi^2 - gS = "<<difg*mf<<"  mgal  (physical) \n";

//  Points 3 & 4 correspond to A' and B'

	gp[3]=gS-W[0][2]*dP[0][0];
	gp[0]=gp[3]-W[2][2]*dP[0][2]; 


	gp[4]=gS-W[1][2]*dP[1][1];
	gp[1]=gp[4]-W[2][2]*dP[1][2]; 

	gp[2]=gS-W[2][2]*dP[2][2];

	for (n=0;n<3;n++)  {
		difg=(gp[n]-gS)*mf;
		res<<"\n Physical system on surface:";
		res<<"\n Point "<<n<<" :  x ="<<setw(10)<<dP[n][0]<<" m ,  y ="<<setw(10)<<dP[n][1]<<" m ,  z ="<<setw(10)<<dP[n][2]<<" m ";
		res<<"\tg = "<<gp[n]*mf<<" mgal \t g("<<n<<") - gS = "<<difg<<" mgal \n";
		gnor<<setprecision(3);
		for(i=0;i<3;i++)
			gnor<<setw(10)<<dP[n][i];
		gnor<<setprecision(13)<<setw(23)<<gp[n]*mf<<endl;
	}
	gnor<<setw(23)<<gp[3]*mf<<endl;
	gnor<<setw(23)<<gp[4]*mf<<endl;
	gnor<<setprecision(6)<<endl;

//-------------------------------------

	for(i=0;i<3;i++) {
		for (j=0;j<3;j++) 
			gnor<<setw(22)<<T[i][j]*eot;  // on surface - physical
		gnor<<endl;
	}
	gnor<<endl;
	for(i=0;i<3;i++) {
		for (j=0;j<3;j++) 
			gnor<<setw(22)<<U[i][j]*eot;  // on surface - physical
		gnor<<endl;
	}
	gnor<<endl;
	for(i=0;i<3;i++) {
		for (j=0;j<3;j++) 
			gnor<<setw(22)<<UG[i][j]*eot;  // on surface - geometric
		gnor<<endl;
	}
	gnor<<endl;
	for(i=0;i<3;i++) {
		for (j=0;j<3;j++) 
			gnor<<setw(22)<<TG[i][j]*eot;  // on surface - geometric
		gnor<<endl;
	}

	gnor<<setw(27)<<Uyy*eot<<endl;  // on surface - geometric
	gnor<<setprecision(6)<<setw(10)<<ksi<<setw(10)<<eta<<endl<<endl;  //on surface
	
  }
	return 0;
}
