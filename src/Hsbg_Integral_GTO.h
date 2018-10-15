#ifndef Hsbg_Integral_GTO_H
#define Hsbg_Integral_GTO_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include <cmath>		// erf
#include <Eigen/Dense>

#include "Hsbg_Const.h"
#include "Hsbg_Tools.h"
#include "Hsbg_Global.h"
#include "Hsbg_Tasker.h"


using namespace std;
using namespace Hsbg;
using namespace Eigen;



//////////////////////////////////////////////////////////////////////////////////////////////////
/*

		the general integral function as follow
        all integral formula refer <<Modern quantum chemistry computational methods>>, Zhizhong Wang.
        and << quantum chenistry >> by Guangxian Xu, Lemin Li.
*/
//////////////////////////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////////////////////////
// declarition of functions
double gto_normfactor(double& alpha, int& l, int& m, int& n);
int parity(const int& n);
int comb(const int& a, const int& b);
int fact( const int& n );
int dfact(const int& n);
double Fn(const int& m, const double& w);
double fun_f(const int& l, const int& l1, const int& l2, const double& a, const double& b);
double fun_G
    (const int& I, const int& l1, const int& l2, const double& a, const double& b, const double& c, const double& w);
double fun_H
    (const int& l, const int& l1, const int& l2, const double& a, const double& b, const double& w);
double fun_D
    (const int& I, const int& l1, const int& l2, const int& l3, const int& l4,
    const double& PQ, const double& a, const double& b, const double& c, const double& d,
    const double& del,  const double& a1, const double& a2, const double& a3, const double& a4 );
double gto_S
(
const double& x1, const double& y1, const double& z1, const double& a1, const int& l1, const int& m1, const int& n1,
const double& x2, const double& y2, const double& z2, const double& a2, const int& l2, const int& m2, const int& n2
);
double gto_T
(
const double& x1, const double& y1, const double& z1, const double& a1, const int& l1, const int& m1, const int& n1,
const double& x2, const double& y2, const double& z2, const double& a2, const int& l2, const int& m2, const int& n2
);
double gto_V
(
const double& x1, const double& y1, const double& z1, const double& a1, const int& l1, const int& m1, const int& n1,
const double& x2, const double& y2, const double& z2, const double& a2, const int& l2, const int& m2, const int& n2,
const double& xc, const double& yc, const double& zc
);
double gto_ERI
(
const double& x1, const double& y1, const double& z1, const double& a1, const int& l1, const int& m1, const int& n1,
const double& x2, const double& y2, const double& z2, const double& a2, const int& l2, const int& m2, const int& n2,
const double& x3, const double& y3, const double& z3, const double& a3, const int& l3, const int& m3, const int& n3,
const double& x4, const double& y4, const double& z4, const double& a4, const int& l4, const int& m4, const int& n4
);
//////////////////////////////////////////////////////////////////////////////////////////////////


double gto_normfactor(double& alpha, int& l, int& m, int& n)
{
    int lmn = l + m + n;
    return sqrt( pow(2,2*lmn) * pow(alpha,lmn+1.5) * pow(2/PI,1.5) /
           ( dfact(2*l-1) * dfact(2*m-1) * dfact(2*n-1) )          );
}

int parity(const int& n)
{
    return (n%2==0 ? 1: -1);
}

int comb(const int& a, const int& b)
{
	return (a <= b || b==0 ) ? 1 : comb(a-1,b) + comb(a-1,b-1) ;
}

int fact( const int& n )
{
    return (n == 1 || n <= 0) ? 1 : fact(n - 1) * n;
}

int dfact(const int& n)
{
	return (n == 1 || n <= 0) ? 1 : dfact(n - 2) * n;
}

double fun_f(const int& l, const int& l1, const int& l2, const double& a, const double& b)
{
    double result = 0;
    for(int i1=0; i1<= l1; i1++)
    {
        for(int i2=0; i2<= l2; i2++)
        {
            if(i1+i2==l)
            {
                result += comb(l1,i1) * comb(l2,i2) * pow(a,l1-i1) * pow(b,l2-i2);
            }
        }
    }
    return result;
}

double fun_G
(const int& I, const int& l1, const int& l2, const double& a, const double& b, const double& c, const double& w)
{
    double result = 0;
    for(int i=0; i<= l1+l2; i++)
    {
        for(int r=0; r<=i/2; r++)
        {
            for(int u=0; u<=(i-2*r)/2; u++)
            {
                if(i-2*r-u==I) // NOT SURE
                {
                    result +=
                    (
                        parity(i) * fun_f(i,l1,l2,a,b)
                      * parity(u) * fact(i) / double( fact(r)*fact(u)*fact(i-2*r-2*u) )
                      * pow(c, i-2*r-2*u) * pow(0.25/w, r+u)
                    );
                }
            }
        }
    }
    return result;
}

double fun_H
(const int& l, const int& l1, const int& l2, const double& a, const double& b, const double& w)
{
    double result = 0;
    for(int i=0; i<=l1+l2; i++)
    {
        for(int r=0; r<= i/2; r++)
        {
            if(i-2*r == l) // NOT SURE
            {
                result +=
                (
                    fact(i)* fun_f(i,l1,l2,a,b)
                  / double( fact(r) * fact(i-2*r) * pow(4*w, i-r) )
                );
            }
        }
    }
    return result;
}

double fun_D
(const int& I, const int& l1, const int& l2, const int& l3, const int& l4,
const double& PQ, const double& a, const double& b, const double& c, const double& d,
const double& del,  const double& a1, const double& a2, const double& a3, const double& a4 )
{
    double result = 0;
    for(int l=0; l<=l1+l2; l++)
    {
        for(int lp=0; lp<=l3+l4; lp++)
        {
            for(int u=0; u<= (l+lp)/2; u++)
            {
                if(l+lp-u == I) // NOT SURE
                {
                    result +=
                    (
                        parity(u+lp) * fact(l+lp) * pow(PQ,l+lp-2*u)
                      / double(fact(u) * fact(l+lp-2*u) * pow(del, l+lp-u) )
                      * fun_H(l,l1,l2,a,b,a1+a2) * fun_H(lp,l3,l4,c,d,a3+a4)
                    );
                }
            }
        }
    }
    return result;
}




double gto_S
(
const double& x1, const double& y1, const double& z1, const double& a1, const int& l1, const int& m1, const int& n1,
const double& x2, const double& y2, const double& z2, const double& a2, const int& l2, const int& m2, const int& n2
)
{	
    double AB2, PAx, PBx, PAy, PBy, PAz, PBz;
    double kcoeff, a12, Sx, Sy, Sz;
    AB2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
    a12 = a1*a2/(a1+a2);
    kcoeff = exp( - a12 * AB2 );

    if(l1<0 || l2<0 || m1<0 || m2<0 || n1<0 || n2 < 0) return 0;
    
    // x direction
    Sx = 0;
    PAx = a2*(x2-x1)/(a1+a2);
    PBx = a1*(x1-x2)/(a1+a2);
    for(int i=0; i<=(l1+l2)/2; i++)
	{
		Sx += fun_f(2*i,l1,l2,PAx,PBx) * dfact(2*i-1) / pow(2*(a1+a2), i);
	}
	// y direction
	Sy = 0;
    PAy = a2*(y2-y1)/(a1+a2);
    PBy = a1*(y1-y2)/(a1+a2);
    for(int i=0; i<=(m1+m2)/2; i++)
	{
		Sy += fun_f(2*i,m1,m2,PAy,PBy) * dfact(2*i-1) / pow(2*(a1+a2), i);
	}
    // z direction
	Sz = 0;
    PAz = a2*(z2-z1)/(a1+a2);
    PBz = a1*(z1-z2)/(a1+a2);
    for(int i=0; i<=(n1+n2)/2; i++)
	{
		Sz += fun_f(2*i,n1,n2,PAz,PBz) * dfact(2*i-1) / pow(2*(a1+a2), i);
	}

	return pow(PI/(a1+a2), 1.5) * kcoeff * Sx * Sy * Sz;
}


double gto_T
(
const double& x1, const double& y1, const double& z1, const double& a1, const int& l1, const int& m1, const int& n1,
const double& x2, const double& y2, const double& z2, const double& a2, const int& l2, const int& m2, const int& n2
)
{
    double S00, Spp, Smm;
    int lmn2 = l2+m2+n2;
    S00 = a2 * (2*lmn2+3) *gto_S(x1,y1,z1,a1,l1,m1,n1, x2,y2,z2,a2,l2,m2,n2);
    Spp = 2*a2*a2*(
                    gto_S(x1,y1,z1,a1,l1,m1,n1, x2,y2,z2,a2,l2+2,m2,n2)
                  + gto_S(x1,y1,z1,a1,l1,m1,n1, x2,y2,z2,a2,l2,m2+2,n2)
                  + gto_S(x1,y1,z1,a1,l1,m1,n1, x2,y2,z2,a2,l2,m2,n2+2) 
                  );
    Smm = 0.500 * (
                    gto_S(x1,y1,z1,a1,l1,m1,n1, x2,y2,z2,a2,l2-2,m2,n2)*l2*(l2-1)
                  + gto_S(x1,y1,z1,a1,l1,m1,n1, x2,y2,z2,a2,l2,m2-2,n2)*m2*(m2-1)
                  + gto_S(x1,y1,z1,a1,l1,m1,n1, x2,y2,z2,a2,l2,m2,n2-2)*n2*(n2-1)
                  );
	return S00 - Spp - Smm;
}


double gto_V
(
const double& x1, const double& y1, const double& z1, const double& a1, const int& l1, const int& m1, const int& n1,
const double& x2, const double& y2, const double& z2, const double& a2, const int& l2, const int& m2, const int& n2,
const double& xc, const double& yc, const double& zc
)
{
    double AB2, PC2, PAx, PBx, PCx, PAy, PBy, PCy, PAz, PBz, PCz;
    double kcoeff, a12, w, V;
    int N;
	
    AB2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
    a12 = a1*a2/(a1+a2);
    w = a1+a2;
    kcoeff = exp( -a12 * AB2 );
    
    N = l1+m1+n1+l2+m2+n2;
    
    PAx = a2*(x2-x1)/(a1+a2);
    PBx = a1*(x1-x2)/(a1+a2);
    PCx = (a1*x1+a2*x2)/(a1+a2) - xc;
    PAy = a2*(y2-y1)/(a1+a2);
    PBy = a1*(y1-y2)/(a1+a2);
    PCy = (a1*y1+a2*y2)/(a1+a2) - yc;
    PAz = a2*(z2-z1)/(a1+a2);
    PBz = a1*(z1-z2)/(a1+a2);
    PCz = (a1*z1+a2*z2)/(a1+a2) - zc;
    PC2 = PCx*PCx + PCy*PCy + PCz*PCz;
    
    V = 0;
    for(int nu=0; nu<=N; nu++)
    {
        double term = 0;
        for(int I=0; I<= l1+l2; I++)
        {
            for(int J=0; J<= m1+m2; J++)
            {
                for(int K=0; K<= n1+n2; K++)
                {
                    if(I+J+K == nu) // NOT SURE
                    {
                        term += (
                               fun_G(I,l1,l2,PAx,PBx,PCx,w)
                             * fun_G(J,m1,m2,PAy,PBy,PCy,w)
                             * fun_G(K,n1,n2,PAz,PBz,PCz,w)
                             );
                    }
                }
            }
        }
        V += term * Fn(nu, (a1+a2)*PC2);
    }
    return 2*PI/(a1+a2) * kcoeff * V;
}

double gto_ERI
(
const double& x1, const double& y1, const double& z1, const double& a1, const int& l1, const int& m1, const int& n1,
const double& x2, const double& y2, const double& z2, const double& a2, const int& l2, const int& m2, const int& n2,
const double& x3, const double& y3, const double& z3, const double& a3, const int& l3, const int& m3, const int& n3,
const double& x4, const double& y4, const double& z4, const double& a4, const int& l4, const int& m4, const int& n4
)
{
    double AB2, CD2, PAx, PBx, PAy, PBy, PAz, PBz, QCx, QDx, QCy, QDy, QCz, QDz;
    double PQ2, PQx, PQy, PQz;
    double a12, a34, d;
    double sum, ERI;
    int N1, N2, Ilim, Jlim, Klim;
	
    a12 = a1*a2/(a1+a2);
    a34 = a3*a4/(a3+a4);
    N1 = l1+m1+n1+l2+m2+n2;
    N2 = l3+m3+n3+l4+m4+n4;
    Ilim = l1+l2+l3+l4;
    Jlim = m1+m2+m3+m4;
    Klim = n1+n2+n3+n4;
    d = (a1+a2+a3+a4)/( 4*(a1+a2)*(a3+a4) );
    
    AB2 = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
    CD2 = (x3-x4)*(x3-x4)+(y3-y4)*(y3-y4)+(z3-z4)*(z3-z4);
    
    PAx = a2*(x2-x1)/(a1+a2);
    PBx = a1*(x1-x2)/(a1+a2);
    PAy = a2*(y2-y1)/(a1+a2);
    PBy = a1*(y1-y2)/(a1+a2);
    PAz = a2*(z2-z1)/(a1+a2);
    PBz = a1*(z1-z2)/(a1+a2);
    
    QCx = a4*(x4-x3)/(a3+a4);
    QDx = a3*(x3-x4)/(a3+a4);
    QCy = a4*(y4-y3)/(a3+a4);
    QDy = a3*(y3-y4)/(a3+a4);
    QCz = a4*(z4-z3)/(a3+a4);
    QDz = a3*(z3-z4)/(a3+a4);
    
    PQx = (a1*x1+a2*x2)/(a1+a2) - (a3*x3+a4*x4)/(a3+a4);
    PQy = (a1*y1+a2*y2)/(a1+a2) - (a3*y3+a4*y4)/(a3+a4);
    PQz = (a1*z1+a2*z2)/(a1+a2) - (a3*z3+a4*z4)/(a3+a4);
    PQ2 = PQx*PQx + PQy*PQy + PQz*PQz;
    
    sum = 0;
    for(int nu=0; nu<= N1+N2; nu++)
    {
        double term = 0;
        for(int I=0; I<=Ilim; I++)
        {
            for(int J=0; J<= Jlim; J++)
            {
                for(int K=0; K<=Klim; K++)
                {
                    if(I+J+K == nu)
                    {
                        term +=
                        (
                            fun_D(I,l1,l2,l3,l4,PQx,PAx,PBx,QCx,QDx,d,a1,a2,a3,a4)
                          * fun_D(J,m1,m2,m3,m4,PQy,PAy,PBy,QCy,QDy,d,a1,a2,a3,a4)
                          * fun_D(K,n1,n2,n3,n4,PQz,PAz,PBz,QCz,QDz,d,a1,a2,a3,a4)
                        );
                    }
                }
            }
        }
        sum += term * Fn(nu, PQ2/(4*d));
    }
    
    ERI = sum * 2*PI*PI/( (a1+a2)*(a3+a4) ) * sqrt(PI/(a1+a2+a3+a4)) * exp(-a12*AB2-a34*CD2);
    return ERI;
}


double Fn(const int& m, const double& w)
{
	double a1, a2, a3, a4, a5, b1, b2, b3, b4, b5, b6;
	double a0 = pow( (double)(2*m+1), -2/((double)(2*m+1)) );
	if(m==0){
		if(w>16.3578)
		{
			return (dfact(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2);
		}else{
			a1=0.213271302431420;
			a2=0.629344460255614E-1;
			a3=0.769838037756759E-2;
			a4=0.758433197127160E-3;
			a5=0.564691197633667E-4;
			b1=0.879937801660182;
			b2=0.338450368470103;
			b3=0.738522953299624E-1;
			b4=0.101431553402629E-1;
			b5=0.955528842975585E-3;
			b6=0.720266520392572E-4;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) + a3*pow(w,3) + a4*pow(w,4) + a5*pow(w,5)) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) + b4*pow(w,4) + b5*pow(w,5) + b6*pow(w,6));
			return pow(temp, m+0.5);
		}
	}
	if(m==1){
		if(w>17.4646)
		{
			return (dfact(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2);
		}else{
			a1=0.295195994716045E-1;
			a2=0.128790985465415E-1;
			a3=0.998165499553218E-3;
			a4=0.970927983276419E-4;
			a5=0.493839847029699E-5;
			b1=0.461403194579124;
			b2=0.108494164372449;
			b3=0.171462934845042E-1;
			b4=0.196918657845508E-2;
			b5=0.160138863265254E-3;
			b6=0.857708713007233E-5;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) + a3*pow(w,3) + a4*pow(w,4) + a5*pow(w,5)) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) + b4*pow(w,4) + b5*pow(w,5) + b6*pow(w,6));
			return pow(temp, m+0.5);
		}
	}
	if(m==2){
		if(w>15.2368)
		{
			return (dfact(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2);
		}else{
			a1=-0.575763488635418E-2;
			a2=0.731474973333076E-2;
			a3=0.251276149443393E-3;
			a4=0.264336244559094E-4;
			a5=0;
			b1=0.274754154712841;
			b2=0.425364830353043E-1;
			b3=0.493902790955943E-2;
			b4=0.437251500927601E-3;
			b5=0.288914662393981E-4;
			b6=0;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) + a3*pow(w,3) + a4*pow(w,4) + a5*pow(w,5)) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) + b4*pow(w,4) + b5*pow(w,5) + b6*pow(w,6));
			return pow(temp, m+0.5);
		}
	}
	if(m==3){
		if(w>16.0419)
		{
			return (dfact(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2);
		}else{
			a1=-0.290110430424666E-1;
			a2=0.561884370781462E-2;
			a3=0.301628267382713E-4;
			a4=0.110671035361856E-4;
			a5=0;
			b1=0.171637608242892;
			b2=0.187571417256877E-1;
			b3=0.178536829675118E-2;
			b4=0.137360778130936E-3;
			b5=0.791915206883054E-5;
			b6=0;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) + a3*pow(w,3) + a4*pow(w,4) + a5*pow(w,5)) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) + b4*pow(w,4) + b5*pow(w,5) + b6*pow(w,6));
			return pow(temp, m+0.5);
		}
	}
	if(m==4){
		if(w>16.8955)
		{
			return (dfact(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
		}else{
			a1=-0.452693111179624E-1;
			a2=0.490070062899003E-2;
			a3=-0.561789719979307E-4;
			a4=0.550814626951998E-5;
			a5=0;
			b1=0.108051989937231;
			b2=0.855924943430755E-2;
			b3=0.724968571389473E-3;
			b4=0.502338223156067E-4;
			b5=0.249107837399141E-5;
			b6=0;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) + a3*pow(w,3) + a4*pow(w,4) + a5*pow(w,5)) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) + b4*pow(w,4) + b5*pow(w,5) + b6*pow(w,6));
			return pow(temp, m+0.5);
		}
	}
	if(m==5){
		if(w>17.7822)
		{
			return (dfact(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
		}else{
			a1=-0.566143259316101E-1;
			a2=0.455916894577203E-2;
			a3=-0.894152721395639E-4;
			a4=0.328096732308082E-5;
			a5=0;
			b1=0.662932958471386E-1;
			b2=0.383724443872493E-2;
			b3=0.327167659811839E-3;
			b4=0.210430437682548E-4;
			b5=0.883562935089333E-6;
			b6=0;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) + a3*pow(w,3) + a4*pow(w,4) + a5*pow(w,5)) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) + b4*pow(w,4) + b5*pow(w,5) + b6*pow(w,6));
			return pow(temp, m+0.5);
		}
	}
	if(m==6){
		if(w>15.8077)
		{
			return (dfact(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
		}else{
			a1=-0.503249167534352E-1;
			a2=0.273135625430953E-2;
			a3=-0.310733624819100E-4;
			b1=0.586609328033372E-1;
			b2=0.194044691497128E-2;
			b3=0.109442742502192E-3;
			b4=0.613406236401726E-5;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) + a3*pow(w,3) ) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) + b4*pow(w,4) );
			return pow(temp, m+0.5);
		}
	}
	if(m==7){
		if(w>16.5903)
		{
			return (dfact(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
		}else{
			a1=-0.548201062615785E-1;
			a2=0.253099908233175E-2;
			a3=-0.333589469427863E-4;
			b1=0.389873128779298E-1;
			b2=0.569890860832729E-3;
			b3=0.422187129333708E-4;
			b4=0.286010059144633E-5;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) + a3*pow(w,3) ) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) + b4*pow(w,4) );
			return pow(temp, m+0.5);
		}
	}
	if(m==8){
		if(w>17.3336)
		{
			return (dfact(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
		}else{
			a1=-0.581618006078160E-1;
			a2=0.238525529084601E-2;
			a3=-0.329989020317093E-4;
			b1=0.240929282666615E-1;
			b2=-0.202677647499956E-3;
			b3=0.119820675974460E-4;
			b4=0.145762086904409E-5;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) + a3*pow(w,3) ) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) + b4*pow(w,4) );
			return pow(temp, m+0.5);
		}
	}
	if(m==9){
		if(w>15.6602)
		{
			return (dfact(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
		}else{
			a1=-0.334843993901400E-1;
			a2=0.846637494147059E-3;
			b1=0.495875606944471E-1;
			b2=0.946642302340943E-3;
			b3=0.108367772249790E-4;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) ) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) );
			return pow(temp, m+0.5);
		}
	}
	if(m==10){
		if(w>16.5258)
		{
			return (dfact(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
		}else{
			a1=-0.335292171805959E-1;
			a2=0.749168957990503E-3;
			b1=0.421492932021683E-1;
			b2=0.582840985360327E-3;
			b3=0.237676790577455E-5;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) ) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) );
			return pow(temp, m+0.5);
		}
	}
	if(m==11){
		if(w>17.5395)
		{
			return (dfact(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
		}else{
			a1=-0.332669773770348E-1;
			a2=0.668720489602687E-3;
			b1=0.363057685289467E-1;
			b2=0.345646100984643E-3;
			b3=-0.190872330373450E-5;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) ) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) );
			return pow(temp, m+0.5);
		}
	}
	if(m==12){
		if(w>18.5783)
		{
			return (dfact(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
		}else{
			a1=-0.326241966410798E-1;
			a2=0.598705175467956E-3;
			b1=0.318680048277695E-1;
			b2=0.202419662347765E-3;
			b3=-0.362095173837973E-5;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) ) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) );
			return pow(temp, m+0.5);
		}
	}
	if(m==13){
		if(w>19.6511)
		{
			return (dfact(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
		}else{
			a1=-0.317754368014894E-1;
			a2=0.537678595933584E-3;
			b1=0.284036027081815E-1;
			b2=0.113673420662576E-3;
			b3=-0.416076810552774E-5;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) ) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) );
			return pow(temp, m+0.5);
		}
	}
	if(m==14){
		if(w>20.7839)
		{
			return (dfact(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
		}else{
			a1=-0.308755854748829E-1;
			a2=0.485046451960769E-3;
			b1=0.255694625434059E-1;
			b2=0.542010192055080E-4;
			b3=-0.424759498527876E-5;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) ) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) );
			return pow(temp, m+0.5);
		}
	}
	if(m==15){
		if(w>21.9998)
		{
			return (dfact(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
		}else{
			a1=-0.3001438066719997E-1;
			a2=0.439983033427912E-3;
			b1=0.231478878674366E-1;
			b2=0.105546581596674E-4;
			b3=-0.418932957034726E-5;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) ) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) );
			return pow(temp, m+0.5);
		}
	}
	if(m==16){
		if(w>20.9225)
		{
			return (dfact(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
		}else{
			a1=-0.288346417991609E-1;
			a2=0.397161796318408E-3;
			b1=0.215021933120724E-1;
			b2=-0.128592457453950E-5;
			b3=-0.362120651688135E-5;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) ) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) );
			return pow(temp, m+0.5);
		}
	}
	return ( (2*m - 1) * Fn(m-1,w) - exp(-w) ) / (2*w);
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// calculate constract GTOs

double IntecGTO_S(
	Orbital_cgto& cgto1,
	Orbital_cgto& cgto2
)
{
	double sum=0;
	double x1,x2,y1,y2,z1,z2,a1,a2;
	int L1,L2,M1,M2,N1,N2;
	
	x1=cgto1.x; y1=cgto1.y; z1=cgto1.z; L1=cgto1.L; M1=cgto1.M; N1=cgto1.N;
    x2=cgto2.x; y2=cgto2.y; z2=cgto2.z; L2=cgto2.L; M2=cgto2.M; N2=cgto2.N;
		    
	for(int i = 0; i < cgto1.cn; i++){
		for(int j = 0; j < cgto2.cn; j++){
		    a1=cgto1.alphas[i];
		    a2=cgto2.alphas[j];
		    
			sum += cgto1.normGTO(i) * cgto2.normGTO(j) *
			cgto1.coeffs[i] * cgto2.coeffs[j] *  (
			gto_S( x1,y1,z1,a1,L1,M1,N1,
			       x2,y2,z2,a2,L2,M2,N2	 ) );
		}
	}
	return sum;
}

double IntecGTO_T(
	Orbital_cgto& cgto1,
	Orbital_cgto& cgto2
)
{
	double sum=0;
	double x1,x2,y1,y2,z1,z2,a1,a2;
	int L1,L2,M1,M2,N1,N2;
	
	x1=cgto1.x; y1=cgto1.y; z1=cgto1.z; L1=cgto1.L; M1=cgto1.M; N1=cgto1.N;
    x2=cgto2.x; y2=cgto2.y; z2=cgto2.z; L2=cgto2.L; M2=cgto2.M; N2=cgto2.N;

	for(int i = 0; i < cgto1.cn; i++){
		for(int j = 0; j < cgto2.cn; j++){
		
		    a1=cgto1.alphas[i];
		    a2=cgto2.alphas[j];
		
			sum +=  cgto1.normGTO(i) * cgto2.normGTO(j) *
			cgto1.coeffs[i] * cgto2.coeffs[j] *  (
			gto_T( x1,y1,z1,a1,L1,M1,N1,
			       x2,y2,z2,a2,L2,M2,N2	 ) );
		}
	}
	return sum;
}

double IntecGTO_V(
	Orbital_cgto& cgto1,
	Orbital_cgto& cgto2,
	Atom& PN
)
{
	double sum=0;
	double x1,x2,y1,y2,z1,z2,a1,a2, xc,yc,zc;
	int L1,L2,M1,M2,N1,N2;
	xc = PN.x; yc = PN.y; zc = PN.z;
	
	x1=cgto1.x; y1=cgto1.y; z1=cgto1.z; L1=cgto1.L; M1=cgto1.M; N1=cgto1.N;
    x2=cgto2.x; y2=cgto2.y; z2=cgto2.z; L2=cgto2.L; M2=cgto2.M; N2=cgto2.N;
    
	for(int i = 0; i < cgto1.cn; i++){
		for(int j = 0; j < cgto2.cn; j++){
		    a1=cgto1.alphas[i];
		    a2=cgto2.alphas[j];
		    
			sum += cgto1.normGTO(i) * cgto2.normGTO(j) *
			cgto1.coeffs[i] * cgto2.coeffs[j] *  (
			gto_V( x1,y1,z1,a1,L1,M1,N1,
			       x2,y2,z2,a2,L2,M2,N2,
			       xc,yc,zc ) );
		}
	}
	return sum * (-1) * PN.znum;
}

double IntecGTO_ERI(
	Orbital_cgto& cgto1,
	Orbital_cgto& cgto2,
	Orbital_cgto& cgto3,
	Orbital_cgto& cgto4
)
{
	double sum=0;
	double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,a1,a2,a3,a4;
	int L1,L2,L3,L4,M1,M2,M3,M4,N1,N2,N3,N4;
	
	x1=cgto1.x; y1=cgto1.y; z1=cgto1.z; L1=cgto1.L; M1=cgto1.M; N1=cgto1.N;
	x2=cgto2.x; y2=cgto2.y; z2=cgto2.z; L2=cgto2.L; M2=cgto2.M; N2=cgto2.N;
	x3=cgto3.x; y3=cgto3.y; z3=cgto3.z; L3=cgto3.L; M3=cgto3.M; N3=cgto3.N;
	x4=cgto4.x; y4=cgto4.y; z4=cgto4.z; L4=cgto4.L; M4=cgto4.M; N4=cgto4.N;
	
	for(int i = 0; i < cgto1.cn; i++)
	{
		for(int j = 0; j < cgto2.cn; j++)
		{
			for(int k=0; k< cgto3.cn; k++)
			{
				for(int m=0; m < cgto4.cn; m++)
				{
				    a1=cgto1.alphas[i];
		            a2=cgto2.alphas[j];
		            a3=cgto3.alphas[k];
		            a4=cgto4.alphas[m]; 
		            				    
					sum += cgto1.normGTO(i) * cgto2.normGTO(j) * cgto3.normGTO(k) * cgto4.normGTO(m) *
					cgto1.coeffs[i] * cgto2.coeffs[j] * cgto3.coeffs[k] * cgto4.coeffs[m]  *(
					gto_ERI(
			               x1,y1,z1,a1,L1,M1,N1,
			               x2,y2,z2,a2,L2,M2,N2,
			               x3,y3,z3,a3,L3,M3,N3,
			               x4,y4,z4,a4,L4,M4,N4
                           ) );
				}
			}
		}
	}
	return sum;
}


#endif




