#ifndef basic_InteG_H
#define basic_InteG_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include <cmath>		// erf
#include <Eigen/Dense>

#include "basic_Const.h"
#include "basic_Tools.h"
#include "basic_Global.h"
#include "basic_Parser.h"


using namespace std;
using namespace Eigen;

//  declarations of basic s-s type integration of contract gauss fucntion
/////////////////////////////////////////////////////////////////////////////////
	double integral_S_sstype(	 Orbital_cgto& cgto1,  Orbital_cgto& cgto2);
	double integral_ERI_sstype(  Orbital_cgto& cgto1, 
								 Orbital_cgto& cgto2, 
								 Orbital_cgto& cgto3, 
								 Orbital_cgto& cgto4);
	double integral_T_sstype( 	 Orbital_cgto& cgto1,  Orbital_cgto& cgto2);
	double integral_V_sstype( 	 Orbital_cgto& cgto1,  Orbital_cgto& cgto2, 
								Atom& P );
//  declarations of general integration of contract gauss fucntion
//  basic mathematics
	int binomial(int a, int b);
	int factorial(int a);
	int doublefactorial(int a);
	double Hermite_normial(double x, int n);
	double Kcoeff(
		Orbital& orb1,
		Orbital& orb2
	);
	double general_Kcoeff(
		Orbital& orb1,
		Orbital& orb2,
		int d1, int d2, int d3, int d4, int d5, int d6
	);
	double general_Kcoeff_ERI(
		Orbital& orb1,
		Orbital& orb2,
		Orbital& orb3,
		Orbital& orb4,
		int d1, int d2, int d3, int d4, int d5, int d6,
		int d7, int d8, int d9, int d10, int d11, int d12
	);
	double Fn(double w, int m);
	double Fn_pade(double w, int m);
	double general_Fn_V(
		Orbital& orb1,
		Orbital& orb2,
		Atom& PN,
		int d1, int d2, int d3, int d4, int d5, int d6
	);
	double general_Fn_ERI(
		Orbital& orb1,
		Orbital& orb2,
		Orbital& orb3,
		Orbital& orb4,
		int d1, int d2, int d3, int d4, int d5, int d6,
		int d7, int d8, int d9, int d10, int d11, int d12
	);
	double IntecGTO_S(
		Orbital_cgto& cgto1,
		Orbital_cgto& cgto2
	);
	double InteGTO_S(
		Orbital& orb1, 
		Orbital& orb2
	);
	double IScoeff(
		Orbital& orb1, 
		Orbital& orb2,
		char flag
	);
	double IntecGTO_T(
		Orbital_cgto& cgto1,
		Orbital_cgto& cgto2
	);
	double InteGTO_T(
		Orbital& orb1, 
		Orbital& orb2
	);
	double IntecGTO_V(
		Orbital_cgto& cgto1,
		Orbital_cgto& cgto2,
		Atom& PN
	);
	double general_InteGTO_V(
		Orbital& orb1, 
		Orbital& orb2,
		Atom& PN,
		int f1, int f2, int f3, int f4, int f5, int f6,
		int d1, int d2, int d3, int d4, int d5, int d6
	);
	double IntecGTO_ERI(
		Orbital_cgto& cgto1,
		Orbital_cgto& cgto2,
		Orbital_cgto& cgto3,
		Orbital_cgto& cgto4
	);
	double general_InteGTO_ERI(
		Orbital& orb1, 
		Orbital& orb2,
		Orbital& orb3, 
		Orbital& orb4,
		int f1,	int f2, int f3, int f4, int f5, int f6, int f7, int f8, int f9, int f10, int f11, int f12,
		int d1, int d2, int d3, int d4, int d5, int d6, int d7, int d8, int d9, int d10, int d11, int d12
	);


//  definations of basic integrations of s-s type
double integral_S_sstype( Orbital_cgto& cgto1,  Orbital_cgto& cgto2)
{
	double AB2, sum;
	
	AB2 = ( cgto1.get_Point() - cgto2.get_Point() ) * ( cgto1.get_Point() - cgto2.get_Point() );
	sum = 0;
	
	double a1, a2, c1, c2, norm;
	
	for(int i = 0; i < cgto1.cn; i++){
		cgto1.get_CA(i, c1, a1);
		for(int j = 0; j < cgto2.cn; j++){
			cgto2.get_CA(j, c2, a2);
						
			norm = pow(4*a1*a2/(PI*PI), 0.75);
			sum += norm * c1 * c2 * 
			exp(- AB2*a1*a2 / (a1+a2)) *
			pow(PI/(a1+a2),1.5);
			
			/* 	note: 	norm = (4*a1*a2/PI**2)**0.75
						expterm = exp( - a1*a2/(a1+a2) * r12**2 )
						factor = (PI**2/(a1+a2)**2)*0.75
						norm*factor
							= (4*a1*a2/(a1+a2)**2)*0.75		, formulate easily
				but to avoid confusion,
						we use three term: norm, exp-term, factor
			*/
		}
	}
	return sum;
}

double integral_ERI_sstype(
	 Orbital_cgto& cgto1, 
	 Orbital_cgto& cgto2, 
	 Orbital_cgto& cgto3, 
	 Orbital_cgto& cgto4
){
	double eps = 1e-10;
	double sum, F0;
	double AB2, CD2;

	AB2 = ( cgto1.get_Point() - cgto2.get_Point() ) * ( cgto1.get_Point() - cgto2.get_Point() );
	CD2 = ( cgto3.get_Point() - cgto4.get_Point() ) * ( cgto3.get_Point() - cgto4.get_Point() );
	
	sum = 0;
	double a1, a2, a3, a4, w, a12, a34, a1234, c1, c2, c3, c4, PQx, PQy, PQz, PQ2, norm;
	for(int i = 0; i < cgto1.cn; i++){
		cgto1.get_CA(i, c1, a1);
		for(int j = 0; j < cgto2.cn; j++){
			cgto2.get_CA(j, c2, a2);
			a12 = a1 + a2;
			for(int k = 0; k < cgto3.cn; k++){
				cgto3.get_CA(k, c3, a3);
				for(int m = 0; m < cgto4.cn; m++){
					cgto4.get_CA(m, c4, a4);
					a34 = a3 + a4;
					PQx = (a1*cgto1.x + a2*cgto2.x)/a12 - (a3*cgto3.x + a4*cgto4.x)/a34;
					PQy = (a1*cgto1.y + a2*cgto2.y)/a12 - (a3*cgto3.y + a4*cgto4.y)/a34;
					PQz = (a1*cgto1.z + a2*cgto2.z)/a12 - (a3*cgto3.z + a4*cgto4.z)/a34;
					PQ2 = PQx*PQx + PQy*PQy + PQz*PQz;
					
					w = PQ2*a12*a34 / (a12 + a34);
					F0 = (abs(w) < eps ? 1.0 : 0.5 * sqrt(PI / w) * erf(sqrt(w)));
					
					norm = pow(16*a1*a2*a3*a4/(PI*PI*PI*PI), 0.75);
					sum += 2*norm*c1*c2*c3*c4*pow(PI, 2.5) * 
					exp(- a1*a2*AB2/a12 - a3*a4*CD2/a34) * F0/
					(a12 * a34 * sqrt(a12 + a34));
					// now soppose it's correct, sure correct
				}
			}
		}
	}
	return sum;
}

double integral_T_sstype(
	 Orbital_cgto& cgto1, 
	 Orbital_cgto& cgto2 
){
	double AB2, sum;
	
	AB2 = ( cgto1.get_Point() - cgto2.get_Point() ) * ( cgto1.get_Point() - cgto2.get_Point() );
	
	sum = 0;
	double a1, a2, c1, c2, temp, norm;
	
	for(int i = 0; i < cgto1.cn; i++){
		cgto1.get_CA(i, c1, a1);
		for(int j = 0; j < cgto2.cn; j++){
			cgto2.get_CA(j, c2, a2);
			temp = a1 * a2 / (a1 + a2);
			norm = pow(4*a1*a2/(PI*PI), 0.75);
			sum += norm * c1 * c2 * temp * (3 - 2*temp*AB2) * 
			exp(- temp*AB2) * pow(PI / (a1+a2),1.5);
			
			/* 	note: 	
						Hpq = (a1*a2)/(a1+a2) * (3-2*a1*a2/(a1+a2)*r12**2) * Spq
			*/
			
		}
	}
	return sum;
}

double integral_V_sstype(
	 Orbital_cgto& cgto1, 
	 Orbital_cgto& cgto2, 
	 Atom& P
){
	double eps = 1e-10;
	double px, py, pz, AB2, sum, F0;
	
	AB2 = ( cgto1.get_Point() - cgto2.get_Point() ) * ( cgto1.get_Point() - cgto2.get_Point() );
	
	sum = 0;
	double a1, a2, c1, c2, a12, w, norm;
	for(int i = 0; i < cgto1.cn; ++i){
		cgto1.get_CA(i, c1, a1);
		for(int j = 0; j < cgto2.cn; ++j){
			cgto2.get_CA(j, c2, a2);
			a12 = a1 + a2;
			px = (cgto1.x*a1 + cgto2.x*a2) / a12;
			py = (cgto1.y*a1 + cgto2.y*a2) / a12;
			pz = (cgto1.z*a1 + cgto2.z*a2) / a12;			
			w = a12 * 
			((px - P.x)*(px - P.x) + (py - P.y)*(py - P.y) + (pz - P.z)*(pz - P.z));
			F0 = (abs(w) < eps ? 1.0 : 0.5 * sqrt(PI / w) * erf(sqrt(w)));
			norm = pow(4*a1*a2/(PI*PI), 0.75);
			sum += norm * c1 * c2 * 2 * PI *
			exp(- AB2*a1*a2 / a12) * F0 / a12;
			
			/* 	note: 	F0:		sqrt(PI/w) * erf(sqrt(w) )	= erf(...)*sqrt(PI/(a1+a2)) / |_rc|
								w = (a1+a2) * |rc|
								_rc = _r' - _rn				// _A stand A is vector
								_r = (_r1*a1 + _r2*a2) / (a1+a2) 
								
						Vpq(r) 	= -Z * erf(...) * Spq / |_rc|
								= -Z * F0/ sqrt(PI/(a1+a1)) *Spq
								= -Z * F0 * (PI/(a1+a2)) *  [ (PI/(a1+a2))**(-1.5) * Spq ]
			*/
		}
	}
	return sum * (-1) * P.znum;		// -1 note for electron, and znum for the nuclues
}




//////////////////////////////////////////////////////////////////////////////////////////////////
/*

				the general integral function as follow

*/
//////////////////////////////////////////////////////////////////////////////////////////////////
//	needed add some general function
/*	1: basic mathematic functions
	binomial(a,b), return C^b_a
	factorial(a), return a!
	doublefactorial(a), return a!!

*/
int binomial(int a, int b)
{
	if(a<b){cerr << "error" << endl; exit(-1);}
	int result = 1;
	int factor = 1;
	for(int i=1; i<=b; i++){
		result*=(a-i+1);
		factor*=i;
	}
	return result/factor;
}

int factorial(int a)
{
	if(a==0) return 1;
	int result=1;
	for(int i=a; i>1; i--){
		result*=i;
	}
	return result;
}

int doublefactorial(int a)
{
	if(a==0 || a==-1) return 1;
	int result=1;
	for(int i=a; i>=2; i-=2){
		result*=i;
	}
	return result;
}

double Hermite_normial(
	double x,
	int n
)
{
	if(n==0) return 1;
	if(n==1) return 2*x;
	if(n==2) return 4*x*x-2;
	if(n==3) return 8*x*x*x -12*x;
	if(n==4) return 16*pow(x,4)-48*pow(x,2)+12;
	if(n==5) return 32*pow(x,5)-160*pow(x,3)+120*x;
	if(n==6) return 64*pow(x,6)-480*pow(x,4)+720*pow(x,2)-120;
	if(n==7) return 128*pow(x,7)-1344*pow(x,5)+3360*pow(x,3)-1680*x;
	if(n==8) return 256*pow(x,8)-3584*pow(x,6)+13440*pow(x,4)-13440*x*x+1680;
	
	double sum = 0;
	for(int m=0; m<=n/2; m++)
	{
		int hermi = 1;
		int factor = 1;
		for(int j=1; j<=m; j++)
		{
			hermi *= (n-2*j+2)*(n-2*j+1);
			factor *= (-2*j);
		}
		sum += (double)(hermi/factor) * pow(x, n-2*m);
	}
	return sum;
}

double Kcoeff(
	Orbital& orb1, 
	Orbital& orb2
)
{
	double AB2 = (orb1.get_Point() - orb2.get_Point()).norm2() ;
	double effalpha = (orb1.alpha*orb2.alpha) / (orb1.alpha+orb2.alpha);
	return exp( - effalpha * AB2 );
}

double general_Kcoeff
(
	Orbital& orb1,
	Orbital& orb2,
	int d1, int d2, int d3, int d4, int d5, int d6
)
{
	int LMN = orb1.L + orb1.M + orb1.N;
	int LMNLMN = orb1.L + orb1.M + orb1.N + orb2.L + orb2.M + orb2.N;
	int LL = orb1.L + orb2.L;
	int MM = orb1.M + orb2.M;
	int NN = orb1.N + orb2.N;
	double effalpha = orb1.alpha*orb2.alpha/(orb1.alpha+orb2.alpha);
	double effABx = sqrt(effalpha) * (orb1.get_Point() - orb2.get_Point()).x;
	double effABy = sqrt(effalpha) * (orb1.get_Point() - orb2.get_Point()).y;
	double effABz = sqrt(effalpha) * (orb1.get_Point() - orb2.get_Point()).z;
	
	return Kcoeff(orb1, orb2) * (LMN%2==0?1:-1) * pow(effalpha,0.5*LMNLMN) * 
	Hermite_normial(effABx,LL) * Hermite_normial(effABy,MM) * Hermite_normial(effABz,NN);
}

double general_Kcoeff_ERI(
	Orbital& orb1,
	Orbital& orb2,
	Orbital& orb3,
	Orbital& orb4,
	int d1, int d2, int d3, int d4, int d5, int d6,
	int d7, int d8, int d9, int d10, int d11, int d12
)
{
	return general_Kcoeff(orb1, orb2, d1, d2, d3, d4, d5, d6) * 
	general_Kcoeff(orb3, orb4, d7, d8, d9, d10, d11, d12);
}


double Fn(
	double w,
	int m
)
{
	double sum=0;
	for(int i=0; i<20; i++)
	{
		int factor = 2*m+1;
		for(int k=1; k<=i; k++)
		{
			factor*=(2*m+1+2*k);
		}
		sum+=pow(2*w,i)/((double)factor);
	}
	return exp(-w)*sum;
}


double general_Fn_V(
	Orbital& orb1,
	Orbital& orb2,
	Atom& PN,
	int d1, int d2, int d3, int d4, int d5, int d6
)
{
	int LMN = orb1.L + orb1.M + orb1.N;
	int LMNLMN = orb1.L + orb1.M + orb1.N + orb2.L + orb2.M + orb2.N;
	int LL = orb1.L + orb2.L;
	int MM = orb1.M + orb2.M;
	int NN = orb1.N + orb2.N;
	double effalpha = orb1.alpha + orb2.alpha;
	Point PC = (orb1.get_Point()*orb1.alpha + orb2.get_Point()*orb2.alpha) / 
			(orb1.alpha + orb2.alpha) - PN.get_Point();
	double w = effalpha * PC.norm2();
	double effPCx = sqrt(effalpha) * PC.x;
	double effPCy = sqrt(effalpha) * PC.y;
	double effPCz = sqrt(effalpha) * PC.z;
	
	return pow(orb1.alpha/effalpha, LMN) * pow(orb2.alpha/effalpha, LMNLMN-LMN) * 
	(LMNLMN%2==0?1:-1) * pow(effalpha, 0.5*LMNLMN) * Fn_pade(w, LMNLMN) *
	Hermite_normial(effPCx,LL) * Hermite_normial(effPCy,MM) * Hermite_normial(effPCz,NN); 
}


double general_Fn_ERI(
	Orbital& orb1,
	Orbital& orb2,
	Orbital& orb3,
	Orbital& orb4,
	int d1, int d2, int d3, int d4, int d5, int d6,
	int d7, int d8, int d9, int d10, int d11, int d12
)
{
	int LMN1 = orb1.L + orb1.M + orb1.N;
	int LMN2 = orb2.L + orb2.M + orb2.N;
	int LMN3 = orb3.L + orb3.M + orb3.N;
	int LMN4 = orb4.L + orb4.M + orb4.N;
	int LMN12 = LMN1+LMN2;
	int LMN34 = LMN3+LMN4;
	int LMN1234 = LMN12 + LMN34;
	int LL = orb1.L + orb2.L + orb3.L + orb4.L;
	int MM = orb1.M + orb2.M + orb3.M + orb4.M;
	int NN = orb1.N + orb2.N + orb3.N + orb4.N;
	double effalpha = (orb1.alpha + orb2.alpha) * (orb3.alpha + orb4.alpha)/
					(orb1.alpha + orb2.alpha + orb3.alpha + orb4.alpha);
	Point P = (orb1.get_Point() * orb1.alpha + orb2.get_Point() * orb2.alpha) / (orb1.alpha + orb2.alpha);
	Point Q = (orb3.get_Point() * orb3.alpha + orb4.get_Point() * orb4.alpha) / (orb3.alpha + orb4.alpha);
	double w = effalpha * (P-Q).norm2();
	double effPQx = sqrt(effalpha) * (P-Q).x;
	double effPQy = sqrt(effalpha) * (P-Q).y;
	double effPQz = sqrt(effalpha) * (P-Q).z;
	
	return pow(orb1.alpha/(orb1.alpha+orb2.alpha), LMN1) * pow(orb2.alpha/(orb1.alpha+orb2.alpha), LMN2) * 
	pow(orb3.alpha/(orb3.alpha+orb4.alpha), LMN3) * pow(orb4.alpha/(orb3.alpha+orb4.alpha), LMN4) * 
	(LMN12%2==0?1:-1) * pow(effalpha, 0.5*LMN1234) * Fn_pade(w, LMN1234) *
	Hermite_normial(effPQx,LL) * Hermite_normial(effPQy,MM) * Hermite_normial(effPQz,NN); 
}

/////////////////////////////////////////////////// general integral function of S matrix
double IntecGTO_S(
	Orbital_cgto& cgto1,
	Orbital_cgto& cgto2
)
{
	double sum=0;
	for(int i = 0; i < cgto1.cn; i++){
		for(int j = 0; j < cgto2.cn; j++){
			sum += cgto1.coeffs[i] * cgto2.coeffs[j] * cgto1.normGTO(i) * cgto2.normGTO(j)*
			InteGTO_S( cgto1.gtos[i], cgto2.gtos[j] );
		}
	}
	return sum;
}


double InteGTO_S(
	Orbital& orb1, 
	Orbital& orb2
)
{
	if(orb1.L<0 || orb2.L <0 || orb1.M<0 || orb2.M <0 || orb1.N<0 || orb2.N <0) return 0;
	return Kcoeff(orb1, orb2) * IScoeff(orb1, orb2, 'x') * IScoeff(orb1, orb2, 'y') * IScoeff(orb1, orb2, 'z');
}

double IScoeff(
	Orbital& orb1, 
	Orbital& orb2,
	char flag
)
{
	double I = 0;
	int L1, L2;
	Point P = (orb1.get_Point()*orb1.alpha + orb2.get_Point()*orb2.alpha) / (orb1.alpha + orb2.alpha) ;
	double PA;
	double PB;
	
	if(flag=='x'){
		L1 = orb1.L; L2=orb2.L;
		PA = (orb1.get_Point() - P).x;
		PB = (orb2.get_Point() - P).x;
	}else if(flag=='y'){
		L1 = orb1.M; L2=orb2.M;
		PA = (orb1.get_Point() - P).y;
		PB = (orb2.get_Point() - P).y;
	}else if(flag=='z'){
		L1 = orb1.N; L2=orb2.N;
		PA = (orb1.get_Point() - P).z;
		PB = (orb2.get_Point() - P).z;
	}else{
		cerr << "error, flag of I of Integral of S is incorrect" << endl; exit(-1);
	}
	
	for(int i=0; i<=(L1+L2)/2; i++)
	{
		double temp=0;
		for(int j=0; j<=L1; j++)
		{
			for(int k=0; k<=L2; k++)
			{
				if(j+k==2*i)
				{
					temp+=binomial(L1,j)*binomial(L2,j)*pow(PA, L1-j)*pow(PB, L2-k);
				}
			}
		}
		I+=temp*doublefactorial(2*i-1) / pow(2*(orb1.alpha+orb2.alpha), i);
	}
	return sqrt(PI/(orb1.alpha+orb2.alpha)) * I;
}



///////////////////////////////////////////////////////////////////////////////////////
double IntecGTO_T(
	Orbital_cgto& cgto1,
	Orbital_cgto& cgto2
)
{
	double sum=0;
	for(int i = 0; i < cgto1.cn; i++){
		for(int j = 0; j < cgto2.cn; j++){
			sum += cgto1.coeffs[i] * cgto2.coeffs[j] * cgto1.normGTO(i) * cgto2.normGTO(j)*
			InteGTO_T( cgto1.gtos[i], cgto2.gtos[j] );
		}
	}
	return sum;
}

double InteGTO_T(
	Orbital& orb1, 
	Orbital& orb2
)
{
	Orbital orb3 = Orbital();
	double result=0;
	
	
	orb3.set_XYZ(orb2.x, orb2.y, orb2.z);
	orb3.set_Alpha(orb2.alpha);
	
	// x direction
	orb3.set_LMN(orb2.L-2, orb2.M, orb2.N);
	result += -(orb2.L*(orb2.L+1))/2 * InteGTO_S(orb1, orb3);
	result += orb2.alpha*(2*orb2.L+1) * InteGTO_S(orb1, orb2);
	orb3.set_LMN(orb2.L+2, orb2.M, orb2.N);
	result += - 2* (orb2.alpha*orb2.alpha) * InteGTO_S(orb1, orb3);
	
	// y direction
	orb3.set_LMN(orb2.L, orb2.M-2, orb2.N);
	result += -(orb2.M*(orb2.M+1))/2 * InteGTO_S(orb1, orb3);
	result += orb2.alpha*(2*orb2.M+1) * InteGTO_S(orb1, orb2);
	orb3.set_LMN(orb2.L, orb2.M+2, orb2.N);
	result += - 2* (orb2.alpha*orb2.alpha) * InteGTO_S(orb1, orb3);
	
	// z direction
	orb3.set_LMN(orb2.L, orb2.M, orb2.N-2);
	result += -(orb2.N*(orb2.N+1))/2 * InteGTO_S(orb1, orb3);
	result += orb2.alpha*(2*orb2.N+1) * InteGTO_S(orb1, orb2);
	orb3.set_LMN(orb2.L, orb2.M, orb2.N+2);
	result += - 2* (orb2.alpha*orb2.alpha) * InteGTO_S(orb1, orb3);
	
	return result;
}


//////////////////////////////////////////////////////////////////////////////////
double IntecGTO_V(
	Orbital_cgto& cgto1,
	Orbital_cgto& cgto2,
	Atom& PN
)
{
	double sum=0;
	for(int i = 0; i < cgto1.cn; i++){
		for(int j = 0; j < cgto2.cn; j++){
			sum += cgto1.coeffs[i] * cgto2.coeffs[j] * cgto1.normGTO(i) * cgto2.normGTO(j)*
			general_InteGTO_V( cgto1.gtos[i], cgto2.gtos[j], PN, cgto1.L, cgto1.M, cgto1.N,
			cgto2.L, cgto2.M, cgto2.N, 0,0,0,0,0,0);
		}
	}
	return sum * (-1) * PN.znum;
}

double general_InteGTO_V(
	Orbital& orb1, 
	Orbital& orb2,
	Atom& PN,
	int f1, int f2, int f3, int f4, int f5, int f6,
	int d1, int d2, int d3, int d4, int d5, int d6
)
{
	if(f1<0 || f2<0 || f3<0 || f4<0 || f5<0 || f6<0 ) return 0;
	if(f1==0 && f2==0 && f3==0 && f4==0 && f5==0 && f6==0 &&
	d1==0 && d2==0 && d3==0 && d4==0 && d5==0 && d6==0)
	{
		return 2*PI/(orb1.alpha + orb2.alpha)*Kcoeff(orb1,orb2) * general_Fn_V(orb1, orb2, PN, 0,0,0,0,0,0);
	}
	
	if(f1>0)
	{
		return 0.5/(orb1.alpha)*
		( general_InteGTO_V(orb1, orb2, PN, f1-1, f2, f3, f4, f5, f6, d1+1, d2, d3, d4, d5, d6) +
		(f1-1)*general_InteGTO_V(orb1, orb2, PN, f1-2, f2, f3, f4, f5, f6, d1, d2, d3, d4, d5, d6) );
	}
	if(f2>0)
	{
		return  0.5/(orb1.alpha)*
		( general_InteGTO_V(orb1, orb2, PN, f1, f2-1, f3, f4, f5, f6, d1, d2+1, d3, d4, d5, d6) +
		(f2-1)*general_InteGTO_V(orb1, orb2, PN, f1, f2-2, f3, f4, f5, f6, d1, d2, d3, d4, d5, d6) );
	}
	if(f3>0)
	{
		return  0.5/(orb1.alpha)*
		( general_InteGTO_V(orb1, orb2, PN, f1, f2, f3-1, f4, f5, f6, d1, d2, d3+1, d4, d5, d6) +
		(f3-1)*general_InteGTO_V(orb1, orb2, PN, f1, f2, f3-2, f4, f5, f6, d1, d2, d3, d4, d5, d6) );
	}
	if(f4>0)
	{
		return 0.5/(orb2.alpha)*
		( general_InteGTO_V(orb1, orb2, PN, f1, f2, f3, f4-1, f5, f6, d1, d2, d3, d4+1, d5, d6) +
		(f4-1)*general_InteGTO_V(orb1, orb2, PN, f1, f2, f3, f4-2, f5, f6, d1, d2, d3, d4, d5, d6) );
	}
	if(f5>0)
	{
		return  0.5/(orb2.alpha)*
		( general_InteGTO_V(orb1, orb2, PN, f1, f2, f3, f4, f5-1, f6, d1, d2, d3, d4, d5+1, d6) +
		(f5-1)*general_InteGTO_V(orb1, orb2, PN, f1, f2, f3, f4, f5-2, f6, d1, d2, d3, d4, d5, d6) );
	}
	if(f6>0)
	{
		return  0.5/(orb2.alpha)*
		( general_InteGTO_V(orb1, orb2, PN, f1, f2, f3, f4, f5, f6-1, d1, d2, d3, d4, d5, d6+1) +
		(f6-1)*general_InteGTO_V(orb1, orb2, PN, f1, f2, f3, f4, f5, f6-2, d1, d2, d3, d4, d5, d6) );
	}
	// else return f1=f2=f3=f4=f5=f6=0, only d1~d6 need to be dealed!
	double sum=0;
	for(int i1=0; i1<=d1; i1++)
	{
		for(int i2=0;i2<=d2;i2++)
		{
			for(int i3=0; i3<=d3; i3++)
			{
				for(int i4=0; i4<=d4; i4++)
				{
					for(int i5=0; i5<=d5; i5++)
					{
						for(int i6=0; i6<=d6; i6++)
						{
							sum+= binomial(d1,i1)* binomial(d2,i2)* binomial(d3,i3)*
							binomial(d4,i4)* binomial(d5,i5)* binomial(d6,i6)*
							general_Kcoeff(orb1, orb2, i1, i2, i3, i4, i5, i6)*
							general_Fn_V(orb1, orb2, PN, d1-i1, d2-i2, d3-i3, d4-i4, d5-i5, d6-i6);
						}
					}
				}
			}
		}
	}
	sum=2*PI/(orb1.alpha + orb2.alpha)*sum;
	return sum;
}


//////////////////////////////////////////////////////////////////////////////////////////
double IntecGTO_ERI(
	Orbital_cgto& cgto1,
	Orbital_cgto& cgto2,
	Orbital_cgto& cgto3,
	Orbital_cgto& cgto4
)
{
	double sum=0;
	for(int i = 0; i < cgto1.cn; i++)
	{
		for(int j = 0; j < cgto2.cn; j++)
		{
			for(int k=0; k< cgto3.cn; k++)
			{
				for(int m=0; m < cgto4.cn; m++)
				{
					sum += cgto1.coeffs[i] * cgto2.coeffs[j] * cgto3.coeffs[k] * cgto4.coeffs[m] 
					* cgto1.normGTO(i) * cgto2.normGTO(j) * cgto3.normGTO(k) * cgto4.normGTO(m) *
					general_InteGTO_ERI( cgto1.gtos[i], cgto2.gtos[j], cgto3.gtos[k], cgto4.gtos[m],
					cgto1.gtos[i].L, cgto1.gtos[i].M, cgto1.gtos[i].N,
					cgto2.gtos[j].L, cgto2.gtos[j].M, cgto2.gtos[j].N,
					cgto3.gtos[k].L, cgto3.gtos[k].M, cgto3.gtos[k].N,
					cgto4.gtos[m].L, cgto4.gtos[m].M, cgto4.gtos[m].N,
					0,0,0,0,0,0,0,0,0,0,0,0 );
				}
			}
		}
	}
	return sum;
}

double general_InteGTO_ERI(
	Orbital& orb1, 
	Orbital& orb2,
	Orbital& orb3, 
	Orbital& orb4,
	int f1,	int f2, int f3, int f4, int f5, int f6, int f7, int f8, int f9, int f10, int f11, int f12,
	int d1, int d2, int d3, int d4, int d5, int d6, int d7, int d8, int d9, int d10, int d11, int d12
)
{
	if(f1<0 || f2<0 || f3<0 || f4<0 || f5<0 || f6<0 || f7<0 || f8<0 || f9<0 || f10<0 || f11<0 || f12<0 ) return 0;
	if(f1==0 && f2==0 && f3==0 && f4==0 && f5==0 && f6==0 && f7<0 && f8<0 && f9<0 && f10<0 && f11<0 && f12<0 && 
	d1==0 && d2==0 && d3==0 && d4==0 && d5==0 && d6==0 && d7==0 && d8==0 && d9==0 && d10==0 && d11==0 && d12==0)
	{
		double temp = (orb1.alpha+orb2.alpha) * (orb3.alpha+orb4.alpha) *
					sqrt(orb1.alpha+orb2.alpha+orb3.alpha+orb4.alpha);
		return 2*pow(PI,2.5) / temp *Kcoeff(orb1,orb2) * Kcoeff(orb3,orb4) * 
		general_Fn_ERI(orb1, orb2, orb3, orb4, 0,0,0,0,0,0,0,0,0,0,0,0);
	}
	
	if(f1>0)
	{
		return 0.5/(orb1.alpha)*
		( general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1-1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12,
		d1+1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12) +
		(f1-1)*general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1-2, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12,
		d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12) );
	}
	if(f2>0)
	{
		return 0.5/(orb1.alpha)*
		( general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1, f2-1, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12,
		d1, d2+1, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12) +
		(f2-1)*general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1, f2-2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12,
		d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12) );
	}
	if(f3>0)
	{
		return 0.5/(orb1.alpha)*
		( general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1, f2, f3-1, f4, f5, f6, f7, f8, f9, f10, f11, f12,
		d1, d2, d3+1, d4, d5, d6, d7, d8, d9, d10, d11, d12) +
		(f3-1)*general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1, f2, f3-2, f4, f5, f6, f7, f8, f9, f10, f11, f12,
		d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12) );
	}
	if(f4>0)
	{
		return 0.5/(orb2.alpha)*
		( general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1, f2, f3, f4-1, f5, f6, f7, f8, f9, f10, f11, f12,
		d1, d2, d3, d4+1, d5, d6, d7, d8, d9, d10, d11, d12) +
		(f4-1)*general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1, f2, f3, f4-2, f5, f6, f7, f8, f9, f10, f11, f12,
		d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12) );
	}
	if(f5>0)
	{
		return 0.5/(orb2.alpha)*
		( general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1, f2, f3, f4, f5-1, f6, f7, f8, f9, f10, f11, f12,
		d1, d2, d3, d4, d5+1, d6, d7, d8, d9, d10, d11, d12) +
		(f5-1)*general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1, f2, f3, f4, f5-2, f6, f7, f8, f9, f10, f11, f12,
		d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12) );
	}
	if(f6>0)
	{
		return 0.5/(orb2.alpha)*
		( general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1, f2, f3, f4, f5, f6-1, f7, f8, f9, f10, f11, f12,
		d1, d2, d3, d4, d5, d6+1, d7, d8, d9, d10, d11, d12) +
		(f6-1)*general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1, f2, f3, f4, f5, f6-2, f7, f8, f9, f10, f11, f12,
		d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12) );
	}
	if(f7>0)
	{
		return 0.5/(orb3.alpha)*
		( general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1, f2, f3, f4, f5, f6, f7-1, f8, f9, f10, f11, f12,
		d1, d2, d3, d4, d5, d6, d7+1, d8, d9, d10, d11, d12) +
		(f7-1)*general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1, f2, f3, f4, f5, f6, f7-2, f8, f9, f10, f11, f12,
		d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12) );
	}
	if(f8>0)
	{
		return 0.5/(orb3.alpha)*
		( general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1, f2, f3, f4, f5, f6, f7, f8-1, f9, f10, f11, f12,
		d1, d2, d3, d4, d5, d6, d7, d8+1, d9, d10, d11, d12) +
		(f8-1)*general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1, f2, f3, f4, f5, f6, f7, f8-2, f9, f10, f11, f12,
		d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12) );
	}
	if(f9>0)
	{
		return 0.5/(orb3.alpha)*
		( general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1, f2, f3, f4, f5, f6, f7, f8, f9-1, f10, f11, f12,
		d1, d2, d3, d4, d5, d6, d7, d8, d9+1, d10, d11, d12) +
		(f9-1)*general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1, f2, f3, f4, f5, f6, f7, f8, f9-2, f10, f11, f12,
		d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12) );
	}
	if(f10>0)
	{
		return 0.5/(orb4.alpha)*
		( general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10-1, f11, f12,
		d1, d2, d3, d4, d5, d6, d7, d8, d9, d10+1, d11, d12) +
		(f10-1)*general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10-2, f11, f12,
		d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12) );
	}
	if(f11>0)
	{
		return 0.5/(orb4.alpha)*
		( general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11-1, f12,
		d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11+1, d12) +
		(f11-1)*general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11-2, f12,
		d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12) );
	}
	if(f12>0)
	{
		return 0.5/(orb4.alpha)*
		( general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12-1,
		d1, d2, d3, d4, d5, d6, d7, d8, d9, d10+1, d11, d12) +
		(f12-1)*general_InteGTO_ERI(orb1, orb2, orb3, orb4, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12-2,
		d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12) );
	}
	// else return f1=f2=f3=f4=f5=f6=0, only d1~d6 need to be dealed!
	double sum=0;
	for(int i1=0; i1<=d1; i1++)
	{
		for(int i2=0;i2<=d2;i2++)
		{
			for(int i3=0; i3<=d3; i3++)
			{
				for(int i4=0; i4<=d4; i4++)
				{
					for(int i5=0; i5<=d5; i5++)
					{
						for(int i6=0; i6<=d6; i6++)
						{
							for(int i7=0; i7<=d7; i7++)
							{
								for(int i8=0; i8<=d8; i8++)
								{
									for(int i9=0; i9<=d9; i9++)
									{
										for(int i10=0; i10<=d10; i10++)
										{
											for(int i11=0; i11<=d11; i11++)
											{
												for(int i12=0; i12<=d12; i12++)
												{
													sum+= binomial(d1,i1)* binomial(d2,i2)* binomial(d3,i3)*
													binomial(d4,i4)* binomial(d5,i5)* binomial(d6,i6)*
													binomial(d7,i7)* binomial(d8,i8)* binomial(d9,i9)*
													binomial(d10,i10)* binomial(d11,i11)* binomial(d12,i12)*
													general_Kcoeff(orb1, orb2, d1, d2, d3, d4, d5, d6) * 
													general_Kcoeff(orb3, orb4, d7, d8, d9, d10, d11, d12) *
													general_Fn_ERI(orb1, orb2, orb3, orb4, d1-i1, d2-i2, d3-i3, 
													d4-i4, d5-i5, d6-i6, d7-i7, d8-i8, d9-i9, d10-i10, d11-i11, 
													d12-i12);
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	sum *= 2*pow(PI,2.5) / ( (orb1.alpha+orb2.alpha) * (orb3.alpha+orb4.alpha) *
				sqrt(orb1.alpha+orb2.alpha+orb3.alpha+orb4.alpha) );
	return sum;
}

/////////////////////////////////////////////////////////////////////////////////////////

double Fn_pade(double w, int m)
{
	double a1, a2, a3, a4, a5, b1, b2, b3, b4, b5, b6;
	double a0 = pow( (double)(2*m+1), -2/((double)(2*m+1)) );
	if(m==0){
		if(w>16.3578)
		{
			return (doublefactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2);
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
			return (doublefactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2);
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
			return (doublefactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2);
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
			return (doublefactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2);
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
			return (doublefactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
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
			return (doublefactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
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
			return (doublefactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
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
			return (doublefactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
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
			return (doublefactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
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
			return (doublefactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
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
			return (doublefactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
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
			return (doublefactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
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
			return (doublefactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
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
			return (doublefactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
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
			return (doublefactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
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
			return (doublefactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
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
			return (doublefactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) ;
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
	return ( (2*m - 1) * Fn_pade(w,m-1) - exp(-w) ) / (2*w);
}


#endif




