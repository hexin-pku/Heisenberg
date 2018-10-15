#ifndef Hsbg_InteG_H
#define Hsbg_InteG_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include <cmath>		// erf
#include <Eigen/Dense>

#include "Hsbg_Const.h"
#include "Hsbg_Tools.h"
#include "Hsbg_Global.h"


using namespace std;
using namespace Hsbg;
using namespace Eigen;


//  Integral of contract gauss fucntion

	double integral_S(	 Orbital_cgto& cgto1,  Orbital_cgto& cgto2);
	double integral_ERI(  Orbital_cgto& cgto1, 
								 Orbital_cgto& cgto2, 
								 Orbital_cgto& cgto3, 
								 Orbital_cgto& cgto4);
	double integral_T( 	 Orbital_cgto& cgto1,  Orbital_cgto& cgto2);
	double integral_V( 	 Orbital_cgto& cgto1,  Orbital_cgto& cgto2, 
								Atom& P );
	
	// s-s integral type
	double integral_S_sstype(	 Orbital_cgto& cgto1,  Orbital_cgto& cgto2);
	double integral_ERI_sstype(  Orbital_cgto& cgto1, 
								 Orbital_cgto& cgto2, 
								 Orbital_cgto& cgto3, 
								 Orbital_cgto& cgto4);
	double integral_T_sstype( 	 Orbital_cgto& cgto1,  Orbital_cgto& cgto2);
	double integral_V_sstype( 	 Orbital_cgto& cgto1,  Orbital_cgto& cgto2, 
								Atom& P );
	// s-p integral type
	double integral_S_sptype(	 Orbital_cgto& cgto1,  Orbital_cgto& cgto2){return 0;};
	double integral_ERI_sptype(  Orbital_cgto& cgto1, 
								 Orbital_cgto& cgto2, 
								 Orbital_cgto& cgto3, 
								 Orbital_cgto& cgto4){return 0;};
	double integral_T_sptype( 	 Orbital_cgto& cgto1,  Orbital_cgto& cgto2){return 0;};
	double integral_V_sptype( 	 Orbital_cgto& cgto1,  Orbital_cgto& cgto2, 
								Atom& P ){return 0;};
	// p-p integral type
	double integral_S_pptype(	 Orbital_cgto& cgto1,  Orbital_cgto& cgto2){return 0;};
	double integral_ERI_pptype(  Orbital_cgto& cgto1, 
								 Orbital_cgto& cgto2, 
								 Orbital_cgto& cgto3, 
								 Orbital_cgto& cgto4){return 0;};
	double integral_T_pptype( 	 Orbital_cgto& cgto1,  Orbital_cgto& cgto2){return 0;};
	double integral_V_pptype( 	 Orbital_cgto& cgto1,  Orbital_cgto& cgto2, 
								Atom& P ){return 0;};
	
	double integral_S(	 Orbital_cgto& cgto1,  Orbital_cgto& cgto2){
		if(cgto1.get_BohrL()+cgto2.get_BohrL() == 0 ){
			return integral_S_sstype(cgto1, cgto2);
		}else if( cgto1.get_BohrL()==0 && cgto2.get_BohrL() == 1 )
		{		
			return integral_S_sptype(cgto1, cgto2);
		}else if( cgto1.get_BohrL()==1 && cgto2.get_BohrL() == 0 )
		{		
			return integral_S_sptype(cgto2, cgto1);
		}
		else if(cgto1.get_BohrL()+cgto1.get_BohrL() == 2)
		{
			return integral_S_pptype(cgto1, cgto2);
		}
		cerr << "integral_S parser error" << endl;
		return 0;
	}
	
	double integral_ERI(  Orbital_cgto& cgto1, 
								 Orbital_cgto& cgto2, 
								 Orbital_cgto& cgto3, 
								 Orbital_cgto& cgto4);
	
	double integral_T( 	 Orbital_cgto& cgto1,  Orbital_cgto& cgto2)
	{
		if(cgto1.get_BohrL()+cgto2.get_BohrL() == 0 ){
			return integral_T_sstype(cgto1, cgto2);
		}else if( cgto1.get_BohrL()==0 && cgto2.get_BohrL() == 1 )
		{		
			return integral_T_sptype(cgto1, cgto2);
		}else if( cgto1.get_BohrL()==1 && cgto2.get_BohrL() == 0 )
		{		
			return integral_T_sptype(cgto2, cgto1);
		}
		else if(cgto1.get_BohrL()+cgto1.get_BohrL() == 2)
		{
			return integral_T_pptype(cgto1, cgto2);
		}
		cerr << "integral_S parser error" << endl;
		return 0;
	}
	
	double integral_V( 	 Orbital_cgto& cgto1,  Orbital_cgto& cgto2, 
								Atom& P )
	{
		if(cgto1.get_BohrL() + cgto2.get_BohrL() == 0 ){
			return integral_V_sstype(cgto1, cgto2, P);
		}else if( cgto1.get_BohrL()==0 && cgto2.get_BohrL() == 1 )
		{		
			return integral_V_sptype(cgto1, cgto2, P);
		}else if( cgto1.get_BohrL()==1 && cgto2.get_BohrL() == 0 )
		{		
			return integral_V_sptype(cgto2, cgto1, P);
		}
		else if(cgto1.get_BohrL()+cgto1.get_BohrL() == 2)
		{
			return integral_V_pptype(cgto1, cgto2, P);
		}
		cerr << "integral_S parser error" << endl;
		return 0;
	}
								

	
/////////////////////////////////////////////////////////////////////////////////////////////////	

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


//	needed add some general function
/*
double Inte_S(
	Orbital& orb1, 
	Orbital& orb2
)
{
	if(orb1.L<0 || orb2.L <1 || orb1.M<0 || orb2.M <0 || orb1.N<0 || orb2.N <0) return 0;
	return Kcoeff(orb1, orb2) * IScoeff(orb1, orb2, 'x') * IScoeff(orb1, orb2, 'y') * IScoeff(orb1, orb2, 'z');
}

double Kcoeff(
	Orbital& orb1, 
	Orbital& orb2
)
{
	double AB2 = (orb1.get_Point() - orb.get_Point()).norm2 ;
	return exp( -(orb1.alpha*orb2.alpha) / (orb1.alpha+orb2.alpha) * AB2 );
}

double IScoeff(
	Orbital& orb1, 
	Orbital& orb2,
	char flag
)
{
	double I = 0;
	int L1, L2;
	Point refPoint = (orb1.get_Point() - orb.get_Point()) ;
	double PA;
	double PB;
	
	if(flag='x'){
		L1 = orb1.L; L2=orb2.L;
		PA = orb2.alpha/(orb1.alpha+orb2.alpha) * abs( refPoint.x );
		PB = orb1.alpha/(orb1.alpha+orb2.alpha) * abs( refPoint.x );
	}else if(flag='y'){
		L1 = orb1.M; L2=orb2.M;
		PA = orb2.alpha/(orb1.alpha+orb2.alpha) * abs( refPoint.y );
		PB = orb1.alpha/(orb1.alpha+orb2.alpha) * abs( refPoint.y );
	}else if(flag='z'){
		L1 = orb1.N; L2=orb2.N;
		PA = orb2.alpha/(orb1.alpha+orb2.alpha) * abs( refPoint.z );
		PB = orb1.alpha/(orb1.alpha+orb2.alpha) * abs( refPoint.z );
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
				if(j+k==2i){
					temp+=Bicoeff(L1,j)*Bicoeff(L2,j)*pow(PA, L1-j)*pow(PB, L2-k);
				}else{
					continue;
				}
			}
		}
		I+=temp*dfactorial(2i-1) * pow(2*(orb1.alpha+orb2.alpha), i);
	}
	return sqrt(PI/(orb1.alpha+orb2.alpha)) * I;
}

int Bicoeff(int a, int b)
{
	if(a<b){cerr << "error" << endl; exit(-1);}
	int result = a;
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

int dfactorial(int a)
{
	if(a==0) return 1;
	int result=1;
	for(int i=a; i>=2; i-=2){
		result*=i;
	}
	return result;
}

double Inte_T(
	Orbital& orb1, 
	Orbital& orb2,
){
	HOribtal orb3 = Orbital();
	double result=0;
	
	
	orb3.set_XYZ(orb2.x, orb2.y, orb2.z);
	orb3.set_Alpha(orb2.alpha);
	
	// x direction
	orb3.set_LMN(orb2.L-2, orb2.M, orb2.N);
	result += -(orb2.L*(orb2.L+1))/2 * Inte_S(orb1, orb3);
	result += orb2.alpha*(2*orb2.L+1) * Inte_S(orb1, orb1);
	orb3.set_LMN(orb2.L+2, orb2.M, orb2.N);
	result += - 2* (orb2.alpha*orb2.alpha) * Inte_S(orb1, orb3);
	
	// y direction
	orb3.set_LMN(orb2.L, orb2.M-2, orb2.N);
	result += -(orb2.M*(orb2.M+1))/2 * Inte_S(orb1, orb3);
	result += orb2.alpha*(2*orb2.M+1) * Inte_S(orb1, orb1);
	orb3.set_LMN(orb2.L, orb2.M+2, orb2.N);
	result += - 2* (orb2.alpha*orb2.alpha) * Inte_S(orb1, orb3);
	
	// z durection
	orb3.set_LMN(orb2.L, orb2.M, orb2.N-2);
	result += -(orb2.N*(orb2.N+1))/2 * Inte_S(orb1, orb3);
	result += orb2.alpha*(2*orb2.N+1) * Inte_S(orb1, orb1);
	orb3.set_LMN(orb2.L, orb2.M, orb2.N+2);
	result += - 2* (orb2.alpha*orb2.alpha) * Inte_S(orb1, orb3);
	
	return result;
}

double Inte_V(
	Orbital& orb1, 
	Orbital& orb2,
	PBoint& PB
)
{
	int L1 = orb1.L; int L2 = orb2.L;
	int M1 = orb1.M; int M2 = orb2.M;
	int N1 = orb1.N; int N2 = orb2.N;
	double sum = 0;
	for(int L1=orb1.L; L1>=0; L1--)
	{
		for(int L2=orb1.L; L1>=0; L1--)
		{
			for(int L1=orb1.L; L1>=0; L1--)
			{
				
			}	
		}
	}
}

double Inte_V_s(
	Orbital& orb1, 
	Orbital& orb2,
	PBoint& PB
)
{
	double w = 
	return Kcoeff(orb1,orb2)*2*PI/(orb1.alpha+orb2.alapha) * Fpade(w,0);
}







double Fpade(double w, int m)
{
	double a1, a2, a3, a4, a5, b1, b2, b3, b4, b5, b6;
	double a0 = pow( (double)(2*m+1), -2/((double)(2*m+1)) );
	if(m==0){
		if(w>16.3578)
		{
			return (bifactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) 
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
			return pow(temp, m+0.5)
		}
	}
	if(m==1){
		if(w>17.4646)
		{
			return (bifactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) 
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
			return pow(temp, m+0.5)
		}
	}
	if(m==2){
		if(w>15.2368)
		{
			return (bifactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) 
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
			return pow(temp, m+0.5)
		}
	}
	if(m==3){
		if(w>16.0419)
		{
			return (bifactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) 
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
			return pow(temp, m+0.5)
		}
	}
	if(m==4){
		if(w>16.8955)
		{
			return (bifactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) 
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
			return pow(temp, m+0.5)
		}
	}
	if(m==5){
		if(w>17.7822)
		{
			return (bifactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) 
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
			return pow(temp, m+0.5)
		}
	}
	if(m==6){
		if(w>15.8077)
		{
			return (bifactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) 
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
			return pow(temp, m+0.5)
		}
	}
	if(m==7){
		if(w>16.5903)
		{
			return (bifactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) 
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
			return pow(temp, m+0.5)
		}
	}
	if(m==8){
		if(w>17.3336)
		{
			return (bifactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) 
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
			return pow(temp, m+0.5)
		}
	}
	if(m==9){
		if(w>15.6602)
		{
			return (bifactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) 
		}else{
			a1=-0.334843993901400E-1;
			a2=0.846637494147059E-3;
			b1=0.495875606944471E-1;
			b2=0.946642302340943E-3;
			b3=0.108367772249790E-4;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) ) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) );
			return pow(temp, m+0.5)
		}
	}
	if(m==10){
		if(w>16.5258)
		{
			return (bifactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) 
		}else{
			a1=-0.335292171805959E-1;
			a2=0.749168957990503E-3;
			b1=0.421492932021683E-1;
			b2=0.582840985360327E-3;
			b3=0.237676790577455E-5;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) ) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) );
			return pow(temp, m+0.5)
		}
	}
	if(m==11){
		if(w>17.5395)
		{
			return (bifactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) 
		}else{
			a1=-0.332669773770348E-1;
			a2=0.668720489602687E-3;
			b1=0.363057685289467E-1;
			b2=0.345646100984643E-3;
			b3=-0.190872330373450E-5;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) ) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) );
			return pow(temp, m+0.5)
		}
	}
	if(m==12){
		if(w>18.5783)
		{
			return (bifactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) 
		}else{
			a1=-0.326241966410798E-1;
			a2=0.598705175467956E-3;
			b1=0.318680048277695E-1;
			b2=0.202419662347765E-3;
			b3=-0.362095173837973E-5;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) ) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) );
			return pow(temp, m+0.5)
		}
	}
	if(m==13){
		if(w>19.6511)
		{
			return (bifactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) 
		}else{
			a1=-0.317754368014894E-1;
			a2=0.537678595933584E-3;
			b1=0.284036027081815E-1;
			b2=0.113673420662576E-3;
			b3=-0.416076810552774E-5;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) ) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) );
			return pow(temp, m+0.5)
		}
	}
	if(m==14){
		if(w>20.7839)
		{
			return (bifactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) 
		}else{
			a1=-0.308755854748829E-1;
			a2=0.485046451960769E-3;
			b1=0.255694625434059E-1;
			b2=0.542010192055080E-4;
			b3=-0.424759498527876E-5;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) ) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) );
			return pow(temp, m+0.5)
		}
	}
	if(m==15){
		if(w>21.9998)
		{
			return (bifactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) 
		}else{
			a1=-0.3001438066719997E-1;
			a2=0.439983033427912E-3;
			b1=0.231478878674366E-1;
			b2=0.105546581596674E-4;
			b3=-0.418932957034726E-5;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) ) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) );
			return pow(temp, m+0.5)
		}
	}
	if(m==16){
		if(w>20.9225)
		{
			return (bifactorial(2*m-1)/pow(2*w, m+0.5))*sqrt(PI/2) 
		}else{
			a1=-0.288346417991609E-1;
			a2=0.397161796318408E-3;
			b1=0.215021933120724E-1;
			b2=-0.128592457453950E-5;
			b3=-0.362120651688135E-5;
			double temp = (a0 + a1*pow(w,1) + a2*pow(w,2) ) /
				(1 + b1*pow(w,1) + b2*pow(w,2) + b3*pow(w,3) );
			return pow(temp, m+0.5)
		}
	}
	cerr << "Fpade is not yet completed!" << endl; exit(-1);
	return 0;
}

*/







#endif




