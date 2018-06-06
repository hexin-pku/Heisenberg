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
#include "Hsbg_Point.h"
#include "Hsbg_Geom.h"
#include "Hsbg_Orbital.h"
#include "Hsbg_Basis.h"
#include "Hsbg_Parser.h"


using namespace std;
using namespace Hsbg;
using namespace Eigen;


//  Integral of contract gauss fucntion

	double integral_S(	 HOrbital_cgto& cgto1,  HOrbital_cgto& cgto2);
	double integral_ERI(  HOrbital_cgto& cgto1, 
								 HOrbital_cgto& cgto2, 
								 HOrbital_cgto& cgto3, 
								 HOrbital_cgto& cgto4);
	double integral_T( 	 HOrbital_cgto& cgto1,  HOrbital_cgto& cgto2);
	double integral_V( 	 HOrbital_cgto& cgto1,  HOrbital_cgto& cgto2, 
								BPoint& P );
	
	// s-s integral type
	double integral_S_sstype(	 HOrbital_cgto& cgto1,  HOrbital_cgto& cgto2);
	double integral_ERI_sstype(  HOrbital_cgto& cgto1, 
								 HOrbital_cgto& cgto2, 
								 HOrbital_cgto& cgto3, 
								 HOrbital_cgto& cgto4);
	double integral_T_sstype( 	 HOrbital_cgto& cgto1,  HOrbital_cgto& cgto2);
	double integral_V_sstype( 	 HOrbital_cgto& cgto1,  HOrbital_cgto& cgto2, 
								BPoint& P );
	// s-p integral type
	double integral_S_sptype(	 HOrbital_cgto& cgto1,  HOrbital_cgto& cgto2){return 0;};
	double integral_ERI_sptype(  HOrbital_cgto& cgto1, 
								 HOrbital_cgto& cgto2, 
								 HOrbital_cgto& cgto3, 
								 HOrbital_cgto& cgto4){return 0;};
	double integral_T_sptype( 	 HOrbital_cgto& cgto1,  HOrbital_cgto& cgto2){return 0;};
	double integral_V_sptype( 	 HOrbital_cgto& cgto1,  HOrbital_cgto& cgto2, 
								BPoint& P ){return 0;};
	// p-p integral type
	double integral_S_pptype(	 HOrbital_cgto& cgto1,  HOrbital_cgto& cgto2){return 0;};
	double integral_ERI_pptype(  HOrbital_cgto& cgto1, 
								 HOrbital_cgto& cgto2, 
								 HOrbital_cgto& cgto3, 
								 HOrbital_cgto& cgto4){return 0;};
	double integral_T_pptype( 	 HOrbital_cgto& cgto1,  HOrbital_cgto& cgto2){return 0;};
	double integral_V_pptype( 	 HOrbital_cgto& cgto1,  HOrbital_cgto& cgto2, 
								BPoint& P ){return 0;};
	
	double integral_S(	 HOrbital_cgto& cgto1,  HOrbital_cgto& cgto2){
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
	
	double integral_ERI(  HOrbital_cgto& cgto1, 
								 HOrbital_cgto& cgto2, 
								 HOrbital_cgto& cgto3, 
								 HOrbital_cgto& cgto4);
	
	double integral_T( 	 HOrbital_cgto& cgto1,  HOrbital_cgto& cgto2)
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
	
	double integral_V( 	 HOrbital_cgto& cgto1,  HOrbital_cgto& cgto2, 
								BPoint& P )
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

double integral_S_sstype( HOrbital_cgto& cgto1,  HOrbital_cgto& cgto2)
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
	 HOrbital_cgto& cgto1, 
	 HOrbital_cgto& cgto2, 
	 HOrbital_cgto& cgto3, 
	 HOrbital_cgto& cgto4
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
	 HOrbital_cgto& cgto1, 
	 HOrbital_cgto& cgto2 
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
	 HOrbital_cgto& cgto1, 
	 HOrbital_cgto& cgto2, 
	 BPoint& P
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

#endif




