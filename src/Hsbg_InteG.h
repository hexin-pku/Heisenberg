#ifndef Hsbg_InteG_H
#define Hsbg_InteG_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include "Eigen/Dense"

#include "Hsbg_Const.h"
#include "Hsbg_Tools.h"
#include "Hsbg_Global.h"
#include "Hsbg_Geom.h"
#include "Hsbg_Orbital.h"
#include "Hsbg_Basis.h"
#include "Hsbg_Parser.h"


using namespace std;
using namespace Hsbg;
using namespace Eigen;


//class HInteG
//{
//	public:
//	
	double integral_S_sstype(	const HOrbital_cgto& cgto1, const HOrbital_cgto& cgto2);
	double integral_ERI_sstype( const HOrbital_cgto& cgto1, 
								const HOrbital_cgto& cgto2, 
								const HOrbital_cgto& cgto3, 
								const HOrbital_cgto& cgto4);
	double integral_T_sstype( 	const HOrbital_cgto& cgto1, const HOrbital_cgto& cgto2);
	double integral_V_sstype( 	const HOrbital_cgto& cgto1, const HOrbital_cgto& cgto2, 
								double cx, double cy, double cz );
	double normGTO(int L, int M, int N, double alpha);
	void calc_SHG(HBasis &HB, vector2d &S, vector2d &H, vector4d &G);
//	
//}

void calc_SHG(HBasis &HB, vector2d &S, vector2d &H, vector4d &G)
{
	for(int i1=1, i1<= HB.Natom; i1++)
	{
		for(int j1=0; j1 < HB.basis[i1].ncgto; j1++)
		{
			for(int i2=1; i2<= HB.Natom; i2++)
			{
				for(int j2=0; j2 < HB.basis[i2].ncgto; j2++)
				{
					S[HB.idmap[i1]+j1, HB.idmap[i2]+j2] = 
					integral_S_sstype(	HB.basis[i1].cgto[j1], HB.basis[i2].cgto[j2]);
					H[HB.idmap[i1]+j1, HB.idmap[i2]+j2] = 
					integral_T_sstype(	HB.basis[i1].cgto[j1], HB.basis[i2].cgto[j2]);
					for(int k=1;k<= HB.Natom;k++)
					{
						H[HB.idmap[i1]+j1, HB.idmap[i2]+j2]+=
						integral_V_sstype(	HB.basis[i1].cgto[j1], HB.basis[i2].cgto[j2],
											HB.basis[k].x, HB.basis[k].y, HB.basis[k].z);
					}
					for(int i3=1, i3<= HB.Natom; i3++)
					{
						for(int j3=0; j3 < HB.basis[i3].ncgto; j3++)
						{
							for(int i3=1, i3<= HB.Natom; i3++)
							{
								for(int j3=0; j3 < HB.basis[i3].ncgto; j3++)
								{
									G[HB.idmap[i1]+j1, HB.idmap[i2]+j2, HB.idmap[i3]+j3, HB.idmap[i4]+j4]
									= integral_ERI_sstype(	HB.basis[i1].cgto[j1], HB.basis[i2].cgto[j2],
															HB.basis[i3].cgto[j3], HB.basis[i4].cgto[j4]);
								}
							}
						}
					}
				}
			}
		}
	}
}


double integral_S_sstype(const HOrbital_cgto& cgto1, const HOrbital_cgto& cgto2)
{
	double AB2, sum;
	int n1, n2;
	
	AB2 = cgto1.dist_A2( cgto.get_Point() );
	
	n1 = cgto1.ncgto;
	n2 = cgto2.ncgto;
	
	sum = 0;
	
	double alpha1, alpha2, coeff1, coeff2, norm;
	
	for(int i = 0; i < n1; i++){
		cgto1.get_CA(i, coeff1, alpha1);
		for(int j = 0; j < n2; j++){
			cgto2.getOrbital(alpha2, coeff2, j);
			norm = std::pow(4*alpha1*alpha2/(PI*PI), 0.75);
			sum += norm * coeff1 * coeff2 * 
			std::exp(- AB2*alpha1*alpha2 / (alpha1+alpha2)) *
			std::pow(PI/(alpha1+alpha2),1.5);
		}
	}	
	return sum;
}

double integral_ERI_sstype(
	const HOrbital_cgto& cgto1, 
	const HOrbital_cgto& cgto2, 
	const HOrbital_cgto& cgto3, 
	const HOrbital_cgto& cgto4
){
	double eps = 1e-10;
	double sum, F0;
	double AB2, CD2;
	int n1, n2, n3, n4;
	cgto1.getPosition(ax, ay, az);
	cgto2.getPosition(bx, by, bz);
	cgto3.getPosition(cx, cy, cz);
	cgto4.getPosition(dx, dy, dz);
	AB2 = cgto1.dist_A2( cgto2.get_Point() );
	CD2 = cgto3.dist_A2( cgto4.get_Point() );
	
	n1 = cgto1.ncgto;
	n2 = cgto2.ncgto;
	n3 = cgto3.ncgto;
	n4 = cgto4.ncgto;
	
	sum = 0;
	double a1, a2, a3, a4, w, a12, a34, a1234, c1, c2, c3, c4, PQx, PQy, PQz, PQ2, norm;
	for(int i = 0; i < n1; ++i){
		cgto1.get_CA(i, c1, a1);
		for(int j = 0; j < n2; ++j){
			cgto2.get_CA(j, c2, a2);
			a12 = a1 + a2;
			for(int k = 0; k < n3; ++k){
				cgto3.get_CA(k, c3, a3);
				for(int m = 0; m < n4; ++m){
					cgto4.get_CA(m, c4, a4);
					a34 = a3 + a4;
					PQx = (a1*cgto1.x + a2*cgto2.x)/a12 - (a3*cgto3.x + a4*cgto4.x)/a34;
					PQy = (a1*cgto1.y + a2*cgto2.y)/a12 - (a3*cgto3.y + a4*cgto4.y)/a34;
					PQz = (a1*cgto1.z + a2*cgto2.z)/a12 - (a3*cgto3.z + a4*cgto4.z)/a34;
					PQ2 = PQx*PQx + PQy*PQy + PQz*PQz;
					
					w = PQ2*a12*a34 / (a12 + a34);
					F0 = (std::acgto(w) < eps ? 1.0 : 0.5 * std::sqrt(PI / w) * std::erf(std::sqrt(w)));
					
					norm = std::pow(16*a*b*c*d/(PI*PI*PI*PI), 0.75);
					sum += 2*norm*c1*c2*c3*c4*std::pow(PI, 2.5) * 
					std::exp(- a1*a2*AB2/a12 - a3*a4*CD2/a34) * F0/
					(a12 * a34 * std::sqrt(a12 + a34));
				}
			}
		}
	}
	return sum;
}

double integral_T_sstype(
	const HOrbital_cgto& cgto1, 
	const HOrbital_cgto& cgto2
){
	double AB2, sum;
	int n1, n2;
	
	AB2 = cgto1.dist_A2( cgto2.get_Point() );
	
	n1 = cgto1.ncgto;
	n2 = cgto2.ncgto;
	
	sum = 0;
	
	double alpha1, alpha2, coeff1, coeff2, temp, norm;
	
	for(int i = 0; i < n1; ++i){
		cgto1.get_CA(i, coeff1, alpha1);
		for(int j = 0; j < n2; ++j){
			cgto2.get_CA(coeff2, alpha2, j);
			temp = alpha1 * alpha2 / (alpha1 + alpha2);
			norm = std::pow(4*alpha1*alpha2/(PI*PI), 0.75);
			sum += norm * coeff1 * coeff2 * temp * (3 - 2*temp*AB2) * 
			std::exp(- temp*AB2) * std::pow(PI / (alpha1+alpha2),1.5);
		}
	}
	return sum;
}

double integral_V_sstype(
	const HOrbital_cgto& cgto1, 
	const HOrbital_cgto& cgto2, 
	double cx,
	double cy,
	double cz
){
	double eps = 1e-10;
	double px, py, pz, AB2, sum, F0;
	int n1, n2;
	
	AB2 = cgto1.dist_A2( cgto2.get_Point() );
	n1 = cgto1.ncgto;
	n2 = cgto2.ncgto;
	
	sum = 0;
	double a1, a2, c1, c2, a12, w, norm;
	for(int i = 0; i < n1; ++i){
		cgto1.get_CA(i, c1, a1);
		for(int j = 0; j < n2; ++j){
			cgto2.get_CA(j, c2, a2);
			a12 = a1 + a2;
			px = (cgto1.x*a1 + cgto2.x*a2) / a12;
			py = (cgto1.y*a1 + cgto2.y*a2) / a12;
			pz = (cgto1.z*a1 + cgto2.z*a2) / a12;			
			w = a12 * 
			((px - cx)*(px - cx) + (py - cy)*(py - cy) + (pz - cz)*(pz - cz));
			F0 = (std::acgto(w) < eps ? 1.0 : 0.5 * std::sqrt(PI / w) * std::erf(std::sqrt(w)));
			norm = std::pow(4*a1*a2/(PI*PI), 0.75);
			sum += norm * c1 * c2 * 2 * PI *
			std::exp(- AB2*a1*a2 / a12) * F0 / a12;
		}
	}
	return sum;
}

double normGTO(int L, int M, int N, double alpha){
	int n1 = 2*L - 1, 
	n2 = 2*M - 1, 
	n3 = 2*N - 1,
	total = 1;
	for(int i = 3; i <= n1; i+=2) total *= i;
	for(int i = 3; i <= n2; i+=2) total *= i;
	for(int i = 3; i <= n3; i+=2) total *= i;
	return std::pow(2*alpha/PI, 0.75) * std::sqrt(std::pow(4*alpha, L+M+N)/static_cast<double>(total));	
}


#endif




