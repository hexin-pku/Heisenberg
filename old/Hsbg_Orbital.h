#ifndef Hsbg_Orbital_H
#define Hsbg_Orbital_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>

#include "Hsbg_Const.h"
#include "Hsbg_Point.h"

using namespace std;
using namespace Hsbg;

class HOrbital
{
	public:
		int L;		// lower nlm, is the type of gauss/slater
		int M;
		int N;
		double coeff;	// coefficient
		double alpha;	// alpha
		double x;
		double y;
		double z;
	
	int init_Orbit(double cc, double ee, double xx, double yy, double zz)
	{
		this->set_CA(cc,ee);
		this->set_XYZ(xx,yy,zz);
		return 0;
	}
	
	int set_CA(double cc, double aa)
	{
		this->coeff = cc;
		this->alpha = aa;
		return 0;
	}
		
	int set_LMN(int ll,int mm, int nn)
	{
		this->L=ll;
		this->M=mm;
		this->N=nn;
		return 0;
	}
	
	int set_XYZ(double xx, double yy, double zz)
	{
		this->x=xx;
		this->y=yy;
		this->z=zz;
		return 0;
	}
	
	int set_XYZ(Point P)
	{
		this->x=P.x;
		this->y=P.y;
		this->z=P.z;
		return 0;
	}
	
	Point& get_Point()
    {
    		Point* P = new Point();
    		P->aname = "NULL";
    		P->x = this->x;
    		P->y = this->y;
    		P->z = this->z;
    		return *P;
    }
	
	int get_CA(double& cc, double& aa)
	{
		cc = this->coeff;
		aa = this->alpha;
		return 0;
	}
	
	int get_LMN(int& ll,int& mm, int& nn)
	{
		ll=this->L;
		mm=this->M;
		nn=this->N;
		return 0;
	}
	
	int get_XYZ(double& xx, double& yy, double& zz)
	{
		xx=this->x;
		yy=this->y;
		zz=this->z;
		return 0;
	}
	
	double normGTO( double a)
	{
		int n1 = 2*this->L - 1, // x-direction
		n2 = 2*this->M - 1, 	// y-direction
		n3 = 2*this->N - 1,		// z-direction
		total = 1;
		
		for(int i = 3; i <= n1; i+=2) total *= i;
		for(int i = 3; i <= n2; i+=2) total *= i;
		for(int i = 3; i <= n3; i+=2) total *= i;
		return pow(2*a/PI, 0.75) * sqrt(pow(4*a, L+M+N)/static_cast<double>(total));	
	}
};


class HOrbital_cgto : public HOrbital
{
	public:
		// note BohrL = L + M + N;
		int cn;	// contraction number
		
		double* coeffs;
		double* alphas;
		HOrbital* orbits;
	
	//HOrbital_cgto(int num);
	~HOrbital_cgto()
	{
		delete[] orbits;
	}

	int set_Cgto(int num)
	{
		this->cn = num;
		this->coeffs = new double[num];
		this->alphas = new double[num];
		this->orbits = new HOrbital[num];
		return 0;
	}
	
	int get_CA(int idx, double& cc, double& aa)
	{
		if(idx>=this->cn) {cerr << "error, segment in cgto" << endl; exit(-1);}
		cc = this->coeffs[idx];
		aa = this->alphas[idx];
		return 0;
	}
	
	int set_Pars(int index, double mycoeff, double myalpha)
	{
		if( index>=0 && index < this->cn){
			this->coeffs[index] = mycoeff;
			this->alphas[index] = myalpha;
			return 0;
		}
		return -1;
	}
	
	int conv_AUnit()
	{
		this->x = this->x / c_bohr2ai;
		this->y = this->y / c_bohr2ai;
		this->z = this->z / c_bohr2ai;
	}
	
	int get_BohrL()
	{
		return this->L + this->M + this->N;
	}
	
	//double normGTO( double a); inherit from father class HOrbital 
	
	double normGTO(int k) 		// return k-th component gauss's norm of this contract-GTO, k from 0 count
	{
		if(k >= this->cn || k<0) {cerr << "error, out bounadry" << endl; exit(-1);}
		int n1 = 2*this->L - 1, // x-direction
		n2 = 2*this->M - 1, 	// y-direction
		n3 = 2*this->N - 1,		// z-direction
		total = 1;
		
		for(int i = 3; i <= n1; i+=2) total *= i;
		for(int i = 3; i <= n2; i+=2) total *= i;
		for(int i = 3; i <= n3; i+=2) total *= i;
		return pow(2*this->alphas[k]/PI, 0.75) * sqrt(pow(4*this->alphas[k], L+M+N)/static_cast<double>(total));	
	}
};

#endif





