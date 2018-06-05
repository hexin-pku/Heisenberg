#ifndef Hsbg_SCF_H
#define Hsbg_SCF_H

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>

#include "Eigen/Dense"
#include <unsupported/Eigen/CXX11/Tensor>

#include "Hsbg_Const.h"
#include "Hsbg_Tools.h"
#include "Hsbg_Point.h"
#include "Hsbg_Geom.h"
#include "Hsbg_Orbital.h"
#include "Hsbg_Basis.h"
#include "Hsbg_Parser.h"
#include "Hsbg_InteG.h"

using namespace std;
using namespace Hsbg;
using namespace Eigen;

typedef Eigen::Tensor< double , 4 > Tensor4D;

class HScf
{	
	private:
		ofstream report; 
	public:
		HTask* tasklink;
		HBasis* HBlink;
		//ifstream readbs;
		//ofstream report;
		
		string Scf_name;
		int nsz;
		int nocc ;
		
		double threshold;
		//double **S;
		//double **H;
		//double **G;
		//double **F;
		//double **C;
		//double **P;
		
		MatrixXd 	eigen_S;
		MatrixXd 	eigen_H;
		MatrixXd	eigen_G; // for test
		MatrixXd	eigen_J; // for test
		MatrixXd	eigen_K; // for test
		Tensor4D	eigen_ERI;
		MatrixXd	eigen_F;
		MatrixXd 	eigen_C;
		MatrixXd 	eigen_P;
		MatrixXd	eigen_X; // X=S^(-1/2)
		MatrixXd	eigen_Y; // Y=S^(1/2), for test
		MatrixXd	eigen_Fp;
		MatrixXd	eigen_Cp;
		
		// for analysis
		bool		do_loop;
		VectorXd	eigen_E;
		double		E;
		double		E_old;
		MatrixXd	eigen_P_old;
		int* 		list;
	
	HScf(const HScf&);
	HScf(string scf_name){
		this->Scf_name = scf_name;
		this->do_loop = true;
	}
	
	~HScf()
	{
		this->report.close();
	}
	
	int set_Threshold(double myeps)
	{
		this->threshold = myeps;
		return 0;
	}
	
	int set_Space(HTask& HT)
	{
		this->tasklink = &HT;
		this->HBlink = &HT.TaskBasis;
		this->nsz = HT.N_set;  // wothout add 1
		cout << "size of the space: " << this->nsz << endl << endl;
				
		eigen_S = MatrixXd(this->nsz, this->nsz);
		eigen_H = MatrixXd(this->nsz, this->nsz);
		eigen_G = MatrixXd(this->nsz, this->nsz);
		eigen_J = MatrixXd(this->nsz, this->nsz);
		eigen_K = MatrixXd(this->nsz, this->nsz);
		eigen_ERI = Tensor4D(this->nsz, this->nsz, this->nsz, this->nsz);
		eigen_F = MatrixXd(this->nsz, this->nsz);
		eigen_C = MatrixXd(this->nsz, this->nsz);
		eigen_P = MatrixXd(this->nsz, this->nsz);
		eigen_X = MatrixXd(this->nsz, this->nsz);
		eigen_Fp = MatrixXd(this->nsz, this->nsz);
		eigen_Cp = MatrixXd(this->nsz, this->nsz);
		//
		eigen_P_old = MatrixXd(this->nsz, this->nsz);
		this->E_old = 0;
		this->nocc = (this->tasklink->Nelec + 1) / 2;
		this->list = new int[this->nocc];
		//
		this->report.open("jj.txt");

		return 0;
	}
	
	int calc_SHERI(HBasis &HB, MatrixXd &S, MatrixXd &H, Tensor4D &G);
	int guess_P();
	int calc_XY();
	int calc_Fock();
	int calc_Cprim();
	int do_Cprim2C();
	int calc_PE();
	int check_Loop(int cnt);
	int report_SCF();
	double get_NE();
	double m_Diff(MatrixXd &M, MatrixXd &N);
	
	int loop_SCF()
	{
		double E_old = this->E + 10;
		int cnt=0;
		
		this->calc_SHERI(this->tasklink->TaskBasis, this->eigen_S, this->eigen_H, this->eigen_ERI);
		this->calc_XY();
		this->guess_P();
		
		while( this->do_loop && cnt<128)
		{	
			E_old = this->E;
			cout << "######################### loop " << cnt+1 << " ##########################" << endl;
			this->report << "######################### loop " << cnt+1 << " ##########################" << endl;
			this->calc_Fock();
			this->calc_Cprim();
			this->do_Cprim2C();
			this->calc_PE();
			
			cnt++;
			this->check_Loop(cnt);
		}
		cout << "#########################  DONE  ##########################" << endl;
			this->report << "#########################  DONE  ##########################" << endl;
		//this->anal_Result();
		return 0;
	}
};

int HScf::calc_SHERI(HBasis &HB, MatrixXd &S, MatrixXd &H, Tensor4D &ERI)
{
	ofstream fS, fH, fERI;
	fS.open("S.dat");
	fH.open("H.dat");
	fERI.open("ERI.dat");
	
	
	for(int i1=1; i1<= HB.Natom; i1++)
	{
		for(int j1=0; j1 < HB.basis[i1].ncgto; j1++)
		{
			int m1 = HB.idmap[i1]+j1;

			for(int i2=i1; i2<= HB.Natom; i2++)
			{	
				for(int j2=(i2==i1?j1:0); j2 < HB.basis[i2].ncgto; j2++)
				{
					int m2 = HB.idmap[i2]+j2;

					// Sij matrix
					S(m1, m2) = 
					integral_S_sstype(	HB.basis[i1].cgto[j1], HB.basis[i2].cgto[j2]);
					// S exchange
					if(m1!=m2) S(m2, m1) = S(m1, m2);
					// S fout, m1 < m2
					fS << m1 << " " << m2 << " "<< S(m1, m2) << endl << endl;		//       fout
					
					// Hij matrix
					H(m1, m2) =  
					integral_T_sstype(	HB.basis[i1].cgto[j1], HB.basis[i2].cgto[j2]);
					for(int k=1;k<= HB.Natom;k++)
					{
						H(HB.idmap[i1]+j1, HB.idmap[i2]+j2) +=
						integral_V_sstype(	HB.basis[i1].cgto[j1], HB.basis[i2].cgto[j2],
											HB.basis[k]);
					}
					// H exchange
					if(m1!=m2) H(m2, m1) = H(m1, m2);
					// H fout, m1 < m2
					fH << m1 << " " << m2 << " "<< H(m1, m2) << endl << endl;		//       fout
					
					for(int i3=1; i3<= HB.Natom; i3++)
					{
						for(int j3=0; j3 < HB.basis[i3].ncgto; j3++)
						{
							int m3 = HB.idmap[i3]+j3;
							
							for(int i4=i3; i4<= HB.Natom; i4++)
							{
								for(int j4=(i4==i3?j3:0); j4 < HB.basis[i4].ncgto; j4++)
								{
									int m4 = HB.idmap[i4]+j4;
									
									// ERI matrix storage format
									// m1<m2, m3<m4, m2*(m2+1)+m1 < m4*(m4+1)+m3
									if(m2*(m2+1)+2*m1 <= m4*(m4+1)+2*m3)
									{
										ERI( m1, m2, m3, m4)
										= integral_ERI_sstype(	HB.basis[i1].cgto[j1], HB.basis[i2].cgto[j2],
															HB.basis[i3].cgto[j3], HB.basis[i4].cgto[j4]);
										// ERI exchange
										ERI( m2, m1, m3, m4) = ERI( m1, m2, m3, m4);
										ERI( m1, m2, m4, m3) = ERI( m1, m2, m3, m4);
										ERI( m2, m1, m4, m3) = ERI( m1, m2, m3, m4);
										ERI( m3, m4, m1, m2) = ERI( m1, m2, m3, m4);
										ERI( m3, m4, m2, m1) = ERI( m1, m2, m3, m4);
										ERI( m4, m3, m1, m2) = ERI( m1, m2, m3, m4);
										ERI( m4, m3, m2, m1) = ERI( m1, m2, m3, m4);
										// ERI fout
										fERI << m1 << " " << m2 << " " << m3 << " " << m4 << " " <<
											ERI(m1, m2, m3, m4) << endl << endl;		//       fout	
									}else{
										continue;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	fS.close();
	fH.close();
	fERI.close();
	this->report << "Inte   S_AO: " << endl << S << endl << endl;	
	this->report << "Inte   H_AO: " << endl << H << endl << endl;
	this->report << "Inte ERI_AO: " << endl << ERI << endl << endl;
	// for test
	/*  
	cout << "Inte   S_AO: " << endl << S << endl << endl;	
	cout << "Inte   H_AO: " << endl << H << endl << endl;
	cout << "Inte ERI_AO: " << endl << ERI << endl << endl;
	*/
	return 0;
}

int HScf::guess_P()
{
	// from H we guess a C
	this->eigen_F = this->eigen_H;
	this->eigen_Fp = this->eigen_X.transpose() * this->eigen_F * this->eigen_X;
	this->calc_Cprim();
	this->do_Cprim2C();
	this->calc_PE();
	return 0;
}


int HScf::calc_XY()
{
	EigenSolver<MatrixXd> eig(this->eigen_S);
	MatrixXd L1(this->nsz, this->nsz);	MatrixXd L2(this->nsz, this->nsz);
	L1.setZero(L1.rows(),L1.cols());	
	L2.setZero(L2.rows(),L2.cols());
	
	VectorXd e = eig.eigenvalues().real(); // it should to be all real
	MatrixXd T = eig.eigenvectors().real();
	for(int i=0;i<this->nsz;i++)
	{
		L1(i,i) = 1/sqrt(e(i));
		L2(i,i) = sqrt(e(i));
	}
	this->eigen_X = T * L1 * T.transpose();
	this->eigen_Y = T * L2 * T.transpose();
					/*	note: 	
								S*Tj=e*Tj,	Ti^*S*Tk = e*delta_ik
								T^*S*T = L
						or		S = T * L * T^
						so
								S^(-1/2) = 	T * L^(-1/2) * T^ === X
								S^(1/2)  =	T * L^(1/2) * T^	=== Y
					*/
	return 0;
}

int HScf::calc_Fock()
{
	for(int i1=1; i1 <= this->HBlink->Natom; i1++)
	{
		for(int j1=0; j1< this->HBlink->basis[i1].ncgto; j1++)
		{
			int m1 = this->HBlink->idmap[i1]+j1;
			
			for(int i2=1; i2 <= this->HBlink->Natom; i2++)
			{
				for(int j2=0; j2< this->HBlink->basis[i2].ncgto; j2++)
				{
					int m2 = this->HBlink->idmap[i2]+j2;
				
					this->eigen_F(m1, m2) = this->eigen_H(m1, m2);
					this->eigen_G(m1, m2) = 0;
					this->eigen_J(m1, m2) = 0;
					this->eigen_K(m1, m2) = 0;
					// revise
					for(int i3=1; i3 <= this->HBlink->Natom; i3++)
					{
						for(int j3=0; j3< this->HBlink->basis[i3].ncgto; j3++)
						{
							int m3 = this->HBlink->idmap[i3]+j3;
							
							for(int i4=1; i4 <= this->HBlink->Natom; i4++)
							{
								for(int j4=0; j4< this->HBlink->basis[i4].ncgto; j4++)
								{
									int m4 = this->HBlink->idmap[i4]+j4;
									
									this->eigen_J(m1, m2) 
													+= this->eigen_ERI(m1, m2, m3, m4)* this->eigen_P(m4, m3);
									this->eigen_K(m1, m2) 
													+= this->eigen_ERI(m1, m4, m3, m2)* this->eigen_P(m4, m3);
									this->eigen_G(m1, m2) += (2* this->eigen_ERI(m1, m2, m3, m4)
									-1* this->eigen_ERI(m1, m4, m3, m2) ) * this->eigen_P(m4, m3);
								}
							}
						}
					}
					this->eigen_F(m1, m2) += this->eigen_G(m1, m2);
				}
				
			}
		}
	}
	// additionally, also calculate the Fprime
	this->eigen_Fp = this->eigen_X.transpose() * this->eigen_F * this->eigen_X;
	
	//test
	/*
	cout << "test G_AO" << endl << this->eigen_G << endl << endl;
	cout << "test F_AO" << endl << this->eigen_F << endl << endl;
	cout << "test F_pm" << endl << this->eigen_Fp << endl << endl;
	cout << "test G_pm" << endl << this->eigen_X.transpose() * this->eigen_G * this->eigen_X << endl << endl;
	cout << "test H_pm" << endl << this->eigen_X.transpose() * this->eigen_H * this->eigen_X << endl << endl;
	*/
	return 0;
}

int HScf::calc_Cprim()
{
	EigenSolver<MatrixXd> eig(this->eigen_Fp);
	this->eigen_E.setZero(this->nsz);
	
	this->eigen_Cp = eig.eigenvectors().real();   // only c might be complex, how to fix it?
	this->eigen_E = eig.eigenvalues().real();
	

	for(int i=0; i< this->nsz; i++)
	{
		if(i < this->nocc){ this->list[i] = i;continue;}
		for(int j=0; j<this->nocc; j++)
		{
			if( this->eigen_E(this->list[j]) > this->eigen_E(i) ) 
			{
				this->list[j]=i;
				break;
			}
		}
	}
	// for test
	/*
	cout << "all E: " << endl << this->eigen_E << endl << endl;
	cout << "index of min : "; for(int i=0;i<this->nocc;i++) cout << "i-th "<< i << " " << this->list[i] << endl ;
	*/
	return 0;
}

int HScf::do_Cprim2C()
{
	this->eigen_C = this->eigen_X * this->eigen_Cp;
	// for test
	/*
	cout << "test C_MO" << endl << this->eigen_Cp << endl << endl;
	cout << "test C_AO" << endl << this->eigen_C << endl << endl;
	cout << "test S_MO : C^ * S * C" << endl << this->eigen_C.transpose() * this->eigen_S * this->eigen_C << endl << endl;
	cout << "test F_MO : C^ * F * C" << endl << this->eigen_C.transpose() * this->eigen_F * this->eigen_C << endl << endl;

	cout << "test G_MO" << endl << this->eigen_C.transpose() * this->eigen_G * this->eigen_C << endl << endl;
	cout << "test J_MO" << endl << this->eigen_C.transpose() * this->eigen_J * this->eigen_C << endl << endl;
	cout << "test K_MO" << endl << this->eigen_C.transpose() * this->eigen_K * this->eigen_C << endl << endl;
	cout << "test H_MO" << endl << this->eigen_C.transpose() * this->eigen_H * this->eigen_C << endl << endl;
	*/
	return 0;
}

int HScf::calc_PE()
{	
	for(int i1=1; i1 <= this->HBlink->Natom; i1++)
	{
		for(int j1=0; j1< this->HBlink->basis[i1].ncgto; j1++)
		{
			int m1 = this->HBlink->idmap[i1]+j1;
			
			for(int i2=1; i2 <= this->HBlink->Natom; i2++)
			{
				for(int j2=0; j2< this->HBlink->basis[i2].ncgto; j2++)
				{
					int m2 = this->HBlink->idmap[i2]+j2;
					
					this->eigen_P(m1, m2)=0;			
					for(int k=0; k<this->nocc; k++)
					{
						this->eigen_P(m1, m2 )
						+= 1.000 *this->eigen_C(m1, this->list[k] )
						* this->eigen_C(m2, this->list[k] ); //?? it is right?, may be there is wrong!
					}
				}
			}
		}
	}
	this->E = get_NE();
	for(int k=0; k<this->nocc; k++)
	{
		this->E += this->eigen_E(list[k]);// + (this->eigen_P*this->eigen_H)(list[k],list[k]);
	}
	this->E += (this->eigen_P*this->eigen_H).trace(); // trace only on ocuppied ??
	this->report_SCF();
	return 0;
}

int HScf::report_SCF()
{	this->report << "test C_MO" << endl << this->eigen_Cp << endl << endl;
	this->report << "test C_AO" << endl << this->eigen_C << endl << endl;
	this->report << "test S_MO : C^ * S * C" 
		<< endl << this->eigen_C.transpose() * this->eigen_S * this->eigen_C << endl << endl;
	this->report << "test F_MO : C^ * F * C" 
		<< endl << this->eigen_C.transpose() * this->eigen_F * this->eigen_C << endl << endl;
	this->report << "test G_MO" << endl << this->eigen_C.transpose() * this->eigen_G * this->eigen_C << endl << endl;
	this->report << "test J_MO" << endl << this->eigen_C.transpose() * this->eigen_J * this->eigen_C << endl << endl;
	this->report << "test K_MO" << endl << this->eigen_C.transpose() * this->eigen_K * this->eigen_C << endl << endl;
	this->report << "test H_MO" << endl << this->eigen_C.transpose() * this->eigen_H * this->eigen_C << endl << endl;
	this->report << "Density P" << endl << this->eigen_P << endl << endl;
	this->report << "Energy " << this->E << "    " << "Convergence: " << E_old-this->E << endl << endl;
	return 0;
}

double HScf::get_NE()
{
	double sum=0;
	for(int i=1;i<=this->HBlink->Natom; i++)
	{
		for(int j=i+1;j<=this->HBlink->Natom; j++)
		{
			sum += this->HBlink->basis[i].znum * this->HBlink->basis[j].znum 	/	sqrt(
			(this->HBlink->basis[i].get_Point() - this->HBlink->basis[j].get_Point() )*
			(this->HBlink->basis[i].get_Point() - this->HBlink->basis[j].get_Point() )	);
		}
	}
	return sum;
}

int HScf::check_Loop(int cnt)
{
	if( abs(this->E - this->E_old) < this->threshold && this->m_Diff(this->eigen_P,this->eigen_P_old) )
	{
		this->do_loop = false;
	}
	cout << "# SCF LOOP RESULT" << cnt << endl << endl;
	cout << "    E: " << this->E << "    " << "Convergence: " << E_old-this->E << endl << endl;
	this->E_old = this->E;
	return 0;
}

double HScf::m_Diff(MatrixXd &M, MatrixXd &N)
{
	double find=0;
	if(M.rows()!=N.rows() || M.cols()!=N.cols()) {cerr<< "Matrix size not match"<< endl << endl; exit(-1);}
	for(int i=0; i< M.rows(); i++)
	{
		for(int j=0; j<M.cols(); j++)
		{
			if(find < abs(M(i,j)-N(i,j))) find = abs(M(i,j)-N(i,j));
		}
	}
	return find;
}

#endif



