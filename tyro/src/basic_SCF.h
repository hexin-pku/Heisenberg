#ifndef basic_SCF_H
#define basic_SCF_H

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>

#include <vector>

#include "Eigen/Dense"
#include <unsupported/Eigen/CXX11/Tensor>

#include "basic_Const.h"
#include "basic_Tools.h"
#include "basic_Global.h"
#include "basic_Parser.h"
#include "basic_Integral.h"

using namespace std;
using namespace Eigen;

typedef Eigen::Tensor< double , 4 > Tensor4D;

class SCFer
{	
	private:
		ofstream report; 
	public:
		Task* tasklink;
		System* SYSlink;
		
		string Scf_name;
		int nsz;
		int nocc ;
		
		double threshold;
		
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
	
	SCFer(const SCFer&);
	SCFer(string scf_name){
		this->Scf_name = scf_name;
		this->do_loop = true;
	}
	
	~SCFer()
	{
		this->report.close();
	}
	
	int set_Threshold(double myeps)
	{
		this->threshold = myeps;
		return 0;
	}
	
	int set_Space(Task& T)
	{
		this->tasklink = &T;
		this->SYSlink = &T.Sys;
		this->nsz = T.Nbasis;  // wothout add 1
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
		this->report.open(T.Logfile.data());

		return 0;
	}
	
	int calc_SHERI(System &SYS, MatrixXd &S, MatrixXd &H, Tensor4D &G);
	int guess_P();
	int calc_XY();
	int calc_Fock();
	int calc_Cprim();
	int tr_Cprim2C();
	int calc_PE();
	int check_Loop(int cnt);
	int report_SCF();
	double get_NE();
	double m_Diff(MatrixXd &M, MatrixXd &N);
	
	int loop_SCF()
	{
		double E_old = this->E + 10;
		int cnt=0;
		
		this->report << *(this->tasklink) << endl;
		this->calc_SHERI(this->tasklink->Sys, this->eigen_S, this->eigen_H, this->eigen_ERI);
		this->calc_XY();
		// cout << this->eigen_X << endl;
		
		this->guess_P();
		
		while( this->do_loop && cnt<32)
		{	
			E_old = this->E;
			cout << "######################### loop " << cnt+1 << " ##########################" << endl;
			this->report << "######################### loop " << cnt+1 << " ##########################" << endl;
			this->calc_Fock();
			this->calc_Cprim();
			this->tr_Cprim2C();
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

int SCFer::calc_SHERI(System &SYS, MatrixXd &S, MatrixXd &H, Tensor4D &ERI)
{
	ofstream fS, fH, fERI;
	fS.open("S.dat");
	fH.open("H.dat");
	fERI.open("ERI.dat");
	
	
	for(int m1=0; m1<SYS.Nbasis; m1++)
	{
		for(int m2=m1; m2<SYS.Nbasis; m2++)
		{
			//S(m1, m2) = integral_S_sstype(SYS[m1], SYS[m2]);
			S(m1, m2) = IntecGTO_S(SYS[m1], SYS[m2]);
			if(m1!=m2) S(m2, m1) = S(m1, m2);
			fS << m1 << " " << m2 << " "<< S(m1, m2) << endl << endl;
			H(m1, m2) = integral_T_sstype( SYS[m1], SYS[m2] );
			//H(m1, m2) = IntecGTO_T(SYS[m1], SYS[m2]);
			
			for(int k=1; k<= SYS.Natom; k++)
			{
				H(m1, m2) += integral_V_sstype(SYS[m1], SYS[m2], SYS.atoms[k]);
				//H(m1, m2) += IntecGTO_V(SYS[m1], SYS[m2], SYS.atoms[k]);
			}
			if(m1!=m2) H(m2, m1) = H(m1, m2);
			fH << m1 << " " << m2 << " "<< H(m1, m2) << endl << endl;
			
			for(int m3=0; m3<SYS.Nbasis; m3++)
			{
				for(int m4=m3; m4<SYS.Nbasis; m4++)
				{
					if(m2*(m2+1)+2*m1 <= m4*(m4+1)+2*m3)
					{
						ERI( m1, m2, m3, m4) = integral_ERI_sstype( SYS[m1], SYS[m2], SYS[m3], SYS[m4] );
						//ERI( m1, m2, m3, m4) = IntecGTO_ERI( SYS[m1], SYS[m2], SYS[m3], SYS[m4] );
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
	
	fS.close();
	fH.close();
	fERI.close();
	this->report << "Inte   S_AO: " << endl << S << endl << endl;	
	this->report << "Inte   H_AO: " << endl << H << endl << endl;
	this->report << "Inte ERI_AO: " << endl << ERI << endl << endl;
	// for test
	 
	cout << "Inte   S_AO: " << endl << S << endl << endl;	
	cout << "Inte   H_AO: " << endl << H << endl << endl;
	cout << "Inte ERI_AO: " << endl << ERI << endl << endl;
	
	return 0;
}

int SCFer::guess_P()
{
	// from H we guess a C
	this->eigen_F = this->eigen_H;
	this->eigen_Fp = this->eigen_X.transpose() * this->eigen_F * this->eigen_X;
	this->calc_Cprim();
	this->tr_Cprim2C();
	this->calc_PE();
	return 0;
}


int SCFer::calc_XY()
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

int SCFer::calc_Fock()
{
	for(int m1=0; m1 < this->SYSlink->Nbasis; m1++)
	{
		for(int m2=0; m2 < this->SYSlink->Nbasis; m2++)
		{
			this->eigen_F(m1, m2) = this->eigen_H(m1, m2);
			this->eigen_G(m1, m2) = 0;
			this->eigen_J(m1, m2) = 0;
			this->eigen_K(m1, m2) = 0;
			for(int m3=0; m3 < this->SYSlink->Nbasis; m3++)
			{
				for(int m4=0; m4 < this->SYSlink->Nbasis; m4++)
				{
					this->eigen_J(m1, m2) += this->eigen_ERI(m1, m2, m3, m4)* this->eigen_P(m4, m3);
					this->eigen_K(m1, m2) += this->eigen_ERI(m1, m4, m3, m2)* this->eigen_P(m4, m3);
					this->eigen_G(m1, m2) += 
							(2*this->eigen_ERI(m1, m2, m3, m4) - this->eigen_ERI(m1, m4, m3, m2) ) 
							* this->eigen_P(m4, m3);
				}
			}
			this->eigen_F(m1, m2) += this->eigen_G(m1, m2);
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

int SCFer::calc_Cprim()
{
	EigenSolver<MatrixXd> eig(this->eigen_Fp);
	this->eigen_E.setZero(this->nsz);
	
	this->eigen_Cp = eig.eigenvectors().real();
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

int SCFer::tr_Cprim2C()
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

int SCFer::calc_PE()
{	
	for(int i1=1; i1 <= this->SYSlink->Natom; i1++)
	{
		for(int j1=0; j1< this->SYSlink->atoms[i1].ncgto; j1++)
		{
			int m1 = this->SYSlink->idmap[i1]+j1;
			
			for(int i2=1; i2 <= this->SYSlink->Natom; i2++)
			{
				for(int j2=0; j2< this->SYSlink->atoms[i2].ncgto; j2++)
				{
					int m2 = this->SYSlink->idmap[i2]+j2;
					
					this->eigen_P(m1, m2)=0;			
					for(int k=0; k<this->nocc; k++)
					{
						this->eigen_P(m1, m2 )
						+= 1.000 *this->eigen_C(m1, this->list[k] )
						* this->eigen_C(m2, this->list[k] ); 
					}
				}
			}
		}
	}
	this->E = get_NE();
	for(int k=0; k<this->nocc; k++)
	{
		this->E += this->eigen_E(this->list[k]) ;//+ (this->eigen_P*this->eigen_H)(list[k],list[k]);
	}
	this->E += (this->eigen_P*this->eigen_H).trace(); // trace only on ocuppied ??
	this->report_SCF();
	return 0;
}

int SCFer::report_SCF()
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

double SCFer::get_NE()
{
	double sum=0;
	for(int i=1;i<=this->SYSlink->Natom; i++)
	{
		for(int j=i+1;j<=this->SYSlink->Natom; j++)
		{
			sum += this->SYSlink->atoms[i].znum * this->SYSlink->atoms[j].znum 	/	sqrt(
			(this->SYSlink->atoms[i].get_Point() - this->SYSlink->atoms[j].get_Point() )*
			(this->SYSlink->atoms[i].get_Point() - this->SYSlink->atoms[j].get_Point() )	);
		}
	}
	return sum;
}

int SCFer::check_Loop(int cnt)
{
	// output the result at each loop
	cout << "# SCF LOOP RESULT" << cnt << endl << endl;
	cout << "    Energy_MO values: " <<  endl << this->eigen_E << endl;
	cout << "    Occupied indexes: "; for(int i=0; i<this->nocc; i++) cout << " " << this->list[i]; cout << endl;
	cout << "    Energy_Total : " << this->E << "    " << "Convergence: " << E_old-this->E << endl << endl;
	//
	if( abs(this->E - this->E_old) < this->threshold && this->m_Diff(this->eigen_P,this->eigen_P_old) )
	{
		this->do_loop = false;
		
		// Load poems to output from file
		std::ifstream infile("attachment.dat");
		vector<string> vPoem;
		
		// Reading and spliting
		string temp;
		string verse;
		while(std::getline(infile, temp)){
			if (temp[0] == '@'){
				vPoem.push_back(verse);

				verse = "";
			} // If one verse ends
			verse += (temp + "\n");
		}
		vPoem.push_back(verse);
		
		// Choose random one and output as bottom lines
		srand((unsigned)time(NULL));
		int nRand =  rand() % vPoem.size();
		cout << vPoem[nRand] << std::endl;
		
	}
	this->E_old = this->E;
	return 0;
}

double SCFer::m_Diff(MatrixXd &M, MatrixXd &N)
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



