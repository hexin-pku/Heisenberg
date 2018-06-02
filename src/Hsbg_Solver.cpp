#include <iostream>
#include <string>
#include <fstream>

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

class HSolver
{
	public:
		HTask* tasklink;
		string solver_name;
		int nspace;
		double threshold = 0.000001;
		Matrix S;
		Matrix H;
		Matrix J;
		Matrix G;
		Matrix K;
		Matrix F;
		Matric C;
		Matrix P;
		
		// for analysis
		double E=0;
		
		
	int set_Solver(HTask& HT)
	{
		this->tasklink = HT;
		this->nspace = HT.TaskBasis.idmap[HT.TaskBasis.Natom] 
			+ HT.TaskBasis.basis[HT.TaskBasis.Natom].ncgto ;  // wothout add 1
		
	}
	
	int guess_P();
	int calc_SHG();
	int calc_Fock();
	int do_Fock2Fockprim();
	int calc_Cprim();
	int do_Cprim2C();
	int calc_PE();
	
	int loop_Solver()
	{
		double E_old = this->E + 10;
		
		this->calc_SHG();
		this->guess_P();
		
		while( std::abs( E - E_old ) > this->threshold )
		{
			this->calc_Fock();
			this->do_Fock2Fockprim();
			this->calc_C();
			this->do_Cprim2C();
			this->calc_PE();
		}
		this->anal_Result()
		return 0;
	}
	
}

