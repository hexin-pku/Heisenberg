#include <iostream>
#include <string>
#include <fstream>

#include "Hsbg_Const.h"
#include "Hsbg_Global.h"
#include "Hsbg_Geom.h"
#include "Hsbg_Basis.h"
#include "Hsbg_Parser.h"

#include "Eigen/Dense"

using namespace std;
using namespace Hsbg;
using namespace Eigen;

class HSolver
{
	public:
		HTask* tasklink;
		string solver_name;
		int nspace;
		Matrix S;
		Matrix J;
		Matrix K;
		Matrix F;
		Matrix H;
		Matrix G;

	int creat_Space(HTask* HT)
	{
		return 0;
	}

}

