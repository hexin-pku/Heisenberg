#ifndef _Hsbg_Basis_H_
#define _Hsbg_Basis_H_

#include <iostream>
#include <sstream>
#include <string>

#include "Hsbg_Const.h"
#include "Hsbg_Global.h"
#include "Hsbg_Orbital.h"

using namespace std;
using namespace Hsbg;

class HBasis
{
    public:
    	string  name; // gto, sto, etc, not same with the basis
    	int  Natom;
        BPoint* basis;
        
    //HBasis();
    //~HBasis();
    
    int set_Basis(int num) // num is the atom position number
    {
    	this->basis = new BPoint[num+1];
    	return 0;
    }    
    
    friend ostream &operator<<( ostream &output, HBasis &HB )
    {
        output << endl << "the Geometry of "<< HB.name << " shows as follow:" << endl;
        for(int i=1; i <= HB.Natom; i++)
        {
            output << setiosflags(ios::fixed) << setprecision(8) << HB.basis[i].name ;
            output << "         ";
            output << HB.basis[i].x;
            output << "    ";
            output << HB.basis[i].y;
            output << "    ";
            output << HB.basis[i].z << endl;
            //
        }
        output << endl;
        for(int i=1; i <= HB.Natom; i++)
        {
        	output << HB.basis[i];
        }
        return output;          
    }
};

#endif
