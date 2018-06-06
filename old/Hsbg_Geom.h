#ifndef Hsbg_Geom_H
#define Hsbg_Geom_H

#include <iostream>
#include <iomanip>			// IO format control
#include <string>
#include <sstream>

#include "Hsbg_Const.h"
#include "Hsbg_Tools.h"
#include "Hsbg_Point.h"

using namespace std;
using namespace Hsbg;

class HGeom: public GPoint	// defination of GPoint is found at "Hsbg_Global.h"
{
    public:
        string 	Gname; 		// system name, default filename
        int 	Natom;		// the number of atoms
        int 	Nfrag;		// the number of fragments
        int 	iatom;
        int 	ifrag;
        GPoint* geom;
    
    //HGeom();
    
    ~HGeom()
    {
        delete[] this->geom;
    }
        
    int set_Geom(int num)
    {
    	if(num <= 0) {cerr << "warning: set geometry with zero atom!"; exit(-1);}
        this->Natom = num;
        this->iatom = 0;
        this->geom = new GPoint[this->Natom+1];
        return 0;
    }
           
    int read_Geom(string line) 			// read from a string line
    {
        string  buff;
        istringstream istr(line);
        if(istr == NULL)
	    {
	    	cerr << "error: istringstream is null, in ( Hsbg_Geom, HGeom::read_Geom)" << endl;
	    	exit(-1);
	    }
	    
	    this->iatom++;
	    if(this->iatom > this->Natom)
	    {
	    	cerr << "error: the number of atom is overload the prepared" << endl;
	    	exit(-1);
	    }
	    istr >> (this->geom[this->iatom]).aname;
	    istr >> (this->geom[this->iatom]).x;
	    istr >> (this->geom[this->iatom]).y;
	    istr >> (this->geom[this->iatom]).z;
	    (this->geom[this->iatom]).indx = this->iatom;
    }
    
    int solve_Top(string line)			// solve the topology or connection of molecules
    {
	    return 0;
    }
    
    friend ostream &operator<<( ostream &output, HGeom &HG )
    {
        output << endl << "the Geometry of "<< HG.Gname << " shows as follow:" << endl;
        for(int i=1; i <= HG.Natom; i++)
        {
            output << HG.geom[i].get_Point();
        }
        return output;          
    }
};

#endif
