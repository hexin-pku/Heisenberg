#ifndef _Hsbg_Geom_H_
#define _Hsbg_Geom_H_

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

#include "Hsbg_Const.h"
#include "Hsbg_Global.h"

using namespace std;

typedef struct Gnode
{
    string name; // atom name (neccessary with fragment ?)
    int indx;
    int znum;
    int mass;
    int mole; // means fragments
    double x;
    double y;
    double z;
    Gnode* next; // for solving topology problem
} *pGnode;

class HGeom
{
    public:
        string name; // system name, might same with filename
        int Natom;
        int Nmole;
        int Nbasis;
        int iatom;
        int imole;
        int ibasis;
        pGnode geom;
        
    HGeom()
    {
        int Natom=0;
        int Nmole=0;
        int Nbasis=0;
        int iatom=0;
        int imole=0;
        int ibasis=0;
        pGnode geom = NULL;
    }
    
    ~HGeom()
    {
        delete[] this->geom;
    }
    
    int set_Geom()
    {
        if(this->Natom==0) cerr << "warning: set geometry with zero atom!";
        this->geom = new Gnode[this->Natom+1]; //?
    }
    
    int set_Geom(int num)
    {
        if(this->Natom==0) cerr << "warning: set geometry with zero atom!";
        this->geom = new Gnode[num+1];
    }
        
    int read_Geom(string line) // called by Hsbg_Parser, read a line
    {
        string  buff;
        istringstream istr(line);
        if(istr == NULL)
	    {
	    	cerr << "error: istringstream is null, in ( Hsbg_Geom, HGeom::read_Geom)" << endl;
	    	exit(-1);
	    }
	    this->iatom++;
	    istr >> (this->geom[this->iatom]).name;
	    istr >> (this->geom[this->iatom]).x;
	    istr >> (this->geom[this->iatom]).y;
	    istr >> (this->geom[this->iatom]).z;
	    (this->geom[this->iatom]).indx = this->iatom;
    }
    
    int solve_Top(string line)
    {
	    return 0;
    }
    
    friend ostream &operator<<( ostream &output, HGeom &HG )
    {
        output << endl << "the Geometry of "<< HG.name << " shows as follow:" << endl;
        for(int i=1; i <= HG.Natom; i++)
        {
            output << setiosflags(ios::fixed) << setprecision(8) << HG.geom[i].name ;
            output << "         ";
            output << HG.geom[i].x;
            output << "    ";
            output << HG.geom[i].y;
            output << "    ";
            output << HG.geom[i].z << endl;
        }
        return output;          
    }
};

#endif
