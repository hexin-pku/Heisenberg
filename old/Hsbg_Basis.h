#ifndef Hsbg_Basis_H
#define Hsbg_Basis_H

#include <iostream>
#include <iomanip>			// IO format control
#include <sstream>
#include <string>

#include "Hsbg_Const.h"
#include "Hsbg_Tools.h"
#include "Hsbg_Point.h"
#include "Hsbg_Orbital.h"

using namespace std;
using namespace Hsbg;

class BPoint: public GPoint
{
	public:
		int 	ncgto;		// all sum of basis needed
		int 	split;		// the split number of valent-shell
		int* 	numcs; 		// contraction numbers of core and shell, such as 6-311
		int* 	idxcs;	 	// the sum of contraction numbers, for 6311 => 6,9,10,11
							// for H atom, the numcs be 0-311
							//			   the indcs be 0,3,4,5
		HOrbital_cgto *cgto;
						// core and shell
	~BPoint()
	{
		delete[] cgto;
	}
	
	int setfrom_GPoint(GPoint GP)
	{
		this->x = GP.x;
		this->y = GP.y;
		this->z = GP.z;
		this->aname = GP.aname; 	// atom name, attribute of [Point]
		
   		this->znum = GP.znum;
   		this->mass = GP.mass;
   		this->perd = GP.perd;
   		this->fmly = GP.fmly;
   		
   		this->indx = GP.indx;
   		this->frag = GP.frag;
   		return 0;
	}
	
	int set_Space(string my_bname)
	{
		this->split=0;
		for(int i = 0; i < my_bname.size(); i++)
       	{
       		if( (int)'0' <= (int)my_bname[i] && (int)my_bname[i] <= (int)'9' ) this->split ++;
       	}
       	this->split--;
       	this->numcs = new int[this->split+1];
       	this->idxcs = new int[this->split+1];
       	for(int i = 0, j=0; i < my_bname.size(); i++)
       	{
       		if( (int)'0' <= (int)my_bname[i] && (int)my_bname[i] <= (int)'9' )
       		{
       			this->numcs[j] = (int)my_bname[i]-(int)'0';
       			if(j==0)
       			{
       				this->idxcs[j] = this->numcs[j];
       			}else{
       				this->idxcs[j] = this->idxcs[j-1] + this->numcs[j];
       			}
       			j++;
       		}
       		if(j > this->split) break;
       	}
		
		this->ncgto = (this->perd==1? this->split : 1*1 + (this->perd-2)*4 + (this->split)*4);
		this->cgto = new HOrbital_cgto[this->ncgto];
		
		if(this->perd == 1)
		{
			for(int i=0; i<this->split; i++)
			{
				this->cgto[i].set_Cgto(this->numcs[i+1]);
				this->cgto[i].set_LMN(0,0,0);
			}
		}else{
			this->cgto[0].set_Cgto(this->numcs[0]);
			this->cgto[0].set_LMN(0,0,0);
			for(int i=1; i < this->perd-1; i++)
			{
				this->cgto[4*i-3].set_Cgto(this->numcs[0]);
				this->cgto[4*i-3].set_LMN(0,0,0);
				this->cgto[4*i-2].set_Cgto(this->numcs[0]);
				this->cgto[4*i-2].set_LMN(1,0,0);
				this->cgto[4*i-1].set_Cgto(this->numcs[0]);
				this->cgto[4*i-1].set_LMN(0,1,0);
				this->cgto[4*i].set_Cgto(this->numcs[0]);
				this->cgto[4*i].set_LMN(0,0,1);
			}
			for(int i=0; i < this->split; i++)
			{
				this->cgto[4*(this->perd-1)+4*i-3].set_Cgto(this->numcs[i+1]);
				this->cgto[4*(this->perd-1)+4*i-3].set_LMN(0,0,0);
				this->cgto[4*(this->perd-1)+4*i-2].set_Cgto(this->numcs[i+1]);
				this->cgto[4*(this->perd-1)+4*i-2].set_LMN(1,0,0);
				this->cgto[4*(this->perd-1)+4*i-1].set_Cgto(this->numcs[i+1]);
				this->cgto[4*(this->perd-1)+4*i-1].set_LMN(0,1,0);
				this->cgto[4*(this->perd-1)+4*i].set_Cgto(this->numcs[i+1]);
				this->cgto[4*(this->perd-1)+4*i].set_LMN(0,0,1);
			}
		}
		return 0;
	}
	
	int set_BPoint(string my_bname, string my_aname)
	{		
		this->aname = my_aname;
		this->link_Info(my_aname);

		this->set_Space(my_bname);
		// read basis-dat file:
		fstream fin;
		string  basisfile="../basis/"+my_bname+".dat";
		string 	line;
		fin.open(basisfile.data(), ios::in);
		if(fin==NULL) {cerr << "basis file " << basisfile << " doesn't exist!" << endl; exit(-1);}
		int cnt=0;
		int ssp; // tmp variable note s/sp, 2,3
		while(getline(fin,line))
    	{
    		line=trim(line);
    		if(line=="\r" || line.size()<=1) continue;
    		if( line[0] == '*' ) cnt++;
    		if( cnt < this->znum ) continue;
			this->read_BPoint(fin);
    		break;
		}
		fin.close();
		
		for(int i=0; i<this->ncgto; i++)
		{
			this->cgto[i].set_XYZ(this->x, this->y, this->z);
		}
		return 0;
	}
	
	int read_BPoint( fstream& input)
   	{
		int ssp=0;		//s type
		string buff;
		
		// check the elements name
		input >> buff;
    	if(buff!=this->aname)
    	{
    		cout << "error, find atom basis error of " << buff << endl;
    		exit(-1);
    	}
    	input >> buff;
    	
		for(int j=0; j < this->ncgto; j++)
		{
			if(this->cgto[j].get_BohrL() ==0){
				for(int i=0; i<3; i++){ input >> buff;}
			}
			for(int k=0; k < this->cgto[j].cn && k<100; k++)
			{
				input >> this->cgto[j].alphas[k];
				input >> this->cgto[j].coeffs[k];
				if(ssp!=0)
				{
					this->cgto[j+1].alphas[k] = this->cgto[j].alphas[k];
					input >> this->cgto[j+1].coeffs[k];
					this->cgto[j+2].alphas[k] = this->cgto[j].alphas[k];
					this->cgto[j+2].coeffs[k] = this->cgto[j+1].coeffs[k];
					this->cgto[j+3].alphas[k] = this->cgto[j].alphas[k];
					this->cgto[j+3].coeffs[k] = this->cgto[j+1].coeffs[k];
				}
			}
			if(j!=0) j+=ssp;
			if(this->perd!=1)ssp=3;
		}
		return 0;	
	}
	
	int conv_AUnit()
	{
		this->x = this->x / c_bohr2ai;
		this->y = this->y / c_bohr2ai;
		this->z = this->z / c_bohr2ai;
		for(int i=0; i<this->ncgto; i++)
		{
			this->cgto[i].conv_AUnit();
		}
	}
	
	friend ostream &operator<<( ostream &output, BPoint &BP )
   	{
		int ssp=0;
		int cnt=1;
		output << setiosflags(ios::left) << setprecision(8);
		output << "the shell of " << BP.aname << " :" << endl;
		
		for(int j=0; j < BP.ncgto; j++)
		{
			if(BP.cgto[j].get_BohrL() ==0){
				output << "    shell of " << cnt << (ssp==0? " S":" SP") << endl;
			}
			for(int k=0; k<BP.cgto[j].cn; k++)
			{
				output << "        " << BP.cgto[j].coeffs[k];
				output << "    " << BP.cgto[j].alphas[k];
				if(ssp!=0)
				{
					output << "    " << BP.cgto[j+1].alphas[k];
				}
				output << endl;
			}
			if(j!=0) j+=ssp;
			if(BP.perd!=1)ssp=3;
			if(cnt < BP.perd) cnt++;
		}
		output << "the AU location is:" << endl << BP.get_Point() << endl;
		output << endl;
		return output;
	}
	
};

class HBasis: public BPoint
{
    public:
    	string  Bname; 		// gto, sto, etc, not same with the basis
    	int  	Natom;
    	int		Allznum;
        BPoint* basis;
        
        int*	idmap;
        //		idmap[i] returns i-th atom's fisrt cgto index in the all array.
    
    int set_Basis(int num) // num is the atom position number
    {
    	this->Natom = num;
    	this->basis = new BPoint[num+1];
    	this->idmap = new int[num+1];
    	return 0;
    }    
    
    int set_Map()
    {
    	int cnt=0;
    	for(int i=1; i <= this->Natom; i++)
    	{
    		this->idmap[i]=cnt;
    		for(int j=0; j<this->basis[i].ncgto; j++) cnt++;
    	}
    }
    
    int count_Allznum()
    {
    	int cnt=0;
    	for(int i=1; i<=this->Natom; i++)
    	{
    		cnt+=this->basis[i].znum;
    	}
    	this->Allznum = cnt;
    }
    
    friend ostream &operator<<( ostream &output, HBasis &HB )
    {
        for(int i=1; i <= HB.Natom; i++)
        {
        	output << HB.basis[i];
        }
        for(int i=1; i<= HB.Natom; i++) output << "idmap: " << i << " "<< HB.idmap[i] << endl;
        return output;          
    }
};

#endif
