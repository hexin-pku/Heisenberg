#ifndef Hsbg_Global_H
#define Hsbg_Global_H

#include <iostream>
#include <iomanip>			// IO format control
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>

#include "Hsbg_Const.h"
#include "Hsbg_Tools.h"

using namespace std;

namespace Hsbg
{
	/* here defines basic class: 
	
		Point, Orbital, Orbital_cgto, Atom, Molecule, System
		relationship as following:
		
		~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		!                                                                         !
		!		[ note: ----- inherit, ===== list, =-=-= inherit & list ]         !
		!                                                                         !
		!		Point ----- Orbital =-=-= Orbital_cgto                            !
		!		   `                        //                                    !
		!			   `                   //                                     !
		!			       `              //                                      !
		!			           `         //                                       !
		!			               `  Atom  =-=-=-=-= Molecule =-=-=-=-= System   !
		!	                                                                      !
		~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	     
			 
	*/
	
	class Point
	{
	  	public:
	  		string name; 	// name of function, atom, molecule or the system
		   	double x;		// position
	   		double y;
	   		double z;
	   	
	   	// constructors and replicators
	   	Point()
		{
			this->x = 0;
			this->y = 0;
			this->z = 0;
		}
	
		Point(const Point &P)
		{
			this->name = P.name;
			this->x = P.x;
			this->y = P.y;
			this->z = P.z;
		}
	
		// methods
		double norm()
		{
			return sqrt( (this->x)*(this->x) + (this->y)*(this->y) + (this->z)*(this->z) );
		}
	
		double norm2()
		{
			return (this->x)*(this->x) + (this->y)*(this->y) + (this->z)*(this->z);
		}
	
		Point ref_Point(Point &A)
		{
			Point P = Point();
			P.x = this->x - A.x;
			P.y = this->y - A.y;
			P.z = this->z - A.z;
			return P;
		}
	
		double ref_norm(Point A)
		{
			return (this->ref_Point(A)).norm();
		}
	
		double ref_norm2(Point A)
		{
			return (this->ref_Point(A)).norm2();
		}
	
		friend double dist_AB(Point A, Point B)
		{
			return sqrt( (A.x-B.x)*(A.x-B.x) + (A.y-B.y)*(A.y-B.y) + (A.z-B.z)*(A.z-B.z) );
		}
	
		friend double dist_AB2(Point A, Point B)
		{
			return (A.x-B.x)*(A.x-B.x) + (A.y-B.y)*(A.y-B.y) + (A.z-B.z)*(A.z-B.z);
		}
	
		//addition
		Point operator+(const Point& Pa)
		{
			Point P = Point();
			P.x = this->x + Pa.x;
			P.y = this->y + Pa.y;
			P.z = this->z + Pa.z;
			P.name = "temp";
			return P;
		};
	
		// subtraction
		Point operator-(const Point& Pa)
		{
			Point P = Point();
			P.x = this->x - Pa.x;
			P.y = this->y - Pa.y;
			P.z = this->z - Pa.z;
			P.name = "temp";
			return P;
		};
		
		// scalar multiplication
		Point operator*(const int& a)
		{
			Point P = Point();
			P.x = this->x*a;
			P.y = this->y*a;
			P.z = this->z*a;
			return  P;
		}
	
		// inner product
		double operator*(const Point& Pa)
		{
			return Pa.x * this->x + Pa.y * this->y + Pa.z * this->z;
		}
	
		friend ostream &operator<<( ostream &output, Point &P )
	   	{
		  	    output << setiosflags(ios::fixed) << setprecision(8) << P.name ;
		  	    output << "         ";
		  	    output << P.x;
		  	    output << "    ";
		  	    output << P.y;
		  	    output << "    ";
		  	    output << P.z << endl;
		   	return output;          
	   	}
	};
	
	class Orbital: public Point
	{
		public:
			/*
			// these attributes are inherited from Point
			string name; 	// name of function, such as sto, gto
			double x;		// position
			double y;
			double z;
			*/
			int L;			// of the type of gauss/slater
			int M;
			int N;
			double alpha;	// alpha of gauss/slater
	
		int set_Alpha(double aa)
		{
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
				P->name = "NULL";
				P->x = this->x;
				P->y = this->y;
				P->z = this->z;
				return *P;
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
	
		int conv_AUnit()
		{
			this->x = this->x / c_bohr2ai;
			this->y = this->y / c_bohr2ai;
			this->z = this->z / c_bohr2ai;
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



	class Orbital_cgto : public Orbital
	{
		public:
			/*
			// these attributes are inherited from Point
			string name; 	// name of function, such as sto, gto
			double x;		// position
			double y;
			double z;
			
			// these attributes are inherited from Orbital
			int L;			// of gauss/slater
			int M;
			int N;
			double alpha;	// alpha of gauss/slater
			*/
			// here we note BohrL = L + M + N;
			int cn;	// contraction number
		
			double* coeffs;
			double* alphas;	// will abandon in the futures
			Orbital* gtos;
	
		//Orbital_cgto(int num);
		~Orbital_cgto()
		{
			delete[] gtos;
			delete[] coeffs;
		}

		int set_Cgto(int num)
		{
			this->cn = num;
			this->coeffs = new double[num];
			this->alphas = new double[num];
			this->gtos = new Orbital[num];
			return 0;
		}
	
		// get coeff and alpha of idx-th of gto
		int get_CA(int idx, double& cc, double& aa)
		{
			if(idx>=this->cn) {cerr << "error, segment in cgto" << endl; exit(-1);}
			cc = this->coeffs[idx];
			aa = this->gtos[idx].alpha;
			return 0;
		}
	
		// overload
		int set_LMN(int ll,int mm, int nn)
		{
			this->L=ll;
			this->M=mm;
			this->N=nn;
			for(int i=0; i<this->cn; i++) this->gtos[i].set_LMN(ll,mm,nn);
			return 0;
		}
	
		// overload
		int set_XYZ(double xx, double yy, double zz)
		{
			this->x=xx;
			this->y=yy;
			this->z=zz;
			for(int i=0; i<this->cn; i++) this->gtos[i].set_XYZ(xx,yy,zz);
			return 0;
		}
	
		// overload
		int set_XYZ(Point P)
		{
			this->x=P.x;
			this->y=P.y;
			this->z=P.z;
			for(int i=0; i<this->cn; i++) this->gtos[i].set_XYZ(P);
			return 0;
		}
	
		// overload
		int conv_AUnit()
		{
			this->x = this->x / c_bohr2ai;
			this->y = this->y / c_bohr2ai;
			this->z = this->z / c_bohr2ai;
			for(int i=0; i<this->cn; i++) (this->gtos[i]).conv_AUnit();
			return 0;
		}
	
		// get Bohr's L
		int get_BohrL()
		{
			return this->L + this->M + this->N;
		}
	
		// overload
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
			return 	pow(2*this->gtos[k].alpha/PI, 0.75) * 
					sqrt(pow(4*this->gtos[k].alpha, L+M+N)/static_cast<double>(total));	
		}
	}; 


	class Atom: public Point
	{	
		public:
			/*
			// these attributes are inherit from Point
			string name; 	// atom name (what if it contents fragment information ?)
			double x;		// atom positions
	   		double y;
	   		double z;
			*/
	
			// elements information
			int 	znum;		// Z-number
			double  mass; 		// Atom mass
			// basis information		
			
			// though public, we don't use the following attributes for inheriting
	   		int 	ncgto;		// the number of contract GTO of an atom ~~~~~~~~~~~~~~~
	   		Orbital_cgto *cgto;                                                       //
								// contract GTO (short for cgto)                      //
								                                                      //
			// elements information plus                                              //
			int 	perd;		// the peroid the atom belong to                      //
			int 	fmly;		// the family the atom belong to                      //
			// geometry information                                                   //
			int 	indx;   	// index in the geometry list                         //
	   		int 	frag; 		// such as 1, 2, 3, ..., stand different fragments~~~~//
	   	
	   	Point& get_Point()
	   	{
	   		Point* P = new Point();
	   		P->name = this->name;
	   		P->x = this->x;
	   		P->y = this->y;
	   		P->z = this->z;
	   		return *P;
	   	}
	   	
	   	int get_Point(Point& P)
	   	{
	   		P.name = this->name;
	   		P.x = this->x;
	   		P.y = this->y;
	   		P.z = this->z;
	   		return 0;
	   	}
	   	
	   	// read atom name and postion from a string line
	   	int read_Atom(string line) 
	   	{
		   	string  buff;
		   	istringstream istr(line);
		   	if(istr == NULL)
			{
				cerr << "error: istringstream is null, in (Atom::read_atom)" << endl;
				exit(-1);
			}
			istr >> this->name;
			istr >> this->x;
			istr >> this->y;
			istr >> this->z;
			this->link_Info(this->name);
			return 0;
	   	} 
	   	
	   	// from atom name (such as H), to find znumber, mass, peroid, family information 
	   	int link_Info(string my_name)
		{
			if(my_name=="" || my_name==" ")
			{
				my_name = this->name;
			}
				for(int i=0; i < MAX_ELEMENT ; i++)
			{
				if( string(elements_label[i]) == my_name )
				{
					this->znum = i;
					this->mass = elements_mass[i];
					if(i<=2)
					{
						this->perd = 1;
						this->fmly = i-0;
					}else if(i<=10){
						this->perd = 2;
						this->fmly = i-2;
					}else if(i<=18){
						this->perd = 3;
						this->fmly = i-10;
					}else if(i<=36){
						this->perd = 4;
						this->fmly = i-18;
					}else if(i<=54){
						this->perd = 5;
						this->fmly = i-36;
					}else if(i<=86){
						this->perd = 6;
						this->fmly = i-54;
					}else if(i<=118){
						this->perd = 7;
						this->fmly = i-86;
					}else{
						this->perd = 0;
						this->fmly = 0;
					}
				}
				if(i==MAX_ELEMENT) this->znum = 0;
			}
			return 0;
		}
	   	
		// total arrangment of setting atoms with basis-list
		int	set_Basis(string my_bname, int mysplit, int mynumcs[], int myidxcs[])
		{
			// set Basis space
			this->set_Basisspace(mysplit, mynumcs, myidxcs);
			
			// constrcut/pre-read basis-dat file:
			fstream fin;
			string  basisfile="../basis/"+my_bname+".dat";
			string 	line;
			fin.open(basisfile.data(), ios::in);
			if(fin==NULL) {cerr << "basis file " << basisfile << " doesn't exist!" << endl; exit(-1);}
			int cnt=0;
			while(getline(fin,line))
	   		{
	   			line=trim(line);
	   			if(line=="\r" || line.size()<=1) continue;
	   			if( line[0] == '*' ) cnt++;
	   			if( cnt < this->znum ) continue;
	   			
	   			// read the basis data
				this->read_Basis(fin);
	   			break;
			}
			fin.close();
			
			for(int i=0; i<this->ncgto; i++)
			{
				this->cgto[i].set_XYZ(this->x, this->y, this->z);
			}
			
			return 0;
		}

		// from basis-index-list (s uch as 6-311g --> [6,3,1,1] ), to allocate CGTO space
		int set_Basisspace(int split, int numcs[], int idxcs[])
		{
			this->ncgto = (this->perd==1? split : 1*1 + (this->perd-2)*4 + split*4);
			this->cgto = new Orbital_cgto[this->ncgto];
		
			if(this->perd == 1)
			{
				for(int i=0; i< split; i++)
				{
					this->cgto[i].set_Cgto(numcs[i+1]);
					this->cgto[i].set_LMN(0,0,0);
				}
			}else{
				this->cgto[0].set_Cgto(numcs[0]);
				this->cgto[0].set_LMN(0,0,0);
				for(int i=1; i < this->perd-1; i++)
				{
					this->cgto[4*i-3].set_Cgto(numcs[0]);
					this->cgto[4*i-3].set_LMN(0,0,0);
					this->cgto[4*i-2].set_Cgto(numcs[0]);
					this->cgto[4*i-2].set_LMN(1,0,0);
					this->cgto[4*i-1].set_Cgto(numcs[0]);
					this->cgto[4*i-1].set_LMN(0,1,0);
					this->cgto[4*i].set_Cgto(numcs[0]);
					this->cgto[4*i].set_LMN(0,0,1);
				}
				for(int i=0; i < split; i++)
				{
					this->cgto[4*(this->perd-1)+4*i-3].set_Cgto(numcs[i+1]);
					this->cgto[4*(this->perd-1)+4*i-3].set_LMN(0,0,0);
					this->cgto[4*(this->perd-1)+4*i-2].set_Cgto(numcs[i+1]);
					this->cgto[4*(this->perd-1)+4*i-2].set_LMN(1,0,0);
					this->cgto[4*(this->perd-1)+4*i-1].set_Cgto(numcs[i+1]);
					this->cgto[4*(this->perd-1)+4*i-1].set_LMN(0,1,0);
					this->cgto[4*(this->perd-1)+4*i].set_Cgto(numcs[i+1]);
					this->cgto[4*(this->perd-1)+4*i].set_LMN(0,0,1);
				}
			}
			return 0;
		}
	
		int read_Basis( fstream& input)
		{
			int ssp=0;		//s type default when ssp=0
			string buff;
		
			// check the elements name
			input >> buff;
	   		if(buff!=this->name)
	   		{
	   			cout << "error, find atom basis error of " << buff << endl;
	   			exit(-1);
	   		}
	   		input >> buff;
	   		
	   		// read basis data
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
		
			// to synchronize alpha from cgto ---> gto
			for(int i=0; i<this->ncgto; i++)
			{
				for(int j=0; j<this->cgto[i].cn; j++)
				{
					this->cgto[i].gtos[j].alpha = this->cgto[i].alphas[j];
				}
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
				(this->cgto[i]).conv_AUnit();
			}
			return 0;
		}

		friend ostream &operator<<( ostream &output, Atom &At )
	  		{
			int ssp=0;
			int cnt=1;
			output << setiosflags(ios::left) << setprecision(8);
			output << "the shell of " << At.name << " :" << endl;
		
			for(int j=0; j < At.ncgto; j++)
			{
				if(At.cgto[j].get_BohrL() ==0){
					output << "    shell of " << cnt << (ssp==0? " S":" SP") << endl;
				}
				for(int k=0; k<At.cgto[j].cn; k++)
				{
					output << "        " << At.cgto[j].coeffs[k];
					output << "    " << At.cgto[j].alphas[k];
					if(ssp!=0)
					{
						output << "    " << At.cgto[j+1].alphas[k];
					}
					output << endl;
				}
				if(j!=0) j+=ssp;
				if(At.perd!=1)ssp=3;
				if(cnt < At.perd) cnt++;
			}
			output << "the AU position is:" << endl << At.get_Point() << endl;
			output << "the Z-number    is:" << "    " << At.znum << endl;
			output << endl;
			return output;
		}	
	};



	class Molecule: public Atom
	{
		public:
			/*
			// these attributes are inherited from Point
			string name; 		// molecule name (what if it contents fragment information ?)
			double x;			// molecule centroid position
			double y;
			double z;
			
			// these attributes are inherited from Atom
			int 	znum;		// molecule Z-number
			double  mass; 		// molecule mass
			*/
			
			int 	Natom;		// the number of atoms
		    int 	iatom;
		    Atom*	atoms;
	};



	class System: public Molecule	// defination of Atom is found at "Hsbg_Global.h"
	{
		public:
			/*
			// these attributes are inherited from Point
			string name; 		// system name (what if it contents fragment information ?)
			double x;			// system centroid position
			double y;
			double z;
			
			// these attributes are inherited from Atom
			int 	znum;		// system Z-number
			double  mass; 		// system mass
			
			// these attributes are inherited from Molecule
			int 	Natom;		// the number of atoms, repite again.
		    int 	iatom;		// Note: the lower level of System is not Molecule, but Atom directly for now!
			Atom*	atoms;
			*/
			
			// molecule setting
		    int		Nmol;		// the number of molecules
		    int		imol;
		    Molecule*	moles;	// after atoms all have been build, then we analyse Molecile by topology analyzation
			
			// basis setting
			string 	Bname;		// basis name, such as 6-311g.
			int		Nbasis;		// the total CGTO basis of the system.
			int 	split;		// the split number of valent-shell of each atom
			int* 	numcs; 		// list of contraction numbers, such as 6-311g => 6,3,1,1
			int* 	idxcs;		// list of sum of contraction numbers, such as 6-311g => 6,9,10,11
			int*	idmap;		// give the an atom first-CGTO index in the All-CGTO space
								// 		or to say, bind the atom-list to basis space.
		//System();
		
		~System()
		{
		    delete[] this->atoms;
		    delete[] this->moles;
		    delete[] this->numcs;
		    delete[] this->idxcs;
		    delete[] this->idmap;
		}
		
		int set_Natom(int num)
		{
			if(num <= 0) {cerr << "warning: set atometry with zero atom!"; exit(-1);}
		    this->Natom = num;
		    this->iatom = 0;
		    this->atoms = new Atom[this->Natom+1];
		    this->idmap = new int[this->Natom+1];
		    return 0;
		}
	
		// read from a string line, iatom record the reading schedule, at i-th atom
		int read_Atom(string line)
		{
			this->iatom++;
			if(this->iatom > this->Natom)
			{
				cerr << "error: the number of atom is overload the prepared" << endl;
				exit(-1);
			}
			this->atoms[this->iatom].read_Atom(line);
			return 0;
		}
		
		// treat with Basis name (string), to get numcs, idxcs, split information
		int set_Basis(string my_bname)
		{
			this->Bname = my_bname;
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
			
			
			// set the basis space for each atoms
			for(int i=1; i<=this->Natom; i++)
			{
				this->atoms[i].set_Basis(my_bname, this->split, this->numcs, this->idxcs);
			}
			return 0;
		}
	
		// set up mapping relations from atom-list to cgto-list, deal with idmap
		// you should first run set_Basis, then set_Map!
		int set_Map()
		{
			int cnt=0;
			for(int i=1; i <= this->Natom; i++)
			{
				this->idmap[i]=cnt;
				for(int j=0; j<this->atoms[i].ncgto; j++) cnt++;
			}
			this->Nbasis = cnt;
			return 0;
		}
		
		Orbital_cgto& OrbC(int &idx)
		{
			for(int i=1; i < this->Natom; i++)
			{
				if(this->idmap[i] <= idx && this->idmap[i+1]>idx)
				{
					return this->atoms[i].cgto[idx-this->idmap[i]];
				}
			}
			if(idx >= this->Nbasis) { cerr << "overflow the basis size" << endl; exit(-1) ;}
			return this->atoms[this->Natom].cgto[idx - this->idmap[this->Natom]];
		}
		
		// count the total znumber of the system. (relating to the total electron number)
		int count_Znum()
		{
			int cnt=0;
			for(int i=1; i<=this->Natom; i++)
			{
				cnt+=this->atoms[i].znum;
			}
			this->znum = cnt;
			return 0;
		}
		
		// conversion Atomic unit, from upper later to low layers
		int conv_AUnit()
		{
			this->x = this->x / c_bohr2ai;
			this->y = this->y / c_bohr2ai;
			this->z = this->z / c_bohr2ai;
			for(int i=1; i<=this->Natom; i++)
			{
				(this->atoms[i]).conv_AUnit();
			}
			return 0;
		}
		
		// divide system to seperate molecules/fragments, and anaylize topology information
		int solve_Top(string line)
		{
			return 0;
		}
		
		
		friend ostream &operator<<( ostream &output, System &S )
		{
		    output << endl << "the Geometry of "<< S.name << " shows as follow:" << endl;
		    for(int i=1; i <= S.Natom; i++)
		    {
		        output << S.atoms[i].get_Point();
		    }
			output << endl;
		   	output << endl << "the Basis of "<< S.name << " shows as follow:" << endl;
		    for(int i=1; i <= S.Natom; i++)
		    {
		        output << S.atoms[i];
		    }
		    return output;
		}
		
		Orbital_cgto& operator[](int &idx)
		{
			return this->OrbC(idx);
		}
	};

};


#endif




