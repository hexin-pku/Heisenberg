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
    class Point
    {
    	public:
    		string aname; 	// atom name (what if it contents fragment information ?)
	    	double x;		// atom positions
    		double y;
    		double z;
    		
    	Point()
		{
			this->x = 0;
			this->y = 0;
			this->z = 0;
		}
		
		double dist_A(Point A)
		{
			return sqrt( (A.x-this->x)*(A.x-this->x) + (A.y-this->y)*(A.y-this->y) + (A.z-this->z)*(A.z-this->z) );
		}
		
		double dist_A2(Point A)
		{
			return (A.x-this->x)*(A.x-this->x) + (A.y-this->y)*(A.y-this->y) + (A.z-this->z)*(A.z-this->z);
		}
		
		friend double dist_AB2(Point A, Point B)
		{
			return (A.x-B.x)*(A.x-B.x) + (A.y-B.y)*(A.y-B.y) + (A.z-B.z)*(A.z-B.z);
		}
		
		friend double dist_AB(Point A, Point B)
		{
			return sqrt( (A.x-B.x)*(A.x-B.x) + (A.y-B.y)*(A.y-B.y) + (A.z-B.z)*(A.z-B.z) );
		}
		
		friend ostream &operator<<( ostream &output, Point &P )
    	{
       	    output << setiosflags(ios::fixed) << setprecision(8) << P.aname ;
       	    output << "         ";
       	    output << P.x;
       	    output << "    ";
       	    output << P.y;
       	    output << "    ";
       	    output << P.z << endl;
        	return output;          
    	}
    };
    
    class GPoint: public Point
	{
		public:
			// elements infofmation
	    	int 	znum;		// Z-number
	    	double  mass; 		// Atom mass
	    	int 	perd;		// the peroid the atom belong to
	    	int 	fmly;		// the family the atom belong to
	    	
	    	// geometry information
	    	int 	indx;   	// index in the geometry list
    		int 	frag; 		// such as 1, 2, 3, ..., stand different fragments
    	
    	Point& get_Point()
    	{
    		Point* P = new Point();
    		P->aname = this->aname;
    		P->x = this->x;
    		P->y = this->y;
    		P->z = this->z;
    		return *P;
    	}
    	
    	int get_Point(Point& P)
    	{
    		P.aname = this->aname;
    		P.x = this->x;
    		P.y = this->y;
    		P.z = this->z;
    		return 0;
    	}
    	
    	int link_Info(string my_name)
		{
			if(my_name=="" || my_name==" ")
			{
				my_name = this->aname;
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
	};
    
    

};


#endif




