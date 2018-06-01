#ifndef _Hsbg_Global_H_
#define _Hsbg_Global_H_

#include <string>
#include <sstream>
#include <fstream>

#include "Hsbg_Const.h"

using namespace std;


template <class T>
int getArrayLen(T& array)
{
return (sizeof(array) / sizeof(array[0]));
};

string& trim(string &s)   
{  
    if (s.empty())   
    {  
        return s;  
    }  
    s.erase(0,s.find_first_not_of(" "));
    s.erase(s.find_last_not_of("\r") + 1);  
    s.erase(s.find_last_not_of(" ") + 1);
    return s;  
}

string& replace_recursive(string& str, const string& old_value, const string& new_value)     
{
    while(true)   {
        string::size_type pos(0);     
        if(   (pos=str.find(old_value))!=string::npos   )     
            str.replace(pos,old_value.length(),new_value);     
        else break;     
    }     
    return str;     
}

string& replace_distinct(string& str, const string& old_value, const string& new_value)     
{
    for(string::size_type pos(0); pos!=string::npos; pos+=new_value.length())
    {     
        if( (pos=str.find(old_value,pos))!=string::npos )     
            str.replace(pos,old_value.length(),new_value);     
        else break;     
    }     
    return str;     
}


namespace Hsbg
{
    class Point
    {
    	public:
	    	double x;
    		double y;
    		double z;
    	Point()
		{
			this->x = 0;
			this->y = 0;
			this->z = 0;
		}
    };
    
    class GPoint: public Point
	{
		public:
	    	string name; // atom name (neccessary with fragment ?)
	    	int znum;
	    	int indx;    // 
    		int mass; // mass
    		int mole; // means fragments
	};
    
    class Cloud
	{
		private:// conflict with GPoint.name??? must be private.
			string namecd;
		public:
			int znumbr; // remove argumrnt flag, by indx
			int peroid;
			int numcs[4]; // core and shell
			int idxcs[4];
			double **c;
			double **s;
		
		Cloud()
		{
			for(int i=0; i<4; i++) this->numcs[i] = 0;
		}
		
		~Cloud()
		{
			for(int i=0; i< getArrayLen(this->s) ; i++) delete this->s[i];
			delete[] this->s;
			if(this->c!=NULL)
			{
				for(int i=0; i< getArrayLen(this->c) ; i++) delete this->c[i];
				delete[] this->c;
			}
		}
		
		int set_Cloud(string basisname, string atomname)
		{		
			this->find_ZnumbrPeroid(atomname);
			
			this->namecd = atomname;			
			for(int i = 0, j=0; i < basisname.size(); i++)
        	{
        		if( (int)'0' <= (int)basisname[i] && (int)basisname[i] <= (int)'9' )
        		{
        			this->numcs[j] = (int)basisname[i] - (int)'0';
        			j++;
        		}
        	}
        	if(this->znumbr <= 2 && this->znumbr > 0) this->numcs[0]=0;
			this->idxcs[0]=this->numcs[0];
			for(int i=1;i<4;i++) this->idxcs[i] = this->idxcs[i-1] + this->numcs[i];
			
			// malloc, throw sometimes only need 2, but for in all, we all used 3!
			if(this->peroid ==1)
			{
				this->c = NULL;
				this->s = new double*[this->idxcs[3]];
				for(int i=0; i < this->idxcs[3]; i++) this->s[i] = new double[3];	
			}
			else
			{
				this->c = new double*[ this->idxcs[0] * (this->peroid-1) ];
				for(int i=1; i < this->idxcs[0]; i++) this->c[i] = new double[3];
				for(int i=this->idxcs[0]; i < this->idxcs[0]*(this->peroid-1); i++) this->c[i] = new double[3];
				
				this->s = new double*[ this->idxcs[3] - this->idxcs[0] ];
				for(int i=0; i < this->idxcs[3] - this->idxcs[0] ; i++) this->s[i] = new double[3];
			}
			
			// read basis-dat file:
			fstream fin;
			string  basisfile="../basis/"+basisname+".dat";
			string 	line;
			fin.open(basisfile.data(), ios::in);
			int cnt=0;
			int ssp; // tmp variable note s/sp, 2,3
			cout << "i jump in" << endl;
			while(getline(fin,line))
	    	{
	    		line=trim(line);
	    		if(line=="\r" || line.size()<=1) continue;
	    		if( line[0] == '*' ) cnt++;
	    		if( cnt < this->znumbr ) continue;
				
				while(getline(fin,line) && trim(line)=="\r") continue;	
	    		istringstream istr(line); string buff; istr>>buff;
	    		if(buff!=atomname) 
	    		{
	    			cout << "error, find atom basis error" << endl;
	    			exit(-1);
	    		}
	    		cout << buff << endl;
	    		
	    		// init core
	    		if(this->numcs[0])
	    		{
	    			while(getline(fin,line) && trim(line)=="\r") continue;
	    			for(int j=0;j < this->numcs[0];j++)
	    			{
	    				while(getline(fin,line) && trim(line)=="\r") continue;
	    				line=trim(line);istringstream istr(line);
	    				//cout << line << endl;
	    				istr >> this->c[j][0] >> this->c[j][1];
	    			}
	    			for(int j=this->numcs[0]; j<this->numcs[0]*(this->peroid-1); j++)
	    			{
	    				if(j%this->numcs[0]==0) getline(fin,line);
	    				while(getline(fin,line) && trim(line)=="\r") continue;
	    				line=trim(line);istringstream istr(line);
	    				//cout << line << endl;
	    				istr >> this->c[j][0] >> this->c[j][1] >> this->c[j][2];
	    			}
	    		}
	    		
	    		// init shell
	    		if(this->numcs[0]==0){
	    			ssp=2;
	    		}else{
	    			ssp=3;
	    		}
	    		
	    		for(int j=0;j < (this->idxcs[3] - this->idxcs[0]);j++)
	    		{
	    			if( j==0 || j == (this->idxcs[1] - this->idxcs[0]) ||
	    				j == (this->idxcs[2] - this->idxcs[0]) )
	    			{
	    				while(getline(fin,line) && trim(line)=="\r") continue;
	    			}
	    			while(getline(fin,line) && trim(line)=="\r") continue;
	    			line=trim(line);istringstream istr(line);
	    			//cout << line << endl;
	    			for(int k=0;k<ssp;k++) istr >> this->s[j][k];
	    		}	    		
	    		break;
			}
			fin.close();
			return 0;
		}
		
		int find_ZnumbrPeroid(string tag) // return peroid
		{
			for(int i=0; i < MAX_ELEMENT ; i++)
			{
				if( string(elements_label[i]) == tag )
				{
					this->znumbr = i;
					if(i<=2)
					{
						this->peroid = 1;
					}else if(i<=10){
						this->peroid = 2;
					}else if(i<=18){
						this->peroid = 3;
					}else if(i<=36){
						this->peroid = 4;
					}else if(i<=54){
						this->peroid = 5;
					}else if(i<=86){
						this->peroid = 6;
					}else if(i<=118){
						this->peroid = 7;
					}else{
						this->peroid = 0;
					}
					return i;
				}
				if(i==MAX_ELEMENT) this->znumbr = 0;
			}
			return 0;
		}
		
		int show_Shells()
    	{
			int ssp;
			cout << "the shell of " << this->namecd << " :" << endl;
			// core
			if((this->numcs[0])!=0)
			{
				for(int j=0; j<this->numcs[0]*(this->peroid-1); j++)
				{
					if(j<this->numcs[0]){ssp=2;} else{ssp=3;}
					if(j%this->numcs[0] == 0)
					{
						cout << "shell " << j/this->numcs[0]+1 << "    type  ";
						if(ssp==2) {
							cout<<"S";
						} else {
							cout<<"SP";
						}
						cout << endl;
					}
				
					for(int k=0;k<ssp;k++) cout << "	" <<this->c[j][k];
					cout << endl;
				}
			}
			// outer-shell
			if(this->numcs[0]==0){ssp=2;} else{ssp=3;}
			for(int j=0;j < (this->numcs[3] - this->numcs[0]); j++)
			{
				for(int k=0; k<4; k++)
				{
					if(j== (this->idxcs[k]-this->idxcs[0]) )
					{
						cout << "shell " << this->peroid << "    type  ";
						if(ssp==2) {
							cout<<"S";
						} else {
							cout<<"SP";
						}
						cout <<"   split: " << k << endl;
					}
				}
				
				for(int k; k<ssp; k++) cout << "	" <<this->s[j][k];
				cout << endl;
			}
			return 0;	
		}
		
		friend ostream &operator<<( ostream &output, Cloud &CD )
    	{
			int ssp;
			output << "the shell of " << CD.namecd << " :" << endl;
			// core
			if((CD.numcs[0])!=0)
			{
				for(int j=0; j<CD.numcs[0]*(CD.peroid-1); j++)
				{
					if(j<CD.numcs[0]){ssp=2;} else{ssp=3;}
					if(j%CD.numcs[0] == 0)
					{
						output << "shell " << j/CD.numcs[0]+1 << "    type  ";
						if(ssp==2) {
							output<<"S";
						} else {
							output<<"SP";
						}
						output << endl;
					}
				
					for(int k=0;k<ssp;k++) output << "	" <<CD.c[j][k];
					output << endl;
				}
			}
			// outer-shell
			if(CD.numcs[0]==0){ssp=2;} else{ssp=3;}
			//cout << "n3 - n0 : " << (CD.idxcs[3] - CD.idxcs[0]) << endl;
			for(int j=0;j < (CD.idxcs[3] - CD.idxcs[0]); j++)
			{
				//cout << "test  " << j << endl;
				for(int k=0; k<4; k++)
				{
					if(j== (CD.idxcs[k]-CD.idxcs[0]) )
					{
						output << "shell " << CD.peroid << "    type  ";
						if(ssp==2) {
							output<<"S";
						} else {
							output<<"SP";
						}
						output <<"   split: " << k << endl;
					}
				}
				
				for(int k=0; k<ssp; k++) output << "	" <<CD.s[j][k];
				output << endl;
			}
			output << endl;
			return output;	
		}
	};
	
	class BPoint : public GPoint, public Cloud
	{
		public:
			string type; // gto or sto
		BPoint()
		{
			this->x = 0;
			this->y = 0;
			this->z = 0;
		};
		int setfrom_GPoint(GPoint GP)
		{
			this->x = GP.x;
			this->y = GP.y;
			this->z = GP.z;
			this->name = GP.name; // atom name (neccessary with fragment ?)
    		this->znum = GP.znum; // charge
    		this->mass = GP.mass; // mass
    		this->mole = GP.mole; // means fragments
    		return 0;
		}
	};
};


#endif




