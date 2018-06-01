// Alias: Hsbg_Tasker.h
// this file is the parser of task, dealing with the geometry and basis

#ifndef _Hsbg_Parser_H_
#define _Hsbg_Parser_H_

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include "Hsbg_Const.h"
#include "Hsbg_Global.h"
#include "Hsbg_Geom.h"
#include "Hsbg_Basis.h"

using namespace std;
using namespace Hsbg;

class HTask
{
    public:
    //  IO/Sys settings
        string  Gjffile;
        //string  Chkfile;
        string  Logfile;
        int     Maxmem;
    
    //  task settings
        string  Job;
        string  Method;
        string  Basis;
        //string  others;
    
    //  Input content
        string  Title;
        int     Charge;
        int     Smulti;
    
    //  others
        HGeom   TaskGeom;
        HBasis  TaskBasis;
        int     Natom;
        int     Nmole; // alias Nfrag
        int     Nbasis;
    // 
        //HTasK();
        //~HTask();
    
    int set_IO(string gjffile, string logfile)
    {
        this->Gjffile = gjffile;
        this->Logfile = logfile;
        return 0;
    }
    
    int set_Job(string job, string method, string basis)
    {
        this->Job = job;
        this->Method = method;
        this->Basis = basis;
        return 0;
    }
    
    int read_Predo() // main count the number of atoms
    {
        int cnt=0;
        string line;
	    fstream fin;
	    
	    fin.open(this->Gjffile.data(), ios::in);
	    if(fin == NULL)
	    {
		    cerr << "error: "<< this->Gjffile <<" doesn't exist" << endl;
		    exit(-1);
	    }
	    
	    while(getline(fin,line))
	    {
	    	line=trim(line);
	    	if(line=="\r" || line.size()<=1) continue; // have been a hidden bug    	
		    if(line[0]=='%' || line[0]=='#' || ( (int)'0' <= (int)line[0] && (int)line[0] <= (int)'9' )) continue;
		    cnt++;
		}
	    this->Natom = cnt-1;
	    this->TaskGeom.Natom = cnt-1;
	    this->TaskGeom.set_Geom(this->TaskGeom.Natom);
	    fin.close();
	    this->TaskGeom.name = this->Gjffile;
	    replace_distinct(this->TaskGeom.name, ".gjf", "");
	    return 0;
    }
    
    int read_Task()
    {
        int cnt=0;
        string line;
        string buff;
	    fstream fin;
	    
	    this->read_Predo();             // pre-do things, relating parser of geometry
	    fin.open(this->Gjffile.data(), ios::in);
	    
	    while(getline(fin,line))
	    {
	    	line=trim(line);
	    	if(line=="\r" || line.size()<=1) continue;	    	
		    if(line[0]=='%') continue;  // here we don't set chkfile
		    if(line[0]=='#')            // main task settings
		    {
		    	line = line.replace(line.find("#"), 1, " ");
		    	line = line.replace(line.find("/"), 1, " ");
		        istringstream istr(line); // error ? not define?
		        while(istr >> buff)
		        {
		            this->Taskparser(buff);
		        }
		        continue;
		    }
		    if(cnt==0)                  // Title of the task
		    {
                this->Title = line;
                cnt++;
		    }
		    else if(cnt==1)             // Charge and Spin-Multiplicity
		    {
		        istringstream istr(line);
		        istr >> this->Charge;
		        istr >> this->Smulti;
		        cnt++;
		    }
		    else if(cnt==2)             // input atoms and its coordinations
		    {
                if( (int)'0' <= (int)line[0] && (int)line[0] <= (int)'9' ) 
                {
                	cnt++;
                }
                else
                {
                	this->TaskGeom.read_Geom(line);
               	}
		    }
		    else if(cnt==3)             // deal with topology, now is unused
		    {
		        this->TaskGeom.solve_Top(line); 
		    }
	    }	    
	    fin.close();
	    
	    this->read_Postdo();
	    return 0;
    }
    
    int read_Postdo()                   // post-do things, relating parser of basis, a little complex
    {        
        this->TaskBasis.name = this->Basis;
        this->TaskBasis.Natom = this->TaskGeom.Natom;
        this->TaskBasis.set_Basis(this->TaskGeom.Natom);
        
        for(int i=1; i <= this->TaskGeom.Natom ; i++)
        {
        	this->TaskBasis.basis[i].setfrom_GPoint(this->TaskGeom.geom[i]);
        }        
        
        for(int i=1; i<= this->TaskBasis.Natom; i++)
        {
        	this->TaskBasis.basis[i].find_ZnumbrPeroid( this->TaskBasis.basis[i].name );
        	this->TaskBasis.basis[i].set_Cloud( this->Basis, this->TaskBasis.basis[i].name );
        }
            
	    return 0;
    }
    
    int Taskparser(string term)
    {
        string dict_job[5] = {"sp", "opt", "freq", "opt+freq", "scan"};
        string dict_method[7] = {"hf", "rhf", "uhf", "dft", "mp2", "ccsd", "cisd"};
        string dict_basis[3] = {"3-21g", "6-31g", "6-311g"};
        
        for(int i=0; i < getArrayLen(dict_job); i++)
        {
            if(dict_job[i]==term)
            {
                this->Job=term;
                return 1;
            }
        }
        for(int i=0; i < getArrayLen(dict_method); i++)
        {
            if(dict_method[i] == term)
            {
                this->Method = term;
                return 2;
            }
        }
        for(int i=0; i < getArrayLen(dict_basis); i++)
        {
            if(dict_basis[i] == term)
            {
                this->Basis = term;
                return 3;
            }
        }
        cerr << "Warning: illegal parameter of "<< term << ", skipped." << endl << endl;
        return -1;
    }
    
    friend ostream &operator<<( ostream &output, HTask &HT )
    {
        output << "========== task ID: ";
        output << HT.Title ;
        output << " ==========" << endl;
        output << "the input  file : " << HT.Gjffile << endl;
        output << "the output file : " << HT.Logfile << endl;
        output << "the settings    : # " << HT.Job << " "<< HT.Method << "/" << HT.Basis << endl;
        output << "(  Q  ,  S  )   : " << HT.Charge << ",     " << HT.Smulti << endl;
        output << HT.TaskGeom;
        return output;
    }
};


#endif
