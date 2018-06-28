// Alias: Hsbg_Parser.h
// this file is the parser of task, dealing with the geometry and basis

#ifndef Hsbg_Parser_H
#define Hsbg_Parser_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include "Hsbg_Const.h"
#include "Hsbg_Tools.h"
#include "Hsbg_Global.h"

using namespace std;
using namespace Hsbg;

class Tasker
{
    public:
    	//  IO/Sys settings
        string  Hiffile;
        //string  Chkfile;
        string  Logfile;
        int     Maxmem;
    
    	//  task settings
        string  Job;
        string  Method;
        string  Basis;
        //string  others;// if needed
    
    	//  Input content
        string  Title;
        int     Charge;
        int     Smulti;
    
    	//  others
        System 	Sys;
        int     Natom;		// number of atom
        int		Nelec;		// number of electron
        int		Nbasis;		// number of contract basis
    
    /*    
    // constructor and disconstructor
    //Tasker();
    Tasker(string myhiffile, string mylogfile)
    {
      	this->Hiffile = myhiffile;
       	this->Logfile = mylogfile;
    }
    //~Tasker();
    */
    
    int set_IO(string Hiffile, string logfile)
    {
        this->Hiffile = Hiffile;
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
    
    // allocate atoms
    int read_Predo()
    {
        int cnt=0;
        string line;
	    fstream fin;
	    
	    fin.open(this->Hiffile.data(), ios::in);
	    if(fin == NULL)
	    {
		    cerr << "error: "<< this->Hiffile <<" doesn't exist" << endl;
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
	    this->Sys.set_Natom(this->Natom);
	    fin.close();
	    this->Sys.name = this->Hiffile;
	    replace_distinct(this->Sys.name, ".gjf", "");
	    return 0;
    }
    
    int read_Task()
    {
        int cnt=0;
        string line;
        string buff;
	    fstream fin;
	    
	    // pre-do task 
	    this->read_Predo();
	    fin.open(this->Hiffile.data(), ios::in);
	    
	    while(getline(fin,line))
	    {
	    	line=trim(line);
	    	if(line=="\r" || line.size()<=1) continue;	    	
		    if(line[0]=='%') continue;  // here we don't set chkfile
		    if(line[0]=='#')            // main task settings
		    {
		    	line = line.replace(line.find("#"), 1, " ");
		    	line = line.replace(line.find("/"), 1, " ");
		        istringstream istr(line);
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
                	this->Sys.read_Atom(line);
               	}
		    }
		    else if(cnt==3)             // deal with topology, now is unused
		    {
		        this->Sys.solve_Top(line); 
		    }
	    }	    
	    fin.close();
	    
	    // post-do task
	    this->Sys.Bname = this->Basis;
        this->Sys.set_Basis(this->Basis);
        this->Sys.conv_AUnit();
        this->Sys.set_Map();
        this->Sys.count_Znum();
        this->Nelec = this->Sys.znum - this->Charge;
        this->Nbasis = this->Sys.idmap[this->Natom] + this->Sys.atoms[this->Natom].ncgto;
        if(this->Nbasis!=this->Sys.Nbasis) { cerr << "Nbasis count error!" << endl; exit(-1); }
	    return 0;
    }
    
    int Taskparser(string term)
    {
        string dict_job[] = {"sp", "opt", "freq", "opt+freq", "scan"};
        string dict_method[] = {"hf", "rhf", "uhf", "dft", "mp2", "ccsd", "cisd"};
        string dict_basis[] = {"3-3g","3-21g", "6-31g", "6-311g"};
        
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
    
    friend ostream &operator<<( ostream &output, Tasker &HT )
    {
        output << "========== task ID: ";
        output << HT.Title ;
        output << " ==========" << endl;
        output << "the input  file : " << HT.Hiffile << endl;
        output << "the output file : " << HT.Logfile << endl;
        output << "the settings    : # " << HT.Job << " "<< HT.Method << "/" << HT.Basis << endl;
        output << "(  Q  ,  S  )   : " << HT.Charge << ",     " << HT.Smulti << endl;
        output << "with Natom	   : " << HT.Natom << endl;
        output << "with Nelec	   : " << HT.Nelec << endl;
        output << "with N_set	   : " << HT.Nbasis << endl << endl;
        for(int i=1; i<=HT.Natom; i++) output << "IDmap:  " << i << "  " << HT.Sys.idmap[i] << endl; 
        output << endl << HT.Sys << endl;
        return output;
    }
};


#endif
