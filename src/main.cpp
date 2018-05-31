#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <string>

#include "Hsbg_Const.h"
#include "Hsbg_Global.h"
#include "Hsbg_Geom.h"
#include "Hsbg_Basis.h"
#include "Hsbg_Parser.h"

using namespace std;
using namespace Hsbg;

int main(int argc, char *argv[])
{
	string gjffile="../test/h2.gjf";
	string logfile="../test/h2.log";
	
	if(argc == 1)
	{
		cout << "need an argument!" << endl;
		cout << "if you want try default settings, you can use argument -d" << endl;
		return -1;
	}
	if(argc == 2 && string(argv[1])!="-d")
	{
		cout << "Unrecognized arguments" << endl;
		return -1;
	}
	if(argc%2==0 && string(argv[1])!="-d")
	{
		cout << "The number of arguments should be pairs" << endl;
		return -1;
	}
	if(argc%2==1)
	{
		for(int i=0; i < argc/2; i++)
		{
			if(string(argv[2*i+1])=="-f") // set IOfiles, can not use switch
			{
				gjffile=argv[2*i+2];
				logfile=gjffile;
				replace_distinct(logfile,".gjf",".log");
				if(logfile==gjffile)
				{
					cout << "The input format (suffix) must be .gjf" << endl;
					return -1;
				}
			}
		}
	}
	
	//cout << gjffile << "\t"<<logfile << endl;
	//return 0;
	
	HTask task = HTask();
	task.set_IO(gjffile,logfile);
	task.set_Job("sp","hf","3-21g"); //default setting
	task.read_Task();
	cout << task << endl;
	//task.solver();
	//task.result();
	return 0;
}
