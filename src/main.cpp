#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <string>

#include "Hsbg_Const.h"
#include "Hsbg_Tools.h"
#include "Hsbg_Global.h"
#include "Hsbg_Tasker.h"
#include "Hsbg_Integral_GTO.h"
#include "Hsbg_SCF.h"

using namespace std;
using namespace Hsbg;

int main(int argc, char *argv[])
{
	string hiffile="../test/HeH.hif";
	string logfile="../test/HeH.hof";
	
	if(argc == 1)
	{
		cout << "Need an argument!" << endl;
		cout << "You can try to use argument -d, this will gives you an example "
			 << "of HeH+ calculation" << endl;
		return -1;
	}
	if(argc == 2 && string(argv[1])!="-d" && string(argv[1])!="-h" )
	{
		cout << "Unrecognized arguments" << endl;
		return -1;
	}
	if(argc == 2 && string(argv[1]) =="-h")
	{
		cout << "Usage:" << endl;
		cout << "    -d                |   run default" << endl;
		cout << "    -f   *.hif/*.gjf  |   run a *.hif/*.gjf file" << endl;
		cout << "    -h                |   help" << endl;
		return -1;
	}
	if(argc%2==0 && string(argv[1])!="-d")
	{
		cout << "The number of arguments should be pairs" << endl;
		return -1;
	}
	if(argc%2==1 && (string(argv[1])=="-d" || string(argv[1])=="-h") )
	{
		cout << "Arguments number error, see by -h" << endl;
	}
	if(argc%2==1)
	{
		for(int i=0; i < argc/2; i++)
		{
			if(string(argv[2*i+1])=="-f")
			{
				hiffile=argv[2*i+2];
				logfile=hiffile;
				replace_distinct(logfile,".hif",".hof");
				if(logfile==hiffile)
				{
					replace_distinct(logfile,".gjf",".log");
				}else{
					continue;
				}
				if(logfile==hiffile)
				{
					cout << "The input format (suffix) must be .hif/.gjf" << endl;
					return -1;
				}
			}
		}
	}

	Tasker task = Tasker();						// build a task	
	task.set_IO(hiffile, logfile);				// set IO settings
	task.set_Job("sp","hf","3-21g"); 			// default setting, if the hiffile doesn't give settings
	task.read_Task();							// read from *.hif file
	
	cout << task ;								// show detials of the task
	SCFer my_scf = SCFer(task.Title);		// build an object of SCF
	my_scf.set_Space(task);						// set SCF basisset space size
	my_scf.set_Threshold(0.00000001);			// set SCF threshold
	my_scf.loop_SCF();							// do SCF job
	return 0;

}



