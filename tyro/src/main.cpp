#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <string>

#include "basic_Const.h"
#include "basic_Tools.h"
#include "basic_Global.h"
#include "basic_Parser.h"
#include "basic_Integral.h"
#include "basic_SCF.h"

using namespace std;

int main(int argc, char *argv[])
{
	string hiffile="../test/HeH.gjf";
	string logfile="../test/HeH.log";
	
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
		cout << "    -f   *.hif/*.gjf  |   run a *.gjf file" << endl;
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
				replace_distinct(logfile,".gjf",".log");
				if(logfile==hiffile)
				{
					cout << "The input format (suffix) must be .gjf" << endl;
					return -1;
				}
			}
		}
	}

	Task task = Task();							// build a task	
	task.set_IO(hiffile, logfile);				// set IO settings
	task.set_Job("sp","hf","3-21g"); 			// default setting, if the hiffile doesn't give settings
	task.read_Task();							// read from *.hif file
	
	cout << task ;								// show detials of the task
	SCFer my_scf = SCFer(task.Title);				// build an object of SCF
	my_scf.set_Space(task);						// set SCF basisset space size
	my_scf.set_Threshold(0.00000001);			// set SCF threshold
	my_scf.loop_SCF();							// do SCF job
	return 0;

}



