#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <stdlib.h>
//#include "Hsbg_unitcovertor"
//#include "Hsbg_Abasis"

using namespace std;

class H_task
{
	public:
		char basis[8];
		char job[4];
		char method[8];
		int charge;
		int spin;
		char chkpath[128];
		char outpath[128];
};

class H_basis
{
	int ncore;
	int nshell;
};

int read_task(char taskfile[16])
{
	string line;
	fstream fin;
	//
	//H_basis hbasis;
	//H_task htask;
	//
	fin.open(taskfile,ios::in);
	if(fin == NULL)
	{
		cerr << "error: taskfile doesn't exist" << endl;
		exit(-1);
	}
	while(getline(fin,line))
	{
		cout << line << endl;
		//cout << line[0] << endl;
		if(line[0]=='%'){
            cout << "% is" << endl;
		}else if(line[0]=='#'){
            cout << "# is" << endl;
		}else{
            cout << "normal" << endl;
		}
	}
	fin.close();
}

int main()
{
    cout << "hello world" << endl;
	read_task("h2.gjf");
}
