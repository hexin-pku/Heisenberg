#include <iostream>
#include <string>
//#include <vector>
#include <fstream>
//#include <stdlib.h>

#include "Hsbg_Const.h"
#include "Hsbg_Global.h"
#include "Hsbg_Geom.h"

HGeom* Read_Geom(HTask task, string geominput)
{
	string line;
	fstream fin;
	int cnt;
    //
	//
	fin.open(task.Gjffile, ios::in);
	if(fin == NULL)
	{
		cerr << "error: "<< geomfile <<" doesn't exist" << endl;
		exit(-1);
	}
	cnt=0;
	while(getline(fin,line))
	{
		//cout << line << endl;
		//cout << line[0] << endl;
		if(line[0]=='!') continue;
		if(line[0]=='*') cnt++; //
		if(cnt==AZ)
		{
            cout << line << endl;
		}
		if(cnt>AZ) break;
	}
	fin.close();
	
	
	
}


