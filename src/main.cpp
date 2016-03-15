/*
	Scalable Functional Bayesian Association --- MCMC (SFBA:MCMC)
    Copyright (C) 2016  Jingjing Yang

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include "param.h"
#include "sfba.h"

using namespace std;

typedef unsigned char uchar;
typedef unsigned short uint16;
typedef unsigned int uint;
typedef short int int16;


int main(int argc, char * argv[])
{ 	
	SFBA cSfba;	
	PARAM cPar;

	if (argc <= 1) {
		cSfba.PrintHeader(); 
		return EXIT_SUCCESS;
	}
	if (argc==2 && argv[1][0] == '-' && argv[1][1] == 'h') {
		cSfba.PrintHelp(0);
		return EXIT_SUCCESS;
	}
	if (argc==3 && argv[1][0] == '-' && argv[1][1] == 'h') {
		string str;
		str.assign(argv[2]);
		cSfba.PrintHelp(atoi(str.c_str()));
		return EXIT_SUCCESS;
	}
	if (argc==2 && argv[1][0] == '-' && argv[1][1] == 'l') {
		cSfba.PrintLicense();
		return EXIT_SUCCESS;
	}
	
	ifstream check_dir("output/");
	if (!check_dir) {
		mkdir("output", S_IRWXU|S_IRGRP|S_IROTH);
	}	
	
	cSfba.Assign(argc, argv, cPar); 
		
	if (cPar.error==true) {return EXIT_FAILURE;}
	     
	if (cPar.mode_silence) {stringstream ss; cout.rdbuf (ss.rdbuf());}
	
	cPar.CheckParam();
	
	if (cPar.error==true) {return EXIT_FAILURE;}
	
	cSfba.BatchRun(cPar);
	
	if (cPar.error==true) {return EXIT_FAILURE;}
	
	cSfba.WriteLog(argc, argv, cPar);
	
    return EXIT_SUCCESS;                                                          
}


 
