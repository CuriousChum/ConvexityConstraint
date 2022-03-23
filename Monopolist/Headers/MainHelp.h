#pragma once
//
//  MainHelp.h
//  Monopolist
//
//  Created by Jean-Marie Mirebeau on 23/03/2022.
//

#include <iostream>
#include <sstream>

///Display a help message, if required. retval : wether the help was displayed
bool MainHelp(int argc,  const char * argv[], const char * help_msg
			  = "This executable is currently untested and in development.\n"){
	std::cout << "Use --help as first argument to learn about usage.\n";
	if(argc<2 || strcmp(argv[1],"--help")) return false;
	else {std::cout << help_msg; return true;}
}

/// Concatenate the given arguments into a string,
std::string ArgsToString(int argc, const char * argv[]){
	std::ostringstream oss;
	for(int i=0; i<argc; ++i) oss << "_" << argv[i];
	return oss.str();
}
