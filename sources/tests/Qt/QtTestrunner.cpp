/**
 * \file QtTestrunner.cpp
 * 2012-01-16 LB Initial implementation
 */

#include "gtest.h"
#include <string>

#include <QApplication>

// An evil global variable
bool g_writeRef = false;
bool g_writeOutput = false;

int main(int argc, char* argv[])
{
	if(argc > 1)
	{
		std::string argument1(argv[1]);
		if(argument1.compare("update") == 0)
			g_writeRef = true;
		else if(argument1.compare("output") == 0)
			g_writeOutput = true;
	}
	
	testing::InitGoogleTest ( &argc, argv );
	
	// Creating a non-gui (console) Qt application
	QApplication app(argc, argv, false);
	
	return RUN_ALL_TESTS();
}