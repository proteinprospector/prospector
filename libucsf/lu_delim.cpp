/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_delim.cpp                                                  *
*                                                                             *
*  Created    : March 25rd 2003                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2003-2008) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lg_io.h>
using std::ostream;
using std::string;

static bool start;
static char delimiter = '\t';

static void printDelimiter ( ostream& os )
{
	if ( start ) start = false;
	else os << delimiter;
}
void delimitedRowStart ( ostream& os )
{
	start = true;
}
void delimitedRowEnd ( ostream& os )
{
	os << '\n';
}
void delimitedHeader ( ostream& os, const string& str )
{
	printDelimiter ( os );
	os << str;
}
void delimitedCellSigFig ( ostream& os, double number, int sigFig )
{
	printDelimiter ( os );
	genPrintSigFig ( os, number, sigFig );
}
void delimitedCell ( ostream& os, double number, int precision )
{
	printDelimiter ( os );
	genPrint ( os, number, precision );
}
void delimitedCellRange ( ostream& os, double number1, double number2, int precision )
{
	printDelimiter ( os );
	genPrint ( os, number1, precision );
	os << " - ";
	genPrint ( os, number2, precision );
}
void delimitedCell ( ostream& os, const string& str )
{
	printDelimiter ( os );
	os << str;
}
void delimitedCell ( ostream& os, char number )
{
	printDelimiter ( os );
	os << number;
}
void delimitedCell ( ostream& os, int number )
{
	printDelimiter ( os );
	os << number;
}
void delimitedEmptyCell ( ostream& os )
{
	printDelimiter ( os );
}
void delimitedEmptyNCells ( ostream& os, int num )
{
	for ( int i = 0 ; i < num ; i++ ) {
		printDelimiter ( os );
	}
}
