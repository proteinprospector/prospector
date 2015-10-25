/******************************************************************************
*                                                                             *
*  Program    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_program.cpp                                                *
*                                                                             *
*  Created    : October 22nd 2001                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lp_frame.h>
#include <lu_html.h>
#include <lu_version.h>
#include <lu_program.h>
#include <lu_xml.h>
#include <lu_param_list.h>
using std::cout;
using std::ostream;
using std::endl;
using std::string;

const ParameterList* MSProgram::paramList = 0;
MSProgram::MSProgram ( const MSProgramParameters& params ) :
	params ( params )
{
	if ( params.getOutputType () != "HTML" && params.getOutputType () != "Tab delimited text" ) {
		FrameIterator::setOutputCharacter ( ' ' );
	}
}
MSProgram::~MSProgram () {}
void MSProgram::outputResults ( bool showParams )
{
	string outputType = params.getOutputType ();
	string searchName = params.getSearchName ();
	if ( outputType == "HTML" ) {
		printAbortFunctions ( cout );
		if ( params.getToFile () ) {
			saveHits ();
			cout << "<p>" << endl;
			params.printOutputFileLink ( cout, "HTML Results", ".htm" );
			cout << "</p>" << endl;
			GenOFStream os ( params.getOutputPath ( ".htm" ) );
			if ( searchName == "msproduct" || searchName == "msisotope" ) {
				init_html ( os, "", "", true );
			}
			else {
				init_html ( os, params.getReportTitle () + " Search Results", "", true );
			}
			printHTML ( os, showParams );
			printProgramInformationHTML ( cout, params.getReportTitle () );
		}
		else
			printHTML ( cout, showParams );
	}
	else if ( outputType == "Tab delimited text" ) {
		if ( params.getToFile () ) {
			printAbortFunctions ( cout );
			saveHits ();
			cout << "<p>" << endl;
			params.printOutputFileLink ( cout, "Tab Delimited Text Results", ".txt" );
			cout << "</p>" << endl;
			GenOFStream os2 ( params.getOutputPath ( ".txt" ) );
			printTabDelimitedText ( os2 );
			printProgramInformationHTML ( cout, params.getReportTitle () );
		}
		else {
			cout << "<pre>" << endl;
				printTabDelimitedText ( cout );
			cout << "</pre>" << endl;
			printProgramInformationHTML ( cout, params.getReportTitle () );
		}
	}
	else if ( outputType == "mgf" ) {
		cout << "<pre>" << endl;
			printMGFSpectrum ( cout );
		cout << "</pre>" << endl;
		printProgramInformationHTML ( cout, params.getReportTitle () );
	}
	else {								// XML
		if ( params.getToFile () ) {
			printAbortFunctions ( cout );
			saveHits ();
			cout << "<p>" << endl;
			params.printOutputFileLink ( cout, "XML Results", ".xml" );
			cout << "</p>" << endl;
			GenOFStream os2 ( params.getOutputPath ( ".xml" ) );
			printXML ( os2 );
			printProgramInformationHTML ( cout, params.getReportTitle () );
		}
		else
			printXML ( cout );
	}
}
void MSProgram::printHTML ( ostream& os, bool showParams )
{
	if ( showParams ) printParamsHTML ( os );
	printBodyHTML ( os );
	printProgramInformationHTML ( os, params.getReportTitle () );
}
void MSProgram::printXML ( ostream& os )
{
	printXMLTop ( os );
	printBodyXML ( os );
	printXMLBottom ( os );
}
void MSProgram::printTabDelimitedText ( ostream& os )
{
	printBodyTabDelimitedText ( os );
}
void MSProgram::printMGFSpectrum ( ostream& os )
{
	printBodyMGFSpectrum ( os );
}
void MSProgram::printXMLTop ( ostream& os )
{
	printXMLHeader ( os );
	printXMLVersion ( os );
	ToFileParameters tfp = params.getToFileParameters ();
	if ( tfp.getScript () ) {
		os << "<?xml-stylesheet type=\"";
		os << tfp.getScriptType ();
		os << "\" href=\"";
		os << tfp.getScriptName ();
		os << "\"?>";
		os << endl;
	}
	os << "<" << params.getSearchName () << "_report>" << endl;
	paramList->XMLParameters ( os );
}
void MSProgram::printXMLBottom ( ostream& os )
{
	os << "</" << params.getSearchName () << "_report>" << endl;
}
void MSProgram::printParamsHTML ( ostream& os ) const
{
	os << "<p>" << endl;
	ExpandableJavascriptBlock ejb ( "Parameters" );
	ejb.printHeader ( os );
	os << "<p>" << endl;
	printParamsBodyHTML ( os );
	os << "</p>" << endl;
	ejb.printFooter ( os );
	os << "</p>" << endl;
}
