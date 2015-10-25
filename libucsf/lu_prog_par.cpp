/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_prog_par.cpp                                               *
*                                                                             *
*  Created    : September 10th 2001                                           *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lg_new.h>
#include <lu_html.h>
#include <lu_prog_par.h>
#include <lu_app_gr.h>
#include <lu_param_list.h>
using std::cout;
using std::endl;

MSProgramParameters::MSProgramParameters ( const ParameterList* params ) :
	searchName			( params->getStringValue( "search_name", "" ) ),
	reportTitle			( params->getStringValue( "report_title", "" ) ),
	outputType			( params->getStringValue( "output_type", "HTML" ) ),
	toFileParameters	( params, searchName ),
	displayGraph		( params->getBoolValue	( "display_graph", false ) )
{
	SpectrumGraph::setDrawGraph ( displayGraph );
}
MSProgramParameters::~MSProgramParameters ()
{
}
void initialiseProspector ()
{
	gen_set_new_handler ();

	cout << "Content-type: text/html" << endl << endl;

	static HTMLErrorHandler h;
}
