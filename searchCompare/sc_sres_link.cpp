/******************************************************************************
*                                                                             *
*  Program    : searchCompare                                                 *
*                                                                             *
*  Filename   : sc_sres_link.cpp                                              *
*                                                                             *
*  Created    : July 15th 2002                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2003-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif
#include <lu_cgi_val.h>
#include <lu_param_list.h>
#include <sc_search_res.h>
#include <sc_sres_link.h>
using std::ostream;
using std::string;

SResLink::SResLink ( bool rawForwarding, bool xlink ) :
	urlProg ( rawForwarding ? "searchCompareRawData" : "searchCompare" ),
	xlink ( xlink )
{
}
void SResLink::write ( ostream& os, const string& accessionNumber, const string& id, const string& rank ) const
{
	StringVector sv;
	sv.push_back ( accessionNumber );
	ProgramLink::openLink ( os, "sresLink", -1 );
	printCGIContainer ( os, "accession_nums", sv );
	if ( !id.empty () && id != SearchResults::getDefaultID () ) {
		StringVector sv2;
		sv2.push_back ( id );
		printCGIContainer ( os, "id_filter_list", sv2 );
	}
	else {
		printCGIContainer ( os, "id_filter_list", idFilterList );
	}
	os << "\\\">";
	os << rank;
	ProgramLink::closeLink ( os );
}
void SResLink::putCGI ( ostream& os ) const
{
	ParameterList p = *(ProgramLink::getParams ());
	idFilterList = p.getPQStringVectorValue ( "id_filter_list" );
	p.removeName ( "accession_nums" );
	p.removeName ( "id_filter_list" );
	p.removeName ( "remove" );
	p.removeName ( "preferred_species" );
	p.setValue ( "report_type", "Peptide" );
	//if ( xlink ) {
	//	p.setValue ( "sort_type", "Start Residue" );
	//	p.setValue ( "sort_type_2", "End Residue" );
	//}
	p.copyToCGI ( os );
}
void SResLink::printHTML ( ostream& os ) const
{
	os << "sresLink" << "=\"";
	os << ProgramLink::getURLStart ( urlProg );
	os << "?";
	putCGI ( os );
	os << "\";\n";
}
