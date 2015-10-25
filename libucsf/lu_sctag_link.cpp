/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_sctag_link.cpp                                             *
*                                                                             *
*  Created    : September 9th 2005                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2005-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef BATCHTAG
#include <cstdlib>
#include <lu_srch_form.h>
#include <lu_spec_id.h>
#include <lu_sctag_link.h>
#include <lu_cgi_val.h>
#include <lu_param_list.h>
#include <lu_prog.h>

using std::string;
using std::ostream;

void SCMSTagLink::putCGI ( ostream& os ) const
{
	const ParameterList* params = ProgramLink::getParams ();
	printCGIString ( os, "form", "mstag" );
	if ( params->getStringValue ( "msms_pk_filter" ) != "Unprocessed MSMS" ) {
		params->copyToCGI ( os, "msms_pk_filter" );
	}
	params->copyToCGI ( os, "msms_max_peaks" );
	params->copyToCGI ( os, "msms_max_reported_hits" );
// Start crosslinking parameters
	string linkSearchType = params->getStringValue ( "link_search_type" );
	if ( linkSearchType != "No Link" ) {
		params->copyToCGI ( os, "link_search_type" );
		params->copyToCGI ( os, "max_saved_tag_hits" );
	}
	if ( linkSearchType == "User Defined Link" ) {
		params->copyToCGI ( os, "link_aa" );
		for ( StringVectorSizeType i = 0 ; i < CrosslinkingForm::getNumUserMods () ; i++ ) {
			string sNum = gen_itoa ( i + 1 );
			params->copyToCGI ( os, "mod_" + sNum + "_label" + sNum );
			params->copyToCGI ( os, "aa_modified_" + sNum );
			params->copyToCGI ( os, "mod_" + sNum + "_composition" );
		}
	}
// End crosslinking parameters
}
void SCMSTagLink::write ( ostream& os, const string& searchKey, const SpecID& specID, const string& str ) const
{
	ProgramLink::openLink ( os, "scMSTagLink", -1 );
	printCGIString ( os, "search_key", searchKey );
	specID.putCGI ( os );
	os << "\\\">";
	os << str;
	ProgramLink::closeLink ( os );
}
void SCMSTagLink::printHTML ( ostream& os ) const
{
	os << "scMSTagLink" << "=\"";
	os << ProgramLink::getURLStart ( "msform" );
	os << "?";
	putCGI ( os );
	os << "\";\n";
}
#endif
