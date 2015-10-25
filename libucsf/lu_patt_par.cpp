/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_patt_par.cpp                                               *
*                                                                             *
*  Created    : June 14th 2001                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_patt_par.h>
#include <lu_param_list.h>
using std::string;
using std::ostream;
using std::endl;

MSPatternParameters::MSPatternParameters ( const ParameterList* params ) :

	MSSearchParameters ( params ),

	maxAASubstitutions	( params->getIntValue		( "max_aa_substitutions", 0 ) ),
	maxPeptideHits	( params->getIntValue	( "max_hits", 100000 ) ),

	massType			( params->getStringValue	( "parent_mass_convert", "monoisotopic" ) ),
	monoisotopicFlag	( massType == "monoisotopic" )
{
	string svalue;
	if ( params->getValue ( "regular_expression", svalue ) )	{
		regularExpression = ( svalue == "" ) ? "." : svalue;
	}
}
void MSPatternParameters::printHTML ( ostream& os ) const
{
	MSSearchParameters::printHTML ( os );
	if ( regularExpression != "" && regularExpression != "." ) {
		ParameterList::printHTML ( os, "Max AA Substitutions", maxAASubstitutions );
		ParameterList::printHTML ( os, "Regular Expression", regularExpression );
	}
	os << "<p />" << endl;
}
