/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_tag_par.cpp                                                *
*                                                                             *
*  Created    : July 3rd 2001                                                 *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lg_string.h>
#include <lu_get_link.h>
#include <lu_inst.h>
#include <lu_tag_par.h>
#include <lu_param_list.h>
using std::string;
using std::ostream;
using std::endl;

MSTagParameters::MSTagParameters ( const ParameterList* params, bool noFilesOK ) :
	MSSearchParameters	( params, "msms_" ),

	msmsPeakFilterOptions	( new MSMSPeakFilterOptions ( params ) ),

	specID				( params ),
	dataSetInfo			( new MSMSDataSetInfo ( params, noFilesOK ) ),

	parentMassTolerance	( "msms_parent_mass", params ),
	productMassTolerance( "fragment_masses", params ),

	systematicError		( params->getDoubleValue ( "msms_parent_mass_systematic_error", false ) ),
	massInfo			( params ),

	consideredAA		( params ),

	compSearchParams	( params ),

	modificationParameters	( params, "msms_" ),
	linkInfo ( new LinkInfo ( params ) ),
	numSavedCrosslinkHits ( params->getIntValue ( "max_saved_tag_hits", 1000 ) ),

	biemannParams	( params ),
	expectationMethod ( params->getStringValue ( "expect_calc_method", "None" ) ),

	scoreHistogramOnly ( params->getBoolValue ( "score_histogram_only" ) )
{
	string svalue;
	regularExpression = ".";
	if ( params->getValue ( "regular_expression", svalue ) )	{
		regularExpression = ( svalue == "" ) ? "." : svalue;
	}
}
MSTagParameters::~MSTagParameters ()
{
	delete dataSetInfo;
}
char MSTagParameters::getAllowNonSpecificType () const
{
	string allowNonSpecific = getAllowNonSpecific ();
	if ( getNoEnzyme () ) return 'E';
	else if ( allowNonSpecific == "N termini-1=D" ) return 'D';
	else return allowNonSpecific [3];
}
unsigned int MSTagParameters::getLinkAminoAcid ( int num ) const
{
	static string aa = "nKCQcADEFGHILMNPRSTVWY";
	string linkAA = ( num == 1 ) ? linkInfo->getLinkAminoAcid1 () : linkInfo->getLinkAminoAcid2 ();
	StringVector sv = genGetSeparatedValues ( linkAA, "," );
	unsigned int mask = 0;
	for ( StringVectorSizeType i = 0 ; i < sv.size () ; i++ ) {
		int idx;
		if ( sv [i] == "Protein N-term" || sv [i] == "Peptide N-term" )		idx = aa.find ( 'n' );
		else if ( sv [i] == "Protein C-term" || sv [i] == "Peptide C-term" )idx = aa.find ( 'c' );
		else																idx = aa.find ( sv [i] );
		mask |= 1 << idx;
	}
	return mask;
}
string MSTagParameters::getBridgeFormula () const
{
	return linkInfo->getBridgeFormula ();
}
bool MSTagParameters::isCrosslinking () const
{
	return !linkInfo->getNoLink ();
}
void MSTagParameters::printHTML ( ostream& os ) const
{
	MSSearchParameters::printHTML ( os );
	if ( regularExpression != "" && regularExpression != "." )
		os << "Peptide contains the regular expression: <b>" << regularExpression << "</b><br />" << endl;
	biemannParams.printHTML ( os );
	modificationParameters.printHTML ( os );
	os << "Peptide Masses are: <b>" << massInfo.getMassType () << "</b><br />" << endl;
	os << "<p />" << endl;
}
int MSTagParameters::getNumSavedHits () const
{
	if ( isCrosslinking () )	return numSavedCrosslinkHits;
	else						return getMaxReportedHits ();
}
void MSTagParameters::setDataSetInfo ( const ParameterList* params, bool noFilesOK )
{
	delete dataSetInfo;
	dataSetInfo = new MSMSDataSetInfo ( params, noFilesOK );
}
