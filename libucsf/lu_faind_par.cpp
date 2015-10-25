/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_faind_par.cpp                                              *
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
*  Copyright (2001-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lg_string.h>
#include <lu_faind_par.h>
#include <lu_param_list.h>
using std::vector;

FAIndexParameters::FAIndexParameters ( const ParameterList* params ) :
	MSProgramParameters	( params ),
	aaInitInfo		( params ),
	database		( params->getStringValue	( "database", "" ) ),
	createSubDatabase	( params->getBoolValue	( "create_sub_database", false ) ),
	createUserDatabase	( params->getBoolValue	( "create_user_database", false ) ),
	createDatabaseIndicies	( params->getBoolValue	( "create_database_indicies", false ) )
{
}
FAIndexNormalParameters::FAIndexNormalParameters ( const ParameterList* params ) :
	FAIndexParameters	( params ),
	dnaToProtein		( params->getBoolValue	( "dna_to_protein", false ) ),
	randomDatabase		( params->getBoolValue	( "random_database", false ) ),
	reverseDatabase		( params->getBoolValue	( "reverse_database", false ) ),
	concatDatabase		( params->getBoolValue	( "concat_database", false ) ),
	deleteDNADatabase	( params->getBoolValue	( "delete_dna_database", false ) )
{
}
FAIndexSubsetDatabaseParameters::FAIndexSubsetDatabaseParameters ( const ParameterList* params ) :
	FAIndexParameters	( params ),
	subDatabaseID		( params->getStringValue( "sub_database_id", ".sub" ) ),
	preSearchInfo		( params )
{
}
IntVector FAIndexSubsetDatabaseParameters::getIndicies ( FastaServer* fsPtr )
{
	vector <FastaServer*> fsv;
	fsv.push_back ( fsPtr );
	preSearchInfo.doSearch ( fsv );
	return preSearchInfo.getIndicies ( 0 );
}
FAIndexDatabaseSummaryParameters::FAIndexDatabaseSummaryParameters ( const ParameterList* params ) :
	FAIndexParameters	( params ),
	dnaReadingFrame		( params->getIntValue	( "dna_reading_frame", 1 ) ),
	startIndexNumber	( params->getIntValue	( "start_index_number", 1 ) ),
	endIndexNumber		( params->getIntValue	( "end_index_number", 1 ) ),
	allIndicies			( params->getBoolValue	( "all_indicies", false ) ),
	hideProteinSequence	( params->getBoolValue	( "hide_protein_sequence", false ) )
{
}
FAIndexCreateOrAppendParameters::FAIndexCreateOrAppendParameters ( const ParameterList* params ) :
	FAIndexParameters	( params ),
	nameField			( params->getStringValue( "name_field", "" ) ),
	accessionNum		( params->getStringValue( "accession_num", "L39370" ) ),
	species				( params->getStringValue( "species", "All" ) )
{
	const char* value;
	if ( params->getValue ( "user_protein_sequence", value ) )	userProteinSequence = gen_strstrip ( value );
}
