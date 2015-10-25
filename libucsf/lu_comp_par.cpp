/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_comp_par.cpp                                               *
*                                                                             *
*  Created    : June 9th 2001                                                 *
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
#include <lu_iso_par.h>
#include <lu_comp_par.h>
#include <lu_param_list.h>
using std::ostream;

MSCompParameters::MSCompParameters ( const ParameterList* params ) :
	MSProgramParameters	( params ),

	aaInitInfo			( params ),
	parentMass			( params->getDoubleValue	( "msms_parent_mass", 1000.0 ) ),
	parentCharge		( params->getIntValue	( "msms_parent_charge", 1 ) ),
	combinationType		( params->getStringValue	( "combination_type", "Amino Acid" ) ),
	consideredAA		( params ),
	maxReportedHits		( params->getIntValue	( "msms_max_reported_hits", 50 ) ),
	parentMassTolerance( "msms_parent_mass", params ),
	massInfo			( params ),
	compSearchParams	( params ),
	requestedIonTypes	( params->getStringVectorValue( "it" ) )
{
}
void MSCompParameters::printHTML ( ostream& os ) const
{
	aaInitInfo.printHTML ( os );
	massInfo.printHTML ( os );
	consideredAA.printHTML ( os );
	compSearchParams.printHTML ( os );
}
