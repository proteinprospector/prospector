/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_prot_srch.cpp                                              *
*                                                                             *
*  Created    : Septmber 5th 2001                                             *
*                                                                             *
*  Purpose    : Database search and hit base classes.                         *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lgen_error.h>
#include <lu_prot_srch.h>
#include <lu_srch_par.h>
using std::ostream;
using std::string;

ProteinSearch::ProteinSearch ( const MSSearchParameters& params ) :
	DatabaseSearch ( params )
{
	ProteinHits* proteinHits = new ProteinHits ();
	for ( int i = 0 ; i < fs.size () ; i++ ) {
		int numIndicies = params.getNumIndicies ( i );
		if ( numIndicies > 2000000 ) {
			string message = "Too many protein hits (" + gen_itoa (numIndicies) + ")";
			if ( fs.size () > 1 ) message += " in database " + fs [i]->getFileName ();
			message += ".\n";
			ErrorHandler::genError ()->error ( message );
		}
	}
	for ( int j = 0 ; j < fs.size () ; j++ ) {
		FastaServer* fsPtr = fs [j];
		ProteinHit::addFS ( fsPtr, j );
		IntVector indicies ( params.getIndicies ( j ) );
		for ( IntVectorSizeType k = 0 ; k < indicies.size () ; k++ ) {
			proteinHits->addHit ( ProteinHit ( fsPtr, indicies [k], 1, 0 ) );
		}
	}
	numHits = proteinHits->size ();
	databaseHits = proteinHits;
}
void ProteinSearch::printParamsBodyHTML ( ostream& os ) const
{
	params.printHTML ( os );
}
