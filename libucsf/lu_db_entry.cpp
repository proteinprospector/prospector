/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_db_entry.cpp                                               *
*                                                                             *
*  Created    :                                                               *
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
#include <lu_db_entry.h>
#include <lu_dig_par.h>
#include <lu_param_list.h>
#include <lu_table.h>
#include <lu_delim.h>
using std::ostream;
using std::endl;

DatabaseEntry::DatabaseEntry ( int index, int dnaReadingFrame, int openReadingFrame ) :
	index ( index ),
	dnaReadingFrame ( dnaReadingFrame ),
	openReadingFrame ( openReadingFrame )
{
}
bool DatabaseEntry::operator== ( const DatabaseEntry& rhs ) const
{
	if ( index == rhs.index ) {
		if ( dnaReadingFrame == rhs.dnaReadingFrame ) {
			if ( openReadingFrame == rhs.openReadingFrame ) {
				return true;
			}
			else return false;
		}
		else return false;
	}
	else return false;
}
void DatabaseEntry::printHTMLHeader ( ostream& os, bool dnaDatabase ) const
{
	tableHeader ( os, "MS-Digest<br />Index #" );
	if ( dnaDatabase ) {
		tableHeader ( os, "DNA Reading<br />Frame" );
		tableHeader ( os, "Open Reading<br />Frame" );
	}
}
void DatabaseEntry::printDelimitedHeader ( ostream& os, bool dnaDatabase ) const
{
	delimitedHeader ( os, "Index" );
	if ( dnaDatabase ) {
		delimitedHeader ( os, "DNA Reading Frame" );
		delimitedHeader ( os, "Open Reading Frame" );
	}
}
void DatabaseEntry::printHTMLHit ( ostream& os, int num, bool dnaEntry ) const
{
	os << "<b>Index: </b>";
	MSDigestLink::write ( os, index, dnaReadingFrame, openReadingFrame, num );
	if ( dnaEntry ) {
		os << " <b>DNA Reading Frame: </b>" << dnaReadingFrame;
		os << " <b>Open Reading Frame: </b>" << openReadingFrame + 1;
	}
	os << " ";
}
void DatabaseEntry::printHTML ( ostream& os, int num, bool dnaDatabase, bool dnaEntry ) const
{
	tableHeaderStart ( os );
		MSDigestLink::write ( os, index, dnaReadingFrame, openReadingFrame, num );
		os << endl;
	tableHeaderEnd ( os );
	if ( dnaDatabase ) {
		if ( dnaEntry ) {
			tableCell ( os, dnaReadingFrame, true );
			tableCell ( os, openReadingFrame + 1, true );
		}
		else tableEmptyNCells ( os, 2 );
	}
}
void DatabaseEntry::printGapHTML ( ostream& os, bool dnaDatabase ) const
{
	tableEmptyCell ( os );
	if ( dnaDatabase ) {
		tableEmptyNCells ( os, 2 );
	}
}
void DatabaseEntry::printDelimited ( ostream& os, int num, bool dnaDatabase, bool dnaEntry ) const
{
	delimitedCell ( os, index );
	if ( dnaDatabase ) {
		if ( dnaEntry ) {
			delimitedCell ( os, dnaReadingFrame );
			delimitedCell ( os, openReadingFrame + 1 );
		}
		else delimitedEmptyNCells ( os, 2 );
	}
}
void DatabaseEntry::printXML ( ostream& os, bool dnaEntry ) const
{
	ParameterList::printXML ( os, "index_number", index );
	if ( dnaEntry ) {
		ParameterList::printXML ( os, "dna_reading_frame", dnaReadingFrame );
		ParameterList::printXML ( os, "open_reading_frame", openReadingFrame + 1 );
	}
}
