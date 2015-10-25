/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_msp.cpp                                                    *
*                                                                             *
*  Created    : March 6th 2014                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2014-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <string>
#include <lg_io.h>
#include <lg_string.h>
#include <lgen_define.h>
#include <lgen_error.h>
#include <lgen_file.h>
#include <lgen_math.h>
#include <lu_aa_info.h>
#include <lu_mass_elem.h>
#include <lu_msp.h>

using std::istringstream;
using std::istream;
using std::ostream;
using std::string;
using std::getline;
using std::endl;

LibraryFileRead::LibraryFileRead ( const string& fname ) :
	fname ( fname )
{
}
void LibraryFileRead::read ( const string& peakListPath, StringVectorVector& header, StringVectorVector& rows )
{
	genCreateDirectory ( peakListPath );
	StringVector hCols;
	hCols.push_back ( "Prev AA" );
	hCols.push_back ( "DB Peptide" );
	hCols.push_back ( "Next AA" );
	hCols.push_back ( "Modifications" );
	hCols.push_back ( "Precursor MZ" );
	hCols.push_back ( "Precursor Charge" );
	hCols.push_back ( "Spectrum" );
	hCols.push_back ( "Num Peaks" );
	hCols.push_back ( "Peptide Type" );
	hCols.push_back ( "Inst" );
	hCols.push_back ( "Protein" );
	header.push_back ( hCols );
	string line;
	string fileName = "f";
	peakListFractionNames.push_back ( fileName );
	peakListCentroidFileNames.push_back ( fileName + ".mgf" );
	peakListCentroidFilePaths.push_back ( peakListPath + SLASH + peakListCentroidFileNames.back () );
	GenOFStream os ( peakListCentroidFilePaths.back () );
	GenIFStream ist ( fname );
	processFile ( ist, os, rows );
}
string LibraryFileRead::getMods ( const string& value ) const
{
//	Mods=0
//	Mods=1/-1,Y,USM_n_230.170762
//	Mods=1/0,C,ICAT-C:13C(9)
//	Mods=1/4,C,ICAT-C:13C(9)
//	Mods=2/5,C,ICAT-C/6,M,Oxidation
	string mods;
	if ( value [0] != '0' ) {
		string::size_type s1 = value.find ( '/' );
		int num = atoi ( value.substr ( 0, s1 ).c_str () );				// Number of modifications
		for ( int i = 0 ; i < num ; i++ ) {
			s1++;
			string::size_type e1 = value.find ( '/', s1 );
			string m = value.substr ( s1, e1-s1 );
			string::size_type e2 = m.find ( ',' );
			int idx = atoi ( m.substr ( 0, e2 ).c_str () ) + 1;
			string::size_type e3 = m.find ( ',', e2+1 );
			string m1 = m.substr ( e3+1 );
			if ( isPrefix ( m1, "Unknown_" ) ) {
				m1 = m1.substr ( 8 );
			}
			else if ( isPrefix ( m1, "Nterm_" ) ) {
				m1 = m1.substr ( 6 );
			}
			/*		This code converts any USM codes into a mass mod. However the mods don't completely agree in mass as the mod
					includes the mass of the amino acid.
			else if ( isPrefix ( m1, "USM_" ) ) {
				char aa = m1 [4];
				m1 = m1.substr ( 6 );
				double mass = atof ( m1.c_str () );
				if ( aa == 'n' ) {
					static double h1 = formula_to_monoisotopic_mass ( "H" );
					mass -= h1;
				}
				else {
					mass -= AAInfo::getInfo ().getMonoisotopicMass ( aa );
				}
				m1 = gen_ftoa ( mass, "%.6f" );
			}
			*/
			mods += m1 + "@";
			if ( idx == 0 ) mods += "N-term";
			else			mods += gen_itoa ( idx );
			if ( i != num - 1 ) mods += ";";
			s1 = e1;
		}
	}
	return mods;
}
void LibraryFileRead::addRow ( StringVectorVector& rows ) const
{
	StringVector cols;
	cols.push_back ( prevAA );
	cols.push_back ( dbPeptide );
	cols.push_back ( nextAA );
	cols.push_back ( mods );
	cols.push_back ( sPrecursorMZ );
	cols.push_back ( sPrecursorZ );
	cols.push_back ( spectrum );
	cols.push_back ( sNumPeaks );
	cols.push_back ( peptideType );
	cols.push_back ( inst );
	cols.push_back ( protein );
	rows.push_back ( cols );
}
void LibraryFileRead::writeSpectrum ( istream& istr, ostream& os, int num, const string& sPrecursorMZ, const string& sPrecursorZ, int numPeaks ) const
{
	os << "BEGIN IONS" << endl;
	os << "TITLE=Scan ";
	os << num;
	os << " ";
	os << "(rt=";
	os << num;
	os << ")";
	os << " ";
	os << "[Prospector Created]";
	os << endl;

	os << "PEPMASS=" << sPrecursorMZ << endl;
	os << "CHARGE=" << sPrecursorZ << "+" << endl;

	string line;
	for ( DoubleVectorSizeType i = 0 ; i < numPeaks ; i++ ) {
		getline ( istr, line );
		istringstream iii ( line );
		string mass;
		string inten;
		iii >> mass;
		iii >> inten;

		os << mass << " " << inten << endl;
	}
	os << "END IONS" << endl;
	getline ( istr, line );
}
MSPRead::MSPRead ( const string& name ) :
	LibraryFileRead ( name )
{
}
void MSPRead::processFile ( istream& ist, ostream& os, StringVectorVector& rows )
{
	string line;
	specNum = 1;
	while ( getline ( ist, line ) ) {
		if ( isPrefix ( line, "Name" ) ) {
			string::size_type start = line.find ( ' ' ) + 1;
			string::size_type end = line.find ( '/' );
			dbPeptide = gen_strstriptags2 ( line.substr ( start, end-start ), '(', ')' );
			end++;
			string::size_type end2 = line.find ( 13, end );
			sPrecursorZ = line.substr ( end, end2-end );
		}
		else if ( isPrefix ( line, "Comment" ) ) {
			string::size_type start = 0;
			string::size_type end;
			string s1 = genNextString ( line, " ", start, end );	// Skip comment
			for ( ; ; ) {
				string name = genNextString ( line, "=", start, end );
				if ( name == "Pep" ) {
					peptideType = genNextString ( line, " ", start, end );
				}
				else if ( name == "Fullname" ) {
					string value = genNextString ( line, " ", start, end );
					prevAA = string ( 1, value [0] );
					string::size_type s = value.find ( '/' ) - 1;
					nextAA = string ( 1, value [s] );
				}
				else if ( name == "Mods" ) {
					mods = getMods ( genNextString ( line, " ", start, end ) );
				}
				else if ( name == "Parent" ) {
					sPrecursorMZ = genNextString ( line, " ", start, end );
				}
				else if ( name == "Inst" ) {
					inst = genNextString ( line, " ", start, end );
				}
				else if ( name == "Protein" ) {
					start++;
					protein = genNextString ( line, "\"", start, end );
					end++;
					break;
				}
				else {
					string value = genNextString ( line, " ", start, end );
				}
			}
		}
		else if ( isPrefix ( line, "Num peaks" ) ) {
			string::size_type start = line.find ( ':' ) + 2;
			string::size_type end = line.find ( 13 );		// Checks for dos 13
			sNumPeaks = line.substr ( start, end-start );
			writeSpectrum ( ist, os, specNum, sPrecursorMZ, sPrecursorZ, atoi ( sNumPeaks.c_str () ) );
			spectrum = gen_itoa ( specNum );
			specNum++;
			addRow ( rows );
		}
		else if ( isPrefix ( line, "NumPeaks" ) ) {			// This is probably a SpectraST file
			ErrorHandler::genError ()->error ( "File format error.\n" );
		}
	}
}
SPTXTRead::SPTXTRead ( const string& name ) :
	LibraryFileRead ( name )
{
}
void SPTXTRead::processFile ( istream& ist, ostream& os, StringVectorVector& rows )
{
	string line;
	specNum = 1;
	while ( getline ( ist, line ) ) {
		if ( isPrefix ( line, "Name" ) ) {
			string::size_type start = line.find ( ' ' ) + 1;
			string::size_type end = line.find ( '/' );
			dbPeptide = gen_strstriptags2 ( line.substr ( start, end-start ), '[', ']' );
			if ( dbPeptide [0] == 'n' ) dbPeptide = dbPeptide.substr ( 1 );					// N-terminal modification
			end++;
			string::size_type end2 = line.find ( 13, end );
			sPrecursorZ = line.substr ( end, end2-end );
		}
		else if ( isPrefix ( line, "PrecursorMZ" ) ) {
			string::size_type start = line.find ( ':' ) + 2;
			string::size_type end = line.find ( 13 );		// Checks for dos 13
			sPrecursorMZ = line.substr ( start, end-start );
		}
		else if ( isPrefix ( line, "FullName" ) ) {
			string::size_type start = line.find ( ' ' ) + 1;
			string::size_type end = line.find ( '/' );
			string value = line.substr ( start, end-start );
			prevAA = string ( 1, value [0] );
			nextAA = string ( 1, value [value.length ()-1] );
		}
		else if ( isPrefix ( line, "Comment" ) ) {
			string::size_type start = 0;
			string::size_type end;
			string s1 = genNextString ( line, " ", start, end );	// Skip comment
			for ( ; ; ) {
				string name = genNextString ( line, "=", start, end );
				if ( name == "Pep" ) {
					peptideType = genNextString ( line, " ", start, end );
				}
				else if ( name == "Mods" ) {
					mods = getMods ( genNextString ( line, " ", start, end ) );
				}
				else if ( name == "Inst" ) {
					inst = genNextString ( line, " ", start, end );
				}
				else if ( name == "Protein" ) {
					protein = genNextString ( line, " ", start, end );
					break;
				}
				else {
					string value = genNextString ( line, " ", start, end );
				}
			}
		}
		else if ( isPrefix ( line, "NumPeaks" ) ) {
			string::size_type start = line.find ( ':' ) + 2;
			string::size_type end = line.find ( 13 );		// Checks for dos 13
			sNumPeaks = line.substr ( start, end-start );
			writeSpectrum ( ist, os, specNum, sPrecursorMZ, sPrecursorZ, atoi ( sNumPeaks.c_str () ) );
			spectrum = gen_itoa ( specNum );
			specNum++;
			addRow ( rows );
		}
	}
}
