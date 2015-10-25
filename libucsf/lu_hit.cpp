/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_hit.cpp                                                    *
*                                                                             *
*  Created    : May 21st 2001                                                 *
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
#include <lu_pi.h>
#include <lu_mass_seq.h>
#include <lu_hit.h>
#include <lu_delim.h>
#include <lu_acc_link.h>
#include <lu_param_list.h>
#include <lu_fasta.h>
#include <lu_table.h>
using std::ostream;
using std::ostringstream;
using std::endl;
using std::string;

bool ProteinHit::uniprot = false;
bool ProteinHit::dnaDatabase = false;
std::map <const FastaServer*, int> ProteinHit::idxMap;

ProteinHit::ProteinHit ( FastaServer* fs, int ind, int drf, int orf ) :
	fs ( fs ),
	databaseEntry ( ind, drf, orf ),
	entrySet ( false )
{
}
void ProteinHit::printHTMLHeader ( ostream& os ) const
{
	databaseEntry.printHTMLHeader ( os, dnaDatabase );
	printHTMLHeader2 ( os );
}
void ProteinHit::printHTMLHeader2 ( ostream& os ) const
{
	tableHeader ( os, "Protein MW<br />(Da)/pI" );
	tableHeader ( os, "Accession<br />#" );
	tableHeader ( os, "Species" );
	tableHeader ( os, "Protein Name", "", "left" );
}
void ProteinHit::printHTMLHeader3 ( ostream& os ) const
{
	databaseEntry.printHTMLHeader ( os, dnaDatabase );
	tableHeader ( os, "Protein MW<br />(Da)/pI" );
	tableHeader ( os, "Accession<br />#" );
	fs->firstLine ( databaseEntry.getIndex () );
	string s = fs->getLineSequenceVersion ();
	if ( !s.empty () ) uniprot = true;
	if ( uniprot ) {
		tableHeader ( os, "UniProt ID" );
		tableHeader ( os, "Organism" );
		tableHeader ( os, "Gene" );
		tableHeader ( os, "Existence" );
		tableHeader ( os, "Version" );
	}
	tableHeader ( os, "Species" );
	tableHeader ( os, "Protein Name", "", "left" );
}
void ProteinHit::printDelimitedHeader ( ostream& os ) const
{
	databaseEntry.printDelimitedHeader ( os, dnaDatabase );
	printDelimitedHeader2 ( os );
}
void ProteinHit::printDelimitedHeader2 ( ostream& os ) const
{
	delimitedHeader ( os, "Protein MW (Da)" );
	delimitedHeader ( os, "Protein pI" );
	delimitedHeader ( os, "Accession #" );
	delimitedHeader ( os, "Species" );
	delimitedHeader ( os, "Protein Name" );
}
void ProteinHit::printHit ( ostream& os ) const
{
	int idx = getDBIndex ();
	for ( fs->firstLine ( databaseEntry.getIndex () ) ; fs->isDoneLine () ; fs->nextLine () ) {
		os << "<b>Acc. #: </b>";
		printAccessionNumber ( os, fs->getLineAccessionNumber (), false, idx );
		os << " ";

		string uniprotID = fs->getLineUniprotID ();
		if ( !uniprotID.empty () )	os << "<b>Uniprot ID: </b>" << uniprotID << " ";

		os << "<b>Species: </b>";
		os << fs->getLineSpecies ();
		os << " ";

		os << "<b>Name: </b>";
		os << fs->getLineName ();
		os << "<br />" << endl;

		string str;

		string organismName = fs->getLineOrganismName ();
		if ( !organismName.empty () )			str += "<b>Organism: </b>"			+ organismName + " ";

		string geneName = fs->getLineGeneName ();
		if ( !geneName.empty () )				str += "<b>Gene: </b>"				+ geneName + " ";

		string proteinExistence = fs->getLineProteinExistence ();
		if ( !proteinExistence.empty () )		str += "<b>Existence: </b>"			+ proteinExistence + " ";

		string sequenceVersion = fs->getLineSequenceVersion ();
		if ( !sequenceVersion.empty () )		str += "<b>Version: </b>"			+ sequenceVersion + " ";

		if ( !str.empty () ) os << str.substr ( 0, str.length () - 1 ) << "<br />" << endl;
	}
	databaseEntry.printHTMLHit ( os, idx, fs->getDNADatabase () );
	os << " ";
	os << "<b>MW: </b>";
	printProteinMW ( os );
	os << " Da";
	os << " ";

	os << "<b>pI: </b>";
	printProteinPI ( os );

	os << "<br />" << endl;
}
void ProteinHit::printDelimitedAccNum ( ostream& os ) const
{
	ostringstream ost;
	if ( idxMap.size () > 1 ) ost << getDBIndex ()+1 << '$';
	ost << getAccessionNumber ();
	if ( fs->getDNADatabase () ) {
		ost << '$' << databaseEntry.getDNAReadingFrame () << '$' << databaseEntry.getOpenReadingFrame () + 1;
	}
	delimitedCell ( os, ost.str () );
}
void ProteinHit::printXMLHit ( ostream& os ) const
{
	ParameterList::printDoubleXMLFixed ( os, "protein_mw", getProteinMW (), ProteinMW::mwPrecision );
	if ( getProteinPI () == PI_NOT_CALCULATED )
		ParameterList::printXML ( os, "protein_pi", "----" );
	else
		ParameterList::printDoubleXMLFixed ( os, "protein_pi", getProteinPI (), ProteinPI::piPrecision );
	if ( idxMap.size () > 1 ) ParameterList::printXML ( os, "database", fs->getFileName () );
	databaseEntry.printXML ( os, fs->getDNADatabase () );
	for ( fs->firstLine ( databaseEntry.getIndex () ) ; fs->isDoneLine () ; fs->nextLine () ) {
		ParameterList::printXML ( os, "accession_number", fs->getLineAccessionNumber () );
		ParameterList::printXML ( os, "species", fs->getLineSpecies () );
		ParameterList::printXML ( os, "name", fs->getLineName () );
	}
}
void ProteinHit::printXML ( ostream& os ) const
{
	os << "<hit>" << endl;
		printXMLHit ( os );
	os << "</hit>" << endl;
}
void ProteinHit::printHTMLHit ( ostream& os ) const
{
	tableHeaderStart ( os, "", "", true );
		printProteinMW ( os );
		os << "/";
		printProteinPI ( os );
		os << endl;
	tableHeaderEnd ( os );
	tableHeaderStart ( os, "", "", true );
		int i = 0;
		string aNum;
		for ( fs->firstLine ( databaseEntry.getIndex () ) ; fs->isDoneLine () ; fs->nextLine (), i++ ) {
			if ( i == 0 ) aNum = fs->getLineAccessionNumber ();
			if ( i > 1 ) break;
		}
		printAccessionNumber ( os, aNum, i > 1, getDBIndex () );
	tableHeaderEnd ( os );
	tableCell ( os, getSpecies () );
	tableCell ( os, getName () );
}
void ProteinHit::printHTMLHit2 ( ostream& os ) const
{
	databaseEntry.printHTML ( os, getDBIndex (), dnaDatabase, fs->getDNADatabase () );
	printHTMLHit ( os );
}
void ProteinHit::printHTMLMultiHit ( ostream& os ) const
{
	bool first = true;
	for ( fs->firstLine ( databaseEntry.getIndex () ) ; fs->isDoneLine () ; fs->nextLine () ) {
		tableRowStart ( os );
			if ( first ) {
				databaseEntry.printHTML ( os, getDBIndex (), dnaDatabase, fs->getDNADatabase () );
				tableHeaderStart ( os, "", "", true );
					printProteinMW ( os );
					os << "/";
					printProteinPI ( os );
					os << endl;
				tableHeaderEnd ( os );
			}
			else {
				databaseEntry.printGapHTML ( os, dnaDatabase );
				tableEmptyCell ( os );
			}
			tableHeaderStart ( os, "", "", true );
				printAccessionNumber ( os, fs->getLineAccessionNumber (), false, getDBIndex () );
			tableHeaderEnd ( os );
			if ( uniprot ) {
				tableCell ( os, fs->getLineUniprotID () );
				tableCell ( os, fs->getLineOrganismName () );
				tableCell ( os, fs->getLineGeneName () );
				tableCell ( os, fs->getLineProteinExistence () );
				tableCell ( os, fs->getLineSequenceVersion () );
			}
			tableCell ( os, fs->getLineSpecies () );
			tableCell ( os, fs->getLineName () );
			first = false;
		tableRowEnd ( os );
	}
}
void ProteinHit::printDelimitedHit ( ostream& os ) const
{
	ostringstream ost;
	printProteinMW ( ost );
	delimitedCell ( os, ost.str () );
	ostringstream ost2;
	printProteinPI ( ost2 );
	delimitedCell ( os, ost2.str () );
	int i = 0;
	string aNum;
	for ( fs->firstLine ( databaseEntry.getIndex () ) ; fs->isDoneLine () ; fs->nextLine (), i++ ) {
		if ( i == 0 ) aNum = fs->getLineAccessionNumber ();
		if ( i > 1 ) break;
	}
	printDelimitedAccessionNumber ( os, aNum, i > 1 );
	delimitedCell ( os, getSpecies () );
	delimitedCell ( os, getName () );
}
void ProteinHit::printDelimitedHit2 ( ostream& os ) const
{
	databaseEntry.printDelimited ( os, getDBIndex (), dnaDatabase, fs->getDNADatabase () );
	printDelimitedHit ( os );
}
void ProteinHit::printAccessionNumber ( ostream& os, const string& an, bool multi, int num ) const
{
	AccessionNumberLinkInfo::write ( os, an, num, true );
	if ( multi ) os << "M";
}
void ProteinHit::printDelimitedAccessionNumber ( ostream& os, const string& an, bool multi ) const
{
	ostringstream ost;
	ost << an;
	if ( multi ) ost << "M";
	delimitedCell ( os, ost.str () );
}
void ProteinHit::printProteinMW ( ostream& os ) const
{
	genPrint ( os, getProteinMW (), ProteinMW::mwPrecision, 8 );
}
void ProteinHit::printProteinPI ( ostream& os ) const
{
	if ( getProteinPI () == PI_NOT_CALCULATED ) os << "----";
	else {
		genPrint ( os, getProteinPI (), ProteinPI::piPrecision );
	}
}
double ProteinHit::getProteinMW () const
{
	if ( entrySet == false ) getEntry ();
	return proteinMW;
}
double ProteinHit::getProteinPI () const
{
	if ( entrySet == false ) getEntry ();
	return proteinPI;
}
string ProteinHit::getSpecies () const
{
	if ( entrySet == false ) getEntry ();
	return species;
}
string ProteinHit::getName () const
{
	if ( entrySet == false ) getEntry ();
	return name;
}
string ProteinHit::getAccessionNumber () const
{
	if ( entrySet == false ) getEntry ();
	return accessionNumber;
}
void ProteinHit::getEntry () const
{
	int nTermLimit = fs->getMaxNTermAA ();
	fs->setMaxNTermAA ( 0 );
	char* frame = fs->getProtein ( databaseEntry );
	ProteinMW pmw ( frame );
	proteinMW = pmw.getMass ();
	ProteinPI ppi ( frame );
	proteinPI = ppi.getProteinPI ();
	int index = databaseEntry.getIndex ();
	accessionNumber = fs->getAccessionNumber ( index );
	species = fs->getSpecies ( index );
	name = fs->getName ( index );
	fs->setMaxNTermAA ( nTermLimit );
	entrySet = true;
}
void ProteinHit::addFS ( const FastaServer* fs, int num )
{
	idxMap [fs] = num;
	if ( fs->getDNADatabase () ) dnaDatabase = true;
}
void ProteinHit::reset ()
{
	idxMap.clear ();
	dnaDatabase = false;
}
bool ProteinHit::isDecoy () const
{
	string a = getAccessionNumber ();
	if ( a [0] == '-' ) return true;
	if ( fs->getDecoy () ) return true;
	return false;
}
