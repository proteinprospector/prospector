/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_species.cpp                                                *
*                                                                             *
*  Created    : July 20th 1996                                                *
*                                                                             *
*  Purpose    : Calculates the masses of fragment ions for a peptide          *
*               sequence.                                                     *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2008-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <stdexcept>
#include <algorithm>
#include <lg_io.h>
#include <lg_string.h>
#include <lu_getfil.h>
#include <lu_species.h>

using std::string;
using std::istringstream;
using std::getline;
using std::reverse;
using std::runtime_error;

class NCBIMergedNodes {
	MapIntToInt mergedNodes;
public:
	NCBIMergedNodes ();
	int getNewNode ( int node ) const;
};

/*
Types of species name:

acronym
anamorph
authority
blast name
common name
equivalent name
genbank acronym
genbank anamorph
genbank common name
genbank synonym
in-part
includes
misnomer
misspelling
scientific name
synonym
teleomorph
unpublished name

Kingdom codes:

'A' for archaea (=archaebacteria),
'B' for bacteria (=prokaryota or eubacteria),
'E' for eukaryota (=eukarya),
'V' for viruses and phages (=viridae).
'X' Root taxonomy node
*/

/*	In order of occurrence in nodes.dmp - July 14th 2008

species	330570
genus	45860
no rank	41529
subspecies	8746
family	6392
varietas	2645
subfamily	1965
tribe	1243
order	1088
subgenus	791
superfamily	703
suborder	343
class	238
species group	210
forma	207
subtribe	196
subclass	120
species subgroup	89
phylum	88
infraorder	80
superorder	59
subphylum	24
infraphylum	11
parvorder	6
superclass	5
superphylum	4
superkingdom	3
kingdom	3
subkingdom	1
infrakingdom	0
microphylum	0
infraclass	0
parvclass	0
magnorder	0
superspecies	0
infraspecies	0
subvarietas	0
subforma	0
*/
namespace {
const char* NCBIRanks [] = {	"No Rank",
								"Superkingdom",
								"Kingdom",
								"Subkingdom",
								"Infrakingdom",
								"Superphylum",
								"Phylum",
								"Subphylum",
								"Infraphylum",
								"Microphylum",
								"Superclass",
								"Class",
								"Subclass",
								"Infraclass",
								"Parvclass",
								"Magnorder",
								"Superorder",
								"Order",
								"Suborder",
								"Infraorder",
								"Parvorder",
								"Superfamily",
								"Family",
								"Subfamily",
								"Tribe",
								"Subtribe",
								"Genus",
								"Subgenus",
								"Species Group",
								"Species Subgroup",
								"Superspecies",
								"Species",
								"Subspecies",
								"Infraspecies",
								"Varietas",
								"Subvarietas",
								"Forma",
								"Subforma",
								"Unknown", 0 };
}

enum NCBIRank {
	NoRank = 0,
	Superkingdom = 1,
	Kingdom = 2,
	Subkingdom = 3,
	Infrakingdom = 4,
	Superphylum = 5,
	Phylum = 6,
	Subphylum = 7,
	Infraphylum = 8,
	Microphylum = 9,
	Superclass = 10,
	Class = 11,
	Subclass = 12,
	Infraclass = 13,
	Parvclass = 14,
	Magnorder = 15,
	Superorder = 16,
	Order = 17,
	Suborder = 18,
	Infraorder = 19,
	Parvorder = 20,
	Superfamily = 21,
	Family = 22,
	Subfamily = 23,
	Tribe = 24,
	Subtribe = 25,
	Genus = 26,
	Subgenus = 27,
	SpeciesGroup = 28,
	SpeciesSubgroup = 29,
	Superspecies = 30,
	Species = 31,
	Subspecies = 32,
	Infraspecies = 33,
	Varietas = 34,
	Subvarietas = 35,
	Forma = 36,
	Subforma = 37,
	Unknown = 38
};

UniProtSpeciesMap::UniProtSpeciesMap ()
{
	GenIFStream ifs ( MsparamsDir::instance ().getParamPath ( "taxonomy/speclist.txt" ) );
	string line;
	int i = 0;
	bool flag = false;
	while ( getline ( ifs, line ) ) {
		if ( gen_strcharcount ( line, '_' ) > 20 ) {
			flag = true;
			break;
		}
	}
	if ( flag ) {	// Format OK
		while ( getline ( ifs, line ) ) {
			if ( line.length () == 0 ) break;
			if ( !isspace ( line [0] ) ) {
				istringstream istr ( line );
				string code;
				char kingdom;
				int node;
				istr >> code;
				istr >> kingdom;
				istr >> node;
				specMap [code] = node;
			}
		}
		int count = 0;
		while ( getline ( ifs, line ) ) {	// Skip past comment lines
			if ( gen_strcharcount ( line, '=' ) > 20 ) count++;
			if ( count == 2 ) break;
		}
		bool flag2 = false;
		while ( getline ( ifs, line ) ) {
			if ( line.length () == 0 ) {
				if ( flag2 ) break;	// Only quit on blank lines after some lines processed
				else continue;
			}
			if ( !isspace ( line [0] ) ) {
				istringstream istr ( line );
				string code;
				char kingdom;
				int node;
				istr >> code;
				istr >> kingdom;
				istr >> node;
				specMap [code] = node;
				flag2 = true;
			}
		}
	}
}
NCBISpeciesMap::NCBISpeciesMap ( bool nodeNames )
{
	char* info = getFileAsCharPtr ( MsparamsDir::instance ().getParamPath ( "taxonomy/names.dmp" ) );
	bool first = true;
	for ( ; ; ) {
		char* ptr = first ? strtok ( info, "\t" ) : strtok ( NULL, "\t" );
		if ( ptr == 0 ) break;
		first = false;
		int node = atoi ( ptr );
		strtok ( NULL, "\t" );									// Skip next tab
		char* spec = strtok ( NULL, "\t" );
		for ( char* p1 = spec ; *p1 ; p1++ ) {
			*p1 = toupper ( *p1 );
		}
		if ( specMap.find ( spec ) == specMap.end () ) {
			specMap [spec] = node;
			if ( nodeNames ) {
				char* c1 = strtok ( NULL, "\t" );
				char* c2 = strtok ( NULL, "\t" );
				if ( c2-c1 == 2 ) strtok ( NULL, "\t" );
				string type = strtok ( NULL, "\t" );
				if ( type == "scientific name" ) nodeNameMap [node] = spec;
			}
		}
		strtok ( NULL, "\n" );									// Skip rest of line
	}
	delete [] info;

}
string NCBISpeciesMap::getScientificName ( int node ) const
{
	MapIntToStringConstIterator cur = nodeNameMap.find ( node );
	if ( cur != nodeNameMap.end () ) {
		return (*cur).second;
	}
	else {
		return gen_itoa ( node );
	}
}
NCBINodes::NCBINodes ( bool parentFlag, bool childFlag, bool rankFlag ) :
	mNodes ( 0 )
{
	char* info = getFileAsCharPtr ( MsparamsDir::instance ().getParamPath ( "taxonomy/nodes.dmp" ) );
	bool first = true;
	for ( ; ; ) {
		char* ptr = first ? strtok ( info, "\n" ) : strtok ( NULL, "\n" );
		if ( ptr == 0 ) break;
		first = false;
		int node = atoi ( ptr );
		while ( *ptr && *ptr != '\t' ) ptr++;					// Skip node
		while ( *ptr && ( *ptr == '\t' || *ptr == '|' ) ) ptr++;// Get to parent node
		int parentNode = atoi ( ptr );
		if ( rankFlag ) {
			while ( *ptr && *ptr != '\t' ) ptr++;					// Skip parent node
			while ( *ptr && ( *ptr == '\t' || *ptr == '|' ) ) ptr++;// Get to rank
			char* ptr2 = ptr;
			while ( *ptr2 && *ptr2 != '\t' ) ptr2++;				// Find end of rank

			parseRank ( string ( ptr, ptr2 - ptr ), node );
		}
		if ( childFlag ) childNodes [parentNode].push_back ( node );
		if ( parentFlag ) parent [node] = parentNode;
	}
	delete [] info;
}
NCBINodes::~NCBINodes ()
{
	delete mNodes;
}
void NCBINodes::parseRank ( const string& r, int node )
{
	if ( r == "species" )				nodeRank [node] = Species;		// By order of occurrence
	else if ( r == "genus" )			nodeRank [node] = Genus;
	else if ( r == "no rank" )			nodeRank [node] = NoRank;
	else if ( r == "subspecies" )		nodeRank [node] = Subspecies;
	else if ( r == "family" )			nodeRank [node] = Family;
	else if ( r == "varietas" )			nodeRank [node] = Varietas;
	else if ( r == "subfamily" )		nodeRank [node] = Subfamily;
	else if ( r == "tribe" )			nodeRank [node] = Tribe;
	else if ( r == "order" )			nodeRank [node] = Order;
	else if ( r == "subgenus" )			nodeRank [node] = Subgenus;
	else if ( r == "superfamily" )		nodeRank [node] = Superfamily;
	else if ( r == "suborder" )			nodeRank [node] = Suborder;
	else if ( r == "class" )			nodeRank [node] = Class;
	else if ( r == "species group" )	nodeRank [node] = SpeciesGroup;
	else if ( r == "forma" )			nodeRank [node] = Forma;
	else if ( r == "subtribe" )			nodeRank [node] = Subtribe;
	else if ( r == "subclass" )			nodeRank [node] = Subclass;
	else if ( r == "species subgroup" )	nodeRank [node] = SpeciesSubgroup;
	else if ( r == "phylum" )			nodeRank [node] = Phylum;
	else if ( r == "infraorder" )		nodeRank [node] = Infraorder;
	else if ( r == "superorder" )		nodeRank [node] = Superorder;
	else if ( r == "subphylum" )		nodeRank [node] = Subphylum;
	else if ( r == "infraphylum" )		nodeRank [node] = Infraphylum;
	else if ( r == "parvorder" )		nodeRank [node] = Parvorder;
	else if ( r == "superclass" )		nodeRank [node] = Superclass;
	else if ( r == "superphylum" )		nodeRank [node] = Superphylum;
	else if ( r == "superkingdom" )		nodeRank [node] = Superkingdom;
	else if ( r == "kingdom" )			nodeRank [node] = Kingdom;
	else if ( r == "subkingdom" )		nodeRank [node] = Subkingdom;

	else if ( r == "infrakingdom" )	nodeRank [node] = Infrakingdom;	// Not used as yet
	else if ( r == "microphylum" )	nodeRank [node] = Microphylum;
	else if ( r == "infraclass" )	nodeRank [node] = Infraphylum;
	else if ( r == "parvclass" )	nodeRank [node] = Parvclass;
	else if ( r == "magnorder" )	nodeRank [node] = Magnorder;
	else if ( r == "superorder" )	nodeRank [node] = Superorder;
	else if ( r == "superspecies" )	nodeRank [node] = Superspecies;
	else if ( r == "infraspecies" )	nodeRank [node] = Infraspecies;
	else if ( r == "subvarietas" )	nodeRank [node] = Subvarietas;
	else if ( r == "subforma" )		nodeRank [node] = Subforma;

	else nodeRank [node] = Unknown;
}
int NCBINodes::getMergedNode ( int node ) const
{
	if ( mNodes == 0 ) mNodes = new NCBIMergedNodes;
	return mNodes->getNewNode ( node );	// Check if the node is a merged node then get it's new node number
}
void NCBINodes::getChildNodes ( int node, IntVector& cn ) const
{
	MapIntToIntVectorConstIterator cur = childNodes.find ( node );
	if ( cur != childNodes.end () ) {
		const IntVector& iv = (*cur).second;
		for ( IntVectorSizeType i = 0 ; i < iv.size () ; i++ ) {
			cn.push_back ( iv [i] );
			getChildNodes ( iv [i], cn );
		}
	}
}
int NCBINodes::getRank ( int node ) const
{
	MapIntToIntConstIterator cur = nodeRank.find ( node );
	if ( cur == nodeRank.end () ) {
		node = getMergedNode ( node );
		if ( node == 0 ) return -1;
		cur = nodeRank.find ( node );
		if ( cur == nodeRank.end () ) return -1;
	}
	return (*cur).second;
}
string NCBINodes::getRankAsString ( int node ) const
{
	MapIntToIntConstIterator cur = nodeRank.find ( node );
	if ( cur == nodeRank.end () ) {
		node = getMergedNode ( node );
		if ( node == 0 ) return "";
		cur = nodeRank.find ( node );
		if ( cur == nodeRank.end () ) return "";
	}
	return NCBIRanks[(*cur).second];
}
IntVector NCBINodes::getTaxonomy ( int node ) const
{
	IntVector iv;
	iv.push_back ( node );	// This is the species itself
	for ( ; ; ) {
		MapIntToIntConstIterator cur = parent.find ( node );
		if ( cur != parent.end () ) {
			int p = (*cur).second;
			if ( p == 1 ) break;
			iv.push_back ( p );
			node = p;
		}
		else {	// Shouldn't really get here
			node = getMergedNode ( node );
			if ( node == 0 ) break;				// Not a merged node - have to leave
		}
	}
	reverse ( iv.begin (), iv.end () );
	return iv;
}
NCBIMergedNodes::NCBIMergedNodes ()
{
	GenIFStream ifs ( MsparamsDir::instance ().getParamPath ( "taxonomy/merged.dmp" ) );
	string line;
	while ( getline ( ifs, line ) ) {
		char* ptr = const_cast <char*> (line.c_str ());
		int oldNode = atoi ( ptr );
		while ( *ptr && *ptr != '\t' ) ptr++;					// Skip node
		while ( *ptr && ( *ptr == '\t' || *ptr == '|' ) ) ptr++;// Get to new node
		int newNode = atoi ( ptr );
		mergedNodes [oldNode] = newNode;
	}
}
int NCBIMergedNodes::getNewNode ( int node ) const
{
	MapIntToIntConstIterator cur = mergedNodes.find ( node );
	if ( cur != mergedNodes.end () ) {
		return (*cur).second;
	}
	else return 0;		// No new node
}

int TaxonomyMatch::count = 0;	// Number of objects
NCBINodes* TaxonomyMatch::nodes = 0;
UniProtSpeciesMap* TaxonomyMatch::usm = 0;
NCBISpeciesMap* TaxonomyMatch::nsm = 0;

TaxonomyMatch::TaxonomyMatch ( const string& preferred, bool parentFlag, bool childFlag, bool rankFlag ) :
	preferredSpeciesStr ( genToUpper ( preferred ) )
{
	if ( count == 0 ) {
		nodes = new NCBINodes ( parentFlag, childFlag, rankFlag );
		usm = new UniProtSpeciesMap;
	}
	preferredNode = usm->getNode ( preferredSpeciesStr );
	if ( preferredNode == 0 ) {
		if ( nsm == 0 ) nsm = new NCBISpeciesMap;
		preferredNode = nsm->getNode ( preferredSpeciesStr );
		if ( preferredNode == 0 ) throw runtime_error ( "Illegal preferred species." );
	}
	preferredMatchTaxonomy = nodes->getTaxonomy ( preferredNode );
	prefLen = preferredMatchTaxonomy.size ();
	count++;
}
TaxonomyMatch::~TaxonomyMatch ()
{
	count--;
	if ( count == 0 ) {
		delete usm;
		delete nsm;
		delete nodes;
		usm = 0;
		nsm = 0;
		nodes = 0;
	}
}
int TaxonomyMatch::getSpeciesMatchLength ( const string& sp )
{
	MapStringToIntConstIterator cur = matchLengths.find ( sp );
	if ( cur != matchLengths.end () ) {
		return (*cur).second;
	}
	else {
		int node = usm->getNode ( sp );
		if ( node == 0 ) {
			if ( nsm == 0 ) nsm = new NCBISpeciesMap;
			node = nsm->getNode ( sp );
			if ( node == 0 ) return 0;
		}
		IntVector iv = nodes->getTaxonomy ( node );
		size_t ivLen = iv.size ();
		size_t len = genMin ( ivLen, prefLen );
		int mLen = 0;
		for ( IntVectorSizeType i = 0 ; i < len ; i++ ) {
			if ( iv [i] == preferredMatchTaxonomy [i] ) mLen++;
			else break;
		}
		if ( mLen == prefLen && ivLen > prefLen ) {	// child node
			matchLengths [sp] = ivLen;
			return ivLen;
		}
		else {
			matchLengths [sp] = mLen;
			return mLen;
		}
	}
}
int TaxonomyMatch::getBestSpeciesIndex ( const StringVector& sv )
{
	if ( sv.empty () ) return -1;			// Error condition
	if ( sv.size () == 1 ) return 0;

	string curSp = sv [0];
	if ( curSp == preferredSpeciesStr ) return 0;

	int bestMatchIndex = 0;
	int bestLenDiff = prefLen - getSpeciesMatchLength ( curSp );

	for ( StringVectorSizeType i = 1 ; i < sv.size () ; i++ ) {
		curSp = sv [i];
		if ( curSp == preferredSpeciesStr ) return i;		// Perfect match
		if ( bestLenDiff != 0 ) {	// if bestLenDiff is zero just looking for perfect match
			int curDiff = prefLen - getSpeciesMatchLength ( curSp );
			bool flag = false;
			if ( curDiff == 0 ) flag = true;	// Found a match node
			else if ( curDiff < 0 ) {	// Child of the preferred node
				if ( bestLenDiff > 0 || curDiff > bestLenDiff ) flag = true;
			}
			else {		// curDiff > 0 : parent of preferred node
				if ( curDiff < bestLenDiff ) flag = true;
			}
			if ( flag ) {
				bestLenDiff = curDiff;
				bestMatchIndex = i;
			}
		}
	}
	return bestMatchIndex;
}
bool TaxonomyMatch::equal ( const string& sp1, const string& sp2 )
{
	if ( sp1 == sp2 ) return true;
	if ( sp1 == preferredSpeciesStr || sp2 == preferredSpeciesStr ) return false;
	return getSpeciesMatchLength ( sp1 ) == getSpeciesMatchLength ( sp2 );
}
bool TaxonomyMatch::greaterThan ( const string& sp1, const string& sp2 )
{
	if ( sp1 == sp2 ) return false;
	if ( sp1 == preferredSpeciesStr ) return true;
	if ( sp2 == preferredSpeciesStr ) return false;
	int sp1LenDiff = prefLen - getSpeciesMatchLength ( sp1 );
	int sp2LenDiff = prefLen - getSpeciesMatchLength ( sp2 );
	if ( sp1LenDiff == 0 ) {
		return sp2LenDiff != 0;
	}
	else if ( sp1LenDiff < 0 ) {
		return sp2LenDiff > 0 || sp1LenDiff > sp2LenDiff;
	}
	else {
		return sp1LenDiff < sp2LenDiff;
	}
}
void TaxonomyMatch::getNodes ( IntVector& n, bool children ) const
{
	n.push_back ( preferredNode );
	if ( children ) nodes->getChildNodes ( preferredNode, n );
}
