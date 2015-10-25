/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_fas_sp.cpp                                                 *
*                                                                             *
*  Created    : July 20th 1996                                                *
*                                                                             *
*  Purpose    : Functions to do species pre-search.                           *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef VIS_C
#include <sys/types.h>
#include <unistd.h>
#include <stdexcept>
#endif
#ifdef VIS_C
#include <process.h>
#endif
#include <algorithm>
#include <lg_string.h>
#include <lgen_error.h>
#include <lgen_file.h>
#include <lu_fas_sp.h>
#include <lu_getfil.h>
#include <lu_species.h>
using std::vector;
using std::string;
using std::istringstream;
using std::unique;
using std::getline;
using std::sort;
using std::stable_sort;
using std::set_difference;
using std::back_inserter;
using std::runtime_error;
using std::endl;
using std::copy;
using std::inserter;
using std::map;

#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif

namespace {
bool getTaxonomyNodesFromCache ( const StringVector& speciesList, IntVector& nodes );
bool updateTaxonomyCache ();
}

TaxonomyNodeList::TaxonomyNodeList ()
{
}
TaxonomyNodeList::~TaxonomyNodeList ()
{
	for ( int i = 0 ; i < tm.size () ; i++ ) {
		delete tm [i];
	}
}
void TaxonomyNodeList::getTaxonomyNodes ( const string& species, IntVector& nodes, bool sort )
{
	StringVector sv;
	if ( TaxonomyGroupNames::instance ().check ( species ) )
		sv = TaxonomyGroupNames::instance ().getAliases ( species );
	else
		sv.push_back ( species );

	for ( int i = 0 ; i < sv.size () ; i++ ) {
		string sp = genToUpper ( sv [i] );
		if ( sp == "UNREADABLE" )
			nodes.push_back ( -1 );
		else if ( sp == "UNKNOWN" )
			nodes.push_back ( 0 );
		else {
			try {
				tm.push_back ( new TaxonomyMatch ( sp, false, true ) );
				tm.back ()->getNodes ( nodes, true );
			}
			catch ( runtime_error e ) {
				//ErrorHandler::genError ()->message ( sp );
			}
		}
	}
	if ( sort ) {
		stable_sort ( nodes.begin (), nodes.end () );
		nodes.erase ( unique ( nodes.begin (), nodes.end () ), nodes.end () );		// This is the list of sorted nodes
	}
}
TaxonomySearch::TaxonomySearch ( const string& databaseName, int numDatabaseEntries, const StringVector& speciesList, bool useCache ) :
	numDatabaseEntries ( numDatabaseEntries )
{
	IntVector nodes;
	if ( useCache ) {
		if ( !getTaxonomyNodesFromCache ( speciesList, nodes ) ) {	// Cache needs updating
			if ( updateTaxonomyCache () ) {							// Try to update the cache
				getTaxonomyNodesFromCache ( speciesList, nodes );
			}
			else {													// Update failed, get the nodes the long way
				useCache = false;
			}
		}
	}
	if ( !useCache ) {
		TaxonomyNodeList* tnl = new TaxonomyNodeList;
		for ( int i = 0 ; i < speciesList.size () ; i++ ) {
			tnl->getTaxonomyNodes ( speciesList [i], nodes, false );
		}
		delete tnl;
	}
	stable_sort ( nodes.begin (), nodes.end () );
	nodes.erase ( unique ( nodes.begin (), nodes.end () ), nodes.end () );		// This is the list of sorted nodes
	if ( !nodes.empty () ) {
		char* info = getFileAsCharPtr ( SeqdbDir::instance ().getSeqdbDir () + databaseName + ".tax" );
		int nodIdx = 0;
		int nextNode = nodes [nodIdx];
		bool first = true;
		for ( ; ; ) {
			char* taxC = first ? strtok ( info, " " ) : strtok ( NULL, " " );
			if ( taxC == 0 ) break;
			int tax = atoi ( taxC );
			first = false;
			int numEntries = atoi ( strtok ( NULL, "\n" ) );
			if ( tax > nextNode ) {
				while ( nodIdx < nodes.size () && tax > nodes [nodIdx] ) {
					nodIdx++;
				}
				if ( nodIdx >= nodes.size () ) break;
				nextNode = nodes [nodIdx];
			}
			if ( tax > nextNode ) break;
			if ( tax == nextNode ) {
				for ( int j = 0 ; j < numEntries ; j++ ) {
					indicies.push_back ( atoi ( strtok ( NULL, "\n" ) ) );
				}
				taxList.push_back ( tax );
				taxListSize.push_back ( numEntries );
			}
			else {
				for ( int j = 0 ; j < numEntries ; j++ ) {
					strtok ( NULL, "\n" );
				}
			}
		}
		delete [] info;

		stable_sort ( indicies.begin (), indicies.end () );
		indicies.erase ( unique ( indicies.begin (), indicies.end () ), indicies.end () );
	}
}

TaxonomyNames::TaxonomyNames ()
{
	char* info = getParamsFileInfo ( "taxonomy.txt" );
	bool first = true;
	for ( ; ; ) {
		char* ptr = first ? strtok ( info, "\n" ) : strtok ( NULL, "\n" );
		first = false;
		if ( ptr == 0 ) break;
		names.push_back ( ptr );
	}
	delete [] info;
}
TaxonomyNames& TaxonomyNames::instance ()
{
	static TaxonomyNames d;
	return d;
}
TaxonomyGroupNames::TaxonomyGroupNames ()
{
	int numEntries;
	char* info = getFileInfo ( MsparamsDir::instance ().getParamPath ( "taxonomy_groups.txt" ), '>', 1, true, &numEntries );
	for ( int i = 0 ; i < numEntries ; i++ ) {
		char* n = ( i == 0 ) ? strtok ( info, "\n" ) : strtok ( NULL, "\n" );
		names.push_back ( n );

		StringVector aliasList;
		for ( ; ; ) {
			string speciesName = strtok ( NULL, "\n" );
			if ( speciesName == ">" ) break;
			aliasList.push_back ( speciesName );
		}
		nameAliases [n] = aliasList;
	}
	delete [] info;
}
TaxonomyGroupNames& TaxonomyGroupNames::instance ()
{
	static TaxonomyGroupNames d;
	return d;
}
bool TaxonomyGroupNames::check ( const string& s ) const
{
	return nameAliases.count ( s ) == 1;
}
StringVector TaxonomyGroupNames::getAliases ( const string& s ) const
{
	MapStringToStringVectorConstIterator cur = nameAliases.find ( s );
	if ( cur != nameAliases.end () )
		return (*cur).second;
	else
		return StringVector ();
}
struct CompareSize {
	bool operator () ( const std::pair<string,int>& p1, const std::pair<string,int>& p2 ) const
	{
		if ( p1.second == p2.second )	return p1.first < p2.first;
		else							return p1.second < p2.second;
	}
};
namespace {
bool updateTaxonomyCache ()
{
	string lockedName = MsparamsDir::instance ().getParamPath ( "taxonomy/taxonomy_cache.txt.wlock" );
	string normalName = MsparamsDir::instance ().getParamPath ( "taxonomy/taxonomy_cache.txt" );
	FileList fList ( MsparamsDir::instance ().getParamPath ( "taxonomy" ), "taxonomy_cache", ".rlock", true );
	if ( fList.size () || genFileExists ( lockedName ) ) {	// Locked for reading or writing can't update
		return false;
	}
	else {		// The file is not locked
		if ( !genFileExists ( normalName ) ) {	// Create the file if it doesn't exist
			GenOFStream ofs ( normalName );
		}
		genRename ( normalName, lockedName );		// Lock it
		GenOFStream ofs ( lockedName );
		TaxonomyNodeList tnl;
		StringVector sv = TaxonomyNames::instance ().getList ();
		std::map <std::pair<std::string,int>, IntVector, CompareSize> msi;
		for ( int i = 0 ; i < sv.size () ; i++ ) {
			IntVector nodes;
			tnl.getTaxonomyNodes ( sv [i], nodes, true );
			msi [std::make_pair (sv [i], nodes.size ())] = nodes;
		}
		StringVector sv2 = TaxonomyGroupNames::instance ().getList ();
		for ( int j = 0 ; j < sv2.size () ; j++ ) {
			IntVector nodes;
			tnl.getTaxonomyNodes ( sv2 [j], nodes, true );
			msi [std::make_pair (sv2 [j], nodes.size ())] = nodes;
		}
		for ( std::map <std::pair<std::string,int>, IntVector, CompareSize>::const_iterator ii = msi.begin () ; ii != msi.end () ; ii++ ) {
			ofs << (*ii).first.first << endl;
			int numEntries = (*ii).first.second;
			const IntVector& iv = (*ii).second;
			ofs << "Entries=" << numEntries << endl;
			for ( int ii = 0 ; ii < numEntries ; ii++ ) {
				ofs << iv [ii] << endl;
			}
		}
		ofs.close ();
		genRename ( lockedName, normalName );		// Unlock it
		return true;
	}
}
bool getTaxonomyNodesFromCache ( const StringVector& speciesList, IntVector& nodes )
{
	if ( speciesList.empty () ) return true;				// Cache not required
	string normalName = MsparamsDir::instance ().getParamPath ( "taxonomy/taxonomy_cache.txt" );
	string rLockedName = MsparamsDir::instance ().getParamPath ( "taxonomy/taxonomy_cache.txt." + gen_itoa ( getpid () ) + ".rlock" );
	string mergedPath = MsparamsDir::instance ().getParamPath ( "taxonomy/merged.dmp" );
	string namesPath = MsparamsDir::instance ().getParamPath ( "taxonomy/names.dmp" );
	string nodesPath = MsparamsDir::instance ().getParamPath ( "taxonomy/nodes.dmp" );
	string speclistPath = MsparamsDir::instance ().getParamPath ( "taxonomy/speclist.txt" );
	string taxPath = MsparamsDir::instance ().getParamPath ( "taxonomy.txt" );
	string taxGroupPath = MsparamsDir::instance ().getParamPath ( "taxonomy_groups.txt" );
	SetString remainingSpecies;
	copy ( speciesList.begin (), speciesList.end (), inserter ( remainingSpecies, remainingSpecies.end () ) );
	if ( !genFileExists ( normalName ) ) return false;	// No cache file
	if ( genLastModifyTime ( normalName ) <= genLastModifyTime ( mergedPath ) )		return false;		// Node files newer than cache
	if ( genLastModifyTime ( normalName ) <= genLastModifyTime ( namesPath ) )		return false;
	if ( genLastModifyTime ( normalName ) <= genLastModifyTime ( nodesPath ) )		return false;
	if ( genLastModifyTime ( normalName ) <= genLastModifyTime ( speclistPath ) )	return false;
	if ( genLastModifyTime ( normalName ) <= genLastModifyTime ( taxPath ) )		return false;
	if ( genLastModifyTime ( normalName ) <= genLastModifyTime ( taxGroupPath ) )	return false;
	GenOFStream ofrl ( rLockedName );						// Create a read lock
	ofrl.close ();
	char* info = getParamsFileInfo ( "taxonomy/taxonomy_cache.txt" );
	genUnlink ( rLockedName );								// Remove the read lock
	bool first = true;
	int entLen = strlen ( "Entries=" );
	for ( ; ; ) {
		char* ptr = first ? strtok ( info, "\n" ) : strtok ( NULL, "\n" );
		first = false;
		if ( ptr == 0 ) break;	// End of file
		int numEntries = atoi ( strtok ( NULL, "\n" ) + entLen );	// Get Number of entries
		SetStringIterator iter = remainingSpecies.find ( ptr );
		if ( iter != remainingSpecies.end () ) {					// This is one of the selected entries
			for ( int i = 0 ; i < numEntries ; i++ ) {
				nodes.push_back ( atoi ( strtok ( NULL, "\n" ) ) );	// Add these entries
			}
			remainingSpecies.erase ( iter );						// Get rid of this entry
			if ( remainingSpecies.empty () ) break;
		}
		else {														// Species not required
			for ( int i = 0 ; i < numEntries ; i++ ) {
				strtok ( NULL, "\n" );								// Skip these entries
			}
		}
	}
	delete [] info;
	if ( remainingSpecies.empty () ) return true;			// All entries have been found
	else {
		nodes.clear ();
		return false;
	}
}
}
