/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_acc_num.cpp                                                *
*                                                                             *
*  Created    : April 4th 2001                                                *
*                                                                             *
*  Purpose    : Functions for doing accession number searches.                *
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
#include <algorithm>
#include <lgen_file.h>
#include <lgen_mmap.h>
#include <lgen_sort.h>
#include <lu_check_db.h>
#include <lu_getfil.h>
#include <lu_acc_num.h>
using std::string;
using std::vector;
using std::lower_bound;
using std::ios;

#define MAX_NON_PRINT 31

class NumericAccessionNumberMap : public AccessionNumberMap {
	VectorIndexedInt accessionNumberMap;
	int numEntries;
	NumericAccessionNumberMap& operator= ( NumericAccessionNumberMap& rhs );
	NumericAccessionNumberMap ( const NumericAccessionNumberMap& rhs );
public:
	NumericAccessionNumberMap ( const string& database );
	~NumericAccessionNumberMap ();
	int getIndexNumber ( const string& accessionNumber ) const;
	bool accessionNumberUnique ( const string& accessionNumber ) const;
	virtual void reset () const {}
};

class GeneralAccessionNumberMap : public AccessionNumberMap {
	UIntVector aNums;
	mutable MMapFile <char>* mm;
	string fullPath;
public:
	GeneralAccessionNumberMap ( const string& database );
	~GeneralAccessionNumberMap ();
	int getIndexNumber ( const string& accessionNumber ) const;
	bool accessionNumberUnique ( const string& accessionNumber ) const;
	virtual void reset () const
	{
		if ( mm != 0 ) {
			delete mm;
			mm = 0;
		}
	}
};

class UserAccessionNumberMap : public AccessionNumberMap {
public:
	UserAccessionNumberMap ( const string& database );
	~UserAccessionNumberMap ();
	int getIndexNumber ( const string& accessionNumber ) const;
	bool accessionNumberUnique ( const string& accessionNumber ) const { return true; }
	virtual void reset () const {}
};

class NoAccessionNumberMap : public AccessionNumberMap {
public:
	NoAccessionNumberMap ( const string& database );
	~NoAccessionNumberMap ();
	int getIndexNumber ( const string& accessionNumber ) const;
	bool accessionNumberUnique ( const string& accessionNumber ) const { return true; }
	virtual void reset () const {}
};

AccessionNumberMap* getAccessionNumberMap ( const string& database )
{
	if ( database.empty () || database == "User Protein Random" || database == "User Protein Reverse" )
		return new NoAccessionNumberMap ( database );
	else if ( is_user_protein_database ( database ) )
		return new UserAccessionNumberMap ( database );
	else {
		if ( is_numeric_acc_number_database ( database ) )
			return new NumericAccessionNumberMap ( database );
		else
			return new GeneralAccessionNumberMap ( database );
	}
}

AccessionNumberMap::~AccessionNumberMap () {}

NumericAccessionNumberMap::NumericAccessionNumberMap ( const string& database )
{
	string fullPath = SeqdbDir::instance ().getSeqdbDir () + database + ".acn";
	if ( genFileExists ( fullPath ) ) {
		numEntries = (int) ( genFileSize ( fullPath ) ) / sizeof (IndexedInt);

		accessionNumberMap.resize ( numEntries );

		if ( numEntries ) {
			GenIFStream ist ( fullPath, ios::binary );
			ist.read ( (char*)&accessionNumberMap[0], numEntries * sizeof (IndexedInt) );
		}
	}
}
NumericAccessionNumberMap::~NumericAccessionNumberMap () {}

int NumericAccessionNumberMap::getIndexNumber ( const string& accessionNumber ) const
{
	if ( accessionNumberMap.empty () ) return -1;
	int number = atoi ( accessionNumber.c_str () );
	VectorIndexedIntConstIterator aNum = lower_bound ( accessionNumberMap.begin (), accessionNumberMap.end (), number, SortIndexedIntAscending () );

	if ( aNum == accessionNumberMap.end () || aNum->number != number )
		return -1;
	else
		return aNum->index;
}
bool NumericAccessionNumberMap::accessionNumberUnique ( const string& accessionNumber ) const
{
	return getIndexNumber ( accessionNumber ) == -1;
}
GeneralAccessionNumberMap::GeneralAccessionNumberMap ( const string& database )
{
	int numEntries;
	fullPath = SeqdbDir::instance ().getSeqdbDir () + database + ".acc";
	if ( genFileExists ( fullPath ) ) {

		GenIFStream idiFile ( SeqdbDir::instance ().getSeqdbDir () + database + ".idi", ios::binary );
		idiFile.read ( (char*) &numEntries, sizeof (unsigned int) );
		idiFile.close ();

		aNums.reserve ( numEntries );

		mm = new MMapFile <char> ( fullPath );
		unsigned int startIndex = 0;
		const char* fpointerStart = mm->getStartPointer ();
		const char* fpointerEnd = mm->getEndPointer ();
		while ( *fpointerEnd != '\n' ) fpointerEnd--;
		const char* pointer = fpointerStart;

		for ( int i = 0 ; i < numEntries ; i++ ) {
			aNums.push_back ( mm->getMapOffset ( pointer ) );
			if ( pointer [0] != ' ' ) {
				while ( *++pointer != ' ' );
			}
			char* end;
			strtol ( pointer+1, &end, 10 );						// Skip integer
			pointer = end;
			while ( pointer <= fpointerEnd && *pointer <= MAX_NON_PRINT ) {
				pointer++;
			}
			if ( pointer > fpointerEnd ) {
				if ( i == numEntries - 1 ) break;
				startIndex += pointer - fpointerStart;
				fpointerStart = mm->getRange ( startIndex, startIndex + 0x7ffff );
				fpointerEnd = mm->getEndPointer ();
				while ( *fpointerEnd != '\n' ) fpointerEnd--;
				pointer = fpointerStart;
			}
		}
		reset ();
	}
}
GeneralAccessionNumberMap::~GeneralAccessionNumberMap ()
{
	reset ();
}
class GeneralAccessionNumberSort {
	mutable MMapFile <char>* mm;
	char delim;
public:
	GeneralAccessionNumberSort ( MMapFile <char>* mm, char delim ) :
		mm ( mm ),
		delim ( delim ) {}

	bool operator () ( const unsigned int& a, const char* b ) const
	{
		const char* first = mm->getRange ( a, a + 0x200 );
		return genStrcasecmp ( first, b, delim ) < 0;
	}
};
int GeneralAccessionNumberMap::getIndexNumber ( const string& accessionNumber ) const
{
	if ( mm == 0 ) mm = new MMapFile <char> ( fullPath );
	const char* fpointerStart = mm->getStartPointer ();
	static char delim = ' ';
	string aann = accessionNumber + delim;
	UIntVectorConstIterator an = lower_bound ( aNums.begin (), aNums.end (), aann.c_str (), GeneralAccessionNumberSort ( mm, delim ) );
	if ( an == aNums.end () ) return -1;

	int idx = an - aNums.begin ();

	unsigned int off = aNums [idx];
	const char* a = mm->getRange ( off, off + 0x200 );
	int ret;
	if ( genStrcasecmp ( a, aann.c_str (), delim ) )
		ret = -1;
	else {
		while ( *++a != delim );
		ret = strtol ( a+1, 0, 10 );
	}
	return ret;
}
bool GeneralAccessionNumberMap::accessionNumberUnique ( const string& accessionNumber ) const
{
	return getIndexNumber ( accessionNumber ) == -1;
}
UserAccessionNumberMap::UserAccessionNumberMap ( const string& database )
{
}
UserAccessionNumberMap::~UserAccessionNumberMap () {}
int UserAccessionNumberMap::getIndexNumber ( const string& accessionNumber ) const
{
	return atoi ( accessionNumber.c_str () );
}
NoAccessionNumberMap::NoAccessionNumberMap ( const string& database )
{
}
NoAccessionNumberMap::~NoAccessionNumberMap () {}
int NoAccessionNumberMap::getIndexNumber ( const string& accessionNumber ) const
{
	return -1;
}
