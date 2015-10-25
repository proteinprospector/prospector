/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_dig_srch.h                                                 *
*                                                                             *
*  Created    : September 27th 2001                                           *
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

#ifndef __lu_dig_srch_h
#define __lu_dig_srch_h

#include <vector>
#include <lu_program.h>
#include <lu_mod_frag.h>

class MSDigestParameters;
class SingleEntry;

class MSSingleSearch : public MSProgram {
protected:
	const MSDigestParameters& digParams;
	std::vector <SingleEntry*> se;
	std::vector <EnzymeFragmentContainer> enzFrags;
public:
	MSSingleSearch ( const MSDigestParameters& digParams );
	std::vector <SingleEntry*> getSingleEntries () const { return se; }
	virtual void printBodyHTML ( std::ostream& os );
	virtual void printBodyXML ( std::ostream& os );

	virtual void printProteinHTML ( std::ostream& os, int searchIndex, int proteinIndex ) const;
	virtual void printProteinCoverage ( std::ostream& os, const CharVector& aaCovered ) const;
	virtual void printResultsHTML ( std::ostream& os, int i ) {}
	virtual void printResultsXML ( std::ostream& os, int i ) {}
	virtual void printParamsBodyHTML ( std::ostream& os ) const {}
	virtual void printHeaderDelimited ( std::ostream& os ) {}
	virtual void printResultsDelimited ( std::ostream& os, int i ) {}
	void printBodyTabDelimitedText ( std::ostream& os );
};

class MSDigestSearch : public MSSingleSearch {
	std::vector <PotentialMSFragmentContainer> potentialMSFragments;
public:
	MSDigestSearch ( const MSDigestParameters& params );
	void printParamsBodyHTML ( std::ostream& os ) const;
	void printResultsHTML ( std::ostream& os, int i );
	void printResultsDelimited ( std::ostream& os, int i );
	void printResultsXML ( std::ostream& os, int i );
	void printProteinCoverage ( std::ostream& os, const CharVector& aaCovered ) const;
	void printHeaderDelimited ( std::ostream& os );
};

#endif /* ! __lu_dig_srch_h */
