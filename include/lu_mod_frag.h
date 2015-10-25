/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mod_frag.h                                                 *
*                                                                             *
*  Created    : July 17th 2001                                                *
*                                                                             *
*  Purpose    : Gets the information on user defined modifications from the   *
*               params/usermod_*.txt files.                                   *
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

#ifndef __lu_mod_frag_h
#define __lu_mod_frag_h

#include <lu_prog.h>
#include <lu_frag_info.h>
#include <lu_prod_par.h>

class ParameterList;
#ifdef CHEM_SCORE
class ChemScore;
#endif

class DigestFragmentParameters {
	double minFragmentMass;
	double maxFragmentMass;
	int minFragmentLength;
public:
	DigestFragmentParameters ( const ParameterList* params );
	double getMinFragmentMass () const { return minFragmentMass; }
	double getMaxFragmentMass () const { return maxFragmentMass; }
	int getMinFragmentLength () const { return minFragmentLength; }
	void printHTML ( std::ostream& os ) const;
};

class PotentialMSFragment {
	EnzymeFragment enzymeFragment;
	MapIntToInt fragMod;
	double miMass;
	double avMass;
	int charge;
	int entryNumber;
#ifdef CHEM_SCORE
	double cScore;
	static ChemScore* chemScore;
#endif
	static bool bullBreeseFlag;
	static bool HPLCFlag;
	static bool monoFlag;
public:
	PotentialMSFragment ( const EnzymeFragment& enzymeFragment, const MapIntToInt& fragMod, int charge, int entryNumber = 0 );
	static void printHeaderHTML ( std::ostream& os, bool modifications );
	static void printHeaderDelimited ( std::ostream& os, bool modifications );
	void printHTML ( std::ostream& os, bool modifications, bool hideLinks, const MSProductLink& productLink, int precision ) const;
	void printDelimited ( std::ostream& os, int precision ) const;
	void printXML ( std::ostream& os, int precision ) const;
	double getMiMass () const { return miMass; }
	int getLength () const { return enzymeFragment.getFragment ().length (); }
	void setAverageMass ();
	friend class SortFragments;
#ifdef CHEM_SCORE
	friend class SortFragmentsByChemScore;
	static void setChemScore ( const std::string& cysMod, double metOxF );
	static void deleteChemScore ();
#endif
	static void setBullBreese ( bool flag );
	static void setHPLC ( bool flag );
};

typedef std::vector <PotentialMSFragment> PotentialMSFragmentVector;
typedef PotentialMSFragmentVector::size_type PotentialMSFragmentVectorSizeType;

#ifdef CHEM_SCORE
class SortFragmentsByChemScore {
	public:
		int operator () ( const PotentialMSFragment& a, const PotentialMSFragment& b ) const
		{
			if ( a.cScore == b.cScore )
				return a.miMass < b.miMass;
			else
				return a.cScore > b.cScore;
		}
};
#endif

class PotentialMSFragmentContainer {
	std::vector <PotentialMSFragment> fragments;
	bool modifications;
	MSProductLink productLink;
	int precision;
	int maxDigestHits;
public:
	PotentialMSFragmentContainer ( const std::vector <EnzymeFragmentContainer>& enzFrags, const DigestFragmentParameters& fragParams, int maxPossCharge, int maxDigestHits );
	PotentialMSFragmentContainer ( const EnzymeFragmentContainer& enzFrags, const DigestFragmentParameters& fragParams, int maxPossCharge, int maxDigestHits, int entryNumber );
	void calculateFragments ( const EnzymeFragmentContainer& enzFrags, int maxPossCharge, int entryNumber, double minMass, double maxMass, int minLength );
	void filterFragments ( const DigestFragmentParameters& fragParams );
	void printHTML ( std::ostream& os, bool hideLinks ) const;
	void printDelimited ( std::ostream& os ) const;
	void printXML ( std::ostream& os ) const;
	void sortFragments ();
};

#endif /* ! __lu_mod_frag_h */
