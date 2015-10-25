/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mut_mtrx.h                                                 *
*                                                                             *
*  Created    : July 21st 1996                                                *
*                                                                             *
*  Purpose    : Functions to deal with mutations.                             *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_mut_mtrx_h
#define __lu_mut_mtrx_h

#include <set>
#include <map>
#include <string>
#include <lu_mass.h>
#include <lu_usermod.h>

struct ExtraMod;
typedef std::vector <ExtraMod*> ExtraModVector;

class PeptideSequence {
	std::string sequence;
	ExtraModVector mods;
	ExtraMod* getMod ( char code ) const;
	std::string getFullSequence () const;
	std::string getModSequence ( const std::string& dbSequence, double err ) const;
	static MapCharToString cMods;
	void printCGI ( std::ostream& os, const std::string& s, double err, int n = 1 ) const;
	std::string getNTerm ( double err = 0.0 ) const;
	std::string getCTerm ( double err = 0.0 ) const;
	std::string getNLoss ( double err = 0.0 ) const;
	bool isNTerm () const;
	bool isCTerm () const;
	bool isNLoss () const;
public:
	PeptideSequence ( const std::string& seq, const ExtraModVector& mods = ExtraModVector ( 0 ) );
	bool operator== ( const PeptideSequence& rhs ) const;
	bool isMassMod () const;
	int getMmod () const;
	unsigned int getMassModAAMask () const;
	void setMassMod ( double mm );
	double getMassMod () const;
	int getMmodIdx () const;
	SetPairIntString getModIndicies () const;
	std::string getSequence () const { return sequence; }
	int getLength () const { return sequence.length (); }
	void applyModifications ( MassType& nTerminusWt, MassType& cTerminusWt ) const;
	int applyModifications ( MassType& nTerminusWt, MassType& cTerminusWt, double diff ) const;
	void putCGI ( std::ostream& os, double err, int n = 1 ) const;
	void printHTML ( std::ostream& os, const std::string& dbSequence, double err ) const;
	void printDelimited ( std::ostream& os, const std::string& dbSequence, double err ) const;
	double getMW () const;
	double getP1 () const;
	void putPeptideStringCGI ( std::ostream& os ) const;
	std::string getPeptideString () const;
	static void addConstMod ( char aa, const std::string& mod );
	double getNLossMass () const;
};
typedef std::vector <PeptideSequence> PeptideSequenceVector;
typedef PeptideSequenceVector::size_type PeptideSequenceVectorSizeType;

class Modification {
	char expectedResidue;
	char modifiedResidue;
	char terminalSpecificity;
	ExtraMod* extraMod;
	double massShift;
	bool mmod;
	friend class sort_modifications;
public:
	Modification ( char expectedResidue, char modifiedResidue, char terminalSpecificity );
	Modification ( const std::string& modAminoAcid, char terminalSpecificity );
	Modification ( const std::string& modAminoAcid, double massShift, char terminalSpecificity, bool mmod );
	~Modification ();
	double getMassShift () const { return massShift; }
	char getExpectedResidue () const { return expectedResidue; }
	char getModifiedResidue () const { return modifiedResidue; }
	ExtraMod* getExtraMod () const { return extraMod; }
	char getTerminalSpecificity () const { return terminalSpecificity; }
	bool getMmod () { return mmod; }
	bool getCation () const;
	bool getMetLoss () const;
	bool getDehydro () const;
	bool getNeutralLoss () { return expectedResidue == '.'; }
	bool getNTerm () { return expectedResidue == 'n'; }
	bool getCTerm () { return expectedResidue == 'c'; }
};

class ModificationParameters {
	std::string prefix;
	StringVector homologyTypes;
	StringVector userMods;
	int maxLevels;
	int maxPeptidePermutations;
	int startNominal;
	int endNominal;
	double defect;
	int maxMMod;
	bool uncleaved;
	StringVector compIon;
	bool modNTermProtein;
	bool modNTerm;
	bool modCTermProtein;
	bool modCTerm;
	bool modNeutralLoss;
	std::string modRangeType;
	int modMaxZ;
public:
	ModificationParameters ( const ParameterList* params, const std::string& prefix );
	void printHTML ( std::ostream& os ) const;
	static void copyToCGI ( std::ostream& os, const ParameterList* params, const std::string& prefix );
	static void copyToHiddenFormJavascriptEntry ( std::ostream& os, const ParameterList* params, const std::string& prefix );
	StringVector getHomologyTypes () const { return homologyTypes; }
	StringVector getUserMods () const { return userMods; }
	bool getAllowErrors () const { return !userMods.empty () || !compIon.empty () || modNTerm || modCTerm || modNeutralLoss || ( !homologyTypes.empty () && homologyTypes [0] != "no errors" && homologyTypes [0] != "parent mass" ); }
	int getMaxLevels () const { return maxLevels; }
	int getMaxPeptidePermutations () const { return maxPeptidePermutations; }
	int getStartNominal () const { return startNominal; }
	int getEndNominal () const { return endNominal; }
	double getDefect () const { return defect; }
	int getMaxMMod () const { return maxMMod; }
	bool getUncleaved () const { return uncleaved; }
	StringVector getCompIon () const { return compIon; }
	bool getModNTermProtein () const { return modNTermProtein; }
	bool getModNTerm () const { return modNTerm; }
	bool getModCTermProtein () const { return modCTermProtein; }
	bool getModCTerm () const { return modCTerm; }
	bool getModNeutralLoss () const { return modNeutralLoss; }
	bool getMassMods () const { return modNTerm || modCTerm || modNeutralLoss || !compIon.empty (); }
	bool getSimpleTagSearch () const;
	std::string getModRangeType () const { return modRangeType; }
	int getModMaxZ () const { return modMaxZ; }
};

typedef std::vector <Modification*> ModificationVector;
typedef ModificationVector::const_iterator ModificationVectorConstIterator;
typedef ModificationVector::size_type ModificationVectorSizeType;
typedef std::pair <ModificationVectorConstIterator,int> ModificationPair;

typedef std::set <Modification*> ModificationSet;
typedef std::vector <ModificationSet> VectorModificationSet;
typedef std::vector <ModificationSet>::size_type VectorModificationSetSizeType;

class MultiModification;
typedef std::vector <MultiModification*> MultiModificationVector;
typedef MultiModificationVector::size_type MultiModificationVectorSizeType;

class ModificationTable {
	double mostNegMassShift;
	double mostPosMassShift;
	double maxParentError;
	ModificationVector modifications;
	MultiModificationVector multiModifications;
	ModificationSet rareMods;
	VectorModificationSet maxMods;
	IntVector maxModsCount;
	IntVector maxModsLimit;
	int maxLevels;
	int maxMMod;
	int numMods;
	int maxSequences;
	bool dehydroFlag;
	bool cationFlag;
	int rare;
	bool maxCh;
	mutable ExtraModVector extraMods;
	mutable int numLevels;
	mutable double nextMultiMod;
	mutable double nextSingleMod;
	double mmodStartExact;
	double mmodEndExact;
	ModificationVector mods;
	void getNextMod ( int start, int level );
	void initialiseHomology ( const StringVector& homologyNames );
	void initialiseUserModHomology ( const StringVector& usermods );
	void initialiseExtraUserModHomology ();
	void initialiseUserDefinedLinkHomology ( const MapPairStringStringExtraUserModInfoConstIterator& iter );
	void initialiseMassOffsetModifications ( const ModificationParameters& modificationParameters );
	void initialiseHomologyMultiMod ();
	void addMultiMutatedSequences ( double startMass, double endMass, const std::string& peptide, bool nTermPeptide, bool cTermPeptide, PeptideSequenceVector& sequenceList, int z ) const;
	void getNextModification ( int modIndex, const CharVector& expectedResidues, const CharVector& modificationResidues, const CharVector& terminalSpecificities, const std::string& peptide, bool nTermPeptide, bool cTermPeptide, PeptideSequenceVector& sequenceList, int z ) const;
	bool checkDehydro ( ModificationVectorSizeType lev );
	bool checkMaxAndRare ( ModificationVectorSizeType lev );
	void addSequences ( const std::string& peptide, const CharVector& modificationResidues, PeptideSequenceVector& sequenceList ) const;

	mutable int numModResidues;
	mutable IntVector indexSet;
	mutable std::string newPeptide;
	mutable SetString sequenceSet;
	mutable int pepLen;

	static bool modMapInitialised;
	static bool massMods;
	static void initialiseModMap ();
	static unsigned int MAX_MODIFICATIONS;
	static StringVector names;
	static MapStringToStringVector modMap;
	static char enzTermSpec;

	int mmod;
	int cat;
	int nT;
	int cT;
	int nTermSpec;
	int cTermSpec;
	int neutralLoss;
	int metLoss;
public:
	ModificationTable ( const ModificationParameters& modificationParameters );
	~ModificationTable ();
	bool checkOffset ( double startMass, double endMass ) const;
	bool checkMultiOffset ( double startMass, double endMass ) const;
	bool checkMatch ( double endMass ) const
	{
		if ( endMass >= nextSingleMod ) return true;
		if ( maxLevels > 1 && endMass >= nextMultiMod ) return true;
		return false;
	}
	ModificationPair getPossibleModifications ( double startMass, double endMass );
	bool getMutatedSequences ( double startMass, double endMass, const std::string& peptide, bool nTermPeptide, bool cTermPeptide, PeptideSequenceVector& sequenceList, int z = 0 ) const;
	double getMostNegMassShift () const { return mostNegMassShift; }
	double getMostPosMassShift () const { return mostPosMassShift; }
	double getMaxParentError () const { return maxParentError; }
	static StringVector getNames ()
	{
		if ( modMapInitialised == false ) initialiseModMap ();
		return names;
	}
	static StringVector getLines ( const std::string& name );
	static bool getMassMods ()
	{
		return massMods;
	}
	void resetNextMods () const;
};

#endif /* ! __lu_mut_mtrx_h */
