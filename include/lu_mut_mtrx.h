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
#include <lgen_reg_exp.h>
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

class Motif {
	std::string motif;
	int offset;
	RegularExpression* re;
	bool on;
	IntVector sites;
	int startIndex;
	int stAA;
public:
	Motif ( const std::string& motif, int offset );
	~Motif ();
	bool operator< ( const Motif* m ) const
	{
		if ( m->motif == motif ) {
			return m->offset < offset;
		}
		else {
			return m->motif < motif;
		}
	}
	bool getOn () const { return on; }
	bool getOn ( int s ) const;
	void setSites ( const char* frame );
	void setFlag ( int startAA, int endAA );
	std::string getMotif () const { return motif; }
	int getOffset () const { return offset; }
};

typedef std::map <char, Motif*> MotifMap;
typedef MotifMap::const_iterator MotifMapConstIterator;

class Modification {
	char expectedResidue;
	char modifiedResidue;
	char terminalSpecificity;
	ExtraMod* extraMod;
	double massShift;
	bool mmod;
	friend class sort_modifications;
	static bool singleLetterMotifFlag;	// Deals with sty and m
	static MotifMap motifMap;
public:
	Modification ( char expectedResidue, char modifiedResidue, char terminalSpecificity, Motif* motif = 0 );
	Modification ( const std::string& modAminoAcid, char terminalSpecificity, Motif* motif = 0 );
	Modification ( const std::string& modAminoAcid, double massShift, char terminalSpecificity, bool mmod, Motif* motif = 0 );
	~Modification ();
	double getMassShift () const { return massShift; }
	char getExpectedResidue () const { return expectedResidue; }
	char getModifiedResidue () const { return modifiedResidue; }
	ExtraMod* getExtraMod () const { return extraMod; }
	char getTerminalSpecificity () const { return terminalSpecificity; }
	bool getMmod () { return mmod; }
	bool getOn () const;
	bool getOn ( int site ) const;
	bool getCation () const;
	bool getMetLoss () const;
	bool getDehydro () const;
	bool getNeutralLoss () { return expectedResidue == '.'; }
	bool getNTerm () { return expectedResidue == 'n'; }
	bool getCTerm () { return expectedResidue == 'c'; }
	std::string getModification () const;
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
	bool modRare;
	std::string modRangeType;
	int modMaxZ;
	StringVector modMotifAA;
	IntVector modMotifOffset;
	StringVector modMotif;
	void setModMotifs ( const ParameterList* params );
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
	bool getModRare () const { return modRare; }
	bool getMassMods () const { return modNTerm || modCTerm || modNeutralLoss || !compIon.empty (); }
	bool getSimpleTagSearch () const;
	std::string getModRangeType () const { return modRangeType; }
	int getModMaxZ () const { return modMaxZ; }
	int getNumModMotif () const { return modMotif.size (); }
	std::string getModMotifAA ( int i ) const { return modMotifAA [i]; }
	int getModMotifOffset ( int i ) const { return modMotifOffset [i]; }
	std::string getModMotif ( int i ) const { return modMotif [i]; }
};

typedef std::vector <Modification*> ModificationVector;
typedef ModificationVector::const_iterator ModificationVectorConstIterator;
typedef ModificationVector::size_type ModificationVectorSizeType;
typedef std::pair <ModificationVectorConstIterator,int> ModificationPair;

typedef std::set <Modification*> ModificationSet;
typedef ModificationSet::const_iterator ModificationSetConstIterator;
typedef std::vector <ModificationSet> VectorModificationSet;
typedef std::vector <ModificationSet>::size_type VectorModificationSetSizeType;

typedef std::map <Modification*, int> MapModificationPtrInt;
typedef MapModificationPtrInt::iterator MapModificationPtrIntIterator;
typedef MapModificationPtrInt::const_iterator MapModificationPtrIntConstIterator;

typedef std::vector <Motif*> MotifPtrVector;

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
	SetPairStringInt motifSet;
	MotifPtrVector motifVector;
	VectorModificationSet labelMods;
	MapModificationPtrInt labelModsMap;
	int maxLevels;
	int maxMMod;
	int numMods;
	int maxSequences;
	bool dehydroFlag;
	bool cationFlag;
	int rare;
	bool modRare;
	bool maxCh;
	bool label;
	mutable ExtraModVector extraMods;
	mutable int numLevels;
	mutable double nextMultiMod;
	mutable double nextSingleMod;
	double mmodStartExact;
	double mmodEndExact;
	ModificationVector mods;
	std::map <std::string, ModificationVector> siteInfoProtein;
	std::map <std::string, ModificationVector> siteInfoPeptide;
	std::map <std::string, ModificationVector> siteInfoSite;
	void getNextMod ( int start, int level );
	void initialiseHomology ( const StringVector& homologyNames );
	void initialiseUserModHomology ( const StringVector& usermods );
	Motif* getMotifPointer ( const std::string& motif, int motifOffset );
	void initialiseExtraUserModHomology ();
	void initialiseUserDefinedLinkHomology ( const MapPairStringStringExtraUserModInfoConstIterator& iter );
	void initialiseMassOffsetModifications ( const ModificationParameters& modificationParameters );
	void initialiseHomologyMultiMod ();
	void addMultiMutatedSequences ( double startMass, double endMass, const std::string& peptide, bool nTermPeptide, bool cTermPeptide, PeptideSequenceVector& sequenceList, int z ) const;
	void getNextModification ( int modIndex, const ModificationVector& mods, const CharVector& modificationResidues, const std::string& peptide, bool nTermPeptide, bool cTermPeptide, PeptideSequenceVector& sequenceList, int z ) const;
	bool checkDehydro ( ModificationVectorSizeType lev );
	bool checkMaxRareAndLabel ( ModificationVectorSizeType lev );
	void addSequence ( PeptideSequenceVector& sequenceList, Modification* m ) const					// Single mod
	{
		if ( label ) {
			MapModificationPtrIntConstIterator cur = labelModsMap.find ( m );					// Find the modification
			if ( cur != labelModsMap.end () && rejectUnmodified ( (*cur).second - 1 ) ) return;
		}
		sequenceList.push_back ( PeptideSequence ( newPeptide, extraMods ) );
	}
	void addSequence ( PeptideSequenceVector& sequenceList, const ModificationVector& mods ) const	// Multiple mods
	{
		if ( label ) {
			for ( ModificationVectorSizeType i = 0 ; i < mods.size () ; i++ ) {
				MapModificationPtrIntConstIterator cur = labelModsMap.find ( mods [i] );		// Find the modification
				if ( cur != labelModsMap.end () && rejectUnmodified ( (*cur).second - 1 ) ) {
					return;
				}
			}
		}
		sequenceList.push_back ( PeptideSequence ( newPeptide, extraMods ) );
	}
	bool rejectUnmodified ( int labIdx ) const
	{
		for ( ModificationSetConstIterator i = labelMods [labIdx].begin () ; i != labelMods [labIdx].end () ; i++ ) {
			if ( newPeptide.find ( (*i)->getExpectedResidue () ) != std::string::npos ) return true;	// Reject if the peptide contains an unmodified residue
		}
		return false;
	}
	void addSequences ( const std::string& peptide, const CharVector& modificationResidues, PeptideSequenceVector& sequenceList, const ModificationVector& mods ) const;

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
	void printModificationCombinations ( std::ostream& os ) const;
	void checkLabelMods ();
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
	void setMotifSites ( const char* frame );
	void setMotifFlags ( int startAA, int endAA );
};

#endif /* ! __lu_mut_mtrx_h */
