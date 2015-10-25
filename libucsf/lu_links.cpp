/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_links.cpp                                                  *
*                                                                             *
*  Created    : February 2nd 2000                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2000-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef VIS_C
#include <stdexcept>
#endif
#include <lgen_error.h>
#ifdef CHEM_SCORE
#include <lu_chem_sc.h>
#endif
#include <lu_mat_score.h>
#include <lu_links.h>
#include <lu_param_list.h>
#include <lu_aa_calc.h>
#include <lu_get_link.h>
#include <lu_table.h>
using std::vector;
using std::string;
using std::ostringstream;
using std::ostream;
using std::endl;
using std::sort;
using std::set;
using std::pair;
using std::find;
using std::make_pair;
using std::runtime_error;

/*
aaFragsArray - the outer array is arranged by the number of potential link AA - passed in
*/
/*
The parameter maxLinkAA is the maximum number of amino acids that can be crosslinked in any one peptide in the proteins
of interest. This would potentially increase if, for example, the number of missed cleavages increased.
*/
/*
Typical links combinations for maxLinkMolecules = 3, maxLinkAA = 2:
0    - 1 peptide,  0 link aa
1    - 1 peptide,  1 link aa
11   - 2 peptides, each 1 link aa
111  - 3 peptides, each 1 link aa
2    - 1 peptide,  2 link aa
21   - 2 peptides
211  - 3 peptides
22   - 2 peptides
221  - 3 peptides
222  - 3 peptides
Here the maximum length of the string is 3 and the maximum number is 2
*/
class Links {
protected:
	PotentialLinkFragmentVectorVector aaFragsArray;
	PotentialLinkFragmentVectorVector fragsArray;
	LinkHitsVector linkHitsContainer;
	PotentialLinkFragmentVector currentHit;
	string arrayNumbers;
	int numMolecules;
	int numBridges;
	const Peak* peak;
	double bridgeMass;
	ElementalFormula bridgeFormula;
	double continueMass;
	AACalculator aaCalc;
	ElementalFormula startElementalFormula;
	bool checkNumLinks () const;
	static int bIndex;
public:
	Links ( PotentialLinkFragmentVectorVector aaFragsArray, const Peak* peak, double bridgeMass, ElementalFormula& bridgeFormula, const AACalculator& aaCalc );
	virtual void initialize_link_array ( const string& arrayNumbers );
	virtual void getLinkMatches ( double mass, int level, int start );
	LinkHitsVector getLinkHits () const { return linkHitsContainer; };
	bool isHit () const { return !linkHitsContainer.empty (); };
	string getMoleculeTypeString ();
	static void setBIndex ( int idx ) { bIndex = idx; }
};

class DelimitedLinks : public Links {
	PotentialLinkFragmentVectorVector aaFragsArrayDelimited;
	void initAAFragsArrayDelimited ( const string& knownSequence );
public:
	DelimitedLinks ( PotentialLinkFragmentVectorVector aaFragsArray, const Peak* peak, double bridgeMass, ElementalFormula& bridgeFormula, const AACalculator& aaCalc, const string& knownSequence );
	void initialize_link_array ( const string& arrayNumbers );
};
DelimitedLinks::DelimitedLinks ( PotentialLinkFragmentVectorVector aaFragsArray, const Peak* peak, double bridgeMass, ElementalFormula& bridgeFormula, const AACalculator& aaCalc, const string& knownSequence ) :
	Links ( aaFragsArray, peak, bridgeMass, bridgeFormula, aaCalc )
{
	initAAFragsArrayDelimited ( knownSequence );
}
void DelimitedLinks::initAAFragsArrayDelimited ( const string& knownSequence )
{
	aaFragsArrayDelimited.resize ( aaFragsArray.size () );
	ParameterList pList ( knownSequence, false, false, false );
	IntVector start = pList.getIntVectorValue ( "start" );
	IntVector end = pList.getIntVectorValue ( "end" );
	string peptide = pList.getStringValue ( "peptide" );

	StringVector mods = pList.getStringVectorValue ( "mod" );
	MapStringToInt msi;
	for ( StringVectorSizeType a = 0 ; a < mods.size () ; a++ ) {
		MapStringToIntIterator iter = msi.find ( mods [a] );
		if ( iter != msi.end () )	(*iter).second++;
		else						msi [mods [a]] = 1;
	}

	StringVector eMods = pList.getStringVectorValue ( "emod" );
	MapStringToInt msi2;
	for ( StringVectorSizeType b = 0 ; b < eMods.size () ; b++ ) {
		MapStringToIntIterator iter = msi2.find ( eMods [b] );
		if ( iter != msi2.end () )	(*iter).second++;
		else						msi2 [eMods [b]] = 1;
	}

	IntVector entryNumber = pList.getIntVectorValue ( "entry" );
	if ( !entryNumber.empty () ) {
		for ( IntVectorSizeType i = 0 ; i < entryNumber.size () ; i++ ) {
			entryNumber [i] -= 1;
		}
	}
	else {
		StringVector accNum = pList.getStringVectorValue ( "acc" );
		if ( !accNum.empty () ) {
			for ( StringVectorSizeType i = 0 ; i < accNum.size () ; i++ ) {
				entryNumber.push_back ( PotentialLinkFragment::getEntryNumber ( accNum [i] ) );
				if ( entryNumber.back () == -1 ) {
					ErrorHandler::genError ()->error ( "accession number " + accNum [i] + " not found in accession number list.\n" );
				}
			}
		}
		else
			if ( PotentialLinkFragment::getNumAccessionNumbers () == 1 ) entryNumber.push_back ( 0 );
	}
	if ( !start.empty () ) {
		if ( start.size () != entryNumber.size () ) {
			ErrorHandler::genError ()->error ( "If you specify a start amino acid for one accession or entry number you need to specify one for the others.\nTo disable the filtering set it to -1.\n" );
		}
	}
	if ( !end.empty () ) {
		if ( end.size () != entryNumber.size () ) {
			ErrorHandler::genError ()->error ( "If you specify an end amino acid for one accession or entry number you need to specify one for the others.To disable the filtering set it to -1.\n\n" );
		}
	}
	for ( int i = 0 ; i < aaFragsArray.size () ; i++ ) {
		for ( int j = 0 ; j < aaFragsArray [i].size () ; j++ ) {
			const PotentialLinkFragment& plf = aaFragsArray [i][j];
			if ( !entryNumber.empty () ) {
				bool flag = false;
				for ( int k = 0 ; k < entryNumber.size () ; k++ ) {
					if ( plf.getEntryNumber () == entryNumber [k] ) {
						if ( start.empty () || start [k] == -1 || plf.getStartAA () >= start [k] ) {
							if ( end.empty () || end [k] == -1 || plf.getEndAA () <= end [k] ) {
								flag = true;
								break;
							}
						}
					}
				}
				if ( flag == false ) continue;
			}
			if ( !peptide.empty () && plf.getFragment () != peptide ) continue;
			if ( !msi.empty () && !plf.containsMods ( msi ) ) continue;
			if ( !msi2.empty () && !plf.isExactMod ( msi2 ) ) continue;
			aaFragsArrayDelimited [plf.getNumLinkAA ()].push_back ( plf );
		}
	}
}
void DelimitedLinks::initialize_link_array ( const string& arrayNumbers )
{
	SetInt si;
	for ( int a = 0 ; a < arrayNumbers.length () ; a++ ) {
		PairSetIntIteratorBool psiib = si.insert ( arrayNumbers [a] - '0' );
		if ( psiib.second ) {
			string an = arrayNumbers;
			int temp = an [0];
			an [0] = an [a];
			an [a] = temp;
			int totalAA = 0;
			numMolecules = an.length ();
			fragsArray.resize ( numMolecules );				// Number of linked peptides
			for ( int i = 0 ; i < numMolecules ; i++ ) {
				int numAA = an [i] - '0';					// Number of linkable aa in this peptide
				fragsArray [i] = (i == 0) ? aaFragsArrayDelimited [numAA] : aaFragsArray [numAA];
				totalAA += numAA;							// Total number of linkable aa in the peptide assembly
			}
			int maxBridges = totalAA / 2;
			int minBridges = numMolecules - 1;
			for ( int numBridges = minBridges ; numBridges <= maxBridges ; numBridges++ ) {
				double startMass = ( - ( numMolecules - 1 ) * ( h1 - ELECTRON_REST_MASS ) ) + ( bridgeMass * numBridges );
				ElementalFormula efTemp;
				startElementalFormula = "";
				efTemp = "H1";
				efTemp.multiply ( - ( numMolecules - 1 ) );
				startElementalFormula += efTemp;
				efTemp = bridgeFormula;
				efTemp.multiply ( numBridges );
				startElementalFormula += efTemp;
				this->arrayNumbers = an;
				this->numBridges = numBridges;
				getLinkMatches ( startMass, 0, 0 );
			}
		}
	}
}
int Links::bIndex = -1;

Links::Links ( PotentialLinkFragmentVectorVector aaFragsArray, const Peak* peak, double bridgeMass, ElementalFormula& bridgeFormula, const AACalculator& aaCalc ) :
	aaFragsArray ( aaFragsArray ),
	peak ( peak ),
	bridgeMass ( bridgeMass ),
	bridgeFormula ( bridgeFormula ),
	continueMass ( peak->getMass () + peak->getTolerance () ),
	aaCalc ( aaCalc )
{
}
void Links::initialize_link_array ( const string& arrayNumbers )
{
	int totalAA = 0;
	numMolecules = arrayNumbers.length ();
	fragsArray.resize ( numMolecules );
	for ( int i = 0 ; i < numMolecules ; i++ ) {
		int numAA = arrayNumbers [i] - '0';
		fragsArray [i] = aaFragsArray [numAA];
		totalAA += numAA;
	}
	int maxBridges = totalAA / 2;
	int minBridges = numMolecules - 1;
	for ( int numBridges = minBridges ; numBridges <= maxBridges ; numBridges++ ) {
		double startMass = ( - ( numMolecules - 1 ) * ( h1 - ELECTRON_REST_MASS ) ) + ( bridgeMass * numBridges );
		ElementalFormula efTemp;
		startElementalFormula = "";
		efTemp = "H1";
		efTemp.multiply ( - ( numMolecules - 1 ) );
		startElementalFormula += efTemp;
		efTemp = bridgeFormula;
		efTemp.multiply ( numBridges );
		startElementalFormula += efTemp;
		this->arrayNumbers = arrayNumbers;
		this->numBridges = numBridges;
		getLinkMatches ( startMass, 0, 0 );
	}
}
void Links::getLinkMatches ( double mass, int level, int start )
{
	const PotentialLinkFragmentVector& frags = fragsArray [level];
	for ( PotentialLinkFragmentVectorSizeType i = start ; i < frags.size () ; i++ ) {
		mass += frags [i].getMass ();
		if ( mass > continueMass ) break;
		level += 1;
		currentHit.push_back ( frags [i] );
		if ( level < numMolecules ) {
			int nextStart = (arrayNumbers [level] == arrayNumbers [level - 1] ) ? i : 0;
			getLinkMatches ( mass, level, nextStart );
		}
		else {
			if ( peak->isMatch ( mass ) && checkNumLinks () ) {
				linkHitsContainer.push_back ( LinkHits ( currentHit, peak, mass, getMoleculeTypeString (), aaCalc, startElementalFormula ) );
			}
		}
		level -= 1;
		currentHit.pop_back ();
		mass -= frags [i].getMass ();
	}
}
bool Links::checkNumLinks () const
{
	if ( bIndex == -1 ) return true;	// both ends of the bridge are equivalent.
	int numALinks = 0;
	int numBLinks = 0;
	int numMolecules = currentHit.size ();
	for ( int i = 0 ; i < numMolecules ; i++ ) {
		int n = currentHit [i].getNumLinkAA ( bIndex );
		numALinks += n;
		numBLinks += currentHit [i].getNumLinkAA () - n;
	}
	if ( numALinks < numBridges || numBLinks < numBridges ) return false;
	return true;
}
string Links::getMoleculeTypeString ()
{
	string typeString;
	int len = arrayNumbers.length ();
	for ( int i = 0 ; i < len ; i++ ) {
		typeString += arrayNumbers [i];
		if ( i != len - 1 ) typeString += '-';
	}
	ostringstream strstream;
	strstream << '(' << numBridges << "B" << ')';
	typeString += strstream.str ();

	return typeString;
}

class LinksCombination {
	StringVector arrayNumbers;
	string array;
	StringVectorSizeType maxLinkMolecules;
	void getNextLinksCombination ( int level, int maxLinkAA );
public:
	LinksCombination ( int maxLinkMolecules, int maxLinkAA );
	StringVector getLinksCombinations () const { return arrayNumbers; };
};

LinksCombination::LinksCombination ( int maxLinkMolecules, int maxLinkAA ) :
	maxLinkMolecules ( maxLinkMolecules )
{
	getNextLinksCombination ( 0, maxLinkAA );
}
void LinksCombination::getNextLinksCombination ( int level, int maxLinkAA )
{
	for ( int i = 1 ; i <= maxLinkAA ; i++ ) {
		array.resize ( level+1 );
		array [level] = i + '0';
		arrayNumbers.push_back ( array );
		if ( array.length () < maxLinkMolecules ) {
			getNextLinksCombination ( level + 1, i );
		}
	}
}
LinkHits::LinkHits ( const PotentialLinkFragmentVector& hitFragments, const Peak* peak, double mass, const string& moleculeType, const AACalculator& aaCalc, ElementalFormula& startElementalFormula ) :
	hitFragments ( hitFragments ),
	peakMatch ( peak, mass ),
	moleculeType ( moleculeType )
{
	elementalFormula = startElementalFormula;
	for ( PotentialLinkFragmentVectorConstIterator i = hitFragments.begin () ; i != hitFragments.end () ; i++ ) {
		elementalFormula += i->getElementalFormula ( aaCalc );
	}
}
bool LinkHits::operator!= ( const LinkHits& rhs ) const
{
	if ( peakMatch.getMatchedMass () != rhs.peakMatch.getMatchedMass () ) return true;
	if ( hitFragments [0].getFragment () != rhs.hitFragments [0].getFragment () ) return true;
	return false;
}
bool LinkHits::isHomologyMatch ( const LinkHits& rhs ) const
{
	if ( hitFragments.size () == 1 && rhs.hitFragments.size () == 1 ) {
		const string& frag = hitFragments [0].getFragment ();
		const string& frag2 = rhs.hitFragments [0].getFragment ();

		if ( frag.length () == frag2.length () ) {
			if ( (double)matrix_score ( frag.c_str (), frag2.c_str () ) / (double)frag.length () > 4.0 ) {
				return true;
			}
		}
	}
	return false;
}
double LinkHits::getError ( const PeakMatchContext& peakMatchContext ) const
{
	return peakMatch.getError ( peakMatchContext );
}
void LinkHits::printHTMLHeader ( ostream& os, const PeakMatchContext& peakMatchContext, bool modifications, bool multiMolecules, bool peptideCombination )
{
	tableRowStart ( os );
		peakMatch.printHeaderHTML ( os, peakMatchContext );
		if ( peptideCombination ) {
			tableHeader ( os, "Peptide<br />Combination" );
			tableHeader ( os, "Elemental<br />Composition" );
		}
		PotentialLinkFragment::printHTMLHeader ( os, modifications, multiMolecules );
	tableRowEnd ( os );
}
void LinkHits::printHTML ( ostream& os, const PeakMatchContext& peakMatchContext, bool modifications, bool multiMolecules, bool peptideCombination, const MSProductLink& productLink, const MSIsotopeLink& isotopeLink, bool hideLinks )
{
	bool first = true;
	for ( PotentialLinkFragmentVectorConstIterator i = hitFragments.begin () ; i != hitFragments.end () ; i++ ) {
		tableRowStart ( os );
		peakMatch.printHTML ( os, peakMatchContext, first );
		if ( peptideCombination ) {
			if ( first ) {
				tableDataStart ( os, "", "left" );
					os << moleculeType << endl;
				tableDataEnd ( os );
			}
			else tableEmptyCell ( os );
			if ( first ) {
				tableDataStart ( os, "", "left" );
					isotopeLink.write ( os, elementalFormula, peakMatch.getCharge () );
				tableDataEnd ( os );
			}
			else tableEmptyCell ( os );
		}
		i->printHTML ( os, modifications, multiMolecules, getMultipleFragments (), productLink, peakMatch.getCharge (), hideLinks );
		tableRowEnd ( os );
		first = false;
	}
}
void LinkHits::printXML ( ostream& os, const PeakMatchContext& peakMatchContext, bool modifications, bool peptideCombination )
{
	os << "<hit>" << endl;
	bool first = true;
	for ( PotentialLinkFragmentVectorConstIterator i = hitFragments.begin () ; i != hitFragments.end () ; i++ ) {
		if ( first ) {
			peakMatch.printXML ( os, peakMatchContext );
		}
		if ( peptideCombination && first ) ParameterList::printXML ( os, "molecule_type", moleculeType );
		i->printXML ( os, modifications );
		first = false;
	}
	os << "</hit>" << endl;
}
LinksSearch::LinksSearch ( const vector <EnzymeFragmentContainer>& enzFragInfo, const PeakContainer& peaks, const LinkInfo* linkInfo, const AACalculator& aaCalc, const StringVector knownSequences ) :
	SingleFitSearch ( peaks )
{
	initLinksSearch ( enzFragInfo, peaks, linkInfo, aaCalc, knownSequences );
}
LinksSearch::LinksSearch ( const EnzymeFragmentContainer& enzFragInfo, const PeakContainer& peaks, const LinkInfo* linkInfo, const AACalculator& aaCalc, const StringVector knownSequences ) :
	SingleFitSearch ( peaks )
{
	vector <EnzymeFragmentContainer> efi;
	efi.push_back ( enzFragInfo );
	initLinksSearch ( efi, peaks, linkInfo, aaCalc, knownSequences );
}
void LinksSearch::initLinksSearch ( const vector <EnzymeFragmentContainer>& enzFragInfo, const PeakContainer& peaks, const LinkInfo* linkInfo, const AACalculator& aaCalc, const StringVector knownSequences )
{
	peptideCombinationFlag = !linkInfo->getNoLink ();
	bool monoisotopicFlag = peaks.getMonoisotopicFlag ();
	if ( monoisotopicFlag ) mass_convert = formula_to_monoisotopic_mass;
	else mass_convert = formula_to_average_mass;
	PotentialLinkFragmentVectorVector aaFragsArray;
	ElementalFormula bridgeFormula;
	double bridgeMass = 0.0;
	multiMolecules = ( enzFragInfo.size () > 1 );
	FragModContainer::setMassType ( monoisotopicFlag );
	int maxLinkMolecules = linkInfo->getMaxLinkMolecules ();
	if ( linkInfo->getNoLink () ) {			// MS-Fit case - no linked molecules.
		maxLinkMolecules = 0;
	}
	else {
		string searchType = linkInfo->getName ();
		bridgeMass = mass_convert ( linkInfo->getBridgeFormula ().c_str () );
		bridgeFormula = linkInfo->getBridgeFormula ();
	}
	double maxMass = peaks.getMaxMassPlusTol ();
	for ( EnzymeFragmentContainerVectorSizeType ii = 0 ; ii < enzFragInfo.size () ; ii++ ) {
		calculatePotentialFragments ( enzFragInfo [ii], ii, maxMass );
	}
	for ( PotentialLinkFragmentVectorConstIterator i = potentialLinkFrags.begin () ; i != potentialLinkFrags.end () ; i++ ) {
		PotentialLinkFragmentVectorVectorSizeType numLinkAA = i->getNumLinkAA ();
		if ( numLinkAA + 1 >= aaFragsArray.size () ) aaFragsArray.resize ( numLinkAA + 1 );
		aaFragsArray [numLinkAA].push_back ( *i );
	}
	int maxLinkAA = aaFragsArray.size () - 1;
/*
The parameter maxLinkAA is the maximum number of amino acids that can be crosslinked in any one peptide in the proteins
of interest. This would potentially increase if, for example, the number of missed cleavages increased.
*/
	for ( int x = 0 ; x <= maxLinkAA ; x++ ) {
		sort ( aaFragsArray[x].begin (), aaFragsArray[x].end (), sort_fragments_by_mass () );
	}
	LinksCombination lc ( maxLinkMolecules, maxLinkAA );
	StringVector arrayNumbers = lc.getLinksCombinations ();
/*
Typical links combinations for maxLinkMolecules = 3, maxLinkAA = 2:
0
1
11
111
2
21
211
22
221
222
Here the maximum length of the string is 3 and the maximum number is 2
*/
	if ( !aaFragsArray.empty () ) {
		if ( linkInfo->getNoLink () ) {			// MS-Fit case - no linked molecules.
			for ( PeakContainerSizeType a = 0 ; a < peaks.size () ; a++ ) {
				bool hitFlag = false;
				const Peak* pk = peaks [a];
				for ( PotentialLinkFragmentVectorSizeType b = 0 ; b < potentialLinkFrags.size () ; b++ ) {
					double mass = potentialLinkFrags [b].getMass ();
					if ( pk->isMatch ( mass ) ) {
						PotentialLinkFragmentVector currentHit;
						currentHit.push_back ( potentialLinkFrags [b] );
						ElementalFormula ef;
						linkHits.push_back ( LinkHits ( currentHit, pk, mass, "", aaCalc, ef ) );
						hitFlag = true;
					}
				}
				peakUsed [a] = hitFlag;
			}
		}
		else {
			Links::setBIndex ( linkInfo->getBIndex () );
			for ( PeakContainerSizeType a = 0 ; a < peaks.size () ; a++ ) {
				Links* links;
				if ( knownSequences.empty () || knownSequences [a].empty () )
					links = new Links ( aaFragsArray, peaks [a], bridgeMass, bridgeFormula, aaCalc );
				else
					links = new DelimitedLinks ( aaFragsArray, peaks [a], bridgeMass, bridgeFormula, aaCalc, knownSequences [a] );
				links->initialize_link_array ( "0" );
				for ( StringVectorSizeType b = 0 ; b < arrayNumbers.size () ; b++ ) {
					links->initialize_link_array ( arrayNumbers [b] );
				}
				LinkHitsVector lh = links->getLinkHits ();
				linkHits.insert ( linkHits.end (), lh.begin (), lh.end () );
				peakUsed [a] = links->isHit ();
				delete links;
			}
		}
	}
	for ( LinkHitsVectorSizeType m = 0 ; m < linkHits.size () ; m++ ) {
		errors.push_back ( linkHits [m].getError ( peakMatchContext ) );
	}
	sort ( linkHits.begin (), linkHits.end (), sortLinkHits () );
	calculateNumUnique ();
	coverageMap.resize ( enzFragInfo.size () );
	for ( EnzymeFragmentContainerVectorSizeType f = 0 ; f < enzFragInfo.size () ; f++ ) {
		calculateAACovered ( f, enzFragInfo [f].getProtLen () );
	}
	calculateTIC ( peaks );
	calculateStats ();
}
LinksSearch::~LinksSearch () {}
void LinksSearch::calculateAACovered ( int index, int protLen )
{
	coverageMap [index].resize ( protLen );
	for ( LinkHitsVectorSizeType i = 0 ; i < linkHits.size () ; i++ ) {
		const PotentialLinkFragmentVector& hitFragments = linkHits [i].getHitFragments ();
		for ( PotentialLinkFragmentVectorConstIterator k = hitFragments.begin () ; k != hitFragments.end () ; k++ ) {
			if ( k->getEntryNumber () == index ) {
				try {
					coverageMap [index].setCoverage ( k->getStartAA (), k->getEndAA () );
				}
				catch ( runtime_error e ) {
					ErrorHandler::genError ()->message ( e );
				}
			}
		}
	}
}
void LinksSearch::calculateNumUnique ()
{
	numUnique = 1;
	for ( LinkHitsVectorSizeType i = 1 ; i < linkHits.size () ; i++ ) {
		if ( linkHits [i] != linkHits [i-1] ) {
			numUnique++;
		}
	}
}
void LinksSearch::calculatePotentialFragments ( const EnzymeFragmentContainer& enzFrags, int entryNumber, double maxMass )
{
	set <PairStringInt> uniqFrag; 
	modifications = false;
	for ( EnzymeFragmentIterator efi ( enzFrags ) ; efi.more () ; efi.advance () ) {
		const EnzymeFragment& ef = efi.getEnzFrag ();
		pair <std::set <PairStringInt>::iterator, bool> flag = uniqFrag.insert ( make_pair ( ef.getFragment (), ef.getStartAA () ) );
		if ( flag.second ) {
			double unmodifiedMass = ef.getMass ();
			if ( unmodifiedMass + FragModContainer::getMinMod () <= maxMass ) {
				FragModContainer fmc ( ef.getFragment (), ef.getFirstFragment (), ef.getLastFragment () );
				for ( fmc.first () ; fmc.isDone () ; fmc.next () ) {
					const MapIntToInt& modList = fmc.getModList ();
					if ( !modifications && !modList.empty () ) modifications = true;
					double mass = unmodifiedMass + fmc.getMass ( modList );
					potentialLinkFrags.push_back ( PotentialLinkFragment ( ef, mass, modList, fmc.countBridgeSites (), entryNumber ) );
				}
			}
		}
	}
}
bool LinksSearch::areHits () const
{
	return !linkHits.empty ();
}
int LinksSearch::calculateNumHomologyMatches ( const LinksSearch* ls ) const
{
	int numMatches = 0;
	if ( areHits () && ls->areHits () ) {
		LinkHits lastLinkHit = linkHits [0];
		for ( LinkHitsVectorSizeType i = 0 ; i < linkHits.size () ; i++ ) {
			const LinkHits& currentLinkHit = linkHits [i];
			if ( i == 0 || currentLinkHit != lastLinkHit ) {
				LinkHits lastHomLinkHit = ls->linkHits [0];
				for ( LinkHitsVectorSizeType j = 0 ; j < ls->linkHits.size () ; j++ ) {
					const LinkHits& homLinkHit = ls->linkHits [j];
					if ( j == 0 || homLinkHit != lastHomLinkHit ) {
						if ( currentLinkHit.isHomologyMatch ( homLinkHit ) ) {
							numMatches++;
							break;
						}
					}
					lastHomLinkHit = homLinkHit;
				}
			}
			lastLinkHit = currentLinkHit;
		}
	}
	return numMatches;
}
void LinksSearch::printHTMLBody ( ostream& os, bool hideLinks )
{
	for ( LinkHitsVectorSizeType i = 0 ; i < linkHits.size () ; i++ ) {
		if ( i == 0 ) linkHits [i].printHTMLHeader ( os, peakMatchContext, modifications, multiMolecules, peptideCombinationFlag );
		linkHits [i].printHTML ( os, peakMatchContext, modifications, multiMolecules, peptideCombinationFlag, productLink, isotopeLink, hideLinks );
	}
}
void LinksSearch::printXMLBody ( ostream& os )
{
	if ( linkHits.size () ) {
		os << "<ms_hits>" << endl;
		for ( LinkHitsVectorSizeType i = 0 ; i < linkHits.size () ; i++ ) {
			linkHits [i].printXML ( os, peakMatchContext, modifications, peptideCombinationFlag );
		}
		os << "</ms_hits>" << endl;
	}
}
StringVector PotentialLinkFragment::accessionNumbers;
PotentialLinkFragment::PotentialLinkFragment ( const EnzymeFragment& enzymeFragment, double mass, const MapIntToInt& modList, const pair <int, IntVector>& numLinkAA, int entryNumber ) :
	enzymeFragment ( enzymeFragment ),
	mass ( mass ),
	modList ( modList ),
	numLinkAA ( numLinkAA ),
	entryNumber ( entryNumber )
{
}
ElementalFormula PotentialLinkFragment::getElementalFormula ( const AACalculator& aaCalc ) const
{
	ElementalFormula ef = enzymeFragment.getElementalFormula ( aaCalc );
	ef += FragModContainer::getElementalFormula ( modList );
	return ef;
}
void PotentialLinkFragment::printHTMLHeader ( ostream& os, bool modifications, bool multiMolecules )
{
#ifdef CHEM_SCORE
	if ( chemScore ) os << "<th>Chem<br />Score";
#endif
	if ( modifications )	tableHeader ( os, "Modifications" );
	if ( multiMolecules )	tableHeader ( os, "Number" );
	EnzymeFragment::printHeaderHTML ( os );
}
void PotentialLinkFragment::printHTML ( ostream& os, bool modifications, bool multiMolecules, bool multipleFragments, const MSProductLink& productLink, const int maxCharge, bool hideLinks ) const
{
#ifdef CHEM_SCORE
	if ( chemScore ) {
		os << "<td>";
		double score = chemScore->getChemScore ( enzymeFragment.getFragment (), FragModContainer::getNumOxidation ( modList ) );
		genPrint ( os, score, 1 );
	}
#endif
	if ( modifications ) FragModContainer::printHTML ( os, modList );
	if ( multiMolecules ) {
		tableDataStart ( os, "", "right" );
			if ( accessionNumbers.empty () )	os << entryNumber+1 << endl;
			else								os << accessionNumbers [entryNumber] << endl;
		tableDataEnd ( os );
	}
	bool hide = hideLinks || multipleFragments || FragModContainer::getModified ( modList );
	enzymeFragment.printHTML ( os, hide, productLink, maxCharge );
}
void PotentialLinkFragment::printXML ( ostream& os, bool modifications ) const
{
	os << "<peptide>" << endl;
	if ( modifications ) {
		FragModContainer::printXML ( os, modList );
	}
	if ( accessionNumbers.empty () )	ParameterList::printXML ( os, "entry_number", entryNumber+1 );
	else								ParameterList::printXML ( os, "entry_number", accessionNumbers [entryNumber] );
	enzymeFragment.printXML ( os );
	os << "</peptide>" << endl;
}
int PotentialLinkFragment::getEntryNumber ( const string& acc )
{
	StringVectorIterator iter = find ( accessionNumbers.begin (), accessionNumbers.end (), acc );
	if ( iter != accessionNumbers.end () ) return iter - accessionNumbers.begin ();
	else return -1;
}
bool PotentialLinkFragment::containsMods ( const MapStringToInt& mods ) const
{
	return FragModContainer::containsMods ( modList, mods );
}
bool PotentialLinkFragment::isExactMod ( const MapStringToInt& mods ) const
{
	return FragModContainer::isExactMod ( modList, mods );
}
#ifdef CHEM_SCORE
ChemScore* PotentialLinkFragment::chemScore = 0;
void PotentialLinkFragment::setChemScore ( const string& cysMod, double metOxF )
{
	chemScore = new ChemScore ( cysMod, metOxF );
}
void PotentialLinkFragment::deleteChemScore ()
{
	delete chemScore;
}
#endif
