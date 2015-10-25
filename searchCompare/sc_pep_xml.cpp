/******************************************************************************
*                                                                             *
*  Program    : searchCompare                                                 *
*                                                                             *
*  Filename   : sc_pep_xml.cpp                                                *
*                                                                             *
*  Created    : December 21st 2010                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2010-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif

#include <lg_time.h>
#include <lu_aa_info.h>
#include <lu_fas_enz.h>
#include <lu_param_list.h>
#include <lu_usermod.h>
#include <sc_pep_xml.h>
#include <sc_search_res.h>
#include <sc_sres_rep.h>

using std::ostream;
using std::string;

SCPepXMLReport::SCPepXMLReport () :
	pList ( PeptidePosition::getParams () ),
	pipeAnal ( "", "" ),
	analSumm ( "database_refresh", genCurrentXSDDateAndTimeString ( ProteinInfo::getDatabaseTime () ) ),
	runSumm ( "", "raw", ".mzXML" ),
	pxse ( 0 ),
	pxsd ( pList.size () ),
	pxesc ( pList.size () ),
	vpxaam ( pList.size () ),
	vpxtm ( pList.size () )
{
	setSearchSummary ();
	setSampleEnzyme ();
	setSearchDatabase ();
	setEnzymaticSearchConstraint ();
	setConstantMods ();
	setVariableMods ();
}
SCPepXMLReport::~SCPepXMLReport ()
{
	for ( int i1 = 0 ; i1 < vpxaam.size () ; i1++ )	{
		for ( int j1 = 0 ; j1 < vpxaam [i1].size () ; j1++ ) delete vpxaam [i1][j1];
	}
	for ( int i2 = 0 ; i2 < vpxtm.size () ; i2++ ) {
		for ( int j2 = 0 ; j2 < vpxtm [i2].size () ; j2++ ) delete vpxtm [i2][j2];
	}
	for ( int i3 = 0 ; i3 < sSumm.size () ; i3++ )		delete sSumm [i3];
	for ( int i4 = 0 ; i4 < pxesc.size () ; i4++ )		delete pxesc [i4];
	for ( int i5 = 0 ; i5 < pxsd.size () ; i5++ )		delete pxsd [i5];
	delete pxse;
	for ( int i6 = 0 ; i6 < enzymeItems.size () ; i6++ )delete enzymeItems [i6];
}
void SCPepXMLReport::setSearchSummary ()
{
	for ( int i = 0 ; i < pList.size () ; i++ ) {
		string pmt = getPrecursorMassType ( pList [i]->getStringValue ( "parent_mass_convert" ) );
		string fmt = getFragmentMassType ( pList [i]->getStringValue ( "parent_mass_convert" ) );
		sSumm.push_back ( new PepXMLSearchSummary ( "", "Protein Prospector Search Compare", pmt, fmt, i+1 ) );
	}
}
void SCPepXMLReport::setSampleEnzyme ()
{
	string enzymeName;
	for ( int i = 0 ; i < pList.size () ; i++ ) {
		string en = pList [i]->getStringValue ( "enzyme" );
		if ( en == "No enzyme" ) continue;
		if ( !enzymeName.empty () ) {
			if ( en != enzymeName ) {	// If more than one enzyme name used (except No enzyme) the enzyme is ambiguous so don't set it
				enzymeName = "";
				break;
			}
		}
		else {				// If empty then set it
			enzymeName = en;
		}
	}
	if ( !enzymeName.empty () ) {
		StringSizeType start = 0;
		StringSizeType end = 0;
		for ( ; ; ) {											// Deals with multiple enzymes
			end = enzymeName.find_first_of ( "/", start );
			string n = enzymeName.substr ( start, end-start );
			string sense = DigestTable::instance ().getSpecificity ( n ) == 'C' ? "C" : "N";
			string cut = DigestTable::instance ().getBreakMask ( n );
			string noCut = DigestTable::instance ().getExcludeMask ( n );
			enzymeItems.push_back ( new PepXMLSpecificity ( sense, -1, cut, noCut ) );
			if ( end == string::npos ) break;
			start = end + 1;
		}
		pxse = new PepXMLSampleEnzyme ( enzymeItems, enzymeName, "", "", enzymeItems.size () == 1 ? -1 : 0 );
	}
}
void SCPepXMLReport::setSearchDatabase ()
{
	for ( int i = 0 ; i < pList.size () ; i++ ) {
		if ( !ProteinInfo::getDatabasePath ().empty () ) {
			pxsd [i] = new PepXMLSearchDatabase ( ProteinInfo::getDatabasePath (), ProteinInfo::getDNADatabase () ? "NA" : "AA", ProteinInfo::getNumEntries () );
		}
		else
			pxsd [i] = 0;
	}
}
void SCPepXMLReport::setEnzymaticSearchConstraint ()
{
	for ( int i = 0 ; i < pList.size () ; i++ ) {
		string enzymeName = pList [i]->getStringValue ( "enzyme" );
		if ( enzymeName != "No enzyme" ) {
			string ans = pList [i]->getStringValue ( "allow_non_specific" );
			int minNumTermini;
			if ( ans == "at 0 termini" || ans == "N termini-1=D" )								minNumTermini = 2;
			else if ( ans == "at N termini" || ans == "at C termini" || ans == "at 1 termini" ) minNumTermini = 1;
			else																				minNumTermini = 0;
			int maxNumInternalCleavages = pList [i]->getIntValue ( "missed_cleavages" );
			pxesc [i] = new PepXMLEnzymaticSearchConstraint ( enzymeName, maxNumInternalCleavages, minNumTermini );
		}
		else
			pxesc [i] = 0;
	}
}
void SCPepXMLReport::setConstantMods ()
{
	for ( int i = 0 ; i < pList.size () ; i++ ) {
		AAInitInfo aaii ( pList [i] );
		MapStringConstModPtr mscmp = aaii.getConstMods ();
		for ( MapStringConstModPtrConstIterator ii = mscmp.begin () ; ii != mscmp.end () ; ii++ ) {
			ConstMod* cMod = (*ii).second;
			string aaList = cMod->getAAList ();
			cMod->setMass ();
			double cMassDiff = cMod->getMass ();
			for ( StringVectorSizeType iii = 0 ; iii < aaList.length () ; iii++ ) {
				string symbol = string ( 1, aaList [iii] );
				bool term = (symbol == "n" || symbol == "c");
				if ( term ) {
					vpxtm [i].push_back ( new PepXMLTerminalModification ( symbol, cMassDiff, "N", "", cMod->getLongName () ) );
				}
				else {
					vpxaam [i].push_back ( new PepXMLAminoAcidModification ( symbol, cMassDiff, "N", "", cMod->getLongName () ) );
				}
			}
		}
	}
}
void SCPepXMLReport::setVariableMods ()
{
	for ( int i = 0 ; i < pList.size () ; i++ ) {
		StringVector usermods = pList [i]->getStringVectorValue	( "msms_mod_AA" );
		Usermod::initialiseUsermodAAInfo ( usermods );
		for ( StringVectorSizeType jj = 0 ; jj < usermods.size () ; jj++ ) {
			Usermod umod ( usermods [jj] );
			string aaList = umod.getAAList ();
			string peptideTerminus;
			char terminalSpecificity = umod.getTerminalSpecificity ();
			ElementalFormula ef = umod.getElementalFormula ();
			double cMassDiff = formula_to_monoisotopic_mass ( ef );
			for ( StringSizeType jjj = 0 ; jjj < aaList.size () ; jjj++ ) {
				string symbol = string ( 1, aaList [jjj] );
				bool term = (symbol == "n" || symbol == "c");
				if ( term ) {
					string proteinTerminus;
					if ( terminalSpecificity == 'N' ) proteinTerminus = "n";
					if ( terminalSpecificity == 'C' ) proteinTerminus = "c";
					vpxtm [i].push_back ( new PepXMLTerminalModification ( symbol, cMassDiff, "Y", proteinTerminus, umod.getOutputString () ) );
				}
				else {
					string peptideTerminus;
					if ( terminalSpecificity == 'n' || terminalSpecificity == 'N' ) peptideTerminus = "n";
					if ( terminalSpecificity == 'c' || terminalSpecificity == 'C' ) peptideTerminus = "c";
					vpxaam [i].push_back ( new PepXMLAminoAcidModification ( symbol, cMassDiff, "Y", peptideTerminus, umod.getOutputString () ) );
				}
			}
		}
	}
}
string SCPepXMLReport::getPrecursorMassType ( const string& pmc )
{
	string pmt = "monoisotopic";
	if ( pmc == "average" || pmc == "Par(av)Frag(mi)" ) pmt = "average";
	return pmt;
}
string SCPepXMLReport::getFragmentMassType ( const string& pmc )
{
	string fmt = "monoisotopic";
	if ( pmc == "average" || pmc == "Par(mi)Frag(av)" ) fmt = "average";
	return fmt;
}
void SCPepXMLReport::updateFractionName ( const string& fractionName )
{
	runSumm.setBaseName ( fractionName );
	sSumm [0]->setBaseName ( fractionName );
}
void SCPepXMLReport::printHeader ( ostream& os ) const
{
	printPepXMLHeader ( os );
	pipeAnal.printOpenTag ( os, 0 );
	analSumm.print ( os, 0 );
	runSumm.printOpenTag ( os, 0 );
	if ( pxse != 0 ) pxse->print ( os, 0 );

	for ( int i = 0 ; i < pList.size () ; i++ ) {
		sSumm [i]->printOpenTag ( os, 0 );
			if ( pxsd [i] != 0 )	pxsd [i]->print ( os, 1 );
			if ( pxesc [i] != 0 )	pxesc [i]->print ( os, 1 );
			for ( int j = 0 ; j < vpxaam [i].size () ; j++ ) {
				vpxaam [i][j]->print ( os, 1 );
			}
			for ( int k = 0 ; k < vpxtm [i].size () ; k++ ) {
				vpxtm [i][k]->print ( os, 1 );
			}
			pList [i]->pepXMLParameters ( os, 1 );
		sSumm [i]->printCloseTag ( os, 0 );
	}
}
void SCPepXMLReport::printFooter ( ostream& os ) const
{
	runSumm.printCloseTag ( os, 0 );
	pipeAnal.printCloseTag ( os, 0 );
}
SCPepXMLReport2::SCPepXMLReport2 ( const SearchResultsPeptideLine* s, int num ) :
	num ( num ),
	cur ( 0 ),
	subItems2 ( 1 ),
	modInfoSubItems ( 1 )
{
	vsrpl.push_back ( s );
}
SCPepXMLReport2::~SCPepXMLReport2 ()
{
	delete pxsq;
	for ( int i1 = 0 ; i1 < pxsh.size () ; i1++ )			delete pxsh [i1];
	for ( int i2 = 0 ; i2 < subItems4.size () ; i2++ )		delete subItems4 [i2];
	for ( int i4 = 0 ; i4 < modInfoSubItems.size () ; i4++ ) {
		for ( int j4 = 0 ; j4 < modInfoSubItems [i4].size () ; j4++ ) {
			delete modInfoSubItems [i4][j4];
		}
	}
	for ( int i5 = 0 ; i5 < subItems2.size () ; i5++ ) {
		for ( int j5 = 0 ; j5 < subItems2 [i5].size () ; j5++ ) {
			delete subItems2 [i5][j5];
		}
	}
}
void SCPepXMLReport2::add ( const SearchResultsPeptideLine* s )
{
	vsrpl.push_back ( s );
}
void SCPepXMLReport2::init ()
{
	string prevPeptide;
	for ( int i = 0 ; i < vsrpl.size () ; i++ ) {
		const SearchResultsPeptideLine* s = vsrpl [i];
		if ( s->getPeptide () == prevPeptide )
			addAdditionalAccessionNumbers ( s );
		else {
			if ( !prevPeptide.empty () ) {
				subItems2 [cur].insert ( subItems2 [cur].end (), subItems2a.begin (), subItems2a.end () );
				subItems2a.clear ();
				static_cast <PepXMLSearchHit*> (pxsh [pxsh.size ()-1])->updateSubItems ( subItems2 [cur] );
				cur++;
				subItems2.resize ( cur+1 );
				modInfoSubItems.resize ( cur+1 );
			}
			prevPeptide = s->getPeptide ();
			createSearchHit ( s );
		}
	}
	subItems2 [cur].insert ( subItems2 [cur].end (), subItems2a.begin (), subItems2a.end () );
	subItems2a.clear ();
	static_cast <PepXMLSearchHit*> (pxsh [pxsh.size ()-1])->updateSubItems ( subItems2 [cur] );
	createSearchResult ( vsrpl [0] );
	createSpectrumQuery ( vsrpl [0], num );
}
void SCPepXMLReport2::createSearchHit ( const SearchResultsPeptideLine* s )
{
	int hitRank				= s->getRank ();
	string peptide			= s->getDBPeptide ();
	string peptidePrevAA	= s->getPrevAA ();
	string peptideNextAA	= s->getNextAA ();
	string accessionInfo = s->getAccessionInfo ();
	string protein;
	if ( !accessionInfo.empty () )	protein = accessionInfo;
	else							protein = s->getAcc ();
	int numTotProteins		= s->getRepeats ();
	int numMatchedIons		= s->getMatched ();
	int totNumIons = -1;			// This the number of possibly matched ions
	double calcNeutralPepMass	= s->getMCalc ();
	double massdiff				= s->getDaError ();
	int numTolTerm = s->getNumTolTerm ();
	int numMissedCleavages			= s->getMissedCleavages () == "-" ? -1 : atoi ( s->getMissedCleavages ().c_str () );
	int isRejected		= 0;		// According to the standard this isn't used yet.
	string proteinDescr	= s->getName ();
	double calcPI		= s->getProteinPI ();
	double proteinMW	= s->getProteinMW ();

	VectorPairIntDouble vpid;
	s->getModMassesAndIndicies ( vpid );
	for ( int i = 0 ; i < vpid.size () ; i++ ) {
		double m =  AAInfo::getInfo ().getMonoisotopicMass ( peptide[vpid [i].first-1] ) + vpid [i].second;
		modInfoSubItems [cur].push_back ( new PepXMLModAminoacidMass ( vpid [i].first, m ) );
	}
	subItems2a.push_back ( new PepXMLModificationInfo ( modInfoSubItems [cur], s->getModNTermMass (), s->getModCTermMass () ) );

	if ( s->getExpectation () != -1.0 ) subItems2a.push_back ( new PepXMLSearchScore ( "expect", s->getExpectation () ) );
	if ( s->getPValue () != -1.0 )		subItems2a.push_back ( new PepXMLSearchScore ( "pvalue", s->getPValue () ) );
	subItems2a.push_back ( new PepXMLSearchScore ( "disc_score", s->getDiscriminantScore () ) );
	subItems2a.push_back ( new PepXMLSearchScore ( "ion_score", s->getScore () ) );
	subItems2a.push_back ( new PepXMLSearchScore ( "ion_score_diff", s->getScoreDifference () ) );

	pxsh.push_back ( new PepXMLSearchHit ( subItems2a, hitRank, peptide, peptidePrevAA, peptideNextAA, protein,
		numTotProteins, numMatchedIons, totNumIons, calcNeutralPepMass, massdiff, numTolTerm,
		numMissedCleavages, isRejected, proteinDescr, calcPI, proteinMW ) );
}
void SCPepXMLReport2::createSearchResult ( const SearchResultsPeptideLine* s )
{
	pxsr = new PepXMLSearchResult ( pxsh, s->getSearchIndex0 () + 1 );
}
void SCPepXMLReport2::createSpectrumQuery ( const SearchResultsPeptideLine* s, int num )
{
	subItems4.push_back ( pxsr );
	int startScan				= 1;
	int endScan					= 1;
	string msmsInfo				= s->getMSMSInfo ();
	if ( !msmsInfo.empty () ) {
		if ( genStringIsInteger ( msmsInfo ) ) {	// Single number
			startScan = atoi ( msmsInfo.c_str () );
			endScan = startScan;
		}
		else {
			int numComma = gen_strcharcount ( msmsInfo, ',' );
			int numColon = gen_strcharcount ( msmsInfo, ':' );
			if ( numComma ) {
				if ( numComma == numColon ) {	// 2:2852,2857
					int s1 = msmsInfo.find ( ":" ) + 1;
					int e1 = msmsInfo.find ( ",", s1 );
					startScan = atoi ( msmsInfo.substr ( s1, e1 ).c_str () );
					endScan = atoi ( msmsInfo.substr ( e1+1 ).c_str () );
				}
				else {
					if ( numColon == 0 ) {
						int e1 = msmsInfo.find ( "," );
						int st1 = atoi ( msmsInfo.substr ( 0, e1 ).c_str () );
						int en1 = atoi ( msmsInfo.substr ( e1+1 ).c_str () );
						if ( en1 > st1 ) {
							startScan = st1; 
							endScan = en1; 
						}
						// else something like 1550,2
					}
					// else 1464,2:1454,3 or 3136,2:3138,2:3147,2
				}
			}
		}
	}
	double precursorNeutralMass	= s->getM ();
	int assumedCharge			= s->getCharge ();
	double retentionTimeSec		= s->getSpotAsNumber ();
	string activationMethod;
	string spectrum				= s->getFractionName () + '.' + gen_itoa ( startScan, "%05d" ) + '.' + gen_itoa ( endScan, "%05d" ) + '.' + gen_itoa ( assumedCharge );
	pxsq = new PepXMLSpectrumQuery ( subItems4, spectrum, startScan, endScan, precursorNeutralMass,
		assumedCharge, "", num, retentionTimeSec, activationMethod );
}
void SCPepXMLReport2::addAdditionalAccessionNumbers ( const SearchResultsPeptideLine* s )
{
	string accessionInfo = s->getAccessionInfo ();
	string protein;
	if ( !accessionInfo.empty () )	protein = accessionInfo;
	else							protein = s->getAcc ();
	string proteinDescr		= s->getName ();
	int numTolTerm			= s->getNumTolTerm ();
	double proteinMW		= s->getProteinMW ();
	string peptidePrevAA	= s->getPrevAA ();
	string peptideNextAA	= s->getNextAA ();
	subItems2 [cur].push_back ( new PepXMLAlternativeProtein ( protein, proteinDescr, numTolTerm, proteinMW, peptidePrevAA, peptideNextAA ) );
}
void SCPepXMLReport2::print ( ostream& os )
{
	init ();
	pxsq->print ( os, 0 );
}
