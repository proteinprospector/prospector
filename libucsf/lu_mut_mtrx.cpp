/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mut_mtrx.cpp                                               *
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
#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif
#include <lg_string.h>
#include <lgen_error.h>
#include <lu_getfil.h>
#include <lu_delim.h>
#include <lu_mut_mtrx.h>
#include <lu_usermod.h>
#include <lu_cgi_val.h>
#include <lu_param_list.h>
#include <lu_aa_info.h>
#include <lu_fas_enz.h>
#include <lu_html.h>
#include <lu_srch_form.h>
using std::ostream;
using std::istringstream;
using std::string;
using std::lower_bound;
using std::sort;
using std::pair;
using std::find;
using std::map;
using std::remove_if;
using std::cout;
using std::make_pair;
using std::fill;

typedef ModificationVector::size_type ModificationVectorSizeType;

class MultiModification {
	CharVector expectedResidues;
	CharVector modifiedResidues;
	CharVector terminalSpecificities;
	ExtraModVector extraMods;
	double massShift;
	friend class sort_multi_modifications;
public:
	MultiModification ( const ModificationVector& mods );
	double getMassShift () const { return massShift; }
	const CharVector& getExpectedResidues () const { return expectedResidues; }
	const CharVector& getModifiedResidues () const { return modifiedResidues; }
	const CharVector& getTerminalSpecificities () const { return terminalSpecificities; }
	const ExtraModVector& getExtraMods () const { return extraMods; }
	bool isExtraMods () const { return !extraMods.empty (); }
};

typedef MultiModificationVector::const_iterator MultiModificationVectorConstIterator;

struct ExtraMod {
	string name;
	char type;
	double mass;
	MassType multiplierMass;
	bool mmod;
	ExtraMod ( char type, double mass, const string& name, bool mmod = false ) :
		type ( type ),
		mass ( mass ),
		name ( name ),
		mmod ( mmod )
	{
		if ( type == 'n' )		this->mass += h1;
		else if ( type == 'c' )	this->mass += o_h;
		multiplierMass = static_cast <MassType> ( this->mass * MASS_TYPE_MULTIPLIER );
	}
	ExtraMod ( char type, const string& nm, double m ) :
		type ( type ),
		mmod ( true )
	{
		mass = m;
		if ( type == 'n' ) {
			mass += h1;
			name = gen_ftoa ( m, "%.4f" );
		}
		else if ( type == 'c' )	{
			mass += o_h;
			name = gen_ftoa ( m, "%.4f" );
		}
		else {
			name = nm.substr ( 0, 2 ) + gen_ftoa ( m, "%.4f" ) + ")";
		}
		multiplierMass = static_cast <MassType> ( mass * MASS_TYPE_MULTIPLIER );
	}
	string getName ( double err = 0.0 ) const
	{
		if ( !mmod ) {
			return name;
		}
		else {
			if ( type == 'u' )
				return name.substr ( 0, 2 ) + gen_ftoa ( atof ( name.substr ( 2, name.length () - 3 ).c_str () ) + err, "%.4f" ) + ")";
			else
				return gen_ftoa ( atof ( name.c_str () ) + err, "%.4f" );
		}
	}
};
typedef ExtraModVector::size_type ExtraModVectorSizeType;
typedef ExtraModVector::iterator ExtraModVectorIterator;

MapCharToString PeptideSequence::cMods;
PeptideSequence::PeptideSequence ( const string& seq, const ExtraModVector& mods ) :
	sequence ( seq ),
	mods ( mods )
{
}
bool PeptideSequence::operator== ( const PeptideSequence& rhs ) const
{
	if ( sequence != rhs.sequence ) return false;
	if ( mods != rhs.mods ) return false;
	return true;
}
unsigned int PeptideSequence::getMassModAAMask () const
{
	string aaMaskArray = "nKCQcADEFGHILMNPRSTVWY";
	for ( ExtraModVectorSizeType i = 0 ; i < mods.size () ; i++ ) {
		const ExtraMod* em = mods [i];
		if ( em->mmod ) {
			int idx;
			if ( em->type == 'n' )		idx = aaMaskArray.find ( 'n' );
			else if ( em->type == 'c' )	idx = aaMaskArray.find ( 'c' );
			else						idx = aaMaskArray.find ( em->name [0] );
			return 1 << idx;
		}
	}
	return 0;
}
void PeptideSequence::setMassMod ( double mm )
{
	for ( ExtraModVectorSizeType i = 0 ; i < mods.size () ; i++ ) {
		if ( mods [i]->mmod ) {
			mods [i] = new ExtraMod ( mods [i]->type, mods [i]->name, mm ); 
			break;
		}
	}
}
bool PeptideSequence::isMassMod () const
{
	for ( ExtraModVectorSizeType i = 0 ; i < mods.size () ; i++ ) {
		if ( mods [i]->mmod ) return true;
	}
	return false;
}
double PeptideSequence::getMassMod () const
{
	double ret = 0.0;
	for ( ExtraModVectorSizeType i = 0 ; i < mods.size () ; i++ ) {
		if ( mods [i]->mmod ) {
			if ( ret ) return 0.0;	// Only 1 mod allowed
			ret = mods [i]->mass;
			if ( mods [i]->type == 'n' ) ret -= h1;
			else if ( mods [i]->type == 'c' ) ret -= o_h; 
		}
	}
	return ret;
}
int PeptideSequence::getMmod () const
{
	int n = 0;
	for ( ExtraModVectorSizeType i = 0 ; i < mods.size () ; i++ ) {
		if ( mods [i]->mmod ) n++;
	}
	return n;
}
int PeptideSequence::getMmodIdx () const
{
	for ( ExtraModVectorSizeType i = 0 ; i < mods.size () ; i++ ) {
		if ( mods [i]->mmod ) {
			char type = mods [i]->type;
			if ( type == 'n' ) return 0;
			if ( type == 'c' ) return getLength ()-1;
			for ( StringSizeType j = 0 ; j < sequence.length () ; j++ ) {
				if ( sequence [j] == type ) return j;
			}
		}
	}
	return -1;
}
void PeptideSequence::applyModifications ( MassType& nTerminusWt, MassType& cTerminusWt ) const
{
	for ( ExtraModVectorSizeType i = 0, j = 0 ; i < mods.size () ; i++ ) {
		char modType = mods [i]->type;
		if ( modType == 'n' )		nTerminusWt = mods [i]->multiplierMass;
		else if ( modType == 'c' )	cTerminusWt = mods [i]->multiplierMass;
		else if ( modType == '.' )	continue;
		else {
			aaArrayMT ['u' + j] = aaArrayMT [mods[i]->name [0]] + mods [i]->multiplierMass;
			j++;
		}
	}
}
int PeptideSequence::applyModifications ( MassType& nTerminusWt, MassType& cTerminusWt, double diff ) const
{
	int mmodIdx;
	int mmodUOffset;
	int nMMods = 0;

	for ( ExtraModVectorSizeType xx = 0, xx2 = 0 ; xx < mods.size () ; xx++ ) {
		if ( mods [xx]->mmod ) {
			mmodIdx = xx;
			mmodUOffset = xx2;
			nMMods++;
		}
		if ( mods [xx]->type == 'u' ) {
			xx2++;
		}
	}
	MassType massMod;
	if ( nMMods != 1 ) mmodIdx = -1;
	else {
		massMod = diff * MASS_TYPE_MULTIPLIER;		// Mass mod is the diff
		char mmodCode = 'u' + mmodUOffset;
		for ( ExtraModVectorSizeType yy = 0 ; yy < mods.size () ; yy++ ) {
			if ( !mods [yy]->mmod ) {
				char modType = mods [yy]->type;
				if ( modType == 'n' ) {
					static MassType adjust = h1 * MASS_TYPE_MULTIPLIER;
					massMod -= mods [yy]->multiplierMass - adjust;
				}
				else if ( modType == 'c' )	{
					static MassType adjust = o_h * MASS_TYPE_MULTIPLIER;
					massMod -= mods [yy]->multiplierMass - adjust;
				}
			}
		}
		for ( StringSizeType zz = 0 ; zz < sequence.length () ; zz++ ) {
			char aa = sequence [zz];
			if ( islower ( aa ) ) {
				if ( aa == 's' || aa == 't' || aa == 'y' ) {
					massMod -= aaArrayMT ['s'] - aaArrayMT ['S'];
				}
				else if ( aa == 'm' ) {
					massMod -= aaArrayMT ['m'] - aaArrayMT ['M'];
				}
				else if ( aa == 'h' ) {
					massMod -= aaArrayMT ['h'] - aaArrayMT ['M'];
				}
				else if ( aa != mmodCode ) {
					massMod -= mods [aa-'u']->multiplierMass;
				}
			}
		}
	}
	for ( ExtraModVectorSizeType i = 0, j = 0 ; i < mods.size () ; i++ ) {
		bool flag = ( i == mmodIdx );
		char modType = mods [i]->type;
		if ( modType == 'n' )		nTerminusWt = flag ? massMod : mods [i]->multiplierMass;
		else if ( modType == 'c' )	cTerminusWt = flag ? massMod : mods [i]->multiplierMass;
		else if ( modType == '.' )	continue;
		else {
			aaArrayMT ['u' + j] = aaArrayMT [mods[i]->name [0]] + ( flag ? massMod : mods [i]->multiplierMass );
			j++;
		}
	}
	return nMMods;
}
/* Old version
void PeptideSequence::applyModifications ( MassType& nTerminusWt, MassType& cTerminusWt, double diff ) const
{
	MassType massMod = diff * MASS_TYPE_MULTIPLIER;
	int mmodIdx;
	int nMMods = 0;
	for ( ExtraModVectorSizeType xx = 0 ; xx < mods.size () ; xx++ ) {
		if ( mods [xx]->mmod ) {
			mmodIdx = xx;
			nMMods++;
		}
		else massMod -= mods [xx]->multiplierMass;
	}
	if ( nMMods != 1 ) mmodIdx = -1;
	for ( ExtraModVectorSizeType i = 0, j = 0 ; i < mods.size () ; i++ ) {
		bool flag = ( i == mmodIdx );
		char modType = mods [i]->type;
		if ( modType == 'n' )		nTerminusWt = flag ? massMod : mods [i]->multiplierMass;
		else if ( modType == 'c' )	cTerminusWt = flag ? massMod : mods [i]->multiplierMass;
		else if ( modType == '.' )	continue;
		else {
			aaArrayMT ['u' + j] = aaArrayMT [mods[i]->name [0]] + ( flag ? massMod : mods [i]->multiplierMass );
			j++;
		}
	}
}
*/
string PeptideSequence::getPeptideString () const
{
	string s;
	string nt = getNTerm ();
	if ( !nt.empty () ) s += nt + '-';
	s += getFullSequence ();
	string ct = getCTerm ();
	if ( !ct.empty () ) s += '-' + ct;
	string nl = getNLoss ();
	if ( !nl.empty () ) s += '+' + nl;
	return s;
}
bool PeptideSequence::isNTerm () const
{
	for ( ExtraModVectorSizeType i = 0 ; i < mods.size () ; i++ ) {
		if ( mods [i]->type == 'n' ) return true;
	}
	return false;
}
bool PeptideSequence::isCTerm () const
{
	for ( ExtraModVectorSizeType i = 0 ; i < mods.size () ; i++ ) {
		if ( mods [i]->type == 'c' ) return true;
	}
	return false;
}
bool PeptideSequence::isNLoss () const
{
	for ( ExtraModVectorSizeType i = 0 ; i < mods.size () ; i++ ) {
		if ( mods [i]->type == '.' ) return true;
	}
	return false;
}
string PeptideSequence::getNTerm ( double err ) const
{
	string s;
	ExtraMod* nTerm = getMod ( 'n' );
	if ( nTerm ) {
		if ( err )	s += nTerm->getName ( err );
		else		s += nTerm->getName ();
	}
	else {							// Check for constant mod
		MapCharToStringConstIterator cur = cMods.find ( 'n' );
		if ( cur != cMods.end () )	s += (*cur).second;
	}
	return s;
}
string PeptideSequence::getCTerm ( double err ) const
{
	string s;
	ExtraMod* cTerm = getMod ( 'c' );
	if ( cTerm ) {
		if ( err )	s += cTerm->getName ( err );
		else		s += cTerm->getName ();
	}
	else {							// Check for constant mod
		MapCharToStringConstIterator cur = cMods.find ( 'c' );
		if ( cur != cMods.end () )	s += (*cur).second;
	}
	return s;
}
string PeptideSequence::getNLoss ( double err ) const
{
	string s;
	ExtraMod* neutralLoss = getMod ( '.' );
	if ( neutralLoss ) {
		if ( err )	s += neutralLoss->getName ( err );
		else		s += neutralLoss->getName ();
	}
	return s;
}
double PeptideSequence::getNLossMass () const
{
	ExtraMod* neutralLoss = getMod ( '.' );
	if ( neutralLoss ) {
		return neutralLoss->mass;
	}
	return 0.0;
}
void PeptideSequence::putPeptideStringCGI ( ostream& os ) const
{
	printCGI ( os, getFullSequence (), 0.0 );
}
void PeptideSequence::printHTML ( ostream& os, const string& dbSequence, double err ) const
{
	err = getMmod () == 1 ? err : 0.0;
	string nt = getNTerm ( err );
	if ( !nt.empty () ) os << nt << '-';
	os << getModSequence ( dbSequence, err );
	string ct = getCTerm ( err );
	if ( !ct.empty () ) os << '-' << ct;
	string nl = getNLoss ( err );
	if ( !nl.empty () ) os << '+' << nl;
}
SetPairIntString PeptideSequence::getModIndicies () const
{
	SetPairIntString si;
	if ( isNTerm () )	si.insert ( make_pair ( -3, string ("0") ) );
	if ( isCTerm () )	si.insert ( make_pair ( -2, string ("0") ) );
	if ( isNLoss () )	si.insert ( make_pair ( -1, string ("0") ) );
	for ( StringSizeType i = 0 ; i < sequence.length () ; i++ ) {
		char aa = sequence [i];
		if ( islower ( sequence [i] ) ) {
			if ( aa >= 'u' && aa <= 'x' ) {
				string s = getMod ( aa )->getName ();
				si.insert ( make_pair ( i, s.substr ( 2, s.length () -3 ) ) );
			}
			else
				si.insert ( make_pair ( i, string ( 1, aa ) ) );
		}
	}
	return si;
}
void PeptideSequence::putCGI ( ostream& os, double err, int n ) const
{
	err = getMmod () == 1 ? err : 0.0;
	string s;
	for ( StringSizeType i = 0 ; i < sequence.length () ; i++ ) {
		char aa = sequence [i];
		if ( aa >= 'u' && aa <= 'x' ) {
			if ( err )	s += getMod ( aa )->getName ( err );
			else		s += getMod ( aa )->getName ();
		}
		else if ( islower ( aa ) )		s += ModCodeMap::instance ().getLowerMod ( aa );
		else {
			MapCharToStringConstIterator cur = cMods.find ( aa );
			if ( cur != cMods.end () )	s += (*cur).second;
			else						s += aa;
		}
	}
	printCGI ( os, s, err, n );
}
void PeptideSequence::printCGI ( ostream& os, const string& s, double err, int n ) const
{
	string num = ( n == 1 ) ? "" : gen_itoa ( n ); 
	printCGIString ( os, "sequence" + num, s );
	string nt = getNTerm ( err );
	if ( !nt.empty () )	printCGIString ( os, "nterm" + num, nt );
	string ct = getCTerm ( err );
	if ( !ct.empty () ) printCGIString ( os, "cterm" + num, ct );
	string nl = getNLoss ( err );
	if ( !nl.empty () ) printCGIString ( os, "nloss" + num, nl );
}
void PeptideSequence::printDelimited ( ostream& os, const string& dbSequence, double err ) const
{
	err = getMmod () == 1 ? err : 0.0;
	string nt = getNTerm ( err );
	if ( !nt.empty () ) delimitedCell ( os, '-' + nt );
	delimitedCell ( os, getModSequence ( dbSequence, err ) );
	string ct = getCTerm ( err );
	if ( !ct.empty () ) delimitedCell ( os, '-' + ct );
	string nl = getNLoss ( err );
	if ( !nl.empty () ) delimitedCell ( os, '+' + nl );
}
string PeptideSequence::getFullSequence () const
{
	string s;
	for ( StringSizeType i = 0 ; i < sequence.length () ; i++ ) {
		char aa = sequence [i];
		if ( aa >= 'u' && aa <= 'x' )		s += aa;
		else if ( islower ( aa ) )			s += ModCodeMap::instance ().getLowerMod ( aa );
		else {
			MapCharToStringConstIterator cur = cMods.find ( aa );
			if ( cur != cMods.end () )	s += (*cur).second;
			else						s += aa;
		}
	}
	return s;
}
string PeptideSequence::getModSequence ( const string& dbSequence, double err ) const
{
	string s;
	for ( StringSizeType i = 0 ; i < dbSequence.length () ; i++ ) {
		char aa = sequence [i];
		if ( aa >= 'u' && aa <= 'x' ) {
			if ( err )	s += getMod ( aa )->getName ( err );
			else		s += getMod ( aa )->getName ();
		}
		else if ( islower ( aa ) )			s += ModCodeMap::instance ().getLowerMod ( aa );
		else if ( dbSequence [i] != aa )	s += ModCodeMap::instance ().getModString ( dbSequence [i], aa );
		else {
			MapCharToStringConstIterator cur = cMods.find ( aa );
			if ( cur != cMods.end () )	s += (*cur).second;
			else						s += aa;
		}
	}
	return s;
}
double PeptideSequence::getMW () const
{
	double nTermWt = n_terminus_wt;
	double cTermWt = c_terminus_wt;
	double neutralLossWt = 0.0;
	for ( ExtraModVectorSizeType i = 0, j = 0 ; i < mods.size () ; i++ ) {
		char modType = mods [i]->type;
		if ( modType == 'n' )		nTermWt += mods [i]->mass - h1;
		else if ( modType == 'c' )	cTermWt += mods [i]->mass - o_h;
		else if ( modType == '.' )	neutralLossWt = mods [i]->mass;
		else {
			amino_acid_wt ['u' + j] = amino_acid_wt [mods[i]->name [0]] + mods [i]->mass;
			j++;
		}
	}
	const char* peptideString = sequence.c_str ();
	double molWt = nTermWt + cTermWt + cation_wt + neutralLossWt;
	while ( *peptideString ) molWt += amino_acid_wt [*peptideString++];
	return molWt;
}
double PeptideSequence::getP1 () const
{
	double nTermWt = n_terminus_wt;
	double cTermWt = c_terminus_wt;
	double neutralLossWt = 0.0;
	for ( ExtraModVectorSizeType i = 0, j = 0 ; i < mods.size () ; i++ ) {
		char modType = mods [i]->type;
		if ( mods [i]->mmod ) {
			if ( modType == 'n' ) {}
			else if ( modType == 'c' ) {}
			else if ( modType == '.' ) {}
			else {
				amino_acid_wt ['u' + j] = amino_acid_wt [mods[i]->name [0]];
				j++;
			}
		}
		else {
			if ( modType == 'n' )		nTermWt += mods [i]->mass - h1;
			else if ( modType == 'c' )	cTermWt += mods [i]->mass - o_h;
			else if ( modType == '.' )	neutralLossWt = mods [i]->mass;
			else {
				amino_acid_wt ['u' + j] = amino_acid_wt [mods[i]->name [0]] + mods [i]->mass;
				j++;
			}
		}
	}
	const char* peptideString = sequence.c_str ();
	double molWt = nTermWt + cTermWt + cation_wt;
	while ( *peptideString ) molWt += amino_acid_wt [*peptideString++];
	return molWt;
}
ExtraMod* PeptideSequence::getMod ( char code ) const
{
	ExtraMod* em = 0;
	for ( ExtraModVectorSizeType i = 0, j = 0 ; i < mods.size () ; i++ ) {
		char modType = mods [i]->type;
		if ( modType == 'u' ) {
			modType += j;
			j++;
		}
		if ( modType == code ) return mods [i];
	}
	return em;
}
void PeptideSequence::addConstMod ( char aa, const string& mod )
{
	string s;
	if ( aa == 'n' || aa == 'c' )	s = mod;
	else							s = string ( 1, aa ) + "(" + mod + ")";
	cMods [aa] = s;
}
Modification::Modification ( char expectedResidue, char modifiedResidue, char terminalSpecificity ) :
	expectedResidue ( expectedResidue ),
	modifiedResidue ( modifiedResidue ),
	terminalSpecificity ( terminalSpecificity ),
	extraMod ( 0 ),
	massShift ( amino_acid_wt [modifiedResidue] - amino_acid_wt [expectedResidue] ),
	mmod ( false )
{
}
Modification::Modification ( const string& modAminoAcid, char terminalSpecificity ) :
	expectedResidue ( modAminoAcid [0] ),
	modifiedResidue ( 'u' ),
	terminalSpecificity ( terminalSpecificity ),
	mmod ( false )
{
	if ( expectedResidue == 'c' ) {
		modifiedResidue = 'c';
		massShift = getCTermMassShift ( modAminoAcid );
		extraMod = new ExtraMod ( 'c', massShift, modAminoAcid.substr ( 2, modAminoAcid.length () - 3 ) );
	}
	else if ( expectedResidue == 'n' ) {
		modifiedResidue = 'n';
		massShift = getNTermMassShift ( modAminoAcid );
		extraMod = new ExtraMod ( 'n', massShift, modAminoAcid.substr ( 2, modAminoAcid.length () - 3 ) );
	}
	else if ( expectedResidue == '.' ) {
		modifiedResidue = '.';
		massShift = getNeutralLossMassShift ( modAminoAcid );
		extraMod = new ExtraMod ( '.', massShift, modAminoAcid.substr ( 2, modAminoAcid.length () - 3 ) );
	}
	else {
		modifyAAInfo ( 'u', modAminoAcid );
		massShift = amino_acid_wt ['u'] - amino_acid_wt [expectedResidue];
		extraMod = new ExtraMod ( 'u', massShift, modAminoAcid );
	}
}
Modification::Modification ( const string& modAminoAcid, double massShift, char terminalSpecificity, bool mmod ) :
	expectedResidue ( modAminoAcid [0] ),
	terminalSpecificity ( terminalSpecificity ),
	massShift ( massShift ),
	mmod ( mmod )
{
	if ( expectedResidue == 'c' ) {
		modifiedResidue = 'c';
		extraMod = new ExtraMod ( 'c', massShift, modAminoAcid.substr ( 2, modAminoAcid.length () - 3 ), mmod );
	}
	else if ( expectedResidue == 'n' ) {
		modifiedResidue = 'n';
		extraMod = new ExtraMod ( 'n', massShift, modAminoAcid.substr ( 2, modAminoAcid.length () - 3 ), mmod );
	}
	else if ( expectedResidue == '.' ) {
		modifiedResidue = '.';
		extraMod = new ExtraMod ( '.', massShift, modAminoAcid.substr ( 2, modAminoAcid.length () - 3 ), mmod );
	}
	else {
		modifiedResidue = 'u';
		extraMod = new ExtraMod ( 'u', massShift, modAminoAcid, mmod );
	}
}
Modification::~Modification ()
{
	delete extraMod;
}
bool Modification::getCation () const
{
	if ( extraMod )	return extraMod->name.find ( "Cation" ) != string::npos;
	else			return false;
}
bool Modification::getMetLoss () const
{
	if ( extraMod )	return extraMod->name == "M(Met-loss)";
	else			return false;
}
bool Modification::getDehydro () const
{
	if ( extraMod )	return extraMod->name.find ( "Dehydro" ) != string::npos;
	else			return false;
}
class sort_modifications {
public:
	bool operator () ( const Modification* a, const Modification* b ) const
	{
		return ( a->massShift < b->massShift );
	}
	bool operator () ( const Modification* a, const double b ) const
	{
		return ( a->massShift < b );
	}
};
class sort_multi_modifications {
public:
	bool operator () ( const MultiModification* a, const MultiModification* b ) const
	{
		return ( a->massShift < b->massShift );
	}
	bool operator () ( const MultiModification* a, const double b ) const
	{
		return ( a->massShift < b );
	}
};
MultiModification::MultiModification ( const ModificationVector& mods )
{
	massShift = 0.0;
	char curMod = 'u';
	int terms = 0;
	for ( ModificationVectorSizeType i = 0 ; i < mods.size () ; i++ ) {
		expectedResidues.push_back ( mods [i]->getExpectedResidue () );
		char mod = mods [i]->getModifiedResidue ();
		ExtraMod* em = mods [i]->getExtraMod ();
		if ( mod == 'u' ) {
			ExtraModVectorIterator emi = find ( extraMods.begin (), extraMods.end (), em );
			if ( emi == extraMods.end () ) {	// New extra mod use a new letter
				modifiedResidues.push_back ( curMod++ );
				extraMods.push_back ( em );
			}
			else {
				int num = emi - extraMods.begin () - terms;		// If there are any terminal groups these don't count
				modifiedResidues.push_back ( 'u' + num );
			}
		}
		else if ( mod == 'n' || mod == 'c' || mod == '.' ) {
			modifiedResidues.push_back ( mod );
			extraMods.push_back ( em );
			terms++;
		}
		else
			modifiedResidues.push_back ( mod );
		terminalSpecificities.push_back ( mods [i]->getTerminalSpecificity () );
		massShift += mods [i]->getMassShift ();
	}
}

ModificationParameters::ModificationParameters ( const ParameterList* params, const string& prefix ) :
	prefix				( prefix ),
	homologyTypes		( params->getStringVectorValue	( prefix + "search_type" ) ),
	userMods			( params->getStringVectorValue	( prefix + "mod_AA" ) ),
	maxLevels			( params->getIntValue			( prefix + "max_modifications", 1 ) ),
	maxPeptidePermutations ( params->getIntValue		( prefix + "max_peptide_permutations", 0 ) ),
	startNominal		( params->getIntValue			( "mod_start_nominal" ) ),
	endNominal			( params->getIntValue			( "mod_end_nominal" ) ),
	defect				( params->getDoubleValue		( "mod_defect" ) ),
	maxMMod				( params->getIntValue			( "mod_max", 1 ) ),
	uncleaved			( params->getBoolValue			( "mod_uncleaved" ) ),
	compIon				( params->getStringVectorValue	( "mod_comp_ion" ) ),
	modNTermProtein		( params->getStringValue		( "mod_n_term_type", "Peptide" ) == "Protein" ),
	modNTerm			( params->getBoolValue			( "mod_n_term", false ) ),
	modCTermProtein		( params->getStringValue		( "mod_c_term_type", "Peptide" ) == "Protein" ),
	modCTerm			( params->getBoolValue			( "mod_c_term", false ) ),
	modNeutralLoss		( params->getBoolValue			( "mod_neutral_loss", false ) ),
	modRangeType		( params->getStringValue		( "mod_range_type", "Da" ) ),
	modMaxZ				( params->getIntValue			( "mod_max_z", 1 ) )
{
	ExtraUserMods::instance ().addUserMods2 ( params, ExtraUserModsForm::getNumUserMods () );
}
void ModificationParameters::printHTML ( ostream& os ) const
{
	ParameterList::printHTMLContainer ( os, "Search Mode", homologyTypes );
	ParameterList::printHTML ( os, "Max Modifications", maxLevels );
}
void ModificationParameters::copyToCGI ( ostream& os, const ParameterList* params, const string& prefix )
{
	params->copyToCGI ( os, prefix + "search_type" );
	params->copyToCGI ( os, prefix + "mod_AA" );
	params->copyToCGI ( os, prefix + "max_modifications" );
}
void ModificationParameters::copyToHiddenFormJavascriptEntry ( ostream& os, const ParameterList* params, const string& prefix )
{
	params->copyToHiddenFormJavascriptEntry ( os, prefix + "search_type" );
	params->copyToHiddenFormJavascriptEntry ( os, prefix + "mod_AA" );
	params->copyToHiddenFormJavascriptEntry ( os, prefix + "max_modifications" );
}
bool ModificationParameters::getSimpleTagSearch () const
{
	return (!getMassMods () && homologyTypes.empty () && userMods.empty () && ExtraUserMods::instance ().empty ()) || maxLevels == 0;
}
class CheckMod {
	string mod;
public:
	CheckMod ( const string& mod ) :
		mod ( mod ) {}
	bool operator () ( const Modification* m )
	{
		ExtraMod* em = m->getExtraMod ();
		if ( em )	return em->name.find ( mod ) != string::npos;
		else		return false;
	}
};

unsigned int ModificationTable::MAX_MODIFICATIONS = 12000000;
bool ModificationTable::modMapInitialised = false;
bool ModificationTable::massMods = false;
map <string, StringVector> ModificationTable::modMap;
StringVector ModificationTable::names;
char ModificationTable::enzTermSpec;

ModificationTable::ModificationTable ( const ModificationParameters& modificationParameters ) :
	maxLevels ( modificationParameters.getMaxLevels () ),
	maxSequences ( modificationParameters.getMaxPeptidePermutations () ),
	maxMMod ( modificationParameters.getMaxMMod () ),
	mostNegMassShift ( 0.0 ),
	mostPosMassShift ( 0.0 ),
	mmodStartExact ( 0.0 ),
	mmodEndExact ( 0.0 )
{
	UpdatingJavascriptMessage ujm;
	ujm.writeMessage ( cout, "Initializing modifications." );
	initialiseHomology ( modificationParameters.getHomologyTypes () );
	initialiseUserModHomology ( modificationParameters.getUserMods () );
	initialiseExtraUserModHomology ();
	initialiseMassOffsetModifications ( modificationParameters );
	sort ( modifications.begin (), modifications.end (), sort_modifications () );
	rare = rareMods.empty () ? 0 : 1;
	maxCh = !maxMods.empty ();
	if ( maxLevels > 1 ) {
		initialiseHomologyMultiMod ();
	}
	int i = remove_if ( modifications.begin (), modifications.end (), CheckMod ( "Dehydro" ) ) - modifications.begin ();	// Dehydros only occur in pairs
	modifications.erase ( modifications.begin () + i, modifications.end () );
	if ( !modifications.empty () ) {
		mostNegMassShift = modifications.front ()->getMassShift ();
		mostPosMassShift = modifications.back ()->getMassShift ();
	}
	if ( !multiModifications.empty () ) {
		mostNegMassShift = genMin ( mostNegMassShift, multiModifications.front ()->getMassShift () );
		mostPosMassShift = genMax ( mostPosMassShift, multiModifications.back ()->getMassShift () );
	}
	maxParentError = genMax ( - mostNegMassShift, mostPosMassShift );
	ujm.deletePreviousMessage ( cout );
}
ModificationTable::~ModificationTable ()
{
	for ( MultiModificationVectorSizeType i = 0 ; i < multiModifications.size () ; i++ ) {
		delete multiModifications [i];
	}
	for ( ModificationVectorSizeType j = 0 ; j < modifications.size () ; j++ ) {
		delete modifications [j];
	}
}
ModificationPair ModificationTable::getPossibleModifications ( double startMass, double endMass )
{
	ModificationPair range;

	range.second = 0;
	if ( startMass >= mostPosMassShift ) {
		return ( range );
	}
	ModificationVectorConstIterator modIter = lower_bound ( modifications.begin (), modifications.end (), startMass, sort_modifications () );
	int number = 0;
	for ( ModificationVectorConstIterator mi = modIter ; mi != modifications.end () ; mi++ ) {
		if ( (*mi)->getMassShift () > endMass ) break;
		number++;
	}
	if ( number ) {
		range.first = modIter;
		range.second = number;
	}
	return ( range );
}
bool ModificationTable::checkOffset ( double startMass, double endMass ) const
{
	if ( startMass < mostPosMassShift ) {
		if ( startMass < 0.0 && endMass > 0.0 ) return true;
		ModificationVectorConstIterator modIter = lower_bound ( modifications.begin (), modifications.end (), startMass, sort_modifications () );
		if ( modIter != modifications.end () && (*modIter)->getMassShift () <= endMass ) {
			return true;
		}
		if ( maxLevels > 1 ) return checkMultiOffset ( startMass, endMass );
		else return false;
	}
	return false;
}
void ModificationTable::resetNextMods () const
{
	if ( !modifications.empty () ) nextSingleMod = modifications.front ()->getMassShift ();
	if ( !multiModifications.empty () ) nextMultiMod = multiModifications.front ()->getMassShift ();
}
struct ModificationTableMaximumSequencesExceeded {};
bool ModificationTable::getMutatedSequences ( double startMass, double endMass, const string& peptide, bool nTermPeptide, bool cTermPeptide, PeptideSequenceVector& sequenceList, int z ) const
{
	ModificationVectorConstIterator modIter = lower_bound ( modifications.begin (), modifications.end (), startMass, sort_modifications () );
	if ( modIter != modifications.end () ) {
		nextSingleMod = (*modIter)->getMassShift ();
		if ( nextSingleMod <= endMass ) {
			newPeptide = peptide;
			pepLen = peptide.length ();
			for ( ModificationVectorConstIterator mi = modIter ; mi != modifications.end () ; mi++ ) {
				if ( (*mi)->getMassShift () > endMass ) break;
				ExtraModVector mal;
				char expectedAA = (*mi)->getExpectedResidue ();
				char terminalSpecificity = (*mi)->getTerminalSpecificity ();
				ExtraMod* extraMod = (*mi)->getExtraMod ();
				if ( extraMod ) mal.push_back ( extraMod );
				if ( expectedAA == '.' ) {			// Neutral loss
					if ( z != 0 && extraMod->mmod ) {
						double mz = (*mi)->getMassShift () / z;
						if ( mz >= mmodStartExact && mz <= mmodEndExact ) {
							sequenceList.push_back ( PeptideSequence ( newPeptide, mal ) );
						}
					}
					else
						sequenceList.push_back ( PeptideSequence ( newPeptide, mal ) );
				}
				else if ( expectedAA == 'n' ) {		// Used for the case of a modified terminal group when you don't
													// care about the amino acid at the terminus.
					if ( terminalSpecificity != 'n' || nTermPeptide ) {		// Covers the case of acetylation
						sequenceList.push_back ( PeptideSequence ( newPeptide, mal ) );
					}
				}
				else if ( expectedAA == 'c' ) {		// Used for the case of a modified terminal group when you don't
													// care about the amino acid at the terminus.
					if ( terminalSpecificity != 'c' || cTermPeptide ) {
						sequenceList.push_back ( PeptideSequence ( newPeptide, mal ) );
					}
				}
				else {
					char modificationAA = (*mi)->getModifiedResidue ();
					int startJ = 0;
					int endJ = pepLen;
					if ( terminalSpecificity != '0' ) {
						if ( terminalSpecificity == 'c' ) startJ = cTermPeptide ? pepLen - 1 : endJ;// Protein C-Terminus
						if ( terminalSpecificity == 'n' ) endJ = nTermPeptide ? 1 : 0;				// Protein N-Terminus
						if ( terminalSpecificity == 'C' ) startJ = pepLen - 1;						// Peptide C-Terminus
						if ( terminalSpecificity == 'N' ) endJ = 1;									// Peptide N-Terminus
						if ( terminalSpecificity == 'e' ) {
							static char enzTermSpec = get_enzyme_terminal_specificity ();
							if ( enzTermSpec == 'C' ) {
								endJ = pepLen - 1;
							}
							if ( enzTermSpec == 'N' ) {
								startJ = 1;
							}
						}
					}
					for ( int j = startJ ; j < endJ ; j++ ) {
						if ( expectedAA == newPeptide [j] ) {
							newPeptide [j] = modificationAA;
							sequenceList.push_back ( PeptideSequence ( newPeptide, mal ) );
							if ( maxSequences && sequenceList.size () > maxSequences ) return false;
							newPeptide [j] = expectedAA; // restores to original sequence
						}
					}
				}
			}
		}
	}
	if ( maxLevels > 1 ) {
		try {
			addMultiMutatedSequences ( startMass, endMass, peptide, nTermPeptide, cTermPeptide, sequenceList, z );
		}
		catch ( ModificationTableMaximumSequencesExceeded ) {
			indexSet.clear ();
			return false;
		}
	}
	return true;
}
bool ModificationTable::checkMultiOffset ( double startMass, double endMass ) const
{
	MultiModificationVectorConstIterator mModIter = lower_bound ( multiModifications.begin (), multiModifications.end (), startMass, sort_multi_modifications () );
	if ( mModIter != multiModifications.end () && (*mModIter)->getMassShift () <= endMass ) {
		return true;
	}
	else {
		return false;
	}
}
void ModificationTable::addMultiMutatedSequences ( double startMass, double endMass, const string& peptide, bool nTermPeptide, bool cTermPeptide, PeptideSequenceVector& sequenceList, int z ) const
{
	MultiModificationVectorConstIterator mModIter = lower_bound ( multiModifications.begin (), multiModifications.end (), startMass, sort_multi_modifications () );
	if ( mModIter != multiModifications.end () ) {
		nextMultiMod = (*mModIter)->getMassShift ();
		if ( nextMultiMod <= endMass ) {
			pepLen = peptide.length ();
			for ( MultiModificationVectorConstIterator mmi = mModIter ; mmi != multiModifications.end () ; mmi++ ) {
				if ( (*mmi)->getMassShift () > endMass ) break;
				numModResidues = (*mmi)->getExpectedResidues ().size ();
				sequenceSet.clear ();
				if ( (*mmi)->isExtraMods () )
					extraMods = (*mmi)->getExtraMods ();
				else
					extraMods.clear ();
				getNextModification ( 0, (*mmi)->getExpectedResidues (), (*mmi)->getModifiedResidues (), (*mmi)->getTerminalSpecificities (), peptide, nTermPeptide, cTermPeptide, sequenceList, z );
			}
		}
	}
}
void ModificationTable::getNextModification ( int modIndex, const CharVector& expectedResidues, const CharVector& modificationResidues, const CharVector& terminalSpecificities, const string& peptide, bool nTermPeptide, bool cTermPeptide, PeptideSequenceVector& sequenceList, int z ) const
{
	char expectedAA = expectedResidues [modIndex];
	char terminalSpecificity = terminalSpecificities [modIndex];
	if ( z != 0 && expectedAA == '.' ) {
		for ( ExtraModVectorSizeType i = 0 ; i < extraMods.size () ; i++ ) {
			const ExtraMod* em = extraMods [i];
			if ( em->type == '.' && em->mmod ) {
				double mz = em->mass / z;
				if ( mz > mmodEndExact || mz < mmodStartExact ) return;
			}
		}
	}
	if ( expectedAA == 'n' ) {
		if ( terminalSpecificity != 'n' || nTermPeptide ) {		// Covers the case of acetylation
			indexSet.push_back ( -1 );
			if ( modIndex + 1 < numModResidues ) getNextModification ( modIndex + 1, expectedResidues, modificationResidues, terminalSpecificities, peptide, nTermPeptide, cTermPeptide, sequenceList, z );
			else {
				addSequences ( peptide, modificationResidues, sequenceList );
			}
			indexSet.pop_back ();
		}
	}
	else if ( expectedAA == 'c' ) {
		if ( terminalSpecificity != 'c' || cTermPeptide ) {
			indexSet.push_back ( -2 );
			if ( modIndex + 1 < numModResidues ) getNextModification ( modIndex + 1, expectedResidues, modificationResidues, terminalSpecificities, peptide, nTermPeptide, cTermPeptide, sequenceList, z );
			else {
				addSequences ( peptide, modificationResidues, sequenceList );
			}
			indexSet.pop_back ();
		}
	}
	else if ( expectedAA == '.' ) {
		indexSet.push_back ( -3 );
		if ( modIndex + 1 < numModResidues ) getNextModification ( modIndex + 1, expectedResidues, modificationResidues, terminalSpecificities, peptide, nTermPeptide, cTermPeptide, sequenceList, z );
		else {
			addSequences ( peptide, modificationResidues, sequenceList );
		}
		indexSet.pop_back ();
	}
	else {
		int startI = 0;
		int endI = pepLen;
		if ( terminalSpecificity != '0' ) {
			if ( terminalSpecificity == 'c' ) startI = cTermPeptide ? pepLen - 1 : endI;
			if ( terminalSpecificity == 'n' ) endI = nTermPeptide ? 1 : 0;
			if ( terminalSpecificity == 'C' ) startI = pepLen - 1;
			if ( terminalSpecificity == 'N' ) endI = 1;
			if ( terminalSpecificity == 'e' ) {
				static char enzTermSpec = get_enzyme_terminal_specificity ();
				if ( enzTermSpec == 'C' ) {
					endI = pepLen - 1;
				}
				if ( enzTermSpec == 'N' ) {
					startI = 1;
				}
			}
		}
		for ( int i = startI ; i < endI ; i++ ) {
			if ( peptide [i] == expectedAA ) {
				bool uniq = true;
				for ( IntVectorSizeType j = 0 ; j < indexSet.size () ; j++ ) {	// Checks if residue already modified.
					if ( indexSet [j] == i ) {
						uniq = false;
						break;
					}
				}
				if ( uniq ) {
					indexSet.push_back ( i );
					if ( modIndex + 1 < numModResidues ) getNextModification ( modIndex + 1, expectedResidues, modificationResidues, terminalSpecificities, peptide, nTermPeptide, cTermPeptide, sequenceList, z );
					else {
						addSequences ( peptide, modificationResidues, sequenceList );
					}
					indexSet.pop_back ();
				}
			}
		}
	}
}
void ModificationTable::addSequences ( const string& peptide, const CharVector& modificationResidues, PeptideSequenceVector& sequenceList ) const
{
	newPeptide = peptide;
	for ( IntVectorSizeType i = 0 ; i < indexSet.size () ; i++ ) {
		if ( indexSet [i] >= 0 ) newPeptide [indexSet[i]] = modificationResidues [i];
	}
	pair <SetStringIterator, bool> flag = sequenceSet.insert ( newPeptide );
	if ( flag.second ) {
		sequenceList.push_back ( PeptideSequence ( newPeptide, extraMods ) );
		if ( maxSequences && sequenceList.size () > maxSequences ) throw ModificationTableMaximumSequencesExceeded ();
	}
}
void ModificationTable::initialiseHomology ( const StringVector& homologyNames )
{
	if ( modMapInitialised == false ) initialiseModMap ();
	for ( StringVectorSizeType i = 0 ; i < homologyNames.size () ; i++ ) {
		StringVector lines = modMap [homologyNames[i]];
		for ( StringVectorSizeType j = 0 ; j < lines.size () ; j++ ) {
			istringstream ist ( lines [j] );
			char aa;
			ist >> aa;
			string value;
			ist >> value;
			StringConstIterator b = value.begin ();
			StringConstIterator e = value.end ();
			while ( b != e ) {
				char modAA = aa;
				char specif = '0';
				char c = *b;
				if ( c != '(' ) {
					modAA = c;
					b++;
				}
				if ( b != e ) {
					if ( *b == '(' ) {
						specif = *(b+1);
						b += 3;
					}
				}
				modifications.push_back ( new Modification ( aa, modAA, specif ) );
			}
		}
	}
}
void ModificationTable::initialiseUserModHomology ( const StringVector& usermods )
{
	Usermod::initialiseUsermodAAInfo ( usermods );
	dehydroFlag = false;
	cationFlag = false;
	for ( StringVectorSizeType i = 0 ; i < usermods.size () ; i++ ) {
		string u = usermods [i];
		bool isRare = false;
		int mmLimit = 0;
		if ( !u.empty () && u [u.length ()-1] != ')' ) {
			int exIdx = u.find_last_of ( ")" );
			string exStr = u.substr ( exIdx + 4 );
			isRare = ( exStr == "Rare" );
			if ( isPrefix ( exStr, "Max" ) ) {
				mmLimit = atoi ( exStr.substr ( 4 ).c_str () );
			}
			u = u.substr ( 0, exIdx + 1 );
		}
		Usermod umod ( u );
		string aaList = umod.getAAList ();
		string outputString = umod.getOutputString ();
		ModificationSet mm;
		if ( !genStrcasecmp ( outputString.c_str (), "Dehydro" ) ) dehydroFlag = true;
		if ( outputString.find ( "Cation" ) != string::npos ) cationFlag = true;
		char terminalSpecificity = umod.getTerminalSpecificity ();
		for ( StringSizeType j = 0 ; j < aaList.size () ; j++ ) {
			char aa = aaList [j];
			string modAA ( 1, aa );
			modAA += '(' + outputString + ')';
			if ( !genStrcasecmp ( outputString.c_str (), "Phospho" ) ) {
				modifications.push_back ( new Modification ( aa, tolower (aa), terminalSpecificity ) );
			}
			else if ( !genStrcasecmp ( u.c_str (), "Oxidation (M)" ) ) {
				modifications.push_back ( new Modification ( 'M', 'm', terminalSpecificity ) );
			}
			else {
				modifications.push_back ( new Modification ( modAA, terminalSpecificity ) );
			}
			if ( isRare )	rareMods.insert ( modifications.back () );
			if ( mmLimit )	mm.insert ( modifications.back () );
		}
		if ( mmLimit ) {
			maxMods.push_back ( mm );
			maxModsLimit.push_back ( mmLimit );
			maxModsCount.push_back ( 0 );
		}
	}
}
void ModificationTable::initialiseExtraUserModHomology ()
{
	MapPairStringStringExtraUserModInfo mseumi = ExtraUserMods::instance ().getExtraUserMods ();
	for ( MapPairStringStringExtraUserModInfoConstIterator i = mseumi.begin () ; i != mseumi.end () ; i++ ) {
		initialiseUserDefinedLinkHomology ( i );
	}
}
void ModificationTable::initialiseUserDefinedLinkHomology ( const MapPairStringStringExtraUserModInfoConstIterator& iter )
{
	string longName		= (*iter).first.first;
	string specificity	= (*iter).first.second;
	string formula		= (*iter).second.getFormula ();
	string userAMass	= (*iter).second.getUserAMass ();
	string limit		= (*iter).second.getLimit ();
	if ( !genStrcasecmp ( longName.c_str (), "Dehydro" ) ) dehydroFlag = true;
	if ( longName.find ( "Cation" ) != string::npos ) cationFlag = true;
	string aaList;
	char terminalSpecificity;
	bool isRare = ( limit == "Rare" );
	bool isMax = isPrefix ( limit, "Max" );
	ModificationSet mm;
	Usermod::parseUsermodSpecificityLine ( specificity, aaList, terminalSpecificity );

	for ( int i = 0 ; i < aaList.length () ; i++ ) {
		char aa = aaList [i];
		string modAA ( 1, aa );
		modAA += '(' + longName + ')';
		AAInfo::getInfo ().addModAminoAcid ( aa, longName, formula );
		if ( !formula.empty () ) {
			modifications.push_back ( new Modification ( modAA, terminalSpecificity ) );
		}
		else {
			double dOffset = atof ( userAMass.c_str () );
			string sOffsetStr = '(' + gen_ftoa ( dOffset, "%.4f" ) + ')';
			double constModAdjuster;
			if ( aa == 'n' )		constModAdjuster = n_terminus_wt - h1;
			else if ( aa == 'c' )	constModAdjuster = c_terminus_wt - o_h;
			else if ( aa == '.' )	constModAdjuster = 0.0;
			else					constModAdjuster = amino_acid_wt [aa] - massConvert ( AAInfo::getInfo ().getElementalString ( aa ).c_str () );
			modifications.push_back ( new Modification ( aa + sOffsetStr, dOffset + constModAdjuster, terminalSpecificity, false ) );
		}
		if ( isRare )	rareMods.insert ( modifications.back () );
		if ( isMax )	mm.insert ( modifications.back () );
	}
	if ( isMax ) {
		maxMods.push_back ( mm );
		maxModsLimit.push_back ( atoi ( limit.substr ( 4 ).c_str () ) );
		maxModsCount.push_back ( 0 );
	}
}
void ModificationTable::initialiseMassOffsetModifications ( const ModificationParameters& modParams )
{
	StringVector aas = modParams.getCompIon ();
	bool uncleaved = modParams.getUncleaved ();
	bool nTermProtein = modParams.getModNTermProtein ();
	bool cTermProtein = modParams.getModCTermProtein ();
	if ( modParams.getModNeutralLoss () )	aas.push_back ( "." );
	if ( modParams.getModNTerm () )			aas.push_back ( "n" );
	if ( modParams.getModCTerm () )			aas.push_back ( "c" );
	if ( !aas.empty () ) {
		DoubleVector constModAdjuster;
		for ( StringVectorSizeType ii = 0 ; ii < aas.size () ; ii++ ) {
			char aaCode = aas [ii][0];
			if ( aaCode == 'n' )		constModAdjuster.push_back ( n_terminus_wt - h1 );
			else if ( aaCode == 'c' )	constModAdjuster.push_back ( c_terminus_wt - o_h );
			else if ( aaCode == '.' )	constModAdjuster.push_back ( 0.0 );
			else						constModAdjuster.push_back ( amino_acid_wt [aaCode] - massConvert ( AAInfo::getInfo ().getElementalString ( aaCode ).c_str () ) );
		}
		int startNominal = modParams.getStartNominal ();
		int endNominal = modParams.getEndNominal ();
		int modMaxZ = 1;
		if ( modParams.getModRangeType () == "m/z" ) {
			modMaxZ = modParams.getModMaxZ ();
			startNominal *= modMaxZ;
			endNominal *= modMaxZ;
		}
		double defect = modParams.getDefect ();
		double offset = 1.0 + defect;
		for ( int i = startNominal ; i <= endNominal ; i++ ) {
			if ( i != 0 ) {
				double dOffset = i * offset;
				string sOffsetStr = '(' + gen_ftoa ( dOffset, "%.4f" ) + ')';
				dOffset = atof ( sOffsetStr.substr ( 1 ).c_str () );	// Makes sure that the offset is the rounded number
				for ( StringVectorSizeType j = 0 ; j < aas.size () ; j++ ) {
					char ts = '0';
					string a = aas [j];
					if ( uncleaved && ( a != "n" && a != "c" && a != "." ) ) ts = 'e';
					if ( nTermProtein && a == "n" ) ts = 'n';
					if ( cTermProtein && a == "c" ) ts = 'c';
					double mOffset = dOffset - constModAdjuster [j];
					modifications.push_back ( new Modification ( a [0] + sOffsetStr, mOffset, ts, true ) );
					if ( i == startNominal )	mmodStartExact = mOffset / modMaxZ;
					if ( i == endNominal )		mmodEndExact = mOffset / modMaxZ;
					massMods = true;
				}
			}
		}
	}
}
void ModificationTable::initialiseHomologyMultiMod ()
{
//
// maxLevels is the maximum number of modifications on a given peptide
// numLevels is the number of modifications currently under consideration
// mods is the combination of modifications currently under consideration (this is continually overwritten)
// multiModifications is a list of all valid modification combinations
// numMods is a list of all valid modifications
// level is the position in the mods array
//
	numMods = modifications.size ();
	mmod = 0;
	cat = 0;
	nT = 0;
	cT = 0;
	nTermSpec = 0;
	cTermSpec = 0;
	neutralLoss = 0;
	metLoss = 0;
	for ( numLevels = 2 ; numLevels <= maxLevels ; numLevels++ ) {
		mods.resize ( numLevels );
		getNextMod ( 0, 0 );
	}
	sort ( multiModifications.begin (), multiModifications.end (), sort_multi_modifications () );
}
void ModificationTable::getNextMod ( int start, int level )
{
	for ( int i = start ; i < numMods ; i++ ) {
		mods [level] = modifications [i];
		Modification* m = mods [level];
		bool mmFlag = m->getMmod ();
		if ( mmFlag ) mmod++;
		if ( mmod <= maxMMod ) {
			if ( m->getNTerm () ) nT++;
			if ( nT <= 1 ) {
				if ( m->getCTerm () ) cT++;
				if ( cT <= 1 ) {
					char tSpec = m->getTerminalSpecificity ();
					if ( tSpec == 'N' )			nTermSpec++;		// Peptide
					else if ( tSpec == 'n' ) {						// Protein
						if ( m->getMetLoss () ) metLoss++;
						if ( metLoss ) {
							if ( !mmFlag ) {
								nTermSpec++;
							}
						}
						else {
							nTermSpec++;
						}
					}
					if ( nTermSpec <= 1 ) {
						if ( tSpec == 'c' || tSpec == 'C' )	cTermSpec++;
						if ( cTermSpec <= 1 ) {
							if ( m->getNeutralLoss () ) neutralLoss++;
							if ( neutralLoss <= 1 ) {
								if ( cationFlag && m->getCation () )cat++;
								if ( ( !cationFlag || cat <= 1 ) ) {
									if ( checkMaxAndRare ( level ) ) {
										if ( level + 1 < numLevels )
											getNextMod ( i, level + 1 );
										else {
											if ( multiModifications.size () < MAX_MODIFICATIONS ) {
												if ( checkDehydro ( level ) ) {
													multiModifications.push_back ( new MultiModification ( mods ) );
												}
											}
											else 
												ErrorHandler::genError ()->error ( "Too many modification combinations selected.\nFunction: initialiseHomologyMultiMod.\n" );
										}
									}
								}
								if ( cationFlag && m->getCation () ) cat--;
							}
							if ( m->getNeutralLoss () ) neutralLoss--;
						}
						if ( tSpec == 'c' || tSpec == 'C' )	cTermSpec--;
					}
					if ( tSpec == 'N' )			nTermSpec--;		// Peptide
					else if ( tSpec == 'n' ) {						// Protein
						if ( metLoss ) {
							if ( !mmFlag ) {
								nTermSpec--;
							}
						}
						else {
							nTermSpec--;
						}
						if ( m->getMetLoss () ) metLoss--;
					}
				}
				if ( m->getCTerm () ) cT--;
			}
			if ( m->getNTerm () ) nT--;
		}
		if ( mmFlag ) mmod--;
	}
}
bool ModificationTable::checkDehydro ( ModificationVectorSizeType lev )	// There must be an even number of dehydros
{
	if ( !dehydroFlag ) return true;
	int n = 0;
	for ( ModificationVectorSizeType i = 0 ; i <= lev ; i++ ) {
		if ( mods [i]->getDehydro () ) n++;
	}
	return genIsEven ( n );
}
bool ModificationTable::checkMaxAndRare ( ModificationVectorSizeType lev )
{
	if ( rare ) {
		int n = 0;
		for ( ModificationVectorSizeType i = 0 ; i <= lev ; i++ ) {
			if ( rareMods.find ( mods [i] ) != rareMods.end () ) {
				n++;
				if ( n > rare ) return false;
			}
		}
	}
	if ( maxCh ) {
		fill ( maxModsCount.begin (), maxModsCount.end (), 0 );
		for ( ModificationVectorSizeType i = 0 ; i <= lev ; i++ ) {
			for ( VectorModificationSetSizeType j = 0 ; j < maxMods.size () ; j++ ) {
				if ( maxMods [j].find ( mods [i] ) != maxMods [j].end () ) {
					maxModsCount [j]++;
					if ( maxModsCount [j] > maxModsLimit [j] ) return false;
				}
			}
		}
	}
	return true;
}
void ModificationTable::initialiseModMap ()
{
	int numEntries;
	char* info = getFileInfo ( MsparamsDir::instance ().getParamPath ( "homology.txt" ), '>', 1, true, &numEntries );
	for ( int i = 0 ; i < numEntries ; i++ ) {
		string name = ( i == 0 ) ? strtok ( info, "\n" ) : strtok ( NULL, "\n" );
		names.push_back ( name );
		StringVector lines;
		string line;
		while ( ( line = strtok ( NULL, "\n" ) ) != ">" ) {
			lines.push_back ( line );
		}
		modMap [name] = lines;
	}
	delete [] info;
	modMapInitialised = true;
}
StringVector ModificationTable::getLines ( const string& name )
{
	StringVector lines = modMap [name];
	return lines;
}
