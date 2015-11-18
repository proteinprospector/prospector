/******************************************************************************
*                                                                             *
*  Program    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_ambiguity.cpp                                              *
*                                                                             *
*  Created    : February 15th 2011                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2011-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif
#include <algorithm>
#include <cmath>
#include <sstream>
#include <lg_io.h>
#include <lg_string.h>
#include <lu_ambiguity.h>
#include <lu_delim.h>
#include <lu_param_list.h>
#include <lu_table.h>
#include <lu_usermod.h>
using std::set;
using std::string;
using std::vector;
using std::make_pair;
using std::ostringstream;
using std::stable_sort;
using std::copy;
using std::ostream;
using std::back_inserter;
using std::sort;

namespace {
string getMod ( int i )
{
	if ( i == -3 )		return "N-term";
	else if ( i == -2 )	return "C-term";
	else if ( i == -1 )	return "Neutral loss";
	else				return gen_itoa ( i );
}
}

class sortModByResidue {
public:
	int operator () ( const PairStringInt& a, const PairStringInt& b ) const
	{
		return a.second < b.second;
	}
};
ModificationAmbiguity::ModificationAmbiguity ( const string& mods, bool chopScores ) :
	cur ( 0 )
{
	string ambiguity;
	string unambiguity;
	SCModInfo::getAmbiguityString ( mods, ambiguity, unambiguity, chopScores );		// Gets an ambiguity string for passing to MS-Product
	if ( !ambiguity.empty () ) initAmbiguity ( ambiguity );
	if ( !unambiguity.empty () ) initUnambiguity ( unambiguity, chopScores );
}
void ModificationAmbiguity::initUnambiguity ( const string& unambiguity, bool chopScores )
{
	string::size_type idx1 = 0;
	for ( ; ; ) {
		string::size_type idx2 = unambiguity.find ( ';', idx1 );		// Find the first ';'
		if ( idx2 == string::npos ) break;
		string ms1 = unambiguity.substr ( idx1, idx2-idx1 );
		string::size_type idx3 = ms1.find ( '@' );
		string mod = ms1.substr ( 0, idx3 );
		string resStr = ms1.substr ( idx3+1 );
		if ( !chopScores ) {
			string::size_type idx4 = resStr.find ( '=' );
			if ( idx4 != string::npos ) {
				slip.push_back ( atoi ( resStr.substr ( idx4+1 ).c_str () ) );
				resStr = resStr.substr ( 0, idx4 );		// Chop off the score element
			}
			else {
				slip.push_back ( -1 );
			}
		}
		int res;
		if ( resStr == "N-term" )			res = -3;
		else if ( resStr == "C-term" )		res = -2;
		else if ( resStr == "Neutral loss" )res = -1;
		else								res = atoi ( resStr.c_str () );
		unambigMods.push_back ( make_pair ( mod, res ) );
		idx1 = idx2 + 1;
	}
}
void ModificationAmbiguity::initAmbiguity ( const string& ambiguity )
{
// a. Phospho@1|2|3
// b. Phospho@11&15|11&17|15&17
// c. Oxidation@6&Phospho@1|Oxidation@7&Phospho@1|Oxidation@6&Phospho@2|Oxidation@7&Phospho@2
// d. Oxidation@6|7;Phospho@1|2
	string::size_type idx1 = 0;
	for ( ; ; ) {
		string::size_type idx2 = ambiguity.find ( ';', idx1 );		// Find the first ';'
		if ( idx2 == string::npos ) break;
		string ms1 = ambiguity.substr ( idx1, idx2-idx1 );
		VectorVectorPairStringInt vvpsi;
		string::size_type idx3 = 0;
		for ( ; ; ) {
			string::size_type idx4 = ms1.find ( "||", idx3 );					// Find the first "||"
			string ms2;
			if ( idx4 == string::npos ) ms2 = ms1.substr ( idx3 );
			else						ms2 = ms1.substr ( idx3, idx4-idx3 );
			if ( ms2.find ( "@(" ) == string::npos )
				parseSimple ( vvpsi, ms2 );
			else
				parseComplex ( vvpsi, ms2 );
			if ( idx4 == string::npos ) break;
			idx3 = idx4 + 2;
		}
		vvvpsi.push_back ( vvpsi );
		idx1 = idx2 + 1;
	}
	getNextMod ( 0 );
}
int ModificationAmbiguity::getResidueIdx ( const string& resStr )
{
	if ( resStr == "N-term" )			return -3;
	else if ( resStr == "C-term" )		return -2;
	else if ( resStr == "Neutral loss" )return -1;
	else								return atoi ( resStr.c_str () );
}
void ModificationAmbiguity::parseSimple ( VectorVectorPairStringInt& vvpsi, const string& ms1 )
{
	string::size_type idx3 = 0;
	string mod;
	for ( ; ; ) {
		string::size_type idx4 = ms1.find ( '|', idx3 );
		bool last1 = ( idx4 == string::npos );
		string ms2 = last1 ? ms1.substr ( idx3 ) : ms1.substr ( idx3, idx4-idx3 );
// a. Phospho@1
// b. Phospho@11&15
// c. Oxidation@6&Phospho@1
		string::size_type idx5 = 0;
		vector <PairStringInt> vpsi;
		for ( ; ; ) {
			string::size_type idx6 = ms2.find ( '&', idx5 );
			bool last2 = ( idx6 == string::npos );
			string ms3 = last2 ? ms2.substr ( idx5 ) : ms2.substr ( idx5, idx6-idx5 );
			string::size_type idx7 = 0;
			string::size_type idx8 = ms3.find ( '@', idx7 );
			string resStr;
			if ( idx8 == string::npos ) {
				resStr = ms3.c_str();
			}
			else {
				mod = ms3.substr ( idx7, idx8-idx7 );
				resStr = ms3.substr ( idx8+1 ).c_str ();
			}
			vpsi.push_back ( make_pair ( mod, getResidueIdx ( resStr ) ) );
			if ( last2 ) break;
			idx5 = idx6 + 1;
		}
		vvpsi.push_back ( vpsi );
		if ( last1 ) break;
		idx3 = idx4 + 1;
	}
}
// Eg
// HexNAc&Phospho&Phospho@(17&13&(20|21|23|24|25|27|28))|(20&(9&24)|(13&(17|21|23|24|25|27|28)))|(21&13&(17|20|23|24|25|27|28))|(23&13&(17|20|21|24|25|27|28))|(25&13&(17|20|21|23|24|27|28))|(27&13&(17|20|21|23|24|25|28))|(28&13&(17|20|21|23|24|25|27))
void ModificationAmbiguity::parseComplex ( VectorVectorPairStringInt& vvpsi, const string& ms1 )
{
	compPairs.clear ();
	string::size_type idx3 = ms1.find ( "@(" );
	string mods = ms1.substr ( 0, idx3 );			// Extract HexNAc&Phospho&Phospho
	string::size_type start = 0;
	string::size_type end;
	for ( ; ; ) {
		string m = genNextString ( mods, "&", start, end );			// m is the mod ie HexNAC or Phospho
		if ( end != string::npos )
			compPairs.push_back ( make_pair ( m, 0 ) );
		else {
			compPairs.push_back ( make_pair ( mods.substr ( start ), 0 ) );
			break;
		}
	}
	complexVar = ms1.substr ( idx3+1 );								// This is everything after the @ symbol
	vIdx = 0;
	parseNextComplex ( 0, vvpsi );
}
// A compPair is a vector of a pair of a modification and a position. It thus represents a single modification state out of multiple possibilities.
// vvpsi then represents a list of these modification state
// vvvpsi is then a list of vvpsi
void ModificationAmbiguity::parseNextComplex ( int start, VectorVectorPairStringInt& vvpsi )
{
	string s;
	for ( int i = start ; ; ) {
		char ch = complexVar [vIdx++];
		if ( vIdx > complexVar.length () ) break;
		if ( ch == '(' ) parseNextComplex ( i, vvpsi );
		else if ( ch == '|' ) {
			compPairs [i].second = getResidueIdx ( s );
			if ( i == compPairs.size () - 1 ) {
				vvpsi.push_back ( compPairs );
			}
			s = "";
		}
		else if ( ch == '&' ) {
			compPairs [i].second = getResidueIdx ( s );
			s = "";
			i++;
		}
		else if ( ch == ')' ) {
			if ( !s.empty () ) {
				compPairs [i].second = getResidueIdx ( s );
				if ( i == compPairs.size () - 1 ) {
					vvpsi.push_back ( compPairs );
				}
				s = "";
			}
			break;
		}
		else s += ch;
	}
}
void ModificationAmbiguity::getNextMod ( int level )
{
	VectorVectorPairStringInt& vvpsi = vvvpsi [level];
	for ( int i = 0 ; i < vvpsi.size () ; i++ ) {
		vector <PairStringInt>& vpsi = vvpsi [i];
		for ( int j = 0 ; j < vpsi.size () ; j++ ) {
			curMod.push_back ( vpsi [j] );
		}
		if ( level+1 < vvvpsi.size () ) getNextMod ( level+1 );
		else {
			modArray.push_back ( curMod );
			stable_sort ( modArray.back ().begin (), modArray.back ().end (), sortModByResidue () );
		}
		for ( int k = 0 ; k < vpsi.size () ; k++ ) {
			curMod.pop_back ();
		}
	}
}
bool ModificationAmbiguity::getNextSequence ( string& s, string& nTerm, string& cTerm, string& nLoss ) const
{
	bool amb = ( cur < modArray.size () );
	vector <PairStringInt> vpsi = unambigMods;
	if ( amb ) copy ( modArray [cur].begin (), modArray [cur].end (), back_inserter ( vpsi ) );		// Append the ambiguous mods on the end of the unambiguous mods
	stable_sort ( vpsi.begin (), vpsi.end (), sortModByResidue () );
	vector <PairStringInt> vpsi2;
	for ( int iii = 0 ; iii < vpsi.size () ; iii++ ) {
		if ( iii != 0 && vpsi [iii].second == vpsi [iii-1].second ) {
			vpsi2.pop_back ();
			vpsi [iii].first = vpsi [iii-1].first + "+" + vpsi [iii].first;
		}
		vpsi2.push_back ( vpsi [iii] ); 
	}
	string t;
	int num = 0;
	string mod = vpsi2 [num].first;
	int residue = vpsi2 [num].second;
	if ( residue == -3 ) {
		nTerm = mod;
		if ( ++num >= vpsi2.size () ) {
			if ( amb ) cur++;
			return cur < modArray.size ();
		}
		mod = vpsi2 [num].first;
		residue = vpsi2 [num].second;
	}
	if ( residue == -2 ) {
		cTerm = mod;
		if ( ++num >= vpsi2.size () ) {
			if ( amb ) cur++;
			return cur < modArray.size ();
		}
		mod = vpsi2 [num].first;
		residue = vpsi2 [num].second;
	}
	if ( residue == -1 ) {
		nLoss = mod;
		if ( ++num >= vpsi2.size () ) {
			if ( amb ) cur++;
			return cur < modArray.size ();
		}
		mod = vpsi2 [num].first;
		residue = vpsi2 [num].second;
	}
	if ( residue == 0 ) {					// Ignore cases where residue == 0, this means one of the variable options is no mod
		if ( ++num >= vpsi2.size () ) {
			if ( amb ) cur++;
			return cur < modArray.size ();
		}
		mod = vpsi2 [num].first;
		residue = vpsi2 [num].second;
	}
	for ( StringSizeType i = 0, idx = 1 ; i < s.length () ; i++ ) {
		bool replace = ( residue == idx );
		char aa = s [i];
		t += aa;
		if ( replace ) {
			t += '(' + mod + ')';
			if ( ++num < vpsi2.size () ) {
				mod = vpsi2 [num].first;
				residue = vpsi2 [num].second;
			}
		}
		idx++;
	}
	s = t;
	if ( amb ) cur++;
	return cur < modArray.size ();
}
void ModificationAmbiguity::getUnambiguousIndexList ( const string& mod, IntVector& sites, IntVector& scores ) const
{
	for ( VectorPairStringIntSizeType i = 0 ; i < unambigMods.size () ; i++ ) {
		if ( unambigMods [i].first == mod ) {
			sites.push_back ( unambigMods [i].second );
			scores.push_back ( slip [i] );
		}
	}
	for ( VectorVectorVectorPairStringIntSizeType j = 0 ; j < vvvpsi.size () ; j++ ) {
		for ( VectorVectorPairStringIntSizeType k = 0 ; k < vvvpsi [j].size () ; k++ ) {
			for ( VectorPairStringIntSizeType l = 0 ; l < vvvpsi [j][k].size () ; l++ ) {
				if ( vvvpsi [j][k][l].first == mod ) {
					sites.push_back ( vvvpsi [j][k][l].second );
					scores.push_back ( 0 );
				}
			}
		}
	}
}

//search # <map <pair <acc#, pair <mod, residue> >, pair <SpecID*, score>>>
//vector <std::map <pair <string, PairStringInt >, pair <const SpecID*, double> > > SCModInfo::bSS;
StringVectorVector SCModInfo::cMods;
StringVector SCModInfo::cNTerm;
StringVector SCModInfo::cCTerm;
double SCModInfo::threshold = 0.0;
vector <MapSpecIDToScoreString> SCModInfo::msiss;
void SCModInfo::init ( int searchNumber, const StringVector& constMods, double t )
{
	threshold = t;
	if ( msiss.size () < searchNumber ) msiss.resize ( searchNumber );

	cNTerm.resize ( searchNumber );
	cCTerm.resize ( searchNumber );
	cMods.resize ( searchNumber );
	int snInx = searchNumber - 1;
	for ( StringVectorSizeType i = 0 ; i < constMods.size () ; i++ ) {
		const string& cm = constMods [i];
		size_t idx1 = cm.rfind ( '(' ) + 1;
		size_t idx2 = cm.length () - 1;
		string aa = cm.substr ( idx1, idx2-idx1 );
		string mod = cm.substr ( 0, idx1-2 );
		if ( aa == "N-term" )		cNTerm [snInx] = mod;
		else if ( aa == "C-term" )	cCTerm [snInx] = mod;
		else {
			for ( int j = 0 ; j < aa.size () ; j++ ) {
				cMods [snInx].push_back ( string ( 1, aa [j] ) + "(" + mod + ")" );
			}
		}
	}
}
string SCModInfo::removeCMods ( const string& peptide )
{
	string pep2 = peptide;
	for ( int i = 0 ; i < cMods.back ().size () ; i++ ) {
		const string& cm = cMods.back () [i];
		pep2 = genReplaceSubstrings ( pep2, cm, string ( 1, cm [0] ) );
	}
	return pep2;
}
void SCModInfo::addHit ( const string& nTerm, const string& peptide, const string& cTerm, const string& neutralLoss, double score, double expectation, int numSpectra )
{	// map < pair <Residue Index, Modification>, score >
	string nT = ( nTerm == cNTerm.back () ) ? "" : nTerm;
	string cT = ( cTerm == cCTerm.back () ) ? "" : cTerm;
	string pep2 = removeCMods ( peptide );
	string dbPeptide = 	gen_strstriptags2 ( pep2, '(', ')' );
	double ps;
	if ( expectation == -1.0 ) {
		ps = score;
		eval = false;
	}
	else {
		ps = -10.0 * log10 ( expectation / numSpectra );
		eval = true;
	}
	bool flag = false;
	for ( int i = 0 ; i < topDBPeptide.size () ; i++ ) {
		if ( dbPeptide == topDBPeptide [i] ) {
			MapPairIntStringToDouble mods = getMods ( nT, pep2, cT, neutralLoss );
			MapPairIntStringToDouble& tm = topMods [i];
			double diff = topDBScore [i] - ps;
			if ( diff <= threshold ) {
				zeroMods [i].insertIfEmpty ( topMods [i] );
				zeroMods [i].insert ( mods );
			}
			for ( MapPairIntStringToDoubleIterator j = tm.begin () ; j != tm.end () ; j++ ) {
				if ( mods.find ( (*j).first ) == mods.end () ) {		// If one of top mods not in current mod
					if ( (*j).second == -1.0 ) {	// Set that top mod's diff score
						(*j).second = diff;
					}
				}
			}
			flag = true;
			break;
		}
	}
	if ( !flag ) {
		topDBPeptide.push_back ( dbPeptide );
		topDBScore.push_back ( ps );
		topMods.push_back ( getMods ( nT, pep2, cT, neutralLoss ) );
		zeroMods.push_back ( ZeroMod () );
	}
}
void SCModInfo::compress ()
{
	for ( int i = 0 ; i < topDBPeptide.size () ; i++ ) {
		const MapPairIntStringToDouble& tm = topMods [i];
		for ( MapPairIntStringToDoubleConstIterator j = tm.begin () ; j != tm.end () ; j++ ) {
			if ( (*j).second > threshold || (*j).second == -1.0 ) {
				zeroMods [i].deleteMod ( (*j).first );
			}
		}
	}
}
MapStringToString SCModInfo::getScoreString ()
{
	MapStringToString mss;
	for ( int i = 0 ; i < topDBPeptide.size () ; i++ ) {
		const MapPairIntStringToDouble& tm = topMods [i];	// map of pair <residue, mod> to score

		ostringstream ostr;
		for ( MapPairIntStringToDoubleConstIterator j = tm.begin () ; j != tm.end () ; j++ ) {
			double score = (*j).second;
			if ( score > threshold || score == -1.0 ) {
				ostr << (*j).first.second << "@";
				ostr << getMod ( (*j).first.first );
				if ( score != -1.0 ) {
					ostr << "=";
					genPrint ( ostr, score, eval ? 0 : 1 );
				}
				ostr << ';';
			}
		}
		zeroMods [i].print ( ostr );
		mss [topDBPeptide [i]] = genStrtrimSemiColon ( ostr.str () );
	}
	return mss;
}
MapPairIntStringToDouble SCModInfo::getMods ( const string& nTerm, const string& peptide, const string& cTerm, const string& neutralLoss )
{
	MapPairIntStringToDouble mcs;
	if ( !nTerm.empty () ) mcs [make_pair ( -3, nTerm )] = -1.0;
	for ( StringSizeType i = 0, idx = 0 ; i < peptide.length () ; i++ ) {
		if ( peptide [i] == '(' ) {
			int bracket = 0;
			StringSizeType start = i+1;
			for ( ; i < peptide.length () ; i++ ) {
				char a = peptide [i];
				if ( a == '(' ) bracket++;
				if ( a == ')' ) bracket--;
				if ( bracket == 0 ) break;
			}
			mcs [make_pair ( idx, peptide.substr ( start, i-start ))] = -1.0;
		}
		else idx++;
	}
	if ( !cTerm.empty () ) mcs [make_pair ( -2, cTerm )] = -1.0;
	if ( !neutralLoss.empty () ) mcs [make_pair ( -1, neutralLoss)] = -1.0;
	return mcs;
}
void SCModInfo::add ( const SpecID* spID )
{
	compress ();
	msiss.back ()[spID] = getScoreString ();
}
string SCModInfo::addStartAAToModString ( const string& modString, int startAA, int len )
{
	string s;
	string::size_type idx1 = 0;
	startAA -= 1;
	for ( ; ; ) {
		string::size_type idx2 = modString.find ( ';', idx1 );
		bool last = ( idx2 == string::npos );
		string ms1;
		if ( last ) ms1 = modString.substr ( idx1 );
		else		ms1 = modString.substr ( idx1, idx2-idx1 );
		string::size_type idx3 = ms1.find ( '@' );
		string mod = ms1.substr ( 0, idx3+1 );
		string ms2 = ms1.substr ( idx3+1 );
		if ( ms2 [0] == '(' ) {				// new compact style
			s += mod;
			string::size_type len = ms2.length ();
			for ( int i = 0 ; i < len ; ) {
				char c = ms2 [i];
				if ( isdigit ( c ) || c == 'N' || c == 'C' ) {	// residue, C-term, N-term, Neutral loss
					string::size_type idx4 = ms2.find_first_of ( "|&)", i );
					string::size_type len2 = idx4-i;
					string ms3 = ms2.substr ( i, len2 );
					int num = atoi ( ms3.c_str () );
					bool nLossFlag = false;
					if ( num == 0 ) {
						if ( ms3 == "N-term" ) {
							num = 1;
						}
						else if ( ms3 == "C-term" ) {
							num = len;
						}
						else if ( ms3 == "Neutral loss" ) {
							nLossFlag = true;
						}
					}
					if ( nLossFlag ) {
						s += gen_itoa ( startAA+1 ) + "-" + gen_itoa ( startAA+len );
					}
					else {
						num += startAA;						// residue number set here
						s += gen_itoa ( num );
					}
					i += len2;
				}
				else {
					s += c;
					i++;
				}
			}
			if ( last ) break;
			s += ';';
		}
		else {								// old style
			string::size_type idx4 = ms2.find ( '=' );	// Defines if there is a score
			string res = ( idx4 == string::npos ) ? ms2 : ms2.substr ( 0, idx4 );
			string ms3 = ( idx4 == string::npos ) ? "" : ms2.substr ( idx4 );
			string::size_type idx5 = 0;
			bool ok = true;
			string output = mod;
			for ( ; ; ) {
				string::size_type idx6 = res.find_first_of ( "|&", idx5 );
				bool last2 = ( idx6 == string::npos );
				string ms4 = last2 ? res.substr ( idx5 ) : res.substr ( idx5, idx6-idx5 );
				string::size_type idx7 = ms4.find ( '@' );
				if ( idx7 != string::npos ) {
					mod = ms4.substr ( 0, idx7+1 );
					output += mod;
					ms4 = ms4.substr ( idx7+1 );
				}
				int num = atoi ( ms4.c_str () );
				bool nLossFlag = false;
				if ( num == 0 ) {
					if ( ms4 == "N-term" ) {
						num = 1;
					}
					else if ( ms4 == "C-term" ) {
						num = len;
					}
					else if ( ms4 == "Neutral loss" ) {
						nLossFlag = true;
					}
					else {
						ok = false;
						break;
					}
				}
				if ( nLossFlag ) {
					output += gen_itoa ( startAA+1 ) + "-" + gen_itoa ( startAA+len );
				}
				else {
					num += startAA;						// residue number set here
					output += gen_itoa ( num );
				}
				if ( last2 ) break;
				else output += res [idx6];

				idx5 = idx6 + 1;
			}
			output += ms3;
			if ( ok ) s += output;
			if ( last ) break;
			if ( ok ) s += ';';
		}
		idx1 = idx2 + 1;
	}
	return s;
}
string SCModInfo::getConstModsString ( const string& dbPeptide, int searchIdx )
{
	string cModString;
	if ( !cNTerm [searchIdx].empty () ) {
		if ( !cModString.empty () ) cModString += ';';
		cModString += cNTerm [searchIdx] + "@N-term";
	}
	if ( !cCTerm [searchIdx].empty () ) {
		if ( !cModString.empty () ) cModString += ';';
		cModString += cCTerm [searchIdx] + "@C-term";
	}
	for ( int i = 0 ; i < cMods [searchIdx].size () ; i++ ) {
		const string& cm = cMods [searchIdx][i];
		int startInd = 0;
		for ( ; ; ) {
			int ind = dbPeptide.find ( cm [0], startInd );
			if ( ind == string::npos ) break;
			else {
				if ( !cModString.empty () ) cModString += ';';
				cModString += cm.substr ( 2, cm.length () - 3 ) + "@" + gen_itoa ( ind+1 );
				startInd = ind + 1;
			}
		}
		//pep2 = genReplaceSubstrings ( pep2, cm, string ( 1, cm [0] ) );
	}
	return cModString;
}
string SCModInfo::getAllModsString ( int searchNumber, const SpecID* spID, const string& dbPeptide, int startAA )
{
	string allMods = getModsString ( searchNumber, spID, dbPeptide, startAA );
	string constMods = getConstModsString ( dbPeptide, searchNumber );
	if ( allMods.empty () ) return constMods;
	else {
		if ( constMods.empty () )	return allMods;
		else						return constMods + ';' + allMods;
	}
}
string SCModInfo::getModsString ( int searchNumber, const SpecID* spID, const string& dbPeptide, int startAA )
{
	MapSpecIDToScoreStringConstIterator cur = msiss [searchNumber].find ( spID );
	if ( cur != msiss [searchNumber].end () ) {
		const MapStringToString& mss = (*cur).second;
		MapStringToStringConstIterator cur2 = mss.find ( dbPeptide );
		if ( cur2 != mss.end () )
			return startAA == 0 ? (*cur2).second : addStartAAToModString ( (*cur2).second, startAA, dbPeptide.length () );
		else
			return "";
	}
	else
		return "";
}
//string SCModInfo::getAmbiguityString ( int searchNumber, const SpecID* spID, const string& dbPeptide )
//{
//	return getAmbiguityString ( getModsString ( searchNumber, spID, dbPeptide ) );
//}
void SCModInfo::getAmbiguityString ( const string& modsString, string& ambString, string& unambString, bool chopScores )		// Gets an ambiguity string for passing to MS-Product
{
	ambString = "";
	unambString = "";
	string::size_type idx1 = 0;
	bool end = false;
	for ( ; ; ) {
		string::size_type idx2 = modsString.find ( ';', idx1 );
		if ( idx2 == string::npos ) end = true;
		string as;
		bool amb = false;
		if ( end )	as = modsString.substr ( idx1 );				// Get each individual mod
		else		as = modsString.substr ( idx1, idx2-idx1 );
		if ( as.find ( '|' ) != string::npos || as.find ( '&' ) != string::npos ) amb = true;			// Is the mod ambiguous
		if ( chopScores ) {
			string::size_type idx3 = as.find ( '=' );
			if ( idx3 != string::npos ) as = as.substr ( 0, idx3 );		// Chop off the score element
		}
		if ( amb )	ambString += as + ';';
		else		unambString += as + ';';
		if ( end ) break;
		idx1 = idx2+1;
		while ( isspace ( modsString [idx1] ) ) idx1++;
	}
}
void ZeroMod::insertIfEmpty ( const MapPairIntStringToDouble& tm )
{
	if ( zm.empty () ) insert ( tm );
}
void ZeroMod::insert ( const MapPairIntStringToDouble& tm )
{
	SetPairIntString spis;
	for ( MapPairIntStringToDoubleConstIterator i = tm.begin () ; i != tm.end () ; i++ ) {
		spis.insert ( (*i).first );
	}
	zm.insert ( spis );
}
void ZeroMod::deleteMod ( const PairIntString& pss )
{
	for ( SetSetPairIntStringIterator i = zm.begin () ; i != zm.end () ; i++ ) {
		SetPairIntString& m = const_cast<SetPairIntString&>((*i));
		m.erase ( pss );
	}
}
void ZeroMod::print2 ( ostream& os )
{
	int x3 = 0;
	bool first = true;
	SetString modSet;
	for ( SetSetPairIntStringConstIterator i = zm.begin () ; i != zm.end () ; i++, x3++ ) {		// This is all the ambiguous mods for a given peptide
		MapStringToString mso;
		const SetPairIntString& zm2 = (*i);
		for ( SetPairIntStringConstIterator j = zm2.begin () ; j != zm2.end () ; j++ ) {
			string mod = (*j).second;
			int res = (*j).first;
			if ( mso.find ( mod ) != mso.end () ) {
				mso [mod] += "&" + getMod ( res );
			}
			else mso [mod] = getMod ( res );
		}
		int n3 = 0;
		for ( MapStringToStringConstIterator k = mso.begin () ; k != mso.end () ; k++, n3++ ) {
			modSet.insert ( (*k).first );
			if ( modSet.size () > 1 || first ) os << (*k).first + "@";
			first = false;
			os << (*k).second;
			if ( n3 != mso.size () - 1 ) os << '&';
		}
		if ( x3 != zm.size () - 1 ) os << '|';
	}
}
typedef std::map <std::string, IntVector> MapStringToIntVector;
typedef MapStringToIntVector::iterator MapStringToIntVectorIterator;
typedef MapStringToIntVector::const_iterator MapStringToIntVectorConstIterator;

class sortIVVByIndex {
public:
	bool operator () ( const IntVector& lhs, const IntVector& rhs ) const
	{
		int siz = lhs.size ();
		for ( int i = 0 ; i < siz ; i++ ) {
			if ( lhs [i] < rhs [i] ) return true;
			if ( lhs [i] > rhs [i] ) break;
		}
		return false;
	}
};

void ZeroMod::print ( ostream& os )
{
	std::map <StringVector, IntVectorVector> zz;
	for ( SetSetPairIntStringConstIterator i = zm.begin () ; i != zm.end () ; i++ ) {
		MapStringToIntVector mso;
		const SetPairIntString& zm2 = (*i);
		for ( SetPairIntStringConstIterator j = zm2.begin () ; j != zm2.end () ; j++ ) {
			string mod = (*j).second;
			int res = (*j).first;
			if ( mso.find ( mod ) != mso.end () ) {
				mso [mod].push_back ( res );
			}
			else mso [mod].push_back ( res );
		}
		StringVector sv;
		IntVector iv;
		for ( MapStringToIntVectorConstIterator k = mso.begin () ; k != mso.end () ; k++ ) {
			string s = (*k).first;
			IntVector iv2 = (*k).second;
			for ( int m = 0 ; m < iv2.size () ; m++ ) {
				sv.push_back ( s );
				iv.push_back ( iv2 [m] );
			}
		}
		zz [sv].push_back ( iv );
	}
	int num = 0;
	for ( std::map <StringVector, IntVectorVector>::iterator n = zz.begin () ; n != zz.end () ; n++ ) {
		if ( num > 0 ) os << "||";
		const StringVector& sv = (*n).first;
		ivv = (*n).second;
		sort ( ivv.begin (), ivv.end (), sortIVVByIndex () ); 
		for ( int x = 0 ; x < sv.size () ; x++ ) {
			os << sv [x];
			if ( x != sv.size () - 1 ) os << '&';
		}
		siz1 = ivv.size ();
		siz2 = ivv [0].size ();
		if ( siz2 ) {
			os << '@';
			writeNext ( os, 0, siz1, 0 );
		}
		num++;
	}
}
void ZeroMod::writeNext ( ostream& os, int start, int end, int level )
{
	int num = 1;
	bool bracket = false;
	bool topLevel = ( level == siz2 - 1 );
	if ( ( topLevel && ( end - start > 1 ) && siz2 > 1 ) || ( level == 0 && siz2 > 1 ) ) os << '(';
	for ( int i = start ; i < end ; i++ ) {
		bool lastI = ( i == end - 1 );
		int cur = ivv [i][level];
		int next = -100;
		if ( !lastI ) next = ivv [i+1][level];
		if ( cur != next ) {
			if ( !lastI ) bracket = true;
			if ( bracket && !topLevel ) os << '(';
			os << getMod ( cur );
			if ( level+1 < siz2 ) {
				os << '&';
				writeNext ( os, start, start+num, level+1 );
				start += num;
				num = 0;
			}
			if ( bracket && !topLevel ) os << ')';
			if ( !lastI ) os << '|';
		}
		num++;
	}
	if ( ( topLevel && ( end - start > 1 ) && siz2 > 1 ) || ( level == 0 && siz2 > 1 ) ) os << ')';
}
ProbabilityAmbiguity::ProbabilityAmbiguity ( const string& probStr, const string& mod, double probLimit )
{
	if ( !probStr.empty () ) {
		double totalProb = 0;
		for ( StringSizeType i = 0, idx = 0 ; i < probStr.length () ; i++ ) {
			if ( probStr [i] == '(' ) {
				i++;
				StringSizeType start = i;
				for ( ; i < probStr.length () ; i++ ) {
					if ( probStr [i] == ')' ) break;
				}
				StringSizeType end = i;
				double prob = atof ( probStr.substr ( start, end-start ).c_str () );
				totalProb += prob;
				if ( prob > probLimit ) {
					siteList.push_back ( idx );
				}
			}
			else idx++;
		}
		numSites = floor ( totalProb + 0.5 );
		siteListLength = siteList.size ();
		getNext ( 0 );
		modStr = mod;
		modStr += '@';
		for ( IntVectorVectorSizeType j = 0 ; j < hits.size () ; j++ ) {
			IntVector& iv = hits [j];
			for ( IntVectorVectorSizeType k = 0 ; k < iv.size () ; k++ ) {
				modStr += gen_itoa ( iv [k] );
				if ( k < iv.size () - 1 ) modStr += '&';
			}
			if ( j < hits.size () - 1 ) modStr += '|';
		}
	}
}
void ProbabilityAmbiguity::getNext ( int level )
{
	for ( int i = level ; i < siteListLength ; i++ ) {
		sequence.push_back ( siteList [i] );
		if ( sequence.size () < numSites ) {
			level += 1;
			getNext ( level );
		}
		else hits.push_back ( sequence );
		sequence.pop_back ();
	}
}
bool SiteInfo::getIndex ( int& idx ) const
{
	if ( index == -2 ) return false;
	else {
		idx = index;
		return true;
	}
}
void SiteInfo::printHTMLHeaderSLIP ( ostream& os, int idx, int rowspan )
{
	string header = "SLIP";
	if ( idx != -1 ) header += " " + gen_itoa ( idx+1 );
	tableHeader ( os, header, "", "", true, 0, rowspan );
}
void SiteInfo::printHTMLSLIP ( ostream& os ) const
{
	if ( slip == -2 ) {
		tableEmptyCell ( os );						// no match
	}
	else {
		if ( slip == -1 )	tableCell ( os, "-" );	// no ambiguity
		else				tableCell ( os, slip );	// slip
	}
}
bool SiteInfo::printHTMLAA ( ostream& os ) const
{
	if ( index == -2 ) return false;
	else {
		tableCell ( os, aa, true );
		return true;
	}
}
void SiteInfo::printDelimitedHeaderSLIP ( ostream& os, int idx )
{
	string header = "SLIP";
	if ( idx != -1 ) header += " " + gen_itoa ( idx+1 );
	delimitedHeader ( os, header );
}
void SiteInfo::printDelimitedSLIP ( ostream& os ) const
{
	if ( slip == -2 ) {
		delimitedEmptyCell ( os );
	}
	else {
		if ( slip == -1 )	delimitedCell ( os, "-" );	// no ambiguity
		else				delimitedCell ( os, slip );	// slip
	}
}
bool SiteInfo::printDelimitedAA ( ostream& os ) const
{
	if ( index == -2 ) return false;
	else {
		delimitedCell ( os, aa );
		return true;
	}
}

class SiteInfoAscending {
public:
	bool operator () ( const SiteInfo& a, const SiteInfo& b ) const
	{
		return a.site < b.site;
	}
};

class GlycoSiteInfoAscending {
public:
	bool operator () ( const GlycoSiteInfo& a, const GlycoSiteInfo& b ) const
	{
		if ( a.site == b.site )
			return a.mod < b.mod;
		else
			return a.site < b.site;
	}
};

SetString SiteScores::glycoMods;
MapStringToMapCharToInt SiteScores::msmci;
int SiteScores::idx = 0;
bool SiteScores::initialised = false;

SiteScores::SiteScores ()
{
	siteScores.resize ( idx );
}
void SiteScores::initMod ( const string& mod, char aa )
{
	MapStringToMapCharToIntConstIterator cur1 = msmci.find ( mod );
	if ( cur1 != msmci.end () ) {									// mod exists
		const MapCharToInt& mci = (*cur1).second;
		MapCharToIntConstIterator cur2 = mci.find ( aa );
		if ( cur2 != mci.end () ) {									// aa exists
			return;
		}
	}
	msmci [mod][aa] = idx++;
}
void SiteScores::init ( const ParameterList* pList )
{
	StringVector mods = pList->getStringVectorValue ( "msms_mod_AA" );
	init ( mods );
}
void SiteScores::init ( const StringVector& mods, bool flag )
{
	if ( mods.empty () ) return;
	for ( StringVectorSizeType i = 0 ; i < mods.size () ; i++ ) {
		string u = mods [i];
		int pos = u.rfind ( '(' );
		string aa = u.substr ( pos+1, u.length()-pos-2 );
		string outputString = u.substr ( 0, pos-1 );
		if ( aa != "N-term" && aa != "C-term" && aa != "Neutral loss" ) {
			initMod ( outputString, aa [0] );
		}
		if ( Usermod::isGlyco ( outputString ) ) glycoMods.insert ( outputString );
	}
	initialised = true;
}
void SiteScores::init ( const StringVector& mods )
{
	if ( initialised == true ) return;
	for ( StringVectorSizeType i = 0 ; i < mods.size () ; i++ ) {
		string u = mods [i];
		if ( !u.empty () && u [u.length ()-1] != ')' ) {
			int exIdx = u.find_last_of ( ")" );
			u = u.substr ( 0, exIdx + 1 );
		}
		Usermod umod ( u );
		string aaList = umod.getAAList ();
		string outputString = umod.getOutputString ();
		for ( int j = 0 ; j < aaList.length () ; j++ ) {
			char c = aaList [j];
			if ( isupper ( c ) ) initMod ( outputString, c );
		}
		if ( Usermod::isGlyco ( outputString ) ) glycoMods.insert ( outputString );
	}
}
void SiteScores::add ( const string& peptide, const string& mods, int start, int index )
{
	ModificationAmbiguity ma ( mods, false );
	for ( MapStringToMapCharToIntConstIterator i = msmci.begin () ; i != msmci.end () ; i++ ) {
		IntVector sites, scores;
		ma.getUnambiguousIndexList ( (*i).first, sites, scores );
		const MapCharToInt& mci = (*i).second;
		for ( int j = 0 ; j < scores.size () ; j++ ) {
			int site = sites [j];
			char aa = peptide [site-start];
			MapCharToIntConstIterator cur2 = mci.find ( aa );
			if ( cur2 != mci.end () ) {
				int ssIdx = (*cur2).second;
				int score = scores [j];
				MapIntToPairIntIntIterator cur = siteScores [ssIdx].find ( site );
				if ( cur == siteScores [ssIdx].end () ) {
					siteScores [ssIdx][site] = make_pair ( score, index );
				}
				else if ( score > (*cur).second.first ) {
					(*cur).second.first = score;
					(*cur).second.second = index;
				}
			}
		}
	}
}
void SiteScores::clear ()
{
	for ( int i = 0 ; i < siteScores.size () ; i++ ) {
		siteScores [i].clear ();
	}
}
void SiteScores::getSiteInfo ( StringVector& vModString, SiteInfoVectorVector& vvSiteInfo, GlycoSiteInfoVector& vGlycoSiteInfo ) const
{
	for ( MapStringToMapCharToIntConstIterator i = msmci.begin () ; i != msmci.end () ; i++ ) {
		string modString = (*i).first;
		bool glycoModFlag = glycoMods.find ( modString ) != glycoMods.end ();
		const MapCharToInt& mci = (*i).second;
		SiteInfoVector vSiteInfo;
		for ( MapCharToIntConstIterator j = mci.begin () ; j != mci.end () ; j++ ) {
			char aa = (*j).first;
			int idx = (*j).second;
			for ( MapIntToPairIntIntConstIterator k = siteScores [idx].begin () ; k != siteScores [idx].end () ; k++ ) {
				int site = (*k).first;
				int slip = (*k).second.first;
				int index = (*k).second.second;
				if ( slip ) {
					if ( glycoModFlag )
						vGlycoSiteInfo.push_back ( GlycoSiteInfo ( site, aa, slip, index, modString ) );	// Don't use ambiguous sites
					else
						vSiteInfo.push_back ( SiteInfo ( site, aa, slip, index ) );							// Don't use ambiguous sites
				}
			}
		}
		if ( !glycoModFlag ) {	// Not a glyco mod
			stable_sort ( vSiteInfo.begin (), vSiteInfo.end (), SiteInfoAscending () );
			vModString.push_back ( modString );
			vvSiteInfo.push_back ( vSiteInfo );
		}
	}
	if ( !vGlycoSiteInfo.empty () ) {
		stable_sort ( vGlycoSiteInfo.begin (), vGlycoSiteInfo.end (), GlycoSiteInfoAscending () );
	}
}
