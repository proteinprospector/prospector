/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mod_frag.cpp                                               *
*                                                                             *
*  Created    : July 17th 2001                                                *
*                                                                             *
*  Purpose    :                                                               *
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
#include <iomanip>
#include <lu_html.h>
#include <lu_indicies.h>
#include <lu_mod_frag.h>
#include <lu_mass_conv.h>
#ifdef CHEM_SCORE
#include <lu_chem_sc.h>
#endif
#include <lu_param_list.h>
#include <lu_inst.h>
#include <lu_usermod.h>
#include <lu_table.h>
#include <lu_delim.h>
using std::ostream;
using std::string;
using std::endl;
using std::fixed;
using std::setprecision;
using std::sort;
using std::vector;
using std::pair;
using std::set;
using std::make_pair;

class SortFragments {
	public:
		int operator () ( const PotentialMSFragment& a, const PotentialMSFragment& b ) const
		{
			if ( a.miMass == b.miMass ) {
				if ( a.fragMod == b.fragMod )
					return a.enzymeFragment.getStartAA () < b.enzymeFragment.getStartAA ();
				else
					return a.fragMod < b.fragMod;
			}
			else
				return a.miMass < b.miMass;
		}
};

DigestFragmentParameters::DigestFragmentParameters ( const ParameterList* params ) :
	minFragmentMass	( params->getDoubleValue	( "min_digest_fragment_mass", 500.0 ) ),
	maxFragmentMass	( params->getDoubleValue	( "max_digest_fragment_mass", 4000.0 ) ),
	minFragmentLength	( params->getIntValue	( "min_digest_fragment_length", 5 ) )
{
}
void DigestFragmentParameters::printHTML ( ostream& os ) const
{
	ParameterList::printHTML ( os, "Minimum Digest Fragment Mass", minFragmentMass );
	ParameterList::printHTML ( os, "Maximum Digest Fragment Mass", maxFragmentMass );
	ParameterList::printHTML ( os, "Minimum Digest Fragment Length", minFragmentLength );
}

bool PotentialMSFragment::monoFlag = false;
PotentialMSFragment::PotentialMSFragment ( const EnzymeFragment& enzymeFragment, const MapIntToInt& fragMod, int charge, int entryNumber ) :
	enzymeFragment ( enzymeFragment ),
	fragMod ( fragMod ),
	charge ( charge ),
	entryNumber ( entryNumber )
{
	if ( !monoFlag ) {
		modified_mass_convert ( true );
		monoFlag = true;
	}
	double mass = enzymeFragment.getMass () + FragModContainer::getMonoisotopicMass ( fragMod );
	miMass = mPlusHToMOverZ ( mass, charge, true );
#ifdef CHEM_SCORE
	if ( chemScore ) {
		cScore = chemScore->getChemScore ( enzymeFragment.getFragment (), FragModContainer::getNumOxidation ( fragMod ) );
	}
#endif
}
void PotentialMSFragment::setAverageMass ()
{
	if ( monoFlag ) {
		modified_mass_convert ( false );
		monoFlag = false;
	}
	double mass = enzymeFragment.getMass () + FragModContainer::getAverageMass ( fragMod );
	avMass = mPlusHToMOverZ ( mass, charge, false );
}
void PotentialMSFragment::printHeaderHTML ( ostream& os, bool modifications )
{
	tableRowStart ( os );
		tableHeader ( os, "Number" );
		tableHeader ( os, "m/z<br />(mi)" );
		tableHeader ( os, "m/z<br />(av)" );
#ifdef CHEM_SCORE
		if ( chemScore ) os << "<th>Chem<br />Score";
#endif
		if ( modifications )	tableHeader ( os, "Modifications" );
		if ( bullBreeseFlag )	tableHeader ( os, "BB<br />(%Hydrophobicity)" );
		if ( HPLCFlag )			tableHeader ( os, "HPLC" );
		EnzymeFragment::printHeaderHTML ( os );
	tableRowEnd ( os );
}
void PotentialMSFragment::printHeaderDelimited ( ostream& os, bool modifications )
{
	delimitedRowStart ( os );
		delimitedHeader ( os, "Number" );
		delimitedHeader ( os, "m/z (mi)" );
		delimitedHeader ( os, "m/z (av)" );
		delimitedHeader ( os, "Charge" );
#ifdef CHEM_SCORE
		if ( chemScore )		delimitedHeader ( os, "Chem Score" );
#endif
		if ( modifications )	delimitedHeader ( os, "Modifications" );
		if ( bullBreeseFlag )	delimitedHeader ( os, "BB (%Hydrophobicity)" );
		if ( HPLCFlag )			delimitedHeader ( os, "HPLC" );
		EnzymeFragment::printHeaderDelimited ( os );
	delimitedRowEnd ( os );
}
void PotentialMSFragment::printHTML ( ostream& os, bool modifications, bool hideLinks, const MSProductLink& productLink, int precision ) const
{
	tableRowStart ( os );
		tableCell ( os, entryNumber+1 );
		tableDataStart ( os, "", "left" );
			genPrint ( os, miMass, precision );
			print_charge_superscript ( os, charge );
			os << endl;
		tableDataEnd ( os );
		tableDataStart ( os, "", "left" );
			genPrint ( os, avMass, precision );
			print_charge_superscript ( os, charge );
			os << endl;
		tableDataEnd ( os );
#ifdef CHEM_SCORE
		if ( chemScore ) {
			tableDataStart ( os, "", "left" );
			genPrint ( os, cScore, 1 );
			tableDataEnd ( os );
		}
#endif
		if ( modifications ) FragModContainer::printHTML ( os, fragMod );

		string fragment = enzymeFragment.getFragment ();
		if ( bullBreeseFlag ) {
			setIndexType ( "Bull Breese" );
			double bb = indexScore ( fragment ) / fragment.length ();
			tableDataStart ( os, "", "middle" );
				os << fixed << setprecision ( 1 ) << ( bb - 970 ) * 100 / ( -1650 - 970 ) << endl;
			tableDataEnd ( os );
		}
		if ( HPLCFlag ) {
			setIndexType ( "HPLC Index" );
			tableDataStart ( os, "", "right" );
				os << fixed << setprecision ( 1 ) << indexScore ( fragment ) << endl;
			tableDataEnd ( os );
		}
		enzymeFragment.printHTML ( os, hideLinks || FragModContainer::getModified ( fragMod ), productLink, charge );
	tableRowEnd ( os );
}
void PotentialMSFragment::printDelimited ( ostream& os, int precision ) const
{
	delimitedRowStart ( os );
		delimitedCell ( os, entryNumber+1 );
		delimitedCell ( os, miMass, precision );
		delimitedCell ( os, avMass, precision );
		delimitedCell ( os, charge );
#ifdef CHEM_SCORE
		if ( chemScore ) delimitedCell ( os, cScore, 1 );
#endif
		FragModContainer::printDelimited ( os, fragMod );

		string fragment = enzymeFragment.getFragment ();
		if ( bullBreeseFlag ) {
			setIndexType ( "Bull Breese" );
			double bb = indexScore ( fragment ) / fragment.length ();
			delimitedCellSigFig ( os, ( bb - 970 ) * 100 / ( -1650 - 970 ), 3 );
		}
		if ( HPLCFlag ) {
			setIndexType ( "HPLC Index" );
			delimitedCellSigFig ( os, indexScore ( fragment ), 3 );
		}
		enzymeFragment.printDelimited ( os );
	delimitedRowEnd ( os );
}
void PotentialMSFragment::printXML ( ostream& os, int precision ) const
{
	os << "<peptide>" << endl;
	ParameterList::printXML ( os, "protein_index", entryNumber+1 );
	ParameterList::printDoubleXMLFixed ( os, "mi_m_over_z", miMass, precision );
	ParameterList::printDoubleXMLFixed ( os, "av_m_over_z", avMass, precision );
	ParameterList::printXML ( os, "charge", charge );
	FragModContainer::printXML ( os, fragMod );

#ifdef CHEM_SCORE
	if ( chemScore ) {
		ParameterList::printDoubleXMLFixed ( os, "chem_score", cScore, 1 );
	}
#endif
	string fragment = enzymeFragment.getFragment ();
	if ( bullBreeseFlag ) {
		setIndexType ( "Bull Breese" );
		double bb = indexScore ( fragment ) / fragment.length ();
		ParameterList::printDoubleXMLSigFig ( os, "bull_breese", ( bb - 970 ) * 100 / ( -1650 - 970 ), 3 );
	}
	if ( HPLCFlag ) {
		setIndexType ( "HPLC Index" );
		ParameterList::printDoubleXMLSigFig ( os, "hplc", indexScore ( fragment ), 3 );
	}
	enzymeFragment.printXML ( os );
	os << "</peptide>" << endl;
}
#ifdef CHEM_SCORE
ChemScore* PotentialMSFragment::chemScore = 0;
void PotentialMSFragment::setChemScore ( const string& cysMod, double metOxF )
{
	chemScore = new ChemScore ( cysMod, metOxF );
}
#endif
bool PotentialMSFragment::bullBreeseFlag = false;
void PotentialMSFragment::setBullBreese ( bool flag )
{
	bullBreeseFlag = flag;
}
bool PotentialMSFragment::HPLCFlag = false;
void PotentialMSFragment::setHPLC ( bool flag )
{
	HPLCFlag = flag;
}
#ifdef CHEM_SCORE
void PotentialMSFragment::deleteChemScore ()
{
	delete chemScore;
}
#endif
PotentialMSFragmentContainer::PotentialMSFragmentContainer ( const vector <EnzymeFragmentContainer>& enzFrags, const DigestFragmentParameters& fragParams, int maxPossCharge, int maxDigestHits ) :
	productLink ( "msdigest" ),
	precision ( instInf->getParentPeakPrecision ().getMassDecimalPlaces () ),
	maxDigestHits ( maxDigestHits )
{
	modifications = false;
	double minFragmentMass = fragParams.getMinFragmentMass ();
	double maxFragmentMass = fragParams.getMaxFragmentMass ();
	int minFragmentLength = fragParams.getMinFragmentLength ();
	for ( EnzymeFragmentContainerVectorSizeType i = 0 ; i < enzFrags.size () ; i++ ) {
		calculateFragments ( enzFrags [i], maxPossCharge, i, minFragmentMass, maxFragmentMass, minFragmentLength );
	}
	sort ( fragments.begin (), fragments.end (), SortFragments () );
	for ( int j = 0 ; j < fragments.size () ; j++ ) {
		fragments [j].setAverageMass ();
	}
}
PotentialMSFragmentContainer::PotentialMSFragmentContainer ( const EnzymeFragmentContainer& enzFrags, const DigestFragmentParameters& fragParams, int maxPossCharge, int maxDigestHits, int entryNumber ) :
	productLink ( "msdigest" ),
	precision ( instInf->getParentPeakPrecision ().getMassDecimalPlaces () ),
	maxDigestHits ( maxDigestHits )
{
	modifications = false;
	double minFragmentMass = fragParams.getMinFragmentMass ();
	double maxFragmentMass = fragParams.getMaxFragmentMass ();
	int minFragmentLength = fragParams.getMinFragmentLength ();
	calculateFragments ( enzFrags, maxPossCharge, entryNumber, minFragmentMass, maxFragmentMass, minFragmentLength );
	sort ( fragments.begin (), fragments.end (), SortFragments () );
	for ( int j = 0 ; j < fragments.size () ; j++ ) {
		fragments [j].setAverageMass ();
	}
}
void PotentialMSFragmentContainer::calculateFragments ( const EnzymeFragmentContainer& enzFrags, int maxPossCharge, int entryNumber, double minMass, double maxMass, int minLength )
{
	set <PairStringInt> uniqFrag; 
	for ( EnzymeFragmentIterator efi ( enzFrags ) ; efi.more () ; efi.advance () ) {
		const EnzymeFragment& ef = efi.getEnzFrag ();
		pair <std::set <PairStringInt>::iterator, bool> flag = uniqFrag.insert ( make_pair ( ef.getFragment (), ef.getStartAA () ) );
		if ( flag.second ) {
			double unmodifiedMass = ef.getMass ();
			if ( maxPossCharge == 0 || ( unmodifiedMass + FragModContainer::getMinModMi () <= maxMass && unmodifiedMass + FragModContainer::getMaxModMi () >= minMass ) ) {
				FragModContainer fmc ( ef.getFragment (), ef.getFirstFragment (), ef.getLastFragment () );
				int maxZ;
				if ( maxPossCharge )	// Impose a maximum charge
					maxZ = genMin ( maxPossCharge, ef.getMaxCharge () );
				else
					maxZ = ef.getMaxCharge ();
				for ( fmc.first () ; fmc.isDone () ; fmc.next () ) {
					const MapIntToInt& modList = fmc.getModList ();
					if ( !modifications && !modList.empty () ) modifications = true;
					for ( int z = 1 ; z <= maxZ ; z++ ) {
						PotentialMSFragment frag ( ef, modList, z, entryNumber );
						double mass = frag.getMiMass ();
						int length = frag.getLength ();
						if ( mass >= minMass && mass <= maxMass && length >= minLength ) {
							fragments.push_back ( frag );
							if ( fragments.size () > maxDigestHits ) {
								ErrorHandler::genError ()->error ( "The maximum number of hits has been exceeded.\n" );
							}
						}
					}
				}
			}
		}
	}
}
void PotentialMSFragmentContainer::printHTML ( ostream& os, bool hideLinks ) const
{
	startJavascript ( os );
	productLink.printHTML ( os );
	endJavascript ( os );

	os << "<table>" << endl;

	for ( PotentialMSFragmentVectorSizeType i = 0 ; i < fragments.size () ; i++ ) {
		if ( i == 0 ) PotentialMSFragment::printHeaderHTML ( os, modifications );
		fragments [i].printHTML ( os, modifications, hideLinks, productLink, precision );
	}
	os << "</table>" << endl;
	os << "<p />" << endl;
}
void PotentialMSFragmentContainer::printDelimited ( ostream& os ) const
{
	for ( PotentialMSFragmentVectorSizeType i = 0 ; i < fragments.size () ; i++ ) {
		fragments [i].printDelimited ( os, precision );
	}
}
void PotentialMSFragmentContainer::printXML ( ostream& os ) const
{
	for ( PotentialMSFragmentVectorSizeType i = 0 ; i < fragments.size () ; i++ ) {
		fragments [i].printXML ( os, precision );
	}
}
