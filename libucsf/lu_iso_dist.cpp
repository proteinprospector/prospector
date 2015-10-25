/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_iso_dist.cpp                                               *
*                                                                             *
*  Created    : September 9th 1997                                            *
*                                                                             *
*  Purpose    : Functions for calculating isotopic distributions.             *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1997-2011) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef VIS_C
#include <stdexcept>
#endif
#include <lg_io.h>
#include <lg_string.h>
#include <lgen_error.h>
#include <lgen_math.h>
#include <lu_iso.h>
#include <lu_mass.h>
#include <lu_mass_conv.h>
#include <lu_mass_elem.h>
using std::vector;
using std::string;
using std::fill;
using std::sort;
using std::ostringstream;
using std::runtime_error;

class SortIsotopicDistributionComponents {
	public:
		bool operator () ( const IsotopicDistributionComponent& a, const IsotopicDistributionComponent& b ) const
		{
			return ( a.massOffset < b.massOffset );
		}
};
bool IsotopicDistribution::initialised = false;

StringVector IsotopicDistribution::elementTable;
int IsotopicDistribution::numValidElements;
DoubleVectorVector IsotopicDistribution::isotopeMass;
DoubleVectorVector IsotopicDistribution::isotopeAbundance;

IsotopicDistribution::IsotopicDistribution ( double mOverZ, int charge, const string& averagineFormulaStr, double probabilityLimit ) :
	charge ( charge ),
	probabilityLimit ( probabilityLimit )
{
	double mass = mOverZToMPlusH ( mOverZ, charge, true );
	Formula <double> averagineFormula ( averagineFormulaStr );
	double averagineMass = 0.0;
	for ( averagineFormula.first () ; averagineFormula.isDone () ; averagineFormula.next () ) {
		averagineMass += averagineFormula.multiplier () * formula_to_monoisotopic_mass ( averagineFormula.element ().c_str () );
	}
	averagineFormula.multiply ( mass / averagineMass );

	ElementalFormula ef;
	for ( averagineFormula.first () ; averagineFormula.isDone () ; averagineFormula.next () ) {
		string temp = averagineFormula.element () + gen_itoa ( (int) ( averagineFormula.multiplier () + 0.5 ) );
		ef += temp; 
	}
	double mass2 = mass - formula_to_monoisotopic_mass ( ef );
	int numHydrogens = mass2 > 0 ? (int) ( mass2 + 0.5 ) : (int) ( mass2 - 0.5 );
	string temp = "H" + gen_itoa ( numHydrogens );
	ef += temp;
	createIsotopicDistribution ( ef );
}
IsotopicDistribution::IsotopicDistribution ( const string& eFormula, int charge, double probabilityLimit ) :
	charge ( charge ),
	probabilityLimit ( probabilityLimit )
{
	ElementalFormula f ( eFormula );
	createIsotopicDistribution ( f );
}
IsotopicDistribution::IsotopicDistribution ( ElementalFormula& ef, int charge, double probabilityLimit ) :
	charge ( charge ),
	probabilityLimit ( probabilityLimit )
{
	createIsotopicDistribution ( ef );
}
void IsotopicDistribution::createIsotopicDistribution ( ElementalFormula& f )
{
	initialise ();

	elementalFormula = f;

	monoisotopicMass = formula_to_monoisotopic_mass ( f.getFormula ().c_str () ) - ELECTRON_REST_MASS;
	idealMonoisotopicMass = formula_to_ideal_monoisotopic_mass ( f.getFormula ().c_str () ) - ELECTRON_REST_MASS;
	idealMonoisotopicMZ = mPlusHToMOverZ ( idealMonoisotopicMass, charge, true );

	for ( f.first () ; f.isDone () ; f.next () ) {
		int j;
		if ( numValidElements == 0 ) {
			ErrorHandler::genError ()->error ( "No valid elements defined.\n" );
		}
		bool present = false;
		for ( j = 0 ; j < numValidElements ; j++ ) {
			if ( f.element () == elementTable [j] ) {
				present = true;
				break;
			}
		}
		if ( !present ) {
			ErrorHandler::genError ()->error ( "Element \"" + f.element () + "\" not currently available.\n" );
		}
		int numIsotopes = isotopeMass [j].size ();
		if ( numIsotopes > 1 ) {
			Element e;
			e.index = j;
			e.multiplier = f.multiplier ();
			if ( e.multiplier < 0 ) throw runtime_error ( "Elemental formula " + f.getFormula () + " has negative value for one or more elements." );
			e.isotopeAtoms = new int [numIsotopes];
			fill ( e.isotopeAtoms, e.isotopeAtoms + numIsotopes, 0 );
			elements.push_back ( e );
		}
	}

	monoisotopicPeakTooSmall = false;
	if ( elements.size () ) {
		calcIsotopeDistributionComponent ( 0, 1, 1.0 );
		noDistribution = false;
		sort ( isotopicDistributionComponent.begin (), isotopicDistributionComponent.end (), SortIsotopicDistributionComponents () );
		if ( isotopicDistributionComponent.size () == 0 )
			monoisotopicPeakTooSmall = true;
	}
	else
		noDistribution = true;

	totalAbundance = 0.0;
	for ( IsotopicDistributionComponentVectorSizeType k = 0 ; k < isotopicDistributionComponent.size () ; k++ )
		totalAbundance += isotopicDistributionComponent [k].probability;
}
IsotopicDistribution::~IsotopicDistribution ()
{
	for ( ElementVectorSizeType i = 0 ; i < elements.size () ; i++ )
		delete [] elements [i].isotopeAtoms;
}
void IsotopicDistribution::calcIsotopeDistributionComponent ( int elementNumber, int isotopeNumber, double probability )
{
	int index = elements [elementNumber].index;
	int multiplier = elements [elementNumber].multiplier;
	int* isotopeAtoms = elements [elementNumber].isotopeAtoms;
	int numClasses = isotopeMass [index].size () - 1;

	for ( int i = 0 ; i <= multiplier ; i++ ) {
		double prob;
		isotopeAtoms [isotopeNumber] = i;

		if ( isotopeNumber == numClasses ) {
			prob = gen_multinomial_probability ( multiplier, &(isotopeAbundance [index][1]), isotopeAtoms + 1, numClasses );
			probability *= prob;
		}
		if ( probability <= probabilityLimit ) break;
		ElementVectorSizeType nextElement = ( isotopeNumber == numClasses ) ? elementNumber + 1 : elementNumber;
		int nextIsotope = ( isotopeNumber == numClasses ) ? 1 : isotopeNumber + 1;
		if ( nextElement < elements.size () ) calcIsotopeDistributionComponent ( nextElement, nextIsotope, probability );
		else addIsotropicDistributionComponent ( probability );
		if ( isotopeNumber == numClasses ) probability /= prob;
	}
}
void IsotopicDistribution::addIsotropicDistributionComponent ( double probability )
{
	ostringstream formula;
	int* isotopeAtoms;
	double massOffset;
	int index;

	massOffset = 0.0;
	ElementVectorSizeType i;
	DoubleVectorVectorSizeType j;
	for ( i = 0 ; i < elements.size () ; i++ ) {
		index = elements [i].index;
		isotopeAtoms = elements [i].isotopeAtoms;
		for ( j = 1 ; j < isotopeMass [index].size () ; j++ ) {
			if ( isotopeAtoms [j] ) {
				genPrint ( formula, isotopeMass [index][j], 0 );
				formula << elementTable [index];
				formula << isotopeAtoms [j];
				formula << " ";
				massOffset += isotopeAtoms [j] * ( isotopeMass [index][j] - isotopeMass [index][0] );
			}
		}
	}
	IsotopicDistributionComponent idc;
	idc.probability = probability;
	idc.massOffset = mPlusHToMOverZ ( massOffset + monoisotopicMass, charge, true );
	string s = formula.str ();
	string::size_type len = s.length ();
	if ( len && s [len-1] == ' ' )	idc.formula = s.substr ( 0, len-1 );
	else							idc.formula = formula.str ();
	isotopicDistributionComponent.push_back ( idc );
}
void IsotopicDistribution::initialise ()
{
	if ( initialised == false ) initialiseIsoDistCalculations ();
}
void IsotopicDistribution::initialiseIsoDistCalculations ()
{
	get_elemental_info ( &elementTable, &isotopeMass, &isotopeAbundance );
	numValidElements = elementTable.size ();
	initialised = true;
}
const double IsotopePeakStats::DEFAULT_PROBABILITY_LIMIT = 0.0000001;
const string IsotopePeakStats::DEFAULT_AVERAGINE_FORMULA = "C4.9384 H7.7583 N1.3577 O1.4773 S0.0417";
// Default Averagine Mass = 111.0543052
IsotopePeakStats::IsotopePeakStats ( double mOverZ, int charge, const string& averagineFormulaStr, double probabilityLimit ) :
	IsotopicDistribution ( mOverZ, charge, averagineFormulaStr, probabilityLimit )
{
	createIsotopePeakStats ();
}
IsotopePeakStats::IsotopePeakStats ( const string& eFormula, int charge, double probabilityLimit ) :
	IsotopicDistribution ( eFormula, charge, probabilityLimit )
{
	createIsotopePeakStats ();
}
IsotopePeakStats::IsotopePeakStats ( ElementalFormula& ef, int charge, double probabilityLimit ) :
	IsotopicDistribution ( ef, charge, probabilityLimit )
{
	createIsotopePeakStats ();
}
void IsotopePeakStats::createIsotopePeakStats ()
{
	double groupSpacing = 0.5 / IsotopicDistribution::getCharge ();
	bool newSet = true;
	double intensity;
	double avMass;
	double maxProbability;
	int n;
	for ( IsotopicDistributionConstIterator i ( this ) ; i.more () ; ) {
		double probability = i.probability ();
		double currentMassOffset = i.massOffset ();
		if ( newSet ) {
			intensity = probability;
			avMass = probability * currentMassOffset;
			maxProbability = probability;
			n = 1;
			newSet = false;
		}
		else {
			intensity += probability;
			avMass += probability * currentMassOffset;
			maxProbability = genMax( probability, maxProbability );
			n++;
		}
		i.advance ();
		if ( !i.more () || i.massOffset () - currentMassOffset >= groupSpacing ) {
			averageMass.push_back ( avMass / intensity );
			totalProbability.push_back ( intensity );
			maximumProbability.push_back ( maxProbability );
			numComponents.push_back ( n );
			newSet = true;
		}
	}
}
int IsotopePeakStats::getProbabilityIndexForMass ( double mass ) const
{
	return genNearestIndex ( averageMass, mass );
}
int IsotopePeakStats::getProbabilityIndexForIdealMonoisotopicMZ () const
{
	return genNearestIndex ( averageMass, getIdealMonoisotopicMZ () );
}
bool IsotopeProfile::detailed = false;
double IsotopeProfile::getStep ( double step, double resolution, double mass )
{
	if ( step == 0.0 ) {
		double width = mass / resolution;
		double increment = width / 5.0;
		if ( increment < 0.01 ) {
			if ( increment < 0.005 ) {
				if ( increment < 0.002 ) {
					if ( increment < 0.001 ) {
						if ( increment < 0.0005 ) {
							if ( increment < 0.0002 ) {
								if ( increment < 0.0001 ) {
									if ( increment < 0.00005 ) {
										if ( increment < 0.00002 ) {
											return 0.00001;
										}
										else return 0.00002;
									}
									else return 0.00005;
								}
								else return 0.0001;
							}
							else return 0.0002;
						}
						else return 0.0005;
					}
					else return 0.001;
				}
				else return 0.002;
			}
			else return 0.005;
		}
		else return 0.01;
	}
	else return step;
}
void IsotopeProfile::getStartAndEndMass ( vector <IsotopePeakStats*>& ips, double& startMass, double& endMass )
{
	for ( int i = 0 ; i < ips.size () ; i++ ) {
		double s = ips [i]->getStartMass ();
		double e = ips [i]->getEndMass ();
		startMass	= ( i == 0 ) ? s : genMin ( s, startMass );
		endMass		= ( i == 0 ) ? e : genMax ( e, endMass );
	}
}
GaussianIsotopeProfile::GaussianIsotopeProfile ( vector <IsotopePeakStats*>& ips, const DoubleVector& intensity, double resolution )
{
	double startMass;
	double endMass;
	getStartAndEndMass ( ips, startMass, endMass );
	double deltaM = startMass / resolution;						// Assumes resolution is that of Mono mass
	double sdMultiplier = 2.0 * sqrt ( - 2.0 * log ( 0.5 ) );	// Resolution at FWHM
	double sd = deltaM / sdMultiplier;
	double twoSDSquared = 2.0 * sd * sd;
	double startX = startMass - ( 3.0 * sd ) - 1.0;
	double endX = endMass + ( 3.0 * sd );

	double step = getStep ( 0.0, ( resolution < 50000 ) ? 50000 : resolution, startMass );
	for ( double x = startX ; x < endX ; x += step ) {
		double y = 0.0;
		for ( int i = 0 ; i < ips.size () ; i++ ) {
			double localY = 0.0;
			if ( detailed )	detailedProfile ( twoSDSquared, ips [i], x, localY );
			else			simpleProfile ( twoSDSquared, ips [i], x, localY );
			localY *= intensity [i];
			y += localY;
		}
		if ( x == startX || y >= 1e-4 ) XYData::add ( x, y );
	}
}
void GaussianIsotopeProfile::simpleProfile ( double twoSDSquared, IsotopePeakStats* isotopePeakStats, const double& x, double& y )
{
	for ( IsotopePeakStatsConstIterator i ( isotopePeakStats ) ; i.more () ; i.advance () ) {
		double top = x - i.averageMass ();
		top *= top;
		y += i.totalProbability () * exp ( - top / twoSDSquared );
	}
}
void GaussianIsotopeProfile::detailedProfile ( double twoSDSquared, IsotopePeakStats* isotopePeakStats, const double& x, double& y )
{
	IsotopicDistributionConstIterator id ( isotopePeakStats );
	IsotopePeakStatsConstIterator ips ( isotopePeakStats );
	for ( int i = 0 ; id.more () ; id.advance (), i++ ) {
		if ( i == ips.numComponents () ) {
			ips.advance ();
			i = 0;
		}
		double top = x - id.massOffset ();
		top *= top;
		y += id.probability () * exp ( - top / twoSDSquared );
	}
}
LorentzianIsotopeProfile::LorentzianIsotopeProfile ( vector <IsotopePeakStats*>& ips, const DoubleVector& intensity, double resolution )
{
	double startMass;
	double endMass;
	getStartAndEndMass ( ips, startMass, endMass );
	double deltaM = startMass / resolution;						// Assumes resolution is that of Mono mass
	double beta = deltaM / 2.0;									// Resolution at FWHM
	double betaSquared = beta * beta;
	double height = 1.0;
	double startX = startMass - ( 6.0 * beta ) - 1.0;
	double endX = endMass + ( 6.0 * beta );

	double step = getStep ( 0.0, resolution, startMass );
	for ( double x = startX ; x < endX ; x += step ) {
		double y = 0.0;
		for ( int i = 0 ; i < ips.size () ; i++ ) {
			double localY = 0.0;
			if ( detailed )	detailedProfile ( beta, betaSquared, ips [i], x, localY );
			else			simpleProfile ( beta, betaSquared, ips [i], x, localY );
			localY *= intensity [i];
			y += localY;
		}
		if ( x == startX || y >= 1e-4 ) XYData::add ( x, y );
	}
}
void LorentzianIsotopeProfile::simpleProfile ( double beta, double betaSquared, IsotopePeakStats* isotopePeakStats, const double& x, double& y )
{
	using nr::pi;
	for ( IsotopePeakStatsConstIterator i ( isotopePeakStats ) ; i.more () ; i.advance () ) {
		double bottom = x - i.averageMass ();
		bottom *= bottom;
		bottom += betaSquared;
		bottom *= pi;
		y += i.totalProbability () * beta / bottom;
	}
}
void LorentzianIsotopeProfile::detailedProfile ( double beta, double betaSquared, IsotopePeakStats* isotopePeakStats, const double& x, double& y )
{
	using nr::pi;
	IsotopicDistributionConstIterator id ( isotopePeakStats );
	IsotopePeakStatsConstIterator ips ( isotopePeakStats );
	for ( int i = 0 ; id.more () ; id.advance (), i++ ) {
		if ( i == ips.numComponents () ) {
			ips.advance ();
			i = 0;
		}
		double bottom = x - id.massOffset ();
		bottom *= bottom;
		bottom += betaSquared;
		bottom *= pi;
		y += id.probability () * beta / bottom;
	}
}
StickIsotopeProfile::StickIsotopeProfile ( vector <IsotopePeakStats*>& isotopePeakStats )
{
	if ( detailed )	detailedProfile ( isotopePeakStats [0] );
	else			simpleProfile ( isotopePeakStats [0] );
}
void StickIsotopeProfile::simpleProfile ( IsotopePeakStats* isotopePeakStats )
{
	for ( IsotopePeakStatsConstIterator i ( isotopePeakStats ) ; i.more () ; i.advance () ) {
		if ( i.totalProbability () >= 1e-5 ) XYData::add ( i.averageMass (), i.totalProbability () );
	}
}
void StickIsotopeProfile::detailedProfile ( IsotopePeakStats* isotopePeakStats )
{
	IsotopicDistributionConstIterator id ( isotopePeakStats );
	IsotopePeakStatsConstIterator ips ( isotopePeakStats );
	for ( int i = 0 ; id.more () ; id.advance (), i++ ) {
		if ( i == ips.numComponents () ) {
			ips.advance ();
			i = 0;
		}
		if ( id.probability () >= 1e-5 ) XYData::add ( id.massOffset (), id.probability () );
	}
}
