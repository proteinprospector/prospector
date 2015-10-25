/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_iso.h                                                      *
*                                                                             *
*  Created    : September 22nd 1997                                           *
*                                                                             *
*  Purpose    : Functions for calculating isotopic distributions.             *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1997-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_iso_h
#define __lu_iso_h

#include <algorithm>
#include <nr.h>
#include <lu_formula.h>

struct Element {
	int index;
	int multiplier;
	int* isotopeAtoms;
};
typedef std::vector <Element> ElementVector;
typedef ElementVector::size_type ElementVectorSizeType;

struct IsotopicDistributionComponent {
	double probability;
	double massOffset;
	std::string formula;
};

typedef std::vector <IsotopicDistributionComponent> IsotopicDistributionComponentVector;
typedef IsotopicDistributionComponentVector::size_type IsotopicDistributionComponentVectorSizeType;

class IsotopicDistribution {
	ElementalFormula elementalFormula;
	std::vector <IsotopicDistributionComponent> isotopicDistributionComponent;
	double monoisotopicMass;
	double idealMonoisotopicMass;
	double idealMonoisotopicMZ;
	int charge;

	bool noDistribution;
	bool monoisotopicPeakTooSmall;
	double probabilityLimit;
	double totalAbundance;

	std::vector <Element> elements;

	static bool initialised;

	static StringVector elementTable;
	static int numValidElements;
	static DoubleVectorVector isotopeMass;
	static DoubleVectorVector isotopeAbundance;

	static void initialise ();
	static void initialiseIsoDistCalculations ();

	void calcIsotopeDistributionComponent ( int elementNumber, int isotopeNumber, double probability );
	void addIsotropicDistributionComponent ( double probability );
public:
	IsotopicDistribution ( ElementalFormula& ef, int charge, double probabilityLimit );
	IsotopicDistribution ( const std::string& eFormula, int charge, double probabilityLimit );
	IsotopicDistribution ( double mOverZ, int charge, const std::string& averagineFormulaStr, double probabilityLimit );
	~IsotopicDistribution ();
	void createIsotopicDistribution ( ElementalFormula& f );
	bool getNoDistribution () const { return noDistribution; }
	bool getMonoisotopicPeakTooSmall () const { return monoisotopicPeakTooSmall; }
	double getTotalAbundance () const { return totalAbundance; }
	ElementalFormula getElementalFormula () const { return elementalFormula; }
	int getCharge () const { return charge; }
	double getIdealMonoisotopicMass () const { return idealMonoisotopicMass; }
	double getIdealMonoisotopicMZ () const { return idealMonoisotopicMZ; }
	friend class IsotopicDistributionConstIterator;
};

class IsotopicDistributionConstIterator {
	const IsotopicDistributionComponentVector& isotopicDistributionComponent;
	IsotopicDistributionComponentVectorSizeType cur;
public:
	IsotopicDistributionConstIterator ( const IsotopicDistribution* isotopicDistribution ) :
		isotopicDistributionComponent ( isotopicDistribution->isotopicDistributionComponent ),
		cur ( 0 ) {}
	const double& probability () const { return isotopicDistributionComponent [cur].probability;}
	const double& massOffset () const { return isotopicDistributionComponent [cur].massOffset;}
	const std::string& formula () const { return isotopicDistributionComponent [cur].formula;}
	bool more () const { return cur < isotopicDistributionComponent.size (); }
	void advance () { cur++; }
};

class IsotopePeakStats : public IsotopicDistribution {
	static const double DEFAULT_PROBABILITY_LIMIT;
	static const std::string DEFAULT_AVERAGINE_FORMULA;
	DoubleVector averageMass;
	DoubleVector totalProbability;
	DoubleVector maximumProbability;
	IntVector numComponents;
	void createIsotopePeakStats ();
public:
	IsotopePeakStats ( ElementalFormula& ef, int charge, double probabilityLimit = DEFAULT_PROBABILITY_LIMIT );
	IsotopePeakStats ( const std::string& eFormula, int charge, double probabilityLimit = DEFAULT_PROBABILITY_LIMIT );
	IsotopePeakStats ( double mOverZ, int charge, const std::string& averagineFormulaStr = DEFAULT_AVERAGINE_FORMULA, double probabilityLimit = DEFAULT_PROBABILITY_LIMIT );
	int size () const { return averageMass.size (); }
	double getStartMass () const { return size () ? averageMass.front () : 0.0; }
	double getEndMass () const { return size () ? averageMass.back () : 0.0; }
	double getGroupMaximumProbability () const { return size () ? *(std::max_element ( totalProbability.begin (), totalProbability.end () )) : 0.0; }
	double getProbability ( int i ) const { return totalProbability [i]; }
	int getProbabilityIndexForMass ( double mass ) const;
	int getProbabilityIndexForIdealMonoisotopicMZ () const;
	DoubleVector getAverageMass () const { return averageMass; }
	DoubleVector getTotalProbability () const { return totalProbability; }
	int getNumPeaks () const { return averageMass.size (); }
	friend class IsotopePeakStatsConstIterator;
};

class IsotopePeakStatsConstIterator {
	const DoubleVector& averageMassList;
	const DoubleVector& totalProbabilityList;
	const DoubleVector& maximumProbabilityList;
	const IntVector& numComponentsList;
	DoubleVectorSizeType cur;
public:
	IsotopePeakStatsConstIterator ( const IsotopePeakStats* ips ) :
		averageMassList ( ips->averageMass ),
		totalProbabilityList ( ips->totalProbability ),
		maximumProbabilityList ( ips->maximumProbability ),
		numComponentsList ( ips->numComponents ),
		cur ( 0 ) {};
	const double& averageMass () const { return averageMassList [cur]; }
	const double& totalProbability () const { return totalProbabilityList [cur]; }
	const double& maximumProbability () const { return maximumProbabilityList [cur]; }
	const int numComponents () const { return numComponentsList [cur]; }
	bool more () const { return cur < averageMassList.size (); }
	void advance () { cur++; }
};

class IsotopeProfile : public XYData {
protected:
	static bool detailed;
public:
	static void getStartAndEndMass ( std::vector <IsotopePeakStats*>& ips, double& startMass, double& endMass );
	static double getStep ( double step, double resolution, double mass );
	static void setDetailed ( bool d ) { detailed = d; }
};
class GaussianIsotopeProfile : public IsotopeProfile {
	void simpleProfile ( double twoSDSquared, IsotopePeakStats* isotopePeakStats, const double& x, double& y );
	void detailedProfile ( double twoSDSquared, IsotopePeakStats* isotopePeakStats, const double& x, double& y );
public:
	GaussianIsotopeProfile ( std::vector <IsotopePeakStats*>& ips, const DoubleVector& intensity, double resolution );
};
class LorentzianIsotopeProfile : public IsotopeProfile {
	void simpleProfile ( double beta, double betaSquared, IsotopePeakStats* isotopePeakStats, const double& x, double& y );
	void detailedProfile ( double beta, double betaSquared, IsotopePeakStats* isotopePeakStats, const double& x, double& y );
public:
	LorentzianIsotopeProfile ( std::vector <IsotopePeakStats*>& ips, const DoubleVector& intensity, double resolution );
};
class StickIsotopeProfile : public IsotopeProfile {
	void simpleProfile ( IsotopePeakStats* isotopePeakStats );
	void detailedProfile ( IsotopePeakStats* isotopePeakStats );
public:
	StickIsotopeProfile ( std::vector <IsotopePeakStats*>& isotopePeakStats );
};
#endif /* ! __lu_iso_h */
