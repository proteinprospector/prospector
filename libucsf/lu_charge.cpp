/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_charge.cpp                                                 *
*                                                                             *
*  Created    : January 28th 1997                                             *
*                                                                             *
*  Purpose    : Functions for dealing with multiply charged data.             *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1997-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef VIS_C
#include <stdexcept>
#endif
#include <cmath>
#include <ln_dot.h>
#include <lgen_define.h>
#include <lu_charge.h>
#include <lu_html.h>
#include <lu_html_form.h>
#include <lu_mass.h>
#include <lu_delim.h>
#include <lu_mass_conv.h>
#include <lu_cgi_val.h>
#include <lu_pk_filter.h>
#include <lu_df_info.h>
#include <lu_pq_vector.h>
#include <lu_param_list.h>
#include <lu_table.h>
#include <lu_distribution.h>
using std::ostream;
using std::string;
using std::fill;
using std::stable_partition;
using std::stable_sort;
using std::sort;
using std::endl;
using std::remove_if;
using std::runtime_error;
using std::min_element;

class sortPeaksByMass {
public:
	bool operator () ( const Peak* a, const Peak* b ) const
	{
		return ( a->mass < b->mass );
	}
	bool operator () ( const Peak& a, const Peak& b ) const
	{
		return ( a.mass < b.mass );
	}
};
class sortPeaksByMOverZ {
public:
	bool operator () ( const Peak* a, const Peak* b ) const
	{
		return ( a->mOverZ < b->mOverZ );
	}
	bool operator () ( const Peak& a, const Peak& b ) const
	{
		return ( a.mOverZ < b.mOverZ );
	}
};

class sortPeaksByDescendingIntensity {
	public:
		int operator () ( const Peak& a, const Peak& b ) const
		{
			return ( a.intensity > b.intensity );
		}
		int operator () ( const Peak* a, const Peak* b ) const
		{
			return ( a->intensity > b->intensity );
		}
};

PeakContainerInfo::PeakContainerInfo ( const ParameterList* params ) :
	toleranceInfo	( "ms_parent_mass", params ),
	massType		( params->getStringValue	( "parent_mass_convert", "monoisotopic" ) ),
	monoisotopicFlag( massType != "average" ),
	systematicError	( params->getDoubleValue ( "ms_parent_mass_systematic_error", 0.0 ) )
{
	const char* value;
	if ( params->getValue ( "average_to_mono_convert", value ) )	getPostQueryVector ( value, averageToMonoConvert, '\n' );
	try {
		if ( params->getValue ( "parent_contaminant_masses", value ) )	getPostQueryVector ( value, parentContaminantMasses, '\n' );
	}
	catch ( runtime_error ) {
		ErrorHandler::genError ()->error ( "The contaminant masses field is incorrectly formatted. There should be a one contaminant mass per line.\n" );
	}
	if ( monoisotopicFlag ) modified_mass_convert ( true );
	else modified_mass_convert ( false );
}
bool PeakContainer::multiChargeAssign = true;
PeakContainer::PeakContainer ( MSDataSetInfo* dsi, const MSPeakFilterOptions& msPeakFilterOptions, const PeakContainerInfo& peakContainerInfo ) :
	spectrumRetained	( true ),
	tolerance			( peakContainerInfo.getTolerance () ),
	monoisotopicFlag	( peakContainerInfo.getMonoisotopicFlag () )
{
	int numDataSets = dsi->getNumDataSets ();
	for ( int i = 0 ; i < numDataSets ; i++ ) {
		initMS ( dsi->getDataSet ( i ), msPeakFilterOptions, peakContainerInfo, numDataSets == 1 ? "" : dsi->getSpecID ( i ) );
	}
	sortPeaks ();
}
PeakContainer::PeakContainer ( MSDataPoint* dataPoint, const MSPeakFilterOptions& msPeakFilterOptions, const PeakContainerInfo& peakContainerInfo ) :
	spectrumRetained	( true ),
	tolerance			( peakContainerInfo.getTolerance () ),
	monoisotopicFlag	( peakContainerInfo.getMonoisotopicFlag () )
{
	initMS ( dataPoint, msPeakFilterOptions, peakContainerInfo, "" );
}
void getPrecursorMs ( DoubleVector& preM, DoubleVector& preM2, const Peak* parentPeak )
{
	int chPar = parentPeak->getCharge ();
	double parMass = parentPeak->getMass ();
	int maxPreRegionZ = genMin ( chPar, 5 );
	for ( int z = maxPreRegionZ ; z >= 1 ; z-- ) {
		preM.push_back ( mPlusHToMOverZ ( parMass - 61, chPar, z, true ) );
		preM2.push_back ( mPlusHToMOverZ ( parMass + 4, chPar, z, true ) );
//std::cout << "m" << z << "=" << preM.back () << " " << preM2.back () << "<br />" << std::endl;
	}
}
typedef std::vector <MapIntToInt> VectorMapIntToInt;
typedef VectorMapIntToInt::size_type VectorMapIntToIntSizeType;

typedef std::vector <MapIntToIntVector> VectorMapIntToIntVector;
typedef VectorMapIntToIntVector::const_iterator VectorMapIntToIntVectorConstIterator;
typedef VectorMapIntToIntVector::size_type VectorMapIntToIntVectorSizeType;
using std::cout;
using std::endl;
void getChargeList ( VectorMapIntToInt& chList, const DataFilePeakVector& pks, const DoubleVector& preM2 )
{
	chList.resize ( pks.size () );
	double maxZ2 = 10000.0;
	double maxZ3 = 10000.0;
	double maxZ4 = 10000.0;
	int last = preM2.size () - 1;
	if ( preM2.size () >= 2 ) maxZ2 = preM2 [last-1];
	if ( preM2.size () >= 3 ) maxZ3 = preM2 [last-2];
	if ( preM2.size () >= 4 ) maxZ4 = preM2 [last-3];
	for ( DataFilePeakVectorSizeType i = 0 ; i < pks.size () ; i++ ) {
		const DataFilePeak& lPeak = pks [i];
		PeakVectorSizeType j = i+1;
		MapIntToInt mii;
		double lowPkMZ = lPeak.getMOverZ ();
		for ( ; j < pks.size () ; j++ ) {
			const DataFilePeak& hPeak = pks [j];
			double z = 1.0025 / ( hPeak.getMOverZ () - lowPkMZ );
			if ( z < 0.95 ) break;	// No more potential isotope peaks
			double ratio = hPeak.getIntensity () / lPeak.getIntensity ();
			int chrg = -1;
			if ( z >= 0.95 && z <= 1.1 ) chrg = 1;
			else if ( z >= 1.92 && z <= 2.08 && lowPkMZ < maxZ2 ) chrg = 2;
			else if ( z >= 2.86 && z <= 3.23 && lowPkMZ < maxZ3 ) chrg = 3;
			else if ( z >= 3.70 && z <= 4.35 && lowPkMZ < maxZ4 ) chrg = 4;
			else if ( z > 4.35 ) chrg = 0;
			if ( chrg >= 1 ) mii [chrg] = j;
		}
		chList [i] = mii;
		//cout << i << '\t' << lPeak.getMOverZ () << '\t' << lPeak.getIntensity () << '\t';
		//for ( MapIntToIntConstIterator k = mii.begin () ; k != mii.end () ; k++ ) {
		//	cout << '[' << (*k).first << ',' << (*k).second << ']';
		//}
		//cout << "\t<br />" << endl;
	}
}
void printChargeTable ( const DataFilePeakVector& pks, const VectorMapIntToInt& chList )
{
	tableStart ( cout, true );
		tableRowStart ( cout );
			tableHeader ( cout, "Pk #" );
			tableHeader ( cout, "m/z" );
			tableHeader ( cout, "1" );
			tableHeader ( cout, "2" );
			tableHeader ( cout, "3" );
			tableHeader ( cout, "4" );
		tableRowEnd ( cout );
		for ( int i = 0 ; i < chList.size () ; i++ ) {
			const MapIntToInt& mii = chList [i];
			tableRowStart ( cout );
				tableCell ( cout, i+1 );
				tableCell ( cout, pks [i].getMOverZ (), 4 );
				int k = 0;
				for ( MapIntToIntConstIterator j = mii.begin () ; k < 4 ; k++ ) {
					if ( j != mii.end () ) {
						tableCell ( cout, (*j).first );
						j++;
					}
					else
						tableCell ( cout, "" );
				}
			tableRowEnd ( cout );
		}
	tableEnd ( cout );
}
void getClusters ( VectorMapIntToIntVector& clusters, VectorMapIntToInt& chList )
{
	clusters.resize ( chList.size () );
	for ( VectorMapIntToIntVectorSizeType aa = 0 ; aa < clusters.size () ; aa++ ) {
		clusters [aa][-2].push_back ( aa );
	}
	for ( VectorMapIntToIntSizeType i = 0 ; i < chList.size () ; i++ ) {
		const MapIntToInt& mii = chList [i];
		for ( MapIntToIntConstIterator j = mii.begin () ; j != mii.end () ; j++ ) {
			int charge = (*j).first;
			int idx = (*j).second;
			IntVector iv;
			if ( charge >= 1 ) {
				iv.push_back ( i );
				iv.push_back ( idx );
				for ( ; ; ) {
					MapIntToInt& mii2 = chList [idx];
					MapIntToIntIterator cur = mii2.find ( charge );
					if ( cur != mii2.end () ) {
						idx = (*cur).second;
						iv.push_back ( idx );
						mii2.erase ( cur );
					}
					else break;
				}
				clusters [i][charge] = iv;
				for ( int k = 0 ; k < iv.size () ; k++ ) clusters [iv[k]].erase ( -2 );
			}
		}
	}
}
void getExperimentalNormDist ( DoubleVector& eDist, const DataFilePeakVector& pks, const IntVector& iv, int offset )
{
	double norm = 0.0;
	for ( IntVectorSizeType i = offset ; i < iv.size () ; i++ ) {
		double d = pks [iv[i]].getIntensity ();
		if ( d ) norm += d * d;
	}
	if ( norm ) {
		norm = sqrt ( norm );
		for ( IntVectorSizeType j = offset ; j < iv.size () ; j++ ) {
			eDist.push_back ( pks [iv[j]].getIntensity () / norm );
		}
	}
}
void getExperimentalNormDist2 ( DoubleVector& eDist, const DataFilePeakVector& pks, const IntVector& iv, int offset )
{
	double norm = 0.0;
	for ( IntVectorSizeType i = /*iv[offset]*/iv[0] ; i <= iv[iv.size ()-1] ; i++ ) {
		double d = pks [i].getIntensity ();
		if ( d ) norm += d * d;
	}
	if ( norm ) {
		norm = sqrt ( norm );
		for ( IntVectorSizeType j = offset ; j < iv.size () ; j++ ) {
			eDist.push_back ( pks [iv[j]].getIntensity () / norm );
		}
	}
}
void zeroPks ( DataFilePeakVector& pks, const IntVector& iv )
{
	for ( IntVectorSizeType i = 0 ; i < iv.size () ; i++ ) {
		pks [iv[i]].setZeroIntensity ();
	}
}
//typedef std::map <PairIntDouble, IntVector> MapPairIntDoubleToIntVector;
//typedef std::vector <MapPairIntDoubleToIntVector> VectorMapPairIntDoubleToIntVector;
//typedef VectorMapPairIntDoubleToIntVector::size_type VectorMapPairIntDoubleToIntVectorSizeType;

class Cluster {
	double mz;
	int ch;
	double monoInt;
	double maxInt;
	IntVector iv;
	DoubleVector dist;
	double mOverZ;
	double intensity;
	static TheoreticalDistribution* td;
public:
	Cluster ( double mz, int z, double monoInt, double maxInt, const IntVector& iv );
	double getDotProduct ( DataFilePeakVector& pks );
	double getMOverZ () const { return mOverZ; }
	int getCharge () const { return ch; }
	double getIntensity () const { return intensity; }
};
TheoreticalDistribution* Cluster::td = 0;
#include <lu_getfil.h>

Cluster::Cluster ( double mz, int z, double monoInt, double maxInt, const IntVector& iv ) :
	mz ( mz ),
	ch ( z == -2 ? 1 : z ),
	monoInt ( monoInt ),
	maxInt ( maxInt ),
	iv ( iv )
{
	if ( !td ) td = new TheoreticalDistribution ( "Averagine" );
	dist = td->getNormDistribution ( mz, ch );
}
double Cluster::getDotProduct ( DataFilePeakVector& pks )
{
	double maxCS = 0.0;
	for ( int i = 0 ; i < iv.size () ; i++ ) {
		DoubleVector eDist;
		getExperimentalNormDist2 ( eDist, pks, iv, i );
		double cs = dotProduct ( dist, eDist );
		if ( cs > maxCS ) {
			maxCS = cs;
			mOverZ = pks[iv[i]].getMOverZ ();
			intensity = pks[iv[i]].getIntensity () / dist [0];
		}
//std::cout << "mz=" << mz << " ch=" << ch << " iso=" << i << " cs=" << cs << "<br />" << std::endl;
	}
	return maxCS;
}
typedef std::vector <Cluster> VectorCluster;
typedef VectorCluster::size_type VectorClusterSizeType;

class CurrentClusterSet {
	DataFilePeakVector& pks;
	TheoreticalDistribution td;
	double minMass;
	double maxMass;
	VectorCluster clusters;
	DoubleVector preM;
	DoubleVector preM2;
	int preIdx;
public:
	CurrentClusterSet ( DataFilePeakVector& pks, const DoubleVector& preM, const DoubleVector& preM2 );
	void add ( double mz, MapIntToIntVector& c );
	void process ( DataFilePeakVector& dfpv );
};
CurrentClusterSet::CurrentClusterSet ( DataFilePeakVector& pks, const DoubleVector& preM, const DoubleVector& preM2  ) :
	pks ( pks ),
	td ( "Averagine" ),
	preM ( preM ),
	preM2 ( preM2 ),
	preIdx ( 0 )
{
}
void CurrentClusterSet::add ( double mz, MapIntToIntVector& c )
{
	if ( clusters.empty () ) minMass = mz;
	maxMass = mz;
	for ( MapIntToIntVector::reverse_iterator i = c.rbegin () ; i != c.rend () ; i++ ) {
		int ch = (*i).first;
		const IntVector& iv = (*i).second;
		double mx = 0.0;
		for ( int j = 0 ; j < iv.size () ; j++ ) {
			mx = genMax ( mx, pks [iv[j]].getIntensity () );
		}
//cout << "mz=" << mz << " ch=" << ch << " siz=" << iv.size () << "<br />" << endl;
		clusters.push_back ( Cluster ( mz, ch, pks [iv[0]].getIntensity (), mx, iv ) );
	}
}
void CurrentClusterSet::process ( DataFilePeakVector& dfpv )
{
static double maxRange = 0.0;
if ( maxMass-minMass > maxRange ) {
	maxRange = maxMass-minMass;
//	pidLogOutput ( gen_ftoa ( maxRange, "%.4f" ) + " " + gen_ftoa ( minMass, "%.4f" ) );
}
//cout << "min mass=" << minMass << " range=" << maxMass - minMass << " num clusters=" << clusters.size () << "<br />" << endl;
	double maxCS = 0.0;
	int maxIdx = -1;
	for ( VectorClusterSizeType a = 0 ; a < clusters.size () ; a++ ) {
		double cs = clusters [a].getDotProduct ( pks );
//pidLogOutput ( gen_ftoa ( clusters [a].getMOverZ (), "%.4f" ) + " " + gen_ftoa ( cs, "%.4f" ) );
		if ( cs > maxCS ) {
			maxIdx = a;
			maxCS = cs;
		}
	}
	if ( maxIdx != -1 ) {
		double mz = clusters [maxIdx].getMOverZ ();
		int ch = clusters [maxIdx].getCharge ();
//pidLogOutput ( gen_ftoa ( mz, "%.4f" ) + " " + gen_itoa ( ch ) );
		double inten = clusters [maxIdx].getIntensity ();
		if ( mz > preM2 [preIdx] ) preIdx++;
		int preChrg = preM2.size () - preIdx;

//pidLogOutput ( gen_ftoa ( mz, "%.4f" ) + " " + gen_ftoa ( maxCS, "%.4f" ) );
		if ( ( mz < preM [preIdx] || ch != preChrg ) && ( maxCS > 0.85 || ch != 1 ) && clusters.size () < 8 ) {
//tableRowStart ( cout );
//	tableCell ( cout, mz, 4 );
//	tableCell ( cout, ch );
//	tableCellSigFig ( cout, inten, 3 );
//	tableCellSigFig ( cout, maxCS, 3 );
//	tableCell ( cout, clusters.size () );
//tableRowEnd ( cout );
			dfpv.push_back ( DataFilePeak ( mz, ch, inten ) );
		}
	}
	clusters.clear ();
}
DataFilePeakVector makeHRETDPeakList ( const Peak* parentPeak, DataFilePeakVector& pks )
{
// 1st isotope maximum until 1770,885, 590, 442,354
// 2nd isotope maximum until 3320,1660,1106,830,664
	DoubleVector preM;
	DoubleVector preM2;
	getPrecursorMs ( preM, preM2, parentPeak );

	VectorMapIntToInt chList;
	getChargeList ( chList, pks, preM2 );
	//printChargeTable ( pks, chList );

	VectorMapIntToIntVector clusters;
	getClusters ( clusters, chList );

	int preIdx = 0;
	DataFilePeakVector dfpv;
	double lastMZ = 0.0;
//tableStart ( cout, true );
//tableRowStart ( cout );
//	tableHeader ( cout, "m/z" );
//	tableHeader ( cout, "z" );
//	tableHeader ( cout, "Int" );
//	tableHeader ( cout, "CS" );
//	tableHeader ( cout, "Num Dist" );
//tableRowEnd ( cout );
	CurrentClusterSet ccs ( pks, preM, preM2 );
	for ( VectorMapIntToIntVectorSizeType a = 0 ; a < clusters.size () ; a++ ) {
		double mz = pks [a].getMOverZ ();
		if ( a > 0 && mz > lastMZ + 1.055 ) ccs.process ( dfpv );
		if ( mz > preM2.back () ) {
			ccs.process ( dfpv );
			break;
		}
		if ( mz > preM2 [preIdx] ) preIdx++;
		/*if ( mz < preM [preIdx] )*/ ccs.add ( mz, clusters [a] );
		if ( a == clusters.size () - 1 ) ccs.process ( dfpv );
		lastMZ = mz;
	}
//tableEnd ( cout );
	return dfpv;
}
PeakContainer::PeakContainer ( MSMSDataPoint* dataPoint, const MSMSPeakFilterOptions* msmsPeakFilterOptions, const Peak* parentPeak, const Tolerance* parTol, const Tolerance* tolerance, bool monoisotopicFlag, bool averageToMonoConvertFlag ) :
	spectrumRetained ( true ),
	tolerance ( tolerance ),
	monoisotopicFlag ( monoisotopicFlag )
{
	if ( parentPeak->getMassPlusTol () >= msmsPeakFilterOptions->getMinPrecursorMass () ) {
		DataFilePeakVector& dataPeaks = dataPoint->getDataPeaks ();
		if ( msmsPeakFilterOptions->getRawSpectrum () ) {
			PeakVector pks;
			IntVector averageToMonoConvertArray ( dataPoint->size () );
			fill ( averageToMonoConvertArray.begin (), averageToMonoConvertArray.end (), averageToMonoConvertFlag );
			getPeaks ( pks, dataPeaks, tolerance, monoisotopicFlag, averageToMonoConvertArray );
			makePeaks ( pks );
		}
		else {
			if ( msmsPeakFilterOptions->getHighResETDDeisotope () ) {
				dataPeaks = makeHRETDPeakList ( parentPeak, dataPeaks );	// For new centroiding
			}
			if ( msmsPeakFilterOptions->getFTPeakExclusion () ) {
				removeFTPeaks ( dataPeaks, parentPeak, tolerance );
			}
			if ( msmsPeakFilterOptions->getECDorETDSideChainExclusion () ) {
				removeECDorETDSideChainPeaks ( dataPeaks, parentPeak, parTol, tolerance );
			}
			filterQuantitationPeaks ( dataPeaks );
			if ( msmsPeakFilterOptions->getPeakExclusion () ) {
				filterIntensity ( dataPeaks, msmsPeakFilterOptions->getMinIntensity () );
			}
			if ( msmsPeakFilterOptions->getJoinPeaks () ) {
				dataPeaks = joinSplitPeaks ( dataPeaks );
			}
			PeakVector pks;
			IntVector averageToMonoConvertArray ( dataPoint->size () );
			fill ( averageToMonoConvertArray.begin (), averageToMonoConvertArray.end (), averageToMonoConvertFlag );
			getPeaks ( pks, dataPeaks, tolerance, monoisotopicFlag, averageToMonoConvertArray );
			if ( msmsPeakFilterOptions->getDeisotopeHiRes () ) {
				pks = deisotopeHighResolution ( pks );	// For new centroiding
			}
			if ( msmsPeakFilterOptions->getMatrixExclusion () ) {
				removeMatrix ( pks, msmsPeakFilterOptions->getMaxMatrixMass () );
			}
			if ( !multiChargeAssign && msmsPeakFilterOptions->getDeisotopeHiRes () ) {	// For high res deisotoping matrix removal needs the charge state
				pks = setSingleCharge ( pks );
			}
			if ( msmsPeakFilterOptions->getDeisotope () ) {
				pks = deisotope ( pks );
			}
			if ( msmsPeakFilterOptions->getMassExclusion () ) {
				filterMSMSPeaks ( pks, msmsPeakFilterOptions->getMinMass (), msmsPeakFilterOptions->getPrecursorExclusion (), parentPeak );
			}
			if ( msmsPeakFilterOptions->getPeakExclusion () ) {
				if ( pks.size () ) {
					if ( msmsPeakFilterOptions->getPer100Da () )
						retainNPeaksPerDaRangeMSMS ( pks, msmsPeakFilterOptions->getMaxPeaks (), 100, msmsPeakFilterOptions->getMinPeaks () );
					else
						retainPeaksMSMS ( pks, msmsPeakFilterOptions->getMaxPeaks (), msmsPeakFilterOptions->getMinPeaks () );
				}
				if ( pks.size () < msmsPeakFilterOptions->getMinPeaks () ) {
					spectrumRetained = false;
				}
			}
			makePeaks ( pks );
		}
	}
	else {
		spectrumRetained = false;
	}
}
PeakContainer::~PeakContainer ()
{
	for ( PeakVectorSizeType i = 0 ; i < peaks.size () ; i++ ) {
		delete peaks [i];
	}
}
void PeakContainer::initMS ( MSDataPoint* dataPoint, const MSPeakFilterOptions& msPeakFilterOptions, const PeakContainerInfo& peakContainerInfo, const string& specNum )
{
	PeakVector pks;
	DataFilePeakVector& dataPeaks = dataPoint->getDataPeaks ();
	if ( msPeakFilterOptions.getPeakExclusion () ) {
		filterIntensity ( dataPeaks, msPeakFilterOptions.getMinIntensity () );
	}
	IntVector averageToMonoConvert	= peakContainerInfo.getAverageToMonoConvert ();
	double systematicError			= peakContainerInfo.getSystematicError ();
	getPeaks ( pks, dataPeaks, tolerance, monoisotopicFlag, averageToMonoConvert, systematicError, specNum );
	if ( msPeakFilterOptions.getMatrixExclusion () ) {
		removeMatrix ( pks, msPeakFilterOptions.getMaxMatrixMass () );
	}
	deleteContaminantPeaks ( pks, peakContainerInfo.getParentContaminantMasses () );
	if ( msPeakFilterOptions.getMassExclusion () ) {
		filterPeaks ( pks, msPeakFilterOptions.getMinMass (), msPeakFilterOptions.getMaxMass () );
	}
	if ( msPeakFilterOptions.getPeakExclusion () ) {
		retainPeaks ( pks, msPeakFilterOptions.getMaxPeaks (), msPeakFilterOptions.getMinPeaks () );
		if ( pks.size () < msPeakFilterOptions.getMinPeaks () ) {
			spectrumRetained = false;
		}
	}
	makePeaks ( pks );
}
void PeakContainer::makePeaks ( PeakVector& pks )
{
	for ( PeakVectorSizeType i = 0 ; i < pks.size () ; i++ ) {
		const Peak& p = pks [i];
		peaks.push_back ( new Peak ( p.mOverZ, p.tolerance, p.charge, p.intensity, p.mass, p.adductMass, p.averageMass, p.specNumber ) );
	}
}
void PeakContainer::truncate ( int size )
{
	peaks.erase ( peaks.begin () + size, peaks.end () );
}
void PeakContainer::sortPeaks ()
{
	sort ( peaks.begin (), peaks.end (), sortPeaksByMass () );
}
void PeakContainer::getPeaks ( PeakVector& pks, const DataFilePeakVector& dataPeaks, const Tolerance* tolerance, bool monoisotopicFlag, const IntVector& averageToMonoConvertArray, double systematicError, const string& specNumber )
{
	double adductMass = getAdductMass ( monoisotopicFlag );

	for ( DataFilePeakVectorSizeType i = 0 ; i < dataPeaks.size () ; i++ ) {
		double mOZ = dataPeaks [i].getMOverZ ();
		int charge = dataPeaks [i].getCharge ();
		if ( systematicError ) mOZ -= tolerance->getCorrection ( mOZ, systematicError );
		double intensity = dataPeaks [i].getIntensity ();
		double tol = tolerance->getTolerance ( mOZ, charge );

		if ( averageToMonoConvertArray.size () )
			pks.push_back ( Peak ( mOZ, tol, charge, intensity, adductMass, averageToMonoConvertArray [i] ) );
		else
			pks.push_back ( Peak ( mOZ, tol, charge, intensity, adductMass ) );
		pks [i].setSpecNumber ( specNumber );
	}
	stable_sort ( pks.begin (), pks.end (), sortPeaksByMass () );
}
bool PeakContainer::getPrecursorLessThanFragments ( double mass ) const
{
	if ( !peaks.empty () ) {
		return mass < peaks [0]->getMass ();
	}
	return false;
}
class CheckContaminant {
	const DoubleVector& contaminantPeaks;
public:
	CheckContaminant ( const DoubleVector& contaminantPeaks ) :
		contaminantPeaks ( contaminantPeaks ) {}
	bool operator () ( const Peak& peak )
	{
		for ( DoubleVectorSizeType i = 0 ; i < contaminantPeaks.size () ; i++ ) {
			if ( peak.isMatch ( contaminantPeaks [i] ) ) {
				return true;
			}
		}
		return false;
	}
};
void PeakContainer::deleteContaminantPeaks ( PeakVector& pks, const DoubleVector& contaminantPeaks )
{
	if ( !contaminantPeaks.empty () ) {
		int i = remove_if ( pks.begin (), pks.end (), CheckContaminant ( contaminantPeaks ) ) - pks.begin ();
		pks.erase ( pks.begin () + i, pks.end () );
	}
}
void PeakContainer::filterMSMSPeaks ( PeakVector& pks, double minMass, double precursorExclusion, const Peak* parentPeak )
{
	double maxMass = parentPeak->getMassPlusTol () - precursorExclusion;
	filterPeaks ( pks, minMass, maxMass );
}
DataFilePeakVector PeakContainer::joinSplitPeaks ( const DataFilePeakVector& pks )
{
	DataFilePeakVector newPks;
	for ( DataFilePeakVectorSizeType i = 0 ; i < pks.size () ; i++ ) {
		const DataFilePeak& p1 = pks [i];
		double interval = splitInterval ( p1.getMOverZ () );
		if ( i < pks.size () - 1 && pks [i+1].getMOverZ () - p1.getMOverZ () < interval ) {
			double sumOfWeightedMasses = p1.getMOverZ () * p1.getIntensity ();
			double sumOfIntensities = p1.getIntensity ();
			double maxIntensity = p1.getIntensity ();
			DataFilePeakVectorSizeType j = i+1;
			for ( ; j < pks.size () ; j++ ) {
				const DataFilePeak& p2 = pks [j];
				if ( p2.getMOverZ () - p1.getMOverZ () < interval ) {
					sumOfWeightedMasses += p2.getMOverZ () * p2.getIntensity ();
					sumOfIntensities += p2.getIntensity ();
					maxIntensity = genMax ( maxIntensity, p2.getIntensity () );
				}
				else {
					break;
				}
			}
			i = j - 1;
			double mOverZ = sumOfWeightedMasses / sumOfIntensities;
			newPks.push_back ( DataFilePeak ( mOverZ, p1.getCharge (), maxIntensity ) );
		}
		else {
			newPks.push_back ( p1 );
		}
	}
	return newPks;
}
PeakVector PeakContainer::deisotope ( const PeakVector& pks )
{
	PeakVector newPks;
	for ( int i = pks.size () ; i-- ; ) {
		double diff;
		int index = i;
		for ( int j = i ; j-- ; ) {
			const Peak& hPeak = pks [index];
			const Peak& lPeak = pks [j];
			if ( hPeak.isInRange ( lPeak.getMOverZ (), 1.0 ) ) {
				if ( hPeak.isAtLeastRatio ( lPeak.getIntensity (), 0.6 ) ) { // lpeak int >= hpeak int * 0.6
					diff = hPeak.getMOverZ () - lPeak.getMOverZ ();
					index = j;
				}
			}
			else break;
		}
		int charge = 1;
		if ( multiChargeAssign ) {
			if ( index != i ) {
				if ( diff < 0.7 ) charge = 2;
			}
		}
		if ( charge == 1 ) newPks.push_back ( pks [index] );
		else newPks.push_back ( Peak ( pks [index].mOverZ, tolerance->getTolerance ( pks [index].mOverZ, charge ), charge, pks [index].intensity, pks [index].adductMass ) );
		i = index;
	}
	return newPks;
}
PeakVector PeakContainer::deisotopeHighResolution ( const PeakVector& pks )
{
	PeakVector newPks;
	for ( int i = pks.size () ; i-- ; ) {
		double diff;
		int index = i;
		for ( int j = i ; j-- ; ) {
			const Peak& hPeak = pks [index];
			const Peak& lPeak = pks [j];
			if ( hPeak.isInRange ( lPeak.getMOverZ (), 1.0 ) ) {
				if ( hPeak.isAtLeastRatio ( lPeak.getIntensity (), 0.2 ) ) { // lpeak int >= hpeak int * 0.6
					diff = hPeak.getMOverZ () - lPeak.getMOverZ ();
					index = j;
				}
			}
			else break;
		}
		int charge = 1;
		if ( index != i ) {
			if ( diff < 0.52 && diff > 0.48 ) charge = 2;
			else if ( diff < 0.35 && diff > 0.31 ) charge = 3;
			else if ( diff < 0.27 && diff > 0.23 ) charge = 4;
		}
		if ( charge == 1 ) newPks.push_back ( pks [index] );
		else newPks.push_back ( Peak ( pks [index].mOverZ, tolerance->getTolerance ( pks [index].mOverZ, charge ), charge, pks [index].intensity, pks [index].adductMass ) );
		i = index;
	}
	return newPks;
}
PeakVector PeakContainer::setSingleCharge ( const PeakVector& pks )
{
	PeakVector newPks;
	for ( int i = pks.size () ; i-- ; ) {
		int charge = pks [i].getCharge ();
		if ( charge == 1 ) newPks.push_back ( pks [i] );
		else newPks.push_back ( Peak ( pks [i].mOverZ, tolerance->getTolerance ( pks [i].mOverZ, 1 ), 1, pks [i].intensity, pks [i].adductMass ) );
	}
	return newPks;
}
class RemoveAboveMOverZ {
	double maxMOverZ;
public:
	RemoveAboveMOverZ ( double maxMOverZ ) :
		maxMOverZ ( maxMOverZ ) {}
	bool operator () ( const DataFilePeak& peak )
		{ return peak.getMOverZ () > maxMOverZ; }
};
void PeakContainer::removeAboveMOverZ ( DataFilePeakVector& pks, double maxMOverZ )
{
	int i = remove_if ( pks.begin (), pks.end (), RemoveAboveMOverZ ( maxMOverZ ) ) - pks.begin ();
	pks.erase ( pks.begin () + i, pks.end () );
}
class RemoveMOverZRange {
	double minMOverZ;
	double maxMOverZ;
public:
	RemoveMOverZRange ( double minMOverZ, double maxMOverZ ) :
		minMOverZ ( minMOverZ ), maxMOverZ ( maxMOverZ ) {}
	bool operator () ( const DataFilePeak& peak )
		{ return peak.getMOverZ () >= minMOverZ && peak.getMOverZ () <= maxMOverZ; }
};
void PeakContainer::removeMOverZRange ( DataFilePeakVector& pks, double minMOverZ, double maxMOverZ )
{
	int i = remove_if ( pks.begin (), pks.end (), RemoveMOverZRange ( minMOverZ, maxMOverZ ) ) - pks.begin ();
	pks.erase ( pks.begin () + i, pks.end () );
}

class CheckMassRange {
	double minM;
	double maxM;
public:
	CheckMassRange ( double minM, double maxM ) :
		minM ( minM ), maxM ( maxM ) {}
	bool operator () ( const Peak& peak )
		{ return peak.getMassPlusTol () < minM || peak.getMassPlusTol () > maxM; }
};
void PeakContainer::filterPeaks ( PeakVector& pks, double minMass, double maxMass )
{
	int i = remove_if ( pks.begin (), pks.end (), CheckMassRange ( minMass, maxMass ) ) - pks.begin ();
	pks.erase ( pks.begin () + i, pks.end () );
}
class CheckQuantitationPeaks {
public:
	bool operator () ( const DataFilePeak& peak )
	{
		return	( peak.getMOverZ () > 112.8 && peak.getMOverZ () < 119.4 ) ||
				( peak.getMOverZ () > 120.8 && peak.getMOverZ () < 121.4 );
	}
};
void PeakContainer::filterQuantitationPeaks ( DataFilePeakVector& pks )
{
	int i = remove_if ( pks.begin (), pks.end (), CheckQuantitationPeaks () ) - pks.begin ();
	pks.erase ( pks.begin () + i, pks.end () );
}
class CheckIntensity {
	double minIntensity;
public:
	CheckIntensity ( double minIntensity ) :
		minIntensity ( minIntensity ) {}
	bool operator () ( const DataFilePeak& peak )
		{ return peak.getIntensity () < minIntensity; }
};
void PeakContainer::filterIntensity ( DataFilePeakVector& pks, double minIntensity )
{
	int i = remove_if ( pks.begin (), pks.end (), CheckIntensity ( minIntensity ) ) - pks.begin ();
	pks.erase ( pks.begin () + i, pks.end () );
}
void PeakContainer::retainPeaks ( PeakVector& pks, PeakVectorSizeType maxPeaks, PeakVectorSizeType minPeaks )
{
	sort ( pks.begin (), pks.end (), sortPeaksByDescendingIntensity () );
	if ( pks.size () > maxPeaks ) {
		pks.erase ( pks.begin () + maxPeaks, pks.end () );
	}
	sort ( pks.begin (), pks.end (), sortPeaksByMass () );
	if ( pks.size () < minPeaks ) {
		pks.erase ( pks.begin (), pks.end () );
	}
}
class LessThanMOverZ {
	double mOverZ;
public:
	LessThanMOverZ ( double mOverZ ) :
		mOverZ ( mOverZ ) {}
	bool operator () ( const Peak& peak )
		{ return peak.getMOverZ () < mOverZ; }
};
void PeakContainer::retainPeaksMSMS ( PeakVector& pks, PeakVectorSizeType maxPeaks, PeakVectorSizeType minPeaks )
{
	sort ( pks.begin (), pks.end (), sortPeaksByMOverZ () );
	double halfwayMOverZ = ( pks.back ().getMOverZ () - pks.front ().getMOverZ () ) / 2.0;
	PeakVectorIterator middle = stable_partition ( pks.begin (), pks.end (), LessThanMOverZ ( halfwayMOverZ ) );
	sort ( pks.begin (), middle, sortPeaksByDescendingIntensity () );
	sort ( middle, pks.end (), sortPeaksByDescendingIntensity () );
	int highRangeMaxPks = static_cast<int> ( maxPeaks * 0.5 );
	int lowRangeMaxPks = maxPeaks - highRangeMaxPks;
	int lowRangeNumPks = middle - pks.begin ();
	int highRangeNumPks = pks.end () - middle;
	if ( highRangeNumPks > highRangeMaxPks ) {
		pks.erase ( middle + highRangeMaxPks, pks.end () );
	}
	if ( lowRangeNumPks > lowRangeMaxPks ) {
		pks.erase ( pks.begin () + lowRangeMaxPks, middle );
	}
	sort ( pks.begin (), pks.end (), sortPeaksByMass () );
	if ( pks.size () < minPeaks ) {
		pks.erase ( pks.begin (), pks.end () );
	}
}
void PeakContainer::retainNPeaksPerDaRangeMSMS ( PeakVector& pks, int n, double daRange, PeakVectorSizeType minPeaks )
{
	stable_sort ( pks.begin (), pks.end (), sortPeaksByMOverZ () );
	double startMOverZ = pks.front ().getMOverZ ();
	double endMOverZ = pks.back ().getMOverZ ();
	double mOverZLimit = startMOverZ;
	PeakVectorIterator start = pks.begin ();
	do {
		mOverZLimit += daRange;
		PeakVectorIterator middle = stable_partition ( start, pks.end (), LessThanMOverZ ( mOverZLimit ) );
		sort ( start, middle, sortPeaksByDescendingIntensity () );
		int lowRangeNumPks = middle - start;
		if ( lowRangeNumPks > n ) {
			pks.erase ( start + n, middle );
			start += n;
		}
		else
			start += lowRangeNumPks;
	} while ( mOverZLimit <= endMOverZ );
	stable_sort ( pks.begin (), pks.end (), sortPeaksByMass () );
	if ( pks.size () < minPeaks ) {
		pks.erase ( pks.begin (), pks.end () );
	}
}
class CheckMatrix {
	double maxMass;
	bool isMatrixPeak ( double mass, double massError );
public:
	CheckMatrix ( double maxMass ) :
		maxMass ( maxMass ) {}
	bool operator () ( const Peak& peak )
	{
		return peak.getMassMinusTol () <= maxMass && peak.getCharge () == 1 && isMatrixPeak ( peak.getMass (), peak.getTolerance () );
	}
};
bool CheckMatrix::isMatrixPeak ( double mass, double massError )
{
	double deltaM = 0.00048 * mass;			// Mann 43rd ASMS p639
	double width = 0.19 + ( 0.0001 * mass );
	double halfWidth = width / 2.0;
	double minOffset = deltaM - halfWidth - massError;
	double maxOffset = deltaM + halfWidth + massError;
	double offset = mass - floor ( mass );
	if ( offset > minOffset && offset < maxOffset ) return false;
	else return true;
}
void PeakContainer::removeMatrix ( PeakVector& pks, double maxMass )
{
	int i = remove_if ( pks.begin (), pks.end (), CheckMatrix ( maxMass ) ) - pks.begin ();
	pks.erase ( pks.begin () + i, pks.end () );
}
class RemoveIsotopeDistribution {
	double curMZ;
	double step;
	double curTol;
	int nSteps;
public:
	RemoveIsotopeDistribution ( double curMZ, int curZ, double curTol ) :
		curMZ ( curMZ ),
		step ( 1.0 / curZ ),
		curTol ( curTol ),
		nSteps ( 0 ) {}
	bool operator () ( const DataFilePeak& pk )
	{
		if ( nSteps > 5 ) return false;
		if ( pk.getMOverZ () > curMZ + curTol ) {
			curMZ += step;
			nSteps++;
		}
		if ( pk.getMOverZ () < curMZ - curTol ) {
			return false;
		}
		else if ( pk.getMOverZ () > curMZ + curTol ) {
			return false;
		}
		else {
			return true;
		}
	}
};
void PeakContainer::removeFTPeaks ( DataFilePeakVector& pks, const Peak* parentPeak, const Tolerance* tol )
{
	double mz = parentPeak->getMOverZ ();
	int z = parentPeak->getCharge ();
	removeIsotopeDistribution ( pks, mz, z, tol->getTolerance ( mz ) );
	removeIsotopeDistribution ( pks, mz / 2.0, z * 2, tol->getTolerance ( mz / 2.0 ) );
	int reducedZ = z - 1;
	double reducedMZ = mPlusHToMOverZ ( parentPeak->getMass (), reducedZ, true );
	removeIsotopeDistribution ( pks, reducedMZ, reducedZ, tol->getTolerance ( reducedMZ ) );
}
void PeakContainer::removeECDorETDSideChainPeaks ( DataFilePeakVector& pks, const Peak* parentPeak, const Tolerance* parTol, const Tolerance* tol )
{
	double mz = parentPeak->getMOverZ ();
	double parMass = parentPeak->getMass ();
	int z = parentPeak->getCharge ();
	static DoubleVector losses = instInf->getLossMasses ();
	static double maxLoss = *(min_element ( losses.begin (), losses.end () ));
	for ( int ch = 2 ; ch <= z ; ch++ ) {
		for ( int i = 0 ; i < losses.size () ; i++ ) {
			double m = mPlusHToMOverZ ( parMass + losses [i], z, ch, true );
			double mTol = genMax ( tol->getTolerance ( m ), parTol->getTolerance ( m ) );
			if ( ch < z || losses [i] >= -2.0 ) {
				removeMOverZRange ( pks, m - mTol, m + mTol );
			}
		}
	}
	double mOZLow = mPlusHToMOverZ ( parentPeak->getMass () + maxLoss, z, 1, true );
	double lowTol = genMax ( tol->getTolerance ( mOZLow ), parTol->getTolerance ( mOZLow ) );
	removeAboveMOverZ ( pks, mOZLow - lowTol );
}
void PeakContainer::removeIsotopeDistribution ( DataFilePeakVector& pks, double mz, int z, double tol )
{
	int i = remove_if ( pks.begin (), pks.end (), RemoveIsotopeDistribution ( mz, z, tol ) ) - pks.begin ();
	pks.erase ( pks.begin () + i, pks.end () );
}
void PeakContainer::setToleranceValue ( double val )
{
	const_cast <Tolerance*> (tolerance)->setValue ( val );
	for ( PeakPtrVectorSizeType i = 0 ; i < peaks.size () ; i++ ) {
		Peak* p = peaks [i];
		double t = tolerance->getTolerance ( p->getMOverZ (), p->getCharge () );
		p->setTolerance ( t );
	}
}
void PeakContainer::calibrate ( double gradient, double offset )
{
	for ( PeakPtrVectorSizeType i = 0 ; i < peaks.size () ; i++ ) {
		peaks [i]->calibrate ( tolerance, gradient, offset );
	}
}
int PeakContainer::getMaxCharge () const
{
	int maxCharge = 1;
	for ( PeakPtrVectorSizeType i = 0 ; i < peaks.size () ; i++ ) {
		maxCharge = genMax( maxCharge, peaks [i]->getCharge () );
	}
	return maxCharge;
}
double PeakContainer::getMaxIntensity () const
{
	double maxIntensity = 0.0;
	for ( PeakPtrVectorSizeType i = 0 ; i < peaks.size () ; i++ ) {
		maxIntensity = genMax( maxIntensity, peaks [i]->getIntensity () );
	}
	return maxIntensity;
}
double PeakContainer::getMinMassMinusTol () const
{
	if ( !peaks.empty () ) {
		return peaks [0]->getMassMinusTol ();
	}
	return 0.0;
}
double PeakContainer::getMaxMassPlusTol () const
{
	if ( !peaks.empty () ) {
		return peaks [peaks.size ()-1]->getMassPlusTol ();
	}
	return 0.0;
}
bool PeakContainer::getNonUnitChargeData () const
{
	for ( PeakPtrVectorSizeType i = 0 ; i < peaks.size () ; i++ ) {
		if ( peaks [i]->getCharge () != 1 ) return true;
	}
	return false;
}
void PeakContainer::printPeaksHTML ( ostream& os, const PeakPrecision& pp ) const
{
	int numPeaks = peaks.size ();

	if ( numPeaks > 0 ) {
		ExpandableJavascriptBlock ejb ( "Fragment Ions" );
		ejb.printHeader ( os );
			if ( numPeaks == 1 ) os << "Ion";
			else os << "<b>" << numPeaks << "</b> Ions";
			os << " used in search: ";
			os << "<b>";
			for ( int i = 0 ; i < numPeaks ; i++ ) {
				peaks [i]->printMOverZHTML ( os, pp.getMassDecimalPlaces () );
				if ( i != numPeaks - 1 ) os << ", ";
			}
			tolerance->print ( os );
			os << "</b>";
			os << " <br />" << endl;
		ejb.printFooter ( os );
	}
}
void PeakContainer::printPeaksXML ( ostream& os, const string& label, const PeakPrecision& pp ) const
{
	os << "<" + label + "_peaks>" << endl;
	int numPeaks = peaks.size ();
	ParameterList::printXML ( os, "num_peaks", numPeaks );

	for ( int i = 0 ; i < numPeaks ; i++ ) {
		peaks [i]->printXML ( os, pp );
	}
	//tolerance->print ( os );
	os << "</" + label + "_peaks>";
	os << endl;
}
void PeakContainer::putCGI ( ostream& os, const string& name, const Peak* parentPeak ) const
{
	os << "data=";
	parentPeak->putParentCGI ( os );
	for ( PeakPtrVectorSizeType i = 0 ; i < peaks.size () ; i++ ) {
		peaks [i]->putCGI ( os );
	}
	os << "&";
	int maxCharge = getMaxCharge ();
	printCGI ( os, "max_charge", maxCharge );
	printCGIString ( os, "data_format", string ( "PP M/Z Intensity Charge" ) );
	tolerance->putCGI ( os, name );
}
void PeakContainer::putCGI ( ostream& os, const string& name, const BoolDeque& peakUsed ) const
{
	os << "data=";
	for ( PeakPtrVectorSizeType i = 0 ; i < peaks.size () ; i++ ) {
		if ( peakUsed [i] == false ) {
			peaks [i]->putCGI ( os );
		}
	}
	os << "&";
	printCGIString ( os, "data_format", string ( "PP M/Z Intensity Charge" ) );
	tolerance->putCGI ( os, name );
}
void PeakContainer::putCGI ( ostream& os, const string& name ) const
{
	os << "data=";
	for ( PeakPtrVectorSizeType i = 0 ; i < peaks.size () ; i++ ) {
		peaks [i]->putCGI ( os );
	}
	os << "&";
	printCGIString ( os, "data_format", string ( "PP M/Z Intensity Charge" ) );
	tolerance->putCGI ( os, name );
}
void PeakContainer::putHiddenFormEntry ( ostream& os, const string& name, const BoolDeque& peakUsed ) const
{
	os << "<input type=\"hidden\" name=\"data\" value=\"";
	for ( PeakPtrVectorSizeType i = 0 ; i < peaks.size () ; i++ ) {
		if ( peakUsed [i] == false ) {
			peaks [i]->putHiddenFormEntry ( os );
		}
	}
	os << "\" />" << endl;

	printHTMLFORMHidden ( os, "data_format", string ( "PP M/Z Intensity Charge" ) );
	tolerance->putHiddenFormEntry ( os, name );
}
void PeakContainer::putHiddenFormJavascriptEntry ( ostream& os, const string& name, const BoolDeque& peakUsed, const string& str ) const
{
	static int idx = 1;
	startJavascript ( os );
	os << "function printUnmatchedPeaks" << str << " () {" << endl;
		os << "\tdocument.writeln ( \"";
		os << "<input";
		os << " ";
		os << "type=\\\"hidden\\\"";
		os << " ";
		os << "name=\\\"data\\\"";
		os << " ";
		os << "value=\\\"";
		PeakPtrVectorSizeType i = 0;
		for ( ; i < peaks.size () ; i++ ) {
			if ( peakUsed [i] == false ) {				// Write out the first number
				peaks [i]->putHiddenFormJavascriptEntry ( os );
				i++;
				break;
			}
		}
		os << "\" );" << endl;
		for ( ; i < peaks.size () ; i++ ) {
			if ( peakUsed [i] == false ) {
				os << "\tdocument.writeln ( \"";
				peaks [i]->putHiddenFormJavascriptEntry ( os );
				os << "\" );" << endl;
			}
		}
		os << "\tdocument.writeln ( \"\\\" />\" );" << endl;

		printHTMLFORMJavascriptHidden ( os, "data_format", string ( "PP M/Z Intensity Charge" ) );
		tolerance->putHiddenFormJavascriptEntry ( os, name );
	os << "}" << endl;
	endJavascript ( os );
	idx++;
}
void PeakContainer::putHiddenFormEntry ( ostream& os, const string& name ) const
{
	os << "<input type=\"hidden\" name=\"data\" value=\"";
	for ( PeakPtrVectorSizeType i = 0 ; i < peaks.size () ; i++ ) {
		peaks [i]->putHiddenFormEntry ( os );
	}
	os << "\" />" << endl;

	printHTMLFORMHidden ( os, "data_format", string ( "PP M/Z Intensity Charge" ) );
	tolerance->putHiddenFormEntry ( os, name );
}
void PeakContainer::setMassAverage () const
{
	for ( PeakPtrVectorSizeType i = 0 ; i < peaks.size () ; i++ ) peaks [i]->setMassAverage ();
}
Peak::Peak ( double mOverZ, double tolerance, int charge, double intensity, double adductMass, int averageToMonoConvertFlag ) :
	mOverZ ( mOverZ ),
	tolerance ( tolerance ),
	charge ( charge ),
	intensity ( intensity ),
	adductMass ( adductMass )
{
	mass = mOverZToMPlusH ( mOverZ, charge, adductMass );

	if ( averageToMonoConvertFlag ) {
		averageMass = mass;
		mass *= AVERAGE_TO_MONOISOTOPIC_FACTOR;
	}
	else averageMass = 0.0;
}
Peak::Peak ( double mOverZ, double tolerance, int charge, double intensity, double mass, double adductMass, double averageMass, const std::string& specNumber ) :
	mOverZ ( mOverZ ),
	tolerance ( tolerance ),
	charge ( charge ),
	intensity ( intensity ),
	mass ( mass ),
	adductMass ( adductMass ),
	averageMass ( averageMass ),
	specNumber ( specNumber )
{
}
void Peak::calibrate ( const Tolerance* t, double gradient, double offset )
{
	mOverZ -= t->getCorrection ( mOverZ, gradient, offset );
	mass = mOverZToMPlusH ( mOverZ, charge, adductMass );
	if ( averageMass != 0.0 ) {
		averageMass = mass;
		mass *= AVERAGE_TO_MONOISOTOPIC_FACTOR;
	}
}
void Peak::putCGI ( ostream& os ) const
{
	genPrint ( os, mOverZ, 4 );
	os << "+";
	genPrint ( os, intensity, 1 );
	if ( charge != 1 ) {
		os << "+" << charge;
	}
	os << "%0D%0A";
}
void Peak::putParentCGI ( ostream& os ) const
{
	genPrint ( os, mOverZ, 4 );
	if ( charge != 1 ) {
		os << "+" << charge;
	}
	os << "%0D%0A";
}
void Peak::putHiddenFormEntry ( ostream& os ) const
{
	genPrint ( os, mOverZ, 4 );
	os << " ";
	genPrint ( os, intensity, 1 );
	if ( charge != 1 ) {
		os << " " << charge;
	}
	os << "\n";
}
void Peak::putHiddenFormJavascriptEntry ( ostream& os ) const
{
	genPrint ( os, mOverZ, 4 );
	os << " ";
	genPrint ( os, intensity, 1 );
	if ( charge != 1 ) {
		os << " " << charge;
	}
}
void Peak::printMOverZHTML ( ostream& os, int precision ) const
{
	genPrint ( os, mOverZ, precision );
	print_charge_superscript ( os, charge );
}
void Peak::printMOverZDelimited ( ostream& os, int precision ) const
{
	delimitedCell ( os, mOverZ, precision );
	delimitedCell ( os, charge );
}
void Peak::printHTML ( ostream& os, const string& label, const PeakPrecision& pp ) const
{
	os << label << ": <b>";
	printMOverZHTML ( os, pp.getMassDecimalPlaces () );
	if ( charge == 1 ) {
		os << " (+/- ";
		genPrint ( os, tolerance, pp.getErrorSigFig () );
		os << " Da)";
	}
	os << "</b>";
	os << "<br />";
	os << endl;
	if ( charge != 1 ) {
		os << "M+H equivalent: <b>";
		genPrint ( os, mass, pp.getMassDecimalPlaces () );
		os << " (+/- ";
		genPrint ( os, tolerance, pp.getErrorSigFig () );
		os << " Da)";
		os << "</b>";
		os << "<br />";
		os << endl;
	}
}
void Peak::printXML ( ostream& os, const string& label, const PeakPrecision& pp ) const
{
	ParameterList::printDoubleXMLFixed ( os, label + "_m_over_z", mOverZ, pp.getMassDecimalPlaces () );
	ParameterList::printXML ( os, label + "_charge", charge );
	if ( charge == 1 ) {
		ParameterList::printDoubleXMLSigFig ( os, label + "_tolerance", tolerance, pp.getErrorSigFig () );
	}
	else {
		ParameterList::printDoubleXMLFixed ( os, label + "_m_plus_h", mass, pp.getMassDecimalPlaces () );
		ParameterList::printDoubleXMLSigFig ( os, label + "_tolerance", tolerance, pp.getErrorSigFig () );
	}
}
void Peak::printXML ( ostream& os, const PeakPrecision& pp ) const
{
	printXML ( os, mOverZ, charge, intensity, pp.getMassDecimalPlaces (), pp.getIntensitySigFig () );
}
void Peak::printDelimited ( ostream& os, const PeakPrecision& pp ) const
{
	printDelimited ( os, mOverZ, charge, intensity, pp.getMassDecimalPlaces (), pp.getIntensitySigFig () );
}
void Peak::printXML ( ostream& os, double mOverZ, int charge, double intensity, int mPrecision, int iPrecision )
{
	os << "<mzi>";
	genPrint ( os, mOverZ, mPrecision );
	os << ",";
	os << charge;
	os << ",";
	genPrintSigFig ( os, intensity, iPrecision );
	os << "</mzi>";
	os << endl;
}
void Peak::printDelimited ( ostream& os, double mOverZ, int charge, double intensity, int mPrecision, int iPrecision )
{
	delimitedCell ( os, mOverZ, mPrecision );
	delimitedCell ( os, charge );
	delimitedCellSigFig ( os, intensity, iPrecision );
}
PeakMatchContext::PeakMatchContext ( const Tolerance* tolerance, const PeakPrecision& pp, bool nonUnitChargeData ) :
	tolerance ( tolerance ),
	peakPrecision ( pp ),
	nonUnitChargeData ( nonUnitChargeData )
{
}

PeakMatch::PeakMatch ( const Peak* dataPeak, double matchedMass ) :
	dataPeak ( dataPeak ),
	matchedMass ( matchedMass )
{
}
void PeakMatch::printHeaderHTML ( ostream& os, const PeakMatchContext& pmc )
{
	tableHeader ( os, "m/z<br />Submitted" );
	if ( pmc.getNonUnitChargeData () ) {
		tableHeaderStart ( os );
			os << mh_plus_html << "<br />Equivalent";
		tableHeaderEnd ( os );
	}
	tableHeaderStart ( os );
		os << mh_plus_html << "<br />Matched";
	tableHeaderEnd ( os );
	tableHeader ( os, "Intensity" );
	tableHeaderStart ( os );
		os << "Delta <br />" << pmc.getToleranceUnits ();
	tableHeaderEnd ( os );
	if ( dataPeak->getSpecNumber () != "" ) {
		tableHeader ( os, "Data<br />Set" );
	}
}
void PeakMatch::printHeaderDelimited ( ostream& os, const string& errUnits, bool multi )
{
	delimitedHeader ( os, "m/z Submitted" );
	delimitedHeader ( os, "Charge" );
	delimitedHeader ( os, "MH+ Equivalent" );
	delimitedHeader ( os, "MH+ Matched" );
	delimitedHeader ( os, "Intensity" );
	delimitedHeader ( os, "Delta " + errUnits );
	if ( multi ) delimitedHeader ( os, "Data Set" );
}
double PeakMatch::getError ( const PeakMatchContext& pmc ) const
{
	return pmc.getError ( dataPeak->getMass (), matchedMass, dataPeak->getCharge () );
}
void PeakMatch::printHTML ( ostream& os, const PeakMatchContext& pmc, bool printLine ) const
{
	bool nonUnitChargeData = pmc.getNonUnitChargeData ();
	if ( printLine ) {
		PeakPrecision pp = pmc.getPeakPrecision ();
		double measuredMass = dataPeak->getMass ();

		tableDataStart ( os, "", "right" );
			dataPeak->printMOverZHTML ( os, pp.getMassDecimalPlaces () );
			os << endl;
		tableDataEnd ( os );
		if ( nonUnitChargeData ) {
			tableDataStart ( os, "", "right" );
				genPrint ( os, measuredMass, pp.getMassDecimalPlaces () );
				os << endl;
			tableDataEnd ( os );
		}
		tableDataStart ( os, "", "right" );
			genPrint ( os, matchedMass, pp.getMassDecimalPlaces () );
			os << endl;
		tableDataEnd ( os );
		tableDataStart ( os, "", "right" );
			genPrint ( os, dataPeak->getIntensity (), 1 );
			os << endl;
		tableDataEnd ( os );
		tableDataStart ( os, "", "right", true );
			genPrintSigFig ( os, pmc.getError ( measuredMass, matchedMass, dataPeak->getCharge () ), pp.getErrorSigFig () );
			os << endl;
		tableDataEnd ( os );
		if ( dataPeak->getSpecNumber () != "" )	tableCell ( os, dataPeak->getSpecNumber (), true );
	}
	else {
		tableEmptyNCells ( os, 4 );
		if ( nonUnitChargeData ) tableEmptyCell ( os );
		if ( dataPeak->getSpecNumber () != "" ) {
			tableEmptyCell ( os );  // spec ID
		}
	}
}
void PeakMatch::printDelimited ( ostream& os, const PeakMatchContext& pmc, bool printLine ) const
{
	if ( printLine ) {
		PeakPrecision pp = pmc.getPeakPrecision ();
		double massPrecision = pp.getMassDecimalPlaces ();
		double measuredMass = dataPeak->getMass ();

		dataPeak->printMOverZDelimited ( os, massPrecision );
		delimitedCell ( os, measuredMass, massPrecision );
		delimitedCell ( os, matchedMass, massPrecision );
		delimitedCell ( os, dataPeak->getIntensity (), 1 );
		delimitedCellSigFig ( os, pmc.getError ( measuredMass, matchedMass, dataPeak->getCharge () ), pp.getErrorSigFig () );
	}
	else delimitedEmptyNCells ( os, 6 );

	if ( !dataPeak->getSpecNumber ().empty () ) {
		if ( printLine )	delimitedCell ( os, dataPeak->getSpecNumber () );
		else				delimitedEmptyCell ( os );
	}
}
void PeakMatch::printXML ( ostream& os, const PeakMatchContext& pmc ) const
{
	bool nonUnitChargeData = pmc.getNonUnitChargeData ();
	PeakPrecision pp = pmc.getPeakPrecision ();
	double measuredMass = dataPeak->getMass ();

	dataPeak->printXML ( os, pp );
	ParameterList::printDoubleXMLFixed ( os, "peptide_mw", matchedMass, pp.getMassDecimalPlaces () );
	ParameterList::printDoubleXMLSigFig ( os, "error", pmc.getError ( measuredMass, matchedMass, dataPeak->getCharge () ), pp.getErrorSigFig () );
}
void PeakMatch::printTagHitHTML ( ostream& os, const PeakMatchContext& pmc, bool zeroErr ) const
{
	PeakPrecision pp = pmc.getPeakPrecision ();
	double measuredMass = dataPeak->getMass ();

	tableHeaderStart ( os, "", "right" );
		genPrint ( os, matchedMass, pp.getMassDecimalPlaces () );
		os << endl;
	tableHeaderEnd ( os );
	tableHeaderStart ( os, "", "right", true );
		if ( zeroErr )
			genPrint ( os, 0.0, 1 );
		else
			genPrintSigFig ( os, pmc.getError ( measuredMass, matchedMass, dataPeak->getCharge () ), pp.getErrorSigFig () );
		os << endl;
	tableHeaderEnd ( os );
}
