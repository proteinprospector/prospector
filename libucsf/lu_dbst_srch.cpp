/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_dbst_srch.cpp                                              *
*                                                                             *
*  Created    : September 5th 2001                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <iomanip>
#include <lu_histogram.h>
#include <lp_frame.h>
#include <lu_dbst_srch.h>
#include <lu_r_plot.h>
#include <lu_app_gr.h>
#include <lu_param_list.h>
#include <lu_table.h>
using std::vector;
using std::ostream;
using std::sort;
using std::streamsize;
using std::setprecision;
using std::endl;
using std::make_pair;

class AminoAcidStats {
	unsigned int a;
	unsigned int b;
	unsigned int c;
	unsigned int d;
	unsigned int e;
	unsigned int f;
	unsigned int g;
	unsigned int h;
	unsigned int i;
	unsigned int j;
	unsigned int k;
	unsigned int l;
	unsigned int m;
	unsigned int n;
	unsigned int o;
	unsigned int p;
	unsigned int q;
	unsigned int r;
	unsigned int s;
	unsigned int t;
	unsigned int u;
	unsigned int v;
	unsigned int w;
	unsigned int x;
	unsigned int y;
	unsigned int z;
	unsigned int total;
	void accumulate ( FastaServer* fs, const IntVector& indicies, int frame );
public:
	AminoAcidStats ( vector <FastaServer*>& fs, const DBStatParameters& params, int frame );
	AminoAcidStats ( FastaServer* fs, const IntVector& indicies, int frame );
	void printHTML ( std::ostream& os ) const;
	void printHTML ( std::ostream& os, char aa, unsigned int number, unsigned int total ) const;
};

class DBStatHistogram : public Histogram {
public:
	DBStatHistogram () : Histogram () {}
	void drawGraph ( ostream& os ) const;
};
void DBStatHistogram::drawGraph ( ostream& os ) const
{
	compute ();

	if ( xyData.size () ) {
		SpectrumGraph s1 ( "dbstat_hist.par.txt" );
		FileGraphData gd ( xyData, true );
		s1.drawGraph ( os, gd, false );
		os << "<p />" << endl;
	}
}

static DBStatHistogram hist;
static DoubleVector masses;
static double numMasses = 0.0;
static unsigned int MAX_NUM_MASSES = 3000000;
DBStatParameters::DBStatParameters ( const ParameterList* params ) :
	MSSearchParameters ( params ),
	showAAStatistics	( params->getBoolValue ( "show_aa_statistics", false ) ),
	minHistogramMass	( params->getDoubleValue ( "min_histogram_mass", 600.0 ) ),
	maxHistogramMass	( params->getDoubleValue ( "max_histogram_mass", 15000.0 ) ),
	bandwidth			( params->getDoubleValue ( "density_bandwidth", 1.0 ) )
{
	PreSearchInfo::setReportTaxonomy ();
}
DBStatSearch::DBStatSearch ( const DBStatParameters& params ) :
	DBSearch ( params ),
	aaStats ( 0 ),
	bandwidth ( params.getBandwidth () ),
	maxMissedCleavages ( params.getMissedCleavages () ),
	outputFlag ( false )
{
	int frame = 1;
	if ( params.getShowAAStatistics () ) aaStats = new AminoAcidStats ( fs, params, frame );
	init_fasta_enzyme_function ( params.getEnzyme () );
	maxNumAA = 0;
	largestNumEnzymeFragments = 0;
	totalNumFragments = 0.0;
	averageProteinMW = 0.0;
	totalAA = 0.0;
	unsigned int totalIndicies = 0;
	for ( int i = 0 ; i < fs.size () ; i++ ) {
		const IntVector& indicies = params.getIndicies ( i );
		if ( !indicies.empty () ) outputFlag = true;
		double minHistMass = params.getMinHistogramMass ();
		double maxHistMass = params.getMaxHistogramMass ();
		if ( !indicies.empty () ) {
			FrameIterator* fi = new FrameIterator ( fs [i], indicies, make_pair ( frame, frame ), params.getTempOverride () );
			char* openReadingFrame;
			while ( ( openReadingFrame = fi->getNextFrame () ) != NULL ) {
				IntVector& cleavageIndex = enzyme_fragmenter ( openReadingFrame );
				int numEnzymeFragments = cleavageIndex.size ();
				calculateDigestHistogram ( openReadingFrame, cleavageIndex, minHistMass, maxHistMass );
				int numEnzymeFragmentWithMissCleave = getNumFragments ( numEnzymeFragments, maxMissedCleavages );
				totalNumFragments += numEnzymeFragmentWithMissCleave;
				int numAA = strlen ( openReadingFrame );
				totalAA += numAA;
				ProteinMW pmw ( openReadingFrame );
				double proteinMW = pmw.getMass ();
				averageProteinMW += proteinMW;
				if ( numAA > maxNumAA ) {
					maxNumAA = numAA;
					longestProteinDBIndex = i;
					longestProteinIndex = fi->getEntry ();
					longestProteinMW = proteinMW;
					longestProteinOrfNumber = fi->getFrame () + 1;
				}
				if ( numEnzymeFragmentWithMissCleave > largestNumEnzymeFragments ) {
					largestNumEnzymeFragments = numEnzymeFragmentWithMissCleave;
					largestNumEnzymeFragmentsDBIndex = i;
					largestNumEnzymeFragmentsIndex = fi->getEntry ();
					largestNumEnzymeFragmentsOrfNumber = fi->getFrame () + 1;
				}
			}
			delete fi;
			FrameIterator::resetElapsedTime ( 1 );
			totalIndicies += indicies.size ();
		}
	}
	if ( totalIndicies ) averageProteinMW /= totalIndicies;
}
DBStatSearch::~DBStatSearch ()
{
	delete aaStats;
}
int DBStatSearch::getNumFragments ( int numEnzymeFragments, int maxMissedCleavages )
{
	int num = 0;
	for ( int i = 0 ; i <= maxMissedCleavages ; i++ ) {
		int n = numEnzymeFragments - i;
		if ( n > 0 ) num += n;
		else break;
	}
	return num;
}
void DBStatSearch::calculateDigestHistogram ( const char* frame, const IntVector& cleavageIndex, double minHistMass, double maxHistMass )
{
	DoubleVector& enzymeFragmentMassArray = get_cleaved_masses ( frame, cleavageIndex );
	int numEnzymeFragments = cleavageIndex.size ();
	int missedCleavageLimit = maxMissedCleavages;
	for ( int i = 0 ; i < numEnzymeFragments ; i++ ) {
		double mol_wt = terminal_wt;
		for ( int j = i ; j <= missedCleavageLimit ; j++ ) {
			if ( j >= numEnzymeFragments ) break;
			mol_wt += enzymeFragmentMassArray [j];
			if ( cnbr_digest && j == missedCleavageLimit && frame [cleavageIndex [j]] == 'M' ) {
				mol_wt += cnbr_homoserine_lactone_mod;
			}
			if ( mol_wt >= minHistMass && mol_wt <= maxHistMass ) {
				if ( hist.size () < MAX_NUM_MASSES ) hist.add ( mol_wt );
				if ( masses.size () < MAX_NUM_MASSES ) masses.push_back ( mol_wt );
				numMasses += 1.0;
			}
		}
		missedCleavageLimit++;
	}
}
void DBStatSearch::printHTMLHits ( ostream& os )
{
	if ( outputFlag ) {
		os << "Longest protein is index number ";
		os << "<b>";
		if ( fs.size () > 1 ) os << fs [longestProteinDBIndex]->getFileName () << ": ";
		os << longestProteinIndex;
		os << "</b> ";
		if ( fs [longestProteinDBIndex]->getDNADatabase () ) {
			os << "(Open Reading Frame <b>" << longestProteinOrfNumber << "</b>) ";
		}
		os << "which has <b>" << maxNumAA << "</b> amino acids";
		os << endl;
		os << "<br />" << endl;

		os << "Mass of longest protein: <b>";
		genPrint ( os, longestProteinMW, 1 );
		os << "</b> Da";
		os << endl;
		os << "<br />" << endl;
		os << "<br />" << endl;

		os << "Largest number of enzyme fragments is index number ";
		os << "<b>";
		if ( fs.size () > 1 ) os << fs [largestNumEnzymeFragmentsDBIndex]->getFileName () << ": ";
		os << largestNumEnzymeFragmentsIndex;
		os << "</b> ";
		if ( fs [largestNumEnzymeFragmentsDBIndex]->getDNADatabase () ) {
			os << "(Open Reading Frame <b>" << largestNumEnzymeFragmentsOrfNumber << "</b>) ";
		}
		os << "which has <b>" << largestNumEnzymeFragments << "</b> fragments";
		os << endl;
		os << "<br />" << endl;

		os << "Total number of enzyme fragments is <b>";
		genPrint ( os, totalNumFragments, 0 );
		os << "</b>";
		os << endl;
		os << "<br />" << endl;

		os << "Average protein MW is <b>";
		genPrint ( os, averageProteinMW, 0 );
		os << "</b> Da";
		os << endl;
		os << "<br />" << endl;

		os << "Total number of amino acids is <b>";
		genPrint ( os, totalAA, 0 );
		os << "</b>";
		os << endl;
		os << "<br />" << endl;

		os << "Total number of enzyme fragments in histogram range is <b>";
		genPrint ( os, numMasses, 0 );
		os << "</b>";
		os << endl;
		os << "<br />" << endl;
		os << "<br />" << endl;

		if ( aaStats ) aaStats->printHTML ( os );
	}
	hist.drawGraph ( os );
	if ( masses.size () <= 300000 ) {
		if ( RPlot::getRFlag () ) {
			RPlot rplot ( "peptideDensity.R" );
			GenOFStream ofs ( rplot.getDataFileFullPath () );
			genPrint ( ofs, bandwidth, 4 );
			ofs << endl;
			sort ( masses.begin (), masses.end () );
			for ( DoubleVectorSizeType i = 0 ; i < masses.size () ; i++ ) {
				genPrint ( ofs, masses [i], 5 );
				ofs << endl;
			}
			ofs.close ();
			os << "<p>" << endl;
				rplot.printImage ( os, ".png" );
			os << "</p>" << endl;
			os << "<p>" << endl;
				rplot.printImage ( os, ".pdf" );
			os << "</p>" << endl;
		}
	}
}
void DBStatSearch::printXMLHits ( ostream& os ) const
{
	if ( outputFlag ) {
		os << "<results>";
		os << "<longest_protein>";
		ParameterList::printXML ( os, "longest_protein_database", fs [longestProteinDBIndex]->getFileName () );
		ParameterList::printXML ( os, "longest_protein_index", longestProteinIndex );
		if ( fs [longestProteinDBIndex]->getDNADatabase () ) {
			ParameterList::printXML ( os, "orf_number", longestProteinOrfNumber );
		}
		ParameterList::printXML ( os, "num_aa", maxNumAA );
		ParameterList::printDoubleXMLFixed ( os, "mass", longestProteinMW, 1 );
		os << "</longest_protein>";

		os << "<digest>";
		ParameterList::printXML ( os, "largest_num_enzyme_fragments_database", fs [largestNumEnzymeFragmentsDBIndex]->getFileName () );
		ParameterList::printXML ( os, "largest_num_enzyme_fragments_index", largestNumEnzymeFragmentsIndex );
		if ( fs [largestNumEnzymeFragmentsDBIndex]->getDNADatabase () ) {
			ParameterList::printXML ( os, "longest_orf_number", largestNumEnzymeFragmentsOrfNumber );
		}
		ParameterList::printXML ( os, "longest_num_fragments", largestNumEnzymeFragments );
		ParameterList::printDoubleXMLFixed ( os, "total_fragments", totalNumFragments, 0 );
		os << "</digest>";

		ParameterList::printDoubleXMLFixed ( os, "average_protein_mw", averageProteinMW, 0 );
		ParameterList::printDoubleXMLFixed ( os, "total_amino_acids", totalAA, 0 );
		os << "</results>";
	}
}
AminoAcidStats::AminoAcidStats ( vector <FastaServer*>& fs, const DBStatParameters& params, int frame ) :
	a(0),b(0),c(0),d(0),e(0),f(0),g(0),h(0),i(0),j(0),k(0),l(0),m(0),
	n(0),o(0),p(0),q(0),r(0),s(0),t(0),u(0),v(0),w(0),x(0),y(0),z(0),
	total(0)
{
	for ( int i = 0 ; i < fs.size () ; i++ ) {
		accumulate ( fs [i], params.getIndicies ( i ), frame );
	}
}
AminoAcidStats::AminoAcidStats ( FastaServer* fs, const IntVector& indicies, int frame ) :
	a(0),b(0),c(0),d(0),e(0),f(0),g(0),h(0),i(0),j(0),k(0),l(0),m(0),
	n(0),o(0),p(0),q(0),r(0),s(0),t(0),u(0),v(0),w(0),x(0),y(0),z(0),
	total(0)
{
	accumulate ( fs, indicies, frame );
}
void AminoAcidStats::accumulate ( FastaServer* fs, const IntVector& indicies, int frame )
{
	int numIndicies = indicies.size ();

	for ( int ii = 0 ; ii < numIndicies ; ii++ ) {
		char* protein = fs->get_fasta_protein ( indicies [ii], frame );
		for ( int jj = 0 ; protein [jj] != '\0' ; jj++ ) {
			total++;
			switch ( protein [jj] ) {
				case 'A':
					a++;
					break;
				case 'B':
					b++;
					break;
				case 'C':
					c++;
					break;
				case 'D':
					d++;
					break;
				case 'E':
					e++;
					break;
				case 'F':
					f++;
					break;
				case 'G':
					g++;
					break;
				case 'H':
					h++;
					break;
				case 'I':
					i++;
					break;
				case 'J':
					j++;
					break;
				case 'K':
					k++;
					break;
				case 'L':
					l++;
					break;
				case 'M':
					m++;
					break;
				case 'N':
					n++;
					break;
				case 'O':
					o++;
					break;
				case 'P':
					p++;
					break;
				case 'Q':
					q++;
					break;
				case 'R':
					r++;
					break;
				case 'S':
					s++;
					break;
				case 'T':
					t++;
					break;
				case 'U':
					u++;
					break;
				case 'V':
					v++;
					break;
				case 'W':
					w++;
					break;
				case 'X':
					x++;
					break;
				case 'Y':
					y++;
					break;
				case 'Z':
					z++;
					break;
			}
		}
	}
}
void AminoAcidStats::printHTML ( ostream& os ) const
{
	os << "<table border=\"border\" cellspacing=\"3\">" << endl;
		tableRowStart (	os );
			tableHeader ( os, "Amino Acid" );
			tableHeader ( os, "Total in database" );
			tableHeader ( os, "Percent of total" );
		tableRowEnd (	os );
		printHTML ( os, 'A', a, total );
		printHTML ( os, 'B', b, total );
		printHTML ( os, 'C', c, total );
		printHTML ( os, 'D', d, total );
		printHTML ( os, 'E', e, total );
		printHTML ( os, 'F', f, total );
		printHTML ( os, 'G', g, total );
		printHTML ( os, 'H', h, total );
		printHTML ( os, 'I', i, total );
		printHTML ( os, 'J', j, total );
		printHTML ( os, 'K', k, total );
		printHTML ( os, 'L', l, total );
		printHTML ( os, 'M', m, total );
		printHTML ( os, 'N', n, total );
		printHTML ( os, 'O', o, total );
		printHTML ( os, 'P', p, total );
		printHTML ( os, 'Q', q, total );
		printHTML ( os, 'R', r, total );
		printHTML ( os, 'S', s, total );
		printHTML ( os, 'T', t, total );
		printHTML ( os, 'U', u, total );
		printHTML ( os, 'V', v, total );
		printHTML ( os, 'W', w, total );
		printHTML ( os, 'X', x, total );
		printHTML ( os, 'Y', y, total );
		printHTML ( os, 'Z', z, total );
		tableRowStart (	os );
			tableCell ( os, "total" );
			tableCell ( os, total );
			tableEmptyCell ( os );
		tableRowEnd (	os );
	os << "</table>";
}
void AminoAcidStats::printHTML ( ostream& os, char aa, unsigned int number, unsigned int total ) const
{
	tableRowStart (	os );
		tableCell ( os, aa );
		tableCell ( os, number );
		tableDataStart ( os );
			double dnumber = (double) number;
			double dtotal = (double) total;
			double percent = dnumber * 100.0 / dtotal;
			streamsize prec = os.precision ();
			os << setprecision ( 4 ) << percent << setprecision ( prec ) << endl;
		tableDataEnd ( os );
	tableRowEnd (	os );
}
