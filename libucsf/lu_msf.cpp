/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_msf.cpp                                                    *
*                                                                             *
*  Created    : October 15th 2013                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2013-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <iostream>
#include <string>
#include <sqlite3.h>
#include <lg_io.h>
#include <lg_string.h>
#include <lgen_define.h>
#include <lgen_file.h>
#include <lgen_uncompress.h>
#include <lgen_xml.h>
#include <lu_getfil.h>
#include <lu_html.h>
#include <lu_mass_conv.h>
#include <lu_msf.h>
#include <lz_app.h>
using std::string;
using std::ios_base;
using std::cout;
using std::endl;
using std::ostream;
using std::pair;
using std::make_pair;

class PPExpatMSFSpectrum : public PPExpat {
	bool peakCentroidFlag;
	bool headerFlag;
	bool monoisotopicPeakCentroidsFlag;
	double precursorMZ;
	int precursorCharge;
	DoubleVector mass;
	DoubleVector intensity;
	double rt;
	int scanNum;
	int fileID;
	void startElement ( const char* name, const char** attributes );
	void characterDataHandler ( const char* str, int len );
	void endElement ( const char* name );
public:
	PPExpatMSFSpectrum ();
	~PPExpatMSFSpectrum ();
	void writeSpectrum ( const string& filePath ) const;
	double getPrecursorMZ () const { return precursorMZ; }
	int getPrecursorCharge () const { return precursorCharge; }
	int getScanNum () const { return scanNum; }
	double getRT () const { return rt; }
	int getFileID () const { return fileID; }
};
PPExpatMSFSpectrum::PPExpatMSFSpectrum () :
	peakCentroidFlag ( false ),
	headerFlag ( false ),
	monoisotopicPeakCentroidsFlag ( false )
{
}
PPExpatMSFSpectrum::~PPExpatMSFSpectrum () {}

void PPExpatMSFSpectrum::startElement ( const char* name, const char** attributes )
{
	if ( peakCentroidFlag ) {
		if ( !strcmp ( name, "Peak" ) ) {
			double x;
			double y;
			getAttributeValue ( attributes, "X", x );
			getAttributeValue ( attributes, "Y", y );
			mass.push_back ( x );
			intensity.push_back ( y );
		}
	}
	else {
		if ( headerFlag ) {
			if ( !strcmp ( name, "SpectrumIdentifier" ) ) {
				getAttributeValue ( attributes, "ScanNumber", scanNum );
				getAttributeValue ( attributes, "RetentionTime", rt );
				getAttributeValue ( attributes, "FileID", fileID );
			}
		}
		else {
			if ( monoisotopicPeakCentroidsFlag ) {
				if ( !strcmp ( name, "Peak" ) ) {
					getAttributeValue ( attributes, "X", precursorMZ );
				}
			}
			else {
				if ( !strcmp ( name, "Header" ) ) {
					headerFlag = true;
				}
				else if ( !strcmp ( name, "PeakCentroids" ) ) {
					peakCentroidFlag = true;
				}
				else if ( !strcmp ( name, "MonoisotopicPeakCentroids" ) ) {
					monoisotopicPeakCentroidsFlag = true;
				}
				else if ( !strcmp ( name, "PrecursorInfo" ) ) {
					getAttributeValue ( attributes, "Charge", precursorCharge );
				}
			}
		}
	}
}
void PPExpatMSFSpectrum::characterDataHandler ( const char* str, int len )
{
}
void PPExpatMSFSpectrum::endElement ( const char* name )
{
	if ( peakCentroidFlag ) {
		if ( !strcmp ( name, "PeakCentroids" ) ) {
			peakCentroidFlag = false;
		}
	}
	else if ( headerFlag ) {
		if ( !strcmp ( name, "Header" ) ) {
			headerFlag = false;
		}
	}
	else if ( monoisotopicPeakCentroidsFlag ) {
		if ( !strcmp ( name, "MonoisotopicPeakCentroids" ) ) {
			monoisotopicPeakCentroidsFlag = false;
		}
	}
}
void PPExpatMSFSpectrum::writeSpectrum ( const string& filePath ) const
{
	GenOFStream os ( filePath, std::ios_base::out | std::ios_base::app );
	os << "BEGIN IONS" << endl;
	os << "TITLE=Scan ";
	os << scanNum;
	os << " ";
	os << "(rt=";
	genPrint ( os, rt, 3 );
	os << ")";
	os << " ";
	os << "[Prospector Created]";
	os << endl;

	os << "PEPMASS=";
	genPrint ( os, precursorMZ, 4 );
	os << endl;
	os << "CHARGE=" << precursorCharge << "+" << endl;

	int numPeaks = mass.size ();
	for ( DoubleVectorSizeType i = 0 ; i < numPeaks ; i++ ) {
		genPrint ( os, mass [i], 4 );
		os << " ";
		genPrintSigFig ( os, intensity [i], 3 );
		os << endl;
	}
	os << "END IONS" << endl;
}

MapIntToPairStringBool MSFRead::scoreInfo;

SetInt MSFRead::setSID;
MapIntToVectorPairIntString MSFRead::peptideScores;

IntVector MSFRead::peptideID;
IntVector MSFRead::spectrumID;
StringVector MSFRead::sequence;

StringVector MSFRead::modifications;
DoubleVector MSFRead::precursorMZ;
IntVector MSFRead::precursorCharge;
StringVector MSFRead::fractions;
IntVector MSFRead::scan;
DoubleVector MSFRead::rt;
MapIntToInt MSFRead::uniqueSpectrumID;

MapIntToVectorPairIntInt MSFRead::mods;

StringVector MSFRead::modificationName;

MapIntToString MSFRead::fileNames;

MapCharToString MSFRead::constMods;
StringVector MSFRead::msfConstModsStr;

int MSFRead::scoreInfoCallback ( void* NotUsed, int argc, char** argv, char** azColName )
{
	int sid = atoi ( argv [0] );
	string fn = argv [1];
	bool ims = atoi ( argv [2] ) == 0 ? false : true;
	scoreInfo [sid] = make_pair ( fn, ims );
	return 0;
}
int MSFRead::peptideScoresCallback ( void* NotUsed, int argc, char** argv, char** azColName )
{
	int pid = atoi ( argv [0] );
	int sid = atoi ( argv [1] );
	string ps =  argv [2];
	setSID.insert ( sid );
	peptideScores [pid].push_back ( make_pair ( sid, ps ) );
	return 0;
}
int MSFRead::peptidesCallback ( void* NotUsed, int argc, char** argv, char** azColName )
{
	peptideID.push_back		( atoi ( argv [0] ) );
	spectrumID.push_back	( atoi ( argv [1] ) );
	sequence.push_back		( argv [2] );
	return 0;
}
int MSFRead::spectrumHeadersCallback ( void* NotUsed, int argc, char** argv, char** azColName )
{
	int sid = atoi ( argv [0] );
	int usid = atoi ( argv [1] );
	uniqueSpectrumID [sid] = usid;
	return 0;
}
int MSFRead::peptidesAminoAcidModificationsCallback ( void* NotUsed, int argc, char** argv, char** azColName )
{
	int pid = atoi ( argv [0] );
	int aamid = atoi ( argv [1] );
	int pos = atoi ( argv [2] );
	mods [pid].push_back ( make_pair ( pos+1, aamid ) );
	return 0;
}
int MSFRead::peptidesTerminalModificationsCallback ( void* NotUsed, int argc, char** argv, char** azColName )
{
	int pid = atoi ( argv [0] );
	int aamid = atoi ( argv [1] );
	int pos = atoi ( argv [2] );
	if ( pos == 1 || pos == 3 ) {							// 1=peptide, 3=protein
		mods [pid].push_back ( make_pair ( -3, aamid ) );	// N-term
	}
	if ( pos == 2 || pos == 4 ) {							// 2=peptide, 4=protein
		mods [pid].push_back ( make_pair ( -2, aamid ) );	// C-term
	}
	return 0;
}
int MSFRead::aminoAcidModificationsCallback ( void* NotUsed, int argc, char** argv, char** azColName )
{
	modificationName.push_back ( argv [0] );
	return 0;
}
int MSFRead::fileInfosCallback ( void* NotUsed, int argc, char** argv, char** azColName )
{
	int fileID = atoi ( argv [0] );
	string physicalFileName = genShortFilenameFromPath ( argv [1] );
	fileNames [fileID] = physicalFileName;
	return 0;
}
int MSFRead::processingNodeParametersCallback ( void* NotUsed, int argc, char** argv, char** azColName )
{
	addConstMod ( argv [0] );
	return 0;
}
MSFRead::MSFRead ( const string& name )
{
	rc = sqlite3_open_v2 ( name.c_str (), &db, SQLITE_OPEN_READONLY, NULL );
}
MSFRead::~MSFRead ()
{
}
void MSFRead::clear ()
{
	scoreInfo.clear ();
	setSID.clear ();
	peptideScores.clear ();
	peptideID.clear ();
	spectrumID.clear ();
	sequence.clear ();
	modifications.clear ();
	precursorMZ.clear ();
	precursorCharge.clear ();
	fractions.clear ();
	scan.clear ();
	rt.clear ();
	uniqueSpectrumID.clear ();
	mods.clear ();
	modificationName.clear ();
	fileNames.clear ();
	constMods.clear ();
	msfConstModsStr.clear ();
}
void MSFRead::read ( const string& peakListPath, StringVectorVector& rows )
{
	clear ();
	genCreateDirectory ( peakListPath );
	readDB ();
	readSpectra ( peakListPath );
	readResults ( rows );
}
void MSFRead::readDB ()
{
	UpdatingJavascriptMessage ujm;
	ujm.writeMessage ( cout, "Reading processing node scores." );
	string sql;
	sql = "SELECT ScoreID, FriendlyName, IsMainScore from ProcessingNodeScores";
	executeSQL ( sql, scoreInfoCallback );

	ujm.writeMessage ( cout, "Reading peptide scores." );
	sql = "SELECT PeptideID, ScoreID, ScoreValue from PeptideScores";
	sql += " ORDER BY PeptideID, ScoreID";
	executeSQL ( sql, peptideScoresCallback );

	ujm.writeMessage ( cout, "Reading peptides." );
	sql = "SELECT PeptideID, SpectrumID, Sequence from Peptides";
	sql += " ORDER BY SpectrumID, PeptideID";
	executeSQL ( sql, peptidesCallback );

	ujm.writeMessage ( cout, "Reading spectrum headers." );
	sql = "SELECT SpectrumID, UniqueSpectrumID from SpectrumHeaders";
	executeSQL ( sql, spectrumHeadersCallback );

	ujm.writeMessage ( cout, "Reading terminal modifications." );
	sql = "SELECT p.PeptideID, a.AminoAcidModificationID, a.PositionType from PeptidesTerminalModifications AS p, AminoAcidModifications AS a";
	sql += " WHERE p.TerminalModificationID = a.AminoAcidModificationID";
	executeSQL ( sql, peptidesTerminalModificationsCallback );

	ujm.writeMessage ( cout, "Reading amino acid modifications." );
	sql = "SELECT PeptideID, AminoAcidModificationID, Position from PeptidesAminoAcidModifications";
	sql += " ORDER BY PeptideID, Position";
	executeSQL ( sql, peptidesAminoAcidModificationsCallback );

	ujm.writeMessage ( cout, "Amino acid modification names." );
	sql = "SELECT ModificationName from AminoAcidModifications";
	sql += " ORDER BY AminoAcidModificationID";
	executeSQL ( sql, aminoAcidModificationsCallback );
	
	ujm.writeMessage ( cout, "File info." );
	sql = "SELECT FileID, PhysicalFileName from FileInfos";
	sql += " ORDER BY FileID";
	executeSQL ( sql, fileInfosCallback );

	ujm.writeMessage ( cout, "Processing node parameters." );
	sql = "SELECT ParameterValue from ProcessingNodeParameters";
	sql += " WHERE substr ( ParameterName, 1, 9 ) = 'StaticMod'";
	executeSQL ( sql, processingNodeParametersCallback );
	ujm.deletePreviousMessage ( cout );
}
void MSFRead::readSpectra ( const string& path )
{
	ZipOpener zo;				// This is required to open the zip files storing the spectra
	beginTransaction ();
	for ( MapIntToStringConstIterator a = fileNames.begin () ; a != fileNames.end () ; a++ ) {
		peakListFractionNames.push_back ( (*a).second );
		peakListCentroidFileNames.push_back ( (*a).second + ".mgf" );
		peakListCentroidFilePaths.push_back ( path + SLASH + peakListCentroidFileNames.back () );
	}
	SetInt idSet;
	int numSpectra = spectrumID.size ();
	string sNumSpectra = gen_itoa ( numSpectra );
	UpdatingJavascriptMessage ujm;
	ujm.writeMessage ( cout, sNumSpectra + " spectra to process." );
	for ( IntVectorSizeType i = 0 ; i < numSpectra ; i++ ) {
		int num = i + 1;
		string sNum = gen_itoa ( num );
		if ( num % 500 == 0 ) ujm.writeMessage ( cout, sNum + "/" + sNumSpectra + " spectra processed." );
		int sid = spectrumID [i];
		PairSetIntIteratorBool flag = idSet.insert ( sid );
		if ( flag.second ) {								// New spectrum
			MapIntToIntConstIterator cur = uniqueSpectrumID.find ( spectrumID [i] );
			int id = -1;
			if ( cur != uniqueSpectrumID.end () ) {
				id = (*cur).second;
			}
			if ( id != -1 ) {
				string sql = "SELECT Spectrum FROM Spectra";
				sql += " WHERE UniqueSpectrumID = '" + gen_itoa ( id ) + "'";

				sqlite3_stmt* pStmt;
				rc = sqlite3_prepare_v2 ( db, sql.c_str (), -1, &pStmt, 0 );

				if ( rc != SQLITE_OK ) {
					error ( "sqlite3_prepare_v2", sql, rc );
				}
				rc = sqlite3_step ( pStmt );
				if ( rc == SQLITE_ROW ) {
					int rsi = readSpectrumRow ( pStmt, path, zo );
				}
				sqlite3_finalize ( pStmt );
			}
			else {
				precursorMZ.push_back ( 0.0 );
				precursorCharge.push_back ( 0 );
				scan.push_back ( 0 );
				rt.push_back ( 0.0 );
				fractions.push_back ( "" );
			}
		}
		else {
			precursorMZ.push_back ( precursorMZ.back () );
			precursorCharge.push_back ( precursorCharge.back () );
			scan.push_back ( scan.back () );
			rt.push_back ( rt.back () );
			fractions.push_back ( fractions.back () );
		}
	}
	ujm.deletePreviousMessage ( cout );
	endTransaction ();
}
int MSFRead::readSpectrumRow ( sqlite3_stmt* pStmt, const string& path, ZipOpener& zo )
{
/*  Old code - uses unzip

	PPTempFile tempFile ( "spec", ".zip" );										// Create a file name for the zipped blob
	string fname = tempFile.getFullPath ();
	GenOFStream os ( fname, ios_base::binary );									// Write out the file
	int len = sqlite3_column_bytes ( pStmt, 0 );								// Currently the blob is not copied as it is written straight to a file.
	unsigned char* spec = (unsigned char*) sqlite3_column_blob ( pStmt, 0 );	// Read the file out as a zipped blob
 	os.write ( (char*) spec, len * sizeof (char) );
	os.close ();
	string dir = genUnzip ( fname, false, true );								// Unzip it
	FileList fList ( dir, "", ".xml", false );									// Parse it for the spectrum
	{
		PPExpatMSFSpectrum ppems;
		ppems.parseXMLFromFile ( dir + SLASH + fList [0] );
		string fract = fileNames [ppems.getFileID ()];
		string fPath = path + SLASH + fract + ".mgf";
		ppems.writeSpectrum ( fPath );

		precursorMZ.push_back ( ppems.getPrecursorMZ () );
		precursorCharge.push_back ( ppems.getPrecursorCharge () );
		scan.push_back ( ppems.getScanNum () );
		rt.push_back ( ppems.getRT () );
		fractions.push_back ( fract );
	}
	genUnlinkDirectory ( dir );
	genUnlink ( fname );
*/
// New code - uses unzip library
	int len = sqlite3_column_bytes ( pStmt, 0 );
	void* spec = (void*)sqlite3_column_blob ( pStmt, 0 );
	zo.getFirstFile ( spec, len );
	PPExpatMSFSpectrum ppems;
	ppems.parseXMLFromString ( zo.getBuf (), zo.getLen (), true );
	string fract = fileNames [ppems.getFileID ()];
	string fPath = path + SLASH + fract + ".mgf";
	ppems.writeSpectrum ( fPath );

	precursorMZ.push_back ( ppems.getPrecursorMZ () );
	precursorCharge.push_back ( ppems.getPrecursorCharge () );
	scan.push_back ( ppems.getScanNum () );
	rt.push_back ( ppems.getRT () );
	fractions.push_back ( fract );
	return 0;
}
void MSFRead::readHeader ( StringVectorVector& header ) const
{
	StringVector hCols;
	hCols.push_back ( "M/Z" );
	hCols.push_back ( "Charge" );
	hCols.push_back ( "DB Peptide" );
	hCols.push_back ( "Modifications" );
	hCols.push_back ( "Fraction" );
	hCols.push_back ( "Retention Time" );
	hCols.push_back ( "Scan" );
	for ( SetIntConstIterator i = setSID.begin () ; i != setSID.end () ; i++ ) {
		hCols.push_back ( scoreInfo [*i].first );
	}
	header.push_back ( hCols );
}
void MSFRead::readResults ( StringVectorVector& rows ) const
{
	bool nTermConst = false;
	bool cTermConst = false;
	SetPairCharString constAA;
	for ( MapCharToStringConstIterator a = constMods.begin () ; a != constMods.end () ; a++ ) {
		char spec = (*a).first;
		if ( spec == 'n' )		nTermConst = true;
		else if ( spec == 'c' ) cTermConst = true;
		else {
			constAA.insert ( make_pair ( spec, (*a).second ) );
		}
	}
	for ( StringVectorSizeType i = 0 ; i < sequence.size () ; i++ ) {
		int pid = peptideID [i];
		MapIntToVectorPairIntIntConstIterator cur = mods.find ( pid );
		if ( cur != mods.end () ) {
			const VectorPairIntInt& vpii = (*cur).second;
			string m;
			int lastPos = -100;
			for ( VectorPairIntIntSizeType j = 0 ; j < vpii.size () ; j++ ) {
				int pos = vpii [j].first;
				if ( pos == -3 ) {
					if ( !nTermConst ) {
						if ( pos == lastPos ) m = m.substr ( 0, m.rfind ( '@' ) ) + '+';
						m += modificationName [(vpii [j].second)-1] + "@N-term;";
					}
				}
				else if ( pos == -2 ) {
					if ( !cTermConst ) {
						if ( pos == lastPos ) m = m.substr ( 0, m.rfind ( '@' ) ) + '+';
						m += modificationName [(vpii [j].second)-1] + "@C-term;";
					}
				}
				else {
					char aa = sequence [i][pos-1];
					string modName = modificationName [(vpii [j].second)-1];
					if ( constAA.find ( make_pair ( aa, modName ) ) == constAA.end () ) {			// Not a constant mod
						if ( pos == lastPos && !m.empty () ) {
							m = m.substr ( 0, m.rfind ( '@' ) ) + '+';			// If the position is the same as the previous position
						}
						m += modName + "@" + gen_itoa ( pos ) + ";";
					}
				}
				lastPos = pos;
			}
			modifications.push_back ( genStrtrimSemiColon ( m ) );
		}
		else {
			modifications.push_back ( "" );
		}
	}
	for ( StringVectorSizeType k = 0 ; k < sequence.size () ; k++ ) {
		StringVector cols;
		cols.push_back ( gen_ftoa ( precursorMZ [k], "%.4f" ) );
		cols.push_back ( gen_itoa ( precursorCharge [k] ) );
		cols.push_back ( sequence [k] );
		cols.push_back ( modifications [k] );
		cols.push_back ( fractions [k] );
		cols.push_back ( gen_ftoa ( rt [k], "%.3f" ) );
		cols.push_back ( gen_itoa ( scan [k] ) );
		const VectorPairIntString& vpis = peptideScores [peptideID [k]];
		VectorPairIntStringSizeType n = 0;
		for ( SetIntConstIterator m = setSID.begin () ; m != setSID.end () ; m++ ) {
			for ( ; n < vpis.size () ; n++ ) {
				if ( *m == vpis [n].first ) {
					cols.push_back ( vpis [n].second );
					break;
				}
			}
		}
		rows.push_back ( cols );
	}
}
void MSFRead::setConstMods ( const MapCharToString& cm )	// Initialises the const mods to what's on the menu 
{															// before reading the file 
	constMods = cm;
}
MapCharToString MSFRead::getConstModMap ()
{
	return constMods;
}
StringVector MSFRead::getMsfConstModsStr ()
{
	return msfConstModsStr;
}
void MSFRead::addConstMod ( const string& cm )
{
	msfConstModsStr.push_back ( cm );
	string::size_type idx1 = cm.rfind ( '(' );
	string::size_type idx2 = cm.rfind ( ')' );
	string spec = cm.substr ( idx1+1, idx2-idx1-1 );
	string label = gen_strtrim ( cm.substr ( 0, idx1-1 ) );
	if ( spec == "N-term" ) constMods ['n'] = label;
	else if ( spec == "C-term" ) constMods ['c'] = label;
	else {
		for ( string::size_type i = 0 ; i < spec.length () ; i++ ) {
			constMods [spec[i]] = label;
		}
	}
}
