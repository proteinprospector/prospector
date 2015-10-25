/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_t2d.cpp                                                    *
*                                                                             *
*  Created    : January 2nd 2004                                              *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2004-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef RAW_DATA
#include <zlib.h>
#include <lg_string.h>
#include <lgen_file.h>
#include <lu_param_list.h>
#include <lu_t2d.h>
using std::ios;
using std::string;
using std::ostream;
using std::endl;
using std::getline;

enum PBSpectrumFileDataFormats {
    PBSpectrumFileUnsigned16 = 1,
    PBSpectrumFileSigned32 = 2,
    PBSpectrumFileFloat32 = 3,
    PBSpectrumFileDouble64 = 4
};

enum PBSpectrumFileCompressionTypes {
    PBSpectrumFileCompressionNone = 1,
    PBSpectrumFileCompressionZLIB = 2
};

enum PBSpectrumFileMode {
    PBSpectrumFileModeLinear = 1,
    PBSpectrumFileModeReflector = 2,
    PBSpectrumFileModeMSMS = 3,
    PBSpectrumFileModePSD = 4
};

enum PBSpectrumFileEquationType {
    PBSpectrumFileEquationType1   = 1,
    PBSpectrumFileEquationType2   = 2,
    PBSpectrumFileEquationType3   = 3,
    PBSpectrumFileEquationType4   = 4
};

class PBSpectrumFileCalibration {
    unsigned short equationType;
    unsigned short nConstants;
	double c0;
	double c1;
	double c2;
	double c3;
	double c4;
	double c5;
	double c6;
	double c7;
	double c8;
	double c9;
	double sqc1c3;
public:
	void read ( GenIFStream& ifs );
	void readFromFile ( const std::string& filename, unsigned short eType );
	void printHTML ( std::ostream& os ) const;
	double timeToMass1 ( double time );
	double timeToMass2 ( double time );
	double timeToMass3 ( double time );
	double timeToMass4 ( double time );
	unsigned short getEquationType () const { return equationType; };
};
double PBSpectrumFileCalibration::timeToMass1 ( double time )
{
	double r1 = sqrt(c1) * (time - c0);
	double q = 1 - 4 * c2 * r1;
	double r2 = (1 - sqrt(q)) / (2 * c2);
	return r2 * r2;
}
double PBSpectrumFileCalibration::timeToMass2 ( double time )
{
	double q = sqrt(c1) * (time - c0);
	double r = q * (1 + c2 * q) * (1 + c3 * q * q);
	return r * r;
}
double PBSpectrumFileCalibration::timeToMass3 ( double time )
{
	double r1 = sqc1c3 * (time - c0);
	double q = 1 - 4 * c2 * r1;
	double r2 = (1 - sqrt(q))/(2 * c2);
	return c3 * r2 * r2;
}
double PBSpectrumFileCalibration::timeToMass4 ( double time )
{
	double r1 = sqc1c3 * (time - c0);
	double q = 1 - 4 * c2 * r1;
	double r2 = (1 - sqrt(q))/(2 * c2);
	double r2sq = r2 * r2;
	double r3 = (1 - r2sq);
	return c3 * r2sq * (1 + c4 * r3 * r3 * r3);
}

void PBSpectrumFileCalibration::read ( GenIFStream& ifs )
{
	ifs.read ( (char*) &equationType, sizeof (unsigned short) );
	ifs.read ( (char*) &nConstants, sizeof (unsigned short) );
	ifs.read ( (char*) &c0, sizeof (double) );
	ifs.read ( (char*) &c1, sizeof (double) );
	ifs.read ( (char*) &c2, sizeof (double) );
	ifs.read ( (char*) &c3, sizeof (double) );
	ifs.read ( (char*) &c4, sizeof (double) );
	ifs.read ( (char*) &c5, sizeof (double) );
	ifs.read ( (char*) &c6, sizeof (double) );
	ifs.read ( (char*) &c7, sizeof (double) );
	ifs.read ( (char*) &c8, sizeof (double) );
	ifs.read ( (char*) &c9, sizeof (double) );
	if ( equationType == 3 || equationType == 4 ) {
		sqc1c3 = sqrt(c1 / c3);
	}
}
void PBSpectrumFileCalibration::readFromFile ( const string& filename, unsigned short eType )
{
	equationType = eType;
	GenIFStream ifs ( filename );
	string line;
	while ( getline ( ifs, line ) ) {
		if ( isPrefix ( line, "IRCCurrentCalibrationCoefficients" ) ) {
			DoubleVector c ( 10, 0.0 );
			int start = line.find ( '=' ) + 1;
			nConstants = 0;
			int end;
			do {
				end = line.find ( ',', start );
				c [nConstants++] = atof ( line.substr ( start, end ).c_str () );
				start = end+1;
			} while ( end != string::npos );
			c0 = c [0];
			c1 = c [1];
			c2 = c [2];
			c3 = c [3];
			c4 = c [4];
			c5 = c [5];
			c6 = c [6];
			c7 = c [7];
			c8 = c [8];
			c9 = c [9];
			break;
		}
	}
	if ( equationType == 3 || equationType == 4 ) {
		sqc1c3 = sqrt(c1 / c3);
	}
}
void PBSpectrumFileCalibration::printHTML ( ostream& os ) const
{
	ParameterList::printHTML ( os, "Equation Type", equationType );
	ParameterList::printHTML ( os, "NConstants", nConstants );
	ParameterList::printHTML ( os, "Constant 0", c0 );
	ParameterList::printHTML ( os, "Constant 1", c1 );
	ParameterList::printHTML ( os, "Constant 2", c2 );
	ParameterList::printHTML ( os, "Constant 3", c3 );
	ParameterList::printHTML ( os, "Constant 4", c4 );
	ParameterList::printHTML ( os, "Constant 5", c5 );
	ParameterList::printHTML ( os, "Constant 6", c6 );
	ParameterList::printHTML ( os, "Constant 7", c7 );
	ParameterList::printHTML ( os, "Constant 8", c8 );
	ParameterList::printHTML ( os, "Constant 9", c9 );
}
enum PBSpectrumFileSpectrumFlags {
    PBSpectrumFlagSaturated = 0x00000001
};

struct PBSpectrumFileLocation {
    double  xMin;
    double  yMin;
    double  xMax;
    double  yMax;
};

struct PBSpectrumFileLaserIntensity {
    double  lMin;
    double  lMax;
};

/*
File Header

Spectrum Header 1
Spectrum 1

Spectrum Header N
Spectrum N

Spectrum File Indicies
File Trailer
*/
/*
Example octal dump
0000000    ñ   ¶   ¤  \0   4 022   Z   ¥  \0  \0 001  \0  \0  \0 001  \0
0000020    @ 001 024  \0   W   h   o   C   a   r   e   s  \0  \0  \0  \0
0000040   \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0
0000060   \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0
0000100   \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0
0000120   \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0
0000140   \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0
0000160   \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0
0000200   \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0
0000220   \0  \0  \0  \0 001  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0  \0
0000240   \0  \0  \0  \0
0000000  f1b6 a400 3412 5aa5 0000 0100 0000 0100
0000020  4001 1400 5768 6f43 6172 6573 0000 0000
0000040  0000 0000 0000 0000 0000 0000 0000 0000
0000060  0000 0000 0000 0000 0000 0000 0000 0000
0000100  0000 0000 0000 0000 0000 0000 0000 0000
0000120  0000 0000 0000 0000 0000 0000 0000 0000
0000140  0000 0000 0000 0000 0000 0000 0000 0000
0000160  0000 0000 0000 0000 0000 0000 0000 0000
0000200  0000 0000 0000 0000 0000 0000 0000 0000
0000220  0000 0000 0100 0000 0000 0000 0000 0000
0000240  0000 0000
*/

class SpectrumHeader {
	unsigned short checksum;
	unsigned short unused;
	PBSpectrumFileSignatures signature;
	DBKEY spectrumID;
	unsigned long instrumentID;
	unsigned short operatingMode;
	unsigned short dataFormat;
	unsigned short compressionType;
	unsigned short dataChecksum;
	unsigned long sizePoints;
	unsigned long sizeBytes;
	double startTime;
	double incrementTime;
	double timeStamp;
	double totalONCount;
	double basePeakTimeS;
	double basePeakIntensity;
	unsigned long totalShots;
	unsigned long totalAccumulations;
	PBSpectrumFileCalibration defaultCalibration;
	PBSpectrumFileCalibration acquisitionCalibration;
	PBSpectrumFileLocation locationBounds;
	PBSpectrumFileLaserIntensity intensityRange;
	unsigned long flags;
public:
	SpectrumHeader ( GenIFStream& ifs );
	unsigned long getSizeBytes () const { return sizeBytes; };
	unsigned long getSizePoints () const { return sizePoints; };
	double getStartTime () const { return startTime; };
	double getIncrementTime () const { return incrementTime; };
	PBSpectrumFileCalibration getCalibration () const { return defaultCalibration; };
	unsigned short getEquationType () const { return defaultCalibration.getEquationType (); };
	void printHTML ( std::ostream& os ) const;
	static unsigned short size ()
	{
		unsigned short s = 0;
		s += 2 * sizeof (unsigned short);
		s += sizeof (PBSpectrumFileSignatures);
		s += sizeof (DBKEY);
		s += sizeof (unsigned long);
		s += 4 * sizeof (unsigned short);
		s += 2 * sizeof (unsigned long);
		s += 6 * sizeof (double);
		s += 2 * sizeof (unsigned long);
		s += 2 * sizeof (PBSpectrumFileCalibration);
		s += sizeof (PBSpectrumFileLocation);
		s += sizeof (PBSpectrumFileLaserIntensity);
		s += sizeof (unsigned long);
		return s;
	}
};

T2DHeader::T2DHeader ()
{
	//checksum = ;
	fileHeaderSize = T2DHeader::size ();
	signature = PBSpectrumFileHeaderSignature;
	softwareMajorVersion = 0;
	softwareMinorVersion = 1;
	fileMajorVersion = 0;
	fileMinorVersion = 1;
	spectrumHeaderSize = SpectrumHeader::size ();
	fileTrailerSize = T2DTrailer::size ();
	strcpy ( instrumentDBProvider.s, "WhoCares" );
	fileID.dlKey = 0;
}
T2DHeader::T2DHeader ( GenIFStream& ifs )
{
	ifs.read ( (char*) &checksum,				sizeof (unsigned short) );
	ifs.read ( (char*) &fileHeaderSize,			sizeof (unsigned short) );
	ifs.read ( (char*) &signature,				sizeof (PBSpectrumFileSignatures) );
	ifs.read ( (char*) &softwareMajorVersion,	sizeof (unsigned short) );
	ifs.read ( (char*) &softwareMinorVersion,	sizeof (unsigned short) );
	ifs.read ( (char*) &fileMajorVersion,		sizeof (unsigned short) );
	ifs.read ( (char*) &fileMinorVersion,		sizeof (unsigned short) );
	ifs.read ( (char*) &spectrumHeaderSize,		sizeof (unsigned short) );
	ifs.read ( (char*) &fileTrailerSize,		sizeof (unsigned short) );
	ifs.read ( (char*) &instrumentDBProvider,	sizeof (DBProvider) );
	ifs.read ( (char*) &fileID,					sizeof (DBKEY) );
}
void T2DHeader::write ( ostream& os ) const
{
	os.write ( (char*) &checksum,				sizeof (unsigned short) );
	os.write ( (char*) &fileHeaderSize,			sizeof (unsigned short) );
	os.write ( (char*) &signature,				sizeof (PBSpectrumFileSignatures) );
	os.write ( (char*) &softwareMajorVersion,	sizeof (unsigned short) );
	os.write ( (char*) &softwareMinorVersion,	sizeof (unsigned short) );
	os.write ( (char*) &fileMajorVersion,		sizeof (unsigned short) );
	os.write ( (char*) &fileMinorVersion,		sizeof (unsigned short) );
	os.write ( (char*) &spectrumHeaderSize,		sizeof (unsigned short) );
	os.write ( (char*) &fileTrailerSize,		sizeof (unsigned short) );
	os.write ( (char*) &instrumentDBProvider,	sizeof (DBProvider) );
	os.write ( (char*) &fileID,					sizeof (DBKEY) );
}
void T2DHeader::printHTML ( ostream& os ) const
{
	ParameterList::printHTML ( os, "Checksum", checksum );
	ParameterList::printHTML ( os, "File Header Size", fileHeaderSize );
	ParameterList::printHTML ( os, "Signature", signature == PBSpectrumFileHeaderSignature );
	ParameterList::printHTML ( os, "Software Major Version", softwareMajorVersion );
	ParameterList::printHTML ( os, "Software Minor Version", softwareMinorVersion );
	ParameterList::printHTML ( os, "File Major Version", fileMajorVersion );
	ParameterList::printHTML ( os, "File Minor Version", fileMinorVersion );
	ParameterList::printHTML ( os, "Spectrum Header Size", spectrumHeaderSize );
	ParameterList::printHTML ( os, "File Trailer Size", fileTrailerSize );
	ParameterList::printHTML ( os, "Instrument DB Provider", instrumentDBProvider.s );
//	ParameterList::printHTML ( os, "File ID", fileID.dlKey );
	os << "<br />" << endl;
}


SpectrumHeader::SpectrumHeader ( GenIFStream& ifs )
{
	ifs.read ( (char*) &checksum,			sizeof (unsigned short) );
	ifs.read ( (char*) &unused,				sizeof (unsigned short) );
	ifs.read ( (char*) &signature,			sizeof (PBSpectrumFileSignatures) );
	ifs.read ( (char*) &spectrumID,			sizeof (DBKEY) );
	ifs.read ( (char*) &instrumentID,		sizeof (unsigned long) );
	ifs.read ( (char*) &operatingMode,		sizeof (unsigned short) );
	ifs.read ( (char*) &dataFormat,			sizeof (unsigned short) );
	ifs.read ( (char*) &compressionType,	sizeof (unsigned short) );
	ifs.read ( (char*) &dataChecksum,		sizeof (unsigned short) );
	ifs.read ( (char*) &sizePoints,			sizeof (unsigned long) );
	ifs.read ( (char*) &sizeBytes,			sizeof (unsigned long) );
	ifs.read ( (char*) &startTime,			sizeof (double) );
	ifs.read ( (char*) &incrementTime,		sizeof (double) );
	ifs.read ( (char*) &timeStamp,			sizeof (double) );
	ifs.read ( (char*) &totalONCount,		sizeof (double) );
	ifs.read ( (char*) &basePeakTimeS,		sizeof (double) );
	ifs.read ( (char*) &basePeakIntensity,	sizeof (double) );
	ifs.read ( (char*) &totalShots,			sizeof (unsigned long) );
	ifs.read ( (char*) &totalAccumulations,	sizeof (unsigned long) );
	defaultCalibration.read ( ifs );
	acquisitionCalibration.read ( ifs );
	ifs.read ( (char*) &locationBounds,	sizeof (PBSpectrumFileLocation) );
	ifs.read ( (char*) &intensityRange,	sizeof (PBSpectrumFileLaserIntensity) );
	ifs.read ( (char*) &flags,	sizeof (unsigned long) );
}
void SpectrumHeader::printHTML ( ostream& os ) const
{
	ParameterList::printHTML ( os, "Checksum", checksum );
	ParameterList::printHTML ( os, "Unused", unused );
	ParameterList::printHTML ( os, "Signature", signature == PBSpectrumFileSpectrumSignature );
//	ParameterList::printHTML ( os, "Spectrum ID", spectrumID );
	ParameterList::printHTML ( os, "Instrument ID", instrumentID );
	ParameterList::printHTML ( os, "Operating Mode", operatingMode );
	ParameterList::printHTML ( os, "Data Format", dataFormat );
	ParameterList::printHTML ( os, "Compression Type", compressionType );
	ParameterList::printHTML ( os, "Data Checksum", dataChecksum );
	ParameterList::printHTML ( os, "Size Points", sizePoints );
	ParameterList::printHTML ( os, "Size Bytes", sizeBytes );
	ParameterList::printHTML ( os, "Start Time", startTime );
	ParameterList::printHTML ( os, "Increment Time", incrementTime );
	ParameterList::printHTML ( os, "Time Stamp", timeStamp );
	ParameterList::printHTML ( os, "Total ON Count", totalONCount );
	ParameterList::printHTML ( os, "Base Peak Time S", basePeakTimeS );
	ParameterList::printHTML ( os, "Base Peak Intensity", basePeakIntensity );
	ParameterList::printHTML ( os, "Total Shots", totalShots );
	ParameterList::printHTML ( os, "Total Accumulations", totalAccumulations );
	os << "Default Calibration<br />" << endl;
	defaultCalibration.printHTML ( os );
	os << "Acquisition Calibration<br />" << endl;
	acquisitionCalibration.printHTML ( os );

	os << "Location Bounds<br />" << endl;
	ParameterList::printHTML ( os, "X Min", locationBounds.xMin );
	ParameterList::printHTML ( os, "Y Min", locationBounds.yMin );
	ParameterList::printHTML ( os, "X Max", locationBounds.xMax );
	ParameterList::printHTML ( os, "Y Max", locationBounds.yMax );
	os << "Intensity Range<br />" << endl;
	ParameterList::printHTML ( os, "L Min", intensityRange.lMin );
	ParameterList::printHTML ( os, "L Max", intensityRange.lMax );
	os << "<br />" << endl;
	ParameterList::printHTML ( os, "Flags", flags );
	os << "<br />" << endl;
}

T2DTrailer::T2DTrailer ( unsigned long indexOffset )
{
	//checksum = ;
	signature = PBSpectrumFileTrailerSignature;
	//indexChecksum = ;
	numSpectra = 1;
	spectrumIndexOffset = indexOffset;
}
T2DTrailer::T2DTrailer ( GenIFStream& ifs )
{
	ifs.read ( (char*) &checksum,			sizeof (unsigned short) );
	ifs.read ( (char*) &signature,			sizeof (PBSpectrumFileSignatures) );
	ifs.read ( (char*) &indexChecksum,		sizeof (unsigned short) );
	ifs.read ( (char*) &numSpectra,			sizeof (unsigned long) );
	ifs.read ( (char*) &spectrumIndexOffset,sizeof (GENINT64) );
}
void T2DTrailer::write ( ostream& os ) const
{
	os.write ( (char*) &checksum,			sizeof (unsigned short) );
	os.write ( (char*) &signature,			sizeof (PBSpectrumFileSignatures) );
	os.write ( (char*) &indexChecksum,		sizeof (unsigned short) );
	os.write ( (char*) &numSpectra,			sizeof (unsigned long) );
	os.write ( (char*) &spectrumIndexOffset,sizeof (GENINT64) );
}
void T2DTrailer::printHTML ( ostream& os ) const
{
	ParameterList::printHTML ( os, "Checksum", checksum );
	ParameterList::printHTML ( os, "Signature", signature == PBSpectrumFileTrailerSignature );
	ParameterList::printHTML ( os, "Index Checksum", indexChecksum );
	ParameterList::printHTML ( os, "Num Spectra", numSpectra );
//	ParameterList::printHTML ( os, "Spectrum Index Offset", spectrumIndexOffset );
	os << "<br />" << endl;
}

T2DFile::T2DFile ( const string& filename, double startMass, double endMass )
{
	GenIFStream ifs ( filename + string ( ".t2d" ), ios::binary );
	T2DHeader t2dHeader ( ifs );
	SpectrumHeader spectrumHeader ( ifs );
	unsigned long numBytes = spectrumHeader.getSizeBytes ();
	unsigned char* data = new unsigned char [numBytes];
	ifs.read ( (char*) data, numBytes * sizeof (unsigned char) );
	unsigned long numPoints = spectrumHeader.getSizePoints ();
	FloatVector intensity ( numPoints );
	unsigned long destLen = numPoints * sizeof (float);
	uncompress ( (unsigned char*)&intensity[0], &destLen, data, numBytes );
	unsigned short equationType = spectrumHeader.getEquationType ();
	string calFilename = filename + string ( "-1.cal" );
	PBSpectrumFileCalibration cal;
	if ( genFileExists ( calFilename ) ) {
		cal.readFromFile ( calFilename, equationType );
	}
	else {
		cal = spectrumHeader.getCalibration ();
	}
	bool massLimits = ( startMass != 0.0 );
	double currentTime = spectrumHeader.getStartTime ();
	double incrementTime = spectrumHeader.getIncrementTime ();
	double m;
	int j = 0;
	for ( unsigned long i = 0 ; i < numPoints ; i++, currentTime += incrementTime ) {
		if ( equationType == 1 ) m = cal.timeToMass1 ( currentTime );
		else if ( equationType == 2 ) m = cal.timeToMass2 ( currentTime );
		else if ( equationType == 3 ) m = cal.timeToMass3 ( currentTime );
		else if ( equationType == 4 ) m = cal.timeToMass4 ( currentTime );
		if ( massLimits ) {
			if ( m >= startMass ) {
				if ( m > endMass ) break;
				xyData.add ( m, intensity [i] );
			}
		}
		else xyData.add ( m, intensity [i] );
	}
}
void T2DFile::printASCII ( ostream& os ) const
{
	for ( xyData.first () ; xyData.isDone () ; xyData.next () ) {
		genPrint ( os, xyData.x (), 4 );
		os << " ";
		genPrintSigFig ( os, xyData.y (), 3 );
		os << endl;
	}
}
int getT2DMSRunNumberIndex ( int msmsRunNumber, const IntVector& msRunNumbers )
{
	if ( msRunNumbers.empty () ) {
		return -1;						// This is the error condition
	}
	else {
		int diff = msmsRunNumber - msRunNumbers [0];
		int idx = 0;
		for ( int i = 1 ; i < msRunNumbers.size () ; i++ ) {
			int currDiff = msmsRunNumber - msRunNumbers [i];
			if ( ( currDiff == 0 ) || ( currDiff < 0 && diff < 0 && currDiff > diff ) || ( currDiff > 0 && ( diff < 0 || currDiff < diff ) ) ) { 
				diff = currDiff;
				idx = i;
			}
		}
		return idx;
	}
}
#endif
