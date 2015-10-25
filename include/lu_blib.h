 /*****************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_blib.h                                                     *
*                                                                             *
*  Created    : July 18th 2013                                                *
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

#ifndef __lu_blib_h
#define __lu_blib_h

#include <lu_data.h>
#include <lu_sqlite.h>

class BlibRetentionTimesEntry {
	int refSpectraID;
	int redundantRefSpectraID;
	int spectrumSourceID;
	double retentionTime;
	int bestSpectrum;
};

class BlibRefSpectraEntryValue {
	std::string peptideSeq;
	VectorPairIntDouble vpid;
	double precursorMZ;
	std::string prevAA;
	std::string nextAA;
	int copies;
	int fileID;
	int numPeaks;
	double retentionTime;
	std::string specIDinFile;
	std::string specID;
	double pValue;
	double eValue;
	int scoreType;

	DoubleVector rts;
	IntVector fileIDs;
	IntVector indicies;
	int bestSpectrumIdx;
public:
	BlibRefSpectraEntryValue ( const std::string& peptideSeq, const VectorPairIntDouble& vpid, double precursorMZ, const std::string& prevAA, const std::string& nextAA, int copies, int numPeaks, int fileID, double retentionTime, int index, const std::string& specIDinFile, const std::string& specID, double pValue, double eValue );
	BlibRefSpectraEntryValue () {}

	std::string getPeptideSeq () const { return peptideSeq; }
	VectorPairIntDouble getVpid () const { return vpid; }
	double getPrecursorMZ () const { return precursorMZ; }
	std::string getPrevAA () const { return prevAA; }
	std::string getNextAA () const { return nextAA; }
	int getCopies () const { return copies; }
	int getNumPeaks () const { return numPeaks; }
	int getFileID () const { return fileID; }
	double getRetentionTime () const { return retentionTime; }
	std::string getSpecIDinFile () const { return specIDinFile; }
	std::string getSpecID () const { return specID; }
	double getPValue () const { return pValue; }
	double getEValue () const { return eValue; }
	int getScoreType () const { return scoreType; }

	IntVector getFileIDs () const { return fileIDs; }
	DoubleVector getRTs () const { return rts; }
	IntVector getIndicies () const { return indicies; }
	int getBestSpectrumIdx () const { return bestSpectrumIdx; }

	void incrementCopies () { copies++; }
	void update ( double precursorMZ, const std::string& prevAA, const std::string& nextAA, int numPeaks, int fileID, double retentionTime, int index, const std::string& specIDinFile, const std::string& specID, double pValue, double eValue );
	void update ( int fileID, double retentionTime, int index );
};

typedef std::map <PairStringInt, BlibRefSpectraEntryValue> MapBlibRefSpectra;
typedef MapBlibRefSpectra::iterator MapBlibRefSpectraIterator;
typedef MapBlibRefSpectra::const_iterator MapBlibRefSpectraConstIterator;

class Blib : public PPSQLite {
protected:
	double* pM;
	float* pI;
	unsigned int pSiz;
	unsigned char* compMZ;
	unsigned long compMZLen;
	unsigned char* compInt;
	unsigned long compIntLen;
	void resizePMPI ( int numPeaks );
	void resizeCompMZ ( unsigned long compMZSiz );
	void resizeCompInt ( unsigned long compIntSiz );
public:
	Blib ();
	virtual ~Blib ();
	void updateQuery ( const std::string& table, const PairStringString& setClause, const PairStringString& whereClause );
};

class BlibWrite : public Blib {
	void populatePMPI ( const DataFilePeakVector& dfpv );
	void updateRefSpectraNumPeaks ( int spectraID, int numPeaks );
public:
	BlibWrite ( const std::string& name, bool append = false );
	~BlibWrite ();
	void create ();
	void insertLibInfo ( int numSpecs, const std::string& server, bool redundant, const std::string& name );
	void insertModifications ( int refSpectraID, const BlibRefSpectraEntryValue& val );
	void insertModifications ( int refSpectraID, const VectorPairIntDouble& vpid );
	void insertRefSpectra ( const PairStringInt& key, const BlibRefSpectraEntryValue& val );
	void insertRetentionTimes ( int refSpectraID, const BlibRefSpectraEntryValue& val );
	void insertRetentionTimes ( int refSpectraID, int redundantRefSpectraID, int spectrumSourceID, double retentionTime, int bestSpectrum );
	void insertRefSpectraPeaks ( int spectraID, const DataFilePeakVector& dfpv );
	void insertRefSpectraPeaks ( int spectraID, unsigned char* mzComp, int mzCompSiz, unsigned char* intComp, int intCompSiz );
	void insertScoreTypes ();
	void insertSpectrumSourceFiles ( const std::string& file );
};

struct sqlite3_stmt;

class BlibRead : public Blib {
	StringVector peakListFractionNames;
	StringVector peakListCentroidFileNames;
	StringVector peakListCentroidFilePaths;
	void writeSpectrum ( std::ostream& os, int rsi, int numPeaks ) const;
	int readSpectrumRow ( sqlite3_stmt* pStmt );
	void readResults ( StringVectorVector& header, StringVectorVector& rows );
	void readSpectra ( const std::string& path );
public:
	BlibRead ( const std::string& name );
	~BlibRead ();
	void read ( const std::string& peakListPath, StringVectorVector& header, StringVectorVector& rows );
	StringVector getPeakListFractionNames () const { return peakListFractionNames; }
	StringVector getPeakListCentroidFileNames () const { return peakListCentroidFileNames; }
	StringVector getPeakListCentroidFilePaths () const { return peakListCentroidFilePaths; }
};

#endif /* ! __lu_blib_h */
