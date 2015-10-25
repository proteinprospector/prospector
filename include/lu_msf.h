 /*****************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_msf.h                                                      *
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

#ifndef __lu_msf_h
#define __lu_msf_h

struct sqlite3;
class ZipOpener;

#include <lu_sqlite.h>

class MSFRead : public PPSQLite {
	static MapIntToPairStringBool scoreInfo;

	static SetInt setSID;
	static MapIntToVectorPairIntString peptideScores;

	static IntVector peptideID;
	static IntVector spectrumID;
	static StringVector sequence;

	static StringVector modifications;
	static DoubleVector precursorMZ;
	static IntVector precursorCharge;
	static StringVector fractions;
	static IntVector scan;
	static DoubleVector rt;

	static MapIntToInt uniqueSpectrumID;

	static MapIntToVectorPairIntInt mods;

	static StringVector modificationName;

	static MapIntToString fileNames;

	static MapCharToString constMods;
	static StringVector msfConstModsStr;

	static void clear ();
	static int scoreInfoCallback ( void* NotUsed, int argc, char** argv, char** azColName );
	static int peptideScoresCallback ( void* NotUsed, int argc, char** argv, char** azColName );
	static int peptidesCallback ( void* NotUsed, int argc, char** argv, char** azColName );
	static int spectrumHeadersCallback ( void* NotUsed, int argc, char** argv, char** azColName );
	static int peptidesAminoAcidModificationsCallback ( void* NotUsed, int argc, char** argv, char** azColName );
	static int peptidesTerminalModificationsCallback ( void* NotUsed, int argc, char** argv, char** azColName );
	static int aminoAcidModificationsCallback ( void* NotUsed, int argc, char** argv, char** azColName );
	static int fileInfosCallback ( void* NotUsed, int argc, char** argv, char** azColName );
	static int processingNodeParametersCallback ( void* NotUsed, int argc, char** argv, char** azColName );
	StringVector peakListFractionNames;
	StringVector peakListCentroidFileNames;
	StringVector peakListCentroidFilePaths;
	int readSpectrumRow ( sqlite3_stmt* pStmt, const std::string& path, ZipOpener& zo );
	void readDB ();
	void readSpectra ( const std::string& path );
	void readResults ( StringVectorVector& rows ) const;
	static void addConstMod ( const std::string& cm );
public:
	MSFRead ( const std::string& name );
	~MSFRead ();
	void readHeader ( StringVectorVector& header ) const;
	void read ( const std::string& peakListPath, StringVectorVector& rows );
	StringVector getPeakListFractionNames () const { return peakListFractionNames; }
	StringVector getPeakListCentroidFileNames () const { return peakListCentroidFileNames; }
	StringVector getPeakListCentroidFilePaths () const { return peakListCentroidFilePaths; }
	static void setConstMods ( const MapCharToString& cm );
	static MapCharToString getConstModMap ();
	static StringVector getMsfConstModsStr ();
};

#endif /* ! __lu_msf_h */
