 /*****************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_msp.h                                                      *
*                                                                             *
*  Created    : March 6th 2014                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2014-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_msp_h
#define __lu_msp_h

class LibraryFileRead {
protected:
	std::string dbPeptide;
	std::string mods;
	std::string sPrecursorMZ;
	std::string sPrecursorZ;
	std::string peptideType;
	std::string prevAA;
	std::string nextAA;
	std::string spectrum;
	std::string sNumPeaks;
	std::string inst;
	std::string protein;
	int specNum;

	const std::string& fname;
	StringVector peakListFractionNames;
	StringVector peakListCentroidFileNames;
	StringVector peakListCentroidFilePaths;
	virtual void processFile ( std::istream& ist, std::ostream& os, StringVectorVector& rows ) = 0;
	void writeSpectrum ( std::istream& istr, std::ostream& os, int num, const std::string& sPrecursorMZ, const std::string& sPrecursorZ, int numPeaks ) const;
	std::string getMods ( const std::string& value ) const;
	void addRow ( StringVectorVector& rows ) const;
public:
	LibraryFileRead ( const std::string& name );
	void read ( const std::string& peakListPath, StringVectorVector& header, StringVectorVector& rows );
	StringVector getPeakListFractionNames () const { return peakListFractionNames; }
	StringVector getPeakListCentroidFileNames () const { return peakListCentroidFileNames; }
	StringVector getPeakListCentroidFilePaths () const { return peakListCentroidFilePaths; }
};

class MSPRead : public LibraryFileRead {
	void processFile ( std::istream& ist, std::ostream& os, StringVectorVector& rows );
public:
	MSPRead ( const std::string& name );
};

class SPTXTRead : public LibraryFileRead {
	void processFile ( std::istream& ist, std::ostream& os, StringVectorVector& rows );
public:
	SPTXTRead ( const std::string& name );
};

#endif /* ! __lu_msp_h */
