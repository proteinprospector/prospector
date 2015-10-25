/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lu_t2d.h                                                      *
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
*  Copyright (2004-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_t2d_h
#define __lu_t2d_h

#include <string>
#include <nr.h>
#include <lgen_define.h>

struct DBProvider {
    enum {
        MaxDBProviderLength = 128
    };
    char s[MaxDBProviderLength];
};

typedef union {
	// Current 64 bit unsigned integer
	GENUINT64 dlKey;
	// Possible GUID or other 128 bit format
	unsigned char gKey[16];
} DBKEY;

enum PBSpectrumFileSignatures {
    PBSpectrumFileHeaderSignature  = 0xA55A1234,
    PBSpectrumFileTrailerSignature = 0x5AA54321,
    PBSpectrumFileSpectrumSignature = 0x5A5A1243
};

class T2DHeader {
	unsigned short checksum;				// 2 bytes
	unsigned short fileHeaderSize;			// 2 bytes
	PBSpectrumFileSignatures signature;		// 4 bytes
	unsigned short softwareMajorVersion;	// 2 bytes
	unsigned short softwareMinorVersion;	// 2 bytes
	unsigned short fileMajorVersion;		// 2 bytes
	unsigned short fileMinorVersion;		// 2 bytes
	unsigned short spectrumHeaderSize;		// 2 bytes
	unsigned short fileTrailerSize;			// 2 bytes
	DBProvider instrumentDBProvider;		// 128 Bytes
	DBKEY fileID;							// 16 bytes
public:
	T2DHeader ();
	T2DHeader ( GenIFStream& ifs );
	static unsigned short size ()
	{
		return ( 8 * sizeof (unsigned short) ) + sizeof (PBSpectrumFileSignatures) + sizeof (DBProvider) + sizeof (DBKEY);
	}
	void write ( std::ostream& os ) const;
	void printHTML ( std::ostream& os ) const;
};

class T2DTrailer {							// Total 20 bytes
	unsigned short checksum;				// 2 bytes
	PBSpectrumFileSignatures signature;		// 4 bytes
	unsigned short indexChecksum;			// 2 bytes
	unsigned long numSpectra;				// 4 bytes
	GENINT64 spectrumIndexOffset;			// 8 bytes
public:
	T2DTrailer ( unsigned long indexOffset );
	T2DTrailer ( GenIFStream& ifs );
	static unsigned short size ()
	{
		return ( 2 * sizeof (unsigned short) ) + sizeof (PBSpectrumFileSignatures) + sizeof (unsigned long) + sizeof (GENINT64);
	}
	void write ( std::ostream& os ) const;
	void printHTML ( std::ostream& os ) const;
};

class T2DFile {
	XYData xyData;
public:
	T2DFile ( const std::string& filename, double startMass = 0.0, double endMass = 0.0 );
	void printASCII ( std::ostream& os ) const;
	const XYData& getData () const { return xyData; }
};

int getT2DMSRunNumberIndex ( int msmsRunNumber, const IntVector& msRunNumbers );

#endif /* ! __lu_t2d_h */
