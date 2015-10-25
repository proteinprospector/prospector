/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lgen_mmap.h                                                   *
*                                                                             *
*  Created    : June 20th 1996                                                *
*                                                                             *
*  Purpose    : Machine independent memory mapping functions.                 *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lgen_mmap_h
#define __lgen_mmap_h

#ifdef VIS_C
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <lg_stdio.h>
#endif

#include <string>

#include <lgen_define.h>
#include <lu_getfil.h>
#include <lgen_error.h>
#include <lgen_file.h>

template <class T> class MMapFile {
	T* startPointer;
	T* endPointer;
	int granularity;
	GENINT64 mapStart;
	GENINT64 mapEnd;
	GENINT64 byteMapSize;
	GENINT64 fileSize;
	GENINT64 mapLimit;
	//void printSystemInfo ();
	void createView ( GENINT64 startOffset, GENINT64 endOffset = 0 );
	void newView ( GENINT64 startOffset, GENINT64 endOffset = 0 );
#ifdef VIS_C
	HANDLE hFile;
	HANDLE hmmf;
#else
	int fd;
#endif
public:
	MMapFile ( const std::string& name, GENINT64 startOffset = 0, GENINT64 mapLimit = 0 );
	~MMapFile ();
	T& subscript ( GENINT64 index )
	{
		if ( index < mapStart || index > mapEnd ) {
			newView ( index );
		}
		return startPointer [index - mapStart];
	}
	T* getRange ( GENINT64 startIndex, GENINT64 endIndex )
	{
		if ( startIndex < mapStart || endIndex > mapEnd ) {
			newView ( startIndex, endIndex );
		}
		return startPointer + ( startIndex - mapStart );
	}
	T* nextRange ()
	{
		newView ( mapEnd + 1 );
		return startPointer;
	}
	T* getStartPointer () const { return startPointer; }
	T* getEndPointer () const { return endPointer; }
	GENINT64 getMapOffset ( const T* pointer ) const { return mapStart + ( pointer - startPointer ); }
	GENINT64 getFileSize () const { return fileSize; }
};
template <class T>
MMapFile<T>::MMapFile ( const std::string& name, GENINT64 startOffset, GENINT64 mapLimit )
{
#ifdef VIS_C
	hFile = CreateFile ( name.c_str (), GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_READONLY, 0 );
	if ( hFile == INVALID_HANDLE_VALUE ) ErrorHandler::genError ()->error ( "CreateFile failed in MMapFile" );
	hmmf = CreateFileMapping ( hFile, NULL, PAGE_READONLY, 0, 0, genFilenameFromPath ( name ).c_str () );
	if ( hmmf == NULL ) ErrorHandler::genError ()->error ( "CreateFileMapping failed in MMapFile: " + genGetLastWindowsError () );
#else
	fd = gen_open ( name, O_RDONLY, "Function name: gen_open_read_only_mmap_file." );
#endif
	fileSize = genFileSize ( name );
	GenNameValueStream nvs ( MsparamsDir::instance ().getParamPath ( "computer.txt" ) );
	nvs.getValue ( "block_size", granularity );
	int numBlocks;
	nvs.getValue ( "num_blocks", numBlocks );

	if ( mapLimit == 0 ) {
		mapLimit = granularity * numBlocks;
	}
	this->mapLimit = mapLimit;
	createView ( startOffset );
}
template <class T>
MMapFile<T>::~MMapFile ()
{
#ifdef VIS_C
	UnmapViewOfFile ( startPointer );
	CloseHandle ( hmmf );
	CloseHandle ( hFile );
#else
	if ( munmap ( startPointer, byteMapSize ) == -1 ) {
		ErrorHandler::genError ()->error ( "Memory page unmapping (munmap) failure.\n" );
	}
	gen_close ( fd, "Function name: close_mmap_file." );
#endif
}
template <class T>
void MMapFile<T>::createView ( GENINT64 startOffset, GENINT64 endOffset )
{
	GENINT64 byteStartOffset = startOffset * sizeof (T);
	GENINT64 byteMapStart = ( byteStartOffset / granularity ) * granularity;
	if ( endOffset != 0 ) {
		GENINT64 byteEndOffset = endOffset * sizeof (T);
		GENINT64 byteMapEnd = ( 1 + ( byteEndOffset / granularity ) ) * granularity;
		byteMapSize = genMax ( mapLimit, byteMapEnd - byteMapStart );
		byteMapSize = genMin ( byteMapSize, fileSize - byteMapStart );
	}
	else {
		byteMapSize = genMin ( mapLimit, fileSize - byteMapStart );
	}
#ifdef VIS_C
	unsigned int high32 = static_cast <unsigned int> ( byteMapStart >> 32 );
	unsigned int low32 = static_cast <unsigned int> ( byteMapStart & 0xFFFFFFFF );
	startPointer = (T*)MapViewOfFile ( hmmf, FILE_MAP_READ, high32, low32, byteMapSize );
	if ( startPointer == NULL ) ErrorHandler::genError ()->error ( "MapViewOfFile failed in MMapFile" );
#else
	ErrorHandler::genError ()->resetErrorNumber ();
	startPointer = (T*) mmap ( 0, byteMapSize, PROT_READ, MAP_SHARED, fd, byteMapStart );	// addr, len, prot, flags, fd, off
	if ( ErrorHandler::genError ()->getErrorNumber () != 0 ) {
		ErrorHandler::genError ()->error ( "Memory page mapping (mmap) failure.\n" );
	}
#endif
	mapStart = byteMapStart / sizeof (T);
	GENINT64 mapSize = byteMapSize /  sizeof (T);
	mapEnd = mapStart + mapSize - 1;
	endPointer = startPointer + ( mapSize - 1 );
}
template <class T>
void MMapFile<T>::newView ( GENINT64 startIndex, GENINT64 endIndex )
{
#ifdef VIS_C
	UnmapViewOfFile ( startPointer );
#else
	if ( munmap ( startPointer, byteMapSize ) == -1 ) {
		ErrorHandler::genError ()->error ( "Memory page unmapping (munmap) failure.\n" );
	}
#endif
	createView ( startIndex, endIndex );
}
#ifdef VIS_C
/*
template <class T>
void MMapFile<T>::printSystemInfo ()
{
	SYSTEM_INFO systemInfo;

	GetSystemInfo ( &systemInfo );
	cout << "Page size and the granularity of page protection and commitment = " << dec << systemInfo.dwPageSize;
	cout << " (0x" << hex << systemInfo.dwPageSize << ")<br />";
	cout << "Pointer to the lowest memory address accessible to applications and DLLs = " << dec << (unsigned int) systemInfo.lpMinimumApplicationAddress;	// cast needed to prevent printing in hex
	cout << " (0x" << hex << systemInfo.lpMinimumApplicationAddress << ")<br />";
	cout << "Pointer to the highest memory address accessible to applications and dynamic-link libraries (DLLs) = " << dec << (unsigned int) systemInfo.lpMaximumApplicationAddress;	// cast needed to prevent printing in hex
	cout << " (0x" << hex << systemInfo.lpMaximumApplicationAddress << ")<br />";
	cout << "Mask representing the set of processors configured into the system (bit 0 is processor 0; bit 31 is processor 31) = " << dec << systemInfo.dwActiveProcessorMask;
	cout << " (0x" << hex << systemInfo.dwActiveProcessorMask << ")<br />";
	cout << "Number of processors in the system = " << dec << systemInfo.dwNumberOfProcessors;
	cout << " (0x" << hex << systemInfo.dwNumberOfProcessors << ")<br />";
	//cout << dec << systemInfo.dwProcessorType "<br />";	No longer relevant.
	cout << "Granularity with which virtual memory is allocated = " << dec << systemInfo.dwAllocationGranularity;
	cout << " (0x" << hex << systemInfo.dwAllocationGranularity << ")<br />";
	cout << "Processor architecture = " << dec << systemInfo.wProcessorArchitecture;
	cout << " (0x" << hex << systemInfo.wProcessorArchitecture << ")<br />";
	cout << "Processor level = " << dec << systemInfo.wProcessorLevel;
	cout << " (0x" << hex << systemInfo.wProcessorLevel << ")<br />";
	cout << "Processor revision = " << dec << systemInfo.wProcessorRevision;
	cout << " (0x" << hex << systemInfo.wProcessorRevision << ")<br />";
	cout << dec;
}
*/
#endif

#endif /* ! __lgen_mmap_h */
