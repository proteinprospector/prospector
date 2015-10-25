/******************************************************************************
*                                                                             *
*  Library    : libxcalibur                                                   *
*                                                                             *
*  Filename   : lz_app.h                                                      *
*                                                                             *
*  Created    : February 19th 2015                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2015-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef __lz_app_h
#define __lz_app_h

typedef struct _MEMFILE MEMFILE;
typedef struct zlib_filefunc_def_s zlib_filefunc_def;

class ZipOpener {
	char* buf;
	int len;
	int bufLen;
	int ret;
	MEMFILE* handle;
	zlib_filefunc_def* api;
public:
	ZipOpener ();
	~ZipOpener ();
	void getFirstFile ( void* buffer, unsigned int bufferLen );
	char* getBuf () { return buf; }
	int getLen () const { return len; }
};

#endif /* ! __lz_app_h */
