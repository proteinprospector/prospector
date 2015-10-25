/******************************************************************************
*                                                                             *
*  Program    : libzip                                                        *
*                                                                             *
*  Filename   : lz_app.cpp                                                    *
*                                                                             *
*  Created    : February 18th 2015                                            *
*                                                                             *
*  Purpose    : This library was written to allow zip files to be read in     *
*               memory.                                                       *
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
#ifndef VIS_C
#include <cstring>
#include <cstdlib>
#endif
#include "fseek.h"
#include "unzipx.h"
#include <lz_app.h>

typedef struct _MEMFILE
{
  void* buffer;    // Base of the region of memory we're using
  ZPOS_T length;   // Size of the region of memory we're using
  ZPOS_T position; // Current offset in the area
} MEMFILE;

static uLong ZCALLBACK mem_read (voidpf opaque, voidpf stream, void* buf, uLong size)
{
   MEMFILE* handle = (MEMFILE*) stream;

   if ( (handle->position + size) > handle->length)
   {
      size = handle->length - handle->position;
   }
   memcpy(buf, ((char*)handle->buffer) + handle->position, size);
   handle->position+=size;
   return size;
}
static uLong ZCALLBACK mem_write (voidpf opaque, voidpf stream, const void* buf, uLong size)
{
   MEMFILE* handle = (MEMFILE*) stream;

   if ((handle->position + size) > handle->length)
   {
      handle->length = handle->position + size;
      handle->buffer = realloc(handle->buffer, handle->length);
   }

   memcpy(((char*)handle->buffer) + handle->position, buf, size);
   handle->position+=size;

   return size;
}
static ZPOS_T ZCALLBACK mem_tell (voidpf opaque, voidpf stream)
{
   MEMFILE *handle = (MEMFILE *)stream;
   return handle->position;
}
static long ZCALLBACK mem_seek (voidpf opaque, voidpf stream, ZOFF_T offset, int origin)
{
   MEMFILE* handle = (MEMFILE*)stream;
   return fseek_calc(offset, origin, &handle->position, handle->length);
}
int ZCALLBACK mem_close (voidpf opaque, voidpf stream)
{
    MEMFILE *handle = (MEMFILE *)stream;
    // Note that once we've written to the buffer we don't tell anyone
    // about it here. Probably the opaque handle could be used to inform
    // some other component of how much data was written.

    // This, and other aspects of writing through this interface, has
    // not been tested.

    //free (handle->buffer);		// This buffer is passed externally so should not be freed - PB
    //free (handle);				// The handle is freed elsewhere
    return 0;
}

int ZCALLBACK mem_error (voidpf opaque, voidpf stream)
{
    MEMFILE *handle = (MEMFILE *)stream;
    // We never return errors
    return 0;
}

ZipOpener::ZipOpener () :
	buf ( 0 ),
	len ( 0 ),
	bufLen ( 0 )
{
	api = new zlib_filefunc_def;
    api->zopen_file  = NULL;
    api->zread_file  = mem_read;
    api->zwrite_file = mem_write;
    api->ztell_file  = mem_tell;
    api->zseek_file  = mem_seek;
    api->zclose_file = mem_close;
    api->zerror_file = mem_error;
    api->opaque      = handle;
    handle = new MEMFILE;
}
#include <iostream>
ZipOpener::~ZipOpener ()
{
	ZCLOSE(*api, handle);
	delete api;
	delete handle;
	delete [] buf;
}
void ZipOpener::getFirstFile ( void* buffer, unsigned int bufferLen )
{
    handle->position = 0;
    handle->buffer   = buffer;
    handle->length   = bufferLen;
	unzFile unz = unzAttach ( handle, api );
	if ( unz ) {
		ret = unzGoToFirstFile ( unz );											// This function assumes there's only one file
		unz_file_info info;
		ret = unzGetCurrentFileInfo ( unz, &info, NULL, 0, NULL, 0, NULL, 0 );	// Get the info for the file length
		ret = unzOpenCurrentFile ( unz );
		len = info.uncompressed_size;
		if ( len+1 > bufLen ) {													// This function is supposed to be called repeatedly to get an uncompressed buffer
			delete [] buf;														// Make sure the buffer is large enough
			buf = new char [len+1];
			bufLen = len+1;
		}
		ret = unzReadCurrentFile  ( unz, buf, len );
		buf [ret] = 0;
		unzCloseCurrentFile ( unz );
		unzDetach(&unz);
	}
}
