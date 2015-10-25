/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_pq_vector.h                                                *
*                                                                             *
*  Created    : July 18th 2000                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2000-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_pq_vector_h
#define __lu_pq_vector_h

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <lgen_define.h>
#include <lg_string.h>

static void getPostQueryVector ( const char* data, StringVector& valVector, char delim )
{
	std::string dataString ( data );
	std::istringstream ist ( dataString );
	std::string line;
	while ( std::getline ( ist, line, delim ) ) {
		if ( line.length () != 0 ) {
			if ( line [0] != '#' ) {
				if ( line [line.length () - 1] == '\r' ) line = line.substr ( 0, line.length () - 1 );	// Delete a trailing dot
				std::string line2 = gen_strtrim ( line );	// Trim white space off start and end of string
				if ( line2.length () ) {
					valVector.push_back ( line2 );
				}
			}
		}
	}
};
template<class T> void getPostQueryVector ( const char* data, typename std::vector<T>& valVector, char delim )
{
	std::string dataString ( data );
	std::istringstream ist ( dataString );
	std::string line;
	while ( std::getline ( ist, line, delim ) ) {
		if ( !line.empty () ) {
			if ( line [line.length () - 1] == '\r' ) line = line.substr ( 0, line.length () - 1 );	// Delete a trailing dot
			line = gen_strtrim ( line );	// Trim white space off start and end of string
			if ( !line.empty () ) {
				if ( line [0] != '#' ) {

					std::istringstream istLine ( line );

					T value;
					if ( istLine >> value ) {
						valVector.push_back ( value );
						std::string temp;
						if ( istLine >> temp ) {
							valVector.clear ();
							throw std::runtime_error ( "Too many fields entered." );
						}
					}
				}
			}
		}
	}
};
template<class T1, class T2> bool getPostQuery2Vectors ( const char* data, typename std::vector<T1>& valVector1, typename std::vector<T2>& valVector2 )
{
	std::string dataString ( data );
	std::istringstream ist ( dataString );
	std::string line;
	bool flag = true;
	while ( std::getline ( ist, line ) ) {
		if ( !line.empty () ) {
			if ( line [line.length () - 1] == '\r' ) line = line.substr ( 0, line.length () - 1 );	// Delete a trailing dot
			line = gen_strtrim ( line );	// Trim white space off start and end of string
			if ( !line.empty () ) {
				if ( line [0] != '#' ) {

					std::istringstream istLine ( line );
					T1 value1;
					if ( istLine >> value1 ) {
						valVector1.push_back ( value1 );
						T2 value2;
						if ( istLine >> value2 ) {
							valVector2.push_back ( value2 );
							std::string temp;
							if ( istLine >> temp ) {
								valVector1.clear ();
								valVector2.clear ();
								throw std::runtime_error ( "Too many fields entered." );
							}
						}
						else {
							flag = false;
							break;
						}
					}
					else {
						flag = false;
						break;
					}
				}
			}
		}
	}
	if ( flag == false ) {
		valVector1.clear ();
		valVector2.clear ();
		return false;
	}
	else {
		return true;
	}
};
template<class T1, class T2, class T3> bool getPostQuery3Vectors ( const char* data, typename std::vector<T1>& valVector1, typename std::vector<T2>& valVector2, typename std::vector<T3>& valVector3 )
{
	std::string dataString ( data );
	std::istringstream ist ( dataString );
	std::string line;
	bool flag = true;
	while ( std::getline ( ist, line ) ) {
		if ( !line.empty () ) {

			if ( line [line.length () - 1] == '\r' ) line = line.substr ( 0, line.length () - 1 );	// Delete a trailing dot
			line = gen_strtrim ( line );	// Trim white space off start and end of string
			if ( !line.empty () ) {
				if ( line [0] != '#' ) {

					std::istringstream istLine ( line );
					T1 value1;
					if ( istLine >> value1 ) {
						valVector1.push_back ( value1 );
						T2 value2;
						if ( istLine >> value2 ) {
							valVector2.push_back ( value2 );

							T3 value3;
							if ( istLine >> value3 ) {
								valVector3.push_back ( value3 );
								std::string temp;
								if ( istLine >> temp ) {
									valVector1.clear ();
									valVector2.clear ();
									valVector3.clear ();
									throw std::runtime_error ( "Too many fields entered." );
								}
							}
							else {
								flag = false;
								break;
							}
						}
						else {
							flag = false;
							break;
						}
					}
					else {
						flag = false;
						break;
					}
				}
			}
		}
	}
	if ( flag == false ) {
		valVector1.clear ();
		valVector2.clear ();
		valVector3.clear ();
		return false;
	}
	else {
		return true;
	}
};

#endif /* ! __lu_pq_vector_h */
