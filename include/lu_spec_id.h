/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_spec_id.h                                                  *
*                                                                             *
*  Created    : December 13th 2003                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2002-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_spec_id_h
#define __lu_spec_id_h

#include <set>
#include <string>
#include <sstream>
class ParameterList;

#ifdef BATCHTAG
class SCMSTagLink;
#endif

class SpecID {
	int fraction;
	std::string spot;
	int run;
	int specNum;
	std::string msmsInfo;
public:
	SpecID ( const std::string& specID );
	SpecID ( const ParameterList* params );
	SpecID () {};
	int getFraction () const { return fraction; }
	void setFraction ( int f ) { fraction = f; }
	std::string getID () const { return spot; }
	double getSpotAsNumber () const { return atof ( spot.c_str () ); }
	int getRun () const { return run; }
	void setRun ( int r ) { run = r; }
	int getSpecNum () const { return specNum; }
	int getCharge () const { return specNum / 10000; }
	int getPureSpecNum100 () const { return specNum % 100; }
	int getPureSpecNum10000 () const { return specNum % 10000; }
	std::string getMSMSInfo () const { return msmsInfo; }
	void setMSMSInfo ( const std::string& s ) { msmsInfo = s; }
	std::string getSpecID () const 
	{
		std::ostringstream ost;
		ost << fraction << '-' << spot << '-' << run << '-' << specNum;
		return ost.str ();
	}
	std::string getFullSpecID () const 
	{
		std::ostringstream ost;
		ost << fraction << '-' << spot << '-' << run << '-' << specNum;
		if ( msmsInfo != "" ) ost << '-' << msmsInfo;
		return ost.str ();
	}
	std::string getSpotID () const 
	{
		std::ostringstream ost;
		ost << fraction << '-' << spot;
		return ost.str ();
	}
	std::string getPureSpecID100 () const
	{
		std::ostringstream ost;
		ost << fraction << '-' << spot << '-' << run << '-' << getPureSpecNum100 ();
		return ost.str ();
	}
	std::string getPureSpecID10000 () const
	{
		std::ostringstream ost;
		ost << fraction << '-' << spot << '-' << run << '-' << getPureSpecNum10000 ();
		return ost.str ();
	}
	bool isMSSpectrumToRead ( int f, const std::string& i ) const
	{
		if ( fraction == -1 ) return true;
		if ( f != fraction ) return false;
		if ( getChoppedSpot ( i ) != getChoppedSpot ( spot ) ) return false;
		return true;
	}
	bool isMSSpectrumToRead ( int f, const std::string& i, int r ) const
	{
		if ( fraction == -1 ) return true;
		if ( f != fraction ) return false;
		if ( getChoppedSpot ( i ) != getChoppedSpot ( spot ) ) return false;
		if ( r != run ) return false;
		return true;
	}
	bool isMSMSSpectrumToRead10000 ( int f, const std::string& i, int r, int s ) const
	{
		if ( fraction == -1 ) return true;
		if ( f != fraction ) return false;
		if ( getChoppedSpot ( i ) != getChoppedSpot ( spot ) ) return false;
		if ( r != run ) return false;
		if ( s != specNum % 10000 ) return false;
		return true;
	}
	bool isMSMSSpectrumToRead100 ( int f, const std::string& i, int r, int s ) const
	{
		if ( fraction == -1 ) return true;
		if ( f != fraction ) return false;
		if ( getChoppedSpot ( i ) != getChoppedSpot ( spot ) ) return false;
		if ( r != run ) return false;
		if ( s != specNum % 100 ) return false;
		return true;
	}
	static std::string getChoppedSpot ( const std::string& s )
	{
		size_t idx1 = s.find ( '.' );			// Find first dot (decimal point)
		if ( idx1 != std::string::npos ) {		// If there is a decimal point
			size_t idx2 = idx1;					// Nothing after decimal point - chop decimal point
			size_t len = s.length ();
			for ( size_t i = idx1+1 ; i < len ; i++ ) {
				if ( s [i] != '0' ) {
					idx2 = i + 1;
					if ( s [i] == '.' ) return s; // multiple dots
				}
			}
			if ( idx2 != len ) return s.substr ( 0, idx2 );
		}
		return s;
	}
	bool isMSMSSpectrumToRead100 ( const std::string& i, int r, int s ) const
	{
		if ( getChoppedSpot ( i ) != getChoppedSpot ( spot ) ) return false;
		if ( r != run ) return false;
		if ( s != specNum % 100 ) return false;
		return true;
	}
	bool isMSMSSpectrumToRead10000 ( const std::string& i, int r, int s ) const
	{
		if ( getChoppedSpot ( i ) != getChoppedSpot ( spot ) ) return false;
		if ( r != run ) return false;
		if ( s != specNum % 10000 ) return false;
		return true;
	}
	bool isFilterSpectrum () const { return fraction != -1; }
	void putCGI ( std::ostream& os ) const;
	void printHTMLHidden ( std::ostream& os ) const;
	std::string getCommandLineNVPair () const;
	static void printDelimitedHeader ( std::ostream& os, bool r, bool s );
	static void printDelimitedMSMSInfoHeader ( std::ostream& os );
	void printDelimitedCell ( std::ostream& os, bool r, bool s ) const;
	void printDelimitedMSMSInfoCell ( std::ostream& os ) const;
	static void printDelimitedEmptyCell ( std::ostream& os, bool r, bool s );
	static void printTableHeader ( std::ostream& os, bool r, bool s, const std::string& styleID, int rowspan = 0 );
	static void printTableMSMSInfoHeader ( std::ostream& os, const std::string& styleID, int rowspan = 0 );
	void printTableCell ( std::ostream& os, bool r, bool s, const std::string& styleID, int colspan = 0, int rowspan = 0 ) const;
	void printTableMSMSInfoCell ( std::ostream& os, const std::string& styleID, int colspan = 0, int rowspan = 0 ) const;
#ifdef BATCHTAG
	void printTableCell2 ( std::ostream& os, bool r, bool s, const std::string& styleID, const SCMSTagLink& smtl, const std::string& searchKey, int colspan = 0, int rowspan = 0 ) const;
#endif
	static void printTableEmptyCell ( std::ostream& os, bool r, bool s, const std::string& styleID );
	static std::string convertTime ( const std::string& oldTime );
};
inline bool operator< ( const SpecID& lhs, const SpecID& rhs )
{
	if ( lhs.getFraction () == rhs.getFraction () ) {
		if ( lhs.getID () == rhs.getID () ) {
			if ( lhs.getRun () == rhs.getRun () ) {
				if ( lhs.getPureSpecNum10000 () == rhs.getPureSpecNum10000 () )
					return lhs.getCharge () < rhs.getCharge ();
				else
					return lhs.getPureSpecNum10000 () < rhs.getPureSpecNum10000 ();
			}
			else return lhs.getRun () < rhs.getRun ();
		}
		else {
			return atof ( lhs.getID ().c_str () ) < atof ( rhs.getID ().c_str () );
		}
	}
	else return lhs.getFraction () < rhs.getFraction ();
}

inline bool operator== ( const SpecID& lhs, const SpecID& rhs )
{
	if ( lhs.getFraction () != rhs.getFraction () ) return false;
	if ( lhs.getID () != rhs.getID () ) return false;
	if ( lhs.getRun () != rhs.getRun () ) return false;
	if ( lhs.getSpecNum () != rhs.getSpecNum () ) return false;
	return true;
}
typedef std::set <SpecID> SetSpecID;
typedef SetSpecID::const_iterator SetSpecIDConstIterator;

#endif /* ! __lu_spec_id_h */
