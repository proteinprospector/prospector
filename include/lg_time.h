/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lg_time.h                                                     *
*                                                                             *
*  Created    : May 22nd 1997                                                 *
*                                                                             *
*  Purpose    : Machine independent interface to time.h.                      *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1997-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lg_time_h
#define __lg_time_h

#include <string>
#include <ctime>
#include <lgen_define.h>

class GenTime {
	std::string dayOfWeek;
	std::string month;
	std::string intMonth;
	std::string dayOfMonth;
	std::string hour;
	std::string minute;
	std::string second;
	std::string year;
public:
	GenTime ( const std::string& timeAndDate );
	std::string getUSDate ( char sep = ' ' ) const;
	std::string getDateAndTime ( char sep = ' ' ) const;
	std::string getDBDateAndTime () const;
	std::string getXSDDateAndTime () const;
	std::string getYear () const { return year; }
	std::string getMonth () const { return month; }
	std::string getIntMonth () const { return intMonth; }
	std::string getDayOfMonth () const { return dayOfMonth; }
	std::string getHour () const { return hour; }
	std::string getMinute () const { return minute; }
	std::string getSecond () const { return second; }
	std::string getMonthAndYear ( char sep = ' ' ) const;
	std::string getYearAndIntMonth ( char sep = ' ' ) const;
	std::string getYearIntMonthAndDay ( char sep = ' ' ) const;
	std::string getHoursMinSec ( char sep = ' ' ) const;
};

inline bool operator< ( const GenTime& lhs, const GenTime& rhs )
{
	int lYear = atoi ( lhs.getYear ().c_str () );
	int rYear = atoi ( rhs.getYear ().c_str () );
	if ( lYear == rYear ) {
		int lMonth = atoi ( lhs.getIntMonth ().c_str () );
		int rMonth = atoi ( rhs.getIntMonth ().c_str () );
		if ( lMonth == rMonth ) {
			int lDayOfMonth = atoi ( lhs.getDayOfMonth ().c_str () );
			int rDayOfMonth = atoi ( rhs.getDayOfMonth ().c_str () );
			if ( lDayOfMonth == rDayOfMonth ) {
				int lHour = atoi ( lhs.getHour ().c_str () );
				int rHour = atoi ( rhs.getHour ().c_str () );
				if ( lHour == rHour ) {
					int lMinute = atoi ( lhs.getMinute ().c_str () );
					int rMinute = atoi ( rhs.getMinute ().c_str () );
					if ( lMinute == rMinute ) {
						int lSecond = atoi ( lhs.getSecond ().c_str () );
						int rSecond = atoi ( rhs.getSecond ().c_str () );
						return lSecond < rSecond;

					}
					else
						return lMinute < rMinute;
				}
				else
					return lHour < rHour;
			}
			else
				return lDayOfMonth < rDayOfMonth;
		}
		else
			return lMonth < rMonth;
	}
	else
		return lYear < rYear;
}
inline bool operator> ( const GenTime& lhs, const GenTime& rhs )
{
	int lYear = atoi ( lhs.getYear ().c_str () );
	int rYear = atoi ( rhs.getYear ().c_str () );
	if ( lYear == rYear ) {
		int lMonth = atoi ( lhs.getIntMonth ().c_str () );
		int rMonth = atoi ( rhs.getIntMonth ().c_str () );
		if ( lMonth == rMonth ) {
			int lDayOfMonth = atoi ( lhs.getDayOfMonth ().c_str () );
			int rDayOfMonth = atoi ( rhs.getDayOfMonth ().c_str () );
			if ( lDayOfMonth == rDayOfMonth ) {
				int lHour = atoi ( lhs.getHour ().c_str () );
				int rHour = atoi ( rhs.getHour ().c_str () );
				if ( lHour == rHour ) {
					int lMinute = atoi ( lhs.getMinute ().c_str () );
					int rMinute = atoi ( rhs.getMinute ().c_str () );
					if ( lMinute == rMinute ) {
						int lSecond = atoi ( lhs.getSecond ().c_str () );
						int rSecond = atoi ( rhs.getSecond ().c_str () );
						return lSecond > rSecond;

					}
					else
						return lMinute > rMinute;
				}
				else
					return lHour > rHour;
			}
			else
				return lDayOfMonth > rDayOfMonth;
		}
		else
			return lMonth > rMonth;
	}
	else
		return lYear > rYear;
}

std::string genTimeAndDateString ( time_t t );
std::string genCurrentTimeAndDateString ();
std::string genCurrentDateString ( char sep = ' ' );
std::string genCurrentDateAndTimeString ( char sep = ' ' );
std::string genCurrentDBDateAndTimeString ();
std::string genCurrentXSDDateAndTimeString ( time_t t );
std::string genCurrentXSDDateAndTimeString ();
std::string genCurrentYear ();
std::string genCurrentMonth ();
std::string genCurrentIntMonth ();
std::string genCurrentMonthAndYear ( char sep = ' ' );
std::string genCurrentYearAndIntMonth ( char sep = ' ' );
std::string genCurrentYearIntMonthAndDay ( char sep = ' ' );
std::string genCurrentHoursMinSec ( char sep = ' ' );

class GenElapsedTime {
	time_t startTime;
public:
	GenElapsedTime ();
	time_t getElapsedTime () const;
	std::string getElapsedTimeString () const;
	void reset ();
	static std::string getTimeString ( time_t seconds );
};
bool genMoreThan2DaysOld ( const std::string& date );
bool genMoreThanNDaysOld ( const std::string& date, int n );

StringVector getYearList ( const std::string& startYear, bool reverseFlag = false );

std::string secToMins ( const std::string& sec );

#endif /* ! __lg_time_h */
