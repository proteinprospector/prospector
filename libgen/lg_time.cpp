/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lg_time.cpp                                                   *
*                                                                             *
*  Created    : May 22nd 1997                                                 *
*                                                                             *
*  Purpose    : Time functions.                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1997-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <lg_io.h>
#include <lg_string.h>
#include <lg_time.h>
using std::string;
using std::istringstream;
using std::ostringstream;
using std::stringstream;
using std::reverse;
using std::setw;
using std::setfill;

class MonthMap {
	static const char* monthStrings [];
	static const char* numMonthStrings [];
	MapStringToString monthMap;
	MonthMap ()
	{
		for ( int i = 0 ; monthStrings [i] != 0 ; i++ ) {
			monthMap [monthStrings [i]] = numMonthStrings [i];
		}
	}
public:
	static MonthMap& instance ()
	{
		static MonthMap d;
		return d;
	}
	string getNumericMonth ( const string& m ) { return monthMap [m]; }
};
const char* MonthMap::monthStrings [] = {
	"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", 0 };
const char* MonthMap::numMonthStrings [] = {
	"01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", 0 };

GenTime::GenTime ( const string& timeAndDate )
{
	istringstream istr ( timeAndDate );
	string t;
	istr >> dayOfWeek;
	istr >> month;
	intMonth = MonthMap::instance ().getNumericMonth ( month );
	istr >> dayOfMonth;
	if ( dayOfMonth.length () == 1 ) dayOfMonth = "0" + dayOfMonth;
	istr >> t;
	hour = t.substr ( 0, 2 );
	minute = t.substr ( 3, 2 );
	second = t.substr ( 6, 2 );
	istr >> year;
}
string GenTime::getUSDate ( char sep ) const
{
	if ( sep == 0 ) return month + dayOfMonth + year;
	else return month + sep + dayOfMonth + sep + year;
}
string GenTime::getDateAndTime ( char sep ) const
{
	if ( sep == 0 ) return year + intMonth + dayOfMonth + hour + minute + second;
	else return year + sep + intMonth + sep + dayOfMonth + sep + hour + sep + minute + sep + second;
}
string GenTime::getDBDateAndTime () const
{
	return year + "-" + intMonth + "-" + dayOfMonth + " " + hour + ":" + minute + ":" + second;
}
string GenTime::getXSDDateAndTime () const
{
	return year + "-" + intMonth + "-" + dayOfMonth + "T" + hour + ":" + minute + ":" + second;
}
string GenTime::getMonthAndYear ( char sep ) const
{
	if ( sep == 0 ) return month + year;
	else return month + sep + year;
}
string GenTime::getYearAndIntMonth ( char sep ) const
{
	if ( sep == 0 ) return year + intMonth;
	else return year + sep + intMonth;
}
string GenTime::getYearIntMonthAndDay ( char sep ) const
{
	if ( sep == 0 ) return year + intMonth + dayOfMonth;
	else return year + sep + intMonth + sep + dayOfMonth;
}
string GenTime::getHoursMinSec ( char sep ) const
{
	if ( sep == 0 ) return hour + minute + second;
	else return hour + sep + minute + sep + second;
}
string genTimeAndDateString ( time_t t )
{
	const char* ct = ctime ( &t );
	return string ( ct, strlen ( ct ) - 1 ); // delete \n character
}
// Fri Jul 19 12:21:10 2013
string genCurrentTimeAndDateString ()
{
	return genTimeAndDateString ( time ( 0 ) );
}
string genCurrentDateString ( char sep )
{
	GenTime t ( genTimeAndDateString ( time ( 0 ) ) );
	return t.getUSDate ( sep );
}
string genCurrentDateAndTimeString ( char sep )
{
	GenTime t ( genTimeAndDateString ( time ( 0 ) ) );
	return t.getDateAndTime ( sep );
}
string genCurrentDBDateAndTimeString ()
{
	GenTime t ( genTimeAndDateString ( time ( 0 ) ) );
	return t.getDBDateAndTime ();
}
string genCurrentXSDDateAndTimeString ( time_t t1 )
{
	GenTime t ( genTimeAndDateString ( t1 ) );
	return t.getXSDDateAndTime ();
}
string genCurrentXSDDateAndTimeString ()
{
	GenTime t ( genTimeAndDateString ( time ( 0 ) ) );
	return t.getXSDDateAndTime ();
}
string genCurrentYear ()
{
	GenTime t ( genTimeAndDateString ( time ( 0 ) ) );
	return t.getYear ();
}
string genCurrentMonth ()
{
	GenTime t ( genTimeAndDateString ( time ( 0 ) ) );
	return t.getMonth ();
}
string genCurrentIntMonth ()
{
	GenTime t ( genTimeAndDateString ( time ( 0 ) ) );
	return t.getIntMonth ();
}
string genCurrentMonthAndYear ( char sep )
{
	GenTime t ( genTimeAndDateString ( time ( 0 ) ) );
	return t.getMonthAndYear ( sep );
}
string genCurrentYearAndIntMonth ( char sep )
{
	GenTime t ( genTimeAndDateString ( time ( 0 ) ) );
	return t.getYearAndIntMonth ( sep );
}
string genCurrentYearIntMonthAndDay ( char sep )
{
	GenTime t ( genTimeAndDateString ( time ( 0 ) ) );
	return t.getYearIntMonthAndDay ( sep );
}
string genCurrentHoursMinSec ( char sep )
{
	GenTime t ( genTimeAndDateString ( time ( 0 ) ) );
	return t.getHoursMinSec ( sep );
}
GenElapsedTime::GenElapsedTime () :
	startTime ( time ( 0 ) )
{
}
time_t GenElapsedTime::getElapsedTime () const
{
	return time ( 0 ) - startTime;
}
string GenElapsedTime::getElapsedTimeString () const
{
	time_t secondsElapsed = time ( 0 ) - startTime;
	return getTimeString ( secondsElapsed );
}
void GenElapsedTime::reset ()
{
	startTime = time ( 0 );
}
string GenElapsedTime::getTimeString ( time_t seconds )
{
	int days = seconds / 86400;
	int rem = seconds % 86400;
	int hours = rem / 3600;
	rem = rem % 3600;
	int minutes = rem / 60;
	int sec = rem % 60;
	ostringstream timeString;
	if ( days ) timeString << days << " dy ";
	if ( days || hours ) timeString << hours << " hr ";
	if ( days || hours || minutes ) timeString << minutes << " min ";
	timeString << sec << " sec";
	return timeString.str ();
}
static time_t genMakeTimeFromDateString ( const string& date )
{
	int st = 0;
	int pos = date.find ( "_" );
	if ( pos == string::npos ) return 0;
	string month = date.substr ( st, pos );
	st = pos+1;
	pos = date.find ( "_", st );
	if ( pos == string::npos ) return 0;
	int day = atoi ( date.substr ( st, pos ).c_str () );
	st = pos+1;
	if ( st >= date.length () - 1 ) return 0;
	int year = atoi ( date.substr ( st ).c_str () );
	time_t t = time ( 0 ) - 1;
	struct tm* timeinfo = localtime ( &t );
	if ( year < 2000 ) return 0;
	timeinfo->tm_year = year - 1900;
	timeinfo->tm_mon = atoi ( MonthMap::instance ().getNumericMonth ( month ).c_str () ) - 1;
	timeinfo->tm_mday = day;
	return mktime ( timeinfo );
}
static time_t genMakeTimeFromDateString2 ( const string& date )
{
	int st = 0;
	int pos = date.find ( "_" );
	if ( pos == string::npos ) return 0;
	int year = atoi ( date.substr ( st, pos ).c_str () );
	st = pos+1;
	pos = date.find ( "_", st );
	if ( pos == string::npos ) return 0;
	int month = atoi ( date.substr ( st, pos ).c_str () );
	st = pos+1;
	if ( st >= date.length () - 1 ) return 0;
	int day = atoi ( date.substr ( st ).c_str () );
	time_t t = time ( 0 ) - 1;
	struct tm* timeinfo = localtime ( &t );
	if ( year < 2000 ) return 0;
	timeinfo->tm_year = year - 1900;
	timeinfo->tm_mon = month - 1;
	timeinfo->tm_mday = day;
	return mktime ( timeinfo );
}
bool genMoreThan2DaysOld ( const string& date )
{
	time_t t = genMakeTimeFromDateString ( date );
	if ( t == 0 ) return false;
	return time ( 0 ) - t > 172800;
}
bool genMoreThanNDaysOld ( const string& date, int n )
{
	time_t t = genMakeTimeFromDateString2 ( date );
	if ( t == 0 ) return false;
	int nsec = n * 86400;
	return time ( 0 ) - t > nsec;
}
StringVector getYearList ( const string& startYear, bool reverseFlag )
{
	int curYear = atoi ( genCurrentYear ().c_str () );
	int stYear = atoi ( startYear.c_str () );
	StringVector sv;
	for ( int i = stYear ; i <= curYear ; i++ ) {
		sv.push_back ( gen_itoa ( i ) );
	}
	if ( reverseFlag ) reverse ( sv.begin (), sv.end () );
	return sv;
}
string secToMins ( const string& sec )
{
	double m = atof ( sec.c_str () );
	m /= 60;
	stringstream sstr;
	genPrint ( sstr, m, 3 );
	return sstr.str ();
}
