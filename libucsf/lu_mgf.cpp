/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mgf.cpp                                                    *
*                                                                             *
*  Created    : April 20th 2007                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif
#include <lg_string.h>
#include <lu_mgf.h>
#include <lu_getfil.h>
#include <lu_xml.h>
using std::string;
using std::make_pair;

/*

Examples of getMSMSInfo:

MGFInstance					- empty string
ExtractMZXMLToMGFInstance	- single scan number
DistillerMGFInstance		- single scan number or string like 2:2852,2857 2=number of scans 2852 to 2857 is start and end of scan range
DtaToMGFInstance			- single scan number or something like 2852,2857 these are the start and end scan
Analyst2MGFInstance			- 1550,2 or 1464,2:1454,3 or 3136,2:3138,2:3147,2 these numbers are (cycle number 1),(experiment number 1):(cycle number 2),(experiment number 2)
AnalystDistillerMGFInstance	- If up to 2 scans are summed it returns the same type of string as above. Otherwise an s is returned.

Examples of getMSMSInfoFromScans:

MGFInstance					- empty string
ExtractMZXMLToMGFInstance	- empty string
DistillerMGFInstance		- empty string
DtaToMGFInstance			- empty string
Analyst2MGFInstance			- empty string
AnalystDistillerMGFInstance	- This function returns an msmsInfo string from the scan numbers.
                              The scan numbers are stored in a SCANS= line. This is invoked if getMSMSInfo returns as s.

Note for TOFTOF data using Protein Prospector format files the msmsInfo string contains the job run item
*/

namespace {
string getMascotParameter ( const string& line, const string& param, char c )
{
	int start = line.find ( param );
	if ( start == string::npos ) return "";
	start += param.length ();
	int end = line.find ( c, start );
	if ( end == string::npos ) return "";
	return line.substr ( start, end - start );
}
string getMascotParameter ( const string& line, const string& param, const string& param2 )
{
	int start = line.find ( param );
	if ( start == string::npos ) return "";
	start += param.length ();
	int end = line.find ( param2, start );
	if ( end == string::npos ) return "";
	return line.substr ( start, end - start );
}
string getMascotParameter ( const string& line, const string& param )
{
	int start = line.find ( param );
	if ( start == string::npos ) return "";
	else return line.substr ( start + param.length () );
}
}
class MGFInstance {
	StringVector start;
	StringVector contains;
	StringVector end;
	string spotStart;
	string spotEnd;
	string runStart;
	string runEnd;
public:
	MGFInstance ( const string& s );
	virtual bool isMGFType ( const string& tLine );
	virtual string getMSMSInfo ( const string& title ) const { return ""; }
	virtual string getMSMSInfoFromScans ( const IntVector& scans ) const { return ""; }
	virtual string getSpot ( const string& title, const string& defaultSpot ) const;
	int getRun ( const string& title ) const;
	bool isSpottingPlate () const;
};
/*
Different Title line formats:
Old Analyst Mascot dll:
TITLE=Elution from: 1.66 to 1.66   period: 0   experiment: 1 cycles:  1  (Charge not auto determined)
TITLE=Elution from: 16.73 to 16.73   period: 0   experiment: 1 cycles:  1  
TITLE=Elution from: 18.94 to 19.7   period: 0   experiments: 1 to 2 cycles:  3  
Agilent 1:
TITLE= +MS2(409.8)
Agilent 2:
TITLE= Cmpd 1, +MSn(531.5), 12.6 min
LCQ_UCSF:
TITLE= (scan=43) (precursorMz=369.1222229) (precursorIntensity=1165.47) (rt=69.8091) (basePeakMz=351.266) (basePeakIntensity=8344.86) (scanType=Full) (msLevel=2) (rawDataFileName=F5081106.RAW)
TOFTOF:
TITLE=Label: 59, Spot_Id: 41307, Peak_List_Id: 3210, MSMS Job_Run_Id: 14369, Comment: 
NIST:
TITLE=Scan:2363 RT:55.55 Filter:+ c d Full ms2 1081.35@35.00 [ 285.00-2000.00]
*/
MGFInstance::MGFInstance ( const string& s ) :
	start ( XMLParser::getStringVectorValue ( s, "start" ) ),
	contains ( XMLParser::getStringVectorValue ( s, "contains" ) ),
	end ( XMLParser::getStringVectorValue ( s, "end" ) ),
	spotStart ( XMLParser::getStringValue ( s, "spot_start", "" ) ),
	spotEnd ( XMLParser::getStringValue ( s, "spot_end", "" ) ),
	runStart ( XMLParser::getStringValue ( s, "run_start", "" ) ),
	runEnd ( XMLParser::getStringValue ( s, "run_end", "" ) )
{
}
bool MGFInstance::isMGFType ( const string& tLine )
{
	if ( !start.empty () ) {
		bool flag1 = false;									// Check the start possibilities
		for ( StringVectorSizeType i = 0 ; i < start.size () ; i++ ) {
			string& s = start [i];
			if ( !tLine.compare ( 0, s.length (), s ) ) {
				flag1 = true;
				break;
			}
		}
		if ( !flag1 ) return false;
	}
	if ( !end.empty () ) {
		bool flag2 = false;									// Check the end possibilities
		string::size_type tlen = tLine.length ();
		for ( StringVectorSizeType j = 0 ; j < end.size () ; j++ ) {
			string& e = end [j];
			string::size_type elen = e.length ();
			if ( tlen >= elen ) {
				if ( !tLine.compare ( tlen - elen, elen, e ) ) {
					flag2 = true;
					break;
				}
			}
		}
		if ( !flag2 ) return false;
	}
	if ( !contains.empty () ) {
		string::size_type idx = 0;
		for ( StringVectorSizeType k = 0 ; k < contains.size () ; k++ ) {
			idx = tLine.find ( contains [k], idx );
			if ( idx == string::npos ) return false;
			idx++;
		}
	}
	return true;
}
string MGFInstance::getSpot ( const string& title, const string& defaultSpot ) const
{
	if ( !spotStart.empty () ) {
		if ( !spotEnd.empty () )	return getMascotParameter ( title, spotStart, spotEnd );
		else						return getMascotParameter ( title, spotStart );
	}
	else
		return defaultSpot;
}
int MGFInstance::getRun ( const string& title ) const
{
	if ( !runStart.empty () )
		return atoi ( ( getMascotParameter ( title, runStart, runEnd ) ).c_str () );
	else
		return 1;
}
bool MGFInstance::isSpottingPlate () const
{
	return !runStart.empty ();
}
class ExtractMZXMLToMGFInstance : public MGFInstance {
public:
	ExtractMZXMLToMGFInstance ( const string& s ) :
		MGFInstance ( s ) {}
	string getMSMSInfo ( const string& title ) const;
};
string ExtractMZXMLToMGFInstance::getMSMSInfo ( const string& title ) const
{
//TITLE= (scan=3) (precursorMz=325.1699219) (precursorIntensity=10573.5) (rt=1.7016) (basePeakMz=304.26) (basePeakIntensity=50625.7) (scanType=Full) (msLevel=2) (rawDataFileName=F5080208.RAW)
	return getMascotParameter ( title, "scan=", ')' );
}
class DataConverter5600Instance : public MGFInstance {
public:
	DataConverter5600Instance ( const string& s ) :
		MGFInstance ( s ) {}
	string getSpot ( const string& title, const string& defaultSpot ) const;
};
string DataConverter5600Instance::getSpot ( const string& title, const string& defaultSpot ) const
{
//TITLE=Locus:1.1.1.320.2
//TITLE=Locus:1.1.1.1624.2 File:"120118ry_201B7-32_1_1.wiff"
	string::size_type idx3 = title.find ( '.' ) + 1;
	string::size_type idx2 = title.find ( '.', idx3 ) + 1;
	string::size_type idx1 = title.find ( '.', idx2 ) + 1;
	string::size_type idx = title.find ( ' ', idx1 );
	if ( idx == string::npos )
		return title.substr ( idx1 );
	else
		return title.substr ( idx1, idx-idx1 );
}
class DataConverter5600mzMLToMGFInstance : public MGFInstance {
public:
	DataConverter5600mzMLToMGFInstance ( const string& s ) :
		MGFInstance ( s ) {}
	string getSpot ( const string& title, const string& defaultSpot ) const;
};
string DataConverter5600mzMLToMGFInstance::getSpot ( const string& title, const string& defaultSpot ) const
{
//TITLE=sample=1 period=1 cycle=327 experiment=2
	string cycle = getMascotParameter ( title, "cycle=", ' ' );
	string expt = getMascotParameter ( title, "experiment=" );
	if ( expt [expt.length ()-1] == '\r' ) expt = expt.substr ( 0, expt.length ()-1 );
	return cycle + "." + expt;
}
class DistillerMGFInstance : public MGFInstance {
public:
	DistillerMGFInstance ( const string& s ) :
		MGFInstance ( s ) {}
	string getMSMSInfo ( const string& title ) const;
};
string DistillerMGFInstance::getMSMSInfo ( const string& title ) const
{
//TITLE=Scan 8678 (rt=63.6507) [C:\Users\robert\serum_PP Mascot compare\C5053103.RAW]
//TITLE=Sum of 2 scans in range 2852 (rt=44.0754) to 2857 (rt=44.1502) [\\Ltqft\RAW\FT_Mar2006\F6032015.RAW]
//TITLE=1: Scan 195 (rt=5.6608) [d:\cygwin\home\lynna\robert_raw\F7010302.RAW]
	int start = title.find ( 'S' );
	if ( start == string::npos ) return "";
	char type = title [start+1];
	if ( type == 'c' ) return getMascotParameter ( title, "Scan " , ' ' );
	if ( type == 'u' ) {
		string specInfo;
		specInfo += getMascotParameter ( title, "of ", ' ' );
		specInfo += ':';
		specInfo += getMascotParameter ( title, "range ", ' ' );
		specInfo += ',';
		specInfo += getMascotParameter ( title, "to ", ' ' );
		return specInfo;
	}
	return "";
}
class DtaToMGFInstance : public MGFInstance {
public:
	DtaToMGFInstance ( const string& s ) :
		MGFInstance ( s ) {}
	string getMSMSInfo ( const string& title ) const;
	string getSpot ( const string& title, const string& defaultSpot ) const;
	bool isMGFType ( const string& tLine );
};
string DtaToMGFInstance::getMSMSInfo ( const string& title ) const
{
//TITLE=F6032010.1093.1095.2.dta
//TITLE=D100930_yeast_SCX10S_rak_ft8E_pc_01.2.2.2
//idx1 is the '.' after the fraction
//idx2 is the '.' after the start scan
//idx3 is the '.' after the end scan
	string::size_type idx4 = title.length ();
	string::size_type idx3 = title.rfind ( '.', idx4-1 );
	string endBit = title.substr ( idx3+1, idx4-idx3-1 );

	if ( !genStringIsInteger ( endBit ) ) {		// Must be something like .dta at the end, skip the charge field
		idx3 = title.rfind ( '.', idx3-1 );
	}
	string::size_type idx2 = title.rfind ( '.', idx3-1 );
	string::size_type idx1 = title.rfind ( '.', idx2-1 );
	string start = title.substr ( idx1+1, idx2-idx1-1 );
	if ( start.length () && start [0] == '0' ) {	// deal with leading zeros
		int s = atoi ( start.c_str () );
		start = gen_itoa ( s );
	}
	string end = title.substr ( idx2+1, idx3-idx2-1 );
	if ( end.length () && end [0] == '0' ) {	// deal with leading zeros
		int e = atoi ( end.c_str () );
		end = gen_itoa ( e );
	}
	if ( start == end ) return start;
	else return start + "," + end;
}
string DtaToMGFInstance::getSpot ( const string& title, const string& defaultSpot ) const
{
	string::size_type idx3 = title.length ();
	string::size_type idx2 = title.rfind ( '.', idx3-1 );
	string endBit = title.substr ( idx2+1, idx3-idx2-1 );
	if ( !genStringIsInteger ( endBit ) ) {		// Must be something like .dta at the end, skip the charge field
		idx2 = title.rfind ( '.', idx2-1 );
	}
	idx2 = title.rfind ( '.', idx2-1 );
	string::size_type idx1 = title.rfind ( '.', idx2-1 );
	return title.substr ( idx1+1, idx2-idx1-1 );
}
//TITLE=F6032010.1093.1095.2.dta
//TITLE=D100930_yeast_SCX10S_rak_ft8E_pc_01.2.2.2
bool DtaToMGFInstance::isMGFType ( const string& tLine )
{
	string::size_type idx1 = tLine.find ( '.' );			// Find first dot
	if ( idx1 == string::npos ) return false;
	string::size_type idx2 = tLine.find ( '.', idx1+1 );	// Find second dot
	if ( idx2 == string::npos ) return false;
	if ( !genStringIsInteger ( tLine.substr ( idx1+1, idx2-idx1-1 ) ) ) return false;
	string::size_type idx3 = tLine.find ( '.', idx2+1 );	// Find third dot
	if ( idx3 == string::npos ) return false;
	if ( !genStringIsInteger ( tLine.substr ( idx2+1, idx3-idx2-1 ) ) ) return false;
	string::size_type idx4 = tLine.find ( '.', idx3+1 );	// Find fourth dot
	string charge;
	if ( idx4 == string::npos ) charge = tLine.substr ( idx3+1 );
	else						charge = tLine.substr ( idx3+1, idx4-idx3-1 );
	if ( !genStringIsInteger ( charge ) ) return false;
	return true;
}
class Analyst2MGFInstance : public MGFInstance {
public:
	Analyst2MGFInstance ( const string& s ) :
		MGFInstance ( s ) {}
	string getMSMSInfo ( const string& title ) const;
};
string Analyst2MGFInstance::getMSMSInfo ( const string& title ) const
{
//TITLE=File: F15uLUCSF.wiff, Sample: F1 25-26_5001 (sample number 1), Elution: 1.655 min, Period: 1, Cycle(s): 97 (Experiment 2)
//TITLE=File: F25uLUCSF.wiff, Sample: F2 26_5-28002 (sample number 1), Elution: 16.781 min, Period: 1, Cycle(s): 976 (Experiment 2) (Charge not auto determined)
//TITLE=File: F25uLUCSF.wiff, Sample: F2 26_5-28002 (sample number 1), Elution: 25.987 to 26.554 min, Period: 1, Cycle(s): 1124 (Experiment 2), 1119 (Experiment 4)
//TITLE=File: F25uLUCSF.wiff, Sample: F2 26_5-28002 (sample number 1), Elution: 26.691 to 28.368 min, Period: 1, Cycle(s): 1125, 1129, 1139, 1150 (Experiment 2)
//TITLE=File: F25uLUCSF.wiff, Sample: F2 26_5-28002 (sample number 1), Elution: 26.813 to 28.437 min, Period: 1, Cycle(s): 1129, 1139, 1150 (Experiment 3), 1125 (Experiment 4)
//TITLE=File: q040907003.wiff.tmp, Sample: Hela SILAC (sample number 1), Elution: 49.85 min, Period: 1, Cycle(s): 2678 (Experiment 3)
	string specInfo;
	string info = getMascotParameter ( title, "Cycle(s): " );
	int start = 0;
	for ( ; ; ) {
		int end = info.find ( "),", start );
		string cur = info.substr ( start, end-start );
		if ( cur [cur.length ()-1] != ')' ) cur += ')';
		string expt = getMascotParameter ( cur, "Experiment ", ')' );
		string cycles = cur.substr ( 0, cur.find ( " (" ) );
		int start2 = 0;
		for ( ; ; ) {
			int end2 = cycles.find ( ",", start2 );
			string cy = cycles.substr ( start2, end2 - start2 );
			int end3 = cy.find ( '-' );
			if ( end3 == string::npos ) {
				specInfo += cy + ',' + expt + ':';
			}
			else {
				int startCycle = atoi ( cy.substr ( 0, end3 ).c_str () );
				int endCycle = atoi ( cy.substr ( end3+1 ).c_str () );
				for ( int i = startCycle ; i <= endCycle ; i++ ) {
					specInfo += gen_itoa ( i ) + ',' + expt + ':';
				}
			}
			if ( end2 == string::npos ) break;
			start2 = end2 + 2;
		}
		if ( end == string::npos ) break;
		start = end + 3;
	}
	return specInfo.substr ( 0, specInfo.length () - 1 );
}
class AnalystDistillerMGFInstance : public MGFInstance {
	mutable int expPerCycle;
	mutable int cycle1;
	mutable int cycle2;
	mutable int expt1;
	mutable int expt2;
public:
	AnalystDistillerMGFInstance ( const string& s ) :
		MGFInstance ( s ) {}
	string getMSMSInfo ( const string& title ) const;
	string getMSMSInfoFromScans ( const IntVector& scans ) const;
};
string AnalystDistillerMGFInstance::getMSMSInfo ( const string& title ) const
{
//TITLE=1: Scan 5 (rt=4.106, p=0, c=1, e=1) [C:\MSDATA\QS20060131_S_18mix_02.wiff]
//TITLE=13: Sum of 2 scans in range 1130 (rt=436.21, p=0, c=376, e=1) to 1133 (rt=439.385, p=0, c=377, e=1) [C:\MSDATA\QS20060131_S_18mix_02.wiff]
//TITLE=109: Sum of 3 scans in range 1602 (rt=859.418, p=0, c=533, e=2) to 1613 (rt=878.393, p=0, c=537, e=1) [C:\MSDATA\QS20060131_S_18mix_02.wiff]
	int start = title.find ( 'S' );
	if ( start == string::npos ) return "";
	char type = title [start+1];
	if ( type == 'c' ) {
		cycle1 = atoi ( getMascotParameter ( title, "c=", ',' ).c_str () ) + 1;
		expt1 = atoi ( getMascotParameter ( title, "e=", ')' ).c_str () ) + 1;
		return gen_itoa ( cycle1 ) + "," + gen_itoa ( expt1 );
	}
	if ( type == 'u' ) {
		cycle1 = atoi ( getMascotParameter ( title, "c=" , ',' ).c_str () ) + 1;
		expt1 = atoi ( getMascotParameter ( title, "e=" , ')' ).c_str () ) + 1;
		string title2 = title.substr ( title.find ( "to" ) );
		cycle2 = atoi ( getMascotParameter ( title2, "c=", ',' ).c_str () ) + 1;
		expt2 = atoi ( getMascotParameter ( title2, "e=", ')' ).c_str () ) + 1;
		if ( !title.compare ( start+7, 2, "2 " ) ) {	// If 2 scans summed
			return gen_itoa ( cycle1 ) + "," + gen_itoa ( expt1 ) + ":" + gen_itoa ( cycle2 ) + "," + gen_itoa ( expt2 );
		}
		else {
			int startScan = atoi ( getMascotParameter ( title, "range ", ' ' ).c_str () );
			int endScan = atoi ( getMascotParameter ( title, "to ", ' ' ).c_str () );
			startScan -= expt1-2;
			endScan -= expt2-2;
			expPerCycle = ( endScan - startScan ) / ( cycle2 - cycle1 );
			return "s";
		}
	}
	return "";
}
string AnalystDistillerMGFInstance::getMSMSInfoFromScans ( const IntVector& scans ) const
{
	string ret = gen_itoa ( cycle1 ) + "," + gen_itoa ( expt1 ) + ":";
	for ( IntVectorSizeType i = 1 ; i < scans.size () - 1 ; i++ ) {
		int diff = scans [i] - scans [0];
		int nCycles = diff / expPerCycle;
		int remainder = diff % expPerCycle;
		int cycle = cycle1 + nCycles;
		int expt;
		if ( remainder + expt1 > expPerCycle ) {
			cycle += 1;
			expt = remainder + expt1 - expPerCycle;
		}
		else
			expt = remainder + expt1;

		ret += gen_itoa ( cycle ) + "," + gen_itoa ( expt ) + ":";
	}
	ret += gen_itoa ( cycle2 ) + "," + gen_itoa ( expt2 );
	return ret;
}
class ProteomeDiscovererMGFInstance : public MGFInstance {
public:
	ProteomeDiscovererMGFInstance ( const string& s ) :
		MGFInstance ( s ) {}
	string getMSMSInfo ( const string& title ) const;
};
string ProteomeDiscovererMGFInstance::getMSMSInfo ( const string& title ) const
{
//TITLE=File: "D:\RAW_DATA\U_Zaman\XL_peptide_T0.raw"; SpectrumID: "1"; scans: "102"
//TITLE=File1694 Spectrum1 scans: 2
	string scan = getMascotParameter ( title, "scans: " );
	if ( scan [scan.length ()-1] == '\r' )	scan = scan.substr ( 0, scan.length ()-1 );
	if ( scan [scan.length ()-1] == '"' )	scan = scan.substr ( 0, scan.length ()-1 );
	if ( scan [0] == '"' )					scan = scan.substr ( 1 );
	return scan;
}
class MSConvertMGFInstance : public MGFInstance {
public:
	MSConvertMGFInstance ( const string& s ) :
		MGFInstance ( s ) {}
	string getMSMSInfo ( const string& title ) const;
};
string MSConvertMGFInstance::getMSMSInfo ( const string& title ) const
{
//TITLE=V20140213-07.7.7.3 File:"V20140213-07.raw", NativeID:"controllerType=0 controllerNumber=1 scan=7"
	return getMascotParameter ( title, "scan=" , '"' );
}
class MSConvert2MGFInstance : public MGFInstance {
public:
	MSConvert2MGFInstance ( const string& s ) :
		MGFInstance ( s ) {}
	string getMSMSInfo ( const string& title ) const;
};
string MSConvert2MGFInstance::getMSMSInfo ( const string& title ) const
{
//TITLE=689RM_Mix_1.345173.345173.3 File:"689RM_Mix_1.d", NativeID:"scanId=345173"
	return getMascotParameter ( title, "scanId=" , '"' );
}
MGFInfo::MGFInfo () :
	currentInstancePtr ( 0 )
{
	char* info = getFileInfo ( MsparamsDir::instance ().getParamPath ( "mgf.xml" ), '\n', 1, false );
	StringVector sv = XMLParser::getStringVectorValue ( info, "mgf_type" );
	delete [] info;
	for ( StringVectorSizeType i = 0 ; i < sv.size () ; i++ ) {
		string& s = sv [i];
		string name = XMLParser::getStringValue ( s, "name" );
		if ( name == "ANALYST_DLL_2" ) 
			mgfi.push_back ( make_pair ( name, new Analyst2MGFInstance ( s ) ) );
		else if ( name == "DISTILLER" )
			mgfi.push_back ( make_pair ( name, new DistillerMGFInstance ( s ) ) );
		else if ( name == "ANALYST_DISTILLER" )
			mgfi.push_back ( make_pair ( name, new AnalystDistillerMGFInstance ( s ) ) );
		else if ( name == "DTA_TO_MGF" )
			mgfi.push_back ( make_pair ( name, new DtaToMGFInstance ( s ) ) );
		else if ( name == "EXTRACTMZXMLTOMGF" )
			mgfi.push_back ( make_pair ( name, new ExtractMZXMLToMGFInstance ( s ) ) );
		else if ( name == "5600_MS_DATA_CONVERTER" )
			mgfi.push_back ( make_pair ( name, new DataConverter5600Instance ( s ) ) );
		else if ( name == "5600_MZML_TO_MGF_MS_DATA_CONVERTER" )
			mgfi.push_back ( make_pair ( name, new DataConverter5600mzMLToMGFInstance ( s ) ) );
		else if ( name == "MS_CONVERT" )
			mgfi.push_back ( make_pair ( name, new MSConvertMGFInstance ( s ) ) );
		else if ( name == "MS_CONVERT_2" )
			mgfi.push_back ( make_pair ( name, new MSConvert2MGFInstance ( s ) ) );
		else if ( name == "PROTEOME_DISCOVERER" )
			mgfi.push_back ( make_pair ( name, new ProteomeDiscovererMGFInstance ( s ) ) );
		else
			mgfi.push_back ( make_pair ( name, new MGFInstance ( s ) ) );
	}
	mgfi.push_back ( make_pair ( string("default"), new MGFInstance ( "" ) ) );
}
MGFInfo::~MGFInfo ()
{
	for ( VectorPairStringMGFInstancePtrSizeType i = 0 ; i < mgfi.size () ; i++ ) {
		delete mgfi [i].second;
	}
}
MGFInfo& MGFInfo::instance ()
{
	static MGFInfo m;
	return m;
}
MGFInstance* MGFInfo::getMGFInstancePtr ( const string& title )
{
	string tLine = gen_strtrim2 ( title.substr ( 6 ) );
	for ( VectorPairStringMGFInstancePtrSizeType i = 0 ; i < mgfi.size () ; i++ ) {
		if ( mgfi [i].second->isMGFType ( tLine ) ) {
			return mgfi [i].second;
		}
	}
	return 0;
}
bool MGFInfo::getTitleParams ( const string& line, string& spot, int& run, string& msmsInfo )
{
	bool comp = ( line [0] == 'T' );
	if ( comp ) {
		if ( currentInstancePtr == 0 ) {
			currentInstancePtr = getMGFInstancePtr ( line );
		}
		run = currentInstancePtr->getRun ( line );
		spot = currentInstancePtr->getSpot ( line, spot );
		msmsInfo = currentInstancePtr->getMSMSInfo ( line );
	}
	return comp;
}
string MGFInfo::getMSMSInfoFromScans ( const IntVector& scans ) const
{
	return currentInstancePtr->getMSMSInfoFromScans ( scans );
}
bool MGFInfo::isSpottingPlate ( const string& line, bool& spottingPlate )
{
	bool comp = ( line [0] == 'T' );
	if ( comp ) {
		if ( currentInstancePtr == 0 ) {
			currentInstancePtr = getMGFInstancePtr ( line );
		}
		spottingPlate = currentInstancePtr->isSpottingPlate ();
	}
	return comp;
}
void MGFInfo::reset ()
{
	currentInstancePtr = 0;
}
