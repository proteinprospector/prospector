/******************************************************************************
*                                                                             *
*  Library    : libxcalibur                                                   *
*                                                                             *
*  Filename   : lx_raw.h                                                      *
*                                                                             *
*  Created    : December 20th 2006                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2006-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef __lx_raw_h
#define __lx_raw_h

#include <afxdisp.h>
#include <string>
#include <nr.h>

typedef struct _datapeak {
	double mz;
	double intensity;
} DataPeak;

class ThermoRawFile : public COleDispatchDriver {
	double timeWindowStart;
	double timeWindowEnd;
	int previousFirstScan;
	int previousLastScan;
	VARIANT varData;
	VARIANT varPeakFlags;
	SAFEARRAY FAR* psa;
	DataPeak* pDataPeaks;
	long nArraySize;

	long Open ( LPCTSTR szFileName );
	long Close ();
	long GetFileName ( BSTR* pbstrFileName );
	long GetFirstSpectrumNumber ( long* pnFirstSpectrum );
	long GetLastSpectrumNumber ( long* pnLastSpectrum );
	long SetCurrentController ( long nControllerType, long nControllerNumber );
	long GetFilterForScanNum ( long nScanNumber, BSTR* pbstrFilter );
	long GetMassListFromScanNum(long* pnScanNumber, LPCTSTR bstrFilter, long nIntensityCutoffType, long nIntensityCutoffValue, long nMaxNumberOfPeaks, long bCentroidResult, double* pdCentroidPeakWidth, VARIANT* pvarMassList, VARIANT* pvarPeakFlags, long* pnArraySize );
	long GetChroData(long nChroType1, long nChroOperator, long nChroType2, LPCTSTR bstrFilter, LPCTSTR bstrMassRanges1, LPCTSTR bstrMassRanges2, double dDelay, double* pdStartTime, double* pdEndTime, long nSmoothingType, long nSmoothingValue, VARIANT* pvarChroData, VARIANT* pvarPeakFlags, long* pnArraySize );
	long GetAverageMassList(long* pnFirstAvgScanNumber, long* pnLastAvgScanNumber, long* pnFirstBkg1ScanNumber, long* pnLastBkg1ScanNumber, long* pnFirstBkg2ScanNumber, long* pnLastBkg2ScanNumber, LPCTSTR bstrFilter, long nIntensityCutoffType, long nIntensityCutoffValue, long nMaxNumberOfPeaks, long bCentroidResult, double* pdCentroidPeakWidth, VARIANT* pvarMassList, VARIANT* pvarPeakFlags, long* pnArraySize);
	long GetAveragedMassSpectrum(long* pnFirstAvgScanNumber, long nScansToAverage, BOOL bCentroidResult, VARIANT* pvarChroData, VARIANT* pvarPeakFlags, long* pnArraySize );
	long GetAveragedLabelData(long* pnFirstAvgScanNumber, long nScansToAverage, VARIANT* pvarChroData, VARIANT* pvarPeakFlags, long* pnArraySize );
	long ScanNumFromRT(double dRT, long* pnScanNumber);
	long RTFromScanNum(long nScanNumber, double* pdRT);
	void makeXYData ( XYData& xyData, const DataPeak* pDataPeaks, long size, double startMass, double endMass );
	void freeResources ();
	void getSpectrum ( XYData& xyData, long scan, double startMass, double endMass );
	void getAveragedSpectrum ( XYData& xyData, IntVector scanList, double startMass, double endMass );
	StringVector getScanLines ( long startScan, long endScan );
	IntVector getScanNumberList ( const std::string& msmsInfo, double mOverZ );
	std::string getScanLine ( int scan );
	int getPrecursorScanNumber ( int scan );
	IntVector getPrecursorScanNumberList ( const IntVector& scanList );
	int getFullScanNumber ( int scan );
	IntVector getFullScanNumberList ( const IntVector& scanList );
public:
	ThermoRawFile () {}
	ThermoRawFile ( const std::string& filename, double timeWindowStart = 0.0, double timeWindowEnd = 0.0 );
	~ThermoRawFile ();
	void getXYData ( std::vector <XYData>& vXYData, const std::string& msmsInfo, double mOverZ, bool ms, bool precursorScan = false, double startMass = 0.0, double endMass = 0.0 );
	int getStartScan ();
	int getEndScan ();
	DoubleVector getTimes ();
	std::string getFilename ();
};
#endif /* ! __lx_raw_h */
