/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_app_gr.h                                                   *
*                                                                             *
*  Created    : October 23rd 1997                                             *
*                                                                             *
*  Purpose    : Java graph applet interface functions.                        *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1997-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_app_gr_h
#define __lu_app_gr_h

#include <string>
#include <ostream>
#include <sstream>
#include <nr.h>

class GraphData : public XYData {
protected:
	static int xPrecision;
	static int yPrecision;
	void writeRangeParams ( std::ostream& os, double xMin, double xMax ) const;
public:
	GraphData () {}
	GraphData ( const XYData& xyData );
	virtual void draw ( std::ostream& os, double xMin = 0.0, double xMax = 0.0 ) const;
	virtual void writeXML ( std::ostream& os ) const;
	virtual void writeTabDelimitedText ( std::ostream& os ) const;
	static void setXPrecision ( int p ) { xPrecision = p; }
	static void setYPrecision ( int p ) { yPrecision = p; }
};

class FileGraphData : public GraphData {
	bool retainZeroPeaks;
public:
	FileGraphData () {}
	FileGraphData ( const XYData& xyData, bool retainZeroPeaks );
	virtual void draw ( std::ostream& os, double xMin = 0.0, double xMax = 0.0 ) const;
};
class ColouredGraphData : public GraphData {
	IntVector colourList;
public:
	ColouredGraphData () {}
	ColouredGraphData ( const XYData& xyData, int colour = 0 );
	void add ( const double& x, const double& y, const int& colour = 0 );
	virtual void draw ( std::ostream& os, double xMin = 0.0, double xMax = 0.0 ) const;
};

class LabelledGraphData : public GraphData {
protected:
	IntVector colourList;
	StringVector labelList;
public:
	LabelledGraphData () {}
	void add ( const double& x, const double& y, const std::string& label = "", const int& colour = 0 );
	virtual void draw ( std::ostream& os, double xMin = 0.0, double xMax = 0.0 ) const;
};

class LabelledCatagorizedGraphData : public LabelledGraphData {
protected:
	StringVector categories;
	BoolDeque showCategory;
	StringVector catagoryList;
	void writeCategories ( std::ostream& os ) const;
public:
	LabelledCatagorizedGraphData ( const StringVector& categories, const BoolDeque& showCategory );
	void add ( const double& x, const double& y, const std::string& label = "", const int& colour = 0, const std::string& category = "" );
	virtual void draw ( std::ostream& os, double xMin = 0.0, double xMax = 0.0 ) const;
	void defineCheckboxes ( std::ostream& os ) const;
};
class LabelledCatagorizedFileGraphData : public LabelledCatagorizedGraphData {
	bool retainZeroPeaks;
	int fileGraphSize;
public:
	LabelledCatagorizedFileGraphData ( const StringVector& categories, const BoolDeque& showCategory, const XYData& xyData, bool retainZeroPeaks );
	void draw ( std::ostream& os, double xMin, double xMax ) const;
};

class SpectrumGraph {
	struct AppletFont {
		std::string family;
		std::string style;
		int points;
	};
	struct AppletColor {
		int red;
		int green;
		int blue;
	};
	std::string xAxisLabel;
	int lineWidth;
	AppletColor appletBackgroundColor;
	AppletColor axesColor;
	AppletColor defaultPeakColor;
	typedef std::vector <AppletColor> AppletColorVector;
	typedef AppletColorVector::size_type AppletColorVectorSizeType;
	AppletColorVector applicationColors;
	int appletWidth;
	int appletHeight;
	AppletFont defaultFont;
	AppletFont peakLabelFont;
	void readParams ( const std::string& fileName );
	static bool drawGraphFlag;
public:
	SpectrumGraph ( const std::string& fileName );
	void drawGraph ( std::ostream& os, GraphData& graphData, bool sorted = true, double xMin = 0.0, double xMax = 0.0, double heightMultipler = 1.0 );
	void drawGraphData ( std::ostream& os, GraphData& graphData );
	static void setDrawGraph ( bool flag );
	static bool getDrawGraph ();
};

#endif /* ! __lu_app_gr_h */
