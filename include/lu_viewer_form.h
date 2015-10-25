/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_viewer_form.h                                              *
*                                                                             *
*  Created    : March 8th 2011                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2011-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_viewer_form_h
#define __lu_viewer_form_h

#include <ostream>
#include <lu_srch_form.h>
#include <lu_form_valid.h>
#include <lu_prod_form.h>

class MSViewerForm : public ProspectorForm {
	bool createDataset;
	FormValidatingJavascript fvj;
	IonTypeForm ionTypeForm;
	SearchForm searchForm;
	MSMSToleranceForm msmsToleranceForm;
	VariableModsForm variableModsForm;
	MassModificationForm massModificationForm;
	ExtraUserModsForm extraUserModsForm;
	ExtraUserModsForm extraUserModsForm2;
	CrosslinkingForm crosslinkingForm;
	MatrixModeForm matrixModeForm;
	bool filterPeaks;
	int numCols;
	void printJavascriptFunctions ( std::ostream& os );
	void printHTMLFirstTable ( std::ostream& os );
	void printHTMLSecondTable ( std::ostream& os );
	void setModVisualizationFlags ( const std::string& val, bool& div_cm, bool& div_cm_col, bool& div_vm_col, bool& div_am_col ) const;
	void setLinkSearchVisualizationFlags ( const std::string& val, bool& div_ls );
	void setResultsFileFormatVisualizationFlags ( const std::string& val, bool& div_info_1, bool& div_info_2, bool& div_info_3, bool& div_info_4, bool& div_info_5, bool& div_info_6, bool& div_info_7, bool& div_info_8, bool& div_info_9, bool& div_info_10, bool& div_info_11 ) const;
	static const int maxSortLevels;
	static const int maxFilterLevels;
	static const int maxRemoveReplicateLevels;
protected:
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	MSViewerForm ( const VectorConstParameterListPtr& params, bool createDataset = true, int numCols = 100 );
	virtual void printHTML ( std::ostream& os );
	void printHTMLForm ( std::ostream& os, bool saveParamsOption, bool initialForm );
	void printTagParamsCGI ( std::ostream& os ) const;
	void printTagParamsHTMLJavascriptHidden ( std::ostream& os ) const;
	static int getMaxSortLevels () { return maxSortLevels; }
	static int getMaxFilterLevels () { return maxFilterLevels; }
	static int getMaxRemoveReplicateLevels () { return maxRemoveReplicateLevels; }
};

struct ViewerFileConverter {
	std::string script;
	int numTitleLines;
	int numHeaderLines;
	std::string delimiter;
	std::string spectrumIdentifier;
	std::string fraction;
	std::string scanID;
	std::string peptide;
	std::string charge;
	std::string modifications;
};
typedef std::map <std::string, ViewerFileConverter> MapStringToViewerFileConverter;
typedef MapStringToViewerFileConverter::iterator MapStringToViewerFileConverterIterator;
typedef MapStringToViewerFileConverter::const_iterator MapStringToViewerFileConverterConstIterator;

class ViewerResultsFileConverter {
	MapStringToViewerFileConverter converters;
	ViewerResultsFileConverter ();
	void init ();
public:
	static ViewerResultsFileConverter& instance ();
	StringVector getConverterList ();
	std::string getScript ( const std::string& name ) const;
	void getScriptParameters ( const std::string& name, int& numTitleLines, int& numHeaderLines, std::string& delimiter ) const;
	void getScriptParameters2 ( const std::string& name, std::string& spectrumIdentifier, std::string& fraction, std::string& scanID, std::string& peptide, std::string& charge, std::string& modifications ) const;
	StringVector getNameList () const;
};

#endif /* ! __lu_viewer_form_h */
