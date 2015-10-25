/******************************************************************************
*                                                                             *
*  Library    : libnrec                                                       *
*                                                                             *
*  Filename   : spline.cpp                                                    *
*                                                                             *
*  Created    :                                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2000-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <limits>
#include <nr.h>
using std::vector;

SplineInterpolator::SplineInterpolator ( const XYData& xyData, double yp1, double ypn ) :
	xa ( xyData.getXList () ), ya ( xyData.getYList () )
{
	na = xyData.size ();
	y2.resize ( na );
	int i,k;
	double p,qn,sig,un;

	DoubleVector u ( na-1 );
	if (yp1 > 0.99e30)
		y2[0]=u[0]=0.0;
	else {
		y2[0] = -0.5;
		u[0]=(3.0/(xa[1]-xa[0]))*((ya[1]-ya[0])/(xa[1]-xa[0])-yp1);
	}
	for (i=1;i<=na-2;i++) {
		sig=(xa[i]-xa[i-1])/(xa[i+1]-xa[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(ya[i+1]-ya[i])/(xa[i+1]-xa[i]) - (ya[i]-ya[i-1])/(xa[i]-xa[i-1]);
		u[i]=(6.0*u[i]/(xa[i+1]-xa[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(xa[na-1]-xa[na-2]))*(ypn-(ya[na-1]-ya[na-2])/(xa[na-1]-xa[na-2]));
	}
	y2[na-1]=(un-qn*u[na-2])/(qn*y2[na-2]+1.0);
	for (k=na-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];

	//for ( double x = xa [0] ; x < xa [na-1] ; x += xInc ) {
	//	XYData::add ( x, getY ( x ) );
	//}
}
double SplineInterpolator::getY ( double x ) const
{
	int klo,khi,k;
	double h,b,a;

	klo=0;
	khi=na-1;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) nrerror("Bad XA input to routine SPLINT");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	return ( a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2[klo]+(b*b*b-b)*y2[khi])*(h*h)/6.0 );
}
XYData averageScans ( const vector <XYData>& vxyd )
{
	XYData xyd;
	double minX = std::numeric_limits<double>::max();
	double maxX = -std::numeric_limits<double>::max();
	for ( std::vector <XYData>::size_type i = 0 ; i < vxyd.size () ; i++ ) {
		std::cout << vxyd [i].minX () << " " << vxyd [i].maxX () << "<br />" << std::endl;
		minX = genMin ( minX, vxyd [i].minX () );
		maxX = genMin ( maxX, vxyd [i].maxX () );
		std::cout << vxyd [i].minX () << " " << vxyd [i].maxX () << "<br />" << std::endl;
	}
	for ( std::vector <XYData>::size_type j = 0 ; j < vxyd.size () ; j++ ) {
		SplineInterpolator si ( vxyd [j] );
	}
	return xyd;
}
