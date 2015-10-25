/******************************************************************************
*                                                                             *
*  Library    : libnrec                                                       *
*                                                                             *
*  Filename   : ln_mrqcof.cpp                                                 *
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
*  Copyright (2000-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <nr.h>
void mrqmin ( double* x, double* y, double* sig, int ndata, double* a, int ma, int* lista, int mfit, double** covar, double** alpha, double* chisq, void (*funcs)(double, double*, double*, double*, int ), double* alamda )
{
	int k,kk,j,ihit;
	static double *da,*atry,**oneda,*beta,ochisq;

	if (*alamda < 0.0) {
		oneda=nrmatrix(1,mfit,1,1);
		atry=nrvector(1,ma);
		da=nrvector(1,ma);
		beta=nrvector(1,ma);
		kk=mfit+1;
		for (j=1;j<=ma;j++) {
			ihit=0;
			for (k=1;k<=mfit;k++)
				if (lista[k] == j) ihit++;
			if (ihit == 0)
				lista[kk++]=j;
			else if (ihit > 1) nrerror("Bad LISTA permutation in MRQMIN-1");
		}
		if (kk != ma+1) nrerror("Bad LISTA permutation in MRQMIN-2");
		*alamda=0.001;
		mrqcof(x,y,sig,ndata,a,ma,lista,mfit,alpha,beta,chisq,funcs);
		ochisq=(*chisq);
	}
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=mfit;k++) covar[j][k]=alpha[j][k];
		covar[j][j]=alpha[j][j]*(1.0+(*alamda));
		oneda[j][1]=beta[j];
	}
	gaussj(covar,mfit,oneda,1);
	for (j=1;j<=mfit;j++)
		da[j]=oneda[j][1];
	if (*alamda == 0.0) {
		covsrt(covar,ma,lista,mfit);
		free_nrvector(beta,1,ma);
		free_nrvector(da,1,ma);
		free_nrvector(atry,1,ma);
		free_nrmatrix(oneda,1,mfit,1,1);
		return;
	}
	for (j=1;j<=ma;j++) atry[j]=a[j];
	for (j=1;j<=mfit;j++)
		atry[lista[j]] = a[lista[j]]+da[j];
	mrqcof(x,y,sig,ndata,atry,ma,lista,mfit,covar,da,chisq,funcs);
	if (*chisq < ochisq) {
		*alamda *= 0.1;
		ochisq=(*chisq);
		for (j=1;j<=mfit;j++) {
			for (k=1;k<=mfit;k++) alpha[j][k]=covar[j][k];
			beta[j]=da[j];
			a[lista[j]]=atry[lista[j]];
		}
	}
	else {
		*alamda *= 10.0;
		*chisq=ochisq;
	}
	return;
}
void mrqcof ( double* x, double* y, double* sig, int ndata, double* a, int ma, int* lista, int mfit, double** alpha, double* beta, double* chisq, void (*funcs)(double, double*, double*, double*, int ) )
{
	int k,j,i;
	double ymod,wt,sig2i,dy,*dyda;

	dyda=nrvector(1,ma);
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=j;k++) alpha[j][k]=0.0;
		beta[j]=0.0;
	}
	*chisq=0.0;
	for (i=1;i<=ndata;i++) {
		(*funcs)(x[i],a,&ymod,dyda,ma);
		sig2i=1.0/(sig[i]*sig[i]);
		dy=y[i]-ymod;
		for (j=1;j<=mfit;j++) {
			wt=dyda[lista[j]]*sig2i;
			for (k=1;k<=j;k++)
				alpha[j][k] += wt*dyda[lista[k]];
			beta[j] += dy*wt;
		}
		(*chisq) += dy*dy*sig2i;
	}
	for (j=2;j<=mfit;j++)
		for (k=1;k<=j-1;k++) alpha[k][j]=alpha[j][k];
	free_nrvector(dyda,1,ma);
}
void fgauss ( double x, double* a, double* y, double* dyda, int na )
{
	int i;
	double fac,ex,arg;
// a [i]................intensity
// a [i+1]..............mean
// a [i+2]..............width [sqrt(2)*stddev]

	*y=0.0;
	for (i=1;i<=na-1;i+=3)
	{
		arg=(x-a[i+1])/a[i+2];
		ex=exp(-arg*arg);
		fac=a[i]*ex*2.0*arg;
		*y += a[i]*ex;
		dyda[i]=ex;
		dyda[i+1]=fac/a[i+2];
		dyda[i+2]=fac*arg/a[i+2];
	}
}
