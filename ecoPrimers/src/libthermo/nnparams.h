/*
 *  nnparams.h
 *  PHunterLib
 *
 * Nearest Neighbor Model Parameters
 *
 *  Created by Tiayyba Riaz on 02/07/09.
 *
 */

#ifndef NNPARAMS_H_
#define NNPARAMS_H_

#include <math.h>
#include <string.h>
//#include "../libecoprimer/ecoprimer.h"

// following defines to simplify coding...
#define ndH(a,b,c,d) nparm->dH[a][b][c][d]
#define ndS(a,b,c,d) nparm->dS[a][b][c][d]
#define forbidden_enthalpy 1000000000000000000.0f
#define	R 1.987f
#define	SALT_METHOD_SANTALUCIA 1
#define	SALT_METHOD_OWCZARZY 2

#define DEF_CONC_PRIMERS 0.0000008
#define DEF_CONC_SEQUENCES 0
#define DEF_SALT 0.05

#define GETNUMCODE(a) bpencoder[a - 'A']
#define GETREVCODE(a) 5-bpencoder[a - 'A']


extern double forbidden_entropy;

static char bpencoder[] = {	1,							// A
		                   	0,                          // b
		                   	2,                          // C
		                   	0,0,0,                      // d, e, f
		                   	3,                          // G
		                   	0,0,0,0,0,0,0,0,0,0,0,0,    // h,i,j,k,l,m,n,o,p,q,r,s
		                   	4,0,                        // T,U
		                   	0,0,0,0,0};                 // v,w,x,y,z


typedef struct CNNParams_st
{
	double Ct1;
	double Ct2;
	double rlogc;
	double kplus;
	double kfac;
	int saltMethod;
	double gcContent;
	double new_TM;
	double dH[6][6][6][6];  // A-C-G-T + gap + initiation (dangling end, $ sign)
	double dS[6][6][6][6];
}CNNParams, * PNNParams;

void nparam_InitParams(PNNParams nparm, double c1, double c2, double kp, int sm);
int nparam_CountGCContent(char * seq );
double nparam_GetEntropy(PNNParams nparm, char x0, char x1, char y0, char y1);
double nparam_GetEnthalpy(PNNParams nparm, char x0, char x1, char y0, char y1);
double nparam_CalcTM(double entropy,double enthalpy);
double nparam_CalcSelfTM(PNNParams nparm, char* seq, int len);
double nparam_CalcTwoTM(PNNParams nparm, char* seq1, char* seq2, int len);

double nparam_GetInitialEntropy(PNNParams nparm) ;
double calculateMeltingTemperatureBasic (char * seq);
//void getThermoProperties (ppair_t* pairs, size_t count, poptions_t options);

#endif
