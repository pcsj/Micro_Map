#define _CRT_SECURE_NO_WARNINGS

#define DEBUG
//#define TEST_OPTICAL_FUNCTIONS
#define CREATE_PNG

#include <cstdio> 
#include <cstdlib>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#if defined(_MSC_VER) || defined(__INTEL_COMPILER)
#pragma warning(disable : 981)
#pragma warning(disable : 383)
#endif


#define CHARGE				1.602176565e-19
#define MP_KG				1.6726231e-27
#define MP_MEV				938.272013
#define SPEED_OF_LIGHT		2.99792458e8
#define FOC					0
#define DEFOC				2



#ifdef __linux
#include <fenv.h>
#endif


using namespace std;



#ifdef DEBUG
void						prod(double *,vector< vector <double> > ,double );
vector< vector <double> >	defocusing (vector< vector <double> > , double , double , FILE * ,int );
vector< vector <double> >	focusing (vector< vector <double> > , double , double , FILE * ,int );
vector< vector <double> >	drift (vector< vector <double> > , double , FILE * ,int );
#else
void						prod(double *,vector< vector <double> > );
vector< vector <double> >	defocusing (vector< vector <double> > , double , double );
vector< vector <double> >	focusing (vector< vector <double> > , double , double );
vector< vector <double> >	drift (vector< vector <double> > , double );
#endif

int							dsMap(double ,double, int);
double *					optics(vector< vector <double> > ,int ,bool);
vector< vector <double> >	prodo(vector< vector <double> > , vector< vector <double> > , int );
void						scrivimatr2D (vector< vector <double> > , FILE * );
void						scrividati (double , double *, double *, FILE * );
void						scrividati_ellissi (double ,double *, double *, FILE * );
void						inizializza3D(double ***,int , int );
void						inizializza2D(double **, int );
void						scrivi_pos_part(FILE * ,double *,double );
vector< vector <double> >	simil(vector< vector <double> > ,vector< vector <double> > , vector< vector <double> >);
double *					assi_ellissi(double *, double);
void						create_gnuplot_file(string , string , double *, int , double , double , double , string *,bool);
double *					optics_T (double * , int , vector< vector <double> > );
void						massimo(double * , double * ,double *);
void						confronto (double *,double *,double ,double,double,FILE * ,bool *);



