#define _CRT_SECURE_NO_WARNINGS

#define DEBUG
//#define TEST_OPTICAL_FUNCTIONS

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


#define EPSILON				1.0
#define ENERGIA				30.0
#define CHARGE				1.602176565e-19
#define MP_KG				1.6726231e-27
#define MP_MEV				938.272013
#define SPEED_OF_LIGHT		2.99792458e8
#define FOC					0
#define DEFOC				2
#define N_STEP				100




#define pos_x 0.0
#define pos_y 0.0
#define imp_x 0.05
#define imp_y 0.05

#ifdef __linux
#include <fenv.h>
#endif


using namespace std;



int dsMap(double ,double, int);

void prod(double *,vector< vector <double> > ,double );

double *optics(vector< vector <double> > ,int ,bool);

vector< vector <double> > prodo(vector< vector <double> > , vector< vector <double> > , int );

vector< vector <double> >  defocusing (vector< vector <double> > , double , double , FILE * ,int );

vector< vector <double> >  focusing (vector< vector <double> > , double , double , FILE * ,int );

vector< vector <double> >  drift (vector< vector <double> > , double , FILE * ,int );

void scrivimatr2D (vector< vector <double> > , FILE * );

void scrividati (double , double * , double *, double *, double *, FILE * );

void inizializza3D(double ***,int , int );

void inizializza2D(double **, int );

void scrivi_pos_part(FILE * ,double *,double );

vector< vector <double> >  simil(vector< vector <double> > ,vector< vector <double> > , vector< vector <double> >);

double * assi_ellissi(double *);

void create_gnuplot_file(string , string , double *, int , double , double , double );

double * optics_T (double * , int , vector< vector <double> > );


