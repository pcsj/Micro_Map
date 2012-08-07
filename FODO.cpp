#include <cstdio> 
#include <cstdlib>
#include <cmath>
#include <iostream>

#ifdef __linux
#include <fenv.h>
#endif

#define DEBUG
#define MOD1
//#define MOD2

using namespace std;



int dsMap(double L)
{
	double ds = 1.0/100;
	int n_p;
	n_p= (int) floor(L/ds);
	return n_p;
}


double *optics(double ** N,int i)
{
	double  oem,s,*A=new double[2];
	oem=acos((N[i][i]+N[i+1][i+1])*0.5);
	s=sin(oem);
	
	if ( s*(N[i][i+1]) < 0.) s=-s;
	
	A[0]= (N[i][i] - cos(oem)) / s;
	A[1]= N[i][i+1] / s;

#ifdef DEBUG
	printf("\nN[%i][%i] = %f\nN[%i][%i] = %f",i,i,N[i][i],i+1,i+1,N[i+1][i+1]);
	printf("\nomega = %f",oem);
	printf("\nsin(omega) = %f",s);
	printf("\nA[0] = %f\nA[1] = %f",A[0],A[1]);
	printf("\n");
#endif	
	
	return A;
}

/*
double **prod(double ** A, double ** B) 
	{
	double ** H = new double*[4];
	for (int kk=0; kk < 4; kk++) H[kk] = new double[4];
	for (int i=0; i<4; i++)
		for (int j=0; j<4; j++) 
			{
	        	H[i][j]=0;
	        	for (int k=0; k<4; k++)             
	        	        H[i][j] += A[i][k] * B[k][j];
			for (int a=0;a<4;a++)
				for (int i=0;i<4;i++)
					A[a][i]=H[a][i];
		        }
	return A;
	}
*/

void scrivimatr (double ** M, FILE * output2)
	{
	for(int k=0; k < 4; k++)
		{
		fprintf(output2,"\n");
		for(int j=0; j < 4; j++) fprintf(output2,"%+12.8f  ", M[k][j]);
		}
	}



int main() 
{ 
#ifdef __linux
	feenableexcept(2);
	feenableexcept(3);
#endif
	double ** F = new double*[4];
	for (int kk=0; kk < 4; kk++) F[kk] = new double[4];
	double ** Dx = new double*[4];
	for (int kk=0; kk < 4; kk++) Dx[kk] = new double[4];
	double ** Fx = new double*[4];
	for (int kk=0; kk < 4; kk++) Fx[kk] = new double[4];
	double ** I = new double*[4];
	for (int kk=0; kk < 4; kk++) I[kk] = new double[4];
	double ** K = new double*[4];
	for (int kk=0; kk < 4; kk++) K[kk] = new double[4];
	double ** O = new double*[4];
	for (int kk=0; kk < 4; kk++) O[kk] = new double[4];
	double ** OI = new double*[4];
	for (int kk=0; kk < 4; kk++) OI[kk] = new double[4];
	double ** FxI = new double*[4];
	for (int kk=0; kk < 4; kk++) FxI[kk] = new double[4];
	double ** DxI = new double*[4];
	for (int kk=0; kk < 4; kk++) DxI[kk] = new double[4];
	for (int a = 0; a < 4; a++)
		for (int i = 0; i < 4; i++)
			K[a][i]=I[a][i]=F[a][i]=O[a][i]=Fx[a][i]=Dx[a][i]=OI[a][i]=FxI[a][i]=DxI[a][i]=0.0;


	FILE * output=fopen("risultatiFINALI.txt","w");
//	FILE * output1=fopen("matrici.txt","w");
	FILE * output2=fopen("Matrici_Prova.txt","w");
/************	
#define EPSILON				1.0
#define DRIFT_LENGTH		0.3
#define CELL_LENGTH			0.2
#define FODO_LENGTH			1.0
#define ENERGIA				30.0
#define GRADIENTE_FOC1		12.15
#define GRADIENTE_DEF1		11.85
#define FOC					0
#define DEFOC				2
************/
	double EPSILON=1.0;
	double DRIFT_LENGTH=0.3;
	double CELL_LENGTH=0.2;
	double FODO_LENGTH=1.0;
	double ENERGIA=30.0;
	double GRADIENTE_FOC1=12.15;
	double GRADIENTE_DEF1=11.85;
	int FOC=0;
	int DEFOC=2;
	int n_step=100;

	double f1 = sqrt(GRADIENTE_FOC1);
	double d1 = sqrt(GRADIENTE_DEF1);
	double e = CELL_LENGTH;

	double S=1.0;
	S = S/n_step;
	double dl=S;

	fprintf(output,"\ngrad. foc.    %+20.10f ", f1*f1);
	fprintf(output,"\ngrad. defoc.  %+20.10f ", d1*d1);

	O[0][0]=O[1][1]=O[2][2]=O[3][3]=1.;
	O[0][1]=O[2][3]=DRIFT_LENGTH;

	Fx[0][0]=Fx[1][1]=cos(f1*e);
	Fx[0][1]=sin(f1*e)/f1;
	Fx[1][0]=-sin(f1*e)*f1;
	Fx[2][2]=Fx[3][3]=cosh(f1*e);
	Fx[2][3]=sinh(f1*e)/f1;
	Fx[3][2]=sinh(f1*e)*f1;

	Dx[0][0]=Dx[1][1]=cosh(d1*e);
	Dx[0][1]=sinh(d1*e)/d1;
	Dx[1][0]=sinh(d1*e)*d1;
	Dx[2][2]=Dx[3][3]=cos(d1*e);
	Dx[2][3]=sin(d1*e)/d1;
	Dx[3][2]=-sin(d1*e)*d1;

#ifdef DEBUG

		fprintf(output2,"\nMATRICE DRIFT");
		scrivimatr(O,output2);
	
		fprintf(output2,"\nMATRICE FOC.");
		scrivimatr(Fx,output2);

		fprintf(output2,"\nMATRICE DEFOC.");
		scrivimatr(Dx,output2);

		fprintf(output2,"\n");

#endif	

/************************************************************************/


//	FODO composizione matrici


	for (int a = 0; a < 4; a++)
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)	
				K[a][i] += Dx[a][j] * O[j][i];
	for (int a = 0; a < 4; a++)
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)	
				I[a][i] += Fx[a][j] * O[j][i];
	for (int a = 0; a < 4; a++)
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)	
				F[a][i] += I[a][j] * K[j][i];

/*
	fprintf(output1,"\nMATRICE FODO x");
	for(int k=0; k < 4; k++)
		{
			fprintf(output1,"\n");
			for(int j=0; j < 4; j++) fprintf(output1,"%+12.8f  ", F[k][j]);
		}
	fprintf(output1,"\n");
*/
#ifdef DEBUG
		fprintf(output2,"\nMATRICE DRIFT * FOC.");
		scrivimatr(I,output2);

		fprintf(output2,"\nMATRICE DRIFT * DEFOC.");
		scrivimatr(K,output2);		

		fprintf(output2,"\nMATRICE 1 * 2");
		scrivimatr(F,output2);
#endif

/************************************************************************/

	// Calcolo Funzioni OTTICHE

	double *A=new double[2], *B=new double[2], AMax, AMin, BMax, BMin;
	
	A = optics(F,FOC);
	B = optics(F,DEFOC);

/************************************************************************/

//ora primi dell'iterazione mi calcolo le micromappe Li di lunghezza S=L/n

	int N;
	N = (int) (2.*dsMap(DRIFT_LENGTH)+2.*dsMap(CELL_LENGTH));

	S=DRIFT_LENGTH/dsMap(DRIFT_LENGTH);
	O[0][0]=O[1][1]=O[2][2]=O[3][3]=1.;
	O[0][1]=O[2][3]=S;
	
	#ifdef	MOD1
		//inversi
		OI[0][0]=OI[1][1]=OI[2][2]=OI[3][3]=1.0;
		OI[0][1]=OI[2][3]=-S;
	#endif	

	S=CELL_LENGTH/dsMap(CELL_LENGTH);	
	Dx[0][0]=Dx[1][1]=cosh(d1*S);
	Dx[0][1]=sinh(d1*S)/d1;
	Dx[1][0]=sinh(d1*S)*d1;
	Dx[2][2]=Dx[3][3]=cos(d1*S);
	Dx[2][3]=sin(d1*S)/d1;
	Dx[3][2]=-sin(d1*S)*d1;

	Fx[0][0]=Fx[1][1]=cos(f1*S);
	Fx[0][1]=sin(f1*S)/f1;
	Fx[1][0]=-sin(f1*S)*f1;
	Fx[2][2]=Fx[3][3]=cosh(f1*S);
	Fx[2][3]=sinh(f1*S)/f1;
	Fx[3][2]=sinh(f1*S)*f1;

	#ifdef	MOD1
		//inversi

		DxI[0][0]=DxI[1][1]=cosh(-d1*S);
		DxI[0][1]=sinh(-d1*S)/d1;
		DxI[1][0]=sinh(-d1*S)*d1;
		DxI[2][2]=DxI[3][3]=cos(-d1*S);
		DxI[2][3]=sin(-d1*S)/d1;
		DxI[3][2]=-sin(-d1*S)*d1;

		FxI[0][0]=FxI[1][1]=cos(-f1*S);
		FxI[0][1]=sin(-f1*S)/f1;
		FxI[1][0]=-sin(-f1*S)*f1;
		FxI[2][2]=FxI[3][3]=cosh(-f1*S);
		FxI[2][3]=sinh(-f1*S)/f1;
		FxI[3][2]=sinh(-f1*S)*f1;
	#endif	

	#ifdef	DEBUG

		fprintf(output2,"\n\nNumero Esatto di Micromappe  %d\n\n",N);
		
		fprintf(output2,"\nMICROMATRICE DRIFT");
		scrivimatr(O,output2);

		fprintf(output2,"\nMICROMATRICE DEF.");
		scrivimatr(Dx,output2);
	
		fprintf(output2,"\nMICROMATRICE FOC.");
		scrivimatr(Fx,output2);

	#endif

	fprintf(output,"\n#   S  ");
	fprintf(output,"      Alpha x  ");
	fprintf(output,"   Beta x ");
	fprintf(output,"  Alpha y  ");
	fprintf(output,"  Beta y ");

	fprintf(output,"\n %+7.4f",S=0.0);
	fprintf(output," %+10.5f",A[1]);
	fprintf(output," %+10.5f",B[1]);
	fprintf(output," %+10.5f",A[0]);
	fprintf(output," %+10.5f",B[0]);

// Primo ciclo per la O

	while (S<DRIFT_LENGTH)
	{
		dl=DRIFT_LENGTH/dsMap(DRIFT_LENGTH);

#ifdef MOD2

		A[1]= A[0] - (O[0][0]*O[1][0] * B[0])  +  (2.*O[1][0]*O[0][1] * A[0])  -  (1./B[0]) * O[0][1]*O[1][1] * (1. - A[0]*A[0]);
		B[1]= pow(O[0][0],2) * B[0]  - 2.*O[0][0]*O[0][1] * A[0]  -  (1./B[0]) * pow(O[0][1],2.) * (1. - (A[0]*A[0]));

		C[1]= C[0] - (O[2][2]*O[3][2] * D[0])  +  (2.*O[3][2]*O[2][3] * C[0])  -  (1./D[0]) * O[2][3]*O[3][3] * (1. - C[0]*C[0]);
		D[1]= pow(O[2][2],2) * D[0]  - 2.*O[2][2]*O[2][3] * C[0]  -  (1./D[0]) * pow(O[2][3],2) * (1. - (C[0]*C[0]));

		A[0]=A[1]; 
		B[0]=B[1];
		C[0]=C[1]; 
		D[0]=D[1];

#endif

#ifdef MOD1

		for (int i=0;i<4;i++)
		{
			for (int j=0;j<4;j++)
			{
				K[i][j]=0.0;
				for (int a=0;a<4;a++)
				{
					for (int l=0;l<4;l++) K[i][j]+=O[i][a]*F[a][l]*OI[l][j];
				}
			}
		}
		for (int i=0;i<4;i++)
		{
			for (int j=0;j<4;j++) F[i][j]=K[i][j];
		}
		fprintf(output2,"\n Num_Step %f\n", dl);
		scrivimatr(F,output2);

	A=optics (F,FOC);
	B=optics (F,DEFOC);


		A=optics (F,FOC);
		B=optics (F,DEFOC);
 
#endif

		//Calcolo assi ellisse		

		AMax= sqrt(A[1]*EPSILON);
		AMin= sqrt((A[0]*A[0]+1)*EPSILON/A[1]);
		BMax= sqrt(B[1]*EPSILON);
		BMin= sqrt((B[0]*B[0]+1)*EPSILON/B[1]);
	
		//--------
	
		S+=dl;
		fprintf(output,"\n %+7.4f",S);
		fprintf(output," %+10.5f",A[1]);
		fprintf(output," %+10.5f",B[1]);
		fprintf(output," %+10.5f",A[0]);
		fprintf(output," %+10.5f",B[0]);
//		fprintf(output," %+10.5f",C[0]);
//		fprintf(output," %+10.5f",D[0]);
		fprintf(output," %+10.5f",AMax);
		fprintf(output," %+10.5f",AMin);
		fprintf(output," %+10.5f",BMax);
		fprintf(output," %+10.5f",BMin);

	}

//Secondo Ciclo per la D in x e la F in y

	while (S<(DRIFT_LENGTH+CELL_LENGTH))
	{
		dl=CELL_LENGTH/dsMap(CELL_LENGTH);
	
#ifdef MOD2
	
		A[1]= A[0] - (Dx[0][0]*Dx[1][0] * B[0])  +  (2.*Dx[1][0]*Dx[0][1] * A[0])  -  (1./B[0]) * Dx[0][1]*Dx[1][1] * (1. - A[0]*A[0]);
		B[1]= pow(Dx[0][0],2) * B[0]  - 2.*Dx[0][0]*Dx[0][1] * A[0]  -  (1./B[0]) * pow(Dx[0][1],2) * (1. - (A[0]*A[0]));
		
		C[1]= C[0] - (Dx[2][2]*Dx[3][2] * D[0])  +  (2.*Dx[3][2]*Dx[2][3] * C[0])  -  (1./D[0]) * Dx[2][3]*Dx[3][3] * (1. - C[0]*C[0]);
		D[1]= pow(Dx[2][2],2) * D[0]  - 2.*Dx[2][2]*Dx[2][3] * C[0]  -  (1./D[0]) * pow(Dx[2][3],2) * (1. - (C[0]*C[0]));

		A[0]=A[1]; 
		B[0]=B[1];
		C[0]=C[1]; 
		D[0]=D[1];

#endif

#ifdef MOD1

		for (int i=0;i<4;i++)
		{
			for (int j=0;j<4;j++)
			{
				K[i][j]=0.0;
				for (int a=0;a<4;a++)
				{
					for (int l=0;l<4;l++) K[i][j]+=Dx[i][a]*F[a][l]*DxI[l][j];
				}
			}
		}
		for (int i=0;i<4;i++)
		{
			for (int j=0;j<4;j++) F[i][j]=K[i][j];
		}
		fprintf(output2,"\n Num_Step %f\n", dl);
		scrivimatr(F,output2);
		
		A=optics (F,FOC);
		B=optics (F,DEFOC);		
			
#endif	

		//Calcolo assi ellisse		

		AMax= sqrt(A[1]*EPSILON);
		AMin= sqrt((A[0]*A[0]+1)*EPSILON/A[1]);
		BMax= sqrt(B[1]*EPSILON);
		BMin= sqrt((B[0]*B[0]+1)*EPSILON/B[1]);
			
		//--------
	
		S+=dl;
		fprintf(output,"\n %+7.4f",S);
		fprintf(output," %+10.5f",A[1]);
		fprintf(output," %+10.5f",B[1]);
		fprintf(output," %+10.5f",A[0]);
		fprintf(output," %+10.5f",B[0]);
//		fprintf(output," %+10.5f",C[0]);
//		fprintf(output," %+10.5f",D[0]);
		fprintf(output," %+10.5f",AMax);
		fprintf(output," %+10.5f",AMin);
		fprintf(output," %+10.5f",BMax);
		fprintf(output," %+10.5f",BMin);

	}

//Terzo ciclo ancora per O

	while (S < (CELL_LENGTH+2.*DRIFT_LENGTH) )
	{
		dl=DRIFT_LENGTH/dsMap(DRIFT_LENGTH);

#ifdef MOD2
		A[1]= A[0] - (O[0][0]*O[1][0] * B[0])  +  (2.*O[1][0]*O[0][1] * A[0])  -  (1./B[0]) * O[0][1]*O[1][1] * (1. - A[0]*A[0]);
		B[1]= pow(O[0][0],2) * B[0]  - 2.*O[0][0]*O[0][1] * A[0]  -  (1./B[0]) * pow(O[0][1],2) * (1. - (A[0]*A[0]));

		C[1]= C[0] - (O[2][2]*O[3][2] * D[0])  +  (2.*O[3][2]*O[2][3] * C[0])  -  (1./D[0]) * O[2][3]*O[3][3] * (1. - C[0]*C[0]);
		D[1]= pow(O[2][2],2) * D[0]  - 2.*O[2][2]*O[2][3] * C[0]  -  (1./D[0]) * pow(O[2][3],2) * (1. - (C[0]*C[0]));

		A[0]=A[1]; 
		B[0]=B[1];
		C[0]=C[1]; 
		D[0]=D[1];
#endif

#ifdef MOD1
		for (int i=0;i<4;i++)
		{
			for (int j=0;j<4;j++)
			{
				K[i][j]=0.0;
				for (int a=0;a<4;a++)
				{
					for (int l=0;l<4;l++) K[i][j]+=O[i][a]*F[a][l]*OI[l][j];
				}
			}
		}
		for (int i=0;i<4;i++)
		{
			for (int j=0;j<4;j++) F[i][j]=K[i][j];
		}
		fprintf(output2,"\n Num_Step %f\n", dl);
		scrivimatr(F,output2);

		A=optics (F,FOC);
		B=optics (F,DEFOC);		

#endif		

		//Calcolo assi ellisse		

		AMax= sqrt(A[1]*EPSILON);
		AMin= sqrt((A[0]*A[0]+1)*EPSILON/A[1]);
		BMax= sqrt(B[1]*EPSILON);
		BMin= sqrt((B[0]*B[0]+1)*EPSILON/B[1]);
	
		//--------

		S+=dl;
		fprintf(output,"\n %+7.4f",S);
		fprintf(output," %+10.5f",A[1]);
		fprintf(output," %+10.5f",B[1]);
		fprintf(output," %+10.5f",A[0]);
		fprintf(output," %+10.5f",B[0]);
//		fprintf(output," %+10.5f",C[0]);
//		fprintf(output," %+10.5f",D[0]);
		fprintf(output," %+10.5f",AMax);
		fprintf(output," %+10.5f",AMin);
		fprintf(output," %+10.5f",BMax);
		fprintf(output," %+10.5f",BMin);

	}

//Quarto Ciclo per la F in x e la D in y

	while (S <= FODO_LENGTH)
	{
		dl=CELL_LENGTH/dsMap(CELL_LENGTH);

#ifdef MOD2
		A[1]= A[0] - (Fx[0][0]*Fx[1][0] * B[0])  +  (2.*Fx[1][0]*Fx[0][1] * A[0])  -  (1./B[0]) * Fx[0][1]*Fx[1][1] * (1. - A[0]*A[0]);
		B[1]= pow(Fx[0][0],2) * B[0]  - 2.*Fx[0][0]*Fx[0][1] * A[0]  -  (1./B[0]) * pow(Fx[0][1],2) * (1. - (A[0]*A[0]));

		C[1]= C[0] - (Fx[2][2]*Fx[3][2] * D[0])  +  (2.*Fx[3][2]*Fx[2][3] * C[0])  -  (1./D[0]) * Fx[2][3]*Fx[3][3] * (1. - C[0]*C[0]);
		D[1]= pow(Fx[2][2],2) * D[0]  - 2.*Fx[2][2]*Fx[2][3] * C[0]  -  (1./D[0]) * pow(Fx[2][3],2) * (1. - (C[0]*C[0]));

		A[0]=A[1]; 
		B[0]=B[1];
		C[0]=C[1]; 
		D[0]=D[1];

#endif

#ifdef MOD1
		for (int i=0;i<4;i++)
		{
			for (int j=0;j<4;j++)
			{
				K[i][j]=0.0;
				for (int a=0;a<4;a++)
				{
					for (int l=0;l<4;l++) K[i][j]+=Fx[i][a]*F[a][l]*FxI[l][j];
				}
			}
		}
		for (int i=0;i<4;i++)
		{
			for (int j=0;j<4;j++) F[i][j]=K[i][j];
		}
		fprintf(output2,"\n Num_Step %f\n", dl);
		scrivimatr(F,output2);

	A=optics (F,FOC);
	B=optics (F,DEFOC);		

	#endif

		//Calcolo assi ellisse		

		AMax= sqrt(A[1]*EPSILON);
		AMin= sqrt((A[0]*A[0]+1)*EPSILON/A[1]);
		BMax= sqrt(B[1]*EPSILON);
		BMin= sqrt((B[0]*B[0]+1)*EPSILON/B[1]);

		//--------
		
		S+=dl;
		fprintf(output,"\n %+7.4f",S);
		fprintf(output," %+10.5f",A[1]);
		fprintf(output," %+10.5f",B[1]);
		fprintf(output," %+10.5f",A[0]);
		fprintf(output," %+10.5f",B[0]);
//		fprintf(output," %+10.5f",C[0]);
//		fprintf(output," %+10.5f",D[0]);
		fprintf(output," %+10.5f",AMax);
		fprintf(output," %+10.5f",AMin);
		fprintf(output," %+10.5f",BMax);
		fprintf(output," %+10.5f",BMin);

	}

	fclose(output);
//	fclose(output1);
	fclose(output2);

}


