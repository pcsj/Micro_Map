
#include "FODO_sing_part.h"

int dsMap(double L)
{
	double ds = 1.0/100;
	int n_p;
	n_p= (int) floor(L/ds);
	return n_p;
}


double *prod(double * A,vector< vector <double> > N,double S)
{
	double *H = new double[4];
	for (int i=0; i<4; i++) H[i] = 0.0;

	for (int i=0; i<4; i++)
		for (int j=0; j<4; j++)
			H[i] += A[j] * N[i][j];
	for (int i=0; i<4; i++)
		A[i]=H[i];
	printf("\nA[0] = %f\nA[1] = %f\nA[2] = %f\nA[3] = %f\nS=%f",A[0],A[1],A[2],A[3],S);
	return A;
}



double *optics(vector< vector <double> > N,int i)
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

void prodo(vector< vector <double> > A, vector< vector <double> > B, int dimensione) 
	{
	vector <vector <double> > H(4,vector<double>(4,0));

	for (int i=0; i<dimensione; i++)
		for (int j=0; j<dimensione; j++) 
        	for (int k=0; k<dimensione; k++)             
	        	        H[i][j] += A[i][k] * B[k][j];
			for (int a=0;a<dimensione;a++)
				for (int i=0;i<dimensione;i++)
					A[a][i]=H[a][i];
	}

		
void  defocusing (vector< vector <double> > M, double d1, double S, FILE * file)
		{
			M[0][0]=M[1][1]=cosh(d1*S);
			M[0][1]=sinh(d1*S)/d1;
			M[1][0]=sinh(d1*S)*d1;
			M[2][2]=M[3][3]=cos(d1*S);
			M[2][3]=sin(d1*S)/d1;
			M[3][2]=-sin(d1*S)*d1;

			for(int k=0; k < 4; k++)
			{
				fprintf(file,"\n");
				for(int j=0; j < 4; j++) fprintf(file,"%+12.8f  ", M[k][j]);
			}
		}

void  focusing (vector< vector <double> > M, double d1, double S, FILE * file)
		{

			M[0][0]=M[1][1]=cos(d1*S);
			M[0][1]=sin(d1*S)/d1;
			M[1][0]=-sin(d1*S)*d1;
			M[2][2]=M[3][3]=cosh(d1*S);
			M[2][3]=sinh(d1*S)/d1;
			M[3][2]=sinh(d1*S)*d1;	
			
			for(int k=0; k < 4; k++)
			{
				fprintf(file,"\n");
				for(int j=0; j < 4; j++) fprintf(file,"%+12.8f  ", M[k][j]);
			}
		}

void  drift (vector< vector <double> > M, double S, FILE * file)
		{
			M[0][0]=M[1][1]=M[2][2]=M[3][3]=1.;
			M[0][1]=M[2][3]=S;

			for(int k=0; k < 4; k++)
			{
				fprintf(file,"\n");
				for(int j=0; j < 4; j++) fprintf(file,"%+12.8f  ", M[k][j]);
			}
		}


void scrivimatr2D (vector< vector <double> > M, FILE * output2)
	{
	for(int k=0; k < 4; k++)
		{
		fprintf(output2,"\n");
		for(int j=0; j < 4; j++) fprintf(output2,"%+12.8f  ", M[k][j]);
		}
	}


void scrividati (double s, double A[], double B[], double AMax, double AMin, double BMax, double BMin, FILE * file)
	{
		fprintf(file,"\n %+7.4f",s);
		fprintf(file," %+10.5f",A[1]);
		fprintf(file," %+10.5f",B[1]);
		fprintf(file," %+10.5f",A[0]);
		fprintf(file," %+10.5f",B[0]);
		fprintf(file," %+10.5f",AMax);
		fprintf(file," %+10.5f",AMin);
		fprintf(file," %+10.5f",BMax);
		fprintf(file," %+10.5f",BMin);
	}

void inizializza3D(double ***M,int a, int dimensione)
	{
	M= new double**[a];
	for (int kk=0; kk < dimensione; kk++)
	{
		M[kk] = new double*[dimensione];
		for (int i=0; i < dimensione; i++) M[kk][i]= new double[dimensione];
	}
	
	for (int k = 0; k < a; k++)
		for (int i = 0; i < dimensione; i++)
			for (int j = 0; j<dimensione; j++)
				M[k][i][j]=0.0;
	}

void inizializza2D(double **M, int dimensione)
	{
	M= new double*[dimensione];
	for (int kk=0; kk < dimensione; kk++)
		M[kk] = new double[dimensione];
	
		for (int i = 0; i < dimensione; i++)
			for (int j = 0; j<dimensione; j++)
				M[i][j]=0.0;
	}


void similitudine_optics(vector< vector <double> > K, vector< vector <double> > OI, vector< vector <double> > O, vector< vector <double> > F, double * vett_i,double L,double S, double lunghezza_elemento, FILE *output2,FILE *posizionePart)
{
	while (S<lunghezza_elemento)
	{
	double * A=new double [2];
	double * B=new double [2];
	double AMax,AMin,BMax,BMin;

	double dl=L/dsMap(L);

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
		fprintf(output2,"\n Num_Step %f\n", S);
		scrivimatr2D(F,output2);

	A=optics (F,FOC);
	B=optics (F,DEFOC);

	vett_i=prod(vett_i,O,S);
	fprintf(posizionePart," %+10.5f",S);
	fprintf(posizionePart," %+10.5f",vett_i[0]);
	fprintf(posizionePart," %+10.5f",vett_i[1]);
	fprintf(posizionePart," %+10.5f",vett_i[2]);
	fprintf(posizionePart," %+10.5f\n",vett_i[3]);

		//Calcolo assi ellisse		

		AMax= sqrt(A[1]*EPSILON);
		AMin= sqrt((A[0]*A[0]+1)*EPSILON/A[1]);
		BMax= sqrt(B[1]*EPSILON);
		BMin= sqrt((B[0]*B[0]+1)*EPSILON/B[1]);
	
		//--------
	
		S+=dl;
		scrividati(S,A,B,AMax,AMin,BMax,BMin,output2);
	}
}



int main() 
{ 
#ifdef __linux
	feenableexcept(2);
	feenableexcept(3);
#endif

	FILE * output=fopen("risultatiFINALI.txt","w");
	FILE * output2=fopen("Matrici_Prova.txt","w");
	FILE * output3=fopen("Matrice.txt","w");

	string utile_per_contare;
	int conta_righe_parametri = 0;
	ifstream parametri("parametri.txt");
//	FILE * parametri=fopen("parametri.txt","r");
	do
	{
		parametri >> utile_per_contare;
		if(parametri.eof()) break;
		parametri.ignore(1000, '\n');
		conta_righe_parametri++;
	}
	while(!parametri.eof());
	parametri.clear();
	parametri.seekg(0,std::ios::beg);

	// qui di leggono tutti i dati
	char *elemento=new char[conta_righe_parametri];
	char *buffer=new char[100];
	double * lunghezza= new double[conta_righe_parametri];
	double * gradiente= new double[conta_righe_parametri];
	int contatore=0;
	for (int i = 0; i < conta_righe_parametri; i++)
	{
		parametri >> elemento[i];
		parametri >> gradiente[i];
		parametri >> lunghezza[i];
		printf("tipo elemento: %c, gradiente: %f, lunghezza: %f\n", elemento[contatore], gradiente[contatore], lunghezza[contatore]);
		contatore++;
	}

	printf("contatore: %d",contatore);
	//Qui inizializziamo le matrici, attenzione: Dx, O, Fx sono dei vettori di matrici e così anche OI,DxI,FxI

	//double ** F = new double*[4];
	//for (int kk=0; kk < 4; kk++) F[kk] = new double[4];
	//double *** Dx, *** Fx, *** O, *** DxI, *** OI, *** FxI;
	//double ** I, ** K; 
	//inizializza3D(Dx,contatore,4);
	//inizializza3D(Fx,contatore,4);
	//inizializza3D(O,contatore,4);
	//inizializza3D(DxI,contatore,4);
	//inizializza3D(OI,contatore,4);
	//inizializza3D(FxI,contatore,4);
	//inizializza2D(K,4);
	//inizializza2D(I,4);


	vector <vector <double> > I(4,vector<double>(4,0));
	vector <vector <double> > K(4,vector<double>(4,0));
	vector <vector <double> > F(4,vector<double>(4,0));
	vector <vector <vector <double> > > Fx(contatore,vector <vector <double> > (4, vector <double> (4,0.0)));
	vector <vector <vector <double> > > Dx(contatore,vector <vector <double> > (4, vector <double> (4,0.0)));
	vector <vector <vector <double> > > OI(contatore,vector <vector <double> > (4, vector <double> (4,0.0)));
	vector <vector <vector <double> > > FxI(contatore,vector <vector <double> > (4, vector <double> (4,0.0)));
	vector <vector <vector <double> > > DxI(contatore,vector <vector <double> > (4, vector <double> (4,0.0)));
	vector <vector <vector <double> > > O(contatore,vector <vector <double> > (4, vector <double> (4,0.0)));

	double gamma_beta=sqrt(2.0*ENERGIA/MP_MEV);
	double gamma_v=gamma_beta*SPEED_OF_LIGHT;
	double *f1 =new double [contatore];
	double *d1 =new double [contatore];
	for (int i=0;i<contatore;i++)
		{
		if (elemento[i]=='F')
		f1[i]=sqrt(gradiente[i]*CHARGE/(MP_KG*gamma_v));
		if (elemento[i]=='D')
		d1[i]=sqrt(gradiente[i]*CHARGE/(MP_KG*gamma_v));
		}

	double S=1.0;
	S = S/n_step;

	for (int i=0;i<contatore;i++)
		{
		fprintf(output,"\ngrad. foc.    %+20.10f ", f1[i]*f1[i]);
		fprintf(output,"\ngrad. defoc.  %+20.10f ", d1[i]*d1[i]);
		}

	for (int i=0;i<contatore;i++)
		{
		if (elemento[i]=='F')
		focusing(Fx[i],f1[i],CELL_LENGTH,output3);
		else if (elemento[i]=='D')
		defocusing(Dx[i],d1[i],CELL_LENGTH,output3);
		else if (elemento[i]=='O')
		drift(O[i],DRIFT_LENGTH,output3);
		else
		printf("Qua si e' sbagliato qualcosa\n");
		}

#ifdef DEBUG

	for (int i=0;i<contatore;i++)
	{
		if (elemento[i]=='O')
		{
//			if (O[i][0][0] == 0.0) continue;
			fprintf(output2,"\nMATRICE DRIFT");
			scrivimatr2D(O[i],output2);
			fprintf(output2,"\n");
		}
		else if (elemento[i]=='F')
		{
//			if (Fx[i][0][0] == 0.0) continue;
			fprintf(output2,"\nMATRICE FOC.");
			scrivimatr2D(Fx[i],output2);
			fprintf(output2,"\n");
		}
		else if (elemento[i]=='D')
		{
//			if (Dx[i][0][0] == 0.0) continue;
			fprintf(output2,"\nMATRICE DEFOC.");
			scrivimatr2D(Dx[i],output2);
			fprintf(output2,"\n");
		}

	}

#endif	

/************************************************************************/


//	FODO composizione matrici

	vector <vector <double> > compos(4,vector<double>(4,0));
	compos[0][0]=compos[1][1]=compos[2][2]=compos[3][3]=1.0;
	
	for (int i=0;i<contatore;i++)
	{
		if (elemento[i]=='O')
		{
			prodo(compos,O[i],4);
		}
		else if (elemento[i]=='F')
		{
			prodo(compos,Fx[i],4);
		}
		else if (elemento[i]=='D')
		{
			prodo(compos,Dx[i],4);
		}

	}
	
/************************************************************************/

	// Calcolo Funzioni OTTICHE

		double *A=new double[2], *B=new double[2], *vett_i=new double[4], AMax, AMin, BMax, BMin;
			
	vett_i[0]=pos_x;
	vett_i[1]=imp_x;
	vett_i[2]=pos_y;
	vett_i[3]=imp_y;

	FILE * posizionePart=fopen("Posizione_Particelle.txt","w");

	fprintf(posizionePart," %+10.5f",S=0.0);
	fprintf(posizionePart," %+10.5f",vett_i[0]);
	fprintf(posizionePart," %+10.5f",vett_i[1]);
	fprintf(posizionePart," %+10.5f",vett_i[2]);
	fprintf(posizionePart," %+10.5f\n",vett_i[3]);

	A = optics(F,FOC);
	B = optics(F,DEFOC);

/************************************************************************/

//ora primi dell'iterazione mi calcolo le micromappe Li di lunghezza S=L/n

	S=DRIFT_LENGTH/dsMap(DRIFT_LENGTH);
//Calcolo MICROMAPPE per il Drift
	for (int i=0; i < contatore; i++)
	{
		if (elemento[i] == 'O')
		{
			drift(O[i],S,output2);
			drift(OI[i],-S,output2);
		}
	}

	S=CELL_LENGTH/dsMap(CELL_LENGTH);
// per il Focusing
	for (int i=0; i < contatore; i++)
	{
		if (elemento[i] == 'F')
		{
			focusing(Fx[i],f1[i],S,output2);
			focusing(FxI[i],f1[i],-S,output2);
		}
	}

	//per il Defocusing
	for (int i=0; i < contatore; i++)
	{
		if (elemento[i] == 'D')
		{
			defocusing(Dx[i],d1[i],S,output2);
			defocusing(DxI[i],d1[i],-S,output2);
		}
	}


	fprintf(output,"\n#   S  ");
	fprintf(output,"      Alpha x  ");
	fprintf(output,"   Beta x ");
	fprintf(output,"  Alpha y  ");
	fprintf(output,"  Beta y ");
	
	AMax= sqrt(A[1]*EPSILON);
	AMin= sqrt((A[0]*A[0]+1)*EPSILON/A[1]);
	BMax= sqrt(B[1]*EPSILON);
	BMin= sqrt((B[0]*B[0]+1)*EPSILON/B[1]);
	
	scrividati(0.0,A,B,AMax,AMin,BMax,BMin,output);

	for (int i=0;i<contatore;i++)
	{
		if (elemento[i]=='O')
		{
			similitudine_optics(K, OI[i], O[i], F, vett_i, DRIFT_LENGTH, S, DRIFT_LENGTH, output2, posizionePart);
		}
		else if (elemento[i]=='F')
		{
			similitudine_optics(K, FxI[i], Fx[i], F, vett_i, CELL_LENGTH, S, FODO_LENGTH, output2, posizionePart);
		}
		else if (elemento[i]=='D')
		{
			similitudine_optics(K, DxI[i], Dx[i], F, vett_i, CELL_LENGTH, S, CELL_LENGTH+DRIFT_LENGTH, output2, posizionePart);
		}

	}
		
	fclose(output);
	fclose(output2);
	fclose(posizionePart);
	fclose(output3);
	parametri.close();
}
