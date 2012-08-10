
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
	printf("S = %f\nx = %f\np_x = %f\ny = %f\np_y = %f\n",S,A[0],A[1],A[2],A[3]);
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

vector< vector <double> > prodo(vector< vector <double> > A, vector< vector <double> > B, int dimensione) 
	{
	vector <vector <double> > H(4,vector<double>(4,0));

	printf("\n");
	printf("prima ciclo prodo\n");
		for(int k=0; k < 4; k++)
		{
			printf("\n");
			for(int j=0; j < 4; j++) printf("%+12.8f  ", A[k][j]);
		}
	printf("\n");
		for(int k=0; k < 4; k++)
		{
			printf("\n");
			for(int j=0; j < 4; j++) printf("%+12.8f  ", B[k][j]);
		}
	printf("\n");

	for (int i=0; i<dimensione; i++)
		for (int j=0; j<dimensione; j++) 
        	for (int k=0; k<dimensione; k++)             
	        	        H[i][j] += A[i][k] * B[k][j];

	for (int a=0;a<dimensione;a++)
		for (int i=0;i<dimensione;i++)
			A[a][i]=H[a][i];
			
	printf("\n");
	printf("dopo ciclo prodo\n");
		for(int k=0; k < 4; k++)
		{
			printf("\n");
			for(int j=0; j < 4; j++) printf("%+12.8f  ", A[k][j]);
		}
	printf("\n");

	return A;
	}

vector< vector <double> >  defocusing (vector< vector <double> > M, double d1, double S, FILE * file,int contatore)
		{
			M[0][0]=M[1][1]=cosh(d1*S);
			M[0][1]=sinh(d1*S)/d1;
			M[1][0]=sinh(d1*S)*d1;
			M[2][2]=M[3][3]=cos(d1*S);
			M[2][3]=sin(d1*S)/d1;
			M[3][2]=-sin(d1*S)*d1;

#ifdef DEBUG
			fprintf(file,"\nMATRICE DEFOC. dentro defocusing N # %d", contatore);
			for(int k=0; k < 4; k++)
			{
				fprintf(file,"\n");
				for(int j=0; j < 4; j++) fprintf(file,"%+12.8f  ", M[k][j]);
			}
#endif
			return M;
		}

vector< vector <double> >  focusing (vector< vector <double> > M, double d1, double S, FILE * file,int contatore)
		{
			M[0][0]=M[1][1]=cos(d1*S);
			M[0][1]=sin(d1*S)/d1;
			M[1][0]=-sin(d1*S)*d1;
			M[2][2]=M[3][3]=cosh(d1*S);
			M[2][3]=sinh(d1*S)/d1;
			M[3][2]=sinh(d1*S)*d1;	

#ifdef DEBUG
			fprintf(file,"\nMATRICE FOC. dentro focusing N # %d", contatore);		
			for(int k=0; k < 4; k++)
			{
				fprintf(file,"\n");
				for(int j=0; j < 4; j++) fprintf(file,"%+12.8f  ", M[k][j]);
			}
#endif
			return M;
		}

vector< vector <double> >  drift (vector< vector <double> > M, double S, FILE * file,int contatore)
		{
			M[0][0]=M[1][1]=M[2][2]=M[3][3]=1.;
			M[0][1]=M[2][3]=S;

#ifdef DEBUG
			fprintf(file,"\nMATRICE DRIFT dentro drift N # %d", contatore);
			for(int k=0; k < 4; k++)
			{
				fprintf(file,"\n");
				for(int j=0; j < 4; j++) fprintf(file,"%+12.8f  ", M[k][j]);
			}
#endif
			return M;
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

double * pos_part(vector< vector <double> > F,FILE * posizionePart,double *vett_i,double S)
{

	vett_i=prod(vett_i,F,S);
	fprintf(posizionePart," %+10.5f",S);
	fprintf(posizionePart," %+10.5f",vett_i[0]);
	fprintf(posizionePart," %+10.5f",vett_i[1]);
	fprintf(posizionePart," %+10.5f",vett_i[2]);
	fprintf(posizionePart," %+10.5f\n",vett_i[3]);
	return vett_i;

}

void assi_ellissi(vector< vector <double> > F,vector< vector <double> > OI, vector< vector <double> > O, double S, FILE * output2)
{

	double * A=new double [2];
	double * B=new double [2];
	double AMax,AMin,BMax,BMin;
	vector <vector <double> > K(4,vector<double>(4,0.0));

	for (int i=0;i<4;i++)
		for (int j=0;j<4;j++)
			for (int a=0;a<4;a++)
				for (int l=0;l<4;l++) K[i][j]+=O[i][a]*F[a][l]*OI[l][j];

	for (int i=0;i<4;i++)
		for (int j=0;j<4;j++) 
			F[i][j]=K[i][j];

	A=optics (F,FOC);
	B=optics (F,DEFOC);

	AMax= sqrt(A[1]*EPSILON);
	AMin= sqrt((A[0]*A[0]+1)*EPSILON/A[1]);
	BMax= sqrt(B[1]*EPSILON);
	BMin= sqrt((B[0]*B[0]+1)*EPSILON/B[1]);
	
	scrividati(S,A,B,AMax,AMin,BMax,BMin,output2);

}


int main() 
{ 
#ifdef __linux
	feenableexcept(2);
	feenableexcept(3);
#endif

	FILE * funzioni_ottiche=fopen("Funzioni_Ottiche.txt","w");
	FILE * matrici_iniziali=fopen("Matrici_Iniziali.txt","w");
	FILE * posizionePart=fopen("Posizione_Particelle.txt","w");

#ifdef DEBUG
	FILE * outputDEBUG=fopen("DEBUG.txt","w");
#endif

	string utile_per_contare;
	int conta_righe_parametri = 0;
	bool fallita_lettura=true;
	ifstream parametri;
	parametri.open("parametri.txt");
	fallita_lettura=parametri.fail();
	if (fallita_lettura) {printf("Impossibile aprire parametri.txt\n");exit(204);}
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


	vector <vector <double> > I(4,vector<double>(4,0.0));
	vector <vector <double> > K(4,vector<double>(4,0.0));
	vector <vector <double> > F(4,vector<double>(4,0.0));
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
		fprintf(outputDEBUG,"\ngrad. foc.    %+20.10f ", f1[i]*f1[i]);
		fprintf(outputDEBUG,"\ngrad. defoc.  %+20.10f ", d1[i]*d1[i]);
		}

	for (int i=0;i<contatore;i++)
		{
		if (elemento[i]=='F')
		Fx[i]=focusing(Fx[i],f1[i],lunghezza[i],matrici_iniziali,i);
		else if (elemento[i]=='D')
		Dx[i]=defocusing(Dx[i],d1[i],lunghezza[i],matrici_iniziali,i);
		else if (elemento[i]=='O')
		O[i]=drift(O[i],lunghezza[i],matrici_iniziali,i);
		else
		fprintf(outputDEBUG,"Elemento[%d]= %c non riconosciuto\n", i,elemento[i]);
		}

#ifdef DEBUG

	for (int i=0;i<contatore;i++)
	{
		if (elemento[i]=='O')
		{
//			if (O[i][0][0] == 0.0) continue;
			fprintf(matrici_iniziali,"\nMATRICE DRIFT");
			scrivimatr2D(O[i],matrici_iniziali);
			fprintf(matrici_iniziali,"\n");
		}
		else if (elemento[i]=='F')
		{
//			if (Fx[i][0][0] == 0.0) continue;
			fprintf(matrici_iniziali,"\nMATRICE FOC.");
			scrivimatr2D(Fx[i],matrici_iniziali);
			fprintf(matrici_iniziali,"\n");
		}
		else if (elemento[i]=='D')
		{
//			if (Dx[i][0][0] == 0.0) continue;
			fprintf(matrici_iniziali,"\nMATRICE DEFOC.");
			scrivimatr2D(Dx[i],matrici_iniziali);
			fprintf(matrici_iniziali,"\n");
		}

	}

#endif	

/************************************************************************/


//	FODO composizione matrici

	vector <vector <double> > compos(4,vector<double>(4,0));
		if (elemento[0]=='O')
		{
			for (int k=0; k<4; k++)
				for (int j=0; j<4; j++)
					compos[k][j]=O[0][k][j];

		}
		else if (elemento[0]=='F')
		{
			for (int k=0; k<4; k++)
				for (int j=0; j<4; j++)
					compos[k][j]=Fx[0][k][j];
		}
		else if (elemento[0]=='D')
		{
			for (int k=0; k<4; k++)
				for (int j=0; j<4; j++)
					compos[k][j]=Dx[0][k][j];
		}

	for (int i=1;i<contatore;i++)
	{
		if (elemento[i]=='O')
		{
			compos=prodo(O[i],compos,4);
		}
		else if (elemento[i]=='F')
		{
			compos=prodo(Fx[i],compos,4);
		}
		else if (elemento[i]=='D')
		{
			compos=prodo(Dx[i],compos,4);
		}

	}
	
for (int i=0;i<4;i++)
	for(int a=0;a<4;a++)
		F[i][a]=compos[i][a];

	/************************************************************************/

	// Calcolo Funzioni OTTICHE

		double *A=new double[2], *B=new double[2], *vett_i=new double[4], AMax, AMin, BMax, BMin;

	vett_i[0]=pos_x;
	vett_i[1]=imp_x;
	vett_i[2]=pos_y;
	vett_i[3]=imp_y;

	fprintf(posizionePart," %+10.5f",S=0.0);
	fprintf(posizionePart," %+10.5f",vett_i[0]);
	fprintf(posizionePart," %+10.5f",vett_i[1]);
	fprintf(posizionePart," %+10.5f",vett_i[2]);
	fprintf(posizionePart," %+10.5f\n",vett_i[3]);
	
	fprintf(outputDEBUG, "\nFODO:");
	scrivimatr2D(F,outputDEBUG);
	A = optics(F,FOC);
	B = optics(F,DEFOC);

/************************************************************************/

//ora primi dell'iterazione mi calcolo le micromappe Li di lunghezza S=L/n


//Calcolo MICROMAPPE per il Drift
	for (int i=0; i < contatore; i++)
	{
		S=lunghezza[i]/dsMap(lunghezza[i]);
		if (elemento[i] == 'O')
		{
			O[i]=drift(O[i],S,matrici_iniziali,i);
			OI[i]=drift(OI[i],-S,matrici_iniziali,i);
		}
		// per il Focusing
		else if (elemento[i] == 'F')
		{
			Fx[i]=focusing(Fx[i],f1[i],S,matrici_iniziali,i);
			FxI[i]=focusing(FxI[i],f1[i],-S,matrici_iniziali,i);
		}
		//per il Defocusing
		else if (elemento[i] == 'D')
		{
			Dx[i]=defocusing(Dx[i],d1[i],S,matrici_iniziali,i);
			DxI[i]=defocusing(DxI[i],d1[i],-S,matrici_iniziali,i);
		}
	}


	fprintf(funzioni_ottiche,"\n#   S  ");
	fprintf(funzioni_ottiche,"      Alpha x  ");
	fprintf(funzioni_ottiche,"   Beta x ");
	fprintf(funzioni_ottiche,"  Alpha y  ");
	fprintf(funzioni_ottiche,"  Beta y ");
	
	AMax= sqrt(A[1]*EPSILON);
	AMin= sqrt((A[0]*A[0]+1)*EPSILON/A[1]);
	BMax= sqrt(B[1]*EPSILON);
	BMin= sqrt((B[0]*B[0]+1)*EPSILON/B[1]);

	scrividati(0.0,A,B,AMax,AMin,BMax,BMin,funzioni_ottiche);

	double dl=0.0;

	for (int i=0;i<contatore;i++)
	{
		if (elemento[i]=='O')
		{
			fprintf(matrici_iniziali,"\n#Drift #%d",i);
			fprintf(funzioni_ottiche,"\n#Drift #%d",i);
			while (S<lunghezza[i])
			{
				dl=lunghezza[i]/dsMap(lunghezza[i]);
				fprintf(matrici_iniziali,"\n\n Num_Step %f", S);
				scrivimatr2D(F,matrici_iniziali);
				
				vett_i=pos_part(F,posizionePart,vett_i,S);
				
				assi_ellissi(F,OI[i],O[i],S,funzioni_ottiche);
				S+=dl;
			}
		}
		else if (elemento[i]=='F')
		{
			fprintf(matrici_iniziali,"\n#Foc. #%d",i);			
			fprintf(funzioni_ottiche,"\n#Foc. #%d",i);
			while (S<lunghezza[i])
			{
				dl=lunghezza[i]/dsMap(lunghezza[i]);
				fprintf(matrici_iniziali,"\n\n Num_Step %f", S);
				scrivimatr2D(F,matrici_iniziali);
				
				vett_i=pos_part(F,posizionePart,vett_i,S);
				
				assi_ellissi(F,FxI[i],Fx[i],S,funzioni_ottiche);
				S+=dl;

			}
		}
		else if (elemento[i]=='D')
		{
			fprintf(matrici_iniziali,"\n#Defoc. #%d",i);
			fprintf(funzioni_ottiche,"\n#Defoc. #%d",i);
			while (S<lunghezza[i])
			{
				dl=lunghezza[i]/dsMap(lunghezza[i]);
				fprintf(matrici_iniziali,"\n\n Num_Step %f", S);
				scrivimatr2D(F,matrici_iniziali);
				
				vett_i=pos_part(F,posizionePart,vett_i,S);
				
				assi_ellissi(F,FxI[i],Fx[i],S,funzioni_ottiche);
				S+=dl;
			}
		}
	}
		
	fclose(funzioni_ottiche);
	fclose(matrici_iniziali);
	fclose(posizionePart);
	fclose(outputDEBUG);
	parametri.close();
}
