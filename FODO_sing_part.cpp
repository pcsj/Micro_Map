
#include "FODO_sing_part.h"


int dsMap(double L, double lunghezzatotale, int n_step)
{
	double ds = lunghezzatotale/n_step;
	int n_p = (int) floor(L/ds);
	return n_p;
}


void prod(double * A, vector< vector <double> > N, double S)
{
	double *H = new double[4];
	for (int i=0; i<4; i++) H[i] = 0.0;

	for (int i=0; i<4; i++)
		for (int j=0; j<4; j++)
			H[i] += A[j] * N[i][j];

	for (int i=0; i<4; i++)
		A[i]=H[i];

#ifdef DEBUG
	printf("S = %f\nx = %f\np_x = %f\ny = %f\np_y = %f\n",S,A[0],A[1],A[2],A[3]);
#endif
}


double *optics(vector< vector <double> > N, int i)
{
	double omega, s, *A=new double[2];
	for (int j = 0; j < 2; j++) A[j] = 0.;
	if (fabs((N[i][i]+N[i+1][i+1])*0.5) > 1.)
	{
		printf("Errore impossibile calcolare le funzioni ottiche!\n");
		return A;		// ritorna zero
	}
	omega = acos((N[i][i]+N[i+1][i+1])*0.5);
	s=sin(omega);
	
	if ( s*(N[i][i+1]) < 0.) s=-s;
	
	A[0]= (N[i][i] - cos(omega)) / s;
	A[1]= N[i][i+1] / s;

#ifdef DEBUG
	printf("\nN[%i][%i] = %f\nN[%i][%i] = %f",i,i,N[i][i],i+1,i+1,N[i+1][i+1]);
	printf("\nomega = %f",omega);
	printf("\nsin(omega) = %f",s);
	printf("\nA[0] = %f\nA[1] = %f",A[0],A[1]);
	printf("\n");
#endif	
	
	return A;
}


vector< vector <double> > prodo(vector< vector <double> > A, vector< vector <double> > B, int dimensione) 
{
	vector <vector <double> > H(4,vector<double>(4,0));

#ifdef DEBUG
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
#endif

	for (int i=0; i<dimensione; i++)
		for (int j=0; j<dimensione; j++) 
        	for (int k=0; k<dimensione; k++)             
	        	        H[i][j] += A[i][k] * B[k][j];

	for (int a=0;a<dimensione;a++)
		for (int i=0;i<dimensione;i++)
			A[a][i]=H[a][i];
	
#ifdef DEBUG
	printf("\n");
	printf("dopo ciclo prodo\n");
		for(int k=0; k < 4; k++)
		{
			printf("\n");
			for(int j=0; j < 4; j++) printf("%+12.8f  ", A[k][j]);
		}
	printf("\n");
#endif

	return A;
}


vector< vector <double> >  defocusing (vector< vector <double> > M, double d1, double S, FILE * file, int contatore)
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


vector< vector <double> >  focusing (vector< vector <double> > M, double d1, double S, FILE * file, int contatore)
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


vector< vector <double> >  drift (vector< vector <double> > M, double S, FILE * file, int contatore)
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


void scrividati (double s, double *A, double *B, double *AMaxMin, double *BMaxMin, FILE * file)
{
	fprintf(file,"\n %+7.4f",s);
	fprintf(file," %+10.5f",A[1]);
	fprintf(file," %+10.5f",B[1]);
	fprintf(file," %+10.5f",A[0]);
	fprintf(file," %+10.5f",B[0]);
	fprintf(file," %+10.5f",AMaxMin[1]);
	fprintf(file," %+10.5f",BMaxMin[1]);
	fprintf(file," %+10.5f",AMaxMin[0]);
	fprintf(file," %+10.5f",BMaxMin[0]);
}


void inizializza3D(double ***M, int a, int dimensione)
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

void scrivi_pos_part(FILE * posizionePart,double *vett_i, double S)
{
	fprintf(posizionePart," %+10.5f",S);
	fprintf(posizionePart," %+10.5f",vett_i[0]);
	fprintf(posizionePart," %+10.5f",vett_i[1]);
	fprintf(posizionePart," %+10.5f",vett_i[2]);
	fprintf(posizionePart," %+10.5f\n",vett_i[3]);
}

vector< vector <double> > simil(vector< vector <double> > F,vector< vector <double> > OI, vector< vector <double> > O)
{
	vector <vector <double> > K(4,vector<double>(4,0.0));

	for (int i=0;i<4;i++)
		for (int j=0;j<4;j++)
			for (int a=0;a<4;a++)
				for (int l=0;l<4;l++) K[i][j]+=O[i][a]*F[a][l]*OI[l][j];

	for (int i=0;i<4;i++)
		for (int j=0;j<4;j++) 
			F[i][j]=K[i][j];

	return F;
}

double * assi_ellissi(double *alpha)
{
	double * minimi_massimi = new double[2];
	for (int i = 0; i < 2; i++) minimi_massimi[i] = 0.0;
	minimi_massimi[0] = sqrt((alpha[1])*EPSILON);
	minimi_massimi[1] = sqrt(((alpha[0])* (alpha[0])+1.)*EPSILON / (alpha[1]));
	return minimi_massimi;
}

void create_gnuplot_file(string gnuplot_filename, string run_name, double *lunghezza, int contatore, double estremo, double zmin, double zmax)
{
	ofstream gnuplot_file;
	double lunghezza_percorsa=0.;
	gnuplot_file.open(gnuplot_filename.c_str());
	gnuplot_file << "#!/gnuplot" << endl;
	gnuplot_file << "FILE=\"" << run_name << ".txt\"" << endl;
	gnuplot_file << "set terminal png enhanced 15" << endl;
	gnuplot_file << "set output \"graph_" << run_name << ".png\"" << endl;
	gnuplot_file << "set xrange[" << zmin << ":" << zmax << "]" << endl;
	gnuplot_file << "#set yrange[" << -estremo << ":" << estremo << "]" << endl;
	gnuplot_file << "set title  \""<< run_name <<" \"" << endl;
	gnuplot_file << "set xlabel \"z (m)\"" << endl;
	if (run_name[0]=='F')
	gnuplot_file << "set ylabel \"Funz. Ottiche x,y\"" << endl;
	else if (run_name[0]=='P')
	gnuplot_file << "set ylabel \"x,y (m)\"" << endl;
	else exit(204);
	for (int i = 0; i < contatore; i++)
	{
		lunghezza_percorsa+=lunghezza[i];
		gnuplot_file << "set arrow from " << lunghezza_percorsa<< "," << -5 << " to "<< lunghezza_percorsa << ","<< 5 << " nohead lc rgb \"black\" lw 1" << endl;
	}
	gnuplot_file << "plot FILE u 1:2 w lines lt 1 lc rgb \"red\" lw 1 t \"x\",\\" << endl;
	gnuplot_file << "FILE u 1:4 w lines lt 1 lc rgb \"blue\" lw 1 t \"y\",\\" << endl;
	gnuplot_file << "FILE u 1:3 w lines lt 1 lc rgb \"orange\" lw 1 t \"p_x\",\\" << endl;
	gnuplot_file << "FILE u 1:5 w lines lt 1 lc rgb \"dark-green\" lw 1 t \"p_y\"" << endl;

	gnuplot_file.close();
}


#ifdef TURK
void optics_T (double * A, int i, vector< vector <double> > O)
{
	double alpha = A[0], beta = A[1];

	A[0] = alpha - ((O[i][i]) * (O[i+1][i]) * beta)  +  (2. * (O[i+1][i]) * (O[i][i+1]) * alpha)  -  (1./beta) * (O[i][i+1]) * (O[i+1][i+1]) * (1. - alpha * alpha);
	A[1]= (O[i][i]) * (O[i][i]) * beta  - 2. * (O[i][i]) * (O[i][i+1]) * alpha  -  (1./beta) * (O[i][i+1]) * (O[i][i+1]) * (1. - (alpha * alpha));
}
#endif


int main() 
{ 
#ifdef __linux
	feenableexcept(2);
	feenableexcept(3);
#endif

	FILE * funzioni_ottiche=fopen("Funzioni_Ottiche.txt","w");
	FILE * matrici_iniziali=fopen("Matrici_Iniziali.txt","w");
	FILE * posizionePart=fopen("Posizione_Particelle.txt","w");

#ifdef TURK
	FILE * funzioni_otticheturk=fopen("Funzioni_Ottiche_T.txt","w");
#endif

#ifdef DEBUG
	FILE * outputDEBUG=fopen("DEBUG.txt","w");
#endif

	string utile_per_contare;
	int conta_righe_parametri = 0;
	bool fallita_lettura=true;
	ifstream parametri;
	parametri.open("parametri.txt");
	fallita_lettura=parametri.fail();
	if (fallita_lettura)
	{
		printf("Impossibile aprire parametri.txt\n");
		exit(204);
	}

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
#ifdef DEBUG
		printf("tipo elemento: %c, gradiente: %f, lunghezza: %f\n", elemento[contatore], gradiente[contatore], lunghezza[contatore]);
#endif
		contatore++;
	}

#ifdef DEBUG
	printf("contatore: %d",contatore);
#endif

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

#ifdef DEBUG
	for (int i=0;i<contatore;i++)
	{
		fprintf(outputDEBUG,"\ngrad. foc.    %+20.10f ", f1[i]*f1[i]);
		fprintf(outputDEBUG,"\ngrad. defoc.  %+20.10f ", d1[i]*d1[i]);
	}
#endif

	for (int i=0;i<contatore;i++)
	{
		if (elemento[i]=='F')
			Fx[i]=focusing(Fx[i],f1[i],lunghezza[i],matrici_iniziali,i);
		else if (elemento[i]=='D')
			Dx[i]=defocusing(Dx[i],d1[i],lunghezza[i],matrici_iniziali,i);
		else if (elemento[i]=='O')
			O[i]=drift(O[i],lunghezza[i],matrici_iniziali,i);
#ifdef DEBUG
		else
			fprintf(outputDEBUG,"Elemento[%d]= %c non riconosciuto\n", i,elemento[i]);
#endif
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
			compos=prodo(O[i],compos,4);
		else if (elemento[i]=='F')
			compos=prodo(Fx[i],compos,4);
		else if (elemento[i]=='D')
			compos=prodo(Dx[i],compos,4);
	}
	
	for (int i=0;i<4;i++)
		for(int a=0;a<4;a++)
			F[i][a]=compos[i][a];

	/************************************************************************/

	// Calcolo Funzioni OTTICHE

	double *vett_i=new double[4];

	double *alpha = new double[2];
	double *beta = new double[2];
	double *aminmax = new double[2];
	double *bminmax = new double[2];
#ifdef TURK
	double *alphaturk = new double[2];
	double *betaturk = new double[2];
	double *aminmaxturk = new double[2];
	double *bminmaxturk = new double[2];
	for (int i = 0; i < 2; i++) alphaturk[i] = betaturk[i] = aminmaxturk[i] = bminmaxturk[i] = 0.;
#endif

	for (int i = 0; i < 2; i++) alpha[i] = beta[i] = aminmax[i] = bminmax[i] = 0.;

	fprintf(funzioni_ottiche,"\n#%7c",'S');
	fprintf(funzioni_ottiche,"%10.8s","Alpha x");
	fprintf(funzioni_ottiche,"%10.7s","Beta x");
	fprintf(funzioni_ottiche,"%12.8s","Alpha y");
	fprintf(funzioni_ottiche,"%10.7s","Beta y");
	fprintf(funzioni_ottiche,"%10s","x");
	fprintf(funzioni_ottiche,"%11s","p_x");
	fprintf(funzioni_ottiche,"%11s","y");
	fprintf(funzioni_ottiche,"%11s","p_y");

	alpha=optics(F,FOC);
	beta=optics(F,DEFOC);
	aminmax = assi_ellissi(alpha);
	bminmax = assi_ellissi(beta);
	scrividati(0.0,alpha,beta,aminmax,bminmax,funzioni_ottiche);

#ifdef TURK
	// non credo sia giusto inizializzare alphaturk e betaturk con l'altra funzione ottica
	// ma le "*turk" sono ricorsive e non permettono un bootstrap per ora...
	alphaturk=optics(F,FOC);
	betaturk=optics(F,DEFOC);
	aminmaxturk = assi_ellissi(alpha);
	bminmaxturk = assi_ellissi(beta);
	scrividati(0.0,alphaturk,betaturk,aminmaxturk,bminmaxturk,funzioni_otticheturk);
#endif

	vett_i[0]=pos_x;
	vett_i[1]=imp_x;
	vett_i[2]=pos_y;
	vett_i[3]=imp_y;

	fprintf(posizionePart," %+10.5f",0.0);
	fprintf(posizionePart," %+10.5f",vett_i[0]);
	fprintf(posizionePart," %+10.5f",vett_i[1]);
	fprintf(posizionePart," %+10.5f",vett_i[2]);
	fprintf(posizionePart," %+10.5f\n",vett_i[3]);
	
#ifdef DEBUG
	fprintf(outputDEBUG, "\nFODO:");
	scrivimatr2D(F,outputDEBUG);
#endif

/************************************************************************/

//ora primi dell'iterazione mi calcolo le micromappe Li di lunghezza S=L/n

	double lunghezzatotale=0.;
	for (int i = 0 ;i < contatore;i++)
		lunghezzatotale+=lunghezza[i];

#ifdef DEBUG
	for (int i=0; i < contatore; i++)
	{
		fprintf(outputDEBUG,"\n#step in elemento %d = %d",i, dsMap(lunghezza[i],lunghezzatotale,N_STEP));
	}
	fprintf(outputDEBUG,"\n");
#endif

	double S = 0.;
	//Calcolo MICROMAPPE per il Drift
	for (int i=0; i < contatore; i++)
	{
		S=lunghezza[i]/dsMap(lunghezza[i],lunghezzatotale,N_STEP);
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

/***********************************************************************/


	double dl=0.;
	double lunghezza_accumulata=0.01;

	S=lunghezza[0]/dsMap(lunghezza[0],lunghezzatotale,N_STEP);
	for (int i=0;i<contatore;i++)
	{
		dl=lunghezza[i]/dsMap(lunghezza[i],lunghezzatotale,N_STEP);
		if (elemento[i]=='O')
		{
			fprintf(matrici_iniziali,"\n#Drift #%d, dl = %f",i,dl);
			fprintf(funzioni_ottiche,"\n#Drift #%d, dl = %f",i,dl);
			while(S<=(lunghezza_accumulata+lunghezza[i]))	
			{
				fprintf(matrici_iniziali,"\n\n Num_Step %f", S);
				scrivimatr2D(F,matrici_iniziali);
				prod(vett_i,O[i],S);
				scrivi_pos_part(posizionePart,vett_i,S);
				F=simil(F,OI[i],O[i]);
				alpha=optics(F,FOC);
				beta=optics(F,DEFOC);
				aminmax = assi_ellissi(alpha);
				bminmax = assi_ellissi(beta);
				scrividati(S,alpha,beta,aminmax,bminmax,funzioni_ottiche);
			#ifdef TURK
//				optics_T(alphaturk,FOC,O[i]);
//				optics_T(betaturk,DEFOC,O[i]);
				optics_T(alphaturk,FOC,F);
				optics_T(betaturk,DEFOC,F);
				aminmaxturk = assi_ellissi(alphaturk);
				bminmaxturk = assi_ellissi(betaturk);
				scrividati(S,alphaturk,betaturk,aminmaxturk,bminmaxturk,funzioni_otticheturk);
			#endif 
				S+=dl;
			}
			lunghezza_accumulata+=lunghezza[i];
		}
		else if (elemento[i]=='F')
		{
			fprintf(matrici_iniziali,"\n#Foc. #%d, dl = %f",i,dl);			
			fprintf(funzioni_ottiche,"\n#Foc. #%d, dl = %f",i,dl);
			while(S<=(lunghezza_accumulata+lunghezza[i]))	
			{
				fprintf(matrici_iniziali,"\n\n Num_Step %f", S);
				scrivimatr2D(F,matrici_iniziali);
				prod(vett_i,Fx[i],S);
				scrivi_pos_part(posizionePart,vett_i,S);
				F=simil(F,FxI[i],Fx[i]);
				alpha=optics(F,FOC);
				beta=optics(F,DEFOC);
				aminmax = assi_ellissi(alpha);
				bminmax = assi_ellissi(beta);
				scrividati(S,alpha,beta,aminmax,bminmax,funzioni_ottiche);
			#ifdef TURK
//				optics_T(alphaturk,FOC,Fx[i]);
//				optics_T(betaturk,DEFOC,Fx[i]);
				optics_T(alphaturk,FOC,F);
				optics_T(betaturk,DEFOC,F);
				aminmaxturk = assi_ellissi(alphaturk);
				bminmaxturk = assi_ellissi(betaturk);
				scrividati(S,alphaturk,betaturk,aminmaxturk,bminmaxturk,funzioni_otticheturk);
			#endif
				S+=dl;
			}
			lunghezza_accumulata+=lunghezza[i];
		}
		else if (elemento[i]=='D')
		{
			fprintf(matrici_iniziali,"\n#Defoc. #%d, dl = %f",i,dl);
			fprintf(funzioni_ottiche,"\n#Defoc. #%d, dl = %f",i,dl);
			while (S<=(lunghezza_accumulata+lunghezza[i]))
			{
				fprintf(matrici_iniziali,"\n\n Num_Step %f", S);
				scrivimatr2D(F,matrici_iniziali);
				prod(vett_i,Dx[i],S);
				scrivi_pos_part(posizionePart,vett_i,S);
				F=simil(F,DxI[i],Dx[i]);
				alpha=optics(F,FOC);
				beta=optics(F,DEFOC);
				aminmax = assi_ellissi(alpha);
				bminmax = assi_ellissi(beta);
				scrividati(S,alpha,beta,aminmax,bminmax,funzioni_ottiche);
			#ifdef TURK
//				optics_T(alphaturk,FOC,Dx[i]);
//				optics_T(betaturk,DEFOC,Dx[i]);
				optics_T(alphaturk,FOC,F);
				optics_T(betaturk,DEFOC,F);
				aminmaxturk = assi_ellissi(alphaturk);
				bminmaxturk = assi_ellissi(betaturk);
				scrividati(S,alphaturk,betaturk,aminmaxturk,bminmaxturk,funzioni_otticheturk);
			#endif 
				S+=dl;
			}
			lunghezza_accumulata+=lunghezza[i];
		}
	}


	create_gnuplot_file( "Posizione.plt", "Posizione_Particelle", lunghezza, contatore, 1 ,0.0, lunghezza_accumulata);
	create_gnuplot_file( "Funzioni_Ottiche.plt", "Funzioni_Ottiche", lunghezza, contatore, 1 ,0.0, lunghezza_accumulata);
	create_gnuplot_file( "Funzioni_Ottiche_Tuchetti.plt", "Funzioni_Ottiche_T", lunghezza, contatore, 1 ,0.0, lunghezza_accumulata);

	fclose(funzioni_ottiche);
	fclose(matrici_iniziali);
	fclose(posizionePart);
	parametri.close();

#ifdef TURK
	fclose(funzioni_otticheturk);
#endif

#ifdef DEBUG
	fclose(outputDEBUG);
#endif

}

