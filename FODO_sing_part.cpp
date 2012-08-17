
#include "FODO_sing_part.h"


int dsMap(double L, double lunghezzatotale, int n_step)
{
	double ds = lunghezzatotale/n_step;
	int n_p = (int) floor(L/ds);
	return n_p;
}


#ifdef DEBUG
void prod(double * A, vector< vector <double> > N, double S)
#else
void prod(double * A, vector< vector <double> > N)
#endif
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


double *optics(vector< vector <double> > N, int i,bool * effettuato_con_successo)
{
	double omega, s, *A=new double[2];
	for (int j = 0; j < 2; j++) A[j] = 0.;
	if (fabs((N[i][i]+N[i+1][i+1])*0.5) > 1.)
	{
		printf("Errore impossibile calcolare le funzioni ottiche!\n");
		* effettuato_con_successo=false;
		return A;		// ritorna zero e basta: al limite dobbiamo lavorarci su in altro modo, fare gli ifndef non funziona perche' sono macro per il preprocessore, non variabili valutate in fase di runtime!
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

#ifdef DEBUG
vector< vector <double> >  defocusing (vector< vector <double> > M, double d1, double S, FILE * file, int contatore)
#else
vector< vector <double> >  defocusing (vector< vector <double> > M, double d1, double S)
#endif
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

#ifdef DEBUG
vector< vector <double> >  focusing (vector< vector <double> > M, double d1, double S, FILE * file, int contatore)
#else
vector< vector <double> >  focusing (vector< vector <double> > M, double d1, double S)
#endif
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


#ifdef DEBUG
vector< vector <double> >  drift (vector< vector <double> > M, double S, FILE * file, int contatore)
#else
vector< vector <double> >  drift (vector< vector <double> > M, double S)
#endif
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
	fprintf(file," %+10.5f",A[1]);			// alpha_x
	fprintf(file," %+10.5f",B[1]);			// beta_x
	fprintf(file," %+10.5f",A[0]);			// alpha_y
	fprintf(file," %+10.5f",B[0]);			// beta_y
	fprintf(file," %+10.5f",AMaxMin[1]);	// proiezione ellisse su asse x
	fprintf(file," %+10.5f",BMaxMin[1]);	// proiezione ellisse su asse p_x
	fprintf(file," %+10.5f",AMaxMin[0]);	// proiezione ellisse su asse y
	fprintf(file," %+10.5f",BMaxMin[0]);	// proiezione ellisse su asse p_y
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

double * assi_ellissi(double *alpha, double emittance)
{
	double * minimi_massimi = new double[2];
	for (int i = 0; i < 2; i++) minimi_massimi[i] = 0.0;
	minimi_massimi[0] = sqrt((alpha[1])*emittance);
	minimi_massimi[1] = sqrt(((alpha[0])* (alpha[0])+1.)*emittance / (alpha[1]));
	return minimi_massimi;
}

void create_gnuplot_file(string gnuplot_filename, string run_name, double *lunghezza, int contatore, double estremo, double zmin, double zmax, string *keys)
{
	ofstream gnuplot_file;
	double lunghezza_percorsa=0.;
	gnuplot_file.open(gnuplot_filename.c_str());
	gnuplot_file << "#!/gnuplot" << endl;
	gnuplot_file << "FILE=\"" << run_name << ".txt\"" << endl;
#if defined (CREATE_PNG)
	gnuplot_file << "set terminal png enhanced 15" << endl;
	gnuplot_file << "set output \"graph_" << run_name << ".png\"" << endl;
#elif defined (CREATE_EPS)
	gnuplot_file << "set terminal postscript eps enhanced colour solid rounded \"Helvetica\" 25" << endl;
	gnuplot_file << "set output \"graph_" << run_name << ".eps\"" << endl;
#endif
	gnuplot_file << "set xrange[" << zmin << ":" << zmax << "]" << endl;
	gnuplot_file << "#set yrange[" << -estremo << ":" << estremo << "]" << endl;
	gnuplot_file << "set title  \""<< run_name <<" \"" << endl;
	gnuplot_file << "set xlabel \" " << keys[0] << "\"" << endl;
	gnuplot_file << "set ylabel \" " << keys[1] << "\"" << endl;
	for (int i = 0; i < contatore; i++)
	{
		lunghezza_percorsa+=lunghezza[i];
		gnuplot_file << "set arrow from " << lunghezza_percorsa<< "," << -5 << " to "<< lunghezza_percorsa << ","<< 5 << " nohead lc rgb \"black\" lw 1" << endl;
	}
	gnuplot_file << "plot FILE u 1:2 w lines lt 1 lc rgb \"red\" lw 1 t \" " << keys[2] << "\",\\" << endl;
	gnuplot_file << "FILE u 1:4 w lines lt 1 lc rgb \"blue\" lw 1 t \" " << keys[3] << "\",\\" << endl;
	gnuplot_file << "FILE u 1:3 w lines lt 1 lc rgb \"orange\" lw 1 t \" " << keys[4] << "\",\\" << endl;
	gnuplot_file << "FILE u 1:5 w lines lt 1 lc rgb \"dark-green\" lw 1 t \" " << keys[5] << "\"" << endl;

	gnuplot_file.close();
}

#ifdef TEST_OPTICAL_FUNCTIONS
double *optics_T (double * A, int i, vector< vector <double> > O)
{
	double alpha=0.;
	double beta=0.;
	alpha = A[0];
	beta = A[1];

	A[0] = alpha - ((O[i][i]) * (O[i+1][i]) * beta)  +  (2. * (O[i+1][i]) * (O[i][i+1]) * alpha)  -  (1./beta) * (O[i][i+1]) * (O[i+1][i+1]) * (1. - alpha * alpha);
	A[1]= (O[i][i]) * (O[i][i]) * beta  - 2. * (O[i][i]) * (O[i][i+1]) * alpha  -  (1./beta) * (O[i][i+1]) * (O[i][i+1]) * (1. - (alpha * alpha));
	return A;
}
#endif



int main(int argc, char *argv[]) 
{ 
#ifdef __linux
	feenableexcept(2);
	feenableexcept(3);
#endif

	bool fallita_lettura_parametri = true;
	bool fallita_lettura_inputdistr = true;
	bool do_transport = false;
	bool do_optics = false;
	ifstream parametri;
	ifstream inputdistr;
	int nstep = 1;

	for (int i = 1; i < argc; i++)
	{
		if (string(argv[i]) == "-p")
		{
			parametri.open(argv[i+1]);
			fallita_lettura_parametri=parametri.fail();
			i++;
		}
		else if (string(argv[i]) == "-i")
		{
			inputdistr.open(argv[i+1]);
			fallita_lettura_inputdistr=inputdistr.fail();
			i++;
		}
		else if (string(argv[i]) == "-transport")
		{
			do_transport=true;
		}
		else if (string(argv[i]) == "-optics")
		{
			do_optics=true;
		}
		else if (string(argv[i]) == "-nstep")
		{
			nstep = atoi(argv[i+1]);
			i++;
		}
		else
		{
			printf("Impossibile riconoscere il parametro %s\n",std::string(argv[i]));
		}
	}

	FILE * matrici_iniziali=fopen("Matrici_Iniziali.txt","w");
	FILE * posizionePart=fopen("Posizione_Particelle.txt","w");

	FILE * funzioni_ottiche=fopen("Funzioni_Ottiche.txt","w");
#ifdef TEST_OPTICAL_FUNCTIONS
	FILE * funzioni_ottiche_t=fopen("Funzioni_Ottiche_T.txt","w");
#endif

#ifdef DEBUG
	FILE * outputDEBUG=fopen("DEBUG.txt","w");
#endif

	string utile_per_contare;
	int conta_righe_parametri = 0;
	if (fallita_lettura_parametri || fallita_lettura_inputdistr)
	{
		printf("Impossibile aprire (o non definito) il file contenente i parametri\no il file contenente la distribuzione/particella iniziale\n");
		exit(204);
	}

	double * dati_iniziali = new double[6];	// emittanza, energia, x, y, px, py
	for (int i = 0; i < 6; i++)
	{
		if(inputdistr.eof())
		{
			printf("Mancano dei dati iniziali!\n");
			exit(123);
		}
		inputdistr >> dati_iniziali[i];
	}
	inputdistr.clear();
	inputdistr.seekg(0,std::ios::beg);

	double emittanza = dati_iniziali[0];
	double energia = dati_iniziali[1];
	double *vett_i=new double[4];
	vett_i[0]=dati_iniziali[2];
	vett_i[1]=dati_iniziali[4];
	vett_i[2]=dati_iniziali[3];
	vett_i[3]=dati_iniziali[5];

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

	double *alpha = new double[2];
	double *beta = new double[2];
	double *aminmax = new double[2];
	double *bminmax = new double[2];
	bool alpha_calcolato_con_successo=true;
	bool beta_calcolato_con_successo=true;

	double gamma_beta=sqrt(2.0*energia/MP_MEV);
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
#ifdef DEBUG
		if (elemento[i]=='F')
			Fx[i]=focusing(Fx[i],f1[i],lunghezza[i],matrici_iniziali,i);
		else if (elemento[i]=='D')
			Dx[i]=defocusing(Dx[i],d1[i],lunghezza[i],matrici_iniziali,i);
		else if (elemento[i]=='O')
			O[i]=drift(O[i],lunghezza[i],matrici_iniziali,i);
		else
			fprintf(outputDEBUG,"Elemento[%d]= %c non riconosciuto\n", i,elemento[i]);
#else
		if (elemento[i]=='F')
			Fx[i]=focusing(Fx[i],f1[i],lunghezza[i]);
		else if (elemento[i]=='D')
			Dx[i]=defocusing(Dx[i],d1[i],lunghezza[i]);
		else if (elemento[i]=='O')
			O[i]=drift(O[i],lunghezza[i]);
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


	if (do_optics)
	{
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


//		Calcolo Funzioni OTTICHE

		for (int i = 0; i < 2; i++) alpha[i] = beta[i] = aminmax[i] = bminmax[i] = 0.;
#ifdef TEST_OPTICAL_FUNCTIONS
		double *alphaturk = new double[2];
		double *betaturk = new double[2];
		double *aminmaxturk = new double[2];
		double *bminmaxturk = new double[2];
		for (int i = 0; i < 2; i++) alphaturk[i] = betaturk[i] = aminmaxturk[i] = bminmaxturk[i] = 0.;
#endif

		alpha=optics(F,FOC,&alpha_calcolato_con_successo);
		beta=optics(F,DEFOC,&beta_calcolato_con_successo);
		aminmax = assi_ellissi(alpha, emittanza);
		bminmax = assi_ellissi(beta, emittanza);

		if (alpha_calcolato_con_successo&&beta_calcolato_con_successo)
		{
			fprintf(funzioni_ottiche,"\n#%7c",'S');
			fprintf(funzioni_ottiche,"%10.8s","Alpha x");
			fprintf(funzioni_ottiche,"%10.7s","Beta x");
			fprintf(funzioni_ottiche,"%12.8s","Alpha y");
			fprintf(funzioni_ottiche,"%10.7s","Beta y");
			fprintf(funzioni_ottiche,"%10s","x");
			fprintf(funzioni_ottiche,"%11s","p_x");
			fprintf(funzioni_ottiche,"%11s","y");
			fprintf(funzioni_ottiche,"%11s","p_y");
		}

		scrividati(0.0,alpha,beta,aminmax,bminmax,funzioni_ottiche);


#ifdef TEST_OPTICAL_FUNCTIONS
		// non credo sia giusto inizializzare alphaturk e betaturk con l'altra funzione ottica
		// ma le "*turk" sono ricorsive e non permettono un bootstrap per ora...
		for (int i=0;i<2;i++)
		{
			alphaturk[i]=alpha[i];
			betaturk[i]=beta[i];
			aminmaxturk[i]=aminmax[i];
			bminmaxturk[i]=bminmax[i];
		}
		scrividati(0.0,alphaturk,betaturk,aminmaxturk,bminmaxturk,funzioni_ottiche_t);
#endif
	}

	if(do_transport)
	{
		fprintf(posizionePart," %+10.5f",0.0);
		fprintf(posizionePart," %+10.5f",vett_i[0]);
		fprintf(posizionePart," %+10.5f",vett_i[1]);
		fprintf(posizionePart," %+10.5f",vett_i[2]);
		fprintf(posizionePart," %+10.5f\n",vett_i[3]);
	}
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
		fprintf(outputDEBUG,"\n#step in elemento %d = %d",i, dsMap(lunghezza[i],lunghezzatotale,nstep));
	}
	fprintf(outputDEBUG,"\n");
#endif

	double S = 0.;
	//Calcolo MICROMAPPE per il Drift
	for (int i=0; i < contatore; i++)
	{
		S=lunghezza[i]/dsMap(lunghezza[i],lunghezzatotale,nstep);
#ifdef DEBUG
		if (elemento[i] == 'O')
		{
			O[i]=drift(O[i],S,matrici_iniziali,i);
			OI[i]=drift(OI[i],-S,matrici_iniziali,i);
		}
		else if (elemento[i] == 'F')
		{
			Fx[i]=focusing(Fx[i],f1[i],S,matrici_iniziali,i);
			FxI[i]=focusing(FxI[i],f1[i],-S,matrici_iniziali,i);
		}
		else if (elemento[i] == 'D')
		{
			Dx[i]=defocusing(Dx[i],d1[i],S,matrici_iniziali,i);
			DxI[i]=defocusing(DxI[i],d1[i],-S,matrici_iniziali,i);
		}
#else
		if (elemento[i] == 'O')
		{
			O[i]=drift(O[i],S);
			OI[i]=drift(OI[i],-S);
		}
		else if (elemento[i] == 'F')
		{
			Fx[i]=focusing(Fx[i],f1[i],S);
			FxI[i]=focusing(FxI[i],f1[i],-S);
		}
		else if (elemento[i] == 'D')
		{
			Dx[i]=defocusing(Dx[i],d1[i],S);
			DxI[i]=defocusing(DxI[i],d1[i],-S);
		}
#endif
	}

/***********************************************************************/


	double dl=0.;
//	double lunghezza_accumulata=0.01;
	double lunghezza_accumulata=lunghezza[0]/dsMap(lunghezza[0],lunghezzatotale,nstep);

	S=lunghezza[0]/dsMap(lunghezza[0],lunghezzatotale,nstep);
	for (int i=0;i<contatore;i++)
	{
		dl=lunghezza[i]/dsMap(lunghezza[i],lunghezzatotale,nstep);
		if (elemento[i]=='O')
		{
			fprintf(matrici_iniziali,"\n#Drift #%d, dl = %f",i,dl);
			fprintf(funzioni_ottiche,"\n#Drift #%d, dl = %f",i,dl);
			while(S<=(lunghezza_accumulata+lunghezza[i]))	
			{
				fprintf(matrici_iniziali,"\n\n Num_Step %f", S);
				scrivimatr2D(F,matrici_iniziali);
				if (do_transport)
				{
#ifdef DEBUG
					prod(vett_i,O[i],S);
#else
					prod(vett_i,O[i]);
#endif
					scrivi_pos_part(posizionePart,vett_i,S);
				}
				if (do_optics)
				{
					F=simil(F,OI[i],O[i]);
					alpha=optics(F,FOC,&alpha_calcolato_con_successo);
					beta=optics(F,DEFOC,&beta_calcolato_con_successo);
					aminmax = assi_ellissi(alpha, emittanza);
					bminmax = assi_ellissi(beta, emittanza);
					scrividati(S,alpha,beta,aminmax,bminmax,funzioni_ottiche);
#ifdef TEST_OPTICAL_FUNCTIONS
					alphaturk=optics_T(alphaturk,FOC,O[i]);
					betaturk=optics_T(betaturk,DEFOC,O[i]);
					aminmaxturk = assi_ellissi(alphaturk, emittanza);
					bminmaxturk = assi_ellissi(betaturk, emittanza);
					scrividati(S,alphaturk,betaturk,aminmaxturk,bminmaxturk,funzioni_ottiche_t);
#endif
				}
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
				if (do_transport)
				{
#ifdef DEBUG
					prod(vett_i,Fx[i],S);
#else
					prod(vett_i,Fx[i]);
#endif
					scrivi_pos_part(posizionePart,vett_i,S);
				}
				if (do_optics)
				{
					F=simil(F,FxI[i],Fx[i]);
					alpha=optics(F,FOC,&alpha_calcolato_con_successo);
					beta=optics(F,DEFOC,&beta_calcolato_con_successo);		
					aminmax = assi_ellissi(alpha, emittanza);
					bminmax = assi_ellissi(beta, emittanza);
					scrividati(S,alpha,beta,aminmax,bminmax,funzioni_ottiche);
#ifdef TEST_OPTICAL_FUNCTIONS
					alphaturk=optics_T(alphaturk,FOC,Fx[i]);
					betaturk=optics_T(betaturk,DEFOC,Fx[i]);
					aminmaxturk = assi_ellissi(alphaturk, emittanza);
					bminmaxturk = assi_ellissi(betaturk, emittanza);
					scrividati(S,alphaturk,betaturk,aminmaxturk,bminmaxturk,funzioni_ottiche_t);
#endif
				}
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
				if (do_transport)
				{
#ifdef DEBUG
					prod(vett_i,Dx[i],S);
#else
					prod(vett_i,Dx[i]);
#endif
					scrivi_pos_part(posizionePart,vett_i,S);
				}
				if (do_optics)
				{
					F=simil(F,DxI[i],Dx[i]);
					alpha=optics(F,FOC,&alpha_calcolato_con_successo);
					beta=optics(F,DEFOC,&beta_calcolato_con_successo);
					aminmax = assi_ellissi(alpha, emittanza);
					bminmax = assi_ellissi(beta, emittanza);
					scrividati(S,alpha,beta,aminmax,bminmax,funzioni_ottiche);
#ifdef TEST_OPTICAL_FUNCTIONS
					alphaturk=optics_T(alphaturk,FOC,Dx[i]);
					betaturk=optics_T(betaturk,DEFOC,Dx[i]);
					aminmaxturk = assi_ellissi(alphaturk, emittanza);
					bminmaxturk = assi_ellissi(betaturk, emittanza);
					scrividati(S,alphaturk,betaturk,aminmaxturk,bminmaxturk,funzioni_ottiche_t);
#endif
				}
				S+=dl;
			}
			lunghezza_accumulata+=lunghezza[i];
		}
	}

	fclose(funzioni_ottiche);
	fclose(matrici_iniziali);
	fclose(posizionePart);
	parametri.close();
	inputdistr.close();

#ifdef TEST_OPTICAL_FUNCTIONS
	fclose(funzioni_ottiche_t);
#endif

	string *etichette_posizione = new string[6];
	string *etichette_ottiche = new string[6];
	etichette_posizione[0] = "z (m)";
	etichette_posizione[1] = "x/y (m), p_x/p_y";
	etichette_posizione[2] = "x";
	etichette_posizione[3] = "y";
	etichette_posizione[4] = "p_x";
	etichette_posizione[5] = "p_y";

	etichette_ottiche[0] = "z (m)";
#if defined (CREATE_EPS)
	etichette_ottiche[1] = "{/Symbol a}, {/Symbol b}";
	etichette_ottiche[2] = "{/Symbol a}_x";
	etichette_ottiche[3] = "{/Symbol a}_y";
	etichette_ottiche[4] = "{/Symbol b}_x";
	etichette_ottiche[5] = "{/Symbol b}_y";
#else
	etichette_ottiche[1] = "Alpha, Beta";
	etichette_ottiche[2] = "Alpha_x";
	etichette_ottiche[3] = "Alpha_y";
	etichette_ottiche[4] = "Beta_x";
	etichette_ottiche[5] = "Beta_y";
#endif

	if (do_transport)
	{
		create_gnuplot_file( "Posizione.plt", "Posizione Particelle", lunghezza, contatore, 1 ,0.0, lunghezza_accumulata, etichette_posizione);
		system ("gnuplot Posizione.plt");
	}
	
	if (do_optics&&alpha_calcolato_con_successo&&beta_calcolato_con_successo)
	{
		create_gnuplot_file( "Funzioni_Ottiche.plt", "Funzioni Ottiche", lunghezza, contatore, 1 ,0.0, lunghezza_accumulata, etichette_ottiche);
		system ("gnuplot Funzioni_Ottiche.plt");
#ifdef TEST_OPTICAL_FUNCTIONS
		create_gnuplot_file( "Funzioni_Ottiche_T.plt", "Funzioni Ottiche T", lunghezza, contatore, 1 ,0.0, lunghezza_accumulata, etichette_ottiche);
		system ("gnuplot Funzioni_Ottiche_T.plt");
#endif
	}
	else	
	{
#if defined (__linux)
		system ("rm Funzioni_Ottiche.txt"); 
#ifdef TEST_OPTICAL_FUNCTIONS
		system ("rm Funzioni_Ottiche_T.txt");
#endif
#elif defined (_WIN32) || defined (_WIN64)
		system ("del Funzioni_Ottiche.txt"); 
#ifdef TEST_OPTICAL_FUNCTIONS
		system ("del Funzioni_Ottiche_T.txt");
#endif
#endif
	} 
	

#ifdef DEBUG
	fclose(outputDEBUG);
#endif

}

