
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
	//if (fabs((N[i][i]+N[i+1][i+1])*0.5) > 1.)
	//{
	//	printf("Errore impossibile calcolare le funzioni ottiche!\n");
	//	* effettuato_con_successo=false;
	//	return A;		// ritorna zero e basta: al limite dobbiamo lavorarci su in altro modo, fare gli ifndef non funziona perche' sono macro per il preprocessore, non variabili valutate in fase di runtime!
	//}
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

	*effettuato_con_successo=true;
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


void scrividati (double s, double *A, double *B, FILE * file)
{
	fprintf(file,"\n %+7.4f",s);
	fprintf(file," %+10.5f",A[1]);			// alpha_x
	fprintf(file," %+10.5f",B[1]);			// beta_x
	fprintf(file," %+10.5f",A[0]);			// alpha_y
	fprintf(file," %+10.5f",B[0]);			// beta_y
}


void scrividati_ellissi (double s,double *AMaxMin, double *BMaxMin, FILE * file)
{
	fprintf(file,"\n %+7.4f",s);
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
	minimi_massimi[1] = (alpha[0] * emittance / (minimi_massimi[0]));
	return minimi_massimi;
}


void create_gnuplot_file(string gnuplot_filename, string data_filename, double *lunghezza, int contatore, double estremo, double zmin, double zmax, string *keys,double *dati_rilevati,int conta_per_confronto)
{
	ofstream gnuplot_file;
	double lunghezza_percorsa=0.;
	gnuplot_file.open(gnuplot_filename.c_str());
	gnuplot_file << "#!/gnuplot" << endl;
	gnuplot_file << "FILE=\"" << data_filename << "\"" << endl;
#if defined (CREATE_PNG)
	gnuplot_file << "set terminal png enhanced 15" << endl;
	gnuplot_file << "set output \"graph_" << keys[0] << ".png\"" << endl;
#elif defined (CREATE_EPS)
	gnuplot_file << "set terminal postscript eps enhanced colour solid rounded \"Helvetica\" 25" << endl;
	gnuplot_file << "set output \"graph_" << run_name << ".eps\"" << endl;
#endif
	gnuplot_file << "set xrange[" << zmin << ":" << zmax << "]" << endl;
	gnuplot_file << "set yrange[" << -estremo << ":" << estremo << "]" << endl;
	gnuplot_file << "set title  \""<< keys[1] <<" \"" << endl;
	gnuplot_file << "set xlabel \" " << keys[2] << "\"" << endl;
	gnuplot_file << "set ylabel \" " << keys[3] << "\"" << endl;
	for (int i = 0; i < contatore; i++)
	{
		lunghezza_percorsa+=lunghezza[i];
		gnuplot_file << "set arrow from " << lunghezza_percorsa<< "," << -estremo << " to "<< lunghezza_percorsa << ","<< estremo << " nohead lc rgb \"black\" lw 1" << endl;
	}
	for (int a = 0; a < conta_per_confronto ; a++)
		gnuplot_file << "set arrow from " << dati_rilevati[a]<< "," << -estremo << " to "<< dati_rilevati[a] << ","<< estremo << " nohead lc rgb \"gray\" lw 1" << endl;
	gnuplot_file << "plot FILE u 1:2 w lines lt 1 lc rgb \"red\" lw 1 t \" " << keys[4] << "\",\\" << endl;
	gnuplot_file << "FILE u 1:4 w lines lt 1 lc rgb \"blue\" lw 1 t \" " << keys[5] << "\",\\" << endl;
	gnuplot_file << "FILE u 1:3 w lines lt 1 lc rgb \"orange\" lw 1 t \" " << keys[6] << "\",\\" << endl;
	gnuplot_file << "FILE u 1:5 w lines lt 1 lc rgb \"dark-green\" lw 1 t \" " << keys[7] << "\"" << endl;
	gnuplot_file.close();
}


#ifdef TEST_OPTICAL_FUNCTIONS
void optics_T (double ** A, int i, vector< vector <double> > O)
{
	double alpha=0.;
	double beta=0.;
	alpha = *A[0];
	beta = *A[1];

	*A[0] = alpha - ((O[i][i]) * (O[i+1][i]) * beta)  +  (2. * (O[i+1][i]) * (O[i][i+1]) * alpha)  -  (1./beta) * (O[i][i+1]) * (O[i+1][i+1]) * (1. - alpha * alpha);
	*A[1]= (O[i][i]) * (O[i][i]) * beta  - 2. * (O[i][i]) * (O[i][i+1]) * alpha  -  (1./beta) * (O[i][i+1]) * (O[i][i+1]) * (1. - (alpha * alpha));
}
//double *optics_T (double * A, int i, vector< vector <double> > O)
//{
//	double alpha=0.;
//	double beta=0.;
//	alpha = A[0];
//	beta = A[1];
//
//	A[0] = alpha - ((O[i][i]) * (O[i+1][i]) * beta)  +  (2. * (O[i+1][i]) * (O[i][i+1]) * alpha)  -  (1./beta) * (O[i][i+1]) * (O[i+1][i+1]) * (1. - alpha * alpha);
//	A[1]= (O[i][i]) * (O[i][i]) * beta  - 2. * (O[i][i]) * (O[i][i+1]) * alpha  -  (1./beta) * (O[i][i+1]) * (O[i][i+1]) * (1. - (alpha * alpha));
//	return A;
//}
//double *optics_T (double * A, int i, vector< vector <double> > O,double lunghezza, double ds)
//{
//	double *alpha=new double[2];
//	double *beta=new double[2];
//	double S=0.;
//	alpha[0] = A[0];
//	beta[0] = A[1];
//	while (lunghezza>=ds)
//	{
//		alpha[1] = alpha[0] - ((O[i][i]) * (O[i+1][i]) * beta[0])  +  (2. * (O[i+1][i]) * (O[i][i+1]) * alpha[0])  -  (1./beta[0]) * (O[i][i+1]) * (O[i+1][i+1]) * (1. - alpha[0] * alpha[0]);
//		beta[1]= (O[i][i]) * (O[i][i]) * beta[0]  - 2. * (O[i][i]) * (O[i][i+1]) * alpha[0]  -  (1./beta[0]) * (O[i][i+1]) * (O[i][i+1]) * (1. - (alpha[0] * alpha[0]));
//		alpha[0]=alpha[1];
//		beta[0]=beta[1];
//		S+=ds;
//	}	
//	A[0]=alpha[0];
//	A[1]=beta[0];
//	return A;
//}
#endif

void massimo_opt(double * optics_x , double * optics_y,double * massimo_temp)
{
	double max= *massimo_temp;
	max>fabs(optics_x[0])?max:max=fabs(optics_x[0]);
	max>fabs(optics_x[1])?max:max=fabs(optics_x[1]);
	max>fabs(optics_y[0])?max:max=fabs(optics_y[0]);
	max>fabs(optics_y[1])?max:max=fabs(optics_y[1]);
	* massimo_temp=max;
}


void massimo_pos(double * vett_i, double * massimo_temp)
{
	double max= *massimo_temp;
	max>fabs(vett_i[0])?max:max=fabs(vett_i[0]);
	max>fabs(vett_i[1])?max:max=fabs(vett_i[1]);
	max>fabs(vett_i[2])?max:max=fabs(vett_i[2]);
	max>fabs(vett_i[3])?max:max=fabs(vett_i[3]);
	* massimo_temp=max;
}

void confronto (double * parametri_dati,double * parametri_ottenuti,double z,/*double gradiente_f,*/double percentuale,FILE * file,bool *confronto_pos,/*string elemento,*/int *conto_per_confronto)
{
	double *diff=new double[2];
	
	//const char *elemento_char;
	//elemento_char = elemento.c_str();	
	diff[0]=fabs(parametri_dati[0]-parametri_ottenuti[0]);
	diff[1]=fabs(parametri_dati[1]-parametri_ottenuti[1]);
	if ((diff[0]<(percentuale*parametri_dati[0]))&&(diff[1]<(percentuale*parametri_dati[1])))
	{
		fprintf(file," \n%+10.5f ",z);
		//fprintf(file,"%s",elemento_char);
		//fprintf(file," %+10.5f ",gradiente_f);
		*confronto_pos=true;
		*conto_per_confronto+=1;
	}
}



int main(int argc, char *argv[]) 
{ 
#ifdef __linux
	feenableexcept(2);
	feenableexcept(3);
#endif

	FILE * matrici_iniziali=fopen("Matrici_Iniziali.txt","w");
	FILE * posizionePart=fopen("Posizione_Particelle.txt","w");
	FILE * ellissi=fopen("Parametri_Ellissi_Funz_Ottiche.txt","w");
	FILE * funzioni_ottiche=fopen("Funzioni_Ottiche.txt","w");
	FILE * confronti=fopen("Math_rilevati.txt","w");
#ifdef TEST_OPTICAL_FUNCTIONS
	FILE * funzioni_ottiche_t=fopen("Funzioni_Ottiche_T.txt","w");
	FILE * ellissi_t=fopen("Parametri_Ellissi_Funz_Ottiche_T.txt","w");
	FILE * confronti_t=fopen("Math_rilevati_T.txt","w");
#endif

#ifdef DEBUG
	FILE * outputDEBUG=fopen("DEBUG.txt","w");
#endif

	bool fallita_lettura_parametri = true;
	bool fallita_lettura_inputdistr = true;
	bool do_transport = false;
	bool do_optics = false;
	bool posso_fare_funzioni_ottiche = false;
	ifstream parametri;
	ifstream inputdistr;
	int nstep = 1;

	double gnuplot_ymax_opt=0.;
	bool calcola_ymax_opt = true;
#ifdef TEST_OPTICAL_FUNCTIONS
	double gnuplot_ymax_opt_T=0.;
	bool calcola_ymax_opt_T = true;
#endif
	double gnuplot_ymax_pos=0.;
	bool calcola_ymax_pos = true;
	double gnuplot_xmax_opt=0.;
	double gnuplot_xmax_pos=0.;
	bool calcola_ymax_ell = true;
	double gnuplot_ymax_ell=0.;
//	double gnuplot_xmax_ell=0.;
//	bool calcola_xmax_ell = true;
	double percentuale=0.03;
	int conto_per_confronto=0;

	double *compare_x=new double[2];
	double *compare_y=new double[2];
	bool confronto_pos_x=false;
	bool confronto_pos_y=false;

	double *paramIniz_X=new double[2];
	double *paramIniz_Y=new double[2];
	bool fai_da_te_x=false;
	bool fai_da_te_y=false;
#ifdef TEST_OPTICAL_FUNCTIONS
	int conto_per_confronto_t_x=0;
	int conto_per_confronto_t_y=0;
	bool confronto_pos_t_y=false;
	bool confronto_pos_t_x=false;
#endif

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
		else if (string(argv[i]) == "-xmax_opt")
		{
			gnuplot_xmax_opt=atof(argv[i+1]);
			i++;
		}
		else if (string(argv[i]) == "-xmax_pos")
		{
			gnuplot_xmax_pos=atof(argv[i+1]);
			i++;
		}
		else if (string(argv[i]) == "-ymax_opt")
		{
			gnuplot_ymax_opt=atof(argv[i+1]);
			calcola_ymax_opt = false;
			i++;
		}
#ifdef TEST_OPTICAL_FUNCTIONS
		else if (string(argv[i]) == "-ymax_opt_T")
		{
			gnuplot_ymax_opt_T=atof(argv[i+1]);
			calcola_ymax_opt_T = false;
			i++;
		}
#endif
		else if (string(argv[i]) == "-ymax_pos")
		{
			gnuplot_ymax_pos=atof(argv[i+1]);
			calcola_ymax_pos = false;
			i++;
		}
		else if (string(argv[i]) == "-ymax_ell")
		{
			gnuplot_ymax_ell=atof(argv[i+1]);
			calcola_ymax_ell = false;
			i++;
		}
		else if (string(argv[i]) == "-compare_X")
		{
			compare_x[0]=atof(argv[i+1]);
			compare_x[1]=atof(argv[i+2]);
			i+=2;
		}
		else if (string(argv[i]) == "-compare_Y")
		{
			compare_y[0]=atof(argv[i+1]);
			compare_y[1]=atof(argv[i+2]);

			i+=2;
		}
		else if (string(argv[i]) == "-perc")
		{
			percentuale=atof(argv[i+1]);
			fprintf(matrici_iniziali,"\n%f\n",percentuale);
			i++;
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
		else if (string(argv[i]) == "-paramIniz_X")
		{
			paramIniz_X[0] = atoi(argv[i+1]);
			paramIniz_X[1] = atoi(argv[i+2]);
			fai_da_te_x=true;
			i+=2;
		}
		else if (string(argv[i]) == "-paramIniz_Y")
		{
			paramIniz_Y[0] = atoi(argv[i+1]);
			paramIniz_Y[1] = atoi(argv[i+2]);
			fai_da_te_y=true;
			i+=2;
		}
		else
		{
			printf("Impossibile riconoscere il parametro %s\n",argv[i]);
		}
	}

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
//	char *elemento=new char[conta_righe_parametri];
	string *elemento=new string[conta_righe_parametri];
	double * lunghezza= new double[conta_righe_parametri];
	double * gradiente= new double[conta_righe_parametri];
	int contatore=0;
	for (int i = 0; i < conta_righe_parametri; i++)
	{
		parametri >> elemento[i];
		parametri >> gradiente[i];
		parametri >> lunghezza[i];
#ifdef DEBUG
		cout << "Tipo elemento: " << elemento[i] << ", gradiente: " << gradiente[i] << ", lunghezza: " << lunghezza[i] << endl;
//		printf("tipo elemento: %c, gradiente: %f, lunghezza: %f\n", elemento[contatore], gradiente[contatore], lunghezza[contatore]);
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
	bool alpha_calcolato_con_successo=false;
	bool beta_calcolato_con_successo=false;

#ifdef TEST_OPTICAL_FUNCTIONS
	double *ottiche_x_t = new double[2];
	double *ottiche_y_t = new double[2];
	double *aminmax_x_t = new double[2];
	double *bminmax_y_t = new double[2];
	for (int i = 0; i < 2; i++) ottiche_x_t[i] = ottiche_y_t[i] = aminmax_x_t[i] = bminmax_y_t[i] = 0.;
#endif


	double gamma_beta=sqrt(2.0*energia/MP_MEV);
	double gamma_v=gamma_beta*SPEED_OF_LIGHT;
	double *f1 =new double [contatore];
	double *d1 =new double [contatore];
	for (int i=0;i<contatore;i++)
	{
		if (elemento[i]=="F")
			f1[i]=sqrt(gradiente[i]*CHARGE/(MP_KG*gamma_v));
		if (elemento[i]=="D")
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
		if (elemento[i]=="F")
			Fx[i]=focusing(Fx[i],f1[i],lunghezza[i],matrici_iniziali,i);
		else if (elemento[i]=="D")
			Dx[i]=defocusing(Dx[i],d1[i],lunghezza[i],matrici_iniziali,i);
		else if (elemento[i]=="O")
			O[i]=drift(O[i],lunghezza[i],matrici_iniziali,i);
		else
			fprintf(outputDEBUG,"Elemento[%d] non riconosciuto\n", i);
#else
		if (elemento[i]=="F")
			Fx[i]=focusing(Fx[i],f1[i],lunghezza[i]);
		else if (elemento[i]=="D")
			Dx[i]=defocusing(Dx[i],d1[i],lunghezza[i]);
		else if (elemento[i]=="O")
			O[i]=drift(O[i],lunghezza[i]);
#endif
	}

#ifdef DEBUG

	for (int i=0;i<contatore;i++)
	{
		if (elemento[i]=="O")
		{
//			if (O[i][0][0] == 0.0) continue;
			fprintf(matrici_iniziali,"\nMATRICE DRIFT");
			scrivimatr2D(O[i],matrici_iniziali);
			fprintf(matrici_iniziali,"\n");
		}
		else if (elemento[i]=="F")
		{
//			if (Fx[i][0][0] == 0.0) continue;
			fprintf(matrici_iniziali,"\nMATRICE FOC.");
			scrivimatr2D(Fx[i],matrici_iniziali);
			fprintf(matrici_iniziali,"\n");
		}
		else if (elemento[i]=="D")
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

		if (elemento[0]=="O")
		{
			for (int k=0; k<4; k++)
				for (int j=0; j<4; j++)
					compos[k][j]=O[0][k][j];
		}
		else if (elemento[0]=="F")
		{
			for (int k=0; k<4; k++)
				for (int j=0; j<4; j++)
					compos[k][j]=Fx[0][k][j];
		}
		else if (elemento[0]=="D")
		{
			for (int k=0; k<4; k++)
				for (int j=0; j<4; j++)
					compos[k][j]=Dx[0][k][j];
		}

		for (int i=1;i<contatore;i++)
		{
			if (elemento[i]=="O")
				compos=prodo(O[i],compos,4);
			else if (elemento[i]=="F")
				compos=prodo(Fx[i],compos,4);
			else if (elemento[i]=="D")
				compos=prodo(Dx[i],compos,4);
		}
	
		for (int i=0;i<4;i++)
			for(int a=0;a<4;a++)
				F[i][a]=compos[i][a];

//		Calcolo Funzioni OTTICHE

		for (int i = 0; i < 2; i++) alpha[i] = beta[i] = aminmax[i] = bminmax[i] = 0.;

		if ( (fabs((F[FOC][FOC]+F[FOC+1][FOC+1])*0.5) <= 1.) && (fabs((F[DEFOC][DEFOC]+F[DEFOC+1][DEFOC+1])*0.5) <= 1.))
			posso_fare_funzioni_ottiche = true;
		//else cout << "Impossibile calcolare le funzioni ottiche!" << endl;
		//cout << "posso_fare_funzioni_ottiche="<<posso_fare_funzioni_ottiche<<endl;
		if (posso_fare_funzioni_ottiche)
		{
			alpha=optics(F,FOC,&alpha_calcolato_con_successo);
			beta=optics(F,DEFOC,&beta_calcolato_con_successo);
			aminmax = assi_ellissi(alpha, emittanza);
			bminmax = assi_ellissi(beta, emittanza);

			if (calcola_ymax_opt) massimo_opt(alpha,beta,&gnuplot_ymax_opt);
			if (calcola_ymax_pos) massimo_pos(vett_i,&gnuplot_ymax_pos);
			if (calcola_ymax_ell) massimo_opt(aminmax,bminmax,&gnuplot_ymax_ell);
//			if (calcola_xmax_ell) massimo_opt(aminmax,bminmax,&gnuplot_xmax_ell);
			fprintf(funzioni_ottiche,"# alpha_successo %d beta_successo %d\n",(int)(alpha_calcolato_con_successo),(int)(beta_calcolato_con_successo));
			if (alpha_calcolato_con_successo&&beta_calcolato_con_successo)
			{
				fprintf(funzioni_ottiche,"\n#%7c",'S');
				fprintf(funzioni_ottiche,"%10.8s","Alpha x");
				fprintf(funzioni_ottiche,"%10.7s","Beta x");
				fprintf(funzioni_ottiche,"%12.8s","Alpha y");
				fprintf(funzioni_ottiche,"%10.7s","Beta y");
				fprintf(ellissi,"%10s","x");
				fprintf(ellissi,"%11s","p_x");
				fprintf(ellissi,"%11s","y");
				fprintf(ellissi,"%11s","p_y");
			}

			scrividati(0.0,alpha,beta,funzioni_ottiche);
			scrividati_ellissi(0.0,aminmax,bminmax,ellissi);


#ifdef TEST_OPTICAL_FUNCTIONS
		// non credo sia giusto inizializzare alphaturk e betaturk con l'altra funzione ottica
		// ma le "*turk" sono ricorsive e non permettono un bootstrap per ora...
			for (int i=0;i<2;i++)
			{
				if (!(fai_da_te_x&&fai_da_te_y))
				{
					ottiche_x_t[i]=alpha[i];
					ottiche_y_t[i]=beta[i];
				}
				else if (fai_da_te_x)
				{
					ottiche_x_t[i]=paramIniz_X[i];
					ottiche_y_t[i]=paramIniz_X[i];

				}
				else if (fai_da_te_y)
				{
					ottiche_x_t[i]=paramIniz_Y[i];
					ottiche_y_t[i]=paramIniz_Y[i];
				}
				else
				{
					ottiche_x_t[i]=paramIniz_X[i];
					ottiche_y_t[i]=paramIniz_Y[i];

				}
			}
			aminmax_x_t = assi_ellissi(ottiche_x_t, emittanza);
			bminmax_y_t = assi_ellissi(ottiche_y_t, emittanza);
			scrividati(0.0,ottiche_x_t,ottiche_y_t,funzioni_ottiche_t);
			scrividati_ellissi(0.0,aminmax_x_t,bminmax_y_t,ellissi_t);
			if (calcola_ymax_opt_T) massimo_opt(ottiche_x_t,ottiche_y_t,&gnuplot_ymax_opt_T);
			//if (calcola_ymax_ell) massimo_opt(aminmaxturk,bminmaxturk,&gnuplot_ymax_ell);
#endif
		}
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
		if (elemento[i] == "O")
		{
			O[i]=drift(O[i],S,matrici_iniziali,i);
			OI[i]=drift(OI[i],-S,matrici_iniziali,i);
		}
		else if (elemento[i] == "F")
		{
			Fx[i]=focusing(Fx[i],f1[i],S,matrici_iniziali,i);
			FxI[i]=focusing(FxI[i],f1[i],-S,matrici_iniziali,i);
		}
		else if (elemento[i] == "D")
		{
			Dx[i]=defocusing(Dx[i],d1[i],S,matrici_iniziali,i);
			DxI[i]=defocusing(DxI[i],d1[i],-S,matrici_iniziali,i);
		}
#else
		if (elemento[i] == "O")
		{
			O[i]=drift(O[i],S);
			OI[i]=drift(OI[i],-S);
		}
		else if (elemento[i] == "F")
		{
			Fx[i]=focusing(Fx[i],f1[i],S);
			FxI[i]=focusing(FxI[i],f1[i],-S);
		}
		else if (elemento[i] == "D")
		{
			Dx[i]=defocusing(Dx[i],d1[i],S);
			DxI[i]=defocusing(DxI[i],d1[i],-S);
		}
#endif
	}

/***********************************************************************/


	double dl=0.;
	double lunghezza_accumulata=0.0;
	S=0.0;

	for (int i=0;i<contatore;i++)
	{
		dl=lunghezza[i]/dsMap(lunghezza[i],lunghezzatotale,nstep);
		if (elemento[i]=="O")
		{
			fprintf(matrici_iniziali,"\n#Drift #%d, dl = %f",i,dl);
			fprintf(funzioni_ottiche,"\n#Drift #%d, dl = %f",i,dl);
			while(S<=(lunghezza_accumulata+lunghezza[i]))
			//for(int l=0;l<=(dsMap(lunghezza[i],lunghezzatotale,nstep));l++)
			{
				S+=dl;
				fprintf(matrici_iniziali,"\n\n Num_Step %f", S);
				scrivimatr2D(F,matrici_iniziali);
				if (do_transport)
				{
#ifdef DEBUG
					prod(vett_i,O[i],S);
#else
					prod(vett_i,O[i]);
#endif
					if (calcola_ymax_pos) massimo_pos(vett_i,&gnuplot_ymax_pos);
					scrivi_pos_part(posizionePart,vett_i,S);
				}
				if (do_optics)
				{
					F=simil(F,OI[i],O[i]);
					if (posso_fare_funzioni_ottiche)
					{	
					alpha=optics(F,FOC,&alpha_calcolato_con_successo);
					beta=optics(F,DEFOC,&beta_calcolato_con_successo);
					aminmax = assi_ellissi(alpha, emittanza);
					bminmax = assi_ellissi(beta, emittanza);
					confronto(compare_x,aminmax,S,percentuale,confronti,&confronto_pos_x,&conto_per_confronto);
					confronto(compare_y,bminmax,S,percentuale,confronti,&confronto_pos_y,&conto_per_confronto);
					scrividati(S,alpha,beta,funzioni_ottiche);
					scrividati_ellissi(S,aminmax,bminmax,ellissi);
					if (calcola_ymax_ell) massimo_opt(aminmax,bminmax,&gnuplot_ymax_ell);
					if (calcola_ymax_opt) massimo_opt(alpha,beta,&gnuplot_ymax_opt);
				
#ifdef TEST_OPTICAL_FUNCTIONS
					optics_T(&ottiche_x_t,FOC,O[i]);
					optics_T(&ottiche_y_t,DEFOC,O[i]);
					aminmax_x_t = assi_ellissi(aminmax_x_t, emittanza);
					bminmax_y_t = assi_ellissi(bminmax_y_t, emittanza);
					confronto(paramIniz_X,aminmax_x_t,S,percentuale,confronti_t,&confronto_pos_t_x,&conto_per_confronto_t_x);
					confronto(paramIniz_Y,bminmax_y_t,S,percentuale,confronti_t,&confronto_pos_t_y,&conto_per_confronto_t_y);
					scrividati(S,ottiche_x_t,ottiche_y_t,funzioni_ottiche_t);
					scrividati_ellissi(S,aminmax_x_t,bminmax_y_t,ellissi_t);	
					//if (calcola_ymax_opt_T) massimo_opt(ottiche_x_t,ottiche_y_t,&gnuplot_ymax_opt_T);
					//if (calcola_ymax_ell_T) massimo_opt(aminmax_x_t,bminmax_y_t,&gnuplot_ymax_ell);
#endif

					}
				}
			}
			lunghezza_accumulata+=lunghezza[i];
		}
		else if (elemento[i]=="F")
		{
			fprintf(matrici_iniziali,"\n#Foc. #%d, dl = %f",i,dl);			
			fprintf(funzioni_ottiche,"\n#Foc. #%d, dl = %f",i,dl);
			while(S<=(lunghezza_accumulata+lunghezza[i]))
			//for(int l=0;l<=(dsMap(lunghezza[i],lunghezzatotale,nstep));l++)
			{
				S+=dl;
				fprintf(matrici_iniziali,"\n\n Num_Step %f", S);
				scrivimatr2D(F,matrici_iniziali);
				if (do_transport)
				{
#ifdef DEBUG
					prod(vett_i,Fx[i],S);
#else
					prod(vett_i,Fx[i]);
#endif
					if (calcola_ymax_pos) massimo_pos(vett_i,&gnuplot_ymax_pos);
					scrivi_pos_part(posizionePart,vett_i,S);
				}
				if (posso_fare_funzioni_ottiche)
				{
					F=simil(F,FxI[i],Fx[i]);
					alpha=optics(F,FOC,&alpha_calcolato_con_successo);
					beta=optics(F,DEFOC,&beta_calcolato_con_successo);		
					aminmax = assi_ellissi(alpha, emittanza);
					bminmax = assi_ellissi(beta, emittanza);
					confronto(compare_x,aminmax,S,percentuale,confronti,&confronto_pos_x,&conto_per_confronto);
					confronto(compare_y,bminmax,S,percentuale,confronti,&confronto_pos_y,&conto_per_confronto);
					scrividati(S,alpha,beta,funzioni_ottiche);
					scrividati_ellissi(S,aminmax,bminmax,ellissi);
					if (calcola_ymax_ell) massimo_opt(aminmax,bminmax,&gnuplot_ymax_ell);
					if (calcola_ymax_opt) massimo_opt(alpha,beta,&gnuplot_ymax_opt);

#ifdef TEST_OPTICAL_FUNCTIONS
					optics_T(&ottiche_x_t,FOC,Fx[i]);
					optics_T(&ottiche_y_t,DEFOC,Fx[i]);
					aminmax_x_t = assi_ellissi(aminmax_x_t, emittanza);
					bminmax_y_t = assi_ellissi(bminmax_y_t, emittanza);
					confronto(paramIniz_X,aminmax_x_t,S,percentuale,confronti_t,&confronto_pos_t_x,&conto_per_confronto_t_x);
					confronto(paramIniz_Y,bminmax_y_t,S,percentuale,confronti_t,&confronto_pos_t_y,&conto_per_confronto_t_y);
					scrividati(S,ottiche_x_t,ottiche_y_t,funzioni_ottiche_t);
					scrividati_ellissi(S,aminmax_x_t,bminmax_y_t,ellissi_t);
					//if (calcola_ymax_opt_T) massimo_opt(ottiche_x_t,ottiche_y_t,&gnuplot_ymax_opt_T);
					//if (calcola_ymax_ell_T) massimo_opt(aminmax_x_t,bminmax_y_t,&gnuplot_ymax_ell);
#endif

				}
			}
			lunghezza_accumulata+=lunghezza[i];
		}
		else if (elemento[i]=="D")
		{
			fprintf(matrici_iniziali,"\n#Defoc. #%d, dl = %f",i,dl);
			fprintf(funzioni_ottiche,"\n#Defoc. #%d, dl = %f",i,dl);
			while (S<=(lunghezza_accumulata+lunghezza[i]))
			//for(int l=0;l<=(dsMap(lunghezza[i],lunghezzatotale,nstep));l++)
			{
				S+=dl;
				fprintf(matrici_iniziali,"\n\n Num_Step %f", S);
				scrivimatr2D(F,matrici_iniziali);
				if (do_transport)
				{
#ifdef DEBUG
					prod(vett_i,Dx[i],S);
#else
					prod(vett_i,Dx[i]);
#endif
					if (calcola_ymax_pos) massimo_pos(vett_i,&gnuplot_ymax_pos);
					scrivi_pos_part(posizionePart,vett_i,S);
				}
				if (posso_fare_funzioni_ottiche)
				{
					F=simil(F,DxI[i],Dx[i]);
					alpha=optics(F,FOC,&alpha_calcolato_con_successo);
					beta=optics(F,DEFOC,&beta_calcolato_con_successo);
					aminmax = assi_ellissi(alpha, emittanza);
					bminmax = assi_ellissi(beta, emittanza);
					confronto(compare_x,aminmax,S,percentuale,confronti,&confronto_pos_x,&conto_per_confronto);
					confronto(compare_y,bminmax,S,percentuale,confronti,&confronto_pos_y,&conto_per_confronto);
					scrividati(S,alpha,beta,funzioni_ottiche);
					scrividati_ellissi(S,aminmax,bminmax,ellissi);
					if (calcola_ymax_ell) massimo_opt(aminmax,bminmax,&gnuplot_ymax_ell);
					if (calcola_ymax_opt) massimo_opt(alpha,beta,&gnuplot_ymax_opt);

#ifdef TEST_OPTICAL_FUNCTIONS
					optics_T(&ottiche_x_t,FOC,Dx[i]);
					optics_T(&ottiche_y_t,DEFOC,Dx[i]);
					aminmax_x_t = assi_ellissi(aminmax_x_t, emittanza);
					bminmax_y_t = assi_ellissi(bminmax_y_t, emittanza);
					confronto(paramIniz_X,aminmax_x_t,S,percentuale,confronti_t,&confronto_pos_t_x,&conto_per_confronto_t_x);
					confronto(paramIniz_Y,bminmax_y_t,S,percentuale,confronti_t,&confronto_pos_t_y,&conto_per_confronto_t_y);
					scrividati(S,ottiche_x_t,ottiche_y_t,funzioni_ottiche_t);
					scrividati_ellissi(S,aminmax_x_t,bminmax_y_t,ellissi_t);
					//if (calcola_ymax_opt_T) massimo_opt(ottiche_x_t,ottiche_y_t,&gnuplot_ymax_opt_T);
					//if (calcola_ymax_ell_T) massimo_opt(aminmax_x_t,bminmax_y_t,&gnuplot_ymax_ell);
#endif

				}
			}
			lunghezza_accumulata+=lunghezza[i];
		}
	}

	fclose(funzioni_ottiche);
	fclose(matrici_iniziali);
	fclose(posizionePart);
	fclose(ellissi);
	fclose(confronti);
	parametri.close();
	inputdistr.close();

#ifdef TEST_OPTICAL_FUNCTIONS
	fclose(funzioni_ottiche_t);
	fclose(ellissi_t);
#endif

	string *etichette_posizione = new string[8];
	string *etichette_ottiche = new string[8];
	string *etichette_ellissi = new string[8];
#ifdef TEST_OPTICAL_FUNCTIONS
	string *etichette_ottiche_T = new string[8];
#endif

	etichette_posizione[0] = "Posizione_Particelle";
	etichette_posizione[1] = "Posizione Particelle";
	etichette_posizione[2] = "z (m)";
	etichette_posizione[3] = "x/y (m), p_x/p_y";
	etichette_posizione[4] = "x";
	etichette_posizione[5] = "y";
	etichette_posizione[6] = "p_x";
	etichette_posizione[7] = "p_y";

	etichette_ottiche[0] = "Funzioni_Ottiche";
	etichette_ottiche[1] = "Funzioni Ottiche";
	etichette_ottiche[2] = "z (m)";
#if defined (CREATE_EPS)
	etichette_ottiche[3] = "{/Symbol a}, {/Symbol b}";
	etichette_ottiche[4] = "{/Symbol a}_x";
	etichette_ottiche[5] = "{/Symbol a}_y";
	etichette_ottiche[6] = "{/Symbol b}_x";
	etichette_ottiche[7] = "{/Symbol b}_y";
#else
	etichette_ottiche[3] = "Alpha, Beta";
	etichette_ottiche[4] = "Alpha_x";
	etichette_ottiche[5] = "Alpha_y";
	etichette_ottiche[6] = "Beta_x";
	etichette_ottiche[7] = "Beta_y";
#endif

	etichette_ellissi[0] = "Parametri_Ellissi_Funz_Ottiche";
	etichette_ellissi[1] = "Parametri Ellissi Funz Ottiche";
	etichette_ellissi[2] = "z (m)";
	etichette_ellissi[3] = "X , P";
	etichette_ellissi[4] = "Xmax";
	etichette_ellissi[5] = "Pmax_x";
	etichette_ellissi[6] = "Ymax";
	etichette_ellissi[7] = "Pmax_y";

#ifdef TEST_OPTICAL_FUNCTIONS
	etichette_ottiche_T[0] = "Funzioni_Ottiche_T";
	etichette_ottiche_T[1] = "Funzioni ottiche test";
	etichette_ottiche_T[2] = "z (m)";
#if defined (CREATE_EPS)
	etichette_ottiche_T[3] = "{/Symbol a}, {/Symbol b}";
	etichette_ottiche_T[4] = "{/Symbol a}_x";
	etichette_ottiche_T[5] = "{/Symbol a}_y";
	etichette_ottiche_T[6] = "{/Symbol b}_x";
	etichette_ottiche_T[7] = "{/Symbol b}_y";
#else
	etichette_ottiche_T[3] = "Alpha, Beta";
	etichette_ottiche_T[4] = "Alpha_x";
	etichette_ottiche_T[5] = "Alpha_y";
	etichette_ottiche_T[6] = "Beta_x";
	etichette_ottiche_T[7] = "Beta_y";
#endif
#endif

/***********************************************************/

	//cout << "conto_per_confronto= "<<conto_per_confronto;
	double *dati_rilevati=new double [conto_per_confronto];
	for (int a=0;a<conto_per_confronto;a++)
		dati_rilevati[a]=0.;
	if (confronto_pos_x||confronto_pos_y)
	{
		ifstream confro;
		confro.open("Math_rilevati.txt");
		for (int i=0;i<conto_per_confronto;i++)
		{
			confro >> dati_rilevati[i];
		}
	}

/***********************************************************/


	if (do_transport&&(confronto_pos_x||confronto_pos_y))
	{
		if (gnuplot_xmax_pos > 0.) create_gnuplot_file( "Posizione.plt", "Posizione_Particelle.txt", lunghezza, contatore, gnuplot_ymax_pos ,0.0, gnuplot_xmax_pos, etichette_posizione,dati_rilevati,conto_per_confronto);
		else create_gnuplot_file( "Posizione.plt", "Posizione_Particelle.txt", lunghezza, contatore, gnuplot_ymax_pos ,0.0, lunghezza_accumulata, etichette_posizione,dati_rilevati,conto_per_confronto);
		system ("gnuplot Posizione.plt");
	}

	
	//cout << "posso_fare_funzioni_ottiche="<<posso_fare_funzioni_ottiche<<endl;
	//cout << "alpha_calcolato_con_successo="<<alpha_calcolato_con_successo<<endl;
	//cout << "beta_calcolato_con_successo="<<beta_calcolato_con_successo<<endl;
	//cout << "do_optics="<<do_optics<<endl;
	//cout << "confronto_pos="<<confronto_pos<<endl;

	if (alpha_calcolato_con_successo&&beta_calcolato_con_successo&&(confronto_pos_x||confronto_pos_y))
	{
		if (gnuplot_xmax_opt > 0.) create_gnuplot_file( "Funzioni_Ottiche.plt", "Funzioni_Ottiche.txt", lunghezza, contatore, gnuplot_ymax_opt ,0.0, gnuplot_xmax_opt, etichette_ottiche,dati_rilevati,conto_per_confronto);
		else create_gnuplot_file( "Funzioni_Ottiche.plt", "Funzioni_Ottiche.txt", lunghezza, contatore, gnuplot_ymax_opt ,0.0, lunghezza_accumulata, etichette_ottiche,dati_rilevati,conto_per_confronto);
		create_gnuplot_file( "Parametri_Ellissi.plt", "Parametri_Ellissi_Funz_Ottiche.txt", lunghezza, contatore, gnuplot_ymax_ell ,0.0, lunghezza_accumulata, etichette_ellissi,dati_rilevati,conto_per_confronto);
		system ("gnuplot Funzioni_Ottiche.plt");
		system ("gnuplot Parametri_Ellissi.plt");

#ifdef TEST_OPTICAL_FUNCTIONS
		if (confronto_pos_t_x||confronto_pos_t_y)
		{
			create_gnuplot_file( "Funzioni_Ottiche_T.plt", "Funzioni_Ottiche_T.txt", lunghezza, contatore, gnuplot_ymax_opt ,0.0, lunghezza_accumulata, etichette_ottiche_T,dati_rilevati,conto_per_confronto);
			create_gnuplot_file( "Parametri_Ellissi_T.plt", "Parametri_Ellissi_Funz_Ottiche_T.txt", lunghezza, contatore, gnuplot_ymax_opt ,0.0, lunghezza_accumulata, etichette_ottiche_T,dati_rilevati,conto_per_confronto);
			system ("gnuplot Funzioni_Ottiche_T.plt");
			system ("gnuplot Parametri_Ellissi_T.plt");
		}
#endif
	}
	else	
	{
#if defined (__linux)
//		system ("rm Funzioni_Ottiche.txt"); 
#ifdef TEST_OPTICAL_FUNCTIONS
//		system ("rm Funzioni_Ottiche_T.txt");
#endif
#elif defined (_WIN32) || defined (_WIN64)
//		system ("del Funzioni_Ottiche.txt"); 
#ifdef TEST_OPTICAL_FUNCTIONS
//		system ("del Funzioni_Ottiche_T.txt");
#endif
#endif
	} 

	//cout << "Confronto: "<<confronto_pos<<endl;
	if ((confronto_pos_x==false)&&(confronto_pos_y==false))
	{
#if defined (__linux)
		system ("rm Math_rilevati.txt"); 
#elif defined (_WIN32) || defined (_WIN64)
		system ("del Math_rilevati.txt"); 
#endif
	}

#ifdef TEST_OPTICAL_FUNCTIONS
//	if ((confronto_pos_t_x==false)&&(confronto_pos_t_y==false))
//	{
//#if defined (__linux)
//		system ("rm Math_rilevati_T.txt"); 
//#elif defined (_WIN32) || defined (_WIN64)
//		system ("del Math_rilevati_T.txt"); 
//#endif
//	}
#endif

#ifdef DEBUG
	fclose(outputDEBUG);
#endif

}

