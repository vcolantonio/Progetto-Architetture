/**************************************************************************************
* 
* CdL Magistrale in Ingegneria Informatica
* Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2020/21
* 
* Progetto dell'algoritmo di Regressione
* in linguaggio assembly x86-32 + SSE
* 
* Fabrizio Angiulli, aprile 2019
* 
**************************************************************************************/

/*
* 
* Software necessario per l'esecuzione:
* 
*    NASM (www.nasm.us)
*    GCC (gcc.gnu.org)
* 
* entrambi sono disponibili come pacchetti software 
* installabili mediante il packaging tool del sistema 
* operativo; per esempio, su Ubuntu, mediante i comandi:
* 
*    sudo apt-get install nasm
*    sudo apt-get install gcc
* 
* potrebbe essere necessario installare le seguenti librerie:
* 
*    sudo apt-get install lib32gcc-4.8-dev (o altra versione)
*    sudo apt-get install libc6-dev-i386
* 
* Per generare il file eseguibile:
* 
* nasm -f elf64 regression64.nasm && gcc -m64 -msse -O0 -no-pie sseutils64.o regression64.o regression46c.c -o regression64c -lm && ./regression64c $pars
* 
* oppure
* 
* ./runregression64
* 
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

#define	type		double
#define	MATRIX		type*
#define	VECTOR		type*

typedef struct {
    MATRIX x; //data set
    VECTOR y; //label set
    MATRIX xast; //data set convertito
    int n; //numero di punti del data set
    int d; //numero di dimensioni del data set    
    int k; //dimensione del batch
    int degree; //grado del polinomio
    type eta; //learning rate
    //STRUTTURE OUTPUT
    VECTOR theta; //vettore dei parametri
    int t; //numero di parametri, dimensione del vettore theta
    int iter; //numero di iterazioni
    int adagrad; //accelerazione adagrad
    int silent; //silenzioso
    int display; //stampa risultati
} params;

/*
* 
*	Le funzioni sono state scritte assumento che le matrici siano memorizzate 
* 	mediante un array (double*), in modo da occupare un unico blocco
* 	di memoria, ma a scelta del candidato possono essere 
* 	memorizzate mediante array di array (double**).
* 
* 	In entrambi i casi il candidato dovrà inoltre scegliere se memorizzare le
* 	matrici per righe (row-major order) o per colonne (column major-order).
*
* 	L'assunzione corrente è che le matrici siano in row-major order.
* 
*/

void* get_block(int size, int elements) { 
    return _mm_malloc(elements*size,32); 
}

void free_block(void* p) { 
    _mm_free(p);
}

MATRIX alloc_matrix(int rows, int cols) {
    return (MATRIX) get_block(sizeof(type),rows*cols);
}

void dealloc_matrix(MATRIX mat) {
    free_block(mat);
}

extern void popola( int tuple,int coldest, int h,  int n, params* input, int* j);

void trasposta(MATRIX data, int rows, int cols){
  MATRIX temp = alloc_matrix(cols,rows);
  for(int j=0; j<cols; j++)
    for(int i=0; i<rows; i++)
      temp[j*rows+i] = data[i*cols+j];
  
  for(int j=0; j<cols; j++)
    for(int i=0; i<rows; i++)
      data[j*rows+i] = temp[j*rows+i];
  dealloc_matrix(temp);
  
}


void trasposta1(MATRIX data, int rows, int cols){
  MATRIX temp = alloc_matrix(cols,rows);
  for(int j=0; j<cols; j++)
    for(int i=0; i<rows; i++)
      temp[j+i*cols] = data[i+j*rows];
  
  for(int j=0; j<cols; j++)
    for(int i=0; i<rows; i++)
      data[j+i*cols] = temp[j+i*cols];
  dealloc_matrix(temp);
}

/*
* 
* 	load_data
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*8 byte: matrix data in row-major order --> numeri floating-point a precisione doppia
* 
*****************************************************************************
*	Se lo si ritiene opportuno, è possibile cambiare la codifica in memoria
* 	della matrice. 
*****************************************************************************
* 
*/
MATRIX load_data(char* filename, int *n, int *k) {
    FILE* fp;
    int rows, cols, status, i;
    
    fp = fopen(filename, "rb");
    
    if (fp == NULL){
        printf("'%s': bad data file name!\n", filename);
        exit(0);
    }
    
    status = fread(&cols, sizeof(int), 1, fp);
    status = fread(&rows, sizeof(int), 1, fp);
    
    MATRIX data = alloc_matrix(rows,cols);
    status = fread(data, sizeof(type), rows*cols, fp);
    fclose(fp);

    trasposta(data, rows, cols);
    
    *n = rows;
    *k = cols;
    
    return data;
}

/*
* 	save_data
* 	=========
* 
*	Salva su file un array lineare in row-major order
*	come matrice di N righe e M colonne
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*8 byte: matrix data in row-major order --> numeri floating-point a precisione doppia
*/
void save_data(char* filename, void* X, int n, int k) {
    FILE* fp;
    int i;
    fp = fopen(filename, "wb");
    if(X != NULL){
        fwrite(&k, 4, 1, fp);
        fwrite(&n, 4, 1, fp);
        for (i = 0; i < n; i++) {
            fwrite(X, sizeof(type), k, fp);
            //printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
            X += sizeof(type)*k;
        }
    }
    fclose(fp);
}



void combinazioni(int *v, int n, int k){

        int num= n+k-1;
        int den=n-1;
        int div=k;

        if(den<div){
          div= den;
          den=k;
        } 

        int parziale=1;
        int fact=1;

        for(int i=num;i>den;i--){
            parziale*=i;
            fact*=div; 
            div--;
        }

        v[k]=parziale/fact;
}


void calcolaJ(int* J, int d, int h){

    for(int k = 0; k<h;k++) J[k]=1;

    int k=h-1;
    int indice = 0;
    while(k>=0){
        int v = J[indice*h + k] + 1;
        if(v>d) k--;
        else{
            //copio <---
            for(int t=0; t<k; t++) J[(indice+1)*h + t] = J[indice*h + t];
            for(int t=k; t<h; t++) J[(indice+1)*h + t] = v;
            k=h-1; indice++;
        }
    }
}

/*
void scrivi_colonne(int tuple, int n, int coldest, int h, int* j, params* input){
  type num;
  int i, oss, k;
  for(i=0;i<tuple;i++){

      for(oss=0;oss<n;oss++){
        num=1;
        for(k=0;k<h;k++){
          num=num* input->x[oss+(j[i*h+k]-1)*n];
          
        }
        input->xast[oss+coldest*n]=num;
      }

      coldest++; //andiamo a scrivere la colonna coldest successiva di xast
  }
}
*/
void convert_data(params* input){

    int d,n,degree,h,oss;
    int coldest; //colonna del dataset di dest in cui inserire il valore calcolato
    type num; //da supporto per calcolare il valore risultato del dataset
    int i,k;
    int* v;
  //  int* j;
    int* vettore_combinazioni;
	int* vettore_col;
    int col=0; //numero colonne del dataset risultato
    d=input->d;
    n=input->n;
    degree= input->degree;


    vettore_combinazioni =(int*) malloc((degree+1)*sizeof(int));
	vettore_col = (int*) malloc((degree+1)*sizeof(int));

  	vettore_col[1]=1;
	vettore_col[0]=0;
	#pragma omp parralel for
    for(h=0;h<=degree;h++){
		combinazioni(vettore_combinazioni,d,h);
    }

    for(h=0;h<=degree;h++){
		
        col+=vettore_combinazioni[h];
		if(h>1){	
			vettore_col[h]=vettore_col[h-1]+vettore_combinazioni[h-1];
		}
		
		
    }

    input->t =col;

    input->xast = alloc_matrix(input->n,col);

	

	//h=3 d=2 degree=4
	// J
	// 1,1,1 -> colonna del dataset risultato x1^3
	// 1,1,2	-> x1^2*x2
	// 1,2,2	-> x1*x2^2
	// 2,2,2	-> x2^3

	//popolamento della prima colonna del dataset convertito
	//Inserisco tutti 1
	
	for(i=0;i<n;i++){
		input->xast[i]=1;
	}

	


	#pragma omp parallel for	
    for(h=1;h<=degree;h++){

    	int tuple=vettore_combinazioni[h]; //corrisponderà al numero di tuple in J
    	int* j= (int*) malloc (tuple*h*sizeof(int)); //h= dimensione della singola tupla

    	//h=3 d=2 tuple=4
    	// -> colonna del dataset risultato x1^3
		// 1,1,2	-> x1^2*x2
		// 1,2,2	-> x1*x2^2
		// 2,2,2	-> x2^3
    	// 1 ---- x1^2*x2, x1*x2^2, x2^3

    	calcolaJ(j,d,h); //il risultato sovrascrive J

		popola(tuple, vettore_col[h], h, n, input, j);
		free(j);
		
		
    }
	
	

	input-> theta= (VECTOR) malloc(col*sizeof(type));
	
}//convert_data

extern void prodottoScalare(type* cost,VECTOR theta,MATRIX xast,int j, int size);
/*
void prodottoScalare(type* cost,VECTOR theta,MATRIX xast,int j, int size){

	type ris=0;

	for(int i=0; i<size; i++){

		ris+= theta[i]* xast[j+ i];

	}

	*cost=ris;
}
*/

extern void sommaVettori(VECTOR daSottrarre, VECTOR temp, int size);
/*
void sommaVettori(VECTOR daSottrarre, VECTOR temp, int size){

	for(int i=0; i<size; i++){
		daSottrarre[i] += temp[i];
	}
}
*/

extern void moltiplicazionePerScalare(VECTOR vett_ris, MATRIX x, int i,int size, type scalare);
/*
void moltiplicazionePerScalare(VECTOR vett_ris, MATRIX x, int i, int size, type scalare){

	for(int k=0; k<size; k++){
		vett_ris[k]= scalare * x[i + k];
	}

}
*/
/*
void moltiplicazionePerScalareADAGRAD(VECTOR g, MATRIX x, int i, int size, type scalare){

	for(int k=0; k<size; k++){
		g[k]= scalare * x[i + k];
	}

}
*/
extern void quadratoVettore(VECTOR g, int riga, MATRIX G, int size);
/*
void quadratoVettore(VECTOR g, int riga, MATRIX G, int size){
  for(int i =0; i<size; i++){
    type n = g[i];
    G[riga + i] +=  n*n;
  }
}
*/
extern void faiSommatoria(VECTOR daSottrarre, int size, VECTOR g, MATRIX G, int riga, type epsilon);
/*
void faiSommatoria(VECTOR daSottrarre, int size, VECTOR g, MATRIX G, int riga, type epsilon){
  type ng, nG;

  for(int h =0; h<size; h++) {
    nG = G[riga + h];
    ng = g[h];
    daSottrarre[h] += ng / ( sqrt(nG)+epsilon )  ;
  }
}
*/

extern void calcolaTheta(VECTOR daSottrarre, VECTOR theta, type eta, int size);

extern void azzeramentoDaSottrarre(VECTOR daSottrarre, int size);

void sgd(params* input){
    // -------------------------------------------------
    // Codificare qui l'algoritmo sgd
    // -------------------------------------------------

    type eta;
    int n,iter,adagrad,size;
    type cost_molt;
	MATRIX xast;
    VECTOR y;
	VECTOR theta;
	VECTOR daSottrarre;


    int i,j,k,it,v;

    eta= input->eta;
    n= input->n;
    iter= input->iter;
    k=input->k;
    adagrad=input->adagrad;
    xast=input->xast;
    y= input->y;
    theta= input->theta;
	size = input->t;



    daSottrarre= (VECTOR) malloc(size*sizeof(type));


    for(i=0;i<size;i++) {
    	theta[i]=0; //inizializzo theta
	   }
	it = 0;
  if(adagrad){
    type epsilon = 1E-8;
    int riga = 0;
    VECTOR g=(VECTOR) malloc(size*sizeof(type));
    MATRIX G =(MATRIX) malloc(size*k*sizeof(type));

    for(int m = 0; m<k*size; m++)
        G[m]=0;
      

    while(it<iter){
      for(i=0; i<n;){
        //calcolo del minimo
          if((n-i)<k) v=n-i;
          else v=k;

          //inizializzo vettore da sottrarre
	      azzeramentoDaSottrarre(daSottrarre,size);
	 
         // for(int index=0;index<size;index++) daSottrarre[index]=0;


          for(j=i; j<i+v; j++){
            riga = j % k ;

            //calcolo riga gj
            prodottoScalare(&cost_molt, theta, xast, j*size, size);
            cost_molt -=  y[j];
            moltiplicazionePerScalare(g, xast, j*size, size, cost_molt);

            // Gj += gj^2
            quadratoVettore(g, riga*size, G, size);

            faiSommatoria(daSottrarre, size, g, G, riga*size, epsilon);


          }

	  calcolaTheta(daSottrarre, theta, eta/v, size);
	  /*
          for(int h =0; h<size; h++) {
            daSottrarre[h] *= eta/v;
            theta[h] -= daSottrarre[h];
          }
	  */

          i=i+v;

        }
      it++;
    }
    free(g); free(G);
  }
  else{
    VECTOR temp= (VECTOR) malloc(size*sizeof(type));
    while(it<iter){
    	for(i=0; i<n; ){

 			//calcolo del minimo
    		if((n-i)<k) v=n-i;
    		else v=k;

    		//inizializzo vettore da sottrarre
		azzeramentoDaSottrarre(daSottrarre,size);
		
		/*
    		for(int index=0;index<size;index++){
          		daSottrarre[index]=0; 
        	}
		*/

    		for(j=i; j<i+v; j++){

    			//calcolo di temp

    		  	prodottoScalare(&cost_molt, theta, xast, j*size, size);


    			cost_molt= cost_molt - y[j];


    			moltiplicazionePerScalare(temp, xast, j*size, size, cost_molt);

    			sommaVettori(daSottrarre, temp, size);
    		}
    		

		 calcolaTheta(daSottrarre, theta, eta/v, size);
		/*
    		for(int indice=0;indice<size;indice++){
    			daSottrarre[indice]= daSottrarre[indice]*(div);
    			theta[indice]= theta[indice]-daSottrarre[indice];
    		}
		*/

        i = i+v;
    	}

    	it++;
    }
   // free(temp);
  }
 // free(xast); free(y); free(theta); free(daSottrarre);
	

}//sgd

int main(int argc, char** argv) {

    char fname[256];
    char* dsname;
    char* filename;
    int i, j, k;
    clock_t t;
    float time;
    int yd = 1;
    
    //
    // Imposta i valori di default dei parametri
    //

    params* input = malloc(sizeof(params));
    
    input->x = NULL;
    input->y = NULL;
    input->xast = NULL;
    input->n = 0;
    input->d = 0;
    input->k = -1;
    input->degree = -1;
    input->eta = -1;
    input->iter = -1;
    input->adagrad = 0;
    input->theta = NULL;
    input->t = 0;
    input->adagrad = 0;
    input->silent = 0;
    input->display = 0;

    //
    // Visualizza la sintassi del passaggio dei parametri da riga comandi
    //

    if(argc <= 1){
        printf("%s D -batch <k> -degree <deg> -eta <eta> -iter <it> [-adagrad]\n", argv[0]);
        printf("\nParameters:\n");
        printf("\tD: il nome del file, estensione .data per i dati x, estensione .labels per le etichette y\n");
        printf("\t-batch <k>: il numero di campini nel batch\n");
        printf("\t-degree <deg>: il grado del polinomio\n");
        printf("\t-eta <eta>: il learning rate\n");
        printf("\t-iter <it>: il numero di iterazioni\n");
        printf("\t-adagrad: l'acceleratore AdaGrad\n");
        exit(0);
    }
    
    //
    // Legge i valori dei parametri da riga comandi
    //
    
    int par = 1;
    while (par < argc) {
        if (par == 1) {
            filename = argv[par];
            par++;
        } else if (strcmp(argv[par],"-s") == 0) {
            input->silent = 1;
            par++;
        } else if (strcmp(argv[par],"-d") == 0) {
            input->display = 1;
            par++;
        } else if (strcmp(argv[par],"-batch") == 0) {
            par++;
            if (par >= argc) {
                printf("Missing batch dimension value!\n");
                exit(1);
            }
            input->k = atoi(argv[par]);
            par++;
        } else if (strcmp(argv[par],"-degree") == 0) {
            par++;
            if (par >= argc) {
                printf("Missing degree value!\n");
                exit(1);
            }
            input->degree = atoi(argv[par]);
            par++;
        } else if (strcmp(argv[par],"-eta") == 0) {
            par++;
            if (par >= argc) {
                printf("Missing eta value!\n");
                exit(1);
            }
            input->eta = atof(argv[par]);
            par++;
        } else if (strcmp(argv[par],"-iter") == 0) {
            par++;
            if (par >= argc) {
                printf("Missing iter value!\n");
                exit(1);
            }
            input->iter = atoi(argv[par]);
            par++;
        } else if (strcmp(argv[par],"-adagrad") == 0) {
            input->adagrad = 1;
            par++;
        } else{
            printf("WARNING: unrecognized parameter '%s'!\n",argv[par]);
            par++;
        }
    }
    
    //
    // Legge i dati e verifica la correttezza dei parametri
    //
    
    if(filename == NULL || strlen(filename) == 0){
        printf("Missing input file name!\n");
        exit(1);
    }

    dsname = basename(strdup(filename));
    sprintf(fname, "%s.data", filename);
    input->x = load_data(fname, &input->n, &input->d);
    sprintf(fname, "%s.labels", filename);
    input->y = load_data(fname, &input->n, &yd);

    if(input->k < 0){
        printf("Invalid value of batch dimension parameter!\n");
        exit(1);
    }
    
    if(input->degree < 0){
        printf("Invalid value of degree parameter!\n");
        exit(1);
    }
    
    if(input->eta < 0){
        printf("Invalid value of eta parameter!\n");
        exit(1);
    }
    
    if(input->iter < 0){
        printf("Invalid value of iter parameter!\n");
        exit(1);
    }
    
    //
    // Visualizza il valore dei parametri
    //
    
    if(!input->silent){
        printf("Input data name: '%s.data'\n", filename);
        printf("Input label name: '%s.labels'\n", filename);
        printf("Data set size [n]: %d\n", input->n);
        printf("Number of dimensions [d]: %d\n", input->d);
        printf("Batch dimension: %d\n", input->k);
        printf("Degree: %d\n", input->degree);
        printf("Eta: %f\n", input->eta);
        if(input->adagrad)
            printf("Adagrad enabled\n");
        else
            printf("Adagrad disabled\n");
    }
    
  
    //
    // Conversione Dati
    //
    
    t = clock();
    convert_data(input);
    t = clock() - t;
    time = ((float)t)/CLOCKS_PER_SEC;
    sprintf(fname, "%s.xast", dsname);
       
    if(!input->silent)
        printf("Conversion time = %.3f secs\n", time);
    else
        printf("%.3f\n", time);
    
    //
    // Regressione
    //

    trasposta1(input->xast, input->n, input->t);
    
    t = clock();
    sgd(input);
    t = clock() - t;
    time = ((float)t)/CLOCKS_PER_SEC;
    
    if(!input->silent)
        printf("Regression time = %.3f secs\n", time);
    else
        printf("%.3f\n", time);

    //
    // Salva il risultato di theta
    //
    if(!input->adagrad)
	    sprintf(fname, "%s.theta.sgdomp", dsname);
    else
	    sprintf(fname, "%s.theta.adagradomp", dsname);
    save_data(fname, input->theta, input->t, 1);
    if(input->display){
        printf("theta: [");
        for(i=0; i<input->t-1; i++)
            printf("%f,", input->theta[i]);
        printf("%f]\n", input->theta[i]);
    }
    
    if(!input->silent)
        printf("\nDone.\n");

    return 0;
}
