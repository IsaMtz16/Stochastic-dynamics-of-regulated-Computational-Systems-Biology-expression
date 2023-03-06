#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

//Definimos el número max de simulaciones
#define N 1000



//Generador de números aleatorios:
//PARISI-RAPUANO
extern float Random(void);
#define NormRANu (2.3283063671E-10)
extern void ini_ran(int SEMILLA);
unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran, ig1, ig2, ig3;


float Random(void)    // generar numero aleatorio
{
    float r;
    ig1=ind_ran-24;
    ig2=ind_ran-55;
    ig3=ind_ran-61;
    irr[ind_ran]=irr[ig1]+irr[ig2];
    ir1=(irr[ind_ran]^irr[ig3]);
    ind_ran++;
    r=ir1*NormRANu;
    return r;
}

void ini_ran(int SEMILLA){
    int INI,FACTOR,SUM,i;
    srand(SEMILLA);
    INI=SEMILLA;
    FACTOR=67397;
    SUM=7364893;
    for(i=0; i<256; i++){
        INI=(INI*FACTOR+SUM);
        irr[i]=INI;
    }
    ind_ran=ig1=ig2=ig3=0;
}



int main(){
    double K=602;
    ini_ran(time(NULL));
    double alpham=10100;
    double deltam=1;
    double alphap=10;
    double deltap=0.1;

    double V=0.6022;

    //Estos coeficientes nos dan probabilidades
    double a1; //generación de RNAm
    double a2; // degeneración de RNAm DEPENDE DEL NUMERO DE RNAM
    double a3; // generacion de proteinas
    double a4; // degeneración de proteínas
    double a0; //sumamos todas las probabilidades

    //Definimos la variable de tiempo
    double t;
    //El tiempo final
    double t_final=20;



   //DEFINIMOS EL NÚMERO DE MÓLECULAS
    double M=1000*V; //RNAm
    double P=10000*V;// Proteinas

    //NUMEROS ALEATORIOS
    double r1;
    double r2;

    //MEDIAS
    double media_M;
    double media_P;

    //PASO DE TIEMPO
    double tau;

    double aux1;
    double aux2;

    //ABRIMOS EL FICHERO
    FILE *f;
    f=fopen("ARN.txt", "w");

    FILE *f1;
    f1=fopen("Proteinas.txt", "w");

    int i;

    //HACEMOS UN BUCLE PARA EL NUMERO DE SIMULACIONES
    int n_iter;
    for(n_iter=0;n_iter<N;n_iter++){

        M=100*V;
        P=10000*V;
        t=0;

        //printf("%d\n",n_iter);
        while(t<t_final){
            a1=V*alpham*K*K/(K*K+P*P);//GENERACION RNAm
            a2=deltam*M; //DEGRADACION RNAm
            a3=alphap*M; //GENERACION PROTEINAS
            a4=deltap*P;//DEGRADACION PROTEINAS

            a0=a1+a2+a3+a4;

            r1=Random();


            tau=-log(r1)/a0;
            t=t+tau;

            r2=a0*Random();

            if(r2<=(a1)){
                M=M+1;
            }
            else{
                if((a1<r2)&&(r2<=(a1+a2))){
                    M=M-1;
                }
                else{
                        if(((a1+a2)<r2)&&(r2<=(a1+a2+a3))){
                            P=P+1;
                        }
                        else{
                            if(((a1+a2+a3)<r2)&&(r2<=(a1+a2+a3+a4))){
                               P=P-1;
                            }
                        }
                }
            }

            i++;
        }
        fprintf(f,"%lf      \n",M/V);
        fprintf(f1,"%lf     \n",P/V);

    }
    fclose(f);
    fclose(f1);

    return 0;
}

