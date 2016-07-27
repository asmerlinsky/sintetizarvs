
/*Realiza envolvente de templados
HAY QUE INTRODUCIRLO POR LINEA DE COMANDO COMO ARGV

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "librerias.c"
#include "nrutil.h"
#include "nrutil.c"
#include "ht.c"
#include "rk4.c"
#include <string.h>

#include"convlv.c"
#include"savgol.c"
#include"realft.c"
#include"twofft.c"
#include"four1.c"
#include"lubksb.c"
#include"ludcmp.c"
#include <stdbool.h>



struct Par {
    double beta; double tau;
} aa;


struct Par1 {
    double pres; double g; double VS;
} bb;



void takens(int n,double v[],double dv[],double t) {
    double x;    
    x=v[0];     
    dv[0]=-(1/aa.tau)*x+fabs(aa.beta);
    return;
}



// Ecuaciones del sistema
void modelo(int n, double v[], double dv[], double t){
  double x,y,z;
  x=v[0];
  y=v[1];
  

  dv[0] = y;
  //dy/dt = -presion*gamma^2 - vs(t)*gamma^2*x   - gamma^2*x^3     - gamma*x^2*y + gamma^2*x^2  - gamma*x*y
  dv[1] = -bb.pres*bb.g*bb.g - bb.VS*bb.g*bb.g*x - bb.g*bb.g*x*x*x - bb.g*x*x*y + bb.g*bb.g*x*x - bb.g*x*y;
  

  return;
}




int nearpow2up(int number){
    int i=0, temp=1;
    while(temp<number) {temp*=2; i++;}
        return(temp);
}


int main(int argc, char *argv[]) {
    int i,k,Ndatos, golord;
    char *nomfile;
    char entrada[100];
    nomfile=argv[1];
    sprintf(entrada,"%s",nomfile);
    //printf("nomfile=%s\n",entrada);
    golord=4;
    
    //CARGA TEMPLATE
    double *temp;
    Ndatos=filesize(entrada,1);

    
    int POT2up=nearpow2up(Ndatos);

    temp=dvector(1,POT2up);
    file_to_vector(entrada,temp,1,Ndatos,1,1);
    printf("templado es %s \n", nomfile);
    printf("Ndatos=%d\n",Ndatos);
    for(i=Ndatos;i<=POT2up;i++){temp[i]=temp[Ndatos];}


   //CALCULA ENVOLVENTE CON HILBERT
    double *hilb;
    hilb=dvector(1,POT2up);
   
    for(i=1;i<=500;i++) hilb[i]=0.0; //guarda si cambia LFILT
   
    hilbert(temp,hilb,POT2up);



    //SUAVIZA ENVOLVENTE CON INTEGRACION
    double v[1],dt, t;
    double *av_sound;


    av_sound=dvector(1,POT2up);
    v[0]=0.0;
    k=0;
    dt=1/10000.;
    
    aa.tau=30/1500.;

    for(i=1;i<=POT2up;i++){
        aa.beta=hilb[i];
        rk4(takens,v,1,t+0.0,dt);
        av_sound[i]=v[0];
    }
    
    
    
    int np,nl,nr,ld,m,index,Nmin;
    char sonido[100],envolvente[100];
    float *c1,*data1,*ans1,dum1;
   
    c1=vector(1,POT2up); data1=vector(1,POT2up); ans1=vector(1,2*POT2up);

    double *sav;
    sav=dvector(1,POT2up);
  
    for(i=1;i<=POT2up;i++) data1[i]=(float) av_sound[i];
    
    savgol(c1,513,256,256,0,golord);
    
    for(index=1;index<=POT2up;index++) data1[index]=fabs(data1[index]);
    
    convlv(data1,POT2up,c1,513,1,ans1);

    for(i=2;i<POT2up-1;i++) {sav[i]=(double) ans1[i];
        
    }
    //LOS DATOS SON sav[i]
    

    double *sintetico;
    sintetico=dvector(1,POT2up);
    
    double envmax=0.0;
    
    //busco maximo*0.05 y ese es umbral de presion
    for(i=1;i<=Ndatos;i++){
        if (envmax<fabs(sav[i])){
            envmax=fabs(sav[i]);
        }

        
    }
    double umbral;
    double rango;
    rango=3.-0.002;
    umbral=0.1*rango;
    umbral+=0.002;
    double renorm;
    renorm=(rango / envmax);
    for(i=1;i<=POT2up;i++){//renormalizo vs entre 0.002 y 2.99
        sav[i]*=renorm;
        sav[i]+=0.002;
        // printf("sav[%d]=%g\n",i,sav[i]);
    }
    
    printf("envmax=%g\n",envmax);
    
    bb.g=10000;
    dt=1./44150.;
    double result[2];
    result[0]=0.0;
    result[1]=0.0;
    printf("umbral=%g\n",umbral);
    for(i=1;i<=POT2up;i++){
        if (sav[i]>umbral){
            bb.pres=0.15;
        }
        else {
        bb.pres=-0.15;
        }
        
        
        bb.VS=sav[i];
        rk4(modelo,result,2,t,dt);
        sintetico[i]=result[0];
    }
    
    strcpy(sonido, "sintetizado.");
    strcat(sonido, nomfile);
    strcat(sonido, ".dat");
    vector_to_file(sonido,sintetico,1,POT2up);
    
    strcpy(envolvente, "envolvente.");
    strcat(envolvente, nomfile);
    strcat(envolvente, ".dat");
    vector_to_file(envolvente,sav,1,POT2up);
    
    
            
    printf("tauintegracion: %g\t",aa.tau);
    printf("orden del filtro: %d\n\n\n",golord);

    
    
    
    free_dvector(sintetico,1,POT2up);
    free_dvector(av_sound,1,POT2up);
    free_dvector(hilb,1,POT2up);
    free_dvector(sav,1,POT2up);

}


