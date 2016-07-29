/*Sintetiza canto
HAY QUE INTRODUCIRLO POR LINEA DE COMANDO LOS ARCHIVOS

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

struct Par2 {
    double PinH; double PinHpunto;
} cc;

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


void helmholtz(int n, double v[], double dv[], double t){
  double a=-540E6, b = -7800, c = 1.8E8, d = 1.2E-2, e = 0.72, f = -0.83E-2, g = -5E2, h = 1E-4;
  double w, x, y, z;
  
  x=v[0];
  y=v[1];
  z=v[2];

  dv[0] = y;
  //dy/dt = -presion*gamma^2 - vs(t)*gamma^2*x   - gamma^2*x^3     - gamma*x^2*y + gamma^2*x^2  - gamma*x*y
  dv[1] = a*x+b*y + c*z + d*cc.PinHpunto + e*cc.PinH;
  dv[2] = f*y + g*z + h*cc.PinH;

  return;
}



int nearpow2up(int number){
    int i=0, temp=1;
    while(temp<number) {temp*=2; i++;}
        return(temp);
}


int main(int argc, char *argv[]) {
    int i,k, Ndatos, golord;
    char *nomfile;
    char entrada[100];
    nomfile=argv[1];
    char canto[100];
    char *entcanto;
    entcanto=argv[2];
    sprintf(canto,"%s",entcanto);
    sprintf(entrada,"%s",nomfile);
 
    int Ndatos1, Ndatos2;
    int POT2up;
    Ndatos1=filesize(entrada,1);
    Ndatos2=filesize(canto,1);
 
    int POT2up1=nearpow2up(Ndatos1);
    int POT2up2=nearpow2up(Ndatos2);
    if (POT2up1>=POT2up2){
        POT2up=POT2up1;
    }
    else {
        POT2up=POT2up2;
    }    
 
 
    golord=4;
    
    //CARGA TEMPLATE
    double *temp;


    
   

    temp=dvector(1,POT2up);
    file_to_vector(entrada,temp,1,Ndatos1,1,1);
    printf("\n\ntemplado es %s \t", nomfile);
    printf("Ndatos=%d\n",Ndatos1);
    for(i=Ndatos1;i<=POT2up;i++){temp[i]=0.0;}
    
    
    
    //CARGA PROXY PRESION
    double *sonido;
 
    sonido=dvector(1,POT2up);
    file_to_vector(canto,sonido,1,Ndatos2,1,1);
    printf("canto es %s \t", canto);
    printf("Ndatos=%d\n",Ndatos);
    for(i=Ndatos2;i<=POT2up;i++){sonido[i]=0.0;}    
    

    
    

   //CALCULA ENVOLVENTE TEMPLATE CON HILBERT
    double *hilb;
    hilb=dvector(1,POT2up);
   
    for(i=1;i<=500;i++) hilb[i]=0.0; //guarda si cambia LFILT
   
    hilbert(temp,hilb,POT2up);


    //SUAVIZA ENVOLVENTE TEMPLATE CON INTEGRACION
    double v[1],dt, t;
    double *av_sound;


    av_sound=dvector(1,POT2up);
    v[0]=0.0;
    k=0;
    dt=1/44150.;
    
    aa.tau=15./1500.;

    for(i=1;i<=POT2up;i++){
        
        aa.beta=hilb[i];
        rk4(takens,v,1,t+0.0,dt);
        av_sound[i]=v[0];
    }
    

    //SUAVIZA CANTO CON INTEGRACION
    dt=1/10000.;
    double *av_pres,*buff;
    char envcanto[100];
    t=0.0;
    buff=dvector(1,POT2up);
    av_pres=dvector(1,POT2up);
    v[0]=0.0;
    k=0;
    aa.tau=30/1500.;

    for(i=1;i<=POT2up;i++){
        aa.beta=sonido[i];
        rk4(takens,v,1,t+0.0,dt);
        buff[i]=v[0];
    }  

    for(i=1;i<=POT2up;i++){
        aa.beta=buff[i];
        rk4(takens,v,1,t+0.0,dt);
        av_pres[i]=v[0];
    }    
    
    strcpy(envcanto, "envolvente.");
    strcat(envcanto, entcanto);
    strcat(envcanto, ".dat");
    vector_to_file(envcanto,av_pres,1,POT2up);


    
    
    //FILTRO CON SAV-GOL
    
    int np,nl,nr,ld,m,index,Nmin;
    char intoutput[100],envolvente[100];
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
    printf("POT2up=%d\n",POT2up);
    //busco maximo de la envolvente vs
    for(i=1;i<=POT2up;i++){
        if (envmax<fabs(sav[i])){
            envmax=fabs(sav[i]);
           
        }
        
    }

    printf("envmax=%g\n",envmax);
    //defino renormalización de la envolvente
    double umbral;
    double rango;
    double min=0.0000001;
    rango=0.1-min;
    umbral=0.000017;//Defino un umbral de presión para prendido y apagado
    //umbral+=0.002;
    double renorm;
    renorm=(rango / envmax);
    for(i=1;i<=POT2up;i++){//renormalizo vs entre 0.002 y 2.99
        sav[i]*=renorm;
        sav[i]+=min;
        // printf("sav[%d]=%g\n",i,sav[i]);
    }
    
    
    
    bb.g=24000;
    int pasosint=40;
    dt=1./(44150.*pasosint);
    double result[2];
    int taux=1;
    result[0]=0.000000001;
    result[1]=0.000000001;
    printf("umbral=%g\n",umbral);
    i=0;
    
    if (av_pres[0]>umbral){bb.pres=0.15;}
    else {bb.pres=-0.15;}       
    bb.VS=sav[0];    
    
    
    
    while (i<=POT2up) {
        if (taux==pasosint){
            sintetico[i]=result[1];
            if (av_pres[i]>umbral){ bb.pres=0.15; } //si la presion supera el umbral prendo el sistema
            else { bb.pres=-0.15; }
                        
            bb.VS=sav[i]-sav[i]*sav[i]*(9E-5);
            taux=0;
            i++;            
        }       
        rk4(modelo,result,2,t,dt); //sintetizo
        taux++;
        t+=dt;
    }
    //guardo el sintetizado y la envolvente de vs
    strcpy(intoutput, "sintetizado.");
    strcat(intoutput, nomfile);
    strcat(intoutput, ".dat");
    vector_to_file(intoutput,sintetico,1,POT2up);
    
    strcpy(envolvente, "envolvente.");
    strcat(envolvente, nomfile);
    strcat(envolvente, ".dat");
    vector_to_file(envolvente,sav,1,POT2up);
    
    //paso por el filtro traqueal
    double *ptraqueal, *ptraqout;
    double ref=0.5;
    ptraqueal=dvector(1,POT2up);
    ptraqout=dvector(1,POT2up);
    
    for (i=0;i<=8;i++){
        ptraqueal[i]=sintetico[i]*av_pres[i];
    }
    for (i=9;i<=POT2up;i++){
        ptraqueal[i]=sintetico[i]*av_pres[i]-ref*(ptraqueal[i-9]);
    }

    for(i=0;i<=4;i++){ ptraqout[i]=0.0;}
    for (i=5;i<=POT2up;i++){
    ptraqout[i]=(1-ref)*ptraqueal[i-5];
    }
    //guardo las señales
    char nomptraqueal[100];
    strcpy(nomptraqueal, "ptraqueal.");
    strcat(nomptraqueal, nomfile);
    strcat(nomptraqueal, ".dat");
    vector_to_file(nomptraqueal,ptraqueal,1,POT2up);


    char nomptraqout[100];
    strcpy(nomptraqout, "ptraqout.");
    strcat(nomptraqout, nomfile);
    strcat(nomptraqout, ".dat");
    vector_to_file(nomptraqout,ptraqout,1,POT2up);
    
    //integro el resonador de helmholtz
    double vv[3]; 
    double *Pfinal;
    Pfinal=dvector(1,POT2up);
    vv[0]=0.000000000001;
    vv[1]=0.000000000001;
    vv[2]=0.000000000001;
    // dt=1./(44150.*pasosint);
    t=0.0;
    Pfinal[0]=0.0;
    i=0;
    taux=0;
    while(i<=POT2up){
    
        if (taux==pasosint){
            Pfinal[i]=vv[2];
            cc.PinHpunto=(ptraqout[i]-ptraqout[i-1])/dt;
            cc.PinH=ptraqout[i];
            taux=0;
            i++;
        }
        rk4(helmholtz,vv,3,t,dt);
        taux++;
        t+=dt;
        
        
        
    }    
    
    Pfinal[POT2up]=0.0;
    
    char salidafinal[100];
    strcpy(salidafinal, "pfinal.");
    strcat(salidafinal, nomfile);
    strcat(salidafinal, ".dat");
    vector_to_file(salidafinal,Pfinal,1,POT2up);    
    
    
    free_dvector(ptraqueal,1,POT2up);
    free_dvector(ptraqout,1,POT2up);
    free_dvector(sintetico,1,POT2up);
    free_dvector(av_sound,1,POT2up);
    free_dvector(hilb,1,POT2up);
    free_dvector(sav,1,POT2up);
    printf("listo\n\n");
    return 0;

}


