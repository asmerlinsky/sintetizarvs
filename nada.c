//ENCUENTRA GTES como MINIMOS LOCALES, MÁXIMOS GLOBALES, INICIOS Y FINALES DE SÍLABAS
//
//p [0:80000]'song.Bird_105.dat' u ($1/1) w l lt 2, 'norma.dat' u ($1) w l lt 1 lw 1, 'gtes1.Bird_105.dat' u 1:(.00) w p ps 1 lt 3, 'gtes2.Bird_105.dat' u 1:(.5) w p ps 1 lt 4, 'gtes3.Bird_105.dat' u 1:(.2) w p ps 1 lt 5, 'gtes4.Bird_105.dat' u 1:(.6) w p lt 7


//El parámetro que selecciona los mínimos locales está en la línea 191, después de dos líneas enteras comentadas


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "librerias.c"
#include "nrutil.h"
#include "nrutil.c"
#include "ht.c"
#include "rk4.c"

#include"convlv.c"
#include"savgol.c"
#include"realft.c"
#include"twofft.c"
#include"four1.c"
#include"lubksb.c"
#include"ludcmp.c"


struct Par {
    double beta; double tau;
} aa;


void takens(int n,double v[],double dv[],double t) {
    double x;    
    x=v[0];     
    dv[0]=-(1/aa.tau)*x+fabs(aa.beta);
    return;
}

int nearpow2(int number){
 	int i=0, temp=1;
	while(temp<number) {temp*=2; i++;}
        return(temp/2);
}



int main() {
   int i,j,k,Ndatos;

    //CARGA SONG
    // re sampleo para tomar 1 de cada 4 
    double *original;
    Ndatos=filesize("Bird_105.Sound.dat",1);
    original=dvector(1,Ndatos); 
    file_to_vector("Bird_105.Sound.dat",original,1,Ndatos,1,1);
    double *sound;
    sound=dvector(1,Ndatos);
    
    //for(i=1;i<=Ndatos/4;i++) sound[i]=original[4*i];       
    //vector_to_file("song.Bird_105.dat",sound,1,Ndatos/4);
    //Ndatos=Ndatos/4; 
    
    for(i=1;i<=Ndatos;i++) sound[i]=original[i];       
    vector_to_file("song.Bird_105.dat",sound,1,Ndatos);


    printf("\n\n\tNdatos %d\n\tPot de 2 más cercana %d\n\n",Ndatos,nearpow2(Ndatos));

   //CALCULA ENVOLVENTE CON HILBERT
   double *hilb;
   hilb=dvector(1,Ndatos);
   hilbert(sound,hilb,Ndatos); 
   vector_to_file("hilbert.Bird_105.dat",hilb,1,Ndatos);


   //SUAVIZA ENVOLVENTE CON INTEGRACION
   double v[1],dt, t;
   double *av_sound;
   av_sound=dvector(1,Ndatos);

   k=0;
   dt=1/10000.;
   aa.tau=20./1500.;
   //aa.tau=.5/1500.;
   for(i=0;i<Ndatos;i++){
        aa.beta=hilb[i];
        rk4(takens,v,1,t+0.0,dt);
        av_sound[i]=v[0];
   }
   vector_to_file("envelope.Bird_105.dat",av_sound,1,Ndatos);


   //SAVITSKY-GOLAY
   int np,nl,nr,ld,m,index;
   int POT=nearpow2(Ndatos);
   float c[POT+1],data[POT+1],ans[2*POT+1],dum;
   double *sav;
   sav=dvector(1,Ndatos);   

   //savgol(c,1025,128,128,0,4);	
   savgol(c,513,256,256,0,4);	
   for(i=1;i<=POT;i++) data[i]=(float) av_sound[i]; 
   //savgol(c,1025,128,128,0,4);		
   //for(i=1;i<=POT;i++) data[i]=(float) hilb[i]; 
   for(index=1;index<=POT;index++) data[index]=fabs(data[index]);

   convlv(data,POT,c,513,1,ans);

   for(i=2;i<POT-1;i++) sav[i]=(double) ans[i]; 
   vector_to_file("env_suave_new.Bird_105.dat",sav,1,Ndatos);
  	 
   //DERIVA Y SUAVIZA con SavGol
   double *deri;
   deri=dvector(1,Ndatos);
   for(i=3;i<POT-1;i++) deri[i]=-sav[i+2]+8.*sav[i+1]-8.*sav[i-1]+sav[i-2];

   //savgol(c,1025,128,128,0,4);
   savgol(c,513,256,256,0,4);		
   for(i=1;i<=POT;i++) data[i]=(float) deri[i]; 
   convlv(data,POT,c,514,1,ans);
   for(i=2;i<POT-1;i++) deri[i]=(double) ans[i]; 
   vector_to_file("deri.Bird_105.dat",deri,1,Ndatos);
   

   //GTEs
   int temp=0, count=0, *inic, *fin;
   double *gte1=dvector(1,1000);
   double *gte2=dvector(1,1000);
   double *gte3a=dvector(1,1000);
   double *gte3b=dvector(1,1000);
   double *gte3=dvector(1,1000);
   double *gte4=dvector(1,1000); 
   
    double *norm; 
   norm=dvector(1,Ndatos);
   normalize(sav, Ndatos, norm);
   inic=ivector(1,100);
   fin=ivector(1,100);
	
   //INICIO Y FIN de silabas        
   for(i=1;i<POT;i++){ 
	if(norm[i]<0.05 && norm[i+1]>0.05) {gte1[++temp]=(double) i; inic[++count]=i; /*printf("%d ",inic[count]);*/}
	if(norm[i]>0.05 && norm[i+1]<0.05) {gte1[++temp]=(double) i; fin[count]=i; /*printf("%d ",fin[count]);*/}
   }	

    temp=0;
   //MAXIMO DE LA SÍLABA
   for(i=1;i<=count;i++) gte2[++temp]=(double) index_maximo(norm,inic[i],fin[i]);



   //MINIMOS INTRASILABICOS
    int comp;        
    temp=0;
    for(i=51;i<POT-52;i++) {        
        for(j=1;j<=count;j++){
            if(i>=inic[j]&&i<=fin[j])
            if(deri[i]<0. && deri[i+1]>0.) 	
	    //if(norm[i]>norm[i-200] && norm[i]>norm[i+200]) gte3[++temp]=(double) i; 
            gte3a[++temp]=(double) i; 
	}
    }
    comp=temp;

   //MAXIMOS INTRASILABICOS       
    temp=0;
    for(i=51;i<POT-52;i++) {        
        for(j=1;j<=count;j++){
            if(i>=inic[j]&&i<=fin[j])
            if(deri[i]>0. && deri[i+1]<0.) 	
	    //if(norm[i]>norm[i-200] && norm[i]>norm[i+200]) gte3[++temp]=(double) i; 
            gte3b[++temp]=(double) i; 
	}
    }    
    if(temp>comp) comp=temp;

    //SELECCIONA MÍNIMOS
    int min_loc, max_ant, max_post;
    temp=0;
    for(i=2;i<=comp;i++) {
	min_loc = (int) gte3a[i];	//indice del minimo

	//busco max consecutivos
        max_ant=0; max_post=0;
	for(j=1;j<=comp;j++){
		if(gte3b[j]<gte3a[i] && gte3b[j+1]>gte3a[i]) {max_ant=(int) gte3b[j];max_post=(int) gte3b[j+1];}
	}
	if(max_ant && max_post) 
    	//----------------------------------------------------------------
    	//----------------------------------------------------------------
		// if(norm[min_loc]<0.8*norm[max_ant] || norm[min_loc]<0.8*norm[max_post]) gte3[++temp]=gte3a[i];
        if(norm[min_loc]<0.8*norm[max_ant]) gte3[++temp]=gte3a[i];

    }

    //ULTIMO MAXIMO
    temp=0;
    int mintemp;
    for(i=1;i<=count;i++){
        mintemp=inic[i];
        for(j=inic[i];j<=fin[i];j++){
            if((deri[j]>0)&&(deri[j+1]<0))mintemp=j;
        }
        gte4[++temp]=(double)mintemp;
    }

    



    vector_to_file("norma.dat",norm,1,Ndatos); 
    vector_to_file("gtes1.Bird_105.dat",gte1,1,1000);
    vector_to_file("gtes2.Bird_105.dat",gte2,1,1000);
    vector_to_file("gtes3a.Bird_105.dat",gte3a,1,1000);
    vector_to_file("gtes4.Bird_105.dat",gte4,1,1000);
    vector_to_file("gtes3b.Bird_105.dat",gte3b,1,1000);
    vector_to_file("gtes3.Bird_105.dat",gte3,1,1000);
    
    free_ivector(inic,1,100);	
    free_ivector(fin,1,100);	
    free_dvector(gte1,1,1000);
    free_dvector(gte2,1,1000);
    free_dvector(gte3a,1,1000);
    free_dvector(gte3b,1,1000);
    free_dvector(gte3,1,1000);
    free_dvector(gte4,1,1000);
   
    free_dvector(norm,1,Ndatos);	
    free_dvector(sound,1,Ndatos);
    free_dvector(av_sound,1,Ndatos);
    free_dvector(hilb,1,Ndatos);
    free_dvector(sav,1,Ndatos);
    free_dvector(deri,1,Ndatos);

}


