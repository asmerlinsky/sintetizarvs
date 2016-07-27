
#include <stdio.h>
#include <stdlib.h>
//#include <unistd.h>
//#include <errno.h>
#include <string.h>
//#include <math.h>


int main(int argc, char *argv[]){

  char *nombre;
  char salida[50];


  //Tomo las variables del input
  nombre = argv[1];

  /// en el caso que le quieras entrar un entero//entero = atoi(argv[1]);


  sprintf(salida,"zaraza--%s.dat", nombre);


  printf("\nel nuevo nombre es %s \n", salida);


  return 1;
}


