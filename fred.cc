#include "string.h"
#include "stdio.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include "stdlib.h"
#include <string.h>


void fred_(int *hea,double *f,double *f1,double *ff,double *ff1,double *fg,double *fg2,double *fg5)
//
// hea[2] - technicalq
// f[192][31] -photon function
// ff[192][31] -hadron function
// f1[192][31] -deriv. photon function
// ff1[192][31] -deriv. hadron function
// fg[3][192][31] -coeff for 3-fit
//  0 -ampl
//  1 -time
//  2 -ped
// fg2[2][24][31] -coeff for 3-fit
//  0 -amp
//  1 -ped
// fg5[5][192][31] -coeff for 3-fit
//  0 -ampl phot
//  1 -time phot
//  2 -ped
//  3 -ampl hadr
//  4 -time hadr
//
//
{
  static FILE *fr=NULL;

  int h,nsiz,nsiz1,j,i;
  char filn[10];
  int k1,k2,ii;
  //  if(*ifl == 0)
    {
      sscanf("arr.dat","%s",filn);

      fr= fopen(filn, "rb");
    }


  if(fr == NULL)printf(" fila net ");
  nsiz=4;

  //	  printf("%d %d %d %d \n",hea[0],hea[1],h,j);
  for(ii=0;ii<192;ii++)
    {
      nsiz=2;
      nsiz1=2;
      h = fread(hea, nsiz, nsiz1, fr);      
      j=hea[0];
      nsiz=8;
      nsiz1=(j)>>3;
      
      //      printf(" f %d %d\n",nsiz1,j);
      h = fread(f+ii*31, nsiz, nsiz1, fr);
      //      printf("%d %d %d %d \n",hea[0],hea[1],h,ii);
      nsiz=2;
      nsiz1=2;
      h = fread(hea, nsiz, nsiz1, fr);
    }
  for(ii=0;ii<192;ii++)
    {
      nsiz=2;
      nsiz1=2;
      h = fread(hea, nsiz, nsiz1, fr);
      j=hea[0];
      nsiz=8;
      nsiz1=(j)>>3;
      h = fread(f1+ii*31, nsiz, nsiz1, fr);
      nsiz=2;
      nsiz1=2;
      h = fread(hea, nsiz, nsiz1, fr);
    }

  for(ii=0;ii<192;ii++)
    {
      nsiz=2;
      nsiz1=2;

      h = fread(hea, nsiz, nsiz1, fr);
      j=hea[0];
       nsiz=8;
       nsiz1=(j)>>3;
       h = fread(ff+ii*31, nsiz, nsiz1, fr);
       nsiz=2;
       nsiz1=2;
       h = fread(hea, nsiz, nsiz1, fr);
    }

  for(ii=0;ii<192;ii++)
    {
      nsiz=2;
      nsiz1=2;

      h = fread(hea, nsiz, nsiz1, fr);
      j=hea[0];
       nsiz=8;
       nsiz1=(j)>>3;

       h = fread(ff1+ii*31, nsiz, nsiz1, fr);
       nsiz=2;
       nsiz1=2;
       h = fread(hea, nsiz, nsiz1, fr);
    }


  for(ii=0;ii<576;ii++)
    {
      nsiz=2;
      nsiz1=2;

      h = fread(hea, nsiz, nsiz1, fr);
      j=hea[0];
       nsiz=8;
       nsiz1=(j)>>3;

       h = fread(fg+ii*31, nsiz, nsiz1, fr);

       nsiz=2;
       nsiz1=2;
       h = fread(hea, nsiz, nsiz1, fr);
    }
  for(ii=0;ii<48;ii++)
    {
      nsiz=2;
      nsiz1=2;

      h = fread(hea, nsiz, nsiz1, fr);
      j=hea[0];
      nsiz=8;
       nsiz1=(j)>>3;

       h = fread(fg2+ii*31, nsiz, nsiz1, fr);
       nsiz=2;
       nsiz1=2;
       h = fread(hea, nsiz, nsiz1, fr);
    }
  for(ii=0;ii<960;ii++)
    {
      nsiz=2;
      nsiz1=2;
      h = fread(hea, nsiz, nsiz1, fr);
      j=hea[0];
       nsiz=8;
       nsiz1=(j)>>3;

       h = fread(fg5+ii*31, nsiz, nsiz1, fr);
       nsiz=2;
       nsiz1=2;
       h = fread(hea, nsiz, nsiz1, fr);
    }


  fclose(fr); 
  return;
}
