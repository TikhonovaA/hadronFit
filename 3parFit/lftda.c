#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include <stdlib.h>
#include <string.h>


void lftda_(short int *id,short int *f,short int *f1,short int *fg41,short int *fg43,short int *fg31,short int *fg32,short int *fg33,
          int *y,int *ttrig, int *n16, int *ch,int *mi5,int *lar,int *ltr,int *lq,int *lch3,int *lch4,int *lcp, int *ask,int *evna,int *ifp)
{
  static long long int k_np[16] = {
    65536,
    32768,
    21845,
    16384,
    13107,
    10923,
    9362,
    8192,
    7282,
    6554,
    5958,
    5461,
    5041,
    4681,
    4369,
    4096};
  long long int b32=4294967295;
  int A0  = (int)*(id+128+*ch)-128;
  int Askip  = (int)*(id+192+*ch)-128;

  int Ahard  = (int)*(id+64+*ch);
  int k_a = (int)*((unsigned char*)id+26);
  int k_b = (int)*((unsigned char*)id+27);
  int k_c = (int)*((unsigned char*)id+28);
  int k_16 = (int)*((unsigned char*)id+29);
  
  int k1_chi =  (int)*((unsigned char*)id+32);
  int k2_chi =  (int)*((unsigned char*)id+33);
  int chi_thres = (int)*(id+15);
  int s1,s2;

 
  long long int z0; 
  int it, it0;   
  long long d_it; 
  int it_h, it_l; 
  long long A1, B1, A2,C1,ch1,ch2,B2,B3,B5 ; 
  int low_ampl,i,T,iter;   


  if(*ifp == 1)printf("v lftda _______________%d %d \n",*evna,*ttrig);
  *ask=Ahard;
  if(k_16+*n16 !=16)
    {
      printf("disagreement in number of the points %d %d \n",k_16,*n16);
    }

  
  int validity_code=0;
  for(i=0, z0=0; i<16; i++)
    z0 += y[i];

  //initial time index
  it0 = 48 + ((23-*ttrig)<<2);
  //  it0 = 96 + ((12-*ttrig)<<3);

  //limits
  it_h=191;
  it_l=0;

  if (it0 < it_l)it0=it_l; 
  if (it0 > it_h)it0=it_h; 
  it=it0;

  //first approximation without time correction


  s1=(*(fg41+*ttrig*16));
  A2 = (s1 * z0);


  for(i=1;i<16;i++){
    s1=(*(fg41+*ttrig*16+i));
    B3 = y[15+i];
    B3 = s1 * B3;
    A2 += B3;
  }

  A2 += (1<<(k_a-1));
  A2 >>= k_a;
  T = it0<<4;


  //too large amplitude
  if(A2>262143){
    A1=A2>>3;
    validity_code=1;
    
    goto ou;
  }


  low_ampl = 0;
  if(A2>=A0){
    for(iter=0, it=it0; iter<3;){
      iter++;
      s1=(*(fg31+it*16));
      s2=(*(fg32+it*16));
      A1 = (s1 * z0);
      B1 = (s2 * z0);

      for(i=1;i<16;i++){
        s1=(*(fg31+i+it*16));
        s2=(*(fg32+i+it*16));
        B5=y[15+i];
        B5=s1 * B5;
        A1 += B5;

        B3=y[15+i];
        B3=s2 * B3;
        B1 += B3;
      }

      A1 += (1<<(k_a-1));

      A1 = A1>>k_a;


      if(A1>262143) 
        goto ou;

      if(A1<A0){
        low_ampl = 1;
        it=it0;

        goto lam;
      }

      if(iter != 3){
        B2 = B1>>(k_b-9);
        B1 = B2>>9;

        B2 += (A1<<9);


        B3=(B2/A1);


        it +=((B3+1)>>1)-256;
        it = it>it_h ? it_h : it;
        it = it<it_l ? it_l : it;
      } 
      else{
        B2 = B1>>(k_b-13);
        B5 = B1>>(k_b-9);
        
        B1 = B2>>13;
        B2 += (A1<<13);
        B3=(B2/A1);
        
        T = ((it)<<4) + ((B3+1)>>1)-4096;
        
        B1=B5>>9;
        B5 += (A1<<9);
        B3=(B5/A1);
        it +=((B3+1)>>1)-256;
        it = it>it_h ? it_h : it; 
        it = it<it_l ? it_l : it;

        T = T > 3071 ?  3071 : T;

        T = T < 0 ? 0 : T;
        C1 = (*(fg33+it*16) * z0);
        for(i=1;i<16;i++)
          C1 +=*(fg33+i+it*16) * y[15+i];
        C1 += (1 << (k_c-1));
        C1 >>= k_c;

      }

    } // for (iter...)
  } // if(A2>A0)
  else
    low_ampl = 1;

  if( low_ampl==1 ){
      
    lam:
      A1 = A2;
      validity_code=0;
      B1=0;
      C1 = (*(fg43+*ttrig*16) * z0);
      for(i=1;i<16;i++){
        B5=y[15+i];
        C1 += *(fg43+i+*ttrig*16) * B5;
      }
      C1 += (1 << (k_c-1));
      C1 >>= k_c;
      //ccf

  }    
  ch2=(A1* *(f+it*16)+B1* *(f1+it*16))>>k1_chi;
  ch2 +=C1;
  ch2 = z0-*n16*ch2;
  ch1=((ch2)*(ch2));
  ch1=ch1*k_np[*n16-1];
  ch1=ch1>>16;
  for(i=1;i<16;i++)
    {
      ch2=A1*(*(f+i+it*16))+B1*(*(f1+i+it*16));
      ch2 >>=k1_chi;
      ch2=(y[i+15]-ch2-C1);

      ch1=ch1 + ch2*ch2;

    }
  B2=(A1>>1)*(A1>>1);
  B2 >>=(k2_chi-2);
  B2 +=chi_thres;
  if (ch1 > B2)validity_code=3;
 ou:
  *mi5=it;
  *lar=A1;
  *ltr=T;

  int ss=(y[20]+y[21]);

  if(ss <= Ahard)validity_code=validity_code+4;


  *lq=validity_code;
  *lch3=B2;

  *lch4=ch1;

  *lcp=C1;

  return ;
}
