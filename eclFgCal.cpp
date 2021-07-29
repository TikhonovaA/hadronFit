#include <iostream>
#include <fstream>
#include <cassert>
#include <stdlib.h>
#include <string.h>
#include <cmath>  
#include <cstdlib>
//#include <TMath.h>
//#include <TF1.h>
#include <vector>
//#include <TCanvas.h>
//#include <TStyle.h>
//#include "TFile.h"
//#include "TChain.h"
//#include "TTree.h"
#include "eclFgCal.h"
//#include "TDirectory.h"
//#include "TDirectoryFile.h"
#include <Eigen/Dense>
#include "funcCalc.h"

using namespace std;


eclFgCal::eclFgCal(/*char* name ,  int cr1*/){
  
  wr_rootf = 0;
  
  
  
  //cr=cr1;
  /*  char histname[100];
  
  sprintf(histname, "Kf_%s_%d",name, cr);
  
  Kf = new TH1F(histname, "ev", 16, 0, 16);
  
  
  sprintf(histname, "Ka_%s_%d",name,cr);
  
  
  Ka = new TH1F(histname, "ev", 16, 0, 16);

  sprintf(histname, "Kb_%s_%d",name,cr);
  
  Kb = new TH1F(histname, "ev", 16, 0, 16);
  
  sprintf(histname, "Kc_%s_%d",name,cr);
  
  Kc = new TH1F(histname, "ev", 16, 0, 16);*/
  
  

  
  //  int nbins = (ampMax - ampMin)/10;
  for ( int ii=0; ii< 16; ii++){
    
   // mask_[ii]=false;
    
    bpf[ii]=35;
    bpf1[ii]=35;
    bhf[ii]=35;
    bhf1[ii]=35;
    bfg51[ii]=35;
    bfg52[ii]=35;
    bfg53[ii]=35;
    bfg54[ii]=35;
    bfg55[ii]=35;
    bfg41[ii]=35;
    bfg43[ii]=35;
    bfg45[ii]=35;
    

    /*
      bf[ii]=25;
      bf1[ii]=25;
      bfg31[ii]=25;
      bfg32[ii]=25;
      bfg33[ii]=25;
      bfg41[ii]=25;
      bfg43[ii]=25;
    */  
   for ( int jj=0; jj< 16; jj++){
     
     if(jj<4){
      mpf[ii][jj]=0.;
      mpf1[ii][jj]=0.;
      mhf[ii][jj]=0.;
      mhf1[ii][jj]=0.;
      mfg51[ii][jj]=0.;
      mfg52[ii][jj]=0.;
      mfg53[ii][jj]=0.;
      mfg54[ii][jj]=0.;
      mfg55[ii][jj]=0.;
      
      mfg41[ii][jj]=0.;
      mfg43[ii][jj]=0.;
      mfg45[ii][jj]=0.;
    }
    for ( int kk = 0; kk < 192; kk++){
      if( kk <16){
        IY[kk][ii][jj]=0.;
      }    
      p_f[ii][kk][jj]=0.;
      p_f1[ii][kk][jj]=0.;
      h_f[ii][kk][jj]=0.;
      h_f1[ii][kk][jj]=0.;

      fg51[ii][kk][jj]=0.;
      fg52[ii][kk][jj]=0.;
      fg53[ii][kk][jj]=0.;
      fg54[ii][kk][jj]=0.;
      fg55[ii][kk][jj]=0.;

      m_pf[ii][kk][jj]=0.;
      m_pf1[ii][kk][jj]=0.;
      m_hf[ii][kk][jj]=0.;
      m_hf1[ii][kk][jj]=0.;
      m_fg51[ii][kk][jj]=0;
      m_fg52[ii][kk][jj]=0;
      m_fg53[ii][kk][jj]=0;
      m_fg54[ii][kk][jj]=0;
      m_fg55[ii][kk][jj]=0;
      
      if( kk < 24 ){
        fg41[ii][kk][jj]=0.;
        fg43[ii][kk][jj]=0.;
        fg45[ii][kk][jj]=0.;
        m_fg41[ii][kk][jj]=0;
        m_fg43[ii][kk][jj]=0;
        m_fg45[ii][kk][jj]=0;
      }
      
     }
   }
   
   
   /*   for ( int jj=0; jj< 192; jj++){
     
     sprintf(histname, "uf_%s_%d_%d_%d",name,ii,jj,cr);

     uf[jj][ii] = new TH1F(histname, "ev", 16, 0, 16);
     
     sprintf(histname, "uf1_%s_%d_%d_%d",name,ii,jj,cr);
     
     uf1[jj][ii] = new TH1F(histname, "ev", 16, 0, 16);
     
     
     sprintf(histname, "uf31_%s_%d_%d_%d",name,ii,jj,cr);
     
     uf31[jj][ii] = new TH1F(histname, "ev", 16, 0, 16);
     
     sprintf(histname, "uf32_%s_%d_%d_%d",name,ii,jj,cr);
     
     uf32[jj][ii] = new TH1F(histname, "ev", 16, 0, 16);
     
     sprintf(histname, "uf33_%s_%d_%d_%d",name,ii,jj,cr);
     
     uf33[jj][ii] = new TH1F(histname, "ev", 16, 0, 16);

     if( jj <24){
       sprintf(histname, "uf41_%s_%d_%d_%d",name,ii,jj,cr);
       
       uf41[jj][ii] = new TH1F(histname, "ev", 16, 0, 16);
       
       sprintf(histname, "uf43_%s_%d_%d_%d",name,ii,jj,cr);
       
       uf43[jj][ii] = new TH1F(histname, "ev", 16, 0, 16);
   
       }
     
  }*/
   
   
  }
  
}


eclFgCal::~eclFgCal(){
//	if (inputFile != NULL) fclose(inputFile);
/*  printf("dtor\n");
  delete Kf;
  delete Ka;
  delete Kb;
  delete Kc;
  printf("dtor 1\n");
  for ( int ii=0; ii< 16; ii++){
    for ( int jj=0; jj< 192; jj++){  
      printf("dtor1 %d %d\n",ii,jj);
  //    Ms[ii]->Delete();
      delete uf[jj][ii];
      printf("dtor2 %d %d\n",ii,jj);
      if (uf1[jj][ii]) {
	printf("dtor2: testing uf1 %d %d\n",ii,jj);
	printf("dtor2 name %s\n",uf1[jj][ii]->GetName());
	delete uf1[jj][ii];
      }
      printf("dtor3 %d %d\n",ii,jj);
      delete uf31[jj][ii];
      printf("dtor4 %d %d\n",ii,jj);
      delete uf32[jj][ii];
      printf("dtor5 %d %d\n",ii,jj);
      delete uf33[jj][ii];
      printf("dtor6 %d %d\n",ii,jj);
      if (jj<24) {
	delete uf41[jj][ii];
	delete uf43[jj][ii];
      }
      printf("dtor7 %d %d\n",ii,jj);
    }
    }*/

}






long int myPow(long int x, int p)
{
  if (p == 0) return 1;
  if (p == 1) return x;
  return x * myPow(x, p - 1);
}


short  int bit_val(double f, int idd, short int bit, int a, int b, long int lim)
{
  short int i = bit;

  long int val;
  long int R;
  int flag=0;

  if (i < 11) {i = 11;}
  do {
    R = myPow(2, i);
    val = (long int)(a * R * f / idd / b + R + 0.5) - R;

    if ((val < lim && val > -lim)) {
      flag = 1;
    } else {
      i = i - 1;

    }

  } while (i > 0 && flag == 0);


  if (i > bit) {i = bit;}
  return i;

}




/*

double  eclFgCal::Sv123(double t, double t01, double tb1, double t02, double tb2, double td1, double ts1)
{

  double sv123 = 0.;
  double  dks0, dks1, dksm,
          dw0, dw1, dwp, dwm, das1, dac1, das0, dac0, dzna, dksm2, ds, dd,
          dcs0, dsn0, dzn0, td, ts, dr,
          dcs0s, dsn0s, dcs0d, dsn0d, dcs1s, dsn1s, dcs1d, dsn1d;


  if (t < 0.) return 0.;

  dr = (ts1 - td1) / td1;
  if (fabs(dr) >= 0.00001) {
    td = td1;
    ts = ts1;
  } else {
    td = td1;
    if (ts1 > td1) {
      ts = td1 * 1.00001;
    } else {
      ts = td1 * 0.99999;
    }
  }



  dr = ((t01 - t02) * (t01 - t02) + (tb1 - tb2) * (tb1 - tb2)) / ((t01) * (t01) + (tb1) * (tb1));
  dks0 = 1.0 / t01;
  dks1 = 1.0 / t02;

  if (dr < 0.0000000001) {

    if (dks0 > dks1) {
      dks0 = dks1 * 1.00001;
    } else {
      dks0 = dks1 * 0.99999;
    }
  }

  //  printf(" ts=%e dr=%e dks0=%e dks1=%e \n",ts,dr,dks0,dks1);

  if (t < 0.) return 0;



  dksm = dks1 - dks0;

  ds = 1. / ts;
  dd = 1. / td;

  dw0 = 1. / tb1;
  dw1 = 1. / tb2;
  dwp = dw0 + dw1;
  dwm = dw1 - dw0;

  dksm2 = dksm * dksm;

  dzna = (dksm2 + dwm * dwm) * (dksm2 + dwp * dwp);


  das0 = dw1 * (dksm2 + dwp * dwm);
  dac0 = -2 * dksm * dw0 * dw1;
  das1 = dw0 * (dksm2 - dwp * dwm);
  dac1 = -dac0;





  dsn0 = (ds - dks0);
  dcs0 = -dw0;
  dzn0 = dcs0 * dcs0 + dsn0 * dsn0;

  dsn0s = (dsn0 * das0 - dcs0 * dac0) / dzn0;
  dcs0s = (dcs0 * das0 + dsn0 * dac0) / dzn0;

  dsn0 = (ds - dks1);
  dcs0 = -dw1;
  dzn0 = dcs0 * dcs0 + dsn0 * dsn0;

  dsn1s = (dsn0 * das1 - dcs0 * dac1) / dzn0;
  dcs1s = (dcs0 * das1 + dsn0 * dac1) / dzn0;


  dsn0 = (dd - dks0);
  dcs0 = -dw0;
  dzn0 = dcs0 * dcs0 + dsn0 * dsn0;

  dsn0d = (dsn0 * das0 - dcs0 * dac0) / dzn0;
  dcs0d = (dcs0 * das0 + dsn0 * dac0) / dzn0;

  dsn0 = (dd - dks1);
  dcs0 = -dw1;
  dzn0 = dcs0 * dcs0 + dsn0 * dsn0;

  dsn1d = (dsn0 * das1 - dcs0 * dac1) / dzn0;
  dcs1d = (dcs0 * das1 + dsn0 * dac1) / dzn0;

  //cppcheck dr = (ts - td) / td;


  sv123 = ((((dsn0s - dsn0d) * sin(dw0 * t)
             + (dcs0s - dcs0d) * cos(dw0 * t)) * exp(-t * dks0)
            - (dcs0s + dcs1s) * exp(-t * ds) + (dcs0d + dcs1d) * exp(-t * dd)
            + ((dsn1s - dsn1d) * sin(dw1 * t)
               + (dcs1s - dcs1d) * cos(dw1 * t)) * exp(-t * dks1)) / dzna / (ts - td));

  //    printf("sv123 %18.15e %18.15e \n",t,sv123);
  
        
  sv123=sv123/(-.109+.919*t01-.261*t01*t01)
    /(-.109+.919*t02-.261*t02*t02)
    /(.262+.174*tb1-.208*tb1*tb1)
    /(.262+.174*tb2-.208*tb2*tb2)
    /(4.56-1.58*td1)/(1.391-0.434*ts1)
    /(1.06-0.578*(t01-tb1)*(t01-tb1))
    /(1.06-0.578*(t02-tb2)*(t02-tb2))
    /(1.2140-0.79645*t01+0.63440*t01*t01)
    /(1.2140-0.79645*t02+0.63440*t02*t02);
  
  //    printf("***sv123 %18.15e %18.15e \n",t,sv123);
  return sv123;


}


double eclFgCal::ShaperDSP_F(double Ti, double* ss)
{
  double svp = 0;
  //  double R =  127.;
  //  double Z =  144.;
  //  double scale =  R/Z;
  double scale =  127./144.;
      double tr1 = Ti * scale;
  //    double tr1 = Ti;

  //    double s[12] = {0, 27.7221, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double s[12] = {0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};


  s[1] =  * (ss + 10);  //
  s[2] =  * (ss + 0);
  s[3] =  * (ss + 1);
  s[4] =  * (ss + 2);
  s[5] =  * (ss + 3);
  s[6] =  * (ss + 4);
  s[7] =  * (ss + 5);
  s[8] =  * (ss + 6);
  s[9] =  * (ss + 7);
  s[10] = * (ss + 8);
  s[11] = * (ss + 9);

                                                                            
  double tr = tr1 - s[2];
  double tr2 = tr + .2;
  double tr3 = tr - .2;
  if (tr2 > 0.) {
    
    svp = (Sv123(tr , s[4], s[5], s[9], s[10], s[3], s[6]) * (1 - s[11])
           + s[11] * .5 * (Sv123(tr2, s[4], s[5], s[9], s[10], s[3], s[6])
                           + Sv123(tr3, s[4], s[5], s[9], s[10], s[3], s[6])));

    //    printf("do svp=  %21.18lf %21.18lf %21.18lf %21.18lf %21.18lf \n ",svp,tr,tr1,tr2,tr3);
    
    //    svp = Sv123(tr , s[4], s[5], s[9], s[10], s[3], s[6]) * (1 - s[11])
    //     + s[11] *Sv123(tr2, s[4], s[5], s[9], s[10], s[3], s[6]);
    double x = tr / s[4];


    svp = s[1] * (svp - s[7] * (exp(-tr / s[8]) *
                                (1 - exp(-x) * (1 + x + x * x / 2 + x * x * x / 6 + x * x * x * x / 24 + x * x * x * x * x / 120))));
    //    printf("posle svp=   %21.18lf \n ",svp);

  } else svp = 0 ;

  //    printf("par %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %e \n ",tr,s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9],s[10],s[11],svp);
  //      printf("par %21.18lf %21.18lf   \n ",tr1,svp);
  return svp;

}





double eclFgCal::ShaperDSP(double Ti)
{
  double svp = 0;
  double tr1 = Ti * 127./144.;

  double s[12] = {0, 27.7221, 0.5, 0.6483, 0.4017, 0.3741, 0.8494, 0.00144547, 4.7071, 0.8156, 0.5556, 0.2752};


  double tr = tr1 - s[2];
  double tr2 = tr + .2;
  double tr3 = tr - .2;
  if (tr2 > 0.) {

    svp = (Sv123(tr , s[4], s[5], s[9], s[10], s[3], s[6]) * (1 - s[11])
           + s[11] * .5 * (Sv123(tr2, s[4], s[5], s[9], s[10], s[3], s[6])
                           + Sv123(tr3, s[4], s[5], s[9], s[10], s[3], s[6])));
    double x = tr / s[4];


    svp = s[1] * (svp - s[7] * (exp(-tr / s[8]) *
                                (1 - exp(-x) * (1 + x + x * x / 2 + x * x * x / 6 + x * x * x * x / 24 + x * x * x * x * x / 120))));
  } else svp = 0 ;
  return svp;

}

*/

void eclFgCal::CalculateFgpar( ){

  double dt;
  double ts0;
  double del;


  int j1, endc, /*j,*/ i;
  double ndt;
  double adt, dt0, t0, ddt, tc1;

  double p_ssssj, p_ssssj1, p_ssssi, p_ssssi1, h_ssssj, h_ssssj1, h_ssssi, h_ssssi1;
  double svp, svm;

  double t[1], tmd[1], tpd[1];

  /*  double MP[10] = {0.5, 0.6483, 0.4017, 0.3741, 0.8494, 0.00144547, 4.7071, 0.8156, 0.5556, 0.2752};*/

  double MZ[11] = {0.75, 0.648324, 0.401711, 0.374167, 0.849417, 0.00144548, 4.70722, 0.815639, 0.555605, 0.2752, 27.7221};
  double spt=0;
  double fm;
  double fc;
 
  // varible for one channel                                                                                                                                                                  

  double  gg[192], gg1[192], g1g1[192], gg2[192], g1g2[192], g2g2[192], gg3[192], g1g3[192], g2g3[192], g3g3[192], gg4[192], g1g4[192], g2g4[192], g3g4[192], g4g4[192]; //, dgg[192];
  double  sg1[16][192], sg[16][192], sg2[16][192], sg3[16][192], sg4[16][192];
  double  dgg1[192];
  int chan;

  double scale =  144./127.;

  for ( chan = 0; chan < 16; chan++){


//for photon component

    /*for (j1 = 0; j1 < 10; j1++) {
      MZ[j1]=p_par[chan][j1]; //par for one chan
    }
    */

    MZ[0]=0.75;
    MZ[1]=0.648324;
    MZ[2]=0.401711;
    MZ[3]=0.374167;
    MZ[4]=0.849417;
    MZ[5]=0.00144548;
    MZ[6]=4.70722;
    MZ[7]=0.815639;
    MZ[8]=0.555605;
    MZ[9]=0.2752;
    MZ[10]=1.0;

    if(MZ[2] >= MZ[3] && MZ[2] >= MZ[4] && MZ[2] >= MZ[5] ){
      spt = MZ[2];
    }
    if(MZ[3] >= MZ[2] && MZ[3] >= MZ[4] && MZ[3] >= MZ[5] ){
      spt = MZ[3];
    }
    if(MZ[4] >= MZ[2] && MZ[4] >= MZ[3] && MZ[4] >= MZ[5] ){
      spt = MZ[4];
    }
    if(MZ[5] >= MZ[2] && MZ[5] >= MZ[3] && MZ[5] >= MZ[4] ){
      spt = MZ[5];
    }

    //    printf("spt=%e \n",spt);

    spt=2.*spt/1000.;
    //        printf("spt=%15.12e \n",spt);
    t[0]=0.;
    fm=0.;
    //    printf("CalculateFGPar p2\n");
    for (j1 = 0; j1 < 4000; j1++) {
      //printf("j== %d %15.12e \n",j1+1,t);
      fc=ShaperDSP_F(t, MZ);
      //fc=ShaperDSP_F(t, MZ);
      if(fc>fm){
        fm=fc;
        //     tm=t;
      }
      t[0]=t[0]+spt;
    }

 //   printf("fm=%lf \n",fm);


    //printf("do chan=%d FM=%21.18e TM=%21.18e \n",chan,fm,tm);

    MZ[10]=1.0/fm;
    t[0]=0.;
    fm=0.;
    for (j1 = 0; j1 < 4000; j1++) {
      fc=ShaperDSP_F(t, MZ);
      // fc=ShaperDSP_F(t, MZ);
      if(fc>fm){
        fm=fc;
      }
      t[0]=t[0]+spt;
    }

  //  printf("FM=%15.12e \n",fm);

    del = 0.;
    ts0 = MZ[0];
    dt = 0.5;


    ndt = 96.;
    adt = 1. / ndt;
    endc = 2 * ndt;


    dt0 = adt * dt;
    t0 = -dt0 * ndt;

    
    //  printf("CalculateFGPar p3\n");
    for (j1 = 0; j1 < endc; j1++) { //from 0 to 191
      t0 = t0 + dt0;
      t[0] = t0 - dt - del;
      
      for (int j = 0; j < 16; j++) {
        t[0] = t[0] + dt;

        if (t[0] > 0) {
          p_f[chan][j1][j] = ShaperDSP_F(t, MZ);

          ddt = 0.005 * dt;
          tc1 = t[0] - ts0;
          tpd[0] = t[0] + ddt;
          tmd[0] = t[0] - ddt;

        
          if (tc1 > ddt) {
            svp = ShaperDSP_F(tpd, MZ);
            svm = ShaperDSP_F(tmd, MZ);
    

            p_f1[chan][j1][j] = (svp - svm) / 2. / ddt;
          } 
          else {
            if (tc1 > 0 ) {
              p_f1[chan][j1][j] = ShaperDSP_F(tpd, MZ) / (ddt + tc1);
            }
            else {
              p_f1[chan][j1][j] = 0.;
            }
          }// else tc1>ddt                                                                                                                                                                    

        } //if t>0                                                                                                                                                                            
        else {
          p_f[chan][j1][j] = 0.;
          p_f1[chan][j1][j] = 0.;
        }
     
      } //for j      

    }//for j1    




//for hadron component

   /* for (j1 = 0; j1 < 10; j1++) {
      MZ[j1]=h_par[chan][j1]; //par for one chan
    }*/

    MZ[0]=0.792623;
    MZ[1]=0.929354;
    MZ[2]=0.556139;
    MZ[3]=0.446967;
    MZ[4]=0.140175;
    MZ[5]=0.0312971;
    MZ[6]=3.12842;
    MZ[7]=0.791012;
    MZ[8]=0.619416;
    MZ[9]=0.385621;
    MZ[10]=29.5092;


    MZ[10]=1.0;

    if(MZ[2] >= MZ[3] && MZ[2] >= MZ[4] && MZ[2] >= MZ[5] ){
      spt = MZ[2];
    }
    if(MZ[3] >= MZ[2] && MZ[3] >= MZ[4] && MZ[3] >= MZ[5] ){
      spt = MZ[3];
    }
    if(MZ[4] >= MZ[2] && MZ[4] >= MZ[3] && MZ[4] >= MZ[5] ){
      spt = MZ[4];
    }
    if(MZ[5] >= MZ[2] && MZ[5] >= MZ[3] && MZ[5] >= MZ[4] ){
      spt = MZ[5];
    }

    //    printf("spt=%e \n",spt);

    spt=2.*spt/1000.;
    //        printf("spt=%15.12e \n",spt);
    t[0]=0.;
    fm=0.;
    //    printf("CalculateFGPar p2\n");
    for (j1 = 0; j1 < 4000; j1++) {
      //printf("j== %d %15.12e \n",j1+1,t);
      fc=ShaperDSP_F(t, MZ);
      //fc=ShaperDSP_F(t, MZ);
      if(fc>fm){
        fm=fc;
        //     tm=t;
      }
      t[0]=t[0]+spt;
    }

  //  printf("fm=%lf \n",fm);

    //printf("do chan=%d FM=%21.18e TM=%21.18e \n",chan,fm,tm);

    MZ[10]=1.0/fm;
    t[0]=0.;
    fm=0.;
    for (j1 = 0; j1 < 4000; j1++) {
      fc=ShaperDSP_F(t, MZ);
      // fc=ShaperDSP_F(t, MZ);
      if(fc>fm){
        fm=fc;
      }
      t[0]=t[0]+spt;
    }

   // printf("FM=%15.12e \n",fm);

    del = 0.;
    ts0 = MZ[0];
    dt = 0.5;


    ndt = 96.;
    adt = 1. / ndt;
    endc = 2 * ndt;


    dt0 = adt * dt;
    t0 = -dt0 * ndt;

    
    //  printf("CalculateFGPar p3\n");
    for (j1 = 0; j1 < endc; j1++) { //from 0 to 191
      t0 = t0 + dt0;
      t[0] = t0 - dt - del;

    
      for (int j = 0; j < 16; j++) {
        t[0] = t[0] + dt;

        if (t[0] > 0) {
          h_f[chan][j1][j] = ShaperDSP_F(t, MZ);
      
          ddt = 0.005 * dt;
          tc1 = t[0] - ts0;
          tpd[0] = t[0] + ddt;
          tmd[0] = t[0] - ddt;

        
          if (tc1 > ddt) {
            svp = ShaperDSP_F(tpd, MZ);
            svm = ShaperDSP_F(tmd, MZ);

            h_f1[chan][j1][j] = (svp - svm) / 2. / ddt;
         } 
          else {
            if (tc1 > 0 ) {
              h_f1[chan][j1][j] = ShaperDSP_F(tpd, MZ) / (ddt + tc1);
            }
            else {
              h_f1[chan][j1][j] = 0.;
            }
          }// else tc1>ddt                                                                                                                                                                    

        } //if t>0                                                                                                                                                                            
        else {
          h_f[chan][j1][j] = 0.;
          h_f1[chan][j1][j] = 0.;
        }

      } //for j                                                                                                                                                                               
     
    }//for j1    


   // TCanvas *mlpa_canvas = new TCanvas("mlpa_canvas","",1000,700);
//TGraph *pilegph = new TGraph(192);
	



    for (j1 = 0; j1 < endc; j1++) { //from 0 to 191

      gg[j1] = 0.;
      gg1[j1] = 0.;
      g1g1[j1] = 0.;

      gg2[j1] = 0.;
      g1g2[j1] = 0.;
      g2g2[j1] = 0.;

      gg3[j1] = 0.;
      g1g3[j1] = 0.;
      g2g3[j1] = 0.;
      g3g3[j1] = 0.;

      gg4[j1] = 0.;
      g1g4[j1] = 0.;
      g2g4[j1] = 0.;
      g3g4[j1] = 0.;
      g4g4[j1] = 0.;

      for (int j = 0; j < 16; j++) {
        sg[j][j1] = 0.;
        sg1[j][j1] = 0.;
        sg2[j][j1] = 0.;
        sg3[j][j1] = 0.;
        sg4[j][j1] = 0.;

        p_ssssj1 = p_f1[chan][j1][j];
        p_ssssj = p_f[chan][j1][j];
        h_ssssj1 = h_f1[chan][j1][j];
        h_ssssj = h_f[chan][j1][j];


        if(j==0){

          if(p_f[chan][j1][j]>mpf[chan][0]){mpf[chan][0]=p_f[chan][j1][j];}
          if(p_f[chan][j1][j]<mpf[chan][1]){mpf[chan][1]=p_f[chan][j1][j];}
          if(p_f1[chan][j1][j]>mpf1[chan][0]){mpf1[chan][0]=p_f1[chan][j1][j];}
          if(p_f1[chan][j1][j]<mpf1[chan][1]){mpf1[chan][1]=p_f1[chan][j1][j];}

          if(h_f[chan][j1][j]>mhf[chan][0]){mhf[chan][0]=h_f[chan][j1][j];}
          if(h_f[chan][j1][j]<mhf[chan][1]){mhf[chan][1]=h_f[chan][j1][j];}
          if(h_f1[chan][j1][j]>mhf1[chan][0]){mhf1[chan][0]=h_f1[chan][j1][j];}
          if(h_f1[chan][j1][j]<mhf1[chan][1]){mhf1[chan][1]=h_f1[chan][j1][j];}

        }

        else{

          if(p_f[chan][j1][j]>mpf[chan][2]){mpf[chan][2]=p_f[chan][j1][j];}
          if(p_f[chan][j1][j]<mpf[chan][3]){mpf[chan][3]=p_f[chan][j1][j];}
          if(p_f1[chan][j1][j]>mpf1[chan][2]){mpf1[chan][2]=p_f1[chan][j1][j];}
          if(p_f1[chan][j1][j]<mpf1[chan][3]){mpf1[chan][3]=p_f1[chan][j1][j];}

          if(h_f[chan][j1][j]>mhf[chan][2]){mhf[chan][2]=h_f[chan][j1][j];}
          if(h_f[chan][j1][j]<mhf[chan][3]){mhf[chan][3]=h_f[chan][j1][j];}
          if(h_f1[chan][j1][j]>mhf1[chan][2]){mhf1[chan][2]=h_f1[chan][j1][j];}
          if(h_f1[chan][j1][j]<mhf1[chan][3]){mhf1[chan][3]=h_f1[chan][j1][j];}

        }


        for (i = 0; i < 16; i++) {

          //	printf("%lf ",IY[chan][j][i]);

          sg[j][j1] = sg[j][j1] + IY[chan][j][i] * p_f[chan][j1][i];
          sg1[j][j1] = sg1[j][j1] + IY[chan][j][i] * p_f1[chan][j1][i];
          sg2[j][j1] = sg2[j][j1] + IY[chan][j][i] * h_f[chan][j1][i];
          sg3[j][j1] = sg3[j][j1] + IY[chan][j][i] * h_f1[chan][j1][i];
          sg4[j][j1] = sg4[j][j1] + IY[chan][j][i];



         // ssssi = f[chan][j1][i];
          //ssssi1 = f1[chan][j1][i];


          p_ssssi1 = p_f1[chan][j1][i];
          p_ssssi = p_f[chan][j1][i];
          h_ssssi1 = h_f1[chan][j1][i];
          h_ssssi = h_f[chan][j1][i];

          gg[j1] = gg[j1] + IY[chan][j][i] * p_ssssj * p_ssssi;
          gg1[j1] = gg1[j1] + IY[chan][j][i] * p_ssssj * p_ssssi1;
          gg2[j1] = gg2[j1] + IY[chan][j][i] * p_ssssj * h_ssssi;
          gg3[j1] = gg3[j1] + IY[chan][j][i] * p_ssssj * h_ssssi1;
          gg4[j1] = gg4[j1] + IY[chan][j][i] * p_ssssj;

          g1g1[j1] = g1g1[j1] + IY[chan][j][i] * p_ssssj1 * p_ssssi1;
          
          g1g2[j1] = g1g2[j1] + IY[chan][j][i] * p_ssssj1 * h_ssssi;
          g2g2[j1] = g2g2[j1] + IY[chan][j][i] * h_ssssj * h_ssssi;

          g1g3[j1] = g1g3[j1] + IY[chan][j][i] * p_ssssj1 * h_ssssi1;
          g2g3[j1] = g2g3[j1] + IY[chan][j][i] * h_ssssj * h_ssssi1;
          g3g3[j1] = g3g3[j1] + IY[chan][j][i] * h_ssssj1 * h_ssssi1;

          g1g4[j1] = g1g4[j1] + IY[chan][j][i] * p_ssssj1;
          g2g4[j1] = g2g4[j1] + IY[chan][j][i] * h_ssssj;
          g3g4[j1] = g3g4[j1] + IY[chan][j][i] * h_ssssj1;
          //if(j1>50&&j1<60&&chan==15) std::cout<<"j1= "<<j1<<", i= "<<i<<", j= "<<j<<", "<<g3g4[j1]<<", "<<IY[chan][j][i]<<", "<<h_ssssj1<<std::endl;
          g4g4[j1] = g4g4[j1] + IY[chan][j][i];


        }   // for i                                                                                                                                                                          
        //	printf("\n ");

      } //for j            
      //  printf("CalculateFGPar p4\n");



      dgg1[j1] = gg[j1] * g2g2[j1] * g4g4[j1] - gg[j1] * g2g4[j1] * g2g4[j1] + 2 * gg2[j1] * g2g4[j1] * gg4[j1] - gg2[j1] * gg2[j1] *
          g4g4[j1] - gg4[j1] * g2g2[j1] * gg4[j1];


    //  double sum=0;
      for (i = 0; i < 16; i++) {

        Eigen::MatrixXd matrix(5, 5);
        matrix(0,0) = gg[j1];
        matrix(1,1) = g1g1[j1];
        matrix(2,2) = g2g2[j1];
        matrix(3,3) = g3g3[j1];
        matrix(4,4) = g4g4[j1];

        matrix(0,1) = matrix(1,0) = gg1[j1];     
        matrix(0,2) = matrix(2,0) = gg2[j1];    
        matrix(0,3) = matrix(3,0) = gg3[j1]; 
        matrix(0,4) = matrix(4,0) = gg4[j1]; 

        matrix(1,2) = matrix(2,1) = g1g2[j1]; 
        matrix(1,3) = matrix(3,1) = g1g3[j1]; 
        matrix(1,4) = matrix(4,1) = g1g4[j1]; 

        matrix(2,3) = matrix(3,2) = g2g3[j1]; 
        matrix(2,4) = matrix(4,2) = g2g4[j1]; 

        matrix(3,4) = matrix(4,3) = g3g4[j1]; 


        Eigen::MatrixXd matrixI = matrix.llt().solve(Eigen::MatrixXd::Identity(5,5));

        //std::cout << matrix*matrixI << std::endl;

        Eigen::VectorXd vector(5);
        vector(0) = sg[i][j1];
        vector(1) = sg1[i][j1];
        vector(2) = sg2[i][j1];
        vector(3) = sg3[i][j1];
        vector(4) = sg4[i][j1];

        Eigen::VectorXd fgpar(5);

        fgpar = matrixI*vector;

        fg51[chan][j1][i] = fgpar(0);
        fg52[chan][j1][i] = fgpar(1);
        fg53[chan][j1][i] = fgpar(2);
        fg54[chan][j1][i] = fgpar(3);
        fg55[chan][j1][i] = fgpar(4);
        
        //if(i==0&&chan==15) pilegph->SetPoint(j1, j1,fg52[chan][j1][i]);

        //sum+=fg55[chan][j1][i]*h_f1[chan][j1][i];

       // std::cout << " fgpar :\n" << fgpar << std::endl;


        if(i==0){
          if(fg51[chan][j1][i]>mfg51[chan][0]){mfg51[chan][0]=fg51[chan][j1][i];}
          if(fg51[chan][j1][i]<mfg51[chan][1]){mfg51[chan][1]=fg51[chan][j1][i];}
          if(fg52[chan][j1][i]>mfg52[chan][0]){mfg52[chan][0]=fg52[chan][j1][i];}
          if(fg52[chan][j1][i]<mfg52[chan][1]){mfg52[chan][1]=fg52[chan][j1][i];}
          if(fg53[chan][j1][i]>mfg53[chan][0]){mfg53[chan][0]=fg53[chan][j1][i];}
          if(fg53[chan][j1][i]<mfg53[chan][1]){mfg53[chan][1]=fg53[chan][j1][i];}
          if(fg54[chan][j1][i]>mfg54[chan][0]){mfg54[chan][0]=fg54[chan][j1][i];}
          if(fg54[chan][j1][i]<mfg54[chan][1]){mfg54[chan][1]=fg54[chan][j1][i];}
          if(fg55[chan][j1][i]>mfg55[chan][0]){mfg55[chan][0]=fg55[chan][j1][i];}
          if(fg55[chan][j1][i]<mfg55[chan][1]){mfg55[chan][1]=fg55[chan][j1][i];}
        }
        else{
          if(fg51[chan][j1][i]>mfg51[chan][2]){mfg51[chan][2]=fg51[chan][j1][i];}
          if(fg51[chan][j1][i]<mfg51[chan][3]){mfg51[chan][3]=fg51[chan][j1][i];}
          if(fg52[chan][j1][i]>mfg52[chan][2]){mfg52[chan][2]=fg52[chan][j1][i];}
          if(fg52[chan][j1][i]<mfg52[chan][3]){mfg52[chan][3]=fg52[chan][j1][i];}
          if(fg53[chan][j1][i]>mfg53[chan][2]){mfg53[chan][2]=fg53[chan][j1][i];}
          if(fg53[chan][j1][i]<mfg53[chan][3]){mfg53[chan][3]=fg53[chan][j1][i];}
          if(fg54[chan][j1][i]>mfg54[chan][2]){mfg54[chan][2]=fg54[chan][j1][i];}
          if(fg54[chan][j1][i]<mfg54[chan][3]){mfg54[chan][3]=fg54[chan][j1][i];}
          if(fg55[chan][j1][i]>mfg55[chan][2]){mfg55[chan][2]=fg55[chan][j1][i];}
          if(fg55[chan][j1][i]<mfg55[chan][3]){mfg55[chan][3]=fg55[chan][j1][i];}

        }

        
    
        /*		uf31[j1][chan]->SetBinContent(i+1,fg31[chan][j1][i]);
        uf32[j1][chan]->SetBinContent(i+1,fg32[chan][j1][i]);
        uf33[j1][chan]->SetBinContent(i+1,fg33[chan][j1][i]);*/
        
      }  // for i 0-16                                                                                                                                                                        
      

  // printf("sum=%lf\n", sum);

      int jk = 23 + ((48 - j1 ) >> 2);
      //  printf( " !jk=%d j1=%d \n",jk,j1);
      //  printf("CalculateFGPar p5\n");
      if (jk >= 0 && jk < 24 && (48 - j1) % 4 == 0) {
        // printf( "jk=%d \n",jk);
       // double sum=0;
        for(i=0;i<16;i++){
          if (dgg1[j1] != 0) {
      

            fg41[chan][jk][i] = ((g2g2[j1] * g4g4[j1] - g2g4[j1] * g2g4[j1]) * sg[i][j1] + (-gg2[j1] * g4g4[j1] + g2g4[j1] * gg4[j1]) * sg2[i][j1] +
                        (gg2[j1] * g2g4[j1] - g2g2[j1] * gg4[j1]) * sg4[i][j1]) / dgg1[j1];


            fg43[chan][jk][i] = ((-gg2[j1] * g4g4[j1] + gg4[j1] * g2g4[j1]) * sg[i][j1] + (gg[j1] * g4g4[j1] - gg4[j1] * gg4[j1]) * sg2[i][j1] +
                          (-gg[j1] * g2g4[j1] + gg2[j1] * gg4[j1]) * sg4[i][j1]) / dgg1[j1];


            fg45[chan][jk][i] = ((gg2[j1] * g2g4[j1] - gg4[j1] * g2g2[j1]) * sg[i][j1] + (-gg[j1] * g2g4[j1] + gg4[j1] * gg2[j1]) * sg2[i][j1] +
                          (gg[j1] * g2g2[j1] - gg2[j1] * gg2[j1]) * sg4[i][j1]) / dgg1[j1];

          // printf("%lf %lf %lf \n", fg41[chan][jk][i], fg43[chan][jk][i], fg45[chan][jk][i]);


          //  sum+=fg45[chan][jk][i]*h_f[chan][j1][i];

            if(i==0){
              if(fg41[chan][jk][i]>mfg41[chan][0]){mfg41[chan][0]=fg41[chan][jk][i];}
              if(fg41[chan][jk][i]<mfg41[chan][1]){mfg41[chan][1]=fg41[chan][jk][i];}
              if(fg43[chan][jk][i]>mfg43[chan][0]){mfg43[chan][0]=fg43[chan][jk][i];}
              if(fg43[chan][jk][i]<mfg43[chan][1]){mfg43[chan][1]=fg43[chan][jk][i];}
              if(fg45[chan][jk][i]>mfg45[chan][0]){mfg45[chan][0]=fg45[chan][jk][i];}
              if(fg45[chan][jk][i]<mfg45[chan][1]){mfg45[chan][1]=fg45[chan][jk][i];}
            }
            else{
              if(fg41[chan][jk][i]>mfg41[chan][2]){mfg41[chan][2]=fg41[chan][jk][i];}
              if(fg41[chan][jk][i]<mfg41[chan][3]){mfg41[chan][3]=fg41[chan][jk][i];}
              if(fg43[chan][jk][i]>mfg43[chan][2]){mfg43[chan][2]=fg43[chan][jk][i];}
              if(fg43[chan][jk][i]<mfg43[chan][3]){mfg43[chan][3]=fg43[chan][jk][i];}
              if(fg45[chan][jk][i]>mfg45[chan][2]){mfg45[chan][2]=fg45[chan][jk][i];}
              if(fg45[chan][jk][i]<mfg45[chan][3]){mfg45[chan][3]=fg45[chan][jk][i];}

            }

          }
          
          //      	uf41[jk][chan]->SetBinContent(i+1,fg41[chan][jk][i]);
          //      	uf43[jk][chan]->SetBinContent(i+1,fg43[chan][jk][i]);
        }  //for i 0-16

       // printf("chan=%d jk=%d sum=%lf \n", chan, jk, sum);
  
      } //                                                              

    } // for j1 <endc           

    //pilegph->Draw();
                                                                             

  } //for chan

//printf("%lf", fg41[15][23][3]);

}




void eclFgCal::CalculateBit( ){

  int /*j1, endc, j,*/ i;
  int ChN;
  long int i16;
  int idd;
  double val_f;


  // double fdd; 
  //  int ibb, iaa, idd, icc;
  //int ia, ib, ic, i16, ilim;
  //  int mxf, mxf1, mxfg31, mxfg32, mxfg33, mxfg41, mxfg43;
  // double xf;
  //  double xf1;
  //double xfg31;
  //double xfg32;
  //double xfg33;
  //double xfg41;
  //double xfg43;

  short int bitpf;                                                                                            
  short int bitpf1;         
  short int bithf;                                                                                            
  short int bithf1;                                                                                       
  short int bitfg51;                                                                                         
  short int bitfg52;                                                                                         
  short int bitfg53;    
  short int bitfg54;
  short int bitfg55;

  short int bitfg41;                                                                                         
  short int bitfg43;
  short int bitfg45;

 // short int bad[16];

 // int badN;

  int cpf;
  int cpf1;
  int chf;
  int chf1;
  int cf51;
  int cf52;
  int cf53;
  int cf54;
  int cf55;
  int cf41;
  int cf43;
  int cf45;
  

  cpf = 1;
  cpf1 = 1;
  chf = 1;
  chf1 = 1;
  cf51 = 1;
  cf52 = 1;
  cf53 = 1;
  cf54 = 1;
  cf55 = 1;
  cf41 = 1;
  cf43 = 1;
  cf45 = 1;


  //  int ipardsp13;
  int ipardsp14;

  int n16;
  int k;


  for (ChN = 0; ChN < 16; ChN++) {
   // if (!mask_[ChN]) continue;
    n16 = 16;

    i16 = myPow(2, 31);

   
    for (i = 0; i < 2; i++) {
      
       if (i == 0) {
        idd = n16;
      } else {
        idd = 1;
      }
    
      for (k = 0; k < 2; k++) {

        //      printf("ch=%d %d %lf \n ",ChN,2*i+k,mfg32[ChN][2*i+k]);
        val_f = mpf[ChN][2*i+k];
        bpf[ChN] = bit_val(val_f, idd, bpf[ChN],  1, 1, i16);
        if (bpf[ChN] < bitpf) {bitpf = bpf[ChN];}
        val_f = mpf1[ChN][2*i+k];
        bpf1[ChN] = bit_val(val_f, idd, bpf1[ChN],  4, 3, i16);
        if (bpf1[ChN] < bitpf1) {bitpf1 = bpf1[ChN];}
        val_f = mhf[ChN][2*i+k];
        bhf[ChN] = bit_val(val_f, idd, bhf[ChN],  1, 1, i16);
        if (bhf[ChN] < bithf) {bithf = bhf[ChN];}
        val_f = mhf1[ChN][2*i+k];
        bhf1[ChN] = bit_val(val_f, idd, bhf1[ChN],  4, 3, i16);
        if (bhf1[ChN] < bithf1) {bithf1 = bhf1[ChN];}
      
        val_f = mfg51[ChN][2*i+k];
        bfg51[ChN] = bit_val(val_f, idd, bfg51[ChN],  1, 1, i16);
        if (bfg51[ChN] < bitfg51) {bitfg51 = bfg51[ChN];}
        
        val_f = mfg52[ChN][2*i+k];
        bfg52[ChN] = bit_val(val_f, idd, bfg52[ChN],  3, 4, i16);
        if (bfg52[ChN] < bitfg52) {bitfg52 = bfg52[ChN];}

        val_f = mfg53[ChN][2*i+k];
        bfg53[ChN] = bit_val(val_f, idd, bfg53[ChN],  1, 1, i16);
        if (bfg53[ChN] < bitfg53) {bitfg53 = bfg53[ChN];}
        
        val_f = mfg54[ChN][2*i+k];
        bfg54[ChN] = bit_val(val_f, idd, bfg54[ChN],  3, 4, i16);
        if (bfg54[ChN] < bitfg54) {bitfg54 = bfg54[ChN];}
        
        val_f = mfg55[ChN][2*i+k];
        bfg55[ChN] = bit_val(val_f, idd, bfg55[ChN],  1, 1, i16);
        if (bfg55[ChN] < bitfg55) {bitfg55 = bfg55[ChN];}
        
        val_f = mfg41[ChN][2*i+k];
        bfg41[ChN] = bit_val(val_f, idd, bfg41[ChN],  1, 1, i16);
        if (bfg41[ChN] < bitfg41) {bitfg41 = bfg41[ChN];}
        val_f = mfg43[ChN][2*i+k];
        bfg43[ChN] = bit_val(val_f, idd, bfg43[ChN],  1, 1, i16);
        if (bfg43[ChN] < bitfg43) {bitfg43 = bfg43[ChN];}
        val_f = mfg45[ChN][2*i+k];
        bfg45[ChN] = bit_val(val_f, idd, bfg45[ChN],  1, 1, i16);
        if (bfg45[ChN] < bitfg45) {bitfg45 = bfg45[ChN];}
        

      } // for k                                                                                                                                                                                                                         



    }  // for i                                                                                                                                                                                
     


    if(bpf1[ChN]<bpf[ChN]){bpf[ChN]=bpf1[ChN];}
    if(bhf1[ChN]<bhf[ChN]){bhf[ChN]=bhf1[ChN];}

    if (bpf[ChN] < cpf) {printf("bpf out of range Ch=%d bpf=%d \n",ChN,bpf[ChN]);}
    if (bpf1[ChN] < cpf1) {printf("bpf1 out of range Ch=%d bpf1=%d \n",ChN,bpf1[ChN]);}
    if (bhf[ChN] < chf) {printf("bhf out of range Ch=%d bhf=%d \n",ChN,bhf[ChN]);}
    if (bhf1[ChN] < chf1) {printf("bhf1 out of range Ch=%d bhf1=%d \n",ChN,bhf1[ChN]);}
    if (bfg51[ChN] < cf51) {printf("bfg51 out of range Ch=%d bf=%d \n",ChN,bfg51[ChN]);}
    if (bfg52[ChN] < cf52) {printf("bfg52 out of range Ch=%d bf=%d \n",ChN,bfg52[ChN]);}
    if (bfg53[ChN] < cf53) {printf("bfg53 out of range Ch=%d bf=%d \n",ChN,bfg53[ChN]);}
    if (bfg54[ChN] < cf54) {printf("bfg54 out of range Ch=%d bf=%d \n",ChN,bfg54[ChN]);}
    if (bfg55[ChN] < cf55) {printf("bfg55 out of range Ch=%d bf=%d \n",ChN,bfg55[ChN]);}
    if (bfg41[ChN] < cf41) {printf("bfg41 out of range Ch=%d bf=%d \n",ChN,bfg41[ChN]);}
    if (bfg43[ChN] < cf43) {printf("bfg43 out of range Ch=%d bf=%d \n",ChN,bfg43[ChN]);}
    if (bfg45[ChN] < cf45) {printf("bfg45 out of range Ch=%d bf=%d \n",ChN,bfg45[ChN]);}


    if (bfg51[ChN] > bfg41[ChN]) {bfg51[ChN]=bfg41[ChN];}
    if (bfg53[ChN] > bfg43[ChN]) {bfg53[ChN]=bfg43[ChN];}
    if (bfg55[ChN] > bfg45[ChN]) {bfg55[ChN]=bfg45[ChN];}


    if (bpf[ChN] < cpf) {bpf[ChN]=cpf;}
    if (bhf[ChN] < chf) {bhf[ChN]=chf;}
   
    if (bfg51[ChN] < cf51) {bfg51[ChN]=cf51;}
    if (bfg52[ChN] < cf52) {bfg52[ChN]=cf52;}
    if (bfg53[ChN] < cf53) {bfg53[ChN]=cf53;}
    if (bfg54[ChN] < cf54) {bfg54[ChN]=cf54;}
    if (bfg55[ChN] < cf55) {bfg55[ChN]=cf55;}


  }//ChN


  short int  bpfc;
  short int  bhfc;
  short int  bf51;
  short int  bf52;
  short int  bf53;
  short int  bf54;
  short int  bf55;
  short int  bf41;
  short int  bf43;
  short int  bf45;

  bpfc = 35; 
  bhfc = 35; 
  bf51 = 35; 
  bf52 = 35; 
  bf53 = 35; 
  bf54 = 35; 
  bf55 = 35; 
  bf41 = 35; 
  bf43 = 35; 
  bf45 = 35; 

/*  bfc = 15; 
  bf31 = 16; 
  bf32 = 16; 
  bf33 = 19; 
  bf41 = 16; 
  bf43 = 19; */

  for (int chn = 0; chn < 16; chn++){

    if(bpf[chn] < bpfc){bpfc=bpf[chn];}
    if(bhf[chn] < bhfc){bhfc=bhf[chn];}
    if(bfg51[chn] < bf51){bf51=bfg51[chn];}
    if(bfg52[chn] < bf52){bf52=bfg52[chn];}
    if(bfg53[chn] < bf53){bf53=bfg53[chn];}
    if(bfg54[chn] < bf54){bf54=bfg54[chn];}
    if(bfg55[chn] < bf55){bf55=bfg55[chn];}
    if(bfg41[chn] < bf41){bf41=bfg41[chn];}
    if(bfg43[chn] < bf43){bf43=bfg43[chn];}
    if(bfg45[chn] < bf45){bf45=bfg45[chn];}

  }

  //   printf("  %d %d  %d  %d  %d  \n",bfc,bf31,bf32,bf33);


  for (int chn = 0; chn < 16; chn++){
    //            printf(" do %d %d  %d  %d  %d  \n",chn,bf[chn],bfg31[chn],bfg32[chn],bfg33[chn]);

    bpf[chn] = bpfc;
    bhf[chn] = bhfc;
    bfg51[chn] = bf51;
    bfg52[chn] = bf52;
    bfg53[chn] = bf53;
    bfg54[chn] = bf54;
    bfg55[chn] = bf55;
    bfg41[chn] = bf41;
    bfg43[chn] = bf43;
    bfg45[chn] = bf45;
    
 // printf("chn=%d, bpf=%d, bpf1=%d, bhf=%d, phf1=%d, bfg51=%d, bfg52=%d, bfg53=%d, bfg54=%d, bfg55=%d, bfg41=%d, bfg43=%d, bfg45=%d  \n",chn, bpf[chn], bpf1[chn], bhf[chn], bhf1[chn], bfg51[chn], bfg52[chn], bfg53[chn], bfg54[chn], bfg55[chn], bfg41[chn], bfg43[chn], bfg45[chn]);

  //printf("chn=%d, mpf=%lf, mpf1=%lf, mhf=%lf, mhf1=%lf, mfg51=%lf, mfg52=%lf, mfg53=%lf, mfg54=%lf, mfg55=%lf, mfg41=%lf, mfg43=%lf, mfg45=%lf  \n",chn, mpf[chn][0], mpf1[chn][0], mhf[chn][0], mhf1[chn][0], mfg51[chn][0], mfg52[chn][0], mfg53[chn][0], mfg54[chn][0], mfg55[chn][0], mfg41[chn][0], mfg43[chn][0], mfg45[chn][0]);

  }

  printf("bpfc=%d, bhfc=%d, bf51=%d, bf52=%d, bf53=%d, bf54=%d, bf55=%d, bf41=%d, bf43=%d, bf45=%d  \n", bpfc, bhfc, bf51, bf52, bf53, bf54, bf55, bf41, bf43, bf45);




}


void eclFgCal::CalculateFginInt( ){

  //  int ChN;
  //  int k;
  //  int i;

  long int iff;
  long int ihf;
  long int ia;
  long int ib;
  long int ic;
  long int id;
  long int ie;

  long int imax1;

  imax1=myPow(2,31);
  double imax;
  imax=(double)imax1;


  //  for (ChN = 9; ChN < 10; ChN++) {


  for (int ChN = 0; ChN < 16; ChN++) {

    //if (!mask_[ChN]) continue;
    //    printf("calc %d %d %d %d  \n",ChN,bf[ChN],bfg31[ChN],bfg32[ChN],bfg33[ChN]);
   
    iff = (long int)1 << bpf[ChN];
    ihf = (long int)1 << bhf[ChN];
    ia = (long int)1 << bfg51[ChN];
    ib = (long int)1 << bfg52[ChN];
    ic = (long int)1 << bfg53[ChN];
    id = (long int)1 << bfg54[ChN];
    ie = (long int)1 << bfg55[ChN];
    
   // printf("")
    
  
    
       //printf("%ld %ld %ld %ld %ld %ld %lf  \n",iff,ia,ib,ic,id,ie,imax);
    
       //printf("%d \n",bfg55[ChN]);
    for (int k = 0; k < 192; k++) {
     
     /*if(ChN==15){

     printf("\n");
      printf("%4d",k);
       for (int i = 0; i < 16; i++) {
         printf("%11lf",fg52[ChN][k][i]);
       }
        printf("\n");
     }*/

      const double isd = 3. / 4., sd = 1 / isd ; // conversion factor (???) 
      
     // long long int sum =0;
     // double dsum =0;
   // if(ChN==15)  printf("%4d",k);
      for (int i = 0; i < 16; i++) {
	double w = i ? 1.0 : 1. / 16.;
	

	m_pf[ChN][k][i] = lrint(p_f[ChN][k][i]  * iff);
 // printf("chn=%d k=%d %lf\n", ChN, k, p_f[ChN][k][0]);
	m_pf1[ChN][k][i] = lrint(p_f1[ChN][k][i] * iff * w * sd);
 
  m_hf[ChN][k][i] = lrint(h_f[ChN][k][i]  * ihf);
	m_hf1[ChN][k][i] = lrint(h_f1[ChN][k][i] * ihf * w * sd);
	
  /*long long int f51 = fg51[ChN][k][i] * ia;
  if(i==0) f51=f51>>4;
  m_fg51[ChN][k][i] = lrint(f51);*/
	m_fg51[ChN][k][i] = lrint(fg51[ChN][k][i] * ia * w);
	m_fg52[ChN][k][i] = lrint(fg52[ChN][k][i] * ib * w * isd);
//if(ChN==15) printf("%11d",m_fg52[ChN][k][i]);

  m_fg53[ChN][k][i] = lrint(fg53[ChN][k][i] * ic * w);
	m_fg54[ChN][k][i] = lrint(fg54[ChN][k][i] * id * w * isd);
	m_fg55[ChN][k][i] = lrint(fg55[ChN][k][i] * ie * w);

 // long long int f52=m_fg52[ChN][k][i];
  //if(i==0) f52=f52<<4;

  //sum+=	f52*m_pf[ChN][k][i]*1000;
  //dsum+=	fg52[ChN][k][i]*p_f[ChN][k][i];//*h_f[ChN][k][i];
	
  //printf("%d, %d  %d, %d, %d\n", 	m_fg51[ChN][k][i], 	m_fg52[ChN][k][i],	m_fg53[ChN][k][i],	m_fg54[ChN][k][i],	m_fg55[ChN][k][i]);
	

	
	if(p_f[ChN][k][i]  * iff * w>imax || p_f[ChN][k][i]  * iff * w<-imax ){
	  printf("alarm ChN=%d k=%d i=%d m_pf=%lf imax=%lf \n",ChN,k,i,p_f[ChN][k][i]  * iff * w,imax);}        
	if(p_f1[ChN][k][i] * iff * w * sd>imax || p_f1[ChN][k][i] * iff * w * sd<-imax ){
	  printf("alarm ChN=%d k=%d i=%d m_pf1=%lf imax=%lf \n",ChN,k,i,p_f1[ChN][k][i] * iff * w * sd,imax);}        
  if(h_f[ChN][k][i]  * ihf * w>imax || h_f[ChN][k][i]  * ihf * w<-imax ){
	  printf("alarm ChN=%d k=%d i=%d m_hf=%lf imax=%lf \n",ChN,k,i,h_f[ChN][k][i]  * ihf * w,imax);}        
	if(h_f1[ChN][k][i] * ihf * w * sd>imax || h_f1[ChN][k][i] * ihf * w * sd<-imax ){
	  printf("alarm ChN=%d k=%d i=%d m_hf1=%lf imax=%lf \n",ChN,k,i,h_f1[ChN][k][i] * ihf * w * sd,imax);}        
	if(fg51[ChN][k][i] * ia * w>imax || fg51[ChN][k][i] * ia * w<-imax ){
	  printf("alarm ChN=%d k=%d i=%d m_fg51=%lf imax=%lf \n",ChN,k,i,fg51[ChN][k][i] * ia * w,imax);}        
	if(fg52[ChN][k][i] * ib * w * isd>imax || fg52[ChN][k][i] * ib * w * isd<-imax ){
	  printf("alarm ChN=%d k=%d i=%d m_fg52=%lf imax=%lf \n",ChN,k,i,fg52[ChN][k][i] * ib * w * isd,imax);}        
	if(fg53[ChN][k][i] * ic * w>imax || fg53[ChN][k][i] * ic * w<-imax ){
	  printf("alarm ChN=%d k=%d i=%d m_fg53=%lf imax=%lf \n",ChN,k,i,fg53[ChN][k][i] * ic * w,imax);}        
	if(fg54[ChN][k][i] * id * w * isd>imax || fg54[ChN][k][i] * id * w * isd<-imax ){
	  printf("alarm ChN=%d k=%d i=%d m_fg54=%lf imax=%lf \n",ChN,k,i,fg54[ChN][k][i] * id * w * isd,imax);}        
	if(fg55[ChN][k][i] * ie * w>imax || fg55[ChN][k][i] * ie * w<-imax ){
	  printf("alarm ChN=%d k=%d i=%d m_fg55=%lf imax=%lf \n",ChN,k,i,fg55[ChN][k][i] * ie * w,imax);}        

  }

 
//	sum=sum/ib;
  //double ddsum=sum;
  //ddsum=ddsum/iff;
  //printf("%lf %lf \n", ddsum, dsum);


	
	int jk = 23 + ((48 - k ) >> 2);
	if (jk >= 0 && jk < 24 && (48 - k) % 4 == 0) {
	//  long long int sum =0;
   // double dsum =0;
	  for (int i1 = 0; i1 < 16; i1++) {
	    double w1 = i1 ? 1.0 : 1. / 16.;
	
	    m_fg41[ChN][jk][i1] = lrint(fg41[ChN][jk][i1] * ia * w1);
	    m_fg43[ChN][jk][i1] = lrint(fg43[ChN][jk][i1] * ic * w1);
      m_fg45[ChN][jk][i1] = lrint(fg45[ChN][jk][i1] * ie * w1);
      //sum+=	m_fg43[ChN][jk][i1]/w1*m_hf[ChN][k][i1]*100;
    //  dsum+=	fg41[ChN][jk][i1]*h_f[ChN][jk][i1];
   //   printf("%d, %d  %d\n", 	m_fg41[ChN][jk][i1], 	m_fg43[ChN][jk][i1],	m_fg45[ChN][jk][i1]);
	

	    if(fg41[ChN][jk][i1] * ia * w1>imax || fg41[ChN][jk][i1] * ia * w1<-imax ){
	      printf("alarm ChN=%d k=%d i=%d m_fg41=%lf imax=%lf \n",ChN,k,i1,fg41[ChN][jk][i1] * ia * w1,imax);}        
	    if(fg43[ChN][jk][i1] * ic * w1>imax || fg43[ChN][jk][i1] * ic * w1<-imax ){
	      printf("alarm ChN=%d k=%d i=%d m_fg43=%lf imax=%lf \n",ChN,k,i1,fg43[ChN][jk][i1] * ic * w1,imax);}        
	    if(fg45[ChN][jk][i1] * ie * w1>imax || fg45[ChN][jk][i1] * ie * w1<-imax ){
	      printf("alarm ChN=%d k=%d i=%d m_fg45=%lf imax=%lf \n",ChN,k,i1,fg45[ChN][jk][i1] * ie * w1,imax);}        
	    
	    
	   //   printf("!!! ch=%d fg41=%lf fg43=%lf fg45=%lf m_fg41=%d m_fg43=%d  m_fg45=%d ia=%ld ic=%ld ie=%ld\n",ChN,fg41[ChN][jk][i1],fg43[ChN][jk][i1], fg45[ChN][jk][i1],m_fg41[ChN][jk][i1],m_fg43[ChN][jk][i1],m_fg45[ChN][jk][i1],ia,ic,ie);
	  }

  //sum=sum/ic;
 //double ddsum=sum;
  //ddsum=ddsum/ihf;
 // printf("%lf\n", ddsum);
	}

      
    }    
    
    
  }


}



void eclFgCal::WriteCalib( const char *dir){

  FILE *fr;
  int h,nsiz,nsiz1; //,j,i;
 
  char fname[256];
  char tname[12];
 
  sprintf(tname, "ECLDSP FILE");

  short int id[256];

  for (int i = 0; i < 256; i++){

    if(i<6){   id[i]=(short int)*(tname+2*i)+256*(short int)*(tname+2*i+1);}
    else{id[i]=0;}
    //	printf("%d %d \n",i,atoi(*(tname+i)));
    //     id=atoi(tname);
  }

  //  id[0]=17221;

  short int  bpfc;
  short int  bhfc;
  short int  bf51;
  short int  bf52;
  short int  bf53;
  short int  bf54;
  short int  bf55;
  short int  bf41;
  short int  bf43;
  short int  bf45;

  bpfc = 35; 
  bhfc = 35; 
  bf51 = 35; 
  bf52 = 35; 
  bf53 = 35; 
  bf54 = 35; 
  bf55 = 35; 
  bf41 = 35; 
  bf43 = 35; 
  bf45 = 35; 


  // fill thresholds                                                                                                                                                                    
  for (int chn = 0; chn < 16; chn++){
    id[128 + chn] = 128 + 40;  // A_low
    id[192 + chn] = 128 - 30;  // A_skip
    id[64  + chn] =  -100;  // A_hard

    if(bpf[chn] < bpfc){bpfc=bpf[chn];}
    if(bhf[chn] < bhfc){bhfc=bhf[chn];}
    if(bfg51[chn] < bf51){bf51=bfg51[chn];}
    if(bfg52[chn] < bf52){bf52=bfg52[chn];}
    if(bfg53[chn] < bf53){bf53=bfg53[chn];}
    if(bfg54[chn] < bf54){bf54=bfg54[chn];}
    if(bfg55[chn] < bf55){bf55=bfg55[chn];}
    if(bfg41[chn] < bf41){bf41=bfg41[chn];}
    if(bfg43[chn] < bf43){bf43=bfg43[chn];}
    if(bfg45[chn] < bf45){bf45=bfg45[chn];}


    
  }

  
  
     
  id[13] = bf51+256*bf52;  // ka kb     
  id[14] =  bf53+256*bf54; // kc kd                               
  id[15] =  bf55+0*256;  // ke  y0                                   
  //  id[15] =  chi+thres;                                  

  id[16] =  30000;                                  
  //  id[16] = 14+k2_chi*256;   //k1_chi k2_chi
  id[17] = bpfc+10*256;   //k1_chi k2_chi
  id[8]=257;
  id[8]=513;  // new version
  id[9]=-1;

  for (int i = 18; i < 64; i++){
    id[i]=0;
  }

  //id[18] = 0;

 

  //for (int i = 0; i < 256; i++){
    //printf("%d ",id[i]);
    //    if(i%16==0){printf(" Q %d\n",i/16);}

  //}

  sprintf(fname,"%s/%s%s",dir,"dsp",".dat");
  printf("eclFgCal::WriteCalib: %s \n",fname);

  fr= fopen(fname, "wb");

  if(fr == NULL){
    std::cout << " can't open  file "
	      << fname << "\n";

    exit(0);


  }
  nsiz=2;
  nsiz1=256;

  h = fwrite(id, nsiz, nsiz1, fr);
  if (h != nsiz1) {
    std::cout << "Error writing id data. Read block size = " << h << '\n';
    exit(0);
  }


  nsiz=4;
  for(int i=0;i<16;i++)
    {

      printf("eclFgCal::WriteCalib: i=%d\n",i);
      nsiz1=384;
      h = fwrite(m_fg41[i][0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
	      std::cout << "Error writing fg41 data. Read block size = " << h << '\n';
        exit(0);
      }

    
      h = fwrite(m_fg43[i][0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
	      std::cout << "Error writing fg43 data. Read block size = " << h << '\n';
        exit(0);
      }
      

      nsiz1=3072;
      h = fwrite(m_fg51[i][0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
	      std::cout << "Error writing fg51 data. Read block size = " << h << '\n';
        exit(0);
      }


      h = fwrite(m_fg52[i][0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
      	std::cout << "Error writing fg52 data. Read block size = " << h << '\n';
        exit(0);
      }

      h = fwrite(m_fg53[i][0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
	      std::cout << "Error writing fg53 data. Read block size = " << h << '\n';
        exit(0);
      }

      h = fwrite(m_fg54[i][0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
	      std::cout << "Error writing fg54 data. Read block size = " << h << '\n';
        exit(0);
      }

      h = fwrite(m_fg55[i][0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
	      std::cout << "Error writing fg55 data. Read block size = " << h << '\n';
        exit(0);
      }

      nsiz1=384;
      h = fwrite(m_fg45[i][0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
	      std::cout << "Error writing fg45 data. Read block size = " << h << '\n';
        exit(0);
      }


      nsiz1=3072;
      h = fwrite(m_pf[i][0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
        std::cout << "Error writing pf data. Read block size = " << h << '\n';
        exit(0);
      }

      h = fwrite(m_pf1[i][0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
	      std::cout << "Error writing pf1 data. Read block size = " << h << '\n';
        exit(0);
      }

      h = fwrite(m_hf[i][0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
        std::cout << "Error writing hf data. Read block size = " << h << '\n';
        exit(0);
      }

      h = fwrite(m_hf1[i][0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
	      std::cout << "Error writing hf1 data. Read block size = " << h << '\n';
        exit(0);
      }

      printf("eclFgCal::WriteCalib1: i=%d\n",i);
    }

  printf("%d", m_fg55[15][23][1]);
  fclose(fr);
}


#ifdef notdef
void eclFgCal::ReadCalib( const char *dir ){


  // FILE *fr;
  //FILE *zr;
  //int h,nsiz,nsiz1,j,i;
 
  //  char fname[256];
  //  char zname[256];
 
  short int id[256];
  short int cf[192][16];
  short int cf1[192][16];
  short int cfg31[192][16];
  short int cfg32[192][16];
  short int cfg33[192][16];
  short int cfg41[24][16];
  short int cfg43[24][16];



  short int idq[256];
  short int cfq[192][16];
  short int cfq1[192][16];
  short int cfgq31[192][16];
  short int cfgq32[192][16];
  short int cfgq33[192][16];
  short int cfgq41[24][16];
  short int cfgq43[24][16];

  int y;

  for ( int u =0 ; u < 16; u ++){
  for ( int k =0 ; k < 192; k ++){

    if(k < 24){
      cfg41[k][u] = 0;
      cfg43[k][u] = 0;
      cfgq41[k][u] = 0;
      cfgq43[k][u] = 0;
    }
   
    cfg31[k][u] = 0;
    cfg32[k][u] = 0;
    cfg33[k][u] = 0;
    cf[k][u]    = 0;
    cf1[k][u]   = 0;

    cfgq31[k][u] = 0;
    cfgq32[k][u] = 0;
    cfgq33[k][u] = 0;
    cfq[k][u]    = 0;
    cfq1[k][u]   = 0;

  }
  }
  /*

  for(y= 1; y <17; y++){
    printf("y=%d \n",y);
     sprintf(fname,"%s/%d/%s_%d%s",dir,cr,"dsp",cr,".dat");
     //     sprintf(zname,"%s/%d/%s_%d%s",dir,cr,"dsp",cr,".dat");
               sprintf(zname,"%s/%d/%s0%d%s",dir,cr,"dsp",1,".dat");
     
     //     sprintf(fname,"%s/%d/%s0%d%s",dir,cr,"dsp",1,".dat");
  //    printf(" read %s \n",fname);

  fr= fopen(fname, "rb");

  if(fr == NULL){
    std::cout << " can't open  file "
	      << fname << "\n";

    exit(0);


  }



  nsiz=2;
  nsiz1=256;

  h = fread(id, nsiz, nsiz1, fr);
  if (h != nsiz1) {
    std::cout << "Error writing id data. Read block size = " << h << '\n';

    exit(0);

  }

  printf("%d %d %d \n ",id[26],id[27],id[28]);



  

  //  for(i=0;i<16;i++){
  for(i=0;i<y;i++){



      nsiz1=384;
      h = fread(cfg41[0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
	std::cout << "Error writing fg41 data. Read block size = " << h << '\n';
      exit(0);

      }





      nsiz1=3072;
            h = fread(cfg31[0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
	std::cout << "Error writing fg31 data. Read block size = " << h << '\n';

      exit(0);

      }





        h = fread(cfg32[0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
	std::cout << "Error writing fg32 data. Read block size = " << h << '\n';

    exit(0);

      }



        h = fread(cfg33[0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
	std::cout << "Error writing fg33 data. Read block size = " << h << '\n';
    exit(0);

      }




      nsiz1=384;
      h = fread(cfg43[0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
	std::cout << "Error writing fg43 data. Read block size = " << h << '\n';
    exit(0);

      }


      nsiz1=3072;
            h = fread(cf[0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
	std::cout << "Error writing f data. Read block size = " << h << '\n';
    exit(0);

      }




        h = fread(cf1[0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
	std::cout << "Error writing f1 data. Read block size = " << h << '\n';
    exit(0);

      }


 

    
    }




  fclose(fr);

  
  zr= fopen(zname, "rb");

  if(zr == NULL){
    std::cout << " can't open  file "
	      << zname << "\n";

    exit(0);


  }

  nsiz=2;
  nsiz1=256;

  h = fread(idq, nsiz, nsiz1, zr);


  //  for(i=0;i<16;i++) {

  for(i=0;i<y;i++){

      nsiz1=384;

      h = fread(cfgq41[0], nsiz, nsiz1, zr);
      if (h != nsiz1) {
	std::cout << "Error writing fg41 data. Read block size = " << h << '\n';
      exit(0);

      }



      nsiz1=3072;


          h = fread(cfgq31[0], nsiz, nsiz1, zr);
      if (h != nsiz1) {
	std::cout << "Error writing fg31 data. Read block size = " << h << '\n';

      exit(0);

      }




        h = fread(cfgq32[0], nsiz, nsiz1, zr);
      if (h != nsiz1) {
	std::cout << "Error writing fg32 data. Read block size = " << h << '\n';

    exit(0);

      }




        h = fread(cfgq33[0], nsiz, nsiz1, zr);
      if (h != nsiz1) {
	std::cout << "Error writing fg33 data. Read block size = " << h << '\n';
    exit(0);

      }


      nsiz1=384;

      h = fread(cfgq43[0], nsiz, nsiz1, zr);
      if (h != nsiz1) {
	std::cout << "Error writing fg41 data. Read block size = " << h << '\n';
      exit(0);

      }

      nsiz1=3072;


           h = fread(cfq[0], nsiz, nsiz1, zr);
      if (h != nsiz1) {
	std::cout << "Error writing f data. Read block size = " << h << '\n';
    exit(0);

      }


        h = fread(cfq1[0], nsiz, nsiz1, zr);
      if (h != nsiz1) {
	std::cout << "Error writing f1 data. Read block size = " << h << '\n';
    exit(0);

      }

    }
  fclose(zr);



   



        
  for (int  k =0 ; k < 192; k ++){
            
      
  for ( j =0 ; j < 16; j ++){

    if(cf[k][j]-cfq[k][j]!=0){printf("f k=%d j=%d d=%d f=%d %d \n", k,j,cf[k][j]-cfq[k][j],cfq[k][j],cf[k][j]);}

    if(cf1[k][j]-cfq1[k][j]!=0){printf("f1 k=%d j=%d d=%d f=%d %d \n", k,j,cf1[k][j]-cfq1[k][j],cfq1[k][j],cf1[k][j]);}


    if(cfg31[k][j]-cfgq31[k][j]!=0){printf("fg31 k=%d j=%d d=%d f=%d %d \n", k,j,cfg31[k][j]-cfgq31[k][j],cfgq31[k][j],cfg31[k][j]);}


    if(cfg32[k][j]-cfgq32[k][j]!=0){printf("fg32 k=%d j=%d d=WriteC%d f=%d %d \n", k,j,cfg41[k][j]-cfgq41[k][j],cfgq41[k][j],cfg41[k][j]);}




    if(cfg43[k][j]-cfgq43[k][j]!=0){printf("fg43 k=%d j=%d d=%d f=%d %d \n", k,j,cfg43[k][j]-cfgq43[k][j],cfgq43[k][j],cfg43[k][j]);}

    }


    
  
    
  }
  }

      
  }    

*/    




 
}

#endif


void eclFgCal::RootWr( ){

  wr_rootf=1;



}
/*

void eclFgCal::RootFill( TFile *rootfile ){
 
  //TDirectory *tdir_crate;
  TDirectory *tdir_shaper;
  TDirectory *tdir_channel;

  TString tdirname;



  if( wr_rootf == 1){


  tdirname.Form("shaper%d",cr);




  tdir_shaper = rootfile->mkdir(tdirname);



  tdir_shaper->cd();


    Kf->Write();
  Ka->Write();
  Kb->Write();
  Kc->Write();
   

	
   for( int icn=0; icn < 16; icn ++){



    tdirname.Form("chan%d",icn);
    tdir_channel = tdir_shaper->mkdir(tdirname);
    tdir_channel->cd();




  for( int i=0; i < 192; i ++){


    
    
    uf[i][icn]->Write();

    uf1[i][icn]->Write();
   
    uf31[i][icn]->Write();
    uf32[i][icn]->Write();
    uf33[i][icn]->Write();
        
    if(i<24){
    uf41[i][icn]->Write();
    uf43[i][icn]->Write();
     }
    	
     }

     }
  }
}



void eclFgCal::RootFuncFill( TFile *rootfile ){
 


  TTree *m_tree;

  
  int           t_sp;
  int           t_nevt;

  double        t_f[192][16];                                                                                                                                                   
  double        t_f1[192][16];                                                                               
  m_tree  = new TTree("function_shape_tree", "function_shape__tree");


  m_tree->Branch("sp",      &t_sp,      "sp/I");
  m_tree->Branch("nevt",       &t_nevt,       "nevt/I");
  m_tree->Branch("f",   t_f,      "f[192][16]/D");
  m_tree->Branch("f1",     t_f1,      "f1[192][16]/D");


  t_nevt=0;
  t_sp=cr;

  for(int i=0; i<16; i++){


 

  for(int k=0; k<192; k++){
  for(int j=0; j<16; j++){

    t_f[k][j]=f[t_nevt][k][j];
    t_f1[k][j]=f1[t_nevt][k][j];

  }

  }

    m_tree->Fill();
    t_nevt++;
  }

  m_tree->Write();
     
 rootfile->Close();


  
 }*/



/*void eclFgCal::RootFuncRead( const char *treeroot ){
 
  TChain fChain("function_shape_tree");

  fChain.Add(treeroot);


  int           t_sp;
  int           t_nevt;

  double        t_f[192][16];                                                                                                                                                   
  double        t_f1[192][16];                                                 
  
  TBranch        *b_sp;   //!
  TBranch        *b_nevt;   //!
  TBranch        *b_f;   //!
  TBranch        *b_f1;   //!



   fChain.SetBranchAddress("sp", &t_sp, &b_sp);
   fChain.SetBranchAddress("nevt", &t_nevt, &b_nevt);
   fChain.SetBranchAddress("f", &t_f, &b_f);
   fChain.SetBranchAddress("f1", &t_f1, &b_f1);





  Int_t nevent = fChain.GetEntries();

  for (Int_t i=0;i<nevent;i++) {


    fChain.GetEntry(i);
    if(i<16){

  for(int k=0; k<192; k++){
  for(int j=0; j<16; j++){

    f[i][k][j]=t_f[k][j];
    f1[i][k][j]=t_f1[k][j];

  }

  }
    

    }

  }  //


  //  fChain.close();
  
  }*/




void eclFgCal::GetResponsePar( const char *fname){

  FILE *fr;
  double a0,a1,a2,a3,a4,a5,a6,a7,a8,a9;
  int u=0;
  int uold;
  fr= fopen(fname, "r");

  if(fr == NULL){
    std::cout << " can't open  file "
      << fname << "\n";
    exit(0);
  }

  fscanf (fr, "%lf   ", &a0);

  while(!feof(fr)){
    uold=u;
    fscanf (fr, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf    ", &u,&a0,&a1,&a2,&a3,&a4,&a5,&a6,&a7,&a8,&a9);
    if(u>0){
      if(u!=uold){
        p_par[u-1][0]=a0+0.25;
        //      par[u-1][0]=a0+0.25-3./192.;
        //      par[u-1][0]=a0;
        p_par[u-1][1]=a1;
        p_par[u-1][2]=a2;
        p_par[u-1][3]=a3;
        p_par[u-1][4]=a4;
        p_par[u-1][5]=a5;
        p_par[u-1][6]=a6;
        p_par[u-1][7]=a7;
        p_par[u-1][8]=a8;
        p_par[u-1][9]=a9;
      }
      else{
        h_par[u-1][0]=a0+0.25;
        //      par[u-1][0]=a0+0.25-3./192.;
        //      par[u-1][0]=a0;
        h_par[u-1][1]=a1;
        h_par[u-1][2]=a2;
        h_par[u-1][3]=a3;
        h_par[u-1][4]=a4;
        h_par[u-1][5]=a5;
        h_par[u-1][6]=a6;
        h_par[u-1][7]=a7;
        h_par[u-1][8]=a8;
        h_par[u-1][9]=a9;
      }
    }
  }

  fclose(fr);
}


void eclFgCal::ReadChiCut( const char *dir  ){


  FILE *fr;

  char fname[256]; 
  int a;
  int b;
    sprintf(fname,"%s/%d/%s%d%s",dir,cr,"wics_",cr,".txt");

    //    sprintf(fname,"%s/%d/%s%d%s","/home/belle/avbobrov/eclshaper_sp/file/tid_x",cr,"wics_",cr,".txt");

  fr= fopen(fname, "r");

  if(fr == NULL){
    std::cout << " can't open  file "
	      << fname << "\n";

    exit(0);


  }




  while(!feof(fr)){
    fscanf (fr, "%d %d   ", &a,&b);

    
    
  }

  fclose(fr);

  k2_chi=b;
  chi_thres=a;


}

void eclFgCal::SaveInverseMatrices(const char *dir){

  char fname[256];

  int icn;
  int i;
  int j;

  double IM[16][16];

  FILE *BMcoIN;

  
  for(icn =0; icn <16; icn++){
    
    for(j =0; j <16; j++){
      for(i =0; i <16; i++){
	      IM[j][i] = IY[icn][j][i];
      }
    }
    

    sprintf(fname,"%s/%s%d%s",dir,"Binmcor_",icn,".dat");
    //    printf("eclFgCal::SaveInverseMatrices %s \n",dir);
    
    printf("eclFgCal::SaveInverseMatrices %s %d \n",fname,icn);
    

    if ( (BMcoIN = fopen(fname, "wb")) == NULL)
      {
        fprintf(stderr, "Error opening file %s",fname);
        exit(1);
      }

    if (fwrite(IM, sizeof(double), 256, BMcoIN) != 256)
      {
        fprintf(stderr, "Error writing to file %s.",fname);
        exit(1);
      }

    fclose(BMcoIN);


  }//endfor icn

  for(icn =0; icn <16; icn++){
    
    for(j =0; j <16; j++){
      for(i =0; i <16; i++){
	      IM[j][i] = IY[icn][j][i];
      }
    }
    
    //    sprintf(fname,"%s/%d/%s_%d%s",dir,cr,"inmcor",icn,".dat");
    sprintf(fname,"%s/%s%d%s",dir,"inmcor_",icn,".dat");



    if ( (BMcoIN = fopen(fname, "w")) == NULL)
      {
        fprintf(stderr, "Error opening file %s",fname);
        exit(1);
      }

    for(i=0;i<16;i++){
      fprintf(BMcoIN,                                                                                                        
      "%.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e\
  \t %.5e \n  ",                                                                                                                  
        IM[0][i], IM[1][i], IM[2][i], IM[3][i], IM[4][i], IM[5][i],IM[6][i], IM[7][i], IM[8][i], IM[9][i], IM[10][i],                                                                                                                IM[11][i], IM[12][i], IM[13][i], IM[14][i], IM[15][i]); 

    }

    fclose(BMcoIN);




  }


}
