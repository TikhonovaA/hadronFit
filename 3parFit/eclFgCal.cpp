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
    
    bf[ii]=15;
    bf1[ii]=15;
    bfg31[ii]=16;
    bfg32[ii]=16;
    bfg33[ii]=19;
    bfg41[ii]=16;
    bfg43[ii]=19;
    

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
      mf[ii][jj]=0.;
      mf1[ii][jj]=0.;
      mfg31[ii][jj]=0.;
      mfg32[ii][jj]=0.;
      mfg33[ii][jj]=0.;
      
      mfg41[ii][jj]=0.;
      mfg43[ii][jj]=0.;
    }
    for ( int kk = 0; kk < 192; kk++){
      if( kk <16){
        IY[kk][ii][jj]=0.;
      }    
      f[ii][kk][jj]=0.;
      f1[ii][kk][jj]=0.;
      fg31[ii][kk][jj]=0.;
      fg32[ii][kk][jj]=0.;
      fg33[ii][kk][jj]=0.;
      m_f[ii][kk][jj]=0.;
      m_f1[ii][kk][jj]=0.;
      m_fg31[ii][kk][jj]=0;
      m_fg32[ii][kk][jj]=0;
      m_fg33[ii][kk][jj]=0;
      
      if( kk < 24 ){
        fg41[ii][kk][jj]=0.;
        fg43[ii][kk][jj]=0.;
        m_fg41[ii][kk][jj]=0;
        m_fg43[ii][kk][jj]=0;
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






int myPow(int x, int p)
{
  if (p == 0) return 1;
  if (p == 1) return x;
  return x * myPow(x, p - 1);
}


short  int bit_val(double f, int idd, short int bit, int a, int b, int lim)
{
  short int i;
  i = bit;

  int val;
  int R;
  int flag;
  flag = 0;
  if (i < 11) {i = 11;}
  do {
    R = myPow(2, i);
    val = (int)(a * R * f / idd / b + R + 0.5) - R;

    if ((val < lim && val > -lim)) {
      flag = 1;
    } else {
      i = i - 1;

    }

  } while (i > 8 && flag == 0);


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

  //    printf("svp=  %e %e20.15 %e20.15 %e %e %e %e %e %e %e %e \n ",s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9],s[10],s[11]);

  // printf("svp=  %21.18lf %21.18lf  %21.18lf %21.18lf %21.18lf %21.18lf %21.18lf %21.18lf %21.18lf %21.18lf %21.18lf  \n ",s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9],s[10],s[11]);

    // printf("svp=  %x %x %x %x %x %x %x %x %x %x %x \n ",s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9],s[10],s[11]);
                                                                                     
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

  double ssssj, ssssj1, ssssi, ssssi1;
  double svp, svm;
  double t[1], tmd[1], tpd[1];

  /*  double MP[10] = {0.5, 0.6483, 0.4017, 0.3741, 0.8494, 0.00144547, 4.7071, 0.8156, 0.5556, 0.2752};*/

  double MZ[11] = {0.75, 0.6483, 0.4017, 0.3741, 0.8494, 0.00144547, 4.7071, 0.8156, 0.5556, 0.2752, 27.7221};
  double spt=0;
  double fm;
  double fc;
 
  // varible for one channel                                                                                                                                                                  

  double  g1g1[192], gg[192], gg1[192]; //, dgg[192];
  double  sg1[16][192], sg[16][192], sg2[16][192], gg2[192], g1g2[192], g2g2[192];
  double  dgg1[192], dgg2[192];
  int chan;

  double scale =  144./127.;

  for ( chan = 0; chan < 16; chan++){

    /*
    for (j1 = 0; j1 < 10; j1++) {
      MZ[j1]=par[chan][j1]; //par for one chan
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
    MZ[10]=27.7221;
    
    //MZ[11]=27.7221;
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

     //   printf("spt=%e \n",spt);

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

   // printf("fm=%lf \n",fm);

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

    //printf("FM=%15.12e \n",fm);

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

      // set but not used!                                                                                                                                                                    
      //tin = t + dt;                                                     

      //    printf("%t=%21.18e t0=%21.18e dt=%21.18e dt0=%21.18e ndt=%21.18e del=%21.18e \n ",t,t0,dt,dt0,ndt,del);                                   
      //    printf("j1=%d t0=%lf \n",j1,t+dt);
      for (int j = 0; j < 16; j++) {
        t[0] = t[0] + dt;

        if (t > 0) {
          f[chan][j1][j] = ShaperDSP_F(t, MZ);
          //	f[chan][j1][j] = ShaperDSP_F(t, MZ);
            //printf("j1=%d j=%d dt=%21.18e f=%21.18e t=%21.18e \n",j1,j,dt,f[chan][j1][j],t);
          ddt = 0.005 * dt;
          tc1 = t[0] - ts0;
          tpd[0] = t[0] + ddt;
          tmd[0] = t[0] - ddt;

        
          if (tc1 > ddt) {
            svp = ShaperDSP_F(tpd, MZ);
            svm = ShaperDSP_F(tmd, MZ);
            //	  svp = ShaperDSP_F(tpd, MZ);
            //	  svm = ShaperDSP_F(tmd, MZ);

            f1[chan][j1][j] = (svp - svm) / 2. / ddt;
            //	  printf("j1=%d j=%d f1=%lf \n",j1,j,f1[chan][j1][j]);

          } 
          else {
            if (tc1 > 0 ) {
              f1[chan][j1][j] = ShaperDSP_F(tpd, MZ) / (ddt + tc1);
                //f1[chan][j1][j] = ShaperDSP_F(tpd, MZ) / (ddt + tc1);
            }
            else {
              f1[chan][j1][j] = 0.;
            }
          }// else tc1>ddt                                                                                                                                                                    

        } //if t>0                                                                                                                                                                            
        else {
          f[chan][j1][j] = 0.;
          f1[chan][j1][j] = 0.;
        }

        //      printf("chan=%d j1=%d j=%d f=%e t=%e \n",chan,j1,j,f[chan][j1][j],f1[chan][j1][j]);

        //      if( chan == 15 ) cout << " cr="<<cr<< " ChN="<< chan <<" Gem k = " << j1<< " i=" <<j<< " f=" << f[chan][j1][j]<< " f1=" << f1[chan][j1][j] << endl;

        //      printf("CalculateFGPar p311\n");
        //      if (!(uf[j1][chan])) printf("No uf for %d %d\n",j1, chan);
        //printf("CalculateFGPar uf name=%s\n",uf[j1][chan]->GetName());
        //      uf[j1][chan]->SetBinContent(j+1,f[chan][j1][j]);
        //      printf("CalculateFGPar p312\n");
        //      uf1[j1][chan]->SetBinContent(j+1,f1[chan][j1][j]);
        //      printf("CalculateFGPar p313\n");

      } //for j                                                                                                                                                                               
      //    if(chan==15){printf("Kluch j=%d i=%d f=%lf f1=%lf \n",j1,j,f[chan][j1][j],f1[chan][j1][j]);}

    }//for j1    



    for (j1 = 0; j1 < endc; j1++) { //from 0 to 191

      gg[j1] = 0.;
      gg1[j1] = 0.;
      g1g1[j1] = 0.;
      gg2[j1] = 0.;
      g1g2[j1] = 0.;
      g2g2[j1] = 0.;
      for (int j = 0; j < 16; j++) {
        sg[j][j1] = 0.;
        sg1[j][j1] = 0.;
        sg2[j][j1] = 0.;

        ssssj1 = f1[chan][j1][j];
        ssssj = f[chan][j1][j];

        if(j==0){

          if(f[chan][j1][j]>mf[chan][0]){mf[chan][0]=f[chan][j1][j];}
          if(f[chan][j1][j]<mf[chan][1]){mf[chan][1]=f[chan][j1][j];}
          if(f1[chan][j1][j]>mf1[chan][0]){mf1[chan][0]=f1[chan][j1][j];}
          if(f1[chan][j1][j]<mf1[chan][1]){mf1[chan][1]=f1[chan][j1][j];}

        }

        else{

          if(f[chan][j1][j]>mf[chan][2]){mf[chan][2]=f[chan][j1][j];}
          if(f[chan][j1][j]<mf[chan][3]){mf[chan][3]=f[chan][j1][j];}
          if(f1[chan][j1][j]>mf1[chan][2]){mf1[chan][2]=f1[chan][j1][j];}
          if(f1[chan][j1][j]<mf1[chan][3]){mf1[chan][3]=f1[chan][j1][j];}

        }

        for (i = 0; i < 16; i++) {

          //	printf("%lf ",IY[chan][j][i]);

          sg[j][j1] = sg[j][j1] + IY[chan][j][i] * f[chan][j1][i];
          sg1[j][j1] = sg1[j][j1] + IY[chan][j][i] * f1[chan][j1][i];

          sg2[j][j1] = sg2[j][j1] + IY[chan][j][i];



          ssssi = f[chan][j1][i];
          ssssi1 = f1[chan][j1][i];

          gg[j1] = gg[j1] + IY[chan][j][i] * ssssj * ssssi;



          gg1[j1] = gg1[j1] + IY[chan][j][i] * ssssj * ssssi1;
          g1g1[j1] = g1g1[j1] + IY[chan][j][i] * ssssi1 * ssssj1;

          gg2[j1] = gg2[j1] + IY[chan][j][i] * ssssj;
          g1g2[j1] = g1g2[j1] + IY[chan][j][i] * ssssj1;
          g2g2[j1] = g2g2[j1] + IY[chan][j][i];


        }   // for i                                                                                                                                                                          
        //	printf("\n ");

      } //for j            
      //  printf("CalculateFGPar p4\n");
      //    dgg[j1] = gg[j1] * g1g1[j1] - gg1[j1] * gg1[j1];
      dgg1[j1] = gg[j1] * g2g2[j1] - gg2[j1] * gg2[j1];
      dgg2[j1] = gg[j1] * g1g1[j1] * g2g2[j1] - gg1[j1] * gg1[j1] * g2g2[j1] + 2 * gg1[j1] * g1g2[j1] * gg2[j1] - gg2[j1] * gg2[j1] *
          g1g1[j1] - g1g2[j1] * g1g2[j1] * gg[j1];


      for (i = 0; i < 16; i++) {
        if (dgg2[j1] != 0) {

          fg31[chan][j1][i] = ((g1g1[j1] * g2g2[j1] - g1g2[j1] * g1g2[j1]) * sg[i][j1] + (g1g2[j1] * gg2[j1] - gg1[j1] * g2g2[j1]) * sg1[i][j1] +
                        (gg1[j1] * g1g2[j1] - g1g1[j1] * gg2[j1]) * sg2[i][j1]) / dgg2[j1];


          fg32[chan][j1][i] = ((g1g2[j1] * gg2[j1] - gg1[j1] * g2g2[j1]) * sg[i][j1] + (gg[j1] * g2g2[j1] - gg2[j1] * gg2[j1]) * sg1[i][j1] +
                        (gg1[j1] * gg2[j1] - gg[j1] * g1g2[j1]) * sg2[i][j1]) / dgg2[j1];


          fg33[chan][j1][i] = ((gg1[j1] * g1g2[j1] - g1g1[j1] * gg2[j1]) * sg[i][j1] + (gg1[j1] * gg2[j1] - gg[j1] * g1g2[j1]) * sg1[i][j1] +
                        (gg[j1] * g1g1[j1] - gg1[j1] * gg1[j1]) * sg2[i][j1]) / dgg2[j1];

          //printf("%lf \n%lf \n%lf \n", fg31[chan][j1][i], fg32[chan][j1][i], fg33[chan][j1][i]);

          if(i==0){
            if(fg31[chan][j1][i]>mfg31[chan][0]){mfg31[chan][0]=fg31[chan][j1][i];}
            if(fg31[chan][j1][i]<mfg31[chan][1]){mfg31[chan][1]=fg31[chan][j1][i];}
            if(fg32[chan][j1][i]>mfg32[chan][0]){mfg32[chan][0]=fg32[chan][j1][i];}
            if(fg32[chan][j1][i]<mfg32[chan][1]){mfg32[chan][1]=fg32[chan][j1][i];}
            if(fg33[chan][j1][i]>mfg33[chan][0]){mfg33[chan][0]=fg33[chan][j1][i];}
            if(fg33[chan][j1][i]<mfg33[chan][1]){mfg33[chan][1]=fg33[chan][j1][i];}
          }
          else{
            if(fg31[chan][j1][i]>mfg31[chan][2]){mfg31[chan][2]=fg31[chan][j1][i];}
            if(fg31[chan][j1][i]<mfg31[chan][3]){mfg31[chan][3]=fg31[chan][j1][i];}
            if(fg32[chan][j1][i]>mfg32[chan][2]){mfg32[chan][2]=fg32[chan][j1][i];}
            if(fg32[chan][j1][i]<mfg32[chan][3]){mfg32[chan][3]=fg32[chan][j1][i];}
            if(fg33[chan][j1][i]>mfg33[chan][2]){mfg33[chan][2]=fg33[chan][j1][i];}
            if(fg33[chan][j1][i]<mfg33[chan][3]){mfg33[chan][3]=fg33[chan][j1][i];}

          }

        }
    
        /*		uf31[j1][chan]->SetBinContent(i+1,fg31[chan][j1][i]);
        uf32[j1][chan]->SetBinContent(i+1,fg32[chan][j1][i]);
        uf33[j1][chan]->SetBinContent(i+1,fg33[chan][j1][i]);*/
        
      }  // for i 0-16                                                                                                                                                                        



      int jk = 23 + ((48 - j1 ) >> 2);
      //  printf( " !jk=%d j1=%d \n",jk,j1);
      //  printf("CalculateFGPar p5\n");
      if (jk >= 0 && jk < 24 && (48 - j1) % 4 == 0) {
        // printf( "jk=%d \n",jk);
        for(i=0;i<16;i++){

          if (dgg1[j1] != 0) {

            fg41[chan][jk][i] = (g2g2[j1] * sg[i][j1] - gg2[j1] * sg2[i][j1]) / dgg1[j1];
            fg43[chan][jk][i] = (gg[j1] * sg2[i][j1] - gg2[j1] * sg[i][j1]) / dgg1[j1];

            if(i==0){
              if(fg41[chan][jk][i]>mfg41[chan][0]){mfg41[chan][0]=fg41[chan][jk][i];}
              if(fg41[chan][jk][i]<mfg41[chan][1]){mfg41[chan][1]=fg41[chan][jk][i];}
              if(fg43[chan][jk][i]>mfg43[chan][0]){mfg43[chan][0]=fg43[chan][jk][i];}
              if(fg43[chan][jk][i]<mfg43[chan][1]){mfg43[chan][1]=fg43[chan][jk][i];}
            }
            else{
              if(fg41[chan][jk][i]>mfg41[chan][2]){mfg41[chan][2]=fg41[chan][jk][i];}
              if(fg41[chan][jk][i]<mfg41[chan][3]){mfg41[chan][3]=fg41[chan][jk][i];}
              if(fg43[chan][jk][i]>mfg43[chan][2]){mfg43[chan][2]=fg43[chan][jk][i];}
              if(fg43[chan][jk][i]<mfg43[chan][3]){mfg43[chan][3]=fg43[chan][jk][i];}

            }

          }
          //      	uf41[jk][chan]->SetBinContent(i+1,fg41[chan][jk][i]);
          //      	uf43[jk][chan]->SetBinContent(i+1,fg43[chan][jk][i]);
        }  //for i 0-16
  
      } //                                                              

    } // for j1 <endc                                                                                        

  } //for chan

}




void eclFgCal::CalculateBit( ){

  int /*j1, endc, j,*/ i;
  int ChN;
  int i16;
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

  short int bitf;                                                                                            
  short int bitf1;                                                                                           
  short int bitfg31;                                                                                         
  short int bitfg32;                                                                                         
  short int bitfg33;                                                                                         
  short int bitfg41;                                                                                         
  short int bitfg43;

  short int bad[16];

  int badN;

  int cf;
  int cf1;
  int cf31;
  int cf32;
  int cf33;
  int cf41;
  int cf43;
  

  cf = 12;
  cf1 = 12;
  cf31 = 9;
  cf32 = 9;
  cf33 = 13;
  cf41 = 9;
  cf43 = 9;


  //  int ipardsp13;
  int ipardsp14;

  int n16;
  int k;


  for (ChN = 0; ChN < 16; ChN++) {
   // if (!mask_[ChN]) continue;
    n16 = 16;

    //  ipardsp13 = 14 + 14 * 256;
    ipardsp14 = 0 * 256 + 17;

    //  ibb = ipardsp13 / 256;
    //  iaa = ipardsp13 - ibb * 256;
    idd = ipardsp14 / 256;
    //  icc = ipardsp14 - idd * 256;

    //  ia = myPow(2, iaa);
    //  ib = myPow(2, ibb);
    //  cout << "ib set but unused. Fix me. ib = " << ib << endl;

    //  ic = myPow(2, icc);
    //  cout << "ic set but unused. Fix me. ic = " << ic << endl;
    i16 = myPow(2, 15);
    //  ilim = myPow(2, 15);
    //  cout << "ilim set but unused. Fix me. ilim = " << ilim << endl;


    //  for (i = 0; i < 16; i++) {
    for (i = 0; i < 2; i++) {
      
       if (i == 0) {
        idd = n16;
        //      fdd = 16.;
      } else {
        idd = 1;
        //      fdd = 1.;
      }
       //       for (k = 0; k < 192; k++) {
      for (k = 0; k < 2; k++) {

        //      printf("ch=%d %d %lf \n ",ChN,2*i+k,mfg32[ChN][2*i+k]);
        val_f = mf[ChN][2*i+k];
        bf[ChN] = bit_val(val_f, idd, bf[ChN],  1, 1, i16);
        if (bf[ChN] < bitf) {bitf = bf[ChN];}
        val_f = mf1[ChN][2*i+k];
        bf1[ChN] = bit_val(val_f, idd, bf1[ChN],  4, 3, i16);
        if (bf1[ChN] < bitf1) {bitf1 = bf1[ChN];}
        val_f = mfg31[ChN][2*i+k];
        bfg31[ChN] = bit_val(val_f, idd, bfg31[ChN],  1, 1, i16);
        if (bfg31[ChN] < bitfg31) {bitfg31 = bfg31[ChN];}
        val_f = mfg32[ChN][2*i+k];
        bfg32[ChN] = bit_val(val_f, idd, bfg32[ChN],  3, 4, i16);
        if (bfg32[ChN] < bitfg32) {bitfg32 = bfg32[ChN];}
        val_f = mfg33[ChN][2*i+k];
        bfg33[ChN] = bit_val(val_f, idd, bfg33[ChN],  1, 1, i16);
        if (bfg33[ChN] < bitfg33) {bitfg33 = bfg33[ChN];}
        
        val_f = mfg41[ChN][2*i+k];
        bfg41[ChN] = bit_val(val_f, idd, bfg41[ChN],  1, 1, i16);
        if (bfg41[ChN] < bitfg41) {bitfg41 = bfg41[ChN];}
        val_f = mfg43[ChN][2*i+k];
        bfg43[ChN] = bit_val(val_f, idd, bfg43[ChN],  1, 1, i16);
        if (bfg43[ChN] < bitfg43) {bitfg43 = bfg43[ChN];}
        



        //      if(ChN==15) cout << " cr="<<cr<< " ChN="<< ChN <<"Pec k = " << k<< " i=" <<i << " f=" << f[ChN][k][i]<< " f1=" << f1[ChN][k][i] << endl;
        //////////////////

        /****
        val_f = f[ChN][k][i];
        bf[ChN] = bit_val(val_f, idd, bf[ChN],  1, 1, i16);
        if (bf[ChN] < bitf) {bitf = bf[ChN];}
        val_f = f1[ChN][k][i];

        bf1[ChN] = bit_val(val_f, idd, bf1[ChN],  4, 3, i16);
        if (bf1[ChN] < bitf1) {bitf1 = bf1[ChN];}
        val_f = fg31[ChN][k][i];
        bfg31[ChN] = bit_val(val_f, idd, bfg31[ChN],  1, 1, i16);
        if (bfg31[ChN] < bitfg31) {bitfg31 = bfg31[ChN];}
        val_f = fg32[ChN][k][i];
        bfg32[ChN] = bit_val(val_f, idd, bfg32[ChN],  3, 4, i16);
        if (bfg32[ChN] < bitfg32) {bitfg32 = bfg32[ChN];}
        val_f = fg33[ChN][k][i];
        bfg33[ChN] = bit_val(val_f, idd, bfg33[ChN],  1, 1, i16);
        if (bfg33[ChN] < bitfg33) {bitfg33 = bfg33[ChN];}
        if (k < 24) {
          val_f = fg41[ChN][k][i];
          bfg41[ChN] = bit_val(val_f, idd, bfg41[ChN],  1, 1, i16);
          if (bfg41[ChN] < bitfg41) {bitfg41 = bfg41[ChN];}
          val_f = fg43[ChN][k][i];
          bfg43[ChN] = bit_val(val_f, idd, bfg43[ChN],  1, 1, i16);
          if (bfg43[ChN] < bitfg43) {bitfg43 = bfg43[ChN];}
        }



        if (fabs(f[ChN][k][i] / fdd)    > xf)   {xf    = fabs(f[ChN][k][i] / fdd);}
        if (fabs(f1[ChN][k][i] / fdd)   > xf1)   {xf1   = fabs(f1[ChN][k][i] / fdd);}
        if (fabs(fg31[ChN][k][i] / fdd) > xfg31) {xfg31 = fabs(fg31[ChN][k][i] / fdd);}
        if (fabs(fg32[ChN][k][i] / fdd) > xfg32) {xfg32 = fabs(fg32[ChN][k][i] / fdd);}
        if (fabs(fg33[ChN][k][i] / fdd) > xfg33) {xfg33 = fabs(fg33[ChN][k][i] / fdd);}
        if (k < 24) {
          if (fabs(fg41[ChN][k][i] / fdd) > xfg41) {xfg41 = fabs(fg41[ChN][k][i] / fdd);}
          if (fabs(fg43[ChN][k][i] / fdd) > xfg43) {xfg43 = fabs(fg43[ChN][k][i] / fdd);}
        }

        if (fabs(f[ChN][k][i] / fdd)    > xf)   {xf    = fabs(f[ChN][k][i] / fdd);}
        if (fabs(f1[ChN][k][i] / fdd)   > xf1)   {xf1   = fabs(f1[ChN][k][i] / fdd);}
        if (fabs(fg31[ChN][k][i] / fdd) > xfg31) {xfg31 = fabs(fg31[ChN][k][i] / fdd);}
        if (fabs(fg32[ChN][k][i] / fdd) > xfg32) {xfg32 = fabs(fg32[ChN][k][i] / fdd);}
        if (fabs(fg33[ChN][k][i] / fdd) > xfg33) {xfg33 = fabs(fg33[ChN][k][i] / fdd);}
        if (k < 24) {
          if (fabs(fg41[ChN][k][i] / fdd) > xfg41) {xfg41 = fabs(fg41[ChN][k][i] / fdd);}
          if (fabs(fg43[ChN][k][i] / fdd) > xfg43) {xfg43 = fabs(fg43[ChN][k][i] / fdd);}
        }

        ////////////////////////////////
        */
      } // for k                                                                                                                                                                                                                         



    }  // for i                                                                                                                                                                                
     /*                                          
    bf[ChN]=14;
    bf1[ChN]=14;
    bfg31[ChN]=15;
    bfg32[ChN]=14;
    bfg33[ChN]=17;
     */

     //     printf("aft %d %d %d %d %d \n",ChN,bf[ChN],bfg31[ChN],bfg32[ChN],bfg33[ChN]);
    if (bf[ChN] < cf || bf1[ChN] < cf1 || bfg31[ChN] < cf31 || bfg32[ChN] < cf32 || bfg33[ChN] < cf33 || bfg41[ChN] < cf41
      || bfg43[ChN] < cf43) {
      badN = badN + 1;

      if (bf[ChN] < cf) {bad[ChN] = bad[ChN] + 1;}
      if (bf1[ChN] < cf1) {bad[ChN] = bad[ChN] + 2;}
      if (bfg31[ChN] < cf31) {bad[ChN] = bad[ChN] + 4;}
      if (bfg32[ChN] < cf32) {bad[ChN] = bad[ChN] + 8;}
      if (bfg33[ChN] < cf33) {bad[ChN] = bad[ChN] + 16;}
      if (bfg41[ChN] < cf41) {bad[ChN] = bad[ChN] + 32;}
      if (bfg43[ChN] < cf43) {bad[ChN] = bad[ChN] + 64;}

    }

    if(bf1[ChN]<bf[ChN]){bf[ChN]=bf1[ChN];}

    if (bf[ChN] < cf) {printf("bf out of range Ch=%d bf=%d \n",ChN,bf[ChN]);}
    if (bf1[ChN] < cf1) {printf("bf1 out of range Ch=%d bf1=%d \n",ChN,bf1[ChN]);}
    if (bfg31[ChN] < cf31) {printf("bfg31 out of range Ch=%d bf=%d \n",ChN,bfg31[ChN]);}
    if (bfg32[ChN] < cf32) {printf("bfg32 out of range Ch=%d bf=%d \n",ChN,bfg32[ChN]);}
    if (bfg33[ChN] < cf33) {printf("bfg33 out of range Ch=%d bf=%d \n",ChN,bfg33[ChN]);}
    if (bfg41[ChN] < cf41) {printf("bfg41 out of range Ch=%d bf=%d \n",ChN,bfg41[ChN]);}
    if (bfg43[ChN] < cf43) {printf("bfg43 out of range Ch=%d bf=%d \n",ChN,bfg43[ChN]);}


    if (bfg31[ChN] > bfg41[ChN]) {bfg31[ChN]=bfg41[ChN];}
    if (bfg33[ChN] > bfg43[ChN]) {bfg33[ChN]=bfg43[ChN];}


    if (bf[ChN] < cf) {bf[ChN]=cf;}
   
    if (bfg31[ChN] < cf31) {bfg31[ChN]=cf31;}
    if (bfg32[ChN] < cf32) {bfg32[ChN]=cf32;}
    if (bfg33[ChN] < cf33) {bfg33[ChN]=cf33;}


    /*
    if (bf[ChN] < 14) {printf("warning bf  Ch=%d bf=%d \n",ChN,bf[ChN]);}
    if (bf1[ChN] < 14) {printf("warning bf1 Ch=%d bf1=%d \n",ChN,bf1[ChN]);}
    if (bfg31[ChN] < 15) {printf("warning  bfg31  Ch=%d bf=%d \n",ChN,bfg31[ChN]);}
    if (bfg32[ChN] < 14) {printf("warning  bfg32 Ch=%d bf=%d \n",ChN,bfg32[ChN]);}
    if (bfg33[ChN] < 17) {printf("warning  bfg33 Ch=%d bf=%d \n",ChN,bfg33[ChN]);}
    */

    //     printf("in %d %d %d %d %d \n",ChN,bf[ChN],bfg31[ChN],bfg32[ChN],bfg33[ChN]);

    /*
      Kf->SetBinContent(ChN+1,bf[ChN]);
      Ka->SetBinContent(ChN+1,bfg31[ChN]);
      Kb->SetBinContent(ChN+1,bfg32[ChN]);
      Kc->SetBinContent(ChN+1,bfg33[ChN]);
    */
  }//ChN


  short int  bfc;
  short int  bf31;
  short int  bf32;
  short int  bf33;
  short int  bf41;
  short int  bf43;

  bfc = 20; 
  bf31 = 20; 
  bf32 = 20; 
  bf33 = 20; 
  bf41 = 20; 
  bf43 = 20; 

  bfc = 15; 
  bf31 = 16; 
  bf32 = 16; 
  bf33 = 19; 
  bf41 = 16; 
  bf43 = 19; 

  for (int chn = 0; chn < 16; chn++){

    if(bf[chn] < bfc){bfc=bf[chn];}
    if(bfg31[chn] < bf31){bf31=bfg31[chn];}
    if(bfg32[chn] < bf32){bf32=bfg32[chn];}
    if(bfg33[chn] < bf33){bf33=bfg33[chn];}
    if(bfg41[chn] < bf41){bf41=bfg41[chn];}
    if(bfg43[chn] < bf43){bf43=bfg43[chn];}
  }

  //   printf("  %d %d  %d  %d  %d  \n",bfc,bf31,bf32,bf33);

  for (int chn = 0; chn < 16; chn++){
    //            printf(" do %d %d  %d  %d  %d  \n",chn,bf[chn],bfg31[chn],bfg32[chn],bfg33[chn]);

    bf[chn] = bfc;
    bfg31[chn] = bf31;
    bfg32[chn] = bf32;
    bfg33[chn] = bf33;
    bfg41[chn] = bf41;
    bfg43[chn] = bf43;

   // printf("chn=%d, mf=%lf, mf1=%lf, mfg31=%lf, mfg32=%lf, mfg33=%lf, mfg41=%lf, mfg43=%lf,\n",chn, mf[chn][0], mf1[chn][0], mfg31[chn][0], mfg32[chn][0], mfg33[chn][0], mfg41[chn][0], mfg43[chn][0]);

    /*Kf->SetBinContent(chn+1,bf[chn]);
  Ka->SetBinContent(chn+1,bfg31[chn]);
  Kb->SetBinContent(chn+1,bfg32[chn]);
  Kc->SetBinContent(chn+1,bfg33[chn]);*/
  //        printf(" posle %d %d  %d  %d  %d  \n",chn,bf[chn],bfg31[chn],bfg32[chn],bfg33[chn]);

  }

  //printf("bfc=%d, bf31=%d, bf32=%d, bf33=%d\n", bfc, bf31, bf32, bf33);



}

void eclFgCal::CalculateFginInt( ){

  //  int ChN;
  //  int k;
  //  int i;

  int iff;
  int ia;
  int ib;
  int ic;

  int imax1;

  imax1=myPow(2,15);
  double imax;
  imax=(double)imax1;


  //  for (ChN = 9; ChN < 10; ChN++) {

  for (int ChN = 0; ChN < 16; ChN++) {

    //if (!mask_[ChN]) continue;
    //    printf("calc %d %d %d %d  \n",ChN,bf[ChN],bfg31[ChN],bfg32[ChN],bfg33[ChN]);
   
    iff = 1 << bf[ChN];
    ia = 1 << bfg31[ChN];
    ib = 1 << bfg32[ChN];
    ic = 1 << bfg33[ChN];
    
    
    
    
    
    //   printf("%d %d %d %d %lf  \n",iff,ia,ib,ic,imax);
    
    //   printf("%d %d %d %d \n",bf[ChN],bfg31[ChN],bfg32[ChN],bfg33[ChN]);
    for (int k = 0; k < 192; k++) {
      

      const double isd = 3. / 4., sd = 1 / isd ; // conversion factor (???) 
      
      for (int i = 0; i < 16; i++) {
	double w = i ? 1.0 : 1. / 16.;
	
	
	
	m_f[ChN][k][i] = lrint(f[ChN][k][i]  * iff * w);
	m_f1[ChN][k][i] = lrint(f1[ChN][k][i] * iff * w * sd);
	
	m_fg31[ChN][k][i] = lrint(fg31[ChN][k][i] * ia * w);
	m_fg32[ChN][k][i] = lrint(fg32[ChN][k][i] * ib * w * isd);
	m_fg33[ChN][k][i] = lrint(fg33[ChN][k][i] * ic * w);
	
	/*
	  m_fg31[ChN][k][i] = 0;
	  m_fg32[ChN][k][i] = 0;
	  m_fg33[ChN][k][i] = 0;
	*/
	
	if(f[ChN][k][i]  * iff * w>imax || f[ChN][k][i]  * iff * w<-imax ){
	  printf("alarm ChN=%d k=%d i=%d m_f=%lf imax=%lf \n",ChN,k,i,f[ChN][k][i]  * iff * w,imax);}        
	if(f1[ChN][k][i] * iff * w * sd>imax || f1[ChN][k][i] * iff * w * sd<-imax ){
	  printf("alarm ChN=%d k=%d i=%d m_f1=%lf imax=%lf \n",ChN,k,i,f1[ChN][k][i] * iff * w * sd,imax);}        
	if(fg31[ChN][k][i] * ia * w>imax || fg31[ChN][k][i] * ia * w<-imax ){
	  printf("alarm ChN=%d k=%d i=%d m_fg31=%lf imax=%lf \n",ChN,k,i,fg31[ChN][k][i] * ia * w,imax);}        
	if(fg32[ChN][k][i] * ib * w * isd>imax || fg32[ChN][k][i] * ib * w * isd<-imax ){
	  printf("alarm ChN=%d k=%d i=%d m_fg32=%lf imax=%lf \n",ChN,k,i,fg32[ChN][k][i] * ib * w * isd,imax);}        
	if(fg33[ChN][k][i] * ic * w>imax || fg33[ChN][k][i] * ic * w<-imax ){
	  printf("alarm ChN=%d k=%d i=%d m_fg33=%lf imax=%lf \n",ChN,k,i,fg33[ChN][k][i] * ic * w,imax);}        

   //  printf("!!! ch=%d i=%d k=%d f=%lf f1=%lf fg31=%lf fg32=%lf fg33=%lf m_f=%d m_f1=%d fm_g31=%d m_fg32=%d m_fg33=%d  ia=%d ib=%d ic=%d w=%lf \n",ChN,i,k,f[ChN][k][i],f1[ChN][k][i],fg31[ChN][k][i],fg32[ChN][k][i],fg33[ChN][k][i],m_f[ChN][k][i],m_f1[ChN][k][i],m_fg31[ChN][k][i],m_fg32[ChN][k][i],m_fg33[ChN][k][i],ia,ib,ic,w);
	  
	/*
	  
	  if(k<24){      printf("!!! ch=%d i=%d k=%d f=%lf f1=%lf fg31=%lf fg32=%lf fg33=%lf fg41=%lf fg43=%lf ia=%d ib=%d ic=%d w=%lf \n",ChN,i,k,f[ChN][k][i],f1[ChN][k][i],fg31[ChN][k][i],fg32[ChN][k][i],fg33[ChN][k][i],fg41[ChN][k][i],fg43[ChN][k][i],ia,ib,ic,w);
	  printf("m ch=%d i=%d k=%d f=%d f1=%d  fg31=%d fg32=%d fg33=%d fg41=%d fg43=%d  ia=%d ib=%d ic=%d w=%lf \n",ChN,i,k,m_f[ChN][k][i],m_f1[ChN][k][i],m_fg31[ChN][k][i],m_fg32[ChN][k][i],m_fg33[ChN][k][i],m_fg41[ChN][k][i],m_fg43[ChN][k][i],ia,ib,ic,w);
	  }
	  else{    
	  printf("!!! ch=%d i=%d k=%d f=%lf f1=%lf fg31=%lf fg32=%lf fg33=%lf  ia=%d ib=%d ic=%d w=%lf \n",ChN,i,k,f[ChN][k][i],f1[ChN][k][i],fg31[ChN][k][i],fg32[ChN][k][i],fg33[ChN][k][i],ia,ib,ic,w);
	  printf("m ch=%d i=%d k=%d f=%d f1=%d  fg31=%d fg32=%d fg33=%d  ia=%d ib=%d ic=%d w=%lf \n",ChN,i,k,m_f[ChN][k][i],m_f1[ChN][k][i],m_fg31[ChN][k][i],m_fg32[ChN][k][i],m_fg33[ChN][k][i],ia,ib,ic,w);
	  
	  }
	*/
	
	int jk = 23 + ((48 - k ) >> 2);
	if (jk >= 0 && jk < 24 && (48 - k) % 4 == 0) {
	   //long long int sum =0;
      //double dsum =0;
	  for (int i1 = 0; i1 < 16; i1++) {
	    double w1 = i1 ? 1.0 : 1. / 16.;
	    
	    m_fg41[ChN][jk][i1] = lrint(fg41[ChN][jk][i1] * ia * w1);
	    m_fg43[ChN][jk][i1] = lrint(fg43[ChN][jk][i1] * ic * w1);
	    
     // sum+=	m_fg41[ChN][jk][i1]/w1*m_f[ChN][jk][i1]/myPow(2,abs(bf[ChN]-bfg31[ChN]));
     // dsum+=	fg41[ChN][jk][i1]*f[ChN][jk][i1];
	    /*
	      m_fg41[ChN][jk][i] = 0;
     m_fg43[ChN][jk][i] = 0;
	    */
	    
	    if(fg41[ChN][jk][i1] * ia * w>imax || fg41[ChN][jk][i1] * ia * w1<-imax ){
	      printf("alarm ChN=%d k=%d i=%d m_fg41=%lf imax=%lf \n",ChN,k,i1,fg41[ChN][jk][i1] * ia * w1,imax);}        
	    if(fg43[ChN][jk][i1] * ic * w1>imax || fg43[ChN][jk][i1] * ic * w1<-imax ){
	      printf("alarm ChN=%d k=%d i=%d m_fg43=%lf imax=%lf \n",ChN,k,i1,fg43[ChN][jk][i1] * ic * w1,imax);}        
	    
	    
	    //  printf("!!! ch=%d i=%d k=%d fg41=%lf fg43=%lf m_fg41=%d m_fg43=%d  ia=%d ib=%d ic=%d w=%lf \n",ChN,i,k,fg41[ChN][k][i],fg43[ChN][k][i],m_fg41[ChN][k][i],m_fg43[ChN][k][i],ia,ib,ic,w);
	  }

  //sum=sum/ia;
  //double ddsum=sum;
  //ddsum=ddsum/ia;
  //printf("%lf %lf %d\n", ddsum, dsum, ia);
	}
	
	
      }
      
      
    }
    
    /*
    for (int k = 0; k < 192; k++) {
      
      for (int i = 0; i < 16; i++) {
	
	if(k<24){     
	  // printf("!!! ch=%d i=%d k=%d f=%18.15e f1=%18.15e fg31=%lf fg32=%lf fg33=%lf fg41=%lf fg43=%lf ia=%d ib=%d ic=%d  \n",ChN,i,k,f[ChN][k][i],f1[ChN][k][i],fg31[ChN][k][i],fg32[ChN][k][i],fg33[ChN][k][i],fg41[ChN][k][i],fg43[ChN][k][i],ia,ib,ic);
	  //      printf("m ch=%d i=%d k=%d f=%d f1=%d  fg31=%d fg32=%d fg33=%d fg41=%d fg43=%d  ia=%d ib=%d ic=%d  \n",ChN,i,k,m_f[ChN][k][i],m_f1[ChN][k][i],m_fg31[ChN][k][i],m_fg32[ChN][k][i],m_fg33[ChN][k][i],m_fg41[ChN][k][i],m_fg43[ChN][k][i],ia,ib,ic);
	}
	else{    
	  //   printf("!!! ch=%d i=%d k=%15.12d f=%15.12e f1=%e fg31=%lf fg32=%lf fg33=%lf  ia=%d ib=%d ic=%d  \n",ChN,i,k,f[ChN][k][i],f1[ChN][k][i],fg31[ChN][k][i],fg32[ChN][k][i],fg33[ChN][k][i],ia,ib,ic);
	  //      printf("m ch=%d i=%d k=%d f=%d f1=%d  fg31=%d fg32=%d fg33=%d  ia=%d ib=%d ic=%d  \n",ChN,i,k,m_f[ChN][k][i],m_f1[ChN][k][i],m_fg31[ChN][k][i],m_fg32[ChN][k][i],m_fg33[ChN][k][i],ia,ib,ic);

	}
	
      }
      
      
    }*/
  //printf("chn=%d, bf=%d, bf1=%d, bfg31=%d, bfg32=%d, bfg33=%d, bfg41=%d, bfg43=%d \n",ChN, bf[ChN], bf1[ChN],bfg31[ChN], bfg32[ChN], bfg33[ChN], bfg41[ChN], bfg43[ChN]);

  //printf("chn=%d, mf=%lf, mf1=%lf, mfg31=%lf, mfg32=%lf, mfg33=%lf, mfg41=%lf, mfg43=%lf \n",ChN, mf[ChN][0], mf1[ChN][0], mfg31[ChN][0], mfg32[ChN][0], mfg33[ChN][0], mfg41[ChN][0], mfg43[ChN][0]);

    
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

  short int  bfc;
  short int  bf31;
  short int  bf32;
  short int  bf33;
  short int  bf41;
  short int  bf43;

  bfc = 20; 
  bf31 = 20; 
  bf32 = 20; 
  bf33 = 20; 
  bf41 = 20; 
  bf43 = 20; 


  // fill thresholds                                                                                                                                                                    
  for (int chn = 0; chn < 16; chn++){
    id[128 + chn] = 128 + 40;  // A_low
    id[192 + chn] = 128 - 30;  // A_skip
    id[64  + chn] =  -100;  // A_hard

    if(bf[chn] < bfc){bfc=bf[chn];}
    if(bfg31[chn] < bf31){bf31=bfg31[chn];}
    if(bfg32[chn] < bf32){bf32=bfg32[chn];}
    if(bfg33[chn] < bf33){bf33=bfg33[chn];}
    if(bfg41[chn] < bf41){bf41=bfg41[chn];}
    if(bfg43[chn] < bf43){bf43=bfg43[chn];}


    /*    
  for (int i = 0; i < 16; i++){
  for (int k = 0; k < 192; k++){

    if(k>23){    printf("ChN=%d k=%d i=%d \t f=%d f1=%d f31=%d f32=%d f33=%d \n",chn,k,i,m_f[chn][k][i],m_f1[chn][k][i],m_fg31[chn][k][i],m_fg32[chn][k][i],m_fg33[chn][k][i]);
  }
    else{      printf("ChN=%d k=%d i=%d \t f=%d f1=%d f31=%d f32=%d f33=%d f41=%d f43=%d \t \n",chn,k,i,m_f[chn][k][i],m_f1[chn][k][i],m_fg31[chn][k][i],m_fg32[chn][k][i],m_fg33[chn][k][i],m_fg41[chn][k][i],m_fg43[chn][k][i]);}

      }
      }
    */
    
  }

    
  // fill some parameters                                               
                                                                        
  /*     
  id[26] = bf31;                                   
  id[27] =  bf32;                                   
  id[28] =  bf33;                                   
  id[29] =  16;                                  
  id[32] =  14;                                  
  id[33] =  10;                                   
  */

  //  printf(" A  %d   %d  %d  %d  \n",bfc,bf31,bf32,bf33);

     
  id[13] = bf31+256*bf32;  // ka kb                                   
  id[14] =  bf33+0*256;  // kc  y0                                   
  //  id[15] =  chi+thres;                                  

  id[15] =  30000;                                  
  //  id[16] = 14+k2_chi*256;   //k1_chi k2_chi
  id[16] = bfc+10*256;   //k1_chi k2_chi
  id[8]=257;
  id[8]=513;  // new version
  id[9]=-1;

  for (int i = 17; i < 64; i++){
    id[i]=0;
  }

  id[17] = 0;

  /*


  k_a  = (int) * ((unsigned char*)id + 26);
  k_b  = (int) * ((unsigned char*)id + 27);
  k_c  = (int) * ((unsigned char*)id + 28);
  k_16 = (int) * ((unsigned char*)id + 29);

  k1_chi = (int) * ((unsigned char*)id + 32); // in sim 32->24 why??                                                                                                                    
  k2_chi = (int) * ((unsigned char*)id + 33); // in sim 33->25 why??                                                                                                                    
  chi_thres = (int) * (id + 15);
  kz_s =  (int)*((unsigned char*)id+34);
  */

  for (int i = 0; i < 256; i++){
    //printf("%d ",id[i]);
    //    if(i%16==0){printf(" Q %d\n",i/16);}

  }

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



  for(int i=0;i<16;i++)
    {

      printf("eclFgCal::WriteCalib: i=%d\n",i);
      nsiz1=384;
      h = fwrite(m_fg41[i][0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
	std::cout << "Error writing fg41 data. Read block size = " << h << '\n';
      exit(0);

      }

      nsiz1=3072;
      h = fwrite(m_fg31[i][0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
	      std::cout << "Error writing fg31 data. Read block size = " << h << '\n';
         exit(0);
      }


      h = fwrite(m_fg32[i][0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
      	std::cout << "Error writing fg32 data. Read block size = " << h << '\n';
        exit(0);
      }

      h = fwrite(m_fg33[i][0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
	      std::cout << "Error writing fg33 data. Read block size = " << h << '\n';
        exit(0);
      }

      nsiz1=384;
      h = fwrite(m_fg43[i][0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
	      std::cout << "Error writing fg43 data. Read block size = " << h << '\n';
        exit(0);
      }


      nsiz1=3072;
      h = fwrite(m_f[i][0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
        std::cout << "Error writing f data. Read block size = " << h << '\n';
        exit(0);
      }

      h = fwrite(m_f1[i][0], nsiz, nsiz1, fr);
      if (h != nsiz1) {
	      std::cout << "Error writing f1 data. Read block size = " << h << '\n';
        exit(0);
      }

      printf("eclFgCal::WriteCalib1: i=%d\n",i);
    }



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
  int u;

  fr= fopen(fname, "r");

  if(fr == NULL){
    std::cout << " can't open  file "
      << fname << "\n";
    exit(0);
  }

  fscanf (fr, "%lf   ", &a0);

  while(!feof(fr)){
    fscanf (fr, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf    ", &u,&a0,&a1,&a2,&a3,&a4,&a5,&a6,&a7,&a8,&a9);
    if(u>0){
      par[u-1][0]=a0+0.25;
      //      par[u-1][0]=a0+0.25-3./192.;
      //      par[u-1][0]=a0;
      par[u-1][1]=a1;
      par[u-1][2]=a2;
      par[u-1][3]=a3;
      par[u-1][4]=a4;
      par[u-1][5]=a5;
      par[u-1][6]=a6;
      par[u-1][7]=a7;
      par[u-1][8]=a8;
      par[u-1][9]=a9;
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
