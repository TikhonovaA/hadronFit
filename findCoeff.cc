#include <TF1.h>
#include <TH1F.h>
#include <TFile.h>
//#include<unuran.h>
#include<iostream> 
#include<iomanip> 
#include<fstream> 
#include<stdlib.h> 
#include<stdio.h> 
#include<math.h> 
#include<fcntl.h> 
#include<string.h> 
#include <cstdlib>
#include <time.h>
#include <sys/timeb.h>
#include "TTree.h"
#include "TRandom.h"
#include <TSystem.h>
#include "TMinuit.h"
#include <Eigen/Dense>
#include "funcCalc.h"

int ttrig;


void minimize(short int *id, int *pf, int *pf1, int *hf, int *hf1, int *fg41, int *fg43, int *fg45, int *fg51, int *fg52, int *fg53, int *fg54, int *fg55,
        int *ttrig, int *ch, double a/*, double b*/, Double_t *lfunc)
{

	int A0  = (int)*(id+128+*ch)-128;
	int Askip  = (int)*(id+192+*ch)-128;

	int Ahard  = (int)*(id+64+*ch);
	int k_a = (int)*((unsigned char*)id+26);
	int k_b = (int)*((unsigned char*)id+27);
	int k_c = (int)*((unsigned char*)id+28);
	int k_d = (int)*((unsigned char*)id+29);
	int k_e = (int)*((unsigned char*)id+30);
	int k_16 = (int)*((unsigned char*)id+31);

	int k1_chi =  (int)*((unsigned char*)id+34);
	int k2_chi =  (int)*((unsigned char*)id+35);
	int chi_thres = (int)*(id+16);
	long long int s1,s2,s3,s4,B1, B2;

	int it, it0;
	double itd, norm; 
	int it_h, it_l; 
	int low_ampl,i,T,iter;   

	double chi2 = 0;
    double hlvl = 0.0; 


	it0 = 48 + ((143-*ttrig)*2/3);

	//printf(" it0 = %d\n", it0);

	//limits
	it_h=191;
	it_l=0;

	if (it0 < it_l)it0=it_l; 
	if (it0 > it_h)it0=it_h; 

    for(int d = -48; d <= 48; d++){
        long long int v51,v52,v53,v54,  h51,h52,h53,h54;
        double w1,w2,w3,w4;
        v51=v52=v53=v54=h51=h52=h53=h54=0;

        it = it0+d;
        
        for(int j=0; j<16; j++){
            s1=(*(fg51+it*16+j));
		    s2=(*(fg52+it*16+j));
		    s3=(*(fg53+it*16+j));
		    s4=(*(fg54+it*16+j));

            B1 = (*(pf+it0*16+j));
            B2 = s1*B1;
            v51 += B2;

            B2 = s2*B1;
            v52 += B2;

            B2 = s3*B1;
            v53 += B2;

            B2 = s4*B1;
            v54 += B2;

            B1 = (*(hf+it0*16+j));
            B2 = s1*B1;
            h51 += B2;

            B2 = s2*B1;
            h52 += B2;

            B2 = s3*B1;
            h53 += B2;

            B2 = s4*B1;
            h54 += B2;
        }

        w1 = ( (1-hlvl) * (double)v51 + hlvl* (double)h51 ) / pow(2,30+k_a);
        w2 = ( (1-hlvl) * (double)v52 + hlvl* (double)h52 ) / pow(2,30+k_b);
        w3 = ( (1-hlvl) * (double)v53 + hlvl* (double)h53 ) / pow(2,30+k_c);
        w4 = ( (1-hlvl) * (double)v54 + hlvl* (double)h54 ) / pow(2,30+k_d);

        chi2 += ( (double)(it - it0)/256 - (w2 + a*w4)/(w1 + w3) ) * ( (double)(it - it0)/256 - (w2 + a*w4)/(w1 + w3) );
    }

    *lfunc = chi2;

   // printf("chi2=%lf\n", chi2);
	
}


void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

	double paq[16][12];
	char in[256];
	strcpy(in, std::getenv("HOME"));
	strcat(in, "/fit/panr_11.dat");

	//fillparam(in, paq);
 
	int chan = 15;

	TRandom *random = new TRandom();
	random->SetSeed((ULong_t)time(NULL));


    //read dsp.dat for one channel
    FILE *fr;
    int h,nsiz,nsiz1;
    char fname[256];
    strcpy(fname, std::getenv("HOME"));
    strcat(fname, "/fit/DSP_exp10/dsp.dat");

    short int id[256];
    int cpf[192][16];
    int cpf1[192][16];
    int chf[192][16];
    int chf1[192][16];
    int cfg51[192][16];
    int cfg52[192][16];
    int cfg53[192][16];
    int cfg54[192][16];
    int cfg55[192][16];
    int cfg41[24][16];
    int cfg43[24][16];
    int cfg45[24][16];

    for (int i = 0; i < 256; i++){
        id[i]=0;
    }

    for ( int u =0 ; u < 16; u ++){
    for ( int k =0 ; k < 192; k ++){
        if(k < 24){
        cfg41[k][u] = 0;
        cfg43[k][u] = 0;
        cfg45[k][u] = 0;
        }
    
        cfg51[k][u] = 0;
        cfg52[k][u] = 0;
        cfg53[k][u] = 0;
        cfg54[k][u] = 0;
        cfg55[k][u] = 0;
        cpf[k][u]    = 0;
        cpf1[k][u]   = 0;
        chf[k][u]    = 0;
        chf1[k][u]   = 0;

    }
    }

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
        std::cout << "Error reading id data. Read block size = " << h << '\n';
        exit(0);
    }


    nsiz=4;

    for(int i=0;i<=chan;i++){

        nsiz1=384;
        h = fread(cfg41[0], nsiz, nsiz1, fr);
        if (h != nsiz1) {
            std::cout << "Error reading fg41 data. Read block size = " << h << '\n';
            exit(0);
        }

        h = fread(cfg43[0], nsiz, nsiz1, fr);
        if (h != nsiz1) {
            std::cout << "Error reading fg43 data. Read block size = " << h << '\n';
            exit(0);
        }

        nsiz1=3072;
        h = fread(cfg51[0], nsiz, nsiz1, fr);
        if (h != nsiz1) {
            std::cout << "Error reading fg51 data. Read block size = " << h << '\n';
            exit(0);
        }

        h = fread(cfg52[0], nsiz, nsiz1, fr);
        if (h != nsiz1) {
            std::cout << "Error reading fg52 data. Read block size = " << h << '\n';
            exit(0);
        }

        h = fread(cfg53[0], nsiz, nsiz1, fr);
        if (h != nsiz1) {
            std::cout << "Error reading fg53 data. Read block size = " << h << '\n';
            exit(0);
        }

        h = fread(cfg54[0], nsiz, nsiz1, fr);
        if (h != nsiz1) {
            std::cout << "Error reading fg54 data. Read block size = " << h << '\n';
            exit(0);
        }

        h = fread(cfg55[0], nsiz, nsiz1, fr);
        if (h != nsiz1) {
            std::cout << "Error reading fg55 data. Read block size = " << h << '\n';
            exit(0);
        }

        nsiz1=384;
        h = fread(cfg45[0], nsiz, nsiz1, fr);
        if (h != nsiz1) {
            std::cout << "Error reading fg45 data. Read block size = " << h << '\n';
            exit(0);
        }

        nsiz1=3072;
        h = fread(cpf[0], nsiz, nsiz1, fr);
        if (h != nsiz1) {
            std::cout << "Error reading pf data. Read block size = " << h << '\n';
            exit(0);
        }

        h = fread(cpf1[0], nsiz, nsiz1, fr);
        if (h != nsiz1) {
            std::cout << "Error reading pf1 data. Read block size = " << h << '\n';
            exit(0);
        }

        h = fread(chf[0], nsiz, nsiz1, fr);
        if (h != nsiz1) {
            std::cout << "Error reading hf data. Read block size = " << h << '\n';
            exit(0);
        }

        h = fread(chf1[0], nsiz, nsiz1, fr);
        if (h != nsiz1) {
            std::cout << "Error reading hf1 data. Read block size = " << h << '\n';
            exit(0);
        }

    }//i

    fclose(fr);

    Double_t lfunc;

    minimize(id, cpf[0], cpf1[0], chf[0], chf1[0], cfg41[0], cfg43[0], cfg45[0], cfg51[0], cfg52[0], cfg53[0], cfg54[0], cfg55[0], &ttrig, &chan, par[0], &lfunc);
    
    f = lfunc;
}


int findCoeff(){

    gSystem->Load("libMinuit");
    const int npar = 1;

    double results[96][3];
    for(int t=48; t<144; t++){
        double ts = (3./2. * (48-(double)t) - 75) / 288.;
        printf("ts = %lf\n", ts);
        ttrig = 218+288*ts;  //129+288*ts; //trigger time

        TMinuit *gMinuit = new TMinuit(npar);  
        gMinuit->SetFCN(fcn);

        Double_t arglist[10];
        Int_t ierflg = 0;

        arglist[0] = 0.5;
        gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

        // Set starting values and step sizes for parameters		
        const Double_t vstart[npar] = {1};
        const Double_t step[npar] =   {0.001};

        gMinuit->mnparm(0, "a", vstart[0], step[0], 0,0,ierflg);
       // gMinuit->mnparm(1, "b", vstart[1], step[1], 0,0,ierflg);
       // gMinuit->mnparm(2, "c", vstart[2], step[2], 0,0,ierflg);

        arglist[0] = 500;
        arglist[1] = 1.;
        gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

        Double_t pval[npar], err[npar];
        for(int i=0; i<npar; i++ ) {
            gMinuit->GetParameter(i, pval[i], err[i]);
        }

        results[t-48][0] = t;
        results[t-48][1] = pval[0];
        results[t-48][2] = err[0];

       // results[t-48][3] = pval[1];
       // results[t-48][4] = err[1];

       // results[t-48][5] = pval[2];
       // results[t-48][6] = err[2];

    }


    for(int i =0; i<96; i++){
        printf("%lf, %lf +- %lf\n", results[i][0], results[i][1], results[i][2]);
      
    }

    return 0;

}