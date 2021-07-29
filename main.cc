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
#include <Eigen/Dense>
#include "funcCalc.h"

/*

double sv123(const double &t, const double &t01, const double &tb1, const double &t02, 
		const double &tb2, const double &td1, const double &ts1){
	//      implicit none
  double sv123 = 0.;

	double dt01,dtb1,dt02,dtb2,dks0,dks1,dksm;
	double dw0,dw1,dwp,dwm,das1,dac1,das0,dac0,dzna,dksm2,ds,dd;
	double dcs0,dsn0,dzn0,td,ts,dr;
	double  dcs0s,dsn0s,dcs0d,dsn0d,dcs1s,dsn1s,dcs1d,dsn1d;

	double rv=0.0;
		if(t < 0) return 0;

	dr=(ts1-td1)/td1;
	//	if(fabs(dr) > 1.0e-5){
    if (fabs(dr) >= 0.0000001) {

		td=td1;
		ts=ts1;
	} else {
		td=td1;
		if(ts1>td1)
			ts=td1*1.00001;
		else
			ts=td1*0.99999;

	}

	dr=((t01-t02)*(t01-t02) + (tb1-tb2)*(tb1-tb2)) / ((t01)*(t01) + (tb1)*(tb1));
	dks0=1.0/t01;
	dks1=1.0/t02;

	if(dr < 1.0e-10){
		if(dks0 > dks1)
			dks0=dks1*1.00001;
		else
			dks0=dks1*0.99999;
	}

	dksm=dks1-dks0;

	ds=1.0/ts;
	dd=1.0/td;

	dw0=1.0/tb1;
	dw1=1.0/tb2;
	dwp=dw0+dw1;
	dwm=dw1-dw0;

	dksm2=dksm*dksm;

	dzna=(dksm2+dwm*dwm)*(dksm2+dwp*dwp);


	das0=dw1*(dksm2+dwp*dwm);
	dac0=-2*dksm*dw0*dw1;
	das1=dw0*(dksm2-dwp*dwm);
	dac1=-dac0;


	dsn0=(ds-dks0);
	dcs0=-dw0;
	dzn0=dcs0*dcs0+dsn0*dsn0;

	dsn0s=(dsn0*das0-dcs0*dac0)/dzn0;
	dcs0s=(dcs0*das0+dsn0*dac0)/dzn0;

	dsn0=(ds-dks1);
	dcs0=-dw1;
	dzn0=dcs0*dcs0+dsn0*dsn0;

	dsn1s=(dsn0*das1-dcs0*dac1)/dzn0;
	dcs1s=(dcs0*das1+dsn0*dac1)/dzn0;


	dsn0=(dd-dks0);
	dcs0=-dw0;
	dzn0=dcs0*dcs0+dsn0*dsn0;

	dsn0d=(dsn0*das0-dcs0*dac0)/dzn0;
	dcs0d=(dcs0*das0+dsn0*dac0)/dzn0;

	dsn0=(dd-dks1);
	dcs0=-dw1;
	dzn0=dcs0*dcs0+dsn0*dsn0;

	dsn1d=(dsn0*das1-dcs0*dac1)/dzn0;
	dcs1d=(dcs0*das1+dsn0*dac1)/dzn0;

	dr=(ts-td)/td;


	
	rv=( ((dsn0s-dsn0d)*sin(dw0*t)	+ (dcs0s-dcs0d)*cos(dw0*t)) 
	* exp(-t*dks0) - (dcs0s+dcs1s)*exp(-t*ds)+(dcs0d+dcs1d)*exp(-t*dd)
	+ ((dsn1s-dsn1d)*sin(dw1*t)	+ (dcs1s-dcs1d)*cos(dw1*t))*exp(-t*dks1) )
	/dzna/(ts-td);
	
	
	rv=rv/(-.109+.919*t01-.261*t01*t01)
	/(-.109+.919*t02-.261*t02*t02)
		/(.262+.174*tb1-.208*tb1*tb1)
		/(.262+.174*tb2-.208*tb2*tb2)
		/(4.56-1.58*td1)/(1.391-0.434*ts1)
		/(1.06-0.578*(t01-tb1)*(t01-tb1))
		/(1.06-0.578*(t02-tb2)*(t02-tb2))
		/(1.2140-0.79645*t01+0.63440*t01*t01)
		/(1.2140-0.79645*t02+0.63440*t02*t02);
	
	
return rv;

}

double sv101(double *xx, double *FITPAD){

	double tr1 = xx[0];
	//double precision FITPAD(12)

	double FITFUN;

	//      common/norm/ped,amp,ts0,td,t0,b1,t1,amm,tmm,t01,b2,a
	double ped,amp,ts0,td,t0,b1,t1,amm,tmm,t01,b2,a;
	double tr,x;
	double tr2, tr3;

	
ped = FITPAD[0];
amp = FITPAD[1];
ts0 = FITPAD[2];
td  = FITPAD[3];
t0  = FITPAD[4];
b1  = FITPAD[5];
t1  = FITPAD[6];

amm = FITPAD[7];
tmm = FITPAD[8];

t01 = FITPAD[9];
b2  = FITPAD[10];
a   = FITPAD[11];


	tr=tr1-ts0;
	tr2=tr+0.2;
 
	tr3=tr-0.2;


	if(tr2 <= 0) return (ped ); 

		FITFUN=( sv123(tr,t0,b1,t01,b2,td,t1)*(1-a) + 
				a*0.5*(sv123(tr2,t0,b1,t01,b2,td,t1) + sv123(tr3,t0,b1,t01,b2,td,t1)) );
			x=tr/t0;

		FITFUN=amp*(FITFUN-amm*(exp(-tr/tmm)*(1-exp(-x)*(1+x+x*x/2+x*x*x/6+x*x*x*x/24+x*x*x*x*x/120))));
	

	FITFUN = FITFUN + ped;

	return FITFUN;
}


*/
double *cholesky(double *A, int n) {
    double *L = (double*)calloc(n * n, sizeof(double));
    if (L == NULL)
        exit(EXIT_FAILURE);
 
    for (int i = 0; i < n; i++){
        for (int j = 0; j < (i+1); j++) {
            double s = 0;
            for (int k = 0; k < j; k++)
                s += L[i * n + k] * L[j * n + k];
				
            L[i * n + j] = (i == j) ?
                           sqrt(A[i * n + i] - s) :
                           (1.0 / L[j * n + j] * (A[i * n + j] - s));
		
        }
	}
 
    return L;
}

void show_matrix(double *A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%2.5f ", A[i * n + j]);
        printf("\n");
    }
}

double distr(double a, double thr, double rand){
	return thr - a*log(rand);
}

void electronic_noise(Double_t *X, char *fname, TRandom *random){

	FILE *PR;
  	PR = fopen(fname, "r");

   	if ( PR == NULL){
		std::cout <<"file "<<fname<< " is not found" <<std::endl;
     	exit (1);
    }


	Int_t iw=0;
	Int_t ih;
	double a_read;
//	double m31[31*31];
	Int_t el = 0;

	Eigen::MatrixXd m31(31,31);
	while(ih!=31){
		fscanf (PR, "%d", &ih);
		while (iw!=31)
		{
			fscanf (PR, "%d %lf", &iw,&a_read);
			//m31[el] = a_read;
			m31(ih-1, el - 31*(ih-1))=a_read;
			el++;
		}
		iw=0;
	}
	fclose(PR);  
	
	
	int npoint = 31;
	//double *Lcosm = cholesky(m31, npoint);
	//show_matrix(Lcosm, 31);

	Eigen::MatrixXd Lcosm  = m31.llt().matrixL();
	//std::cout<<Lcosm<<std::endl;

	// 31 yi
	Eigen::VectorXd Y(31);
	Double_t min=-20;
	Double_t max=20;

		
	for(int i=0;i<npoint;i++)
		Y(i)=random->Gaus();

	for(int j=0;j<npoint;j++)
			X[j]=0;

	//electronics noise
	for(int i=0;i<npoint;i++){
		for(int j=0;j<npoint;j++){
			X[i]+=Lcosm(i,j)*Y(j);
	
	}}

}


double pile_up_noise(double t, TF1 *waveform, TRandom *random, double coeff){
	//Fsr = 1.76MHz
	// pile-up frequency 8 событий в мкс
	Double_t pileupf = 0.8*coeff;
	//Double_t pileupf = coeff;

	Double_t Tmax = 8.0;
	Double_t Tmin = -32.0;
	Double_t Fnoise=0.0;

	Double_t N =(Tmax-Tmin)*pileupf; 
//	TH1F *h1 = new TH1F("pile-up", "pile-up", 100, 0, 150);

		Int_t n = (Int_t)random->Gaus(N,sqrt(N));
		Double_t ti[n];
		Double_t A[n];
		double a = 14.0;
		double thr = 1.0;

		
		for(int i = 0; i<n;i++){
			ti[i] = random->Uniform(Tmin, Tmax);
			A[i] = distr(a, thr, random->Uniform());
			//h1->Fill(A[i]);
			Fnoise+=A[i]*waveform->Eval((Double_t)t-ti[i]);
			//printf("%lf,  %lf\n", Fnoise, t);
		}


	//h1->Draw();

	return Fnoise;

}

void fillparam(char *fname, double paq[][12]){
	Int_t u;
	Double_t c0;
	Double_t c1;
	Double_t c2;
	Double_t c3;
	Double_t c4;
	Double_t c5;
	Double_t c6;
	Double_t c7;
	Double_t c8;
	Double_t c9;
	Double_t Cs[10][16];

	FILE *PR;
  	PR = fopen(fname, "r");

   	if ( PR == NULL){
		std::cout <<"file "<<fname<< " is not found" <<std::endl;
     	exit (1);
    }

	fscanf(PR,"%lf  \n",&c0);

	while(!feof(PR)){
		fscanf (PR, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf    ", &u,&c0,&c1,&c2,&c3,&c4,&c5,&c6,&c7,&c8,&c9);
	//	printf ("%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n   ",u,c0,c1,c2,c3,c4,c5,c6,c7,c8,c9);
		if(u>0){
		      Cs[0][u-1]=c0;
		      Cs[1][u-1]=c1;
		      Cs[2][u-1]=c2;
		      Cs[3][u-1]=c3;
		      Cs[4][u-1]=c4;
		      Cs[5][u-1]=c5;
		      Cs[6][u-1]=c6;
		      Cs[7][u-1]=c7;
		      Cs[8][u-1]=c8;
		      Cs[9][u-1]=c9;
		    }
	}

	fclose(PR);  

	Int_t Chst =1;
    Int_t Ched =16;
	Int_t ih;

	for(ih=Chst;ih<=Ched;ih++){  

	//	printf("ih=%d \n ",ih);

		paq[ih-1][0]  = -0.000; 
		paq[ih-1][1]  = 1.00; 
		//paq[2]  = ((double)poq/72.)-1.-2.07;
	/*	paq[ih-1][2]  = 0.0;
		paq[ih-1][3]  = Cs[1][ih-1]; 
		paq[ih-1][4]  = Cs[2][ih-1]; 
		paq[ih-1][5]  = Cs[3][ih-1]; 
		paq[ih-1][6]  = Cs[4][ih-1]; 
		paq[ih-1][7]  = Cs[5][ih-1]; 
		paq[ih-1][8]  = Cs[6][ih-1]; 
		paq[ih-1][9]  = Cs[7][ih-1]; 
		paq[ih-1][10] = Cs[8][ih-1]; 
		paq[ih-1][11] = 0; */

		paq[ih-1][2]  = 0.75;
		paq[ih-1][3]  = 0.648324; 
		paq[ih-1][4]  = 0.401711; 
		paq[ih-1][5]  = 0.374167; 
		paq[ih-1][6]  = 0.849417; 
		paq[ih-1][7]  = 0.00144548; 
		paq[ih-1][8]  = 4.70722; 
		paq[ih-1][9]  = 0.815639; 
		paq[ih-1][10] = 0.555605;
		paq[ih-1][11] = 0.2752; 
		
	}
}

void lftda_(short int *id, int *pf, int *pf1, int *hf, int *hf1, int *fg41, int *fg43, int *fg45, int *fg51, int *fg52, int *fg53, int *fg54, int *fg55,
          int *y,int *ttrig, int *n16, int *ch,int *mi5,int *p_lar,int *ltr,int *h_lar,int *lq,int *lch3,int *lch4,int *lcp, int *ask, int *yfit, int *dif)
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
  int k_d = (int)*((unsigned char*)id+29);
  int k_e = (int)*((unsigned char*)id+30);
  int k_16 = (int)*((unsigned char*)id+31);
 // printf("ka=%d, kb=%d, kc=%d, kd=%d, ke=%d \n", k_a, k_b, k_c, k_d, k_e);
  
  int k1_chi =  (int)*((unsigned char*)id+34);
  int k2_chi =  (int)*((unsigned char*)id+35);
  int chi_thres = (int)*(id+16);
  long long int s1,s2,s3,s4;

 
  long long int z0; 
  int it, it0;
  double itd, norm; 
 // long long d_it; 
  int it_h, it_l; 
  long long int A1, B1, A2,A3,A4,C1,ch1,ch2,B2h,B3h,B5h,B2,B3,B5,B6; 
  int low_ampl,i,T,iter;   

  *ask=Ahard;
  if(k_16+*n16 !=16)
    {
      printf("disagreement in number of the points %d %d \n",k_16,*n16);
    }


  int validity_code=0;
  for(i=0, z0=0; i<16; i++)
    z0 += y[i];

 /* printf("ttrig = %d\n", *ttrig);

	TCanvas *mlpa_canvas = new TCanvas("mlpa_canvas","",1000,700);
	TGraph *handp = new TGraph(192);
	TGraph *hComp = new TGraph(192);
	TGraph *pComp = new TGraph(192);*/


 // for(int  d =-10; d<=10; d++){
	 // std::cout<<"-----------------------------------"<<std::endl;
	it0 = 48 + ((143-*ttrig)*2/3);//+d;

//	printf(" it0 = %d\n", it0);

  //limits
  it_h=191;
  it_l=0;

  if (it0 < it_l)it0=it_l; 
  if (it0 > it_h)it0=it_h; 
  it=it0;

  T = it0<<4;

//first approximation without time correction

  int ttrig0=*ttrig/6;

  s1=(*(fg41+ttrig0*16));
  s2=(*(fg43+ttrig0*16));
  A3 = (s1 * z0);
  A4 = (s2 * z0);


  for(i=1;i<16;i++){
    s1=(*(fg41+ttrig0*16+i));
	s2=(*(fg43+ttrig0*16+i));
	
    B3 = y[15+i];
    B3 = s1 * B3;
    A3 += B3;

	B3 = y[15+i];
    B3 = s2 * B3;
    A4 += B3;

  }

  A3 += (1<<(k_a-1));
  A3 >>= k_a;

  A4 += (1<<(k_c-1));
  A4 >>= k_c;


  s1 = (*(fg45+ttrig0*16));
  C1=s1*z0;

for(i=1;i<16;i++){
	s1=*(fg45+i+ttrig0*16);
	B5=y[15+i];
	B5=s1*B5;
	C1 += B5;
}
C1 += (1 << (k_e-1));
C1 >>= k_e;

//for(it=0; it<192;it++){
for(iter=0, it=it0; iter<3;){
    iter++;


	s1=(*(fg51+it*16));
	s2=(*(fg53+it*16));
	s3=(*(fg52+it*16));
	s4=(*(fg54+it*16));
	//printf("%4d, ", it);
	

	A1 = (s1 * z0);
	A2 = (s2 * z0);
	B1 = (s3 * z0);
	B6 = (s4 * z0);

	for(i=1;i<16;i++){
		s1=(*(fg51+it*16+i));
		s2=(*(fg53+it*16+i));
		s3=(*(fg52+it*16+i));
		s4=(*(fg54+it*16+i));

		B3 = y[15+i];
		B3 = s1 * B3;
		A1 += B3;

		B3 = y[15+i];
		B3 = s2 * B3;
		A2 += B3;

		B3=y[15+i];
        B3= s3 * B3;
        B1 += B3;

		B3=y[15+i];
        B3= s4 * B3;
        B6 += B3;
		
	}

	A1 += (1<<(k_a-1));
	A1 = A1>>k_a;

	A2 += (1<<(k_c-1));
	A2 = A2>>k_c;

	if(iter != 3){
		//printf("iter=%d \n",iter);
        B2 = B1>>(k_b-9);
        B1 = B2>>9;
        B2 += (A1<<9);
        B3=(B2/A1);
		B3=((B3+1)>>1)-256;
	
		B2h = B6>>(k_d-9);
        B6 = B2h>>9;
        B2h += (A2<<9);
        if (A2 !=0){
			B3h=(B2h/A2);
			B3h=((B3h+1)>>1)-256;
		}
		//std::cout<<"p_A="<<A1<<", p_t="<< B3<<std::endl;
		//std::cout<<"h_A="<<A2<<", h_t="<<B3h<<std::endl;
	
		itd = (double) (B3*A1 + 1.3*B3h*A2);
		norm = (double) (A1 + A2);
		itd = itd / norm;
	
		B3=itd;
		it += B3;
	
        it = it>it_h ? it_h : it;
        it = it<it_l ? it_l : it;
    } 
	else{
		//printf("iter=%d \n",iter);
		B2 =  B1>>(k_b-13);
        B5 =  B1>>(k_b-9);
        B2h = B6>>(k_d-13);
        B5h = B6>>(k_d-9);
        
        B1 = B2>>13;
        B2 += (A1<<13);
        B3=(B2/A1);
		B3=((B3+1)>>1)-4096;

		B6 = B2h>>13;
        B2h += (A2<<13);
		if(A2!=0){
        	B3h=(B2h/A2);
			B3h=((B3h+1)>>1)-4096;
		}
		else B3h=0;
		itd = (double) (B3*A1 + 1.3*B3h*A2);
		norm = (double) (A1 + A2);
		itd = itd / norm;
	
		B3=itd;
        T = ((it)<<4) + B3;
        T = T > 3071 ?  3071 : T;
        T = T < 0 ? 0 : T;

        B1=B5>>9;
		//hComp->SetPoint(it, it, B1);
		//printf("B=%lld, ", B1);
        B5 += (A1<<9);
        B3=(B5/A1);
		B3=((B3+1)>>1)-256;
		//printf("A=%lld, 256*B/A=%lld, ", A1, B3);
		//pComp->SetPoint(it, it, A1);
		//handp->SetPoint(it, it, B3<<8);
		//printf("dt=%lld\n", B3);

		B6=B5h>>9;
        B5h += (A2<<9);
		if(B3h!=0){
       	 	B3h=(B5h/A2);
			B3h=((B3h+1)>>1)-256;
		} 
		else B3h=0;
		
		itd = (double) (B3*A1 + 1.3*B3h*A2);
		norm = (double) (A1 + A2);
		itd = itd / norm;
		
		B3=itd;
		//printf("dt=%lld\n", B3);
		//handp->SetPoint(it, it, B3);

		*dif=B3;

        it +=B3;
		//printf("dt=%lld, it=%d\n", B3, it);
        it = it>it_h ? it_h : it; 
        it = it<it_l ? it_l : it;


        s1 = (*(fg55+it*16));
		C1=s1*z0;
		for(i=1;i<16;i++){
			s1=*(fg55+i+it*16);
			B5=y[15+i];
			B5=s1*B5;
			C1 += B5;
		}
		C1 += (1 << (k_e-1));
		C1 >>= k_e;
    }
}

/*handp->Draw();
hComp->SetLineColor(kRed);
hComp->Draw("same");
pComp->SetLineColor(kBlue);
pComp->Draw("same");

TLegend *legend = new TLegend(.75, .80, .95, .95); 
legend->AddEntry(handp, "B/A * 2^16");
legend->AddEntry(hComp, "B");
legend->AddEntry(pComp, "A");
legend->Draw();
*/
	ou:
	*mi5=it;
	*p_lar=A1;
	*h_lar=A2;
	*ltr=T;

	int ss=(y[20]+y[21]);

	if(ss <= Ahard)validity_code=validity_code+4;


	*lq=validity_code;
	*lch3=B2;

	*lch4=ch1;

	*lcp=C1;

	return ;
}





int main(){

	double paq[16][12];
	char in[256];
	strcpy(in, std::getenv("HOME"));
	strcat(in, "/fit/panr_11.dat");

	//fillparam(in, paq);
 
	int chan = 15;
	TF1 *photonFun = new TF1("fitfun3",&ShaperDSP_F,0.,9.0,12); 


	double MZ[11] = {0.75, 0.648324, 0.401711, 0.374167, 0.849417, 0.00144548, 4.70722, 0.815639, 0.555605, 0.2752, 1.0};
 	double spt=0;
 	double fm, fc;
	double t1[1]={0};

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

	spt=2.*spt/1000.;
    t1[0]=0.;
    fm=0.;

    for (int j1 = 0; j1 < 4000; j1++) {
      fc=ShaperDSP_F(t1, MZ);
      if(fc>fm){
        fm=fc;
      }
      t1[0]=t1[0]+spt;
    }

    MZ[10]=1.0/fm;

	photonFun->SetParameters(MZ);



	TF1 *hadronFun = new TF1("fitfun3",&ShaperDSP_F,0.,9.0,12); 


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
    MZ[10]=1;
 	spt=0;
 
	t1[0]=0.;
    fm=0.;

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

	spt=2.*spt/1000.;
    t1[0]=0.;
    fm=0.;

    for (int j1 = 0; j1 < 4000; j1++) {
      fc=ShaperDSP_F(t1, MZ);
      if(fc>fm){
        fm=fc;
      }
      t1[0]=t1[0]+spt;
    }

    MZ[10]=1.0/fm;

	hadronFun->SetParameters(MZ);

	//fitfun3->Draw()

	TRandom *random = new TRandom();
	random->SetSeed((ULong_t)time(NULL));


double F;
double Fk;
double Ffulk;


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



//пьедистал 4000
//18 разрядов

TTree tree("tree1","waveforms1");
TTree ftree("ftree","editeddata");

struct waveformData{ 
        Int_t y[31]; 
		//Int_t yfit[31]; 
        Float_t time[31]; 
		Int_t Aset;
		Int_t hAset;
		Int_t p_Afit;
		Int_t h_Afit;
		Float_t ts;
		Int_t tfit;
		Float_t mi5;
		Int_t Pset;
		Int_t Pfit;
		Int_t dif;
};


waveformData mywaveform;

tree.Branch("waveform_ampl", &mywaveform, "y[31]/I:time[31]/F:Aset/I:hAset/I:p_Afit/I:h_Afit/I:ts/F:tfit/I:mi5/F:Pset/I:Pfit/I:dif/I");

ftree.Branch("waveform_ampl", &mywaveform, "y[31]/I:time[31]/F:Aset/I:hAset/I:p_Afit/I:h_Afit/I:ts/F:tfit/I:mi5/F:Pset/I:Pfit/I:dif/I");

sprintf(in,"matr31.dat");

//TH1F *h1 = new TH1F("elnoise", "elnoise", 100, -100, 100);
double mean = 0;

for(int k=0; k<5000; k++){
	double ts = random->Uniform(-0.75, -0.26);//*192;
	//ts = (int)ts;
	//ts=ts/192.;
	//std::cout<<ts<<std::endl;
	double Aset = 5000;//random->Uniform(1000, 100000);
	double Aset_h = Aset*0.3;//random->Uniform(0, 0.3);
	Double_t X[31]; //X[i] electronic noise
	electronic_noise(X, in, random); 
	
	//double ts = -0.1;
	//double Aset = 10000;  //amplitude

	//TGraph *pilegph = new TGraph(31);
	double y31[31];   //31 signal points                              
	int points = 15;
	int Pset=2000; //pedestal
	for(double t = -16*0.5, i=0; t<points*0.5; t+=0.5,i++){
		F = pile_up_noise(t, photonFun, random, 1.0); //pile-up noise amplitude at i point
		//std::cout<<fitfun3->Eval(t-ts)<<",   "<<fitfun1->Eval(t-ts)<<std::endl;
		//pilegph->SetPoint(s-1, t,F);
		Fk = X[(int)i] + F;  //full noise amplitude
		//h1->Fill(F);
		Ffulk = Fk + Aset*photonFun->Eval(t-ts) + Aset_h*hadronFun->Eval(t-ts);//+Pset;
		//pilegph->SetPoint((int)i, t, Ffulk);
		y31[(int)i]=Ffulk;
		mywaveform.time[(int)i] = t;
	}
	//pilegph->Draw("AL*");


	int yint31[31];
	for (int i = 0; i<31; i++){
		yint31[i]=lrint(y31[i]);
		mywaveform.y[i] = yint31[i];
	}

    //    48+(143-ttrig)*2/3 = -1 - 192*ts

	int ttrig = 218+288*ts;  //129+288*ts; //trigger time
	int n16 = 16;
	int p_lar, ltr, mi5, h_lar, lq, lch3, lch4, lcp, ask, yfit[31], dif;

	lftda_(id, cpf[0], cpf1[0], chf[0], chf1[0], cfg41[0], cfg43[0], cfg45[0], cfg51[0], cfg52[0], cfg53[0], cfg54[0], cfg55[0], yint31, 
		&ttrig, &n16, &chan, &mi5, &p_lar, &ltr, &h_lar, &lq, &lch3, &lch4, &lcp, &ask, yfit, &dif);


	//printf("mi5=%d, p_lar=%d, ltr=%d, lcp=%d\n", mi5, p_lar, ltr, lcp);

//	for(int i = 0; i<31; i++) mywaveform.yfit[i]=yfit[i]; 
    
	mywaveform.Aset=Aset;
	mywaveform.hAset=Aset_h;
	mywaveform.p_Afit=p_lar;
	mywaveform.h_Afit=h_lar;
	mywaveform.ts=ts;
	mywaveform.tfit=ltr;
	mywaveform.mi5=mi5;
	mywaveform.Pset=Pset;
	mywaveform.Pfit=lcp;
	mywaveform.dif=dif;

	tree.Fill();

	float sigma2 = 46;
	float coeff = 0.8699;

	if(p_lar < -(int)sigma2){
		float sum = (float)p_lar + coeff*(float)h_lar;
		sum = (sum + sigma2) / coeff;
		mywaveform.h_Afit = (int)sum;
		mywaveform.p_Afit = -(int)sigma2;
	}

	if(h_lar < -(int)sigma2){
		float sum = (float)p_lar + coeff*(float)h_lar;
		//printf("sum=%f, p_lar=%d, h_lar=%d, ", sum, p_lar, h_lar);
		sum = sum + coeff*sigma2;
		mywaveform.p_Afit = (int)sum;
		mywaveform.h_Afit = -(int)sigma2;
		//printf("newsum=%f, h_Afit=%d, p_Afit=%d\n", (float)mywaveform.p_Afit+coeff*(float)mywaveform.h_Afit, mywaveform.h_Afit, mywaveform.p_Afit);
	}

	ftree.Fill();
}

TFile fout("output1.root","recreate");
tree.Write();
fout.Close();

TFile fout_ed("fitout.root","recreate");
ftree.Write();
fout_ed.Close();

//h1->Draw();
}
