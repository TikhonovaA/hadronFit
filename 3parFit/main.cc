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
	//double m31[31*31];
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

	Eigen::MatrixXd Lcosm  = m31.llt().matrixL();
	// 31 yi
	//Double_t Y[31];
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
			//X[i]+=Lcosm[i*npoint+j]*Y[j];
	
	}}

}


double pile_up_noise(double t, TF1 *waveform, TRandom *random, double coeff){
	//Fsr = 1.76MHz
	// pile-up frequency 8 ?????????????? ?? ??????
	Double_t pileupf = 0.8*coeff;
	//Double_t pileupf = coeff;

	Double_t Tmax = 8.0;
	Double_t Tmin = -32.0;
	Double_t Fnoise=0.0;

	Double_t N =(Tmax-Tmin)*pileupf; 
	//TH1F *h1 = new TH1F("pile-up", "pile-up", 100, 0, 150);


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

		paq[ih-1][0]  = -0.000; 
		paq[ih-1][1]  = 1.00; 
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


void lftda_(short int *id,short int *f,short int *f1,short int *fg41,short int *fg43,short int *fg31,short int *fg32,short int *fg33,
          int *y,int *ttrig, int *n16, int *ch,int *mi5,int *lar,int *ltr,int *lq,int *lch3,int *lch4,int *lcp, int *ask/*,int *evna,int *ifp*/, int *yfit)
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
  //printf("ka=%d, kb=%d, kc=%d\n", k_a, k_b, k_c);
  
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


 /* if(*ifp == 1)printf("v lftda _______________%d %d \n",*evna,*ttrig);*/
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
	//printf("%d %d %d\n", s1, k_a, s1>>k_a);
    B3 = y[15+i];
    B3 = s1 * B3;
    A2 += B3;
  }

  A2 += (1<<(k_a-1));
  A2 >>= k_a;

  T = it0<<4;

//printf("A2 = %lld\n", A2);

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

	//printf("\n *** ITERATION #%d ***\n\n", iter );

    //printf(" it = %d\n", it);

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
//double D = B1/pow(2,k_b);
//double DD = D/(double)A1*256;
//	 std::cout<<DD<<std::endl;
	

      if(iter != 3){
        B2 = B1>>(k_b-9);
        B1 = B2>>9;

        B2 += (A1<<9);


        B3=(B2/A1);


	//	printf("iter=%d %d\n", iter, it);
        it +=((B3+1)>>1)-256;
		//printf("%lld, %d\n", ((B3+1)>>1)-256, it);
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
      //  printf("iter=%d %d\n", iter, it);
        it +=((B3+1)>>1)-256;
		//printf("%lld, %d\n", ((B3+1)>>1)-256, it);
        it = it>it_h ? it_h : it; 
        it = it<it_l ? it_l : it;

        T = T > 3071 ?  3071 : T;

        T = T < 0 ? 0 : T;
        C1 = (*(fg33+it*16) * z0);
        for(i=1;i<16;i++)
          C1 +=*(fg33+i+it*16) * y[15+i];
        C1 += (1 << (k_c-1));
        C1 >>= k_c;
		//printf("C1 (high amplitude) = %lld\n", C1);
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
     // printf("C1 (low amplitude) = %lld\n", C1);

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
	  yfit[i-1]= C1;
	  yfit[i+15]= ch2+C1;
      ch2=(y[i+15]-ch2-C1);

      ch1=ch1 + ch2*ch2;

    }
	for(i=0; i<16;i++){
		yfit[i]= C1;
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


int main(){

	/*double paq[16][12];
	char in[256];
	strcpy(in, std::getenv("HOME"));
	strcat(in, "/fit/1fit/panr_11.dat");

	fillparam(in, paq);
 
	int chan = 15;
	TF1 *fitfun3 = new TF1("fitfun3",&sv101,0.,9.0,12); 
	fitfun3->SetParameters(paq[chan]);


	TRandom *random = new TRandom();
	random->SetSeed((ULong_t)time(NULL));
*/

double paq[16][12];
char in[256];
strcpy(in, std::getenv("HOME"));
strcat(in, "/fit/3parFit/panr_11.dat");

//fillparam(in, paq);

int chan = 15;
TF1 *fitfun3 = new TF1("fitfun3",&ShaperDSP_F,0.,9.0,12); 


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

fitfun3->SetParameters(MZ);



	TF1 *fitfun1 = new TF1("fitfun3",&ShaperDSP_F,0.,9.0,12); 


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

	fitfun1->SetParameters(MZ);


TRandom *random = new TRandom();
random->SetSeed((ULong_t)time(NULL));


//TH1F *h3 = new TH1F("pile-up", "pile-up", 100, 0, 200);

double F;
double Fk;
double Ffulk;




//read dsp.dat for one channel
FILE *fr;
int h,nsiz,nsiz1;
char fname[256];
strcpy(fname, std::getenv("HOME"));
strcat(fname, "/fit/3parFit/DSP_exp10/dsp.dat");

short int id[256];
short int cf[192][16];
short int cf1[192][16];
short int cfg31[192][16];
short int cfg32[192][16];
short int cfg33[192][16];
short int cfg41[24][16];
short int cfg43[24][16];

for (int i = 0; i < 256; i++){
	id[i]=0;
}

for ( int u =0 ; u < 16; u ++){
  for ( int k =0 ; k < 192; k ++){
    if(k < 24){
      cfg41[k][u] = 0;
      cfg43[k][u] = 0;
    }
   
    cfg31[k][u] = 0;
    cfg32[k][u] = 0;
    cfg33[k][u] = 0;
    cf[k][u]    = 0;
    cf1[k][u]   = 0;

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
	std::cout << "Error writing id data. Read block size = " << h << '\n';
	exit(0);
}


for(int i=0;i<=chan;i++){

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

}//i

fclose(fr);

//?????????????????? 4000
//18 ????????????????


TTree tree("tree1","waveforms1");

struct waveformData{ 
        Int_t y[31]; 
		Int_t yfit[31]; 
        Float_t time[31]; 
		Int_t Aset;
		Int_t hAset;
		Int_t Afit;
		Float_t ts;
		Int_t tfit;
		Float_t mi5;
		Int_t Pset;
		Int_t Pfit;
		Int_t lch4;
		Int_t lq;
};


waveformData mywaveform;

tree.Branch("waveform_ampl", &mywaveform, "y[31]/I:yfit[31]/I:time[31]/F:Aset/I:hAset/I:Afit/I:ts/F:tfit/I:mi5/F:Pset/I:Pfit/I:lch4/I");

sprintf(in,"matr31.dat");

//TH1F *h1 = new TH1F("elnoise", "elnoise", 100, -100, 100);
double mean = 0;
//double Aset = 200;

for(int k=0; k<5000; k++){
	double ts = random->Uniform(-0.75, -0.26);
	double lvl = 0.0;//random->Uniform(0, 0.3);;
	double Aset = random->Uniform(200, 20000);
	double Aset_h = Aset*lvl;
	Aset = Aset*(1-lvl);
	Double_t X[31];
	electronic_noise(X, in, random);

	//TGraph *pilegph = new TGraph(31);
	double y31[31];                                  //array for lftda func
	int points = 15;
	int Pset=4000;
	for(double t = -16*0.5, i=0; t<points*0.5; t+=0.5,i++){
		F=pile_up_noise(t, fitfun3, random, 1.0);
		//pilegph->SetPoint(s-1, t,F);
		Fk = X[(int)i] + F;
		//h1->Fill(Fk);
		//std::cout<<X[(int)i]<<std::endl;
		Ffulk = Fk + Aset*fitfun3->Eval(t-ts) + Aset_h*fitfun1->Eval(t-ts);// +Pset;
		//pilegph->SetPoint((int)i, t, Ffulk);
		y31[(int)i]=Ffulk;
		mywaveform.time[(int)i] = t;
		//printf("%lf\n", X[m]);
	}
	//pilegph->Draw("AL*");


	int yint31[31];
	for (int i = 0; i<31; i++){
		yint31[i]=lrint(y31[i]);
		mywaveform.y[i] = yint31[i];
	}


	int ttrig = 0;
	int n16 = 16;
	int mi5, lar, ltr, lq, lch3, lch4, lcp, ask, yfit[31];

	lftda_(id, cf[0], cf1[0], cfg41[0], cfg43[0], cfg31[0], cfg32[0], cfg33[0], yint31, 
		&ttrig, &n16, &chan, &mi5, &lar, &ltr, &lq, &lch3, &lch4, &lcp, &ask, yfit);

	//printf("mi5=%d, lar=%d, ltr=%d, lq=%d, lch3=%d, lch4=%d, lcp=%d\n", mi5, lar, ltr, lq, lch3, lch4, lcp);

	for(int i = 0; i<31; i++) mywaveform.yfit[i]=yfit[i]; 
    
	mywaveform.Aset=Aset;
	mywaveform.hAset=Aset_h;
	mywaveform.Afit=lar;
	mywaveform.ts=ts;
	mywaveform.tfit=ltr;
	mywaveform.mi5=mi5;
	mywaveform.Pset=Pset;
	mywaveform.Pfit=lcp;
	mywaveform.lch4=lch4;
	mywaveform.lq=lq;


	//mean+=ltr+3072*ts;

	tree.Fill();

}

//mean=mean/5000.0;

//char meanstr[256];
//sprintf (meanstr, "Aset*(tfit+3072*ts-%lf)", mean);
//  printf ("[%s]", meanstr);

//h1->Draw();

//tree.Draw(meanstr);


TFile fout("output1.root","recreate");
tree.Write();
fout.Close();


}

//???????????? ????????, ?????????? +-0.5 ??????, ?? ???????????? ???????????????????? ?????????? 31, ??????????, ???????????? ??????????, ????????????????, ????????, ????????, ?????????????????????? ???????? ???? ??????, 
//???? ??????????????, ?????????????????? ???????????????????? ?????????? - ????????????, ?????????????????? ????????????????.
//16:30