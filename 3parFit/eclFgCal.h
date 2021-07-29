// This class emulates algorithm to fit the signal shape  implimented in ECL shaperDSP boards
//
//#include "TH1F.h"
//#include "TFile.h"
//#define NCHANNELS 16

/*
struct ecldspConfig{
	int A0; 
	int Askip;

	int Ahard;
	int k_a; 
	int k_b;
	int k_c;
	int k_16;
	int k1_chi;
	int k2_chi;
}
*/
#pragma once

class eclFgCal{
	private:

	/*	TH1F *Kf;
		TH1F *Ka;
		TH1F *Kb;
		TH1F *Kc;
		TH1F *uf[192][16];
		TH1F *uf1[192][16];
		TH1F *uf31[192][16];
		TH1F *uf32[192][16];
		TH1F *uf33[192][16];
		TH1F *uf41[24][16];
		TH1F *uf43[24][16]; 
		*/


                double IY[16][16][16];  //M_{ij}^{-1}  [channel][i][j] element of inverse matrix 

                double par[16][10];  // response function parametrization

                double f[16][192][16];  // fitting function parametrization
                double f1[16][192][16];  //first derivative of function parametrization
                double fg31[16][192][16];  // array for amplitude reconstrustion [channel][i][j]
                double fg32[16][192][16];  // array for time reconstrustion [channel][i][j]
                double fg33[16][192][16];  // array for pedestal reconstrustion [channel][i][j]

                double fg41[16][24][16];  // array for amplitude reconstrustion [channel][i][j]

                double fg43[16][24][16];  // array for pedestal reconstrustion [channel][i][j]



                double mf[16][4];  // maximum value of the fitting function parametrization
                double mf1[16][4];  //maximum value of the first derivative of function parametrization
                double mfg31[16][4];  // maximum value of the array for amplitude reconstrustion [channel][i][j]
                double mfg32[16][4];  // maximum value of the array for time reconstrustion [channel][i][j]
                double mfg33[16][4];  // maximum value of the array for pedestal reconstrustion [channel][i][j]

                double mfg41[16][4];  // maximum value of the array for amplitude reconstrustion [channel][i][j]

                double mfg43[16][4];  // maximum value of the array for pedestal reconstrustion [channel][i][j]







                short  int m_f[16][192][16];  // fitting function parametrization
                short  int m_f1[16][192][16];  //first derivative of function parametrization
                short  int m_fg31[16][192][16];  // int  array for amplitude reconstrustion [channel][i][j]
                short  int m_fg32[16][192][16];  // int array for time reconstrustion [channel][i][j]
                short  int m_fg33[16][192][16];  // int array for pedestal reconstrustion [channel][i][j]

                short  int m_fg41[16][24][16];  // int array for amplitude reconstrustion [channel][i][j]
                short  int m_fg43[16][24][16];  // int array for pedestal reconstrustion [channel][i][j]

               short int    bf[16];  //  number of bit for f parametrizition
               short int   bf1[16];  //  number of bit for f1 parametrizition
               short int bfg31[16];  //  number of bit for fg31 parametrizition
               short int bfg32[16];  //  number of bit for fg32 parametrizition
               short int bfg33[16];  //  number of bit for fg33 parametrizition
               short int bfg41[16];  //  number of bit for fg41 parametrizition
               short int bfg43[16];  //  number of bit for fg43 parametrizition
               short int k2_chi;  //  bit quality calculation parameter
               short int chi_thres;  //  chi-squared threshold


                char hsnm[256]; // name historgamm
                int cr; //
                int wr_rootf; //


	//	bool mask_[16];

	public:
		
		eclFgCal(/*char* name, int cr*/ ); // constructor
		~eclFgCal(); //destructor
		
		void GetResponsePar(const char *dir);  // read parameters for response function

		//		void ReadInverseMatrices(const char *dir);  // read inverse covariance matrices
		void CalculateFgpar(); // calculate parameters for time and enegry reconctrustion f f1 fg31 fg32 fg33 fg41 fg43
		//		void CalculateFgparF(); // calculate parameters for time and enegry reconctrustion fg31 fg32 fg33 fg41 fg43	
		//        	void CalculateShape(); // calculate function shape f f1 	
        void CalculateBit(); // calculate bit size for fg

        void ReadChiCut(const char *dir); // read k2_chi and threshold

		void CalculateFginInt(); // calculate  fg in integer

        void WriteCalib(const char *dir); // write calibration in rec.dsp

		//          	void ReadCalib(const char *dir); //  read calibration from rec.dsp

        void RootWr(); //  write calibration from rec.dsp

      //  void RootFill( TFile *rootfile); //  write calibration from rec.dsp

		//          	void RootFuncFill( TFile *rootfile); //  write function shape f f1 in root file
		//          	void RootFuncRead( const char *treeroot); //  read function shape f f1 from root file

		bool setIMElement(int chan, int i, int j, double val) {
		  if ((0<chan) && (chan<=16) && (0<=i) && (i<16) && (0<=j) && (j<16)) {
		    IY[chan-1][i][j]=val;
		   // mask_[chan-1]=true;
		    return true;
		  }
		  return false;
		};

		/** Shaper-DSP Output  Map */
	//	double ShaperDSP(double Ti);

		/** Shaper-DSP Output  Map */
	//	double ShaperDSP_F(double Ti, double* ss);

		/** Shaper-DSP Output  basic */
	//	double  Sv123(double t, double t01, double tb1, double t02, double tb2, double td1, double ts1);

		void SaveInverseMatrices(const char *dir) ;

};

/*
class shaperDspCalculator{
	private:
		TMatrix *FG41;
		// ...

	public:
		makeCovarianceMatrix(char *pedestalFileName);
		calculateCoefficients(TF1 *signalShape, int *par);
		saveMatricesToFile(char *coefsFileName);
};
*/
