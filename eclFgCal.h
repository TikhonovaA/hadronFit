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

                double p_par[16][10];  // response photon function parametrization
				double h_par[16][10];  // response hadron function parametrization

                double p_f[16][192][16];   // fitting photon function parametrization
                double p_f1[16][192][16];  // first derivative of photon function parametrization
				double h_f[16][192][16];   // fitting hadron function parametrization
                double h_f1[16][192][16];  // first derivative of hadron function parametrization
                
				double fg51[16][192][16];  // array for photon amplitude reconstrustion [channel][i][j]
                double fg52[16][192][16];  // array for photon time reconstrustion [channel][i][j]
                double fg53[16][192][16];  // array for hadron amplitude reconstrustion [channel][i][j]
				double fg54[16][192][16];  // array for hadron time reconstrustion [channel][i][j]
 				double fg55[16][192][16];  // array for pedestal reconstrustion [channel][i][j]


                double fg41[16][24][16];  // array for photon amplitude reconstrustion [channel][i][j]
                double fg43[16][24][16];  // array for hadron amplitude reconstrustion [channel][i][j]
                double fg45[16][24][16];  // array for pedestal reconstrustion [channel][i][j]

            

                double mpf[16][4];   // maximum value of the fitting function parametrization
                double mpf1[16][4];  // maximum value of the first derivative of function parametrization
				double mhf[16][4];   // maximum value of the fitting function parametrization
                double mhf1[16][4];  // maximum value of the first derivative of function parametrization
                double mfg51[16][4];  // maximum value of the array for amplitude reconstrustion [channel][i][j]
                double mfg52[16][4];  // maximum value of the array for time reconstrustion [channel][i][j]
                double mfg53[16][4];  // maximum value of the array for pedestal reconstrustion [channel][i][j]
				double mfg54[16][4];
				double mfg55[16][4];

                double mfg41[16][4];  // maximum value of the array for amplitude reconstrustion [channel][i][j]
                double mfg43[16][4];  // maximum value of the array for pedestal reconstrustion [channel][i][j]
				double mfg45[16][4];  // maximum value of the array for pedestal reconstrustion [channel][i][j]







            	int m_pf[16][192][16];    // fitting photon function parametrization
                int m_pf1[16][192][16];   // first derivative of photon function parametrization
				int m_hf[16][192][16];    // fitting hadron function parametrization
                int m_hf1[16][192][16];   // first derivative of hadron function parametrization
                int m_fg51[16][192][16];  // int  array for photon amplitude reconstrustion [channel][i][j]
                int m_fg52[16][192][16];  // int array for photon time reconstrustion [channel][i][j]
                int m_fg53[16][192][16];  // int  array for hadron amplitude reconstrustion [channel][i][j]
                int m_fg54[16][192][16];  // int array for hadron time reconstrustion [channel][i][j]
				int m_fg55[16][192][16];  // int array for pedestal reconstrustion [channel][i][j]

                int m_fg41[16][24][16];  // int array for photon amplitude reconstrustion [channel][i][j]
				int m_fg43[16][24][16];  // int array for hadron amplitude reconstrustion [channel][i][j]
                int m_fg45[16][24][16];  // int array for pedestal reconstrustion [channel][i][j]

               short int bpf[16];    //  number of bit for p_f parametrizition
               short int bpf1[16];   //  number of bit for p_f1 parametrizition
			   short int bhf[16];    //  number of bit for h_f parametrizition
               short int bhf1[16];   //  number of bit for h_f1 parametrizition
               short int bfg51[16];  //  number of bit for fg51 parametrizition
               short int bfg52[16];  //  number of bit for fg52 parametrizition
               short int bfg53[16];  //  number of bit for fg53 parametrizition
			   short int bfg54[16];  //  number of bit for fg53 parametrizition
			   short int bfg55[16];  //  number of bit for fg53 parametrizition

               short int bfg41[16];  //  number of bit for fg41 parametrizition
               short int bfg43[16];  //  number of bit for fg43 parametrizition
			   short int bfg45[16];  //  number of bit for fg45 parametrizition
               short int k2_chi;     //  bit quality calculation parameter
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
		//double ShaperDSP_F(double Ti, double* ss);

		/** Shaper-DSP Output  basic */
		//double  Sv123(double t, double t01, double tb1, double t02, double tb2, double td1, double ts1);

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
