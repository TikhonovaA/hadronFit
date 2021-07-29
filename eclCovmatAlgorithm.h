/**************************************************************************
 * BASF2 (Belle Analysis Framework 2)                                     *
 * Copyright(C) 2018 - Belle II Collaboration                             *
 *                                                                        *
 * Calculate time offsets for cosmics using chi2 minimization algorigthm. *
 *                                                                        *
 * Author: The Belle II Collaboration                                     *
 * Contributors: Mikhail Remnev                                           *
 *                                                                        *
 * This software is provided "as is" without any warranty.                *
 **************************************************************************/

#pragma once
//#include <calibration/CalibrationAlgorithm.h>
//#include "TH1F.h"
#include "eclFgCal.h"
#include <string>
//class eclFgCal;


namespace Belle2 {

  /** Calibrate ecl crystals using muon pair events */
  class eclCovmatAlgorithm{
  public:

    /**..Constructor */
    eclCovmatAlgorithm();

    void calibrate();
    /**..Destructor */
    virtual ~eclCovmatAlgorithm() {}

    /*** Parameters ***/

    int cellIDLo; /**< Fit crystals with cellID0 in the inclusive range [cellIDLo,cellIDHi] */
    int cellIDHi; /**< Fit crystals with cellID0 in the inclusive range [cellIDLo,cellIDHi] */
    int maxIterations; /**< Fit is always stopped after maxIterations */

    bool debugOutput; /**< Save every histogram and fitted function to debugFilename */
    /** Name of file with debug output, eclBhabhaTAlgorithm.root by default */
    std::string debugFilename;


  private:
    //int eclCrateNumber(int cpr, int fin) ;
    //int eclShapersNumber(int cpr_id, int fin) ;
    //    void GetResponsePar( const char *dir/*, const char *dir1*/  );
    //double  Sv123(double t, double t01, double tb1, double t02, double tb2, double td1, double ts1) ;
    //double ShaperDSP_F(double Ti, double* ss) ;
    // void CalculateFgpar( ) ;
    //void CalculateFginInt( );
    //bool channelIsActive(int chan);
    //bool shaperIsActive(int crate, int shaper);
    //    void WriteCalib( Int_t crate, Int_t shaper, const char *dir ) ;
    int myPow(int x, int p)
    {
      if (p == 0) return 1;
      if (p == 1) return x;
      return x * myPow(x, p - 1);
    };
    
  //  static const int m_nchannels=8736;
    //    std::vector<ECLCovmatElements*> m_covmatel;
    
  //  std::string outFilename;

    double par[16][10];  // response function parametrization
/*
    TH1F *Kf;
    TH1F *Ka;
    TH1F *Kb;
    TH1F *Kc;
    TH1F *uf[192][m_nchannels];
    TH1F *uf1[192][m_nchannels];
    TH1F *uf31[192][m_nchannels];
    TH1F *uf32[192][m_nchannels];
    TH1F *uf33[192][m_nchannels];
    TH1F *uf41[24][m_nchannels];
    TH1F *uf43[24][m_nchannels];
    
*/
    
    


    // this    
        double f[16][192][16];  // fitting function parametrization
        double f1[16][192][16];  //first derivative of function parametrization

      double fg31[16][192][16];  // array for amplitude reconstrustion [channel][i][j]
      double fg32[16][192][16];  // array for time reconstrustion [channel][i][j]
    
    double fg33[16][192][16];  // array for pedestal reconstrustion [channel][i][j]

    double fg41[16][24][16];  // array for amplitude reconstrustion [channel][i][j]
    
    double fg43[16][24][16];  // array for pedestal reconstrustion [channel][i][j]
    
    
/*
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

*/
    //    double evt[m_nchannels];  //number of entries for each channel
    //TH1F *Ms[m_nchannels];
    //double Mu[m_nchannels][16][16];  //M_{ij} [channel][i][j] element of covariance matrix
    
    double IY[16][16][16];  //M_{ij}^{-1}  [channel][i][j] element of inverse matrix 
    char hsnm[256]; // name historgamm
    //    int cr; //
    int wr_rootf; // 
    
    eclFgCal* fgcal; // [crate][shaper]
    //    Double_t dt[m_nchannels];
    
    //int crate_;
    //int shaper_;

  };
} // namespace Belle2

