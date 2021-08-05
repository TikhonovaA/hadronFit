#include <eclCovmatAlgorithm.h>
//#include <ecl/dbobjects/ECLCrystalCalib.h>
//#include <ecl/utility/ECLChannelMapper.h>
//#include "../../utility/include/ECLChannelMapper.h"
//#include "ECLCovmatElements.h"


//#include <TFile.h>
#include<fstream> 
#include<iostream> 
#include <cstdlib>
//#include <TF1.h>
//#include <TROOT.h>
//#include <TH2D.h>
//#include <TProfile.h>
//#include <TCL1.h>


#include <Eigen/Dense>

using namespace std;
using namespace Belle2;

eclCovmatAlgorithm::eclCovmatAlgorithm(){

  fgcal=0; //new eclFgCal(name, shaper); ???

  wr_rootf = 0;
  
}



void eclCovmatAlgorithm::calibrate()
{
  /** Put root into batch mode so that we don't try to open a graphics window */
  //gROOT->SetBatch();

  char in[256];
  FILE *PR;
  strcpy(in, std::getenv("HOME"));
	strcat(in, "/fit/3parFit/newmatr31.dat");
  PR = fopen(in, "r");

  if ( PR == NULL){
  std::cout <<"file "<<in<< " is not found" <<std::endl;
    exit (1);
  }

	int i=0;
	int ih;
	double a_read;
	double m31[31*31];
	int el = 0;
  double S1=0;
  double S2;


	while(ih!=31){
		fscanf (PR, "%d", &ih);
		while (i!=31)
		{
			fscanf (PR, "%d %lf", &i,&a_read);
			m31[el] = a_read;
			el++;
		}
		i=0;
	}
	fclose(PR);  

  Eigen::MatrixXd covmatrix(16, 16);  

  for(int i=0; i<16;i++){
    for(int j=0; j<16;j++){
      S1+=m31[i*31+j];
    }
  }

  covmatrix(0,0) = S1/256.;

  for(int k=0; k<16; k++){
    for(int m=1; m<16; m++){
      if(k==0){
        S2=0;
        for(int i=0; i<16;i++){
          S2+=m31[i*31+m+15];
        }
        covmatrix(0, m)=S2/16.;
        covmatrix(m, 0)=covmatrix(0, m);
      }
      else
      {
        covmatrix(k, m)=m31[(k+15)*31+m+15];
      }
    }
  }


  Eigen::MatrixXd covmatrixI = covmatrix.llt().solve(Eigen::MatrixXd::Identity(16,16)); //covmatrix.inverse();

  
  for (int chan=1; chan<=16; chan++) {
    for (int i=0; i<16; i++) {
      for (int j=0; j<16; j++) {
	      IY[chan-1][i][j]=covmatrixI(i,j);
      }
    }
  }
    

	eclFgCal* fg = new eclFgCal(); 
	
	if (fg) {
	  printf("Setting IY\n");
	  for (int chan=1; chan<=16; chan++ ) {
	   
      for (int i=0; i<16; i++) {
        for (int j=0; j<16; j++) {
          bool set=fg->setIMElement(chan,i,j,IY[chan-1][i][j]);
          if (!set) printf("error setting IM for chan=%d\n",chan);
        }
      }
	    
	  }
	  printf("Setting IY finished\n");


	  char invmat_dir_name[200];
    strcpy(invmat_dir_name, std::getenv("HOME"));
	  strcat(invmat_dir_name, "/fit/3parFit/invmat");
	  printf("SaveInverseMatrices\n");
	  fg->SaveInverseMatrices(invmat_dir_name);

	  char fname[200];
    strcpy(fname, std::getenv("HOME"));
	  strcat(fname, "/fit/3parFit/panr_11.dat");
	  printf("GetResponsePar for fname=%s\n",fname);
	  fg->GetResponsePar(fname); 
	  printf("CalculateFgPar\n");
	  fg->CalculateFgpar();
	  
	  /////////////////
	  printf("CalculateBit\n");
	  fg->CalculateBit();
	  
	  printf("CalculateFginInt\n");
	  fg->CalculateFginInt();
	  char dir_name[200];
    strcpy(dir_name, std::getenv("HOME"));
	  strcat(dir_name, "/fit/3parFit/DSP_exp10");
	  //sprintf(dir_name,"/home/belle/shtol/belle2/ECL/development/ecl/modules/eclCovmatCollector/DSP_exp10/%d",crate);
	  
	  //	fg->ReadChiCut(dir_name1);
	  if (fg) {
	    printf("WriteCalib\n");
	    fg->WriteCalib(dir_name);
	    printf("WriteCalib finished\n");
	  } else {
	    printf("No FG object\n");
	  }

	  //fg->ReadCalib(dir_name);
	  //	  printf("RootWr for crate %d shaper %d\n",crate, shaper);
	  ///	  fg->RootWr();
	  //if (rootout) fg->RootFill(rootout);
	  printf("deleting FG\n");
	  delete fg;
	  printf("FG deleted\n");
	}
    
    
  
 // if (rootout)  rootout->Close();
  //return c_OK;
  
}

