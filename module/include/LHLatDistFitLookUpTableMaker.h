#ifndef LH_LAT_DIST_FIT_LOOK_UP_TABLE_MAKER_H_INCLUDED
#define LH_LAT_DIST_FIT_LOOK_UP_TABLE_MAKER_H_INCLUDED
#include <string>
#include <hawc-reco/RecoBase.h>
#include <data-structures/event/SimEvent.h>

#include <hawcnest/HAWCNest.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"


/*!
 * @class LHLatDistFitLookUpTableMaker
 * @author Vikas Joshi 
 Module for generating the lookuptables for the LHLatDistFit option of offline-reconstructor. 
 */

/////lookuptable binning definition/////
const int numXmaxBins_ = 12; const float xmaxBinSize_ = 50.0; // in Xmax in g/cm^2
const int numEnergyBins_ = 37; const float energyBinSize_ = 0.1; // in log10(Energy in GeV)
const int numZenithBins_ = 6; const float ZenithBinSize_ = 0.06; // in cos(ZenithAngle)
const float xmaxInitialValue_ = 150.0; const float xmaxEndValue_ = 750.0;
const float energyInitialValue_ = 250.0; const float energyEndValue_ = 1250000.0; // in GeV
const int numRadialDistBins_ = 250; const int numLogPEBins_ = 120;
//////////////////////////////////////////////////////////////////////


class LHLatDistFitLookUpTableMaker : public virtual RecoBase {

  public:

    LHLatDistFitLookUpTableMaker() { }

    const char* GetName() const {return "LHLatDistFitLookUpTableMaker";}

    Configuration DefaultConfiguration();

    void Initialize(const Configuration& config);

    Module::Result Process(BagPtr bag);
    
    void Finish() {
      fileout_->cd();
      tree_->Write();

      for(int i=0; i<numXmaxBins_; i++){ // Xmax loop
	for(int j=0; j<numEnergyBins_; j++){ // energy loop
	  for(int k=0; k<numZenithBins_; k++){ // Zenith angle loop
	    
	    TH2F* histLogrho = (TH2F*)logrhoDist[i][j][k]->Clone("histLogrho");

	    for(int l = 1 ; l<=numRadialDistBins_; l++){ //radial distance loop
	      float totEntProjY = 0.0;
	      totEntProjY = (float)histLogrho->ProjectionY("py_r0",l,l)->Integral();	      
	      for(int m = 1 ; m<=numLogPEBins_; m++){ //npe loop 
		double prob = -12.0; // is the minimum logprob set
		float binContent = 0.0;
		binContent = histLogrho->GetBinContent(l,m);
		if(binContent > 0 && totEntProjY > 0){
		  prob = log(binContent/totEntProjY);
		}
		logprob[i][j][k]->SetBinContent(l,m,prob);
	      }
	    }

	    delete histLogrho;
	 
	    if(hasInputFile_){
	      TH2F* histNewLogrho = (TH2F*)newLogrhoDist[i][j][k]->Clone("histNewLogrho");

	      for(int l = 1 ; l<=numRadialDistBins_; l++){ //radial distance loop
		float totEntProjY = 0.0;

		totEntProjY = (float)histNewLogrho->ProjectionY("py_r0",l,l)->Integral();
	      
		for(int m = 1 ; m<=100; m++){ //npe loop
		  double binContent = 0.0;
		  double prob = 0.0;
		  binContent = histNewLogrho->GetBinContent(l,m);
		  if(binContent > 0 && totEntProjY > 0){
		    prob = binContent/totEntProjY;
		  }
		  newLogprob[i][j][k]->SetBinContent(l,m,prob);
		}
	      }

	      delete histNewLogrho;
	    }
	    
	    char name_hist[200];	     
	    sprintf(name_hist,"Xmax %d to %d g/cm^{2}, Energy %d to %d GeV and Zenith angle %d to %d deg ;r[m];log10(Number of PE)",int (xmaxInitialValue_+(i*xmaxBinSize_)), int(xmaxInitialValue_+(i+1)*xmaxBinSize_),int (pow(10,log10(energyInitialValue_)+j*energyBinSize_)),int (pow(10,(log10(energyInitialValue_)+(j+1)*energyBinSize_))), int (acos(1-k*ZenithBinSize_)/D2R_), int (acos(1-(k+1)*ZenithBinSize_)/D2R_));
	    logprob[i][j][k]->SetTitle(name_hist);
	    logrhoDist[i][j][k]->SetTitle(name_hist);
	    
	    if(hasInputFile_){
	      newLogprob[i][j][k]->SetTitle(name_hist);
	      newLogrhoDist[i][j][k]->SetTitle(name_hist);
	      meanLogAR[i][j][k]->SetTitle(name_hist);
	      sigmaLogAR[i][j][k]->SetTitle(name_hist);
	    }
	    
	    sprintf(name_hist,"logprob_%d_%d_%d",i,j,k);
	    logprob[i][j][k]->Write(name_hist);
	    delete  logprob[i][j][k];
	    
	    sprintf(name_hist,"logrho_dist_%d_%d_%d",i,j,k); 
	    logrhoDist[i][j][k]->Write(name_hist);	     
	    delete logrhoDist[i][j][k];
	    
	    if(hasInputFile_){
	      sprintf(name_hist,"newlogrho_dist_%d_%d_%d",i,j,k);
	      newLogrhoDist[i][j][k]->Write(name_hist);
	      delete newLogrhoDist[i][j][k];

	      sprintf(name_hist,"newlogprob_%d_%d_%d",i,j,k);
	      newLogprob[i][j][k]->Write(name_hist);
	      delete newLogprob[i][j][k];
	    
	      sprintf(name_hist,"MeanLogA_r_%d_%d_%d",i,j,k);
	      meanLogAR[i][j][k]->Write(name_hist);
	      delete meanLogAR[i][j][k];
	 
	      sprintf(name_hist,"SigmaLogA_r_%d_%d_%d",i,j,k);
	      sigmaLogAR[i][j][k]->Write(name_hist);
	      delete sigmaLogAR[i][j][k];
	    }
	  }
	}
      }      
      fileout_->Close();
      RecoBase::Finish();
    }

    TTree* tree_; 
    TFile* filein_;
    TFile* fileout_;

    TH2F* logrhoDist[numXmaxBins_][numEnergyBins_][numZenithBins_];//// array indices : [Xmax][Energy][ZenithAng]
    TH2F* logprob[numXmaxBins_][numEnergyBins_][numZenithBins_];   ///storing the logprob for a given bin combination

    TH2F* newLogrhoDist[numXmaxBins_][numEnergyBins_][numZenithBins_];
    TH2F* newLogprob[numXmaxBins_][numEnergyBins_][numZenithBins_];
    TH1F* meanLogAR[numXmaxBins_][numEnergyBins_][numZenithBins_];
    TH1F* sigmaLogAR[numXmaxBins_][numEnergyBins_][numZenithBins_];
  
    void InitGeo(const det::Detector& det,
		 const ChannelStatusMap& channelStatusMap, float simZang, float simAang, float simXCore, float Ysimcore );

 protected:
    std::string simEventName_; 
    std::string outputRootFile_;
    std::string inputRootFile_;
    bool hasInputFile_;
    bool outriggers_;
    double costh_,sinth_,cosph_,sinph_;
    double D2R_;                           // pi/180 

    /////////////////for data set information for the tree_////////////////////
    double simEnergy_,simXmax_,simZang_,simAang_,simXCore_,simYCore_;
    ///////////////////////////////////////////////////////////////////////////

    class ChannelStatus {
      
    public:
    ChannelStatus() : x_(0), y_(0), pe_(0),
        nGoodPMTs_(0), nUsedPMTs_(0) { }
      double x_, y_, pe_;
      bool ifMA_;
      int nGoodPMTs_, nUsedPMTs_;
      double time_;
    }; 
    
    std::map<int, ChannelStatus> statusMap_;

    /////////////////for data set information for the tree_////////////////////
    std::vector<double> xTank_;
    std::vector<double> yTank_;
    std::vector<double> peTank_;
    std::vector<double> timeTank_;
    std::vector<bool> typeTank_;
    std::vector<double> dist_;
    ///////////////////////////////////////////////////////////////////////////
};


#endif //LH_LAT_DIST_FIT_LOOK_UP_TABLE_MAKER_H_INCLUDED
