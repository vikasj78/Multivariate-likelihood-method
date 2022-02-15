#ifndef LH_LAT_DIST_FIT_H_INCLUDED
#define LH_LAT_DIST_FIT_H_INCLUDED
#include <string>
#include <hawc-reco/RecoBase.h>

#include <hawcnest/HAWCNest.h>
#include <data-structures/detector/Detector.h>
#include <data-structures/event/SimEvent.h>
#include <data-structures/reconstruction/track-fitter/AngleFitResult.h> 
#include <cmath>
#include <map>

#include "TMinuit.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "TMinuitMinimizer.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TEllipse.h"
#include "TMarker.h"
#include "TLine.h"
#include "TArrow.h"


/*!
 * @class LHLatDistFit
 * @author Vikas Joshi
 * @ingroup core_fitters
 * @By fitting shower lateral charge distribution fits core, energy and Xmax.
 
 * This is a template based Likehood fit method, which uses templates generated 
 * by using MC simulations. 
 * As a fit result it gives core, Energy ,Xmax and Goodness of Fit of the shower.
 */

/* ///////initialize the binnings and ranges for the parameters//////////////// */
const double arrayCenterX_ = 35.0; const double arrayCenterY_ = 245.0;

const int numXmaxBins_ = 12; const double xmaxBinSize_ = 50.0; // in Xmax in g/cm^2
const int numEnergyBins_ = 33; const double energyBinSize_ = 0.1; // in log10(Energy in GeV)
const int numZenithBins_ = 6; const double zenithBinSize_ = 0.06; // in cos(ZenithAngle)
const int numRDistBins_ = 250; const double rDistBinSize_ = 2.0; //in meters
const int numNPEBins_ = 120; const double nPEBinSize_ = 0.1; //in log(nPE_)

const double xmaxInitialValue_ = 150.0; const double xmaxEndValue_ = 750.0;
const double energyInitialValue_ = 2.39; const double energyEndValue_ = 5.69; // log10(250GeV) and log10(500000GeV) */
const double rDistInitialValue_ = 0.0; const double rDistEndValue_ = 500.0; //in meters
const double nPEInitialValue_ = -4.0; const double nPEEndValue_ = 8.0; //in log(nPE_)

const double xmaxVsEnergyNorm_ = 1.17923e+02; const double xmaxVsEnergySlope_ = 7.89186e+01; //comes from the Xmax Energy relation from MC
const double maxRadialDist_ = 495.0; //m lookups are filled till 500 m
//////////////////////////////////////////////////////////////////////

class LHLatDistFitResult;

class LHLatDistFit : public RecoBase {

  public:
  
  LHLatDistFit() { };
  
  const char* GetName() const {return "LHLatDistFit";}
  
  Configuration DefaultConfiguration();
  
  void Initialize(const Configuration& config);
  
  Module::Result Process(BagPtr bag);
    
  void Finish() {
    for(int i=0; i<numXmaxBins_; i++){  //Xmax loop
      for(int j=0; j<numEnergyBins_; j++){ // energy loop
	for(int k=0; k<numZenithBins_; k++){ // Zenith angle loop
	  delete histLogprobMainArray[i][j][k];
	  delete histLogprobOutriggers[i][j][k];
	}
      }
    }
    delete histLogprobXmaxEnergy;
    delete histNhitREnergyGuess;
    delete histLikelihoodVsNhitTanks;

    if(drawLHSurface_ && display_){
      delete lhSurface;  delete lhSurfaceZoom; 
      delete lhEnergy; delete lhXmax;
      delete grMATanks; delete grORTanks;
      delete LDFMA; delete LDFOR;
      
      delete mSeed;  delete mReco;
      delete leg;
      delete ellipse; delete line;
      
      delete padSurfaceXY; delete padSurfaceXYZoom;
      delete padLeg; delete padSurfaceEnergy;
      delete padSurfaceXmax;
      delete canvas;
      delete pad; delete zeroPad; delete padOR; delete zeroPadOR; 
      delete canvasPDF;
    }
    RecoBase::Finish();
  }
  
  TFile* probFileMainArray;
  TFile* probFileOutriggers;
  TFile* probFileXmaxEnergy;
  TFile* guessFileNhitREnergy;
  TFile* lookupFileLikelihoodVsNhitTanks;

  TH2F* histLogprobMainArray[numXmaxBins_][numEnergyBins_][numZenithBins_];
  TH2F* histLogprobOutriggers[numXmaxBins_][numEnergyBins_][numZenithBins_];
  TH2F* histLogprobXmaxEnergy;
  TH2F* histNhitREnergyGuess;
  TH1F* histLikelihoodVsNhitTanks;
  
  ROOT::Math::Minimizer* min_;
  ROOT::Math::Functor* functor_;
  double LikelihoodFunction(const double* P);
  static void LikelihoodToMinuit(int &npar, double* gin, double& L, double* P, int iflag);
  double Prob(float radialDist,  double rho, double xmax, double energy, int* bins, double* vals, bool tankType);
  double ProbHist(float radialDist, double rho, int histXmaxNum, int histEnergyNum, bool tankType);
  double InterpolateInsideHist( TH2F* hist, double x, double y);
  double GetProbFromHist(TH2F* hist, double x, double y);
  double InterpolateXmaxEnergy( TH2F* hist, int xBin, int yBin, double x, double y);
  double InterpolateTwoPoints(double x0, double y0, double x1, double y1, double x);
  double GetEnergyGuess(double coreX, double coreY, int nHitMA);
  double CalculateGoF(double minLikelihood, int nHitTanksMA, int nZeroTanksMA, int nGoodTanksMA);
  void GetEnergyXmaxHistsandValues(double xmax, double energy, int* bins, double* vals);
  void DrawLikelihoodSurfaceLDFandPDF(double xGuess, double yGuess);
  void InitGeo(const det::Detector& det,
	       const ChannelStatusMap& channelStatusMap);
  double FitModule(double* P, bool* PFix);
  void FitProcedure(double xc, double yc);

  ////DrawLikelihoodSurfaceLDFandPDF objects
  TCanvas *canvas, *canvasPDF;
  TPad *padSurfaceXY, *padSurfaceXYZoom, *padLeg, *padSurfaceEnergy, *padSurfaceXmax;
  TPad *pad, *zeroPad, *padOR, *zeroPadOR;
  TH2F *lhSurface, *lhSurfaceZoom, *logprobMA, *logprobOR;
  TH1F *lhEnergy, *lhXmax, *zeroLike, *zeroLikeOR;
  TGraph *grMATanks, *grORTanks, *LDFMA, *LDFOR;
  TMarker *mSeed, *mReco;
  TArrow *radialLine;
  TLegend *leg;
  TEllipse *ellipse;
  TLine *line;

 protected:
  std::string simEventName_;
  std::string angleFitName_;
  std::string probFileMainArray_;
  std::string probFileOutriggers_;
  std::string probFileXmaxEnergy_;
  std::string guessFileNhitREnergy_;
  std::string lookupFileLikelihoodVsNhitTanks_;

  double costh_,sinth_,cosph_,sinph_;
  double minLogPE_;
  int minLogPEBin_;
  double gBestXmax_,gBestEnergy_,gBestXc_,gBestYc_,gMinLikelihood_,gGoF_;
  double gMinL_, gZeroMinL_, gHitMinL_;
  unsigned int nHitTanksMA_,nHitTanksOR_,nZeroTanksMA_,nZeroTanksOR_ ;
  int zenithBin_;
  
  //used for taking into account only nth zero while nHitTanksMA_ is smaller than 
  //to speed up the fitting in real data 
  unsigned int minNHitTanksMA_, everyNthZero_;

  int minimumTanks_;
  double tolerance_;
  int maxIterations_;
  int strategy_;
  bool fitXmax_;
  bool drawLHSurface_;
  bool display_;
 
  class ChannelStatus {
    
  public:
  ChannelStatus() : x_(0), y_(0), pe_(0),
      nGoodPMTs_(0), nUsedPMTs_(0), probWeight_(1) { }
    double x_, y_, pe_;
    bool ifMA_;
    int nGoodPMTs_, nUsedPMTs_, probWeight_;
  }; 
  
  int statusRunID_;
  int statusTimeSliceID_;
  std::map<int, ChannelStatus> statusMap_;
  std::vector<ChannelStatus> data_;
};


#endif // LH_LAT_DIST_FIT_H_INCLUDED

