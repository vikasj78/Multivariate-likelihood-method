#include <core-fitter/LHLatDistFit.h>
#include <core-fitter/CoreChannelStatus.h>

#include <hawc-reco/hit-selection/HitSelector.h>
#include <hawc-reco/hit-selection/HitSelectorOR.h>

#include <data-structures/reconstruction/track-fitter/AngleFitResult.h> 
#include <data-structures/reconstruction/core-fitter/LHLatDistFitResult.h>
#include <data-structures/event/Event.h>
#include <data-structures/event/ChannelEvent.h>
#include <data-structures/detector/Detector.h>

#include <detector-service/DetectorService.h>
#include <detector-service/ChannelStatusMap.h>

#include <hawcnest/HAWCNest.h>
#include <hawcnest/HAWCUnits.h>
#include <hawcnest/ConfigurationUtil.h>
#include <hawcnest/RegisterService.h>

#include <iterator>
#include <algorithm>
#include <map>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"


using namespace std;
using namespace HAWCUnits;
REGISTER_SERVICE(LHLatDistFit);

Configuration LHLatDistFit::DefaultConfiguration() {
  
  Configuration config = RecoBase::DefaultConfiguration();
  config.Parameter<std::string>("angleFitName");
  config.Parameter<bool>("fitXmax"); 
  config.Parameter<bool>("drawLHSurface");
  config.Parameter<bool>("display");
  config.Parameter<std::string>("probFileMainArray");
  config.Parameter<std::string>("probFileOutriggers");
  config.Parameter<std::string>("probFileXmaxEnergy");
  config.Parameter<std::string>("guessFileNhitREnergy");
  config.Parameter<std::string>("lookupFileLikelihoodVsNhitTanks");

  /// Maximum number of minimizer iterations
  config.Parameter<int>("maxIterations", 400);
  
  /// Acceptable tolerance to stop fit
  config.Parameter<double>("tolerance", 0.02);
  
  /// Minimum number of hits required to perform the fit
  config.Parameter<int>("minimumTanks", 2);
  
  /// Search strategy
  config.Parameter<int>("strategy", 0); //0 =fast, 1=default, 2=slow but more accurate
  
  return config;
}

void LHLatDistFit::Initialize(const Configuration& config) {
  log_info("Initializing LHLatDistFit...");
  RecoBase::Initialize(config);
  config.GetParameter("angleFitName",angleFitName_);
  config.GetParameter("fitXmax",fitXmax_);
  config.GetParameter("drawLHSurface",drawLHSurface_);
  config.GetParameter("display",display_);
  config.GetParameter("probFileMainArray",probFileMainArray_);
  config.GetParameter("probFileOutriggers",probFileOutriggers_);
  config.GetParameter("probFileXmaxEnergy",probFileXmaxEnergy_);
  config.GetParameter("guessFileNhitREnergy",guessFileNhitREnergy_);
  config.GetParameter("lookupFileLikelihoodVsNhitTanks",lookupFileLikelihoodVsNhitTanks_);
  config.GetParameter("maxIterations", maxIterations_);
  config.GetParameter("tolerance", tolerance_);
  config.GetParameter("minimumTanks", minimumTanks_);
  config.GetParameter("strategy", strategy_);
  
  statusRunID_ = -1;
  statusTimeSliceID_ = -1;
  costh_ = 1.0, sinth_ = 0.0, cosph_ = 1.0, sinph_ = 0.0;
  minLogPE_ = -3.0;
  zenithBin_ = 50;
  minNHitTanksMA_ = 100;
  everyNthZero_ = 5;
  TDirectory *cur_dir = gDirectory;
  
  probFileMainArray = new TFile(probFileMainArray_.c_str());
  cur_dir->cd();
  if (probFileMainArray->IsZombie() || !probFileMainArray->IsOpen()) {
    log_fatal("The file [" << probFileMainArray_ << "] does not exist ==> BIG PROBLEM !");
  }
  
  probFileOutriggers = new TFile(probFileOutriggers_.c_str());
  cur_dir->cd();
  if (probFileOutriggers->IsZombie() || !probFileOutriggers->IsOpen()) {
    log_fatal("The file [" << probFileOutriggers_ << "] does not exist ==> BIG PROBLEM !");
  }
  
  probFileXmaxEnergy = new TFile(probFileXmaxEnergy_.c_str());
  cur_dir->cd();
  if (probFileMainArray->IsZombie() || !probFileMainArray->IsOpen()) {
    log_fatal("The file [" << probFileXmaxEnergy_ << "] does not exist ==> BIG PROBLEM !");
  }
  
  histLogprobXmaxEnergy = (TH2F*)probFileXmaxEnergy->Get("h_prob");
  if (!histLogprobXmaxEnergy) {
    log_fatal("Can't find the histogram \"h_prob\" in file [" << probFileXmaxEnergy_ << "]");
  }
  histLogprobXmaxEnergy->SetDirectory(0);
  
  guessFileNhitREnergy = new TFile(guessFileNhitREnergy_.c_str());
  cur_dir->cd();
  if (guessFileNhitREnergy->IsZombie() || !guessFileNhitREnergy->IsOpen()) {
    log_fatal("The file [" << guessFileNhitREnergy_ << "] does not exist ==> BIG PROBLEM !");
  }
  
  histNhitREnergyGuess = (TH2F*)guessFileNhitREnergy->Get("nHit_r_energy_guess");
  if (!histNhitREnergyGuess) {
    log_fatal("Can't find the histogram \"nHit_r_energy_guess\" in file [" << guessFileNhitREnergy_ << "]");
  }
  histNhitREnergyGuess->SetDirectory(0);
  
  lookupFileLikelihoodVsNhitTanks = new TFile(lookupFileLikelihoodVsNhitTanks_.c_str());
  cur_dir->cd();
  if (lookupFileLikelihoodVsNhitTanks->IsZombie() || !lookupFileLikelihoodVsNhitTanks->IsOpen()) {
    log_fatal("The file [" << lookupFileLikelihoodVsNhitTanks_ << "] does not exist ==> BIG PROBLEM !");
  }
  
  histLikelihoodVsNhitTanks = (TH1F*)lookupFileLikelihoodVsNhitTanks->Get("likelihood_vs_nHitTanks");
  if (!histLikelihoodVsNhitTanks) {
    log_fatal("Can't find the histogram \"likelihood_vs_nHitTanks\" in file [" << lookupFileLikelihoodVsNhitTanks_ << "]");
  }
  histLikelihoodVsNhitTanks->SetDirectory(0);
  
  bool isMinLogPEBinSet_ = false;
  for(int i=0; i<numXmaxBins_; i++){  //Xmax loop
    for(int j=0; j<numEnergyBins_; j++){ // energy loop
      for(int k=0; k<numZenithBins_; k++){ // Zenith angle loop
        
        std::ostringstream oss_histname;
        oss_histname << "logprob_" << i << "_" << j << "_" << k;
        std::string histname = oss_histname.str();
        
        //getting logprob_hist for main array
        TH2F *tmpMainArrayHist = (TH2F*)probFileMainArray->Get(histname.c_str());
        if (!tmpMainArrayHist) {
          log_fatal("Can't find the histogram [" << histname << "] in file [" << probFileMainArray_ << "]");
        }
        histLogprobMainArray[i][j][k] =  tmpMainArrayHist;
        histLogprobMainArray[i][j][k]->SetDirectory(0);
        if(!isMinLogPEBinSet_) minLogPEBin_ = (int)histLogprobMainArray[i][j][k]->GetYaxis()->FindBin(minLogPE_);
        //getting logprob_hist for outriggers array
        TH2F *tmpOutriggersHist = (TH2F*)probFileOutriggers->Get(histname.c_str());
        if (!tmpOutriggersHist) {
          log_fatal("Can't find the histogram [" << histname << "] in file [" << probFileOutriggers_ << "]");
        }
        histLogprobOutriggers[i][j][k] =  tmpOutriggersHist;
        histLogprobOutriggers[i][j][k]->SetDirectory(0);
      }
    }
  }
  
  delete probFileMainArray;
  delete probFileOutriggers;
  delete probFileXmaxEnergy;
  delete guessFileNhitREnergy;
  delete lookupFileLikelihoodVsNhitTanks;
  
  min_ = ROOT::Math::Factory::CreateMinimizer("Minuit", "Simplex");
  functor_ = new ROOT::Math::Functor(this,&LHLatDistFit::LikelihoodFunction,4);
 
  //initiallize only if drawLHSurface
  if(drawLHSurface_ && display_){
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10);
    gStyle->SetOptStat("");

    double margin = 0.15;

    ///canvas and pads for likelihood surface
    canvas = new TCanvas("canvas", "", 900, 1400);
    padSurfaceXY = new TPad("padSurfaceXY", "", 0.0, 0.4, 0.5, 0.9);
    padSurfaceXYZoom = new TPad("padSurfaceXYZoom", "", 0.5, 0.4, 1.0, 0.9);
    padLeg = new TPad("padLeg","",0.0, 0.9, 1.0, 1.0);
    padSurfaceEnergy = new TPad("padSurfaceEnergy","",0.0,0.0,0.5,0.4);
    padSurfaceXmax = new TPad("padSurfaceXmax", "",0.5,0.0,1.0,0.4);
    canvas->Draw();
    padSurfaceXY->SetTopMargin(0.01);
    padSurfaceXY->SetTicks(1,1);
    padSurfaceXY->Draw();
    padSurfaceXYZoom->SetTopMargin(0.01);
    padSurfaceXYZoom->SetTicks(1,1);
    padSurfaceXYZoom->Draw();
    padLeg->SetBottomMargin(0.01);
    padLeg->Draw();
    padSurfaceEnergy->SetTopMargin(0.01);
    padSurfaceEnergy->Draw();
    padSurfaceXmax->SetTopMargin(0.01);
    padSurfaceXmax->Draw();

    ////canvas and pads for LDF and PDF plot
    canvasPDF = new TCanvas("canvasPDF","LDF and PDF",1200,700);
    pad = new TPad("pad", "", 0.0, 0.25, 0.5, 1.0);
    zeroPad = new TPad("zeroPad", "", 0.0, 0.0, 0.5, 0.368);
    padOR = new TPad("padOR", "", 0.5, 0.25, 1.0, 1.0);
    zeroPadOR = new TPad("zeroPadOR", "", 0.5, 0.0, 1.0, 0.368);
    canvasPDF->Draw();
    pad->SetRightMargin(margin);
    pad->Draw();
    zeroPad->SetTicks(1,1);
    zeroPad->SetRightMargin(margin);
    zeroPad->SetTopMargin(0);
    zeroPad->SetBottomMargin(0.33);
    zeroPad->Draw();
    padOR->SetRightMargin(margin);
    padOR->Draw();
    zeroPadOR->SetTicks(1,1);
    zeroPadOR->SetRightMargin(margin);
    zeroPadOR->SetTopMargin(0);
    zeroPadOR->SetBottomMargin(0.33);
    zeroPadOR->Draw();
    
    //initialize all draw pointers to null
    lhSurface = lhSurfaceZoom = logprobMA = logprobOR = NULL;
    lhEnergy = lhXmax = zeroLike = zeroLikeOR = NULL;
    grMATanks = grORTanks = LDFMA = LDFOR = NULL;
    mSeed = mReco = NULL;
    leg = NULL;
    ellipse = NULL;
    line = NULL;
  }
  log_info("...Finished initializing LHLatDistFit");
}

void LHLatDistFit::InitGeo(const det::Detector& det,
                           const ChannelStatusMap& channelStatusMap) {
  
  // Cycle through Channels
  statusMap_.clear();
  ChannelStatus status;
  for (det::Detector::ConstTankIterator it = det.TanksBegin();
       it != det.TanksEnd(); ++it) {

    status.nGoodPMTs_ = 0;
    // Ignore channels marked as bad
    for (det::Detector::ConstChannelIterator
         channel = it->ChannelsBegin();
         channel != it->ChannelsEnd(); ++channel) {
 
      if (channelStatusMap.IsGood(channel->GetChannelId())) {
        ++status.nGoodPMTs_;	
      }
    }
    
    status.pe_      = 0.;
    status.x_       = it->GetPosition().GetX();
    status.y_       = it->GetPosition().GetY();

    if(it->GetTankId() < 301){ //main array tanks = 300
      status.ifMA_      = true; //flag main array tanks
    }
    else {
      status.ifMA_ = false; //flag outriggers
      status.nUsedPMTs_ = 1;
      status.nGoodPMTs_ = 1;//for now take all outrigger PMTs as used if outriggers are used because there is no good channel list available for ORs.
    }
    statusMap_[it->GetTankId()] = status;
  }
}

double LHLatDistFit::InterpolateTwoPoints(double x0, double y0, double x1, double y1, double x){
  // If sampling points are too close, assume they are the same point
  // and return
  if(fabs(x1 - x0) < 1e-14) return y0;
  else {
    double result = y0 + ((x - x0) * (y1 - y0) / (x1 - x0));
    return result;
  }
}


double LHLatDistFit::GetProbFromHist(TH2F* hist, double x, double y){
 
  int xBin, yBin;
  double prob = 0.;
  xBin = int((x - rDistInitialValue_)/rDistBinSize_) + 1;
  
  if(y > minLogPE_) yBin = int((y - nPEInitialValue_)/nPEBinSize_) + 1;
  else yBin = minLogPEBin_; // zeros stored at the log scale at minLogPE_

  prob = (double)hist->GetBinContent(xBin , yBin);
  return prob;
}


double LHLatDistFit::ProbHist(float radialDist, double logrho, int histXmaxNum, int histEnergyNum, bool tankType){
  
  TH2F* logprobHistTemp;
  
  if(tankType) logprobHistTemp = histLogprobMainArray[histXmaxNum][histEnergyNum][zenithBin_]; //for MA
  else logprobHistTemp = histLogprobOutriggers[histXmaxNum][histEnergyNum][zenithBin_]; // for OR
  
  double probHist =  GetProbFromHist(logprobHistTemp,radialDist,logrho);
  //double probHist = InterpolateInsideHist(logprobHistTemp,radialDist,logrho);
  return probHist;
}

double LHLatDistFit::Prob(float radialDist,  double logrho, double xmax, double energy, int* bins, double* vals, bool tankType){
  
  /////here bins and vals have been defined as follows
  ///bins = {histXmax0, histEnergy0, histXmax1, histEnergy1}
  ///vals = {Xmax0, Energy0, Xmax1, Energy1}
  
  double probHist00 = 0.0, probHist01 = 0.0, probHist10 = 0.0, probHist11 = 0.0;
  double prob = 0.0;
  
  probHist00 = ProbHist(radialDist, logrho, bins[0], bins[1], tankType);
  
  if(bins[0] != bins[2] && bins[1] != bins[3]){//two Xmax and two energy, do interpolation in both
    double probEnergy0 = 0.0, probEnergy1 = 0.0;
    probHist01 = ProbHist(radialDist, logrho, bins[2], bins[1], tankType);
    probHist10 = ProbHist(radialDist, logrho,  bins[0], bins[3], tankType);
    probHist11 = ProbHist(radialDist, logrho, bins[2], bins[3], tankType);
    probEnergy0 = InterpolateTwoPoints(vals[0] , probHist00, vals[2], probHist01, xmax); //interpolate in Xmax
    probEnergy1 = InterpolateTwoPoints(vals[0], probHist10, vals[2], probHist11, xmax); //interpolate in Xmax
    prob = InterpolateTwoPoints(vals[1], probEnergy0, vals[3], probEnergy1, energy); //interpolate in Energy
  }
  
  else if(bins[0] != bins[2] && bins[1] == bins[3]){ //same energy but two different Xmax, therefore only interplote in Xmax
    probHist01 = ProbHist(radialDist, logrho, bins[2], bins[1], tankType);
    prob = InterpolateTwoPoints(vals[0], probHist00, vals[2], probHist01, xmax);
  }
  
  else if(bins[0] == bins[2] && bins[1] != bins[3]){ //same xmax but two different energy, therefore only interpolate in energy
    probHist10 = ProbHist(radialDist, logrho,  bins[0], bins[3], tankType);
    prob = InterpolateTwoPoints(vals[1], probHist00, vals[3], probHist10, energy);
  }
    
  else{ //nothing to interpolate only one value of energy and Xmax 
    prob = probHist00;
  }
  return prob;
}

void LHLatDistFit::GetEnergyXmaxHistsandValues(double xmax, double energy, int* bins, double* vals){
  
  float xmax0 = 0.0, xmax1 = 0.0;
  float energy0 = 0.0, energy1 = 0.0;
  int  histXmaxNum = 50, histEnergyNum = 50;
  
  histXmaxNum = (xmax - xmaxInitialValue_)/xmaxBinSize_;
  histEnergyNum = (energy - energyInitialValue_)/energyBinSize_;
  
  xmax0 = xmaxInitialValue_+((2*histXmaxNum+1)*xmaxBinSize_)/2.0;
  energy0 = energyInitialValue_+((2*histEnergyNum+1)*energyBinSize_)/2.0;
  
  int histXmaxNumTemp = histXmaxNum;
  int histEnergyNumTemp = histEnergyNum;
  
  if (energy0 > energy && histEnergyNum > 0) histEnergyNum--;
  else if(energy0 < energy && histEnergyNum < numEnergyBins_ - 1) histEnergyNum++;
  
  if (xmax0 > xmax && histXmaxNum > 0) histXmaxNum--;
  else if(xmax0 < xmax && histXmaxNum < numXmaxBins_ - 1) histXmaxNum++;
 
  energy1 = energyInitialValue_+((2*histEnergyNum + 1)*energyBinSize_)/2.0;
  xmax1 = xmaxInitialValue_+((2*histXmaxNum + 1)*xmaxBinSize_)/2.0;
  
  bins[0] = histXmaxNumTemp;
  bins[1] = histEnergyNumTemp;
  bins[2] = histXmaxNum;
  bins[3] = histEnergyNum;
  vals[0] = xmax0;
  vals[1] = energy0;
  vals[2] = xmax1;
  vals[3] = energy1;
}

double LHLatDistFit::LikelihoodFunction(const double* P){
  
  bool verbose = true;
  double L = 0.0;
  double LHit = 0.0;
  double LZero = 0.0;
  float radialDist;
  
  double xmax = P[0];
  double energy = P[1];
   
  //if xmax is fixed calculate the value from energy
  if (!fitXmax_) {
    xmax = xmaxVsEnergyNorm_  + xmaxVsEnergySlope_ * energy;
    if(xmax < xmaxInitialValue_) xmax = xmaxInitialValue_;
    if(xmax > xmaxEndValue_) xmax = xmaxEndValue_;
  }

  //get the four bin numbers and values corresponding to xmax0, xmax1, energy0, energy1 where the linear interpolation will be done
  int bins[4] = {50,50,50,50};
  double vals[4] = {0.0,0.0,0.0,0.0};
  GetEnergyXmaxHistsandValues(xmax,energy,bins,vals);
  
  //if we are at the edge of the phase space in distance from the shower core the PDFs become flat therefore no need to do the interpolation between energy and xmax
  //and to get the value of probability at far distances for hits we take the value at maxRadialDist_ m and for zero hits there is no need to reapeat that for all zero therefore we do it for only one case and store it.
  TH2F* histMA = (TH2F*)histLogprobMainArray[bins[0]][bins[1]][zenithBin_];
  TH2F* histOR = (TH2F*)histLogprobOutriggers[bins[0]][bins[1]][zenithBin_];

  if (!histMA) log_fatal("Can't find the histogram histMA");
  if (!histOR) log_fatal("Can't find the histogram histOR");
  
  bool zeroPEDistFlagMA = false;
  bool zeroPEDistFlagOR = false;
  
  double prob = -12.0;
  double probZeroPEEdgeMA = -12.0;
  double probZeroPEEdgeOR = -12.0;
  
  //loop over all the selected tanks 
  for (std::vector<ChannelStatus>::iterator d = data_.begin();  d < data_.end(); ++d) {
    radialDist = sqrt((d->x_-P[2])*(d->x_-P[2])+(d->y_-P[3])*(d->y_-P[3])-sinth_*sinth_*((d->x_-P[2])*cosph_+(d->y_-P[3])*sinph_)*((d->x_-P[2])*cosph_+(d->y_-P[3])*sinph_));
    
    //check for the large impact distances for MA
    if(d->ifMA_ && radialDist > maxRadialDist_){
      if(d->pe_ > minLogPE_){
        int xBin = (int)histMA->GetXaxis()->FindBin(maxRadialDist_);
        int yBin = (int)histMA->GetYaxis()->FindBin(d->pe_);
        prob = histMA->GetBinContent(xBin, yBin);
      }
      
      else{
        if (!zeroPEDistFlagMA){
          zeroPEDistFlagMA = true;
          int xBin = (int)histMA->GetXaxis()->FindBin(maxRadialDist_);
          int yBin = minLogPEBin_;
          probZeroPEEdgeMA = histMA->GetBinContent(xBin, yBin);
	  prob = probZeroPEEdgeMA;
        }
        else prob = probZeroPEEdgeMA;	
	//for the events with low nHitTanksMA_ we only take every nth zero into account to speed up the fitting 
	//there are lot of those in real data and this is where all the computation time goes
	if(nHitTanksMA_ < minNHitTanksMA_) prob = d->probWeight_ * prob;	       
      }
    }

    //check for large impact distances for OR
    else if(!d->ifMA_ && radialDist > maxRadialDist_){
      if(d->pe_ > minLogPE_){
        int xBin = (int)histOR->GetXaxis()->FindBin(maxRadialDist_);
        int yBin = (int)histOR->GetYaxis()->FindBin(d->pe_);
        prob = histOR->GetBinContent(xBin, yBin);
      }
      
      else{
        if (!zeroPEDistFlagOR){
          zeroPEDistFlagOR = true;
          int xBin = (int)histOR->GetXaxis()->FindBin(maxRadialDist_);
          int yBin = minLogPEBin_;
          probZeroPEEdgeOR = histOR->GetBinContent(xBin, yBin);
          prob = probZeroPEEdgeOR;
        }
        else prob = probZeroPEEdgeOR;
	//for the events with low nHitTanksMA_ we only take every nth zero into account to speed up the fitting                                              
        //there are lot of those in real data and this is where all the computation time goes 
	//this we do same for OR because number of hit MA tanks decide the size of the event for now maybe later will change it for OR specifically
        if(nHitTanksMA_ < minNHitTanksMA_) prob = d->probWeight_ * prob;
      } 
    }

    // otherwise go for normal interpolation in all dimensions
    else {
      prob = Prob(radialDist,d->pe_,xmax,energy,bins,vals,d->ifMA_);
      //for the events with low nHitTanksMA_ we only take every nth zero into account to speed up the fitting
      //there are lot of those in real data and this is where all the computation time goes
      if(nHitTanksMA_ <  minNHitTanksMA_ && !(d->pe_ > minLogPE_))  prob = d->probWeight_ * prob; //multiply the probability of nth zero by n      
    }
    
    L -= prob;
    
    if(d->pe_ > minLogPE_) LHit -= prob;
    else LZero -= prob;

    if (L != L && verbose) {
      log_fatal("radialDist: " << radialDist << " prob: " << prob << " Likelihood: " << L);
      verbose = false;
    }
  }
  if (L != L)  log_fatal("Likelihood: " << L << " Xmax: " << P[0] << " Energy: "<< P[1] << " Xcore: "<< P[2] << " Ycore: "<< P[3]);
  
  log_debug("Likelihood: "<< L << " Xmax: "<< P[0] <<" Energy: "<< P[1] << " Xcore:  "<< P[2] <<" Ycore: "<< P[3]);
  
  if(L < gMinL_){
    gMinL_ = L;
    gZeroMinL_ = LZero;
    gHitMinL_ = LHit;
  }

  return L;
}

double LHLatDistFit::GetEnergyGuess(double coreX, double coreY, int nHitMA){
  
  double distFromCenter = sqrt((coreX-arrayCenterX_)*(coreX-arrayCenterX_)+(coreY-arrayCenterY_)*(coreY-arrayCenterY_));
  
  int distBin, tankHitBin;
  distBin = tankHitBin = 1;
  
  if(distFromCenter > 980.) distFromCenter = 980.;
  if(nHitMA > 299) nHitMA = 299;
  
  distBin = histNhitREnergyGuess->GetXaxis()->FindBin(distFromCenter);
  tankHitBin = histNhitREnergyGuess->GetYaxis()->FindBin(nHitMA);
  
  double eGuess = histNhitREnergyGuess->GetBinContent(distBin,tankHitBin);
  
  return eGuess;
}

double LHLatDistFit::FitModule(double* P, bool* PFix){
  min_->Clear();
  min_->SetFunction(*functor_);
  min_->SetMaxFunctionCalls(maxIterations_);
  min_->SetStrategy(strategy_); // 0 = fast
  min_->SetErrorDef(0.5); //for negative log likelihood
  min_->SetPrintLevel(0);

  if(nHitTanksMA_ < minNHitTanksMA_) min_->SetTolerance(10*tolerance_);
  else min_->SetTolerance(tolerance_);

  double distFromCore = sqrt((P[2]-arrayCenterX_)*(P[2]-arrayCenterX_) + (P[3]-arrayCenterY_)*(P[3]-arrayCenterY_));
  double distStep = 3.0 + distFromCore*0.05; /// cos(45)/4/10 = 0.176/10
  double energyStep = 0.3;
  double stepSize[4] = {80.0,energyStep,distStep,distStep};
  
  // Set the free variables to be minimized!
  min_->SetVariable(0,"Xmax", P[0], stepSize[0]);
  min_->SetVariable(1,"Energy", P[1], stepSize[1]);
  min_->SetVariable(2,"xFitCore", P[2], stepSize[2]);
  min_->SetVariable(3,"yFitCore", P[3], stepSize[3]);
  // Set the variable limits! tables are filled until certain value
  min_->SetVariableLimits(0, xmaxInitialValue_, xmaxEndValue_);
  min_->SetVariableLimits(1, energyInitialValue_, energyEndValue_);
  min_->SetVariableLimits(2,  arrayCenterX_ - maxRadialDist_, arrayCenterX_  + maxRadialDist_);
  min_->SetVariableLimits(3, arrayCenterY_ - maxRadialDist_, arrayCenterY_ + maxRadialDist_);  
  //setting bounds on X and Y from the array center,  using lookuptable distance till maxRadialDist_ m   // after this part of the phase sapce we do not have any fitting power so no reason to wander around that area..                                                                                               
  ////fix variable if needed
  if(PFix[0]) min_->FixVariable(0);
  if(PFix[1]) min_->FixVariable(1);
  if(PFix[2]) min_->FixVariable(2);
  if(PFix[3]) min_->FixVariable(3);
  min_->Minimize();

  const double *POut = min_->X();
  
  //returning the fit value of the variables and the minimum likelihood value
  P[0] = POut[0];
  if(PFix[0]) P[0] = xmaxVsEnergyNorm_ + xmaxVsEnergySlope_*POut[1]; //this is because we used this value in the likelihood calculation although it is fixed in the minimizer
  P[1] = POut[1];
  P[2] = POut[2];
  P[3] = POut[3];

  double LMin = min_->MinValue();
  return LMin;
}

void LHLatDistFit::FitProcedure(double xc, double yc){
  
  double PBest[4];
  bool PFix[4] = {false, false, false, false};

  double energy = GetEnergyGuess(xc, yc, nHitTanksMA_); //getting it from GetEnergyGuess function.                                                  
  PBest[0] = xmaxVsEnergyNorm_ + xmaxVsEnergySlope_*energy;
  PBest[1] = energy;
  PBest[2] = xc;
  PBest[3] = yc;
  PFix[0] = !fitXmax_; //fix xmax if it is calculated using energy                                                                                                                               
  gMinLikelihood_ = FitModule(PBest, PFix);
  gBestXmax_ = PBest[0];
  gBestEnergy_ = PBest[1];
  gBestXc_ = PBest[2];
  gBestYc_ = PBest[3];
}

double LHLatDistFit::CalculateGoF(double minLikelihood, int nHitTanksMA, int nZeroTanksMA, int nGoodTanksMA){
  int tanks = nHitTanksMA + (nGoodTanksMA - (nZeroTanksMA + nHitTanksMA));
  double meanL = histLikelihoodVsNhitTanks->GetBinContent(tanks);
  double rmsL = histLikelihoodVsNhitTanks->GetBinError(tanks);
  double gof = -1.e6;

  if(rmsL > 0)  gof =  (minLikelihood - meanL)/rmsL;
  else gof = -1.e6;
  
  return gof;
}

Module::Result LHLatDistFit::Process(BagPtr bag) {
  
  const evt::Event& event  = bag->Get<evt::Event>(event_);
  
  //Load seed
  if (!bag->Exists(seed_)) {
    log_fatal("No seed value found, need to call the COM/SFCF service before, requested seed:" << seed_ );
  }

  const CoreFitResult& seed = bag->Get<CoreFitResult>(seed_);
  const AngleFitResult& angleFitResult = bag->Get<AngleFitResult>(angleFitName_);
  
  float ZangRec = angleFitResult.GetTheta();
  float AangRec = angleFitResult.GetPhi();
  costh_ = cos(ZangRec);
  sinth_ = sin(ZangRec);
  cosph_ = cos(AangRec);
  sinph_ = sin(AangRec);
  
  if(ZangRec < 0.872665){
    zenithBin_ = (1-cos(ZangRec))/zenithBinSize_;// 50 deg to radian = 0.872665
  }
  else{
    zenithBin_ = numZenithBins_ - 1;
  }
  // Load geometry
  DetectorService& detSv   = GetService<DetectorService>(detectorService_);
  const det::Detector& det = detSv.GetDetector(event.GetTime());
  
  // Load bad channel map
  const ChannelStatusMap& channelStatusMap = detSv.GetChannelStatusMap(event.GetTime());
  if (!(event.GetRunID() == statusRunID_ &&
        event.GetTimeSliceID() == statusTimeSliceID_)) {
    InitGeo(det, channelStatusMap);
    statusRunID_ = event.GetRunID();
    statusTimeSliceID_ = event.GetTimeSliceID();
  }
  
  // Set all tanks that have good channels as good
  int nGood = 0;  //stores good tanks before seeing the event.
  for (std::map<int, ChannelStatus>::iterator it = statusMap_.begin();
       it != statusMap_.end(); ++it) {
    it->second.pe_ = 0.;
    it->second.nUsedPMTs_ = it->second.nGoodPMTs_;
    if(it->second.nUsedPMTs_ > 0 && it->second.ifMA_) {
      nGood++; //used to caluculate GoF for now we calculate it using only main array tanks
    }
  }
  
  // Cycle through all hit channels; set sigma=0 so they're initially unused.
  // This interprets unselected channels as no measurement instead of zero
  for (evt::Event::ConstChannelIterator it = event.ChannelsBegin();
       it != event.ChannelsEnd(); ++it) {
    // Ignore channels marked as bad
    if (channelStatusMap.IsGood(it->GetChannelId())){
      --(statusMap_[it->GetTankId()].nUsedPMTs_); 
    }
  }
  
  // Cycle through selected hits, update map.  Make sure not to include
  // multiple hits from a channel.  All hits are presented in the order
  // Tank->Channel->Hit, so multiple hits in a channel must be consecutive.
  int prevChannel = -1;
  
  ConstHitSelector sel(event.HitsBegin(),
                    event.HitsEnd(), GetHitSelection(bag));
  
  for (ConstHitSelectionIterator it = sel.HitSelectionBegin();
    it != sel.HitSelectionEnd(); ++it) {
    
    if(it->channelId_< 1201 && channelStatusMap.IsGood(it->channelId_) &&
       it->channelId_ != prevChannel){
      // Ignore channels marked as bad
      prevChannel = it->channelId_;
      ChannelStatus& status = statusMap_[it->tankId_];
      ++(status.nUsedPMTs_);
      status.pe_ += GetEffectiveCharge(det, *it);
    }
  }
  
  int count = 0;
  ConstORHitSelector selOR(event.ORHitsBegin(),
                           event.ORHitsEnd(), GetORHitSelection(bag));

  for (ConstORHitSelectionIterator it = selOR.HitSelectionBegin();
                                   it != selOR.HitSelectionEnd(); ++it){

      ChannelStatus& status = statusMap_[it->tankId_];
      //++(status.nUsedPMTs_);
      status.pe_ += it->nPE_;
      count++;    
    }
  log_debug("There are " << count << " outrigger hits used for this event in the LHLatDistFit" << std::endl);
  //}

  ////////////////////////////////////////////////////////////////////
  
  // Cycle through the tanks; calculate nFit
  LHLatDistFitResultPtr result = boost::make_shared<LHLatDistFitResult>();
  int nFit = 0;
  for (std::map<int, ChannelStatus>::iterator it = statusMap_.begin();
       it != statusMap_.end(); ++it) {
    if (it->second.nUsedPMTs_ > 0) {
      ++nFit;
      it->second.pe_ /= it->second.nUsedPMTs_;
    }
  }
  
  std::vector<ChannelStatus> data(nFit);
  int idx = 0;
  for (std::map<int, ChannelStatus>::iterator it = statusMap_.begin();
       it != statusMap_.end(); ++it) {
    
    if (it->second.nUsedPMTs_ > 0) {
      data[idx] = it->second;
      ++idx;
    }
  }
  
  // counting number of hit and zeros tanks.
  nHitTanksMA_ = 0;
  nHitTanksOR_ = 0; 
  nZeroTanksMA_ = 0;
  nZeroTanksOR_ = 0;

  for (std::vector<ChannelStatus>::iterator d=data.begin();  d != data.end(); ++d) {
    if(d->pe_> 0.){                     
      d->pe_ = log10(d->pe_);
      if(d->ifMA_) nHitTanksMA_++;
      else if(!d->ifMA_)  nHitTanksOR_++;      
    }
    else{
      d->pe_ = minLogPE_;
      if(d->ifMA_) nZeroTanksMA_++;
      else if(!d->ifMA_) nZeroTanksOR_++;
    }
  }
  
  //trick to use the full data (hits + zeros) when we have more than certain number of hits in the main array
  //otherwise we limit the numer of zeros taken into account in the fit to speed up the fitting procedure
  //this is because we have lot of small nhit showers in data where taking all zeros in the fit makes it as time consuming
  
  data_.clear();
  double xGuess = seed.GetX(); //we will use this SFCF core as seed
  double yGuess = seed.GetY();
  unsigned int counter = 0;

  if(nHitTanksMA_ < minNHitTanksMA_){
    for (std::vector<ChannelStatus>::iterator d=data.begin();  d != data.end(); ++d) {
      if(d->pe_ > minLogPE_) {
	data_.push_back(*d);
      }
      else{
	counter++;
	if(counter == everyNthZero_){
	  d->probWeight_ = everyNthZero_;
	  data_.push_back(*d);
	  counter = 0; 
	}
	else continue;
      }	  
    }
  }
  
  else data_ = data;
  
  result->SetType(CoreFitTypes::LHLatDistFit);
  result->SetXYZ(-1, -1, seed.GetZ());
  result->SetStatus(RECO_FAIL);
  result->SetUncertaintiesCalculated(false);
  result->SetNFit(nFit);
  result->SetNdof(nFit-3);
  result->SetLHLatDistFitXmax(-1.e6);
  result->SetLHLatDistFitEnergy(-1.e6);
  result->SetLHLatDistFitGoF(-1.e6);

#ifdef storeLHLatDistFitExtraVaribales 
  result->SetLHLatDistFitMinLikelihood(1.e6);  
  result->SetLHLatDistFitHitMinLikelihood(1.e6);
  result->SetLHLatDistFitZeroMinLikelihood(1.e6);

  result->SetLHLatDistFitNGoodTanksMA(0); //Good is only defined for MA now now selection for ORs
  result->SetLHLatDistFitNHitTanksMA(0);
  result->SetLHLatDistFitNHitTanksOR(0);
  result->SetLHLatDistFitNZeroTanksMA(0);
  result->SetLHLatDistFitNZeroTanksOR(0);
#endif

  if (nHitTanksMA_ < minimumTanks_ || seed.GetStatus() == RECO_FAIL) {
    bag->Put(resultName_, result);
    return HandleResult(result->GetStatus());
  }

  gMinL_ = gZeroMinL_ = gHitMinL_ = 1e6; //get separate contribution to the likelihood from zeros and hits
  //result variables
  gMinLikelihood_ = 1.e6;
  gBestXc_ = gBestYc_ = gBestXmax_ = gBestEnergy_ = gGoF_ = -1.e6;  
 
  FitProcedure(xGuess,yGuess);
  
  //GoF only works for Main array the code will work with outriggers as well but the model for minLikelihood which we calculate the GoF (coming from simulation is only defined for MA only)
  gGoF_ = CalculateGoF(gMinLikelihood_,nHitTanksMA_,nZeroTanksMA_,nGood);

  result->SetStatus(RECO_SUCCESS);
  result->SetXYZ(gBestXc_, gBestYc_, result->GetZ());
  result->SetLHLatDistFitXmax(gBestXmax_);
  result->SetLHLatDistFitEnergy(gBestEnergy_);
  result->SetLHLatDistFitGoF(gGoF_);
  result->SetLHLatDistFitMinLikelihood(gMinLikelihood_);  
  result->SetLHLatDistFitHitMinLikelihood(gHitMinL_);
  result->SetLHLatDistFitZeroMinLikelihood(gZeroMinL_);
  result->SetLHLatDistFitNGoodTanksMA(nGood);
  result->SetLHLatDistFitNHitTanksMA(nHitTanksMA_);
  result->SetLHLatDistFitNHitTanksOR(nHitTanksOR_);
  result->SetLHLatDistFitNZeroTanksMA(nZeroTanksMA_);
  result->SetLHLatDistFitNZeroTanksOR(nZeroTanksOR_);
  bag->Put(resultName_, result);
  
  // Draw the likelihood surface and LDF PDF of the event when debugging mode is on
 
  if(drawLHSurface_ && display_){
    log_info("LHLatDistFit Resutls: ");
    log_info("Estimated Core X: " << gBestXc_ << " m");
    log_info("Estimated Core Y: " << gBestYc_ << " m");
    log_info("Estimated Energy = " << pow(10,gBestEnergy_)/1.e3 << " TeV");
    log_info("Estimated Xmax = " << gBestXmax_ << " g/cm^2");
    log_info("Goodness of Fit = " << gGoF_ << " for gamma should be about < 3 (doesn't work with outriggers yet)");
    log_info("Drawing Likelihood surface and LDF,PDF for the current event...");
    DrawLikelihoodSurfaceLDFandPDF(xGuess, yGuess);
  }

  return HandleResult(result->GetStatus());
}

//// Draw likelihood surface for reconstructed energy and Xmax and draw the reconstructed core

void LHLatDistFit::DrawLikelihoodSurfaceLDFandPDF(double xGuess, double yGuess){

  //// Plot Full likelihood surface
  canvas->cd(); 
  padSurfaceXY->cd();
  
  if(lhSurface != NULL) delete lhSurface;
  lhSurface = new TH2F("lhSurface","Likelihood surface",39,arrayCenterX_ - 195.,arrayCenterX_ + 195.,39,arrayCenterY_ -195.,arrayCenterY_ + 195.);
  lhSurface->SetDirectory(0);
  for(int xbin = 1; xbin <= lhSurface->GetNbinsX(); xbin++){
    double xpos = lhSurface->GetXaxis()->GetBinCenter(xbin);
    for(int ybin = 1; ybin <= lhSurface->GetNbinsY(); ybin++){
      double ypos = lhSurface->GetYaxis()->GetBinCenter(ybin);
      double P[4]= {gBestXmax_,gBestEnergy_,xpos,ypos};
      double L = LikelihoodFunction(P);
      lhSurface->SetBinContent(xbin,ybin,L);
    }
  }

  lhSurface->SetTitle(";x[m];y[m]");
  lhSurface->SetContour(50);
  lhSurface->Draw("CONT1");

  if(grMATanks != NULL) delete grMATanks;
  grMATanks = new TGraph();
  if(grORTanks != NULL) delete grORTanks;
  grORTanks = new TGraph();
  if(LDFMA != NULL) delete LDFMA;
  LDFMA = new TGraph();
  if(LDFOR != NULL) delete LDFOR;
  LDFOR = new TGraph();

  vector<TEllipse*> el;
  vector<TLine*> zeroLineMA;
  vector<TLine*> zeroLineOR;
  bool hasORs = false;

  for(std::vector<ChannelStatus>::iterator d = data_.begin();  d < data_.end(); ++d){
    double x = d->x_;
    double y = d->y_;
    double dist = sqrt((x-gBestXc_)*(x-gBestXc_)+(y-gBestYc_)*(y-gBestYc_));
    double pe = d->pe_;

    //// Draw charge ellipses on the likelihood surface
    if(pe > minLogPE_){
      float r = 3.0*log10(pow(10,pe)+1.5);
      ellipse = new TEllipse(x, y, r);
      el.push_back(ellipse);
    }

    if(d->ifMA_){
      grMATanks->SetPoint(grMATanks->GetN(),x,y);
      if(pe > minLogPE_){
	LDFMA->SetPoint(LDFMA->GetN(),dist,pe);
      }
      else {
	line = new TLine(dist,-10,dist,1);
	zeroLineMA.push_back(line);
      }
    }
    
    else{
      hasORs = true;
      grORTanks->SetPoint(grORTanks->GetN(),x,y);
      if(pe > minLogPE_){
	LDFOR->SetPoint(LDFOR->GetN(),dist,pe);
      }
      else{
	line = new TLine(dist,-10,dist,1);
	zeroLineOR.push_back(line);
      } 
    }    
  }
    
  grMATanks->SetMarkerStyle(20);
  grMATanks->Draw("p");
  if(hasORs){
    grORTanks->SetMarkerStyle(20);
    grORTanks->SetMarkerColor(kGray+2);
    grORTanks->Draw("p");
  }
  
  //// Draw charge ellipses on top of the tank locations
  for (unsigned int i = 0; i < el.size(); ++i){
    el[i]->SetFillStyle(0);
    el[i]->SetLineWidth(2);
    el[i]->SetLineColor(kMagenta);
    el[i]->Draw();
  }

  if(mSeed != NULL) {
    delete mSeed;
  }
  mSeed = new TMarker(xGuess,yGuess,21);
  mSeed->SetMarkerColor(28);
  mSeed->SetMarkerSize(2);
  mSeed->Draw();
  
  if(mReco != NULL) delete mReco;
  mReco = new TMarker(gBestXc_,gBestYc_,29);
  mReco->SetMarkerColor(3);
  mReco->SetMarkerSize(2.5);
  mReco->Draw();
  
  ////Draw zoomed in likelihood surface around the minimum
  padSurfaceXYZoom->cd();
  if(lhSurfaceZoom != NULL) delete lhSurfaceZoom;
  lhSurfaceZoom = new TH2F("lhSurfaceZoom","Likelihood surface",32,gBestXc_ - 80.,gBestXc_ + 80.,32,gBestYc_ - 80.,gBestYc_ + 80.);
  for(int xbin = 1; xbin <= lhSurfaceZoom->GetNbinsX(); xbin++){
    double xpos = lhSurfaceZoom->GetXaxis()->GetBinCenter(xbin);
    for(int ybin = 1; ybin <= lhSurfaceZoom->GetNbinsY(); ybin++){
      double ypos = lhSurfaceZoom->GetYaxis()->GetBinCenter(ybin);
      double P[4]= {gBestXmax_,gBestEnergy_,xpos,ypos};
      double L = LikelihoodFunction(P);
      lhSurfaceZoom->SetBinContent(xbin,ybin,L);
    }
  }
  lhSurfaceZoom->SetTitle(";x[m];y[m]");
  lhSurfaceZoom->SetContour(70);
  lhSurfaceZoom->Draw("CONT1");
  grMATanks->Draw("p");
  if(hasORs){
    grORTanks->Draw("p");
  }
  mReco->Draw();

  ////Draw legend
  padLeg->cd();
  if(leg != NULL) delete leg;
  leg = new TLegend(0.0,0.0,1.0,1.0);
  leg->AddEntry(grMATanks,"Main array tanks","P");
  if(hasORs){
    leg->AddEntry(grORTanks,"Outrigger tanks","P");
  }
  leg->AddEntry(mReco,"Reconstructed core","P");
  leg->AddEntry(mSeed,"SFCF core","P");
  leg->SetNColumns(4);
  leg->Draw();
  
  ////Draw likelihood surface Energy
  padSurfaceEnergy->cd();
  padSurfaceEnergy->SetLogy();
  if(lhEnergy != NULL) delete lhEnergy;
  lhEnergy = new TH1F("lhEnergy","",330,energyInitialValue_,energyEndValue_);
  for(int enBin = 1; enBin <= lhEnergy->GetNbinsX(); enBin++){
    double energy = lhEnergy->GetBinCenter(enBin);
    double P[4]= {gBestXmax_,energy,gBestXc_,gBestYc_};
    double L = LikelihoodFunction(P);
    lhEnergy->SetBinContent(enBin, L);
  }
  lhEnergy->SetMarkerStyle(20);
  lhEnergy->SetTitle(";log_{10}(E/GeV);Likelihood");
  lhEnergy->GetYaxis()->SetTitleOffset(1.15);
  lhEnergy->GetXaxis()->SetTitleOffset(1.15);
  lhEnergy->Draw("P");

  ////Draw likelihood surface Xmax
  padSurfaceXmax->cd();
  padSurfaceXmax->SetLogy();
  if(lhXmax != NULL) delete lhXmax;
  lhXmax = new TH1F("lhXmax","",600,300.0,xmaxEndValue_);
  for(int xmBin = 1; xmBin <= lhXmax->GetNbinsX(); xmBin++){
    double Xmax = lhXmax->GetBinCenter(xmBin);
    double P[4]= {Xmax,gBestEnergy_,gBestXc_,gBestYc_};
    double L = LikelihoodFunction(P);
    lhXmax->SetBinContent(xmBin, L);
  }
  lhXmax->SetMarkerStyle(20);
  lhXmax->SetTitle(";Xmax[g/cm^{2}];Likelihood");
  lhXmax->GetYaxis()->SetTitleOffset(1.15);
  lhXmax->Draw("P");
  canvas->SaveAs("LikelihoodSurface.pdf"); 

  ////////////////////////////////////////////////////////////////////
  /////Here we draw LDF and PDF

  ////Get the best fit PDF bins
  int xmaxHistNum = (gBestXmax_ - xmaxInitialValue_)/xmaxBinSize_;
  int energyHistNum = (gBestEnergy_ -  energyInitialValue_)/energyBinSize_;

  LDFMA->SetMarkerStyle(20);
  LDFOR->SetMarkerStyle(20);

  ////Here we draw the MA LDF PDF
  canvasPDF->cd();
  pad->cd();
  pad->SetTicks(1,1);
  
  logprobMA = (TH2F*)histLogprobMainArray[xmaxHistNum][energyHistNum][zenithBin_]->Clone();
  logprobMA->UseCurrentStyle();
  logprobMA->GetYaxis()->SetRangeUser(-1.7,5.2);
  logprobMA->GetXaxis()->SetRangeUser(0,350);
  logprobMA->GetYaxis()->SetTitleOffset(0.8);
  logprobMA->GetYaxis()->SetTitle("log_{10}(N_{pe})");
  logprobMA->GetZaxis()->SetTitleOffset(1.0);
  logprobMA->GetZaxis()->SetTitle("log(P)");
  logprobMA->Draw("colz");
  LDFMA->Draw("p");

  zeroPad->cd();
  int zeroBin = minLogPEBin_;
  zeroLike = (TH1F*)logprobMA->ProjectionX("px",zeroBin,zeroBin);
  zeroLike->SetTitle("");
  zeroLike->GetYaxis()->SetRangeUser(-10,1);
  zeroLike->GetYaxis()->SetLabelSize(0.08);
  zeroLike->GetYaxis()->SetTitle("log(P_{0})");
  zeroLike->GetYaxis()->SetTitleOffset(0.43);
  zeroLike->GetYaxis()->SetTitleSize(0.08);

  zeroLike->GetXaxis()->SetTitleSize(0.08);
  zeroLike->GetXaxis()->SetTitleOffset(1.0);
  zeroLike->GetXaxis()->SetTitle("r[m]");
  zeroLike->GetXaxis()->SetRangeUser(0,350);
  zeroLike->GetXaxis()->SetTickSize(0.06);
  zeroLike->GetXaxis()->SetLabelSize(0.08);

  zeroLike->Draw();
  
  for(unsigned int i = 0; i < zeroLineMA.size(); i++){
    zeroLineMA[i]->SetLineStyle(2);
    zeroLineMA[i]->Draw();
  }

  ////Here we draw the OR LDF and PDF
  padOR->cd();
  padOR->SetTicks(1,1);

  logprobOR = (TH2F*)histLogprobOutriggers[xmaxHistNum][energyHistNum][zenithBin_]->Clone();
  logprobOR->UseCurrentStyle();
  logprobOR->GetYaxis()->SetRangeUser(-1.7,5.2);
  logprobOR->GetXaxis()->SetRangeUser(0,350);
  logprobOR->GetYaxis()->SetTitleOffset(0.8);
  logprobOR->GetYaxis()->SetTitle("log_{10}(N_{pe})");
  logprobOR->GetZaxis()->SetTitleOffset(1.0);
  logprobOR->GetZaxis()->SetTitle("log(P)");
  logprobOR->Draw("colz");
  if(hasORs){
    LDFOR->Draw("p");
  }
  zeroPadOR->cd();
  zeroLikeOR = (TH1F*)logprobOR->ProjectionX("px",zeroBin,zeroBin);
  zeroLikeOR->SetTitle("");
  zeroLikeOR->GetYaxis()->SetRangeUser(-10,1);
  zeroLikeOR->GetYaxis()->SetLabelSize(0.08);
  zeroLikeOR->GetYaxis()->SetTitleSize(0.08);
  zeroLikeOR->GetYaxis()->SetTitle("log(P_{0})");
  zeroLikeOR->GetYaxis()->SetTitleOffset(0.43);
  
  zeroLikeOR->GetXaxis()->SetTitleSize(0.08);
  zeroLikeOR->GetXaxis()->SetTitle("r[m]"); 
  zeroLikeOR->GetXaxis()->SetTitleOffset(1.0);
  zeroLikeOR->GetXaxis()->SetTickSize(0.06);
  zeroLikeOR->GetXaxis()->SetRangeUser(0,350);
  zeroLikeOR->GetXaxis()->SetLabelSize(0.08);
 
  zeroLikeOR->Draw();
  if(hasORs){
    for(unsigned int i = 0; i < zeroLineOR.size(); i++){
      zeroLineOR[i]->SetLineStyle(2);
      zeroLineOR[i]->Draw();
    }
  }
  canvasPDF->SaveAs("LDF_and_PDF.pdf");
  return;
}

double LHLatDistFit::InterpolateXmaxEnergy( TH2F* hist, int xBin, int yBin, double x, double y){ //this function was used for the prior between Xmax and energy not used anymore!
  
  double prob00,prob01,prob10,prob11;
  prob00 = prob01 = prob10 = prob11 = 0.0;
  int xBin0, yBin0, xBin1, yBin1;
  xBin0 = yBin0 = xBin1 = yBin1 = 0;
  
  xBin0 = xBin;
  yBin0 = yBin;
  prob00 = (double)hist->GetBinContent(xBin0,yBin0);
  
  double x0 = (double)hist->GetXaxis()->GetBinCenter(xBin0);
  double y0 = (double)hist->GetYaxis()->GetBinCenter(yBin0);
  
  if(x0 > x && xBin0 > 1){
    xBin1 = xBin0 - 1;
    prob01 = (double)hist->GetBinContent(xBin1,yBin0);
  }
  
  else if(x0 < x && xBin0 < hist->GetNbinsX()){
    xBin1 = xBin0 + 1;
    prob01 = (double)hist->GetBinContent(xBin1,yBin0);
  }
  
  else if(xBin0 == 1 || xBin0 == hist->GetNbinsX()){
    xBin1 = xBin0;
    prob01 = (double)hist->GetBinContent(xBin1,yBin0);
  }
  
  if(y0 > y && yBin0 > 1 ){
    yBin1 = yBin0 - 1;
    prob10 = (double)hist->GetBinContent(xBin0,yBin1);
    prob11 = (double)hist->GetBinContent(xBin1,yBin1);
  }
  
  else if(y0 < y && yBin0 < hist->GetNbinsY()){
    yBin1 = yBin0 + 1;
    prob10 = (double)hist->GetBinContent(xBin0,yBin1);
    prob11 = (double)hist->GetBinContent(xBin1,yBin1);
  }
  
  else if(yBin0 == 1 || yBin0 == hist->GetNbinsY()){
    yBin1 = yBin0;
    prob10 = (double)hist->GetBinContent(xBin0,yBin1);
    prob11 = (double)hist->GetBinContent(xBin1,yBin1);
  }
  
  double x1 = (double)hist->GetXaxis()->GetBinCenter(xBin1);
  double y1 = (double)hist->GetYaxis()->GetBinCenter(yBin1);
  
  double prob_temp0 = InterpolateTwoPoints(x0,prob00,x1,prob01,x);
  double prob_temp1 = InterpolateTwoPoints(x0,prob10,x1,prob11,x);
  double prob = InterpolateTwoPoints(y0,prob_temp0,y1,prob_temp1,y);
  return prob;
}

double LHLatDistFit::InterpolateInsideHist( TH2F* hist, double x, double y){ //function to do interpolation within a histogram which we are not doing for pass5 for speed reasons
  // this function is replaced now with GetProbFromHist
  int xBin, yBin;
  double prob00, prob01;
  prob00 = prob01 = 0.0;
  
  xBin = (int)hist->GetXaxis()->FindBin(x);
  
  if(y > minLogPE_) yBin = (int)hist->GetYaxis()->FindBin(y);
  else yBin = minLogPEBin_; // zeros stored at the log scale at minLogPE_
  
  prob00 = (double)hist->GetBinContent(xBin , yBin);
  return prob00;

  double x0 = (double)hist->GetXaxis()->GetBinCenter(xBin);
  double x1 = x0;
  
  
  if(xBin > 1 && xBin < hist->GetNbinsX()){
    if(x0 > x){
      int xBinNum = xBin - 1;
      prob01 = (double)hist->GetBinContent(xBinNum, yBin);
      x1 = (double)hist->GetXaxis()->GetBinCenter(xBinNum);
    }
    
    else {
      int xBinNum = xBin + 1;
      prob01 = (double)hist->GetBinContent(xBinNum, yBin);
      x1 = (double)hist->GetXaxis()->GetBinCenter(xBinNum);
    }
    double prob = InterpolateTwoPoints(x0,prob00,x1,prob01,x);

    return prob;
  }
  
  else{
    double prob = prob00;
    return prob;
  }
}
