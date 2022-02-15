#include <core-fitter/LHLatDistFitLookUpTableMaker.h>
#include <hawc-reco/hit-selection/HitSelector.h>
#include <hawc-reco/hit-selection/HitSelectorOR.h>
#include <data-structures/reconstruction/Reco.h>
#include <data-structures/reconstruction/core-fitter/CoreFitResult.h>
#include <data-structures/event/SimEvent.h>
#include <data-structures/event/Event.h>
#include <data-structures/event/ChannelEvent.h>
#include <data-structures/detector/Detector.h>
#include <detector-service/DetectorService.h>

#include <hawcnest/HAWCNest.h>
#include <hawcnest/HAWCUnits.h>
#include <hawcnest/RegisterService.h>

#include <cstring>
#include <string>
#include <iterator>

using namespace std;
using namespace HAWCUnits;
REGISTER_SERVICE(LHLatDistFitLookUpTableMaker);

Configuration LHLatDistFitLookUpTableMaker::DefaultConfiguration() {

  Configuration config = RecoBase::DefaultConfiguration(); 
  config.Parameter<std::string>("simEventName");
  config.Parameter<std::string>("inputRootFile");
  config.Parameter<std::string>("outputRootFile");
  config.Parameter<bool>("outriggers");
  config.Parameter<bool>("hasInputFile",false);
  
  return config;
}

void LHLatDistFitLookUpTableMaker::Initialize(const Configuration& config) {
  
  log_info("Initializing LookupTableGenerator.....");

  RecoBase::Initialize(config);
  config.GetParameter("simEventName",simEventName_);
  config.GetParameter("inputRootFile",inputRootFile_);
  config.GetParameter("outputRootFile",outputRootFile_);
  config.GetParameter("outriggers",outriggers_); 
  costh_ = 1.0, sinth_ = 0.0, cosph_ = 1.0, sinph_ = 0.0;
  D2R_ = 0.01745329;                           // pi/180             
  
  filein_ = new TFile(inputRootFile_.c_str());
  if(filein_->IsZombie() || !filein_->IsOpen()){
    log_info("The file [" << filein_ << "] does not exist ==> Is this the first Iteration?");
  }
  else hasInputFile_ = true;

  /////////////////for data set information for the tree_////////////////////
  tree_ = new TTree("tree","tree"); 
  tree_->Branch("simEnergy", &simEnergy_, "simEnergy/D");
  tree_->Branch("simZang", &simZang_, "simZang/D");
  tree_->Branch("simAang", &simAang_, "simAang/D");
  tree_->Branch("simXmax", &simXmax_, "simXmax/D");
  tree_->Branch("simXCore", &simXCore_, "simXCore/D");
  tree_->Branch("simYCore", &simYCore_, "simYCore/D");
  tree_->Branch("xTank", &xTank_);
  tree_->Branch("yTank", &yTank_);
  tree_->Branch("peTank", &peTank_);
  tree_->Branch("timeTank",&timeTank_); 
  tree_->Branch("typeTank",&typeTank_); 
  tree_->Branch("dist",&dist_);
  /////////////////for data set information for the tree_////////////////////

  fileout_ = new TFile(outputRootFile_.c_str(),"recreate");
  tree_->SetDirectory(fileout_);
  
  for(int i=0; i<numXmaxBins_; i++){ //Xmax loop
    for(int j=0; j<numEnergyBins_; j++){ //Energy loop
      for(int k=0; k<numZenithBins_; k++){ // Zenith angle loop
	char name_hist[100];
	
	sprintf(name_hist,"logprob_%d_%d_%d",i,j,k);
	logprob[i][j][k]= new TH2F(name_hist,"",numRadialDistBins_,0.0,500.0,numLogPEBins_,-4.0,8.0);

	sprintf(name_hist,"logrho_dist_%d_%d_%d",i,j,k);
	logrhoDist[i][j][k]= new TH2F(name_hist,"",numRadialDistBins_,0.0,500.0,numLogPEBins_,-4.0,8.0);
	if(hasInputFile_){
	  sprintf(name_hist,"newlogprob_%d_%d_%d",i,j,k);
	  newLogprob[i][j][k]= new TH2F(name_hist,"",numRadialDistBins_,0.0,500.0,100.0,-5.0,5.0);

	  sprintf(name_hist,"MeanLogA_r_%d_%d_%d",i,j,k);
	  meanLogAR[i][j][k] =  (TH1F*)filein_->Get(name_hist);
	  meanLogAR[i][j][k]->SetDirectory(0);

	  sprintf(name_hist,"SigmaLogA_r_%d_%d_%d",i,j,k);
	  sigmaLogAR[i][j][k] =  (TH1F*)filein_->Get(name_hist);
	  sigmaLogAR[i][j][k]->SetDirectory(0);

	  sprintf(name_hist,"newlogrho_dist_%d_%d_%d",i,j,k);
	  newLogrhoDist[i][j][k]= new TH2F(name_hist,"",numRadialDistBins_,0.0,500.0,100,-5.0,5.0);
	}
      }
    }
  }
  filein_->Close();
  delete filein_;
}

void LHLatDistFitLookUpTableMaker::InitGeo(const det::Detector& det,
			const ChannelStatusMap& channelStatusMap, float simZang, float simAang, float simXCore, float simYCore) {

  costh_ = cos(simZang);
  sinth_ = sin(simZang);
  cosph_ = cos(simAang);
  sinph_ = sin(simAang);
  
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
    status.x_ =  it->GetPosition().GetX() - simXCore;
    status.y_ =  it->GetPosition().GetY() - simYCore;
    status.time_ = 0.; //vikas
    if(it->GetTankId() < 301){ //main array tanks = 300
      status.ifMA_      = true; //flag main array tanks
    }
    else status.ifMA_ = false; //flag outriggers
    statusMap_[it->GetTankId()] = status;
  }
}

Module::Result LHLatDistFitLookUpTableMaker::Process(BagPtr bag) {

  const SimEvent& simEvent  = bag->Get<SimEvent>(simEventName_);
  float simZang = simEvent.GetEventHeader().theta_; //values are in radians
  float simAang = simEvent.GetEventHeader().phi_;  //values are in radians
  float simXCore = simEvent.GetEventHeader().xcoreDet_;
  float simYCore = simEvent.GetEventHeader().ycoreDet_;
  float simEnergy = simEvent.GetEventHeader().energy_/GeV;
  float simXmax = (simEvent.GetEventHeader().xmax_*cm*cm/g)/(cos(simZang)); //getting the true Xmax 
 
  // ///////for data set information////////////////////
  simEnergy_ = simEnergy;
  simXmax_ = simXmax;
  simAang_ = simAang;
  simZang_ = simZang;
  simXCore_ = simXCore;
  simYCore_ = simYCore;
  /////////////////////////////////////////////////////////  


  const evt::Event& event  = bag->Get<evt::Event>(event_);
  
  //  Load geometry
  DetectorService& detSv   = GetService<DetectorService>(detectorService_);
  const det::Detector& det = detSv.GetDetector(event.GetTime());
  
  // Load bad channel map 
  const ChannelStatusMap& channelStatusMap = detSv.GetChannelStatusMap(event.GetTime());
  InitGeo(det, channelStatusMap, simZang, simAang, simXCore, simYCore);
  
  // Set all tanks that have good channels as good
  for (std::map<int, ChannelStatus>::iterator it = statusMap_.begin();
       it != statusMap_.end(); ++it) {
    it->second.pe_ = 0.;
    it->second.nUsedPMTs_ = it->second.nGoodPMTs_;
  }
  
  // Cycle through all hit channels; set sigma=0 so they're initially unused.
  // This interprets unselected channels as no measurement instead of zero
  for (evt::Event::ConstChannelIterator it = event.ChannelsBegin();
       it != event.ChannelsEnd(); ++it) {
    // Ignore channels marked as bad
    if (channelStatusMap.IsGood(it->GetChannelId())) {
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
    
    //Ignore channels marked as bad
    if (it->channelId_< 1201 && channelStatusMap.IsGood(it->channelId_) &&
        it->channelId_ != prevChannel) {
      
      prevChannel = it->channelId_;
      
      log_debug("channel: " << it->channelId_ << " time: "<< it->calibData_.time_ << " charge: " << it->calibData_.PEs_);
      ChannelStatus& status = statusMap_[it->tankId_];
      ++(status.nUsedPMTs_);
      status.pe_ += GetEffectiveCharge(det, *it);
      status.time_ += it->calibData_.time_;
    }   
  }

  int count = 0;  
  ConstORHitSelector selOR(event.ORHitsBegin(),
			   event.ORHitsEnd(), GetORHitSelection(bag));
  
  for (ConstORHitSelectionIterator it = selOR.HitSelectionBegin();
       it != selOR.HitSelectionEnd(); ++it){
    
    ChannelStatus& status = statusMap_[it->tankId_];
    status.pe_ += it->nPE_;
    count++;
  }
  log_debug("There are " << count << " outrigger hits used for this event in the LHLatDistFit" << std::endl);
  
  //  Cycle through the tanks; calculate nFit
  
  int nFit = 0;
  for (std::map<int, ChannelStatus>::iterator it = statusMap_.begin();
       it != statusMap_.end(); ++it) {
    
    if (it->second.nUsedPMTs_ > 0) {
      ++nFit;
      it->second.pe_ /= it->second.nUsedPMTs_;
      it->second.time_ /= it->second.nUsedPMTs_;
    }
  }
  
  std::vector<ChannelStatus> data(nFit);
  
  int idx = 0;

  /////////////////for data set information for the tree_////////////////////
  xTank_.clear(); 
  yTank_.clear(); 
  peTank_.clear();
  timeTank_.clear();
  typeTank_.clear();
  dist_.clear();
  /////////////////for data set information for the tree_////////////////////

  for (std::map<int, ChannelStatus>::iterator it = statusMap_.begin();
       it != statusMap_.end(); ++it) {
   
    if (it->second.nUsedPMTs_ > 0) {

      /////////////////for data set information for the tree_////////////////////
      xTank_.push_back(it->second.x_);
      yTank_.push_back(it->second.y_);
      /////////////////for data set information for the tree_////////////////////

      data[idx] = it->second;
      ++idx;
    }
  }

  if (simXmax < xmaxEndValue_ && simXmax > xmaxInitialValue_ && simEnergy < energyEndValue_ && simEnergy > energyInitialValue_ && simZang < 0.872665){ // 50deg to radian  = 0.872665
    int xmaxBin = (simXmax-xmaxInitialValue_)/xmaxBinSize_;
    int energyBin = log10(simEnergy/energyInitialValue_)/energyBinSize_;
    int zAngBin = (1-cos(simZang))/ZenithBinSize_; 

    dist_.clear();
    for (std::vector<ChannelStatus>::iterator d=data.begin();  d < data.end(); ++d) {
      float distShPlane = sqrt((d->x_)*(d->x_)+(d->y_)*(d->y_)-sinth_*sinth_*(d->x_*cosph_+d->y_*sinph_)*(d->x_*cosph_+d->y_*sinph_));

      peTank_.push_back(d->pe_);
      timeTank_.push_back(d->time_);
      typeTank_.push_back(d->ifMA_);
      dist_.push_back(distShPlane);
      if(!outriggers_ && d->ifMA_){
	if(d->pe_ > 0){
	  logrhoDist[xmaxBin][energyBin][zAngBin]->Fill(distShPlane,log10(d->pe_));
	  if(hasInputFile_){
	    int binNumber = meanLogAR[xmaxBin][energyBin][zAngBin]->FindBin(distShPlane);
	    double meanLogA = meanLogAR[xmaxBin][energyBin][zAngBin]->GetBinContent(binNumber);
	    double sigma = sigmaLogAR[xmaxBin][energyBin][zAngBin]->GetBinContent(binNumber);
	    
	    double sigmaDist = 0.0;
	    if( sigma > 0.00001){
	      sigmaDist = (log10(d->pe_)-meanLogA)/sigma;
	    }
	    newLogrhoDist[xmaxBin][energyBin][zAngBin]->Fill(distShPlane,sigmaDist); 
	  }
	}
	else logrhoDist[xmaxBin][energyBin][zAngBin]->Fill(distShPlane,-3.0);	
      }
      else if(outriggers_ && !d->ifMA_){
	if(d->pe_ > 0){
	  logrhoDist[xmaxBin][energyBin][zAngBin]->Fill(distShPlane,log10(d->pe_));
	  if(hasInputFile_){
	    int binNumber = meanLogAR[xmaxBin][energyBin][zAngBin]->FindBin(distShPlane);
	    double meanLogA = meanLogAR[xmaxBin][energyBin][zAngBin]->GetBinContent(binNumber);
	    double sigma = sigmaLogAR[xmaxBin][energyBin][zAngBin]->GetBinContent(binNumber);
	    
	    double sigmaDist = 0.0;
	    if( sigma > 0.00001){
	      sigmaDist = (log10(d->pe_)-meanLogA)/sigma;
	    }
	    newLogrhoDist[xmaxBin][energyBin][zAngBin]->Fill(distShPlane,sigmaDist); 
	  }
	}
	else logrhoDist[xmaxBin][energyBin][zAngBin]->Fill(distShPlane,-3.0);			
      }
      else{
	log_debug("This should not happen!!!!,Two cases: either outirggers or main array templates will be made. Something doesn't look right check if you have used correct flag in lookuptable-generator --outrigger-templates");
      }
    }
  }
  fileout_->cd();
  /////////////////for data set information for the tree_////////////////////
  tree_->Fill();
  /////////////////for data set information for the tree_////////////////////
  data.clear();
  
  CoreFitResultPtr result = boost::make_shared<CoreFitResult>();
  result->SetType(CoreFitTypes::LHLatDistFitLookUpTableMaker);
  result->SetStatus(RECO_SUCCESS);
  bag->Put(resultName_, result);
  return HandleResult(result->GetStatus());
}
