/*!
 * @file offline-reconstructor.cc
 * @brief Reconstruct triggered XCD files with TDC data.
 * @author Jim Braun
 * @date 8 Nov 2012
 * @version $Id$
 */
#include <hawc-calibration/calibrator/ChargeCalibrationService.h>
#include <hawc-calibration/calibrator/TimingCalibrationService.h>
#include <hawc-calibration/calibrator/Calibrator.h>
#include <or-calibration/ORTimingCalibrator.h>
#include <or-calibration/ORSignalReconstructor.h>

#include <aerie-io/AERIEFileType.h>
#include <aerie-io/reader/Reader.h>
#include <aerie-io/triggered/TriggeredInputSelector.h>
#include <aerie-io/mc-analysis/MCAnalysisDataSelector.h>
#include <aerie-io/mc-analysis/MCAnalysisInputSelector.h>
#include <aerie-io/writer/BinaryWriter.h>
#include <aerie-io/serializer/DynamicSerializer.h>
#include <aerie-io/hawcsim/HAWCSimInputSelector.h>
#include <aerie-io/outrigger/OutriggerInputSelector.h>
#include <aerie-io/reader/MergerORReader.h>

#include <sim-pmt-modeler/SimulatedPMTModeler.h>
#include <sim-pmt-modeler/DAQSim.h>

#include <trig-mc/SummarizeMCFields.h>
#include <trig-mc/SummarizeSimFields.h>

#include <hawc-reco/hit-selection/HitSelectionRoutines.h>
#include <hawc-reco/hit-selection/HitSelectionRoutinesOR.h>
#include <hawc-reco/CalibrationEventFilter.h>
#include <hawc-reco/ThresholdFilter.h>
#include <hawc-reco/NhitCalculator.h>
#include <hawc-reco/EdgeRefiner.h>
#include <core-fitter/COMCoreFit.h>
#include <core-fitter/GaussCoreFit.h>
#include <core-fitter/TankGaussCoreFit.h>
#include <core-fitter/SFCF.h>
#include <core-fitter/NKGCoreFit.h>
#include <core-fitter/TankLMCoreFit.h>
#include <core-fitter/ChargeResidualGenerator.h>
#include <core-fitter/PMTChi2.h>
#include <core-fitter/LatDist.h>
#include <gamma-filter/CompactnessCalculator.h>
#include <gamma-filter/PairCompactnessCalculator.h>
#include <gamma-filter/PINCCalculator.h>
#include <gamma-filter/TankGHSep.h>
#include <gamma-filter/SOSCalculator.h>
#include <gamma-filter/FeaturesNNCalculator.h>
#include <gamma-filter/addgbt.h>
#include <hawc-reco/curvature-model/TweakedCurvatureService.h>
#include <hawc-reco/curvature-model/ConfigurableCurvatureService.h>
#include <hawc-reco/curvature-model/HawcCurvatureService.h>
#include <hawc-reco/curvature-model/CrabCurvatureService.h>
#include <track-fitter/GaussPlaneFit.h>
#include <track-fitter/LHAngleFit.h>
#include <track-fitter/PDFAngleFit.h>
#include <track-fitter/ZenithAlignment.h>
#include <track-fitter/TimeResidualGenerator.h>
#include <track-fitter/MPFEventSplitter.h>
#include <detector-service/ConfigDirDetectorService.h>
#include <detector-service/StdDetectorService.h>
#include <astro-service/StdAstroService.h>
#include <fiducial-charge-calc/FiducialCharge.h>
#include <fiducial-core-calc/FiducialCore.h> 

#include <summarize/SummarizeRec.h>
#include <summarize/SummarizeMC.h>

#include <energy-estimator/LHEnergyEstimator.h>
#include <rng-service/StdRNGService.h>

#include <energy-estimator/NeuralNetEnergyEstimator.h>

#include <data-structures/event/Event.h>
#include <data-structures/detector/Detector.h>

#include <hawcnest/HAWCUnits.h>
#include <hawcnest/HAWCNest.h>
#include <hawcnest/Logging.h>
#include <hawcnest/CommandLineConfigurator.h>
#include <hawcnest/SoftwareVersion.h>
#include <hawcnest/processing/SequentialMainLoop.h>
#include <hawcnest/xml/XMLReader.h>
#include <hawcnest/NestIniConfig.h>

#include <diagnostics/TimeResidualDiagnostic.h>
#include <diagnostics/TankLightDistribution.h>
#include <diagnostics/MuonFinder.h>
#include <diagnostics/ORDiagnostics.h>

#include <hawc-display/root-display/HAWCDisplay.h>
#include <TApplication.h>
#include <core-fitter/LHLatDistFit.h>

#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <dlfcn.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>

using namespace HAWCUnits;
using namespace std;

typedef vector<string> FileList;


void parseParamTweak(CommandLineConfigurator& cl,
                     HAWCNest& nest){
#if 0
  vector<string> parsedParameters =
  cl.GetArgument<vector<string> >("param-tweak");
  for(unsigned i = 0 ; i < parsedParameters.size() ; i++){
    vector<string> tokens =
    ConfigurationUtil::tokenize(parsedParameters[i],':');
    if(tokens.size() != 3){
      log_fatal("error parsing param-tweak: '"<<parsedParameters[i]<<"'");
    }
    log_info("tweaking service "<<parsedParameters[0]<<"'s parameter "
             <<parsedParameters[1]<<" to "<<parsedParameters[2]);
    nest.SetParameter_decode(parsedParameters[0],
                             parsedParameters[1],
                             parsedParameters[2]);
  }
#endif
}

int main(int argc, char** argv)
{

  CommandLineConfigurator cl("Reconstruct real and simulated events.");
  
  // Input parameters common to data and MC
  cl.AddPositionalOption<FileList>(
                                   "input",
                                   "Input file name(s)");
  cl.AddOption<string>(
                       "output,o",
                       "offlineOutput.xcd",
                       "Output file (XCD format)");
  cl.AddOption<string>(
                       "config-summary",
                       "",
                       "summary of hawcnest configuration");
  cl.AddOption<string>(
      "usermod",
      "Config file to dynamically introduce user modules");
  cl.AddOption<int>(
                    "maxevents",
                    -1,
                    "Maximum number of events to reconstruct, <=0 for all.");
  
  // Core fit options
  OptionGroup& cf = cl.AddOptionGroup("Core Fit Options");
  cf.AddOption<string>(
                       "core-fit", "sfcf",
                       "Set the core fit type. Options: "
                       "sfcf, gauss, gauss-tank, nkg, nkg-tank, nkg-log.");
  cf.AddOption<int>(
                    "sfcf-strategy", 1,
                    "Set strategy when using SFCF core fit (0=bad, 1=default)");
  cf.AddFlag(
             "sfcf-disable-or",
             "Disable use of outrigger tanks in SFCF core fit");
  cf.AddFlag(
             "sfcf-enable-true-core",
             "Use true core location in SFCF core fit (sim only)");
  cf.AddFlag(
             "no-core-fit-in-plane",
             "Use a single step core fit, instead of doing the core fit with only "
             "hits in time with the plane fit");
  cf.AddFlag(
      "pmtChi2",
      "Run individual PMT chi2 calculation after core fit");
  cf.AddOption<int>(
      "LatDistQErr",
       3,
      "LatDist fitter Charge error method version. Available versions same as PINC, 2 or 3 are valid settings, 3 is for Pass Five, 2 for Pass 4");
  cf.AddOption<int>(
      "GPVersion",
       3,
      "Ground parameter energy version to calculate (version 2 and 3 are valid settings");

  // Angle fit options
  OptionGroup& af = cl.AddOptionGroup("Angle Fit Options");
  af.AddOption<int>(
                    "exclude-channel-from-track-fit,x", -1,
                    "Exclude a channel from the track fit");
  af.AddOption<string>(
                       "exclude-channels-from-track-fit-file", "",
                       "File containing PMT Id list to exclude in track fit");
  af.AddOption<double>(
                       "track-fit-min-npe", -1.,
                       "Minimum charge to be used for track fit, in PEs");
  af.AddFlag(
             "no-remove-max-charge",
             "Do not discard hits larger than the max charge in the calibration");
  af.AddOption<int>(
                    "curvature-model",1,
                    "Curvature: 0:vanilla HAWC, 1:tweaked HAWC, 2:Crab based curvature, 3:Configured by file");
  af.AddOption<double>(
                       "tweaked-curvature-lq1", 3.02,
                       "Linear log(q) corection");
  af.AddOption<double>(
                       "tweaked-curvature-lq2", -0.63,
                       "Quadratic log(q) corection");
  af.AddOption<double>(
                       "tweaked-curvature-r1", 0.,
                       "linear r correction");
  af.AddOption<double>(
                       "tweaked-curvature-r2", 0.,
                       "quadratic r correction");
  af.AddOption<double>(
                       "tweaked-curvature-rlq", 0.,
                       "quadratic 'cross term' correction");
  af.AddOption<int>(
                       "gauss-plane-curv-recompute-iter", 0,
                       "Iteration of the Gaussian plane fitter after which the curvature is recomputed in the shower plane. 0 = not done. (default=0)");
  af.AddOption<string>(
                       "Configurable-Curvature-File", "reconstruction/curvature/Curve_v1.txt",
                       "Filename of Curvature definition file in the CONIG_HAWC directory");
  af.AddOption<double>(
                       "time-residual-cut",
                       -1.,
                       "Cut on time residuals in the shower axis fit [ns]");
  af.AddFlag(
             "evenodd",
             "Do even/odd plane fit");
  af.AddFlag(
             "no-axis",
             "Turn off the calculation of lateral distribution with distance from "
             "axis instead of core. "
             "Performed after core fit, saves fit values to rec output");
  af.AddFlag(
             "enable-likelihood",
             "Perform the LH or PDF angle fit.");
  af.AddFlag(
             "pdffit",
             "Use PDFFit angle fit");
  af.AddFlag(
             "no-splitter",
             "Do not use MPFEventSplitter");
  
  // Data diagnostics options
  OptionGroup& di = cl.AddOptionGroup("Data Diagnostics");
  di.AddOption<string>(
                       "diagnostics",
                       "",
                       "Diagnostics output file "
                       "Only activated if set.");
  di.AddOption<double>(
                       "diagnostics-fhit",
                       0.5,
                       "NHit threshold for some diagnostics, denoting a large event");
  // Muon finder
  di.AddOption<string>(
                       "muonfinder-file",
                       "",
                       "Output XCD file for muonlike tank hits");
  
  
  // Gamma-hadron settings
  OptionGroup& gh = cl.AddOptionGroup("Gamma/Hadron Calculations");
  gh.AddOption<int>(
      "pincness",
       3,
      "PINCCalculator switch: 0: do not run, 1: Pretz, 2: AJS 2015, 3: AJS 2017.");
  gh.AddOption<double>(
                       "pincRadius",
                       0,
                       "PINCCalculator exclusion radius (m)");
  gh.AddFlag(
             "pairCompactness",
             "Run PairCompactnessCalculator.");
  //turn on extra CxPE calculation
  gh.AddFlag(
      "allCxPERadii",
      "Calculate CxPE for 20, 30, and 50m.");
  gh.AddFlag(
      "no-gbt",
      "Do not calculate gbt probability (This requires SOS and lh-lat-dist-fit)");
  gh.AddOption<string>("Pinput-file","",
                       "Proton/Data histograms of tank light distribution -- default in $CONFIG_HAWC");
  gh.AddOption<string>("Ginput-file","",
                       "Gamma histograms of tank light distribution -- default in $CONFIG_HAWC");
  gh.AddOption<string>("GHoutput-file","",
                       "Debug root output file (containing a TTree) of the TankGHSep Service");
  gh.AddOption<string>("GHHistogram-file","",
                       "Debug root output, giving renormalized light distributions from the TankGHSep service");
  gh.AddMultiOption< std::vector<double> >("PE-Cut","12 members: First 10 members are sequentially for bins 0-9. Member 11 is same as bin 9 but limit of the bin is 101% instead of 100%. Member 12 is all events without binning. The cuts must be entered separated by space");
  //SOS options
  gh.AddOption<string>("sos-table-file","","generate output for SOS module");
  gh.AddOption<string>("sos-template-file","","input for SOS module");
  gh.AddOption<string>("sos-debug-file","","debug output for likelihood values SOS module");
  gh.AddFlag("disable-sos","enable the SOS (or not)");
  
  // Output settings
  OptionGroup& ws = cl.AddOptionGroup("Output settings");
  ws.AddFlag(
             "extended",
             "Write both extended output and summary rec");
  ws.AddFlag(
             "no-correct-curv-tr",
             "Don't correct curvature-Time Residuals in Dynamic Serializer ");
  
  // Display settings
  OptionGroup& ds = cl.AddOptionGroup("Display settings");
  ds.AddFlag(
             "display",
             "Display events in viewer");
  ds.AddOption<double>(
                       "displaysleep",
                       -2.,
                       "If --display, sets the time delay between events, in seconds. "
                       "If >0, auto update is active at program start. Else, inactive.");
  
  // Energy estimator settings
  OptionGroup& es = cl.AddOptionGroup("Energy estimators");
  es.AddFlag(
             "no-neural-net-energy,n",
             "Disable neural-net energy estimator.");
  es.AddFlag(
             "no-lh-energy",
             "Disable likelihood energy estimator.");
  es.AddOption<string>(
                       "table-config-file,E",
                       "",
                       "Likelihood Energy Estimator XML configuration file name");
  es.AddFlag(
             "no-features-nn",
             "Disable FeaturesNN calculator");
  
  // Event cuts
  OptionGroup& cs = cl.AddOptionGroup("Event cuts");
  cs.AddOption<int>(
      "threshold",
      0,
      "Minimum number of channels to include event");
  cs.AddOption<int>(
      "wthreshold",
      0,
      "Minimum number of hits within a 150ns window");

  // Input parameters specific to data
  OptionGroup& da = cl.AddOptionGroup("Data Configs");
  
  da.AddOption<string>("config-dir,c",
                       "",
                       "Detector + calib config folder. "
                       "Default: $CONFIG_HAWC");
  da.AddOption<string>("charge-calfit-file,f",
                       "",
                       "PE Calibration file or directory name");
  da.AddOption<string>("slew-calfit-file,l",
                       "",
                       "Slewing Calibration file or directory name");
  da.AddOption<string>("showerPDF,P",
                       "",
                       "Shower PDF file for LH fits. "
                       "Default in $CONFIG_HAWC");
  da.AddOption<string>("resi-calfit-file,t",
                       "",
                       "Time Residual Calibration file or directory name");
  da.AddOption<string>("zenith-alignment-file,z",
                       "",
                       "Zenith alignment file name");
  da.AddFlag("NoLinearInterpolation",
             "Do not use linear interpolation method for charge calibration");
  
  // Input parameters specific to MC
  OptionGroup& so = cl.AddOptionGroup("Simulation Configs");
  so.AddOption<int>("seed",-1,
                    "Seed for the random number generator.");
  so.AddOption<string>("config-file", "",
                       "Detector XML configuration file name");
  so.AddFlag("daqsim","Use if you want to use daqsim instead of "
             "SimPMTModeler. If enabled, the mc-parameter-file option "
             "is supplied to daqsim for configuration");
  so.AddOption<string>("mc-parameter-file,m",
                       "",
                       "XML config for SimPMTModeler. (or DAQSim ini config) "
                       "Default in $CONFIG_HAWC");
  so.AddOption<string>("peCurves,p",
                       "",
                       "PE smearing file for SimPMTModeler. "
                       "Default in $CONFIG_HAWC");
  so.AddOption<string>("hitDropProb",
                       "",
                       "Hit dropping probability file for SimPMTModeler. "
                       "Default in $CONFIG_HAWC");
  so.AddOption<string>("broadPulse8inch",
                       "",
                       "Pulse broadening file for 8 inch PMTs");
  so.AddOption<string>("broadPulse10inch",
                       "",
                       "Pulse broadening file for 10 inch PMTs");
  so.AddOption<string>("broadPulse8inchCharge",
                       "",
                       "Charge broadening file for 8 inch PMTs");
  so.AddOption<string>("broadPulse10inchCharge",
                       "",
                       "Charge broadening file for 10 inch PMTs");
  
  so.AddOption<string>("hitMaxChargeDrop",
                       "",
                       "Hit dropping due to max charge file for SimPMTModeler. "
                       "Default in $CONFIG_HAWC");
  so.AddFlag("no-hits", "No hit data will be written to full output file.");
  
  // Parameters for time residuals
  OptionGroup& tr = cl.AddOptionGroup("Time Residual Calculation");
  tr.AddFlag("timeresiduals", "Compute time residuals.");
  tr.AddOption<string>("timeresidualsoutput", "residuals.xcd",
                       "Time residuals output XCD file name.");
  tr.AddOption<int>("timeresidualsmaxevents", 200000,
                    "Max events used for time residuals.");
  tr.AddOption<int>("timeresidualsminhits", 1000,
                    "Min hits per channel for time residuals.");
  tr.AddOption<string>("timeresidualshistfilename", "tresiduals.root",
                       "Output root file for time residuals histograms. "
                       "Leave empty to not save them.");
  tr.AddOption<int>("timeresidualsminimumnchannels", 200,
                    "Minimum number of channels to compute time residuals.");
  tr.AddOption<float>("timeresidualsminzenith", -1e3,
                      "Minimum zenith angle to compute time residuals (deg).");
  tr.AddOption<float>("timeresidualsmaxzenith", 30.,
                      "Maximum zenith angle to compute time residuals (deg).");
  tr.AddOption<float>("timeresidualsazimuth", 0.,
                      "Reference azimuth angle to compute time residuals (deg).");
  tr.AddOption<float>("timeresidualsdeltaazimuth", 1e3,
                      "Maximum delta azimuth from reference azimuth to compute time residuals (deg).");
  tr.AddFlag("tResi-PE", "Do tResi-PE calibration");
  tr.AddOption<string>("tResiPE_CalFile","",
                       "The Time Residual vs PE Calibration File");
  
  // Parameters for charge residuals
  OptionGroup& qr = cl.AddOptionGroup("Charge Residual Calculation");
  qr.AddFlag("chargeresiduals", "Compute charge residuals.");
  qr.AddOption<int>("chargeresidualsmaxevents", -1,
                    "Max events used for charge residuals.");
  qr.AddOption<int>("chargeresidualsminhits", 1000,
                    "Min hits per channel for charge residuals");
  qr.AddOption<string>("chargeresidualshistfilename", "qresiduals.root",
                       "Output root file for charge residuals histograms. "
                       "Leave empty to not save them.");
  qr.AddOption<int>("chargeresidualsminimumnchannels", 200,
                    "Minimum number of channels to compute charge residuals.");
  qr.AddOption<double>("chargeresidualsminimumncharge", 10.,
                       "Minimum number of channels to compute charge residuals.");
  
  //Parameters for Likelihood Lateral Distribution Fit method
  OptionGroup& lhldf = cl.AddOptionGroup("Likelihood Lateral Distribution Fit");
  
  lhldf.AddOption<string>("prob-file-xmax-energy","","path to prob file for xmax and energy");
  lhldf.AddOption<string>("prob-file-main-array","","path to prob file for main array");
  lhldf.AddOption<string>("prob-file-outriggers","","path to prob file for outriggers");
  lhldf.AddOption<string>("guess-file-nHit-r-energy","","path to file for energy guess as a function of r and nHit ");
  lhldf.AddOption<string>("file-likelihood-vs-nHitTanks","","path to file for likelihood vs nHitTanks ");
  lhldf.AddFlag(
                "disable-lh-lat-dist-fit",
                "Do not perform the LH lateral distribution fit.");
  lhldf.AddFlag(
                "fit-Xmax",
                "Fit Xmax, otherwise by default use Xmax and energy correlation from MC instead of fitting Xmax");
  lhldf.AddFlag(
                "draw-lh-surface",
                "Draw likelihood surface and LDF, PDF of evnets. Must be used with --display NOT ALONE");
  
  //Parameters with respect to outriggers options (oro)
  OptionGroup& oro = cl.AddOptionGroup("Options for running with outiggers");
  oro.AddFlag("disable-or","Disables the automatic reading of outrigger files based on the trig_ file names");
  
  oro.AddOption<string>("or-diag-file","","output root file for or - diagnostics");
  oro.AddOption<string>("or-tres-file","","time residual file for outriggers");
  oro.AddFlag("disable-tres-or","Disable timing residuals for outriggers (needed when generating new ones)");
  oro.AddOption<bool>("hit-reco-wave",true,"Reconstruct hits from the FADC waveforms, it will overwrite the online found hits.");
  oro.AddOption<string>("wav-cal-path","","If specified, use this path for the calibration data instead of CONFIG_HAWC");
  
  
  // Parsing command lines
  if (!cl.ParseCommandLine(argc, argv)) {
    return 1;
  }
  
  // Check if even/odd fit selected, then extended output is selected
  if (cl.HasFlag("evenodd") && !cl.HasFlag("extended"))
    log_fatal("If using --evenodd, you must also specify --extended to output "
              "the fits.");
  
  // Check first input file for the file type to determine if we're dealing
  // with simulation or data
  bool isMC = false;
  bool hasOR = false;
  bool isTrigMC = false;
  FileList inputs = cl.GetArgument<FileList>("input");
  aerie_io::AERIEFileType ftype = aerie_io::GetFileType(inputs.front());
  switch (ftype) {
    case aerie_io::XCDF_TRIG:
      isMC = false;
      log_info("Reading in HAWC data files for input.");
      break;
    case aerie_io::XCDF_HAWCSIM:
      isMC = true;
      log_info("Reading in HAWCSIM simulated data files for input.");
      break;
    case aerie_io::XCDF_TRIG_OR:
      isMC = false;
      hasOR = true;
      log_info("Reading in HAWC TRIG_OR data files for input.");
      break;
    case aerie_io::XCDF_OR:
      isMC = false;
      hasOR = true;
      log_info("Reading in HAWC OR data files for input.");
      break;
    case aerie_io::XCDF_TRIG_MC:
      isTrigMC = true;
      log_info("Reading in HAWC TRIG_MC data files for input");
      break;
    default:
      log_error("File type " << ftype << " not supported.");
      return 1;
  }
  
  bool readOutrigger = !(cl.HasFlag("disable-or"));
  FileList inputs_or = inputs;
  unsigned int orFileCount = 0;
  if (ftype == aerie_io::XCDF_TRIG && readOutrigger) {
    for (unsigned int i = 0; i < inputs_or.size(); i++) {
      std::string::size_type pos = inputs_or[i].find("trig_");
      if (pos == std::string::npos) {
        //        log_fatal( "Cannot find trig_run in file name " <<inputs_or[i] );
        log_warn("Will not derive outrigger data names, could not find \"trig_\" part of file name");
        break;
      }
      inputs_or[i].replace(pos, 5, "or_coin_");
      //check if file exists
      if (FILE *file = fopen(inputs_or[i].c_str(), "r")) {
        fclose(file);
        orFileCount++;
      }
      log_info("derived or_coin file-name:" << inputs_or[i]);
    }  
    if (orFileCount == 0) {
      log_info("No outrigger files found to read, will run without");
      readOutrigger = false;
    } else {
      if (orFileCount != inputs.size()) {
        log_fatal(" not same number of OR / main array files found, OR: "
                  << orFileCount << ", Main: " << inputs.size());
      } else {
        log_info ("Reading outrigger data ");
        hasOR = true;
      }
    }
  }
  
  
  
  // Define the prefix for the TriggeredInputSelector
  string trig_prefix = "trig";
  if (isTrigMC) {
    trig_prefix = "trigMC";
  }
  
  bool doDiagnostics = false;
  string diagnosticsDirectory;
  string diagnosticsTag;
  if (cl.GetArgument<string>("diagnostics") != ""){
    doDiagnostics = true;
  }
  
  bool doMuons = false;
  string muonfinderFile;
  if(cl.GetArgument<string>("muonfinder-file") != ""){
    doMuons = true;
    muonfinderFile = cl.GetArgument<string>("muonfinder-file");
  }
  
  bool hasZenithAlign =
       cl.GetArgument<string>("zenith-alignment-file").compare("");

  string config_dir = cl.GetArgument<string>("config-dir");
  string config_file = cl.GetArgument<string>("config-file");
  string pe_curves = cl.GetArgument<string>("peCurves");
  string hitDropProbFile = cl.GetArgument<string>("hitDropProb");
  string hitMaxChargeDropFile = cl.GetArgument<string>("hitMaxChargeDrop");
  string mc_params = cl.GetArgument<string>("mc-parameter-file");
  string broadPulse8inch = cl.GetArgument<string>("broadPulse8inch");
  string broadPulse10inch = cl.GetArgument<string>("broadPulse10inch");
  string broadPulse8inchCharge = cl.GetArgument<string>("broadPulse8inchCharge");
  string broadPulse10inchCharge = cl.GetArgument<string>("broadPulse10inchCharge");
  
  log_info(""<< broadPulse8inch);
  
  log_info("Aerie version: " << SoftwareVersion(AERIE_VERSION_CODE));
  log_info("Aerie trunk revision: " << AERIE_REPOSITORY_REVISION);
  log_info("Aerie build type: " << AERIE_BUILD_TYPENAME);
  log_info("Config file ... " << config_file);
  log_info("ShowerPDF file ... " << cl.GetArgument<string>("showerPDF"));
  
  int random_seed = cl.GetArgument<int>("seed");
  

  HAWCNest nest;
  
  /////////////////////////////////////////
  // Defining the analysis cuts          //
  // (early but needed for data in       //
  // edge refiner)                       //
  //////////////////////////////////

  vector<string> hitSelections;
  vector<string> hitSelectionsOR;
  
  //Propagation plane cut defined below other cuts
  
  // Basic cuts, keep only reasonable hits
  vector<string> stdCuts;
  std::string stdCutsName = "stdCuts";
  hitSelections.push_back(stdCutsName);

  vector<string> stdCutsOR;
  std::string stdCutsNameOR = "stdCutsOR";
  hitSelectionsOR.push_back(stdCutsNameOR);
  
  // Cut on allowed PE values
  nest.Service<PECut>("peCut")
    ("minPE", 0.)
    ("maxPE", 10000.);
  stdCuts.push_back("peCut");

  if (!cl.HasFlag("no-remove-max-charge")) {
    // Filter out hits at larger than the maximum calibration
    nest.Service<MaxChargeCut>("maxChargeCut");
    stdCuts.push_back("maxChargeCut");
  }
  
  // Cut around trigger time
  // add time due to how sim-calibrator works
  nest.Service<TriggerTimeCut>("triggerTimeCut")
    ("minTime", -150*ns )
    ("maxTime",  400*ns );
  stdCuts.push_back("triggerTimeCut");

  // Cut on small hits below 1 p.e.
  nest.Service<TOTCut>("totCut")
    ("minLoTOT", 500);
  stdCuts.push_back("totCut");

  // Cut to remove hits flagged as "bad"
  nest.Service<BadHitCut>("badHitCut");
  stdCuts.push_back("badHitCut");

  // Cut to remove hits flagged as PE promotions
  nest.Service<PEPromHitCut>("pePromHitCut");
  stdCuts.push_back("pePromHitCut");

  // Reco cuts are the same as standard cuts unless MPF is used
  vector<string> recoCuts;
  recoCuts.push_back(stdCutsName);
  std::string recoCutsName = "recoCuts";

  if (!cl.HasFlag("no-splitter")) {
    recoCuts.push_back("MPFMainPlaneCut"); // Defined below
  }

  // In case we want to exclude channels from the track fit
  vector<string> trackFitCuts;
  std::string trackFitCutsName = "stdCutsTrackFit";
  hitSelections.push_back(trackFitCutsName);
  trackFitCuts=recoCuts;

  // Remove low charges from track fit if requested
  if (cl.GetArgument<double>("track-fit-min-npe")>0) {
    // Cut on more than n PEs
    nest.Service<PECut>("npeCut")
      ("minPE", cl.GetArgument<double>("track-fit-min-npe"))
      ("maxPE", 10000.);
    trackFitCuts.push_back("npeCut");
  }
  // Single channel?
  if (cl.GetArgument<int>("exclude-channel-from-track-fit")>=0) {
    nest.Service<ExcludeChannel>("excludeOneChannel")
      ("channel", cl.GetArgument<int>("exclude-channel-from-track-fit"));
    trackFitCuts.push_back("excludeOneChannel");
  }

  // Multiple channels?
  string excludePMTFileName = cl.GetArgument<string>("exclude-channels-from-track-fit-file");
  if (excludePMTFileName.compare("")) {
    ifstream excludePMTFile(excludePMTFileName.c_str());
    if (excludePMTFile.fail())
      log_fatal("Can not read file "<< excludePMTFileName);
    string PMTId;
    while (excludePMTFile>>PMTId){
      string excludeName = "excludeChannel-" + PMTId;
      nest.Service<ExcludeChannel>(excludeName)
        ("channel", atoi(PMTId.c_str()));
      trackFitCuts.push_back(excludeName);
    }
  }

  nest.Service<HitSelectionAND>(trackFitCutsName)
    ("selections", trackFitCuts);

  // Hits compatible with the propagation plane
  vector<string> ppCuts;
  std::string ppCutsName = "cxpeCuts";
  hitSelections.push_back(ppCutsName);
  ppCuts=recoCuts;
  ppCuts.push_back("propagationPlaneCut"); // Defined below

  nest.Service<HitSelectionAND>(ppCutsName)
    ("selections", ppCuts);

  /////////////////////////////////////////
  // Reading input files and calibrating //
  // Different for data and simulation   //
  /////////////////////////////////////////
  // main loop service
  vector<string> chain;
  std::string eventName = "event";
 
  nest.Service<SelectOutriggerChannelsOR>("stdCutsOR")
    ("eventName", eventName);
  stdCutsOR.push_back("stdCutsOR");

  //Input/output services
  string readerName = "reader";
  string serializerName = "serializer";
  
  // Output extended results vector
  vector<string> resultsExtended;
  
  // Output summary results vector
  vector<string> resultsSummary;
  
  nest.Service<StdRNGService>("random")
    ("seed", random_seed);
  
  // Astro transformation service
  nest.Service<StdAstroService>("astroService");
  
  // Check validity of cammand line parameter ranges -
  // Here we check if command line params have sensible ranges and if 
  // not issue a warning or a fatal error as appropriate.
  int gpcri = cl.GetArgument<int>("gauss-plane-curv-recompute-iter");
  if (gpcri<0) {
    log_fatal("Invalid value for 'gauss-plane-curv-recompute-iter'. Value must be >=0. Called with value: "<<gpcri);
  }

  // Curvature Correction
  int curvatureModel = cl.GetArgument<int>("curvature-model");
  if (curvatureModel == 0){
    nest.Service<HawcCurvatureService>("curvatureCorrection");
  }
  if (curvatureModel == 1){
    nest.Service<TweakedCurvatureService>("curvatureCorrection")
      ("lq1",cl.GetArgument<double>("tweaked-curvature-lq1"))
      ("lq2",cl.GetArgument<double>("tweaked-curvature-lq2"))
      ("r1",cl.GetArgument<double>("tweaked-curvature-r1"))
      ("r2",cl.GetArgument<double>("tweaked-curvature-r2"))
      ("rlq",cl.GetArgument<double>("tweaked-curvature-rlq"));
  }
  if (curvatureModel == 2){
    nest.Service<CrabCurvatureService>("curvatureCorrection");
  }
  if (curvatureModel == 3){
    nest.Service<ConfigurableCurvatureService>("curvatureCorrection")
      ("lq1",cl.GetArgument<double>("tweaked-curvature-lq1"))
      ("lq2",cl.GetArgument<double>("tweaked-curvature-lq2"))
      ("r1",cl.GetArgument<double>("tweaked-curvature-r1"))
      ("r2",cl.GetArgument<double>("tweaked-curvature-r2"))
      ("rlq",cl.GetArgument<double>("tweaked-curvature-rlq"))
      ("curveFile",cl.GetArgument<string>("Configurable-Curvature-File"));
  }
  if (curvatureModel < 0 or curvatureModel > 3){
    log_fatal("Invalid curvature-model setting: "<<curvatureModel);
  }
  
  if (!isMC) { // reading data files
    log_info("Assuming data");
    
    // Detector
    nest.Service<ConfigDirDetectorService>("det")
      ("addOR",hasOR)
      ("configDir", config_dir);
    
    // Calibration
    std::string interMethod = "LinearInterpolation";
    if (cl.HasFlag("NoLinearInterpolation")) interMethod = "MilagroFit";
    
    nest.Service<ChargeCalibrationService>("chargeCalibrator")
      ("chargeCalibXCDF",cl.GetArgument<string>("charge-calfit-file"))
      ("interpolationMethod",interMethod)
      ("configDir", config_dir);
    
    string tResiPE_file=cl.GetArgument<string>("tResiPE_CalFile");

    if (tResiPE_file !="" and !cl.HasFlag("tResi-PE"))
      log_fatal("To enable tResi-PE calibration, you must add the flag");

    nest.Service<TimingCalibrationService>("timingCalibrator")
      ("slewingXCDF",cl.GetArgument<string>("slew-calfit-file"))
      ("timeResidualXCDF",cl.GetArgument<string>("resi-calfit-file"))
      ("interpolationMethod","6ParFit")
      ("configDir", config_dir)
      ("tResi", cl.HasFlag("tResi-PE"))
      ("WidePulse", tResiPE_file);
    
    // Binary Event Source
    if (!readOutrigger) {
      nest.Service<Reader>(readerName)
        ("dataName",  "TriggeredData")
        ("inputFiles", inputs);
    } else {
      nest.Service<MergerORReader>(readerName)
        ("dataNameOR",  "ORData")
        ("inputFilesOR", inputs_or)
        ("dataNameMain",  "TriggeredData")
        ("inputFilesMain", inputs);
    }
    // Triggered Input Manager
    nest.Service<TriggeredInputSelector>("selector")
      ("dataMapName",  "TriggeredData")
      ("eventName",  "rawEvent")
      ("detectorService", "det")
      ("prefix", trig_prefix);
    chain.push_back("selector");

    // Filter based on number of hits
    nest.Service<ThresholdFilter>("thresholdFilter")
      ("eventName",   "rawEvent")
      ("resultName",  "thresholdResult")
      ("threshold",   cl.GetArgument<int>("threshold"))
      ("wThreshold",  cl.GetArgument<int>("wthreshold"));
    chain.push_back("thresholdFilter");
    
    // Filter calibration events
    nest.Service<CalibrationEventFilter>("calFilter")
      ("eventName",  "rawEvent")
      ("passPhysicsEvents", true)
      ("passCalibrationEvents", false);
    chain.push_back("calFilter");
    
    // Do edge refining
    nest.Service<EdgeRefiner>("edgeRefiner")
      ("eventName",  "rawEvent")
      ("refinedEventName",  eventName)
      ("trigMinTime", (nest.GetParameter<double>("triggerTimeCut", "minTime")))
      ("trigMaxTime", (nest.GetParameter<double>("triggerTimeCut", "maxTime")));
    chain.push_back("edgeRefiner");
    
    // Calibration
    nest.Service<Calibrator>("calibrator")
      ("event",  eventName)
      ("chargeCalibrationService", "chargeCalibrator")
      ("timingCalibrationService", "timingCalibrator");
    chain.push_back("calibrator");
    
    if (ftype == aerie_io::XCDF_TRIG_OR || readOutrigger) {
      //Outrigger Input Manager
      string orDataName = "ORData";
      if (ftype == aerie_io::XCDF_TRIG_OR) {
        orDataName = "TriggeredData";
      }

      bool recalibrate = cl.GetArgument<bool>("hit-reco-wave");
      nest.Service<OutriggerInputSelector>("selectorOR")
        ("dataMapName", orDataName)
        ("eventName", eventName)
        ("recalibrate",recalibrate)
        ("prefix", "or")
        ("detectorService", "det");
      chain.push_back("selectorOR");
      
      if (recalibrate) {
        
        string sPath =cl.GetArgument<string>("wav-cal-path");
        //get string from config HAWC (default)
        if (sPath == "") {
          const char* const env = getenv("CONFIG_HAWC");
          if (!env) {
            log_fatal("CONFIG_HAWC undefined.");
          }
          sPath = env;
          sPath += "/calibrations/outrigger/signal-reco/";
        }
        
        nest.Service<ORSignalReconstructor>("ORSignalReconstructor")
          ("eventName", eventName)
          ("configDir",sPath)
          ("detectorService", "det");
        chain.push_back("ORSignalReconstructor");
      }
      
      if (!cl.HasFlag("disable-tres-or")) {
        string orTresFile = cl.GetArgument<string>("or-tres-file");
        string configDir ="";
        //no tres file specified by user
        if ( orTresFile == "") {
          //user specified config dir
          if (config_dir!="") {
            configDir=config_dir;
          } else {//use CONFIG_HAWC (Default)
            const char* const env = getenv("CONFIG_HAWC");
            if (!env) {
              log_fatal("CONFIG_HAWC undefined.");
            }
            configDir = env;
          }
          configDir +=
          "/calibrations/outrigger/timing-residuals/";
        }
      
        nest.Service<ORTimingCalibrator>("ORTimingCalib")
        ("eventName", eventName)
        ("configDir", configDir)
        ("timeResFile",orTresFile);
        chain.push_back("ORTimingCalib");
      }
    }
    
  } else { // reading MC files
    log_info("Assuming simulated data");
    // Give user some useful info if running simulation
    log_info("Using mc parameter file: " << mc_params);
    // No mc-parameter file given, but trying to reconstruct simulation, this will fail, stop it now
    // 0 if they are equal
    if (mc_params.compare("")==0)
      {
	log_fatal("You must provide use the -m/--mc-parameter-file flag when reconstructing simulation");
      }
    
    // Detector
    if (config_file.compare("")) {
      nest.Service<StdDetectorService>("det")
        ("configFile", config_file)
        ("validateConfigs", false);
    }
    else{
      log_fatal("--config-file needed for simulated data.");
    }
    
    log_info("added detector service");
    
    // XCDF Reader
    nest.Service<Reader>(readerName)
      ("dataName", "hawcsimData")
      ("inputFiles", inputs);
    
    // MC Event Creator
    nest.Service<HAWCSimInputSelector>("eventSource")
      ("dataMapName","hawcsimData")
      ("eventName","SimEvent")
      ("detectorService", "det")
      ("loadWaterHits",false);
    
    chain.push_back("eventSource");
    resultsExtended.push_back("SimEvent");
    
    log_info("added source");

    if(!cl.HasFlag("daqsim")){
      // Calibrate the MC Event
      nest.Service<SimulatedPMTModeler>("cal")
        ("eventName", eventName)
        ("simEventName","SimEvent")
        ("detectorService",  "det")
        ("rngService", "random")
        ("configDir", config_dir)
        ("peCurves", pe_curves)           // optional override of CONFIG_HAWC
        ("hitDropProb", hitDropProbFile)  // optional override of CONFIG_HAWC
        ("hitMaxChargeDrop", hitMaxChargeDropFile)
        ("mcParameterFile", mc_params)    // optional override of CONFIG_HAWC
        ("broadPulse8inch", broadPulse8inch)
        ("broadPulse10inch", broadPulse10inch)
        ("broadPulse8inchCharge", broadPulse8inchCharge)
        ("broadPulse10inchCharge", broadPulse10inchCharge);
    }
    else{
      NestIniConfig(nest,mc_params);
    }
    chain.push_back("cal");
    log_info("added calibration");

    // Filter based on number of hits
    nest.Service<ThresholdFilter>("thresholdFilter")
      ("eventName",   eventName)
      ("resultName",  "thresholdResult")
      ("threshold",   cl.GetArgument<int>("threshold"))
      ("wThreshold",  cl.GetArgument<int>("wthreshold"));
      chain.push_back("thresholdFilter");
    
      //reconstruct signals from waveforms
      string sPath =cl.GetArgument<string>("wav-cal-path");
      //get string from config HAWC (default)
      if (sPath == "") {
        const char* const env = getenv("CONFIG_HAWC");
        if (!env) {
          log_fatal("CONFIG_HAWC undefined.");
        }
        sPath = env;
        sPath += "/calibrations/outrigger/signal-reco/";
      }
      
      nest.Service<ORSignalReconstructor>("ORSignalReconstructor")
        ("eventName", eventName)
        ("configDir",sPath)
        ("detectorService", "det");
      chain.push_back("ORSignalReconstructor");
    
  }

  /////////////////////////////////////////
  // Physics modules                     //
  /////////////////////////////////////////

  // Standard basic cuts
  nest.Service<HitSelectionAND>(stdCutsName)
    ("detectorServiceName",  "det")
    ("eventName",  eventName)
    ("selections", stdCuts)
    ("useTriggerTimeCut", true)
    ("trigMinTime",       (nest.GetParameter<double>("triggerTimeCut", "minTime")))
    ("trigMaxTime",       (nest.GetParameter<double>("triggerTimeCut", "maxTime")));

  // Standard basics cuts plus MPF, if used
  nest.Service<HitSelectionAND>(recoCutsName)
    ("detectorServiceName",  "det")
    ("eventName",  eventName)
    ("selections", recoCuts)
    ("useTriggerTimeCut", true)
    ("trigMinTime",       (nest.GetParameter<double>("triggerTimeCut", "minTime")))
    ("trigMaxTime",       (nest.GetParameter<double>("triggerTimeCut", "maxTime")));

  // the multi-plane fit module
  string mpfResult = "";
    
  if (!cl.HasFlag("no-splitter")) {

    mpfResult = "MPFplanes";
    
    // MPF splitter: computes and stores the main plane
    nest.Service<MPFEventSplitter>("MPFEventSplitter")
      ("event",  eventName)
      ("detectorService",  "det")
      ("resultName", mpfResult)
      ("planeWidthMaxNs", 30)
      ("hitSelection", stdCutsName)
      ("minHitsPerPlane", 5)
      ("numPlanes", 3);

    chain.push_back("MPFEventSplitter");

    // Hit selector: lets the other modules know if a hit
    // belongs to the main plane
    nest.Service<MPFMainPlaneCut>("MPFMainPlaneCut")
      ("mpfResult", mpfResult);

  }  
  
  // Charge Filter
  // Determine frac of charge in Fiducial Area
  string fiducialChargeName = "fiduCharge";
  nest.Service<FiducialCharge>("fiducialCharge")
    ("detectorService",  "det")
    ("detectorScaleFactor", -1)
    ("event",  eventName)
    ("resultName",  fiducialChargeName);
  chain.push_back("fiducialCharge");
  
  // Reconstruction
  // First-Guess Core Reconstruction
  string coreGuessName = "coreGuess";
  string hitSelectionForCoreName = recoCutsName;
  nest.Service<COMCoreFit>("comCoreFit")
    ("detectorService",  "det")
    ("event",  eventName)
    ("resultName",  coreGuessName)
    ("hitSelection", hitSelectionForCoreName);
  chain.push_back("comCoreFit");
  
  if (!cl.HasFlag("no-core-fit-in-plane")) {
    // Do an extra step of core/plane fit, to define the hit selection
    // compatible with the plane for the next core fit.
    
    // SFCF without plane cut
    string coreSecondGuessName = "coreSecondGuess";
    nest.Service<SFCF>("coreFitGuess")
      ("detectorService", "det")
      ("event", eventName)
      ("resultName", coreSecondGuessName)
      ("seed", coreGuessName)
      ("strategy", cl.GetArgument<int>("sfcf-strategy"))
      ("maxIterations", 100)
      ("useOR", !cl.HasFlag("sfcf-disable-or"))
      ("useTrueCore", cl.HasFlag("sfcf-enable-true-core"))
      ("hitSelection", hitSelectionForCoreName)
      ("hitSelectionOutrigger", "stdCutsOR");
    chain.push_back("coreFitGuess");
    coreGuessName = coreSecondGuessName;
    
    // Plane fit using core fit without plane cut
    nest.Service<GaussPlaneFit>("planeFitGuess")
      ("detectorService",  "det")
      ("event",  eventName)
      ("resultName",  "planeGuess")
      ("coreFit",  coreGuessName)
      ("maxIterations",  5)
      ("peCutThresholds", "0.5,0.5,0.5,0.5,0.5")
      ("chiThresholds", "5.,4.,3.,2.")
      ("hitSelection", trackFitCutsName)
      ("curvatureModel","curvatureCorrection")
      ("maxShowerPlaneIter",cl.GetArgument<int>("gauss-plane-curv-recompute-iter"));
    chain.push_back("planeFitGuess");
    
    // Hit selection at +/- 50 ns of PF0
    nest.Service<PropagationPlaneCut>("propagationPlaneCutGuess")
      ("minTimeResidual", -50*ns)
      ("maxTimeResidual", 50*ns)
      ("angleFitName", "planeGuess") // right now we only use plane cut for cxpe
      ("detectorServiceName", "det")
      ("eventName", eventName)
      ("coreFitName", coreGuessName)
      ("curvatureCorrectTimeResiduals", ! cl.HasFlag("no-correct-curv-tr"))
      ("curvatureModel","curvatureCorrection");

    // Reco hit selection cuts
    vector<string> stdCutsPlane = recoCuts;
    stdCutsPlane.push_back("propagationPlaneCutGuess");
    nest.Service<HitSelectionAND>("stdCutsPlane")
      ("selections", stdCutsPlane);
    hitSelectionForCoreName = "stdCutsPlane";

  }
  resultsExtended.push_back(coreGuessName);
  
  //////////Check if likelihood lateral distribution fit is enabled/////////////
  string coreName = "coreFit"; ////made these two variables globle because of
  string pmtChi2Name = ""; ///// the if statement to check the lhlatdist flag
  
  if(cl.HasFlag("disable-lh-lat-dist-fit")){
    
    // Refined core reconstruction
    string coreType = cl.GetArgument<string>("core-fit");
    
    if (coreType == "sfcf") {
      // SFCF Reconstruction
      nest.Service<SFCF>("SFCF")
        ("detectorService", "det")
        ("event", eventName)
        ("resultName", coreName)
        ("seed", coreGuessName)
        ("strategy", cl.GetArgument<int>("sfcf-strategy"))
        ("maxIterations", 100)
        ("useOR", !cl.HasFlag("sfcf-disable-or"))
        ("useTrueCore", cl.HasFlag("sfcf-enable-true-core"))
        ("hitSelection", hitSelectionForCoreName)
        ("hitSelectionOutrigger", "stdCutsOR")
        ("storeFitData", ( cl.HasFlag("display")  )); // if display, save all fit data
      chain.push_back("SFCF");
    } else if (coreType == "gauss") {
      // PMT Gauss Core Reconstruction
      nest.Service<GaussCoreFit>("gaussCoreFit")
        ("detectorService",  "det")
        ("event",  eventName)
        ("resultName",  coreName)
        ("seed",  coreGuessName)
        ("sigmaLimit", 150*HAWCUnits::meter)
        ("maxIterations",  100)
        ("hitSelection", hitSelectionForCoreName);
      chain.push_back("gaussCoreFit");
    } else if (coreType == "gauss-tank") {
      // Tank Gauss Core Reconstruction
      nest.Service<TankGaussCoreFit>("gaussCoreFit")
        ("detectorService",  "det")
        ("event",  eventName)
        ("resultName",  coreName)
        ("seed",  coreGuessName)
        ("sigmaLimit", 150*HAWCUnits::meter)
        ("maxIterations",  100)
        ("hitSelection", hitSelectionForCoreName);
      chain.push_back("gaussCoreFit");
    } else if (coreType == "nkg") {
      // NKG Core Reconstruction
      nest.Service<NKGCoreFit>("NKGCoreFit")
        ("detectorService",  "det")
        ("event",  eventName)
        ("resultName",  coreName)
        ("seed",  coreGuessName)
        ("sigmaLimit", 150*HAWCUnits::meter)
        ("maxIterations",  100)
        ("hitSelection", hitSelectionForCoreName)
        ("useMINUIT", cl.HasFlag("lnkg"))
        ("useLogErrs", cl.HasFlag("lnkg"));
      chain.push_back("NKGCoreFit");
    } else if (coreType == "nkg-tank") {
      // Tank NKG/R Reconstruction
      nest.Service<TankLMCoreFit>("TankLMFit")
        ("detectorService", "det")
        ("event", eventName)
        ("resultName", coreName)
        ("seed", coreGuessName)
        ("sigmaLimit", 150*HAWCUnits::meter)
        ("maxIterations", 100)
        ("function",2)
        ("initialS",0.9)
        ("radialCut",10.0)
        ("smearValue",0.1225)
        ("hitSelection", hitSelectionForCoreName);
      chain.push_back("TankLMFit");
    } else {
      log_fatal("Invalid core type (" << coreType << "), please check the "
                "--core-fit argument");
    }
    resultsExtended.push_back(coreName);
    
    if (cl.HasFlag("pmtChi2")) {
      pmtChi2Name = "PMTChi2";
      if (coreType == "nkg-tank") {
        nest.Service<PMTChi2>("Chi2NKGR")
          ("detectorService", "det")
          ("event",eventName)
          ("hitSelection", recoCutsName)
          ("coreFit",coreName)
          ("resultName",  pmtChi2Name)
          ("function","nkg/r")
          ("radius",10*HAWCUnits::meter);
        chain.push_back("Chi2NKGR");
      } else if (coreType == "sfcf") {
        nest.Service<PMTChi2>("Chi2SFCF")
          ("detectorService", "det")
          ("event",eventName)
          ("hitSelection", recoCutsName)
          ("coreFit",coreName)
          ("resultName",  pmtChi2Name)
          ("function","sfcf");
        chain.push_back("Chi2SFCF");
      } else if (coreType == "nkg") {
        nest.Service<PMTChi2>("Chi2NKG")
          ("detectorService", "det")
          ("event",eventName)
          ("hitSelection", recoCutsName)
          ("coreFit",coreName)
          ("resultName",  pmtChi2Name)
          ("function","nkg");
        chain.push_back("Chi2NKG");
      } else {
        log_warn("Core type " << coreType << "not supported by --pmtChi2, "
                 "turning off PMTChi2.");
        pmtChi2Name = "";
      }
    }
  }
  
  //////////this will do Likelihood Lateral distribution fit////////////////
  
  if (!cl.HasFlag("disable-lh-lat-dist-fit")){
    
    string probFileMainArray = cl.GetArgument<string>("prob-file-main-array");
    string probFileOutriggers = cl.GetArgument<string>("prob-file-outriggers");
    string probFileXmaxEnergy = cl.GetArgument<string>("prob-file-xmax-energy");
    string guessFilenHitrEnergy = cl.GetArgument<string>("guess-file-nHit-r-energy");
    string lookupFileLikelihoodVsNhitTanks = cl.GetArgument<string>("file-likelihood-vs-nHitTanks");
    
    if(probFileMainArray.empty()){
      const char* const env = getenv("CONFIG_HAWC");
      if (!env)
        log_fatal("CONFIG_HAWC undefined.");
      
      string configHAWCDir = env;
      string LHLatDistFitSubDir = "/reconstruction/LHLatDistFit";
      
      probFileMainArray = configHAWCDir + LHLatDistFitSubDir + "/logprob_main_array.root";
    }
    
    if(probFileOutriggers.empty()){
      const char* const env = getenv("CONFIG_HAWC");
      if (!env)
        log_fatal("CONFIG_HAWC undefined.");
      
      string configHAWCDir = env;
      string LHLatDistFitSubDir = "/reconstruction/LHLatDistFit";
      
      probFileOutriggers = configHAWCDir + LHLatDistFitSubDir + "/logprob_outriggers.root";
    }
    
    if(probFileXmaxEnergy.empty()){
      const char* const env = getenv("CONFIG_HAWC");
      if (!env)
        log_fatal("CONFIG_HAWC undefined.");
      
      string configHAWCDir = env;
      string LHLatDistFitSubDir = "/reconstruction/LHLatDistFit";
      
      probFileXmaxEnergy = configHAWCDir + LHLatDistFitSubDir + "/logprob_xmax_energy.root";
    }
    
    if(guessFilenHitrEnergy.empty()){
      const char* const env = getenv("CONFIG_HAWC");
      if (!env)
        log_fatal("CONFIG_HAWC undefined.");
      
      string configHAWCDir = env;
      string LHLatDistFitSubDir = "/reconstruction/LHLatDistFit";
      
      guessFilenHitrEnergy = configHAWCDir + LHLatDistFitSubDir + "/Distribution_nHit_R_Energy.root";
    }
    
    if(lookupFileLikelihoodVsNhitTanks.empty()){
      const char* const env = getenv("CONFIG_HAWC");
      if (!env)
        log_fatal("CONFIG_HAWC undefined.");
      
      string configHAWCDir = env;
      string LHLatDistFitSubDir = "/reconstruction/LHLatDistFit";
      
      lookupFileLikelihoodVsNhitTanks = configHAWCDir + LHLatDistFitSubDir + "/likelihood_vs_nHitTanks.root";
    }

    //// only doing one iteration with sfcf as the first guess 
    //// to speed things up 
    nest.Service<LHLatDistFit>("LHLatDistFit")
      ("detectorService",  "det")
      ("event", eventName)
      ("resultName",  coreName)
      ("seed",coreGuessName)
      ("angleFitName","planeGuess")
      ("hitSelection", hitSelectionForCoreName)
      ("hitSelectionOutrigger", "stdCutsOR")
      ("fitXmax", cl.HasFlag("fit-Xmax"))
      ("drawLHSurface", cl.HasFlag("draw-lh-surface"))
      ("display", cl.HasFlag("display"))
      ("probFileMainArray",probFileMainArray)
      ("probFileOutriggers",probFileOutriggers)
      ("probFileXmaxEnergy",probFileXmaxEnergy)
      ("guessFileNhitREnergy",guessFilenHitrEnergy)
      ("lookupFileLikelihoodVsNhitTanks",lookupFileLikelihoodVsNhitTanks);
    chain.push_back("LHLatDistFit");
    resultsExtended.push_back(coreName);
    
  }
  ///////////////////////////////////////////////////////////////////////////////////
  
  // Note: the LH/PDF angle fits may not be run. In this case, force all
  // downstream modules to use the results of the plane fit.
  string angleName = "angleFit";
  string planeName = cl.HasFlag("enable-likelihood") ? "planeFit" : "angleFit";
  string cxpeAngle = planeName;
  
  // Definition of propagation plane cut
  nest.Service<PropagationPlaneCut>("propagationPlaneCut")
    ("minTimeResidual", -10*ns)
    ("maxTimeResidual", 10*ns)
    ("angleFitName", cxpeAngle) // right now we only use plane cut for cxpe
    ("detectorServiceName", "det")
    ("eventName", eventName)
    ("coreFitName", coreName)
    ("curvatureCorrectTimeResiduals", ! cl.HasFlag("no-correct-curv-tr"))
    ("curvatureModel","curvatureCorrection");


  // Milagro plane fit -- no PE cuts
  nest.Service<GaussPlaneFit>("gaussPlaneFit")
    ("detectorService",  "det")
    ("event", eventName)
    ("resultName",  planeName)
    ("coreFit",  coreName)
    ("maxIterations",  5)
    ("peCutThresholds", "0.5,0.5,0.5,0.5,0.5")
    ("chiThresholds", "5.,4.,3.,2.")
    ("hitSelection", trackFitCutsName)
    ("curvatureModel","curvatureCorrection")
    ("maxShowerPlaneIter",cl.GetArgument<int>("gauss-plane-curv-recompute-iter"));
  chain.push_back("gaussPlaneFit");
  resultsExtended.push_back(planeName);
  
  // Likelihood Angle Reconstruction
  if (cl.HasFlag("enable-likelihood")) {
    if (cl.HasFlag("pdffit")) {
      // New PDFAngleFit
      nest.Service<PDFAngleFit>("pdfAngleFit")
        ("detectorService",  "det")
        ("event", eventName)
        ("resultName",  angleName)
        ("coreFit",  coreName)
        ("seed", planeName)
        ("minimumHits",  5)
        ("showerPDFTable",  cl.GetArgument<string>("showerPDF"))
        ("hitSelection", trackFitCutsName);
      chain.push_back("pdfAngleFit");
    } else {
      // Old LHAngleFit
      nest.Service<LHAngleFit>("lhAngleFit")
        ("detectorService",  "det")
        ("event", eventName)
        ("resultName",  angleName)
        ("coreFit",  coreName)
        ("seed", planeName)
        ("minimumHits",  5)
        ("showerPDFTable",  cl.GetArgument<string>("showerPDF"))
        ("hitSelection", trackFitCutsName);
      chain.push_back("lhAngleFit");
    }
    resultsExtended.push_back(angleName);
  }
  // Nhit calculation
  std::string nhitname="nhit";
  nest.Service<NhitCalculator>("NhitCal")
    ("detectorService", "det")
    ("event", eventName)
    ("coreFitName", coreName)
    ("angleFitName", angleName)
    ("hitSelection", stdCutsName)
    ("resultName", nhitname)
    ("curvatureModel","curvatureCorrection");
  chain.push_back("NhitCal");



  // Compactness calculation
  
  if (cl.HasFlag("allCxPERadii")) {
    
    nest.Service<CompactnessCalculator>("compCalc20")
      ("detectorService",  "det")
      ("event",  eventName)
      ("radius",  20.*meter)
      ("resultName",  "compactness20")
      ("coreFit",  coreName)
      ("hitSelection", ppCutsName)
      ("angleFitName", cxpeAngle)
      ("curvatureModel","curvatureCorrection");
    chain.push_back("compCalc20");
    resultsExtended.push_back("compactness20");
    
    nest.Service<CompactnessCalculator>("compCalc30")
      ("detectorService",  "det")
      ("event",  eventName)
      ("radius",  30.*meter)
      ("resultName",  "compactness30")
      ("coreFit",  coreName)
      ("hitSelection", ppCutsName)
      ("angleFitName", cxpeAngle)
      ("curvatureModel","curvatureCorrection");
    chain.push_back("compCalc30");
    resultsExtended.push_back("compactness30");
    
  }
  //allways calculate 40m
  nest.Service<CompactnessCalculator>("compCalc40")
    ("detectorService",  "det")
    ("event",  eventName)
    ("radius",  40.*meter)
    ("resultName",  "compactness40")
    ("coreFit",  coreName)
    ("hitSelection", ppCutsName)
    ("angleFitName", cxpeAngle)
    ("curvatureModel","curvatureCorrection");
  chain.push_back("compCalc40");
  resultsExtended.push_back("compactness40");

  
  if (cl.HasFlag("allCxPERadii")) {
    nest.Service<CompactnessCalculator>("compCalc50")
      ("detectorService",  "det")
      ("event",  eventName)
      ("radius",  50.*meter)
      ("resultName",  "compactness50")
      ("coreFit",  coreName)
      ("hitSelection", ppCutsName)
      ("angleFitName", cxpeAngle)
      ("curvatureModel","curvatureCorrection");
    chain.push_back("compCalc50");
    resultsExtended.push_back("compactness50");
    
  }
  
  // Do PINC and pair compactness calculations if extended output is requested
  string pincnessName = "";
  
  if (cl.GetArgument<int>("pincness") != 0) {
    pincnessName = "pincness";
    nest.Service<PINCCalculator>("pincCalc")
      ("detectorService",  "det")
      ("event",  eventName)
      ("radius",  cl.GetArgument<double>("pincRadius")*meter)
      ("resultName",  pincnessName)
      ("coreFit",  coreName)
      ("angleFit", planeName)
      ("hitSelection", ppCutsName)
      ("version", cl.GetArgument<int>("pincness"));
    chain.push_back("pincCalc");
    resultsExtended.push_back("pincness");
  }
  
  // Features neural networks
  if (!cl.HasFlag("no-features-nn")) {
    nest.Service<FeaturesNN>("featuresnn")
      ("detectorService",  "det")
      ("event",  eventName)
      ("resultName",  "FeaturesNN")
      ("coreFit",  coreName)
      ("angleFit", planeName)
      ("hitSelection", recoCutsName);
      chain.push_back("featuresnn");
      resultsExtended.push_back("FeaturesNN");
  }
  
  string pairCompactnessName = "";
  
  if (cl.HasFlag("pairCompactness")) {
    pairCompactnessName = "pairCompactness";
    nest.Service<PairCompactnessCalculator>("pairCompactnessCalc")
      ("detectorService",  "det")
      ("event",  eventName)
      ("resultName",  pairCompactnessName)
      ("hitSelection", ppCutsName);
    chain.push_back("pairCompactnessCalc");
    resultsExtended.push_back("pairCompactness");
  }
  
  //Tank Gamma-Hadron Separation
  string OutputFile = cl.GetArgument<string>("GHoutput-file");
  string HistogramFileName = cl.GetArgument<string>("GHHistogram-file");

  string TGHSName = "TGHS";

  string PInputFile = cl.GetArgument<string>("Pinput-file");
  string GInputFile = cl.GetArgument<string>("Ginput-file");

  if (PInputFile.empty() || GInputFile.empty()) {
    const char* const env = getenv("CONFIG_HAWC");
    if (!env)
      log_fatal("Unable to do GH-Separation, CONFIG_HAWC undefined.");

    string configHAWCDir = env;
    string GHSSubDir = "/reconstruction/TankLightDist-Lookup";
    GInputFile = configHAWCDir + GHSSubDir + "/TankGammaHistogram.root";
    PInputFile = configHAWCDir + GHSSubDir + "/TankDataHistogram.root";
  }

  std::vector<double> veccuts;
  if (cl.HasFlag("PE-Cut")) {
    veccuts = cl.GetArgument< std::vector<double> >("PE-Cut");
  }

  nest.Service<TankGHSep>("tGHSep")
    ("LikelihoodFile",OutputFile)
    ("HistogramFile",HistogramFileName)
    ("ChargeCut",veccuts)
    ("PHistoFile",PInputFile)
    ("GHistoFile",GInputFile)
    ("resultName",  TGHSName)
    ("event", eventName)
    ("hitSelection",recoCutsName)
    ("detectorService",  "det")
    ("simEventName","SimEvent");
  chain.push_back("tGHSep");
  resultsExtended.push_back(TGHSName);
  log_info("added TankGHSep");
          
  const bool runLHEnergy = !cl.HasFlag("no-lh-energy");
  const string lheneRecName = runLHEnergy ? "lheneRec" : "";
  if (runLHEnergy) {
    string LHEConfigFile = cl.GetArgument<string>("table-config-file");
    string energyTable = "";
    if ( LHEConfigFile.empty() ) {
      if ( isMC ){
        const char* aerie_config_char = std::getenv("HAWC_CONFIG");
        string aerie_config = string(aerie_config_char);
        LHEConfigFile = aerie_config + "/LHEnergyConfig.xml";
      }
      else {
        const char* const env = getenv("CONFIG_HAWC");
        if (!env)
          log_fatal("CONFIG_HAWC undefined.");
        
        string configHAWCDir = env;
        
        string LHESubDir = "/reconstruction/energy/LHEnergyConfigDir";
        LHEConfigFile = configHAWCDir + LHESubDir + "/LHEnergyConfig.xml";
        energyTable = configHAWCDir + LHESubDir + "/energytables.fits";
      }
    }
    
    nest.Service<LHEnergyEstimator>("lhEneCalc")
      ("detectorService", "det")
      ("rngService", "random")
      ("event", eventName)
      ("coreFit", coreName)
      ("angleFit",planeName)
      ("resultName", lheneRecName)
      ("tableConfigFile", LHEConfigFile)
      ("energyTable",energyTable)
      ("hitSelection", recoCutsName);
    chain.push_back("lhEneCalc");
  }

  string fiducialCoreName = "fiduCore";
  nest.Service<FiducialCore>("fiducialCore")
    ("detectorService", "det")
    ("event", eventName)
    ("coreFit", coreName)
    ("resultName", fiducialCoreName);
  chain.push_back("fiducialCore");


  string latResult = "";

  if (!cl.HasFlag("no-axis"))
    {
      latResult = "latDistFit";

      nest.Service<LatDist>("LatDist")
        ("detectorService", "det")
        ("event", eventName)
        ("coreFit", coreName)
        ("angleFit", angleName)
        ("hitSelection", recoCutsName)
        ("configDir",  config_dir)
        ("fiducialCore", fiducialCoreName)
        ("resultName", latResult)
        ("verbose", 0)
        ("QErrVersion", cl.GetArgument<int>("LatDistQErr"))
        ("groundParameterVersion", cl.GetArgument<int>("GPVersion"));
      chain.push_back("LatDist");
    }

  // Neural-net energy calculation. Optionally enable LH energy estimator.
  const bool runNnEnergy = !cl.HasFlag("no-neural-net-energy");
  const string nnEneRecName = runNnEnergy ? "nnEneRec" : "";
  if (runNnEnergy) {
    nest.Service<NeuralNetEnergyEstimator>("nnEneCalc")
      ("detectorService",  "det")
      ("event",  eventName)
      ("coreFit",  coreName)
      ("angleFit",  angleName)
      ("latDistFit", latResult)
      ("configDir",  config_dir)
      ("nhit", nhitname)
      ("resultName",  nnEneRecName)
      ("hitSelection",  recoCutsName);
    chain.push_back("nnEneCalc");
  }

  if (doDiagnostics){
    nest.Service<TimeResidualDiagnostic>("tresdiagnostic")
      ("event",eventName)
      ("detectorService","det")
      ("coreFitName",coreName)
      ("angleFitName",planeName)
      ("hitSelection",stdCutsName)
      ("curvatureModel","curvatureCorrection")
      ("fHitThreshold",cl.GetArgument<double>("diagnostics-fhit"))
      ("outfile",cl.GetArgument<string>("diagnostics"));
    chain.push_back("tresdiagnostic");
  }
  if(doMuons){
    nest.Service<MuonFinder>("muonfinder")
      ("outfile",muonfinderFile)
      ("detectorService", "det")
      ("event", eventName)
      ("coreFit", coreName)
      ("angleFit", planeName)
      ("hitSelection", stdCutsName);
    chain.push_back("muonfinder");
  }
  
  
  // Even/Odd reconstruction
  
  if (cl.HasFlag("evenodd")) {
    // Cuts for even/odd PMTs
    vector<string> cutsEven;
    cutsEven=trackFitCuts;
    nest.Service<SelectEvenChannels>("selectEven");
    cutsEven.push_back("selectEven");
    nest.Service<HitSelectionAND>("stdCutsEven")
      ("selections", cutsEven);
    
    vector<string> cutsOdd;
    cutsOdd=trackFitCuts;
    nest.Service<SelectOddChannels>("selectOdd");
    cutsOdd.push_back("selectOdd");
    nest.Service<HitSelectionAND>("stdCutsOdd")
      ("selections", cutsOdd);

    hitSelections.push_back("stdCutsEven");
    hitSelections.push_back("stdCutsOdd");
    
    // Even/odd reconstructions -- no PE cuts
    nest.Service<GaussPlaneFit>("gaussPlaneFitEven")
      ("detectorService",  "det")
      ("event",  eventName)
      ("resultName",  "planeFitEven")
      ("coreFit",  coreName)
      ("maxIterations",  5)
      ("peCutThresholds", "0.5,0.5,0.5,0.5,0.5")
      ("chiThresholds", "5.,4.,3.,2.")
      ("hitSelection", "stdCutsEven")
      ("curvatureModel","curvatureCorrection")
      ("maxShowerPlaneIter",cl.GetArgument<int>("gauss-plane-curv-recompute-iter"));
    chain.push_back("gaussPlaneFitEven");
    resultsExtended.push_back("planeFitEven");
    
    nest.Service<GaussPlaneFit>("gaussPlaneFitOdd")
      ("detectorService",  "det")
      ("event",  eventName)
      ("resultName",  "planeFitOdd")
      ("coreFit",  coreName)
      ("maxIterations",  5)
      ("peCutThresholds", "0.5,0.5,0.5,0.5,0.5")
      ("chiThresholds", "5.,4.,3.,2.")
      ("hitSelection", "stdCutsOdd")
      ("curvatureModel","curvatureCorrection")
      ("maxShowerPlaneIter",cl.GetArgument<int>("gauss-plane-curv-recompute-iter"));
    chain.push_back("gaussPlaneFitOdd");
    resultsExtended.push_back("planeFitOdd");

    if (cl.HasFlag("enable-likelihood")) {

      string evenFitName = "";
      string oddFitName = "";
  
      if (cl.HasFlag("pdffit")) {
        // New PDFAngleFit
        evenFitName = "pdfFitEven";
        nest.Service<PDFAngleFit>("pdfAngleFitEven")
          ("detectorService",  "det")
          ("event",  eventName)
          ("resultName", evenFitName)
          ("coreFit",  coreName)
          ("seed", "planeFitEven")
          ("minimumHits",  5)
          ("showerPDFTable",  cl.GetArgument<string>("showerPDF"))
          ("hitSelection", "stdCutsEven");
        chain.push_back("pdfAngleFitEven");
        
        evenFitName = "pdfFitOdd";
        nest.Service<PDFAngleFit>("pdfAngleFitOdd")
          ("detectorService",  "det")
          ("event",  eventName)
          ("resultName", oddFitName)
          ("coreFit",  coreName)
          ("seed", "planeFitOdd")
          ("minimumHits",  5)
          ("showerPDFTable",  cl.GetArgument<string>("showerPDF"))
          ("hitSelection", "stdCutsOdd");
        chain.push_back("pdfAngleFitOdd");
      }
      else {
        // Old LHAngleFit
        evenFitName = "lhFitEven";
        nest.Service<LHAngleFit>("lhAngleFitEven")
          ("detectorService",  "det")
          ("event",  eventName)
          ("resultName",  evenFitName)
          ("coreFit",  coreName)
          ("seed", "planeFitEven")
          ("minimumHits",  5)
          ("verbose",  0)
          ("showerPDFTable",  cl.GetArgument<string>("showerPDF"))
          ("hitSelection", "stdCutsEven");
        chain.push_back("lhAngleFitEven");

        oddFitName = "lhFitOdd";
        nest.Service<LHAngleFit>("lhAngleFitOdd")
          ("detectorService",  "det")
          ("event",  eventName)
          ("resultName",  oddFitName)
          ("coreFit",  coreName)
          ("seed", "planeFitOdd")
          ("minimumHits",  5)
          ("verbose",  0)
          ("showerPDFTable",  cl.GetArgument<string>("showerPDF"))
          ("hitSelection", "stdCutsOdd");
        chain.push_back("lhAngleFitOdd");
      }
      resultsExtended.push_back(evenFitName);
      resultsExtended.push_back(oddFitName);
    }
  } // end of if (evenodd)

  if (cl.HasFlag("timeresiduals")) {
    nest.Service<TimeResidualGenerator>("residualGenerator")
      ("eventName",  eventName)
      ("detectorService", "det")
      ("curvatureModel","curvatureCorrection")
      ("output", cl.GetArgument<string>("timeresidualsoutput"))
      ("angleFitName", planeName)
      ("coreFitName", coreName)
      ("timeResidualFileName",cl.GetArgument<string>("resi-calfit-file"))
      ("histfilename", cl.GetArgument<string>("timeresidualshistfilename"))
      ("maxEventsPerRun", cl.GetArgument<int>("timeresidualsmaxevents"))
      ("minHitsPerChannel", cl.GetArgument<int>("timeresidualsminhits"))
      ("thetaMin", cl.GetArgument<float>("timeresidualsminzenith") * HAWCUnits::degree)
      ("thetaMax", cl.GetArgument<float>("timeresidualsmaxzenith") * HAWCUnits::degree)
      ("phiRef", cl.GetArgument<float>("timeresidualsazimuth") * HAWCUnits::degree)
      ("phiDelta", cl.GetArgument<float>("timeresidualsdeltaazimuth") * HAWCUnits::degree)
      ("nChanMin", cl.GetArgument<int>("timeresidualsminimumnchannels"));
    chain.push_back("residualGenerator");
  }
  
  if (cl.HasFlag("chargeresiduals")) {
    nest.Service<ChargeResidualGenerator>("chargeResidualGenerator")
      ("event", eventName)
      ("detectorService", "det")
      ("coreFitName", coreName)
      ("angleFitName", planeName)
      ("histfilename", cl.GetArgument<string>("chargeresidualshistfilename"))
      ("nChanMin", cl.GetArgument<int>("chargeresidualsminimumnchannels"))
      ("minCharge", cl.GetArgument<double>("chargeresidualsminimumncharge"))
      ("maxEventsPerRun",  cl.GetArgument<int>("chargeresidualsmaxevents"))
      ("minHitsPerChannel",  cl.GetArgument<int>("chargeresidualsminhits"))
      ("hitSelection", stdCutsName);
    chain.push_back("chargeResidualGenerator");
  }
  
  //Set up OR diagnostics if file is specified.
  string orDiagFile = cl.GetArgument<string>("or-diag-file");
  if (orDiagFile != "") {
    nest.Service<ORDiagnostics>("orDiagnostics")
      ("event", eventName)
      ("outfile",orDiagFile)
      ("detectorService", "det");
    chain.push_back("orDiagnostics");
  }
  
  
  
  // Set up Zenith alignment if requested.
  string zenAlign = "";
  if (hasZenithAlign) {
    // Zenith alignment service
    nest.Service<ZenithAlignment>("zenithAlign")
      ("dataLocation", cl.GetArgument<string>("zenith-alignment-file"));
    zenAlign = "zenithAlign";
  }
  
  //lenght scales in meters
  string SOSName = "SOS";
  if (!cl.HasFlag("disable-sos")) {
    
    string sosFile = cl.GetArgument<std::string>("sos-template-file");
    if (sosFile.empty()) {
      log_info("no template file given, will try to get one from CONFIG_HAWC");
      const char* const env = getenv("CONFIG_HAWC");
      if (!env)
        log_fatal("CONFIG_HAWC undefined.");
      string configHAWCDir = env;
      sosFile = env + string("/reconstruction/sos/sos_templates.root");
      log_info("Using from CONFIG_HAWC: " << sosFile );
    }
    
    ifstream ifsFileCheck(sosFile.c_str());
    if (!ifsFileCheck.good()) {
      log_fatal("SOS template file not found: " << sosFile.c_str());
    }
    
    log_info("template file: " << sosFile );
    
    nest.Service<SOSCalculator>("sosCalc")
      ("detectorService","det")
      ("hitSelection", ppCutsName)
      ("lengthScale",10.)
      ("detectorService","det")
      ("event", "event")
      ("resultName", SOSName)
      ("tableFile",cl.GetArgument<std::string>("sos-table-file"))
      ("templateFile",sosFile)
      ("lhDebugFile",cl.GetArgument<std::string>("sos-debug-file"))
      ("angleFitName", angleName)
      ("coreFit",coreName);
    chain.push_back("sosCalc");
    resultsExtended.push_back(SOSName);
    log_info("added SOSCalculator");
  }
  //GBT calculation
  if (!cl.HasFlag("no-gbt")) {
    if(!cl.HasFlag("disable-sos") && !cl.HasFlag("disable-lh-lat-dist-fit")) {
      nest.Service<addgbt>("gdtcalc")
      ("detectorService",  "det")
      ("event", eventName)
      ("coreFitName", coreName)
      ("LHLatDistFlag", !cl.HasFlag("disable-lh-lat-dist-fit"))
      ("angleFitName", angleName)
      ("latDistFitName", latResult)
      ("hitSelectionName", stdCutsName)
      ("CxPE40Name","compactness40")
      ("PINCName", pincnessName)
      ("GHSName", TGHSName)
      ("sosName",SOSName)
      ("lhEnergyResultName", lheneRecName)
      ("fiducialCoreResultName", fiducialCoreName)
      ("resultName", "GBT")
      ("FeaturesNNName","FeaturesNN");
      chain.push_back("gdtcalc");  
      log_info("added GBT");
    }
    else
      log_fatal("Not enable sos and LHLatdist but enable gbt");
  }

  // Add summary calculator
  nest.Service<SummarizeRec>("summarize")
    ("event", eventName)
    ("detectorService", "det")
    ("astroService", "astroService")
    ("thresholdResultName", "thresholdResult")
    ("mpfSplitterName", mpfResult)
    ("coreFitName", coreName)
    ("LHLatDistFlag", !cl.HasFlag("disable-lh-lat-dist-fit"))
    ("angleFitName", planeName)
    ("planeFitName",planeName)
    ("latDistFitName", latResult)
    ("nnEnergyResultName", nnEneRecName)
    ("hitSelectionName", stdCutsName)
    ("zenithAlignmentName", zenAlign)
    ("nhitName", nhitname)
    ("CxPE20Name","compactness20")
    ("CxPE30Name","compactness30")
    ("CxPE40Name","compactness40")
    ("CxPE50Name","compactness50")
    ("pairCompactnessResultName", pairCompactnessName)
    ("PINCName", pincnessName)
    ("GHSName", TGHSName)
    ("sosName",SOSName)
    ("PMTChi2Name", pmtChi2Name)
    ("lhEnergyResultName", lheneRecName)
    ("fiducialChargeResultName", fiducialChargeName)
    ("fiducialCoreResultName", fiducialCoreName)
    ("TriggerWindow",150*HAWCUnits::ns)
    ("resultName", "rec")
    ("curvatureModel","curvatureCorrection")
    ("FeaturesNNName","FeaturesNN")
    ("GBTName","GBT");

  chain.push_back("summarize");
  resultsSummary.push_back("rec");
  
  if (isMC) {
    // Add MC summary calculator
    nest.Service<SummarizeMC>("mcsummarize")
      ("event", eventName)
      ("detectorService", "det")
      ("coreFit",  coreName)
      ("angleFit",  planeName)
      ("resultName", "mc");
    
    chain.push_back("mcsummarize");
    resultsSummary.push_back("mc");
  }
  
  
  // Binary XCDF output
  nest.Service<BinaryWriter>("writer")
    ("outFile", cl.GetArgument<string>("output"))
    ("blockSize", 10000);
  
  if (isTrigMC) {
    nest.Service<SummarizeMCFields>("summarizeMCFields")
      ("dataMapName","TriggeredData")
      ("mcEventName","mc")
      ("prefix","mc")
      ("resultName","mc");
    chain.push_back("summarizeMCFields");
    resultsSummary.push_back("mc");
    
    nest.Service<SummarizeSimFields>("summarizeSimFields")
      ("dataMapName","TriggeredData")
      ("simEventName","SimEvent")
      ("prefix","SimEvent")
      ("resultName","SimEvent");
    chain.push_back("summarizeSimFields");
    resultsSummary.push_back("SimEvent");
  }
  
  // Extended output
  bool writeHitsFlag = ! cl.HasFlag("no-hits");
  if (cl.HasFlag("extended")) {
    resultsSummary.insert(resultsSummary.end(),
                          resultsExtended.begin(),
                          resultsExtended.end());
  } else {
    // Outputing only summary data
    // No hit
    // No hit selections
    writeHitsFlag = false;
    hitSelections.clear();
    hitSelectionsOR.clear();
  }

  //Serializer
  nest.Service<DynamicSerializer>(serializerName)
    ("readerName", readerName)
    ("writerName", "writer")
    ("dataNames", resultsSummary)
    ("hitSelections", hitSelections)
    ("hitSelectionsOR", hitSelectionsOR)
    ("writeHitData",  writeHitsFlag)
    ("writeRawData",  true)
    ("eventName",  eventName)
    ("detectorService", "det")
    ("coreFitName", coreName)
    ("writeSimComments", isMC)
    ("copySimComments",isTrigMC)
    ("writeRecoDeviations",isMC)
    ("curvatureCorrectTimeResiduals", ! cl.HasFlag("no-correct-curv-tr"))
    ("curvatureModel","curvatureCorrection");

  chain.push_back(serializerName);


  // Insert the user modules specified by the user.
  if (cl.HasFlag("usermod")){

    boost::property_tree::ptree params;
    boost::property_tree::ini_parser::read_ini(cl.GetArgument<string>("usermod"), params);

    string chainLocationStr = params.get<string>("usermod.chainLocation");
    string modulesStr = params.get<string>("usermod.modules");
    string libsStr = params.get<string>("usermod.libs");

    boost::trim(chainLocationStr);
    boost::trim(modulesStr);
    boost::trim(libsStr);
      
    vector<string> chainLocations;
    vector<string> modules;
    vector<string> libs;

    boost::split(chainLocations, chainLocationStr,
                 boost::is_any_of(" "), boost::token_compress_on);
    boost::split(modules, modulesStr,
                 boost::is_any_of(" "), boost::token_compress_on);
    boost::split(libs, libsStr, 
                 boost::is_any_of(" "), boost::token_compress_on);
    
    if (chainLocations.size() != modules.size() ||
        modules.size() != libs.size() ||
        libs.size() != chainLocations.size())
    {
      log_fatal("Wrong user mode config file");
    }

    for (int i = 0; i < modules.size(); ++i ) {

      void* v = dlopen(libs[i].c_str(), RTLD_NOW | RTLD_GLOBAL);
      if (!v){
        log_fatal("Couldn't load dynamic library " << libs[i].c_str());
      }
      
      nest.Service(modules[i], modules[i]);

      vector<string>::iterator loc = find(chain.begin(), chain.end(), chainLocations[i]);

      if (loc == chain.end())
      {
        log_fatal("Service " << chainLocations[i] << " not found in chain.");
      }


      loc++; //Insert after
      
      chain.insert(loc, modules[i]);

      log_info("Added user module '" << modules[i]
               << "' after module '" << chainLocations[i]
               << "' from the library '" <<  libs[i] << "'");
      
    }
  }
  
  //Display or regular loop

  if (cl.HasFlag("display")) {
    parseParamTweak(cl,nest);
    TApplication app("EventDisplay", &argc, argv);
    HAWCDisplay display(gClient->GetRoot(), 1.2*950, 1.2*800,
                        &nest, &chain, cl.GetArgument<double>("displaysleep"));
    app.Run();
  }
  else {
    nest.Service<SequentialMainLoop>("mainloop")
      ("modulechain", chain)
      ("source", readerName)
      ("terminationLimit", cl.GetArgument<int>("maxevents"));
    
    parseParamTweak(cl,nest);
    nest.Configure();
    
    string configSummary = cl.GetArgument<string>("config-summary");
    if(configSummary != ""){
      ofstream fout(configSummary.c_str());
      nest.dump_ini(fout);
    }
    
    MainLoop& main = GetService<MainLoop>("mainloop");
    main.Execute();
  }

  //Save AERIE version
  Writer& writer = GetService<Writer>("writer");
  stringstream versionComment;
  versionComment << "Aerie version: " << SoftwareVersion(AERIE_VERSION_CODE);
  writer.AddComment(versionComment.str());
  
  stringstream revComment;
  revComment << "Aerie trunk revision: " << AERIE_REPOSITORY_REVISION;
  writer.AddComment(revComment.str());
  
  nest.Finish();
  
  log_info("Finished!");
  
  return 0;
}
