//--------------------------------Script description------------------------------------------------
//Date: 7th Jan 2021 Author: Chris Pliatskas
// This is a script will be used to calculate the energy flow fluctuations among jets of increasing jet radius. This is meant as an MC comparison for the Energy flow analysis task, using JEWEL jets as input. Previous implementations were using trees as a data format and were focusing on only hte leading jet for each Rjet wthout any matching conditions. A major requirement was that events would only be saved if a suitable jet was found for every Rjet ("full chain restriction"). The purpose of this script is to loose the latter restriction and simulate the analysis process of the analysis task as best as possible. For that reason,histograms are chosen as a first and simpler approach with the option to extend to a tree implementation which provides more options.

//Including libraries


#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"

#include "fastjet/Selector.hh" //.......... Background Sutraction event by event
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"//.......... Background Sutraction event by event
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequenceAreaBase.hh"
#include "fastjet/contrib/ConstituentSubtractor.hh"

#include "TPDGCode.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "THn.h"

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include "TTree.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TList.h"
#include "TVector3.h"
#include "TMath.h"
#include "THnSparse.h"
#include "TNtuple.h"
#include "TString.h"
#include "TRandom3.h"
#include "TH1D.h"

//Setting flags and using namespaces

using std::cout;
using std::cerr;
using std::endl;
using std::vector;

static const int debug = 0;
static const int do_bkg = 0;//Enable when working with the Recoil sample
static const int charged_jets = 1;

//Secondary function to check if the particle is stable

int is_stable(const HepMC::GenParticle *part) {
// copied from AliStack::IsStable()
int pdg = abs(part->pdg_id());
  if(pdg>1000000000)return kTRUE;

  const Int_t kNstable = 18;
  Int_t i;


  Int_t pdgStable[kNstable] = {
    kGamma,             // Photon
    kElectron,          // Electron
    kMuonPlus,          // Muon 
    kPiPlus,            // Pion
    kKPlus,             // Kaon
    kK0Short,           // K0s
    kK0Long,            // K0l
    kProton,            // Proton 
    kNeutron,           // Neutron
    kLambda0,           // Lambda_0
    kSigmaMinus,        // Sigma Minus
    kSigmaPlus,         // Sigma Plus
    3312,               // Xsi Minus 
    3322,               // Xsi 
    3334,               // Omega
    kNuE,               // Electron Neutrino 
    kNuMu,              // Muon Neutrino
    kNuTau              // Tau Neutrino
  };

  Bool_t isStable = kFALSE;
  for (i = 0; i < kNstable; i++) {
    if (pdg == abs(pdgStable[i])) {
      isStable = kTRUE;
      break;
    }
  }

  return isStable;
}
//Secondary function to check if the particle is stable & charged

int is_stable_charged(const HepMC::GenParticle *part) {
  // copied from AliStack::IsStable() 
 int pdg = abs(part->pdg_id());
  if(pdg>1000000000)return kTRUE;

  const Int_t kNstableCharged = 9;
  Int_t i;


  Int_t pdgStableCharged[kNstableCharged] = {
    kElectron,          // Electron
    kMuonPlus,          // Muon 
    kPiPlus,            // Pion
    kKPlus,             // Kaon
    kProton,            // Proton 
    kSigmaMinus,        // Sigma Minus
    kSigmaPlus,         // Sigma Plus
    3312,               // Xsi Minus 
    3334                // Omega
  };

  Bool_t isStable = kFALSE;
  for (i = 0; i < kNstableCharged; i++) {
    if (pdg == abs(pdgStableCharged[i])) {
      isStable = kTRUE;
      break;
    }
  }
  return isStable;
}

//Secondary function to check if the particle is charged

int is_charged(const HepMC::GenParticle *part) {
  int abs_kf = abs(part->pdg_id());

  if (abs_kf==211 || abs_kf==321 || abs_kf==2212 || abs_kf==11 || abs_kf==13)
    return 1;
  else if (abs_kf != 22 && abs_kf!=111 && abs_kf!=130 && abs_kf!=2112 && abs_kf!=311 && abs_kf!=12 && abs_kf !=14 && abs_kf!=16)
    cout << " Unexpected particle: kf=" << abs_kf << endl;
  return 0;
}

//Secondary function to calculate azimuthal angle difference
float dphi(float phi1, float phi2) {
  float dphi=phi1-phi2;
  float pi = 3.14159;
  if (dphi < -pi)
    dphi+=2*pi;
  if (dphi > pi)
    dphi-=2*pi;
  return dphi;
}

void JetMatcher( std::vector <fastjet::PseudoJet> *genJets,const Int_t &kGenJets,
                 std::vector <fastjet::PseudoJet> *recJets,const Int_t &kRecJets,
                TArrayI &iGenIndex,TArrayI &iRecIndex,Int_t iDebug, Float_t maxDist){

 // Get the implementation from the analysis task

   // Size indepnedendentt Implemenation of jet matching
   // Thepassed TArrayI should be static in the user function an only increased if needed
   //
   // Relate the two input jet Arrays
   //
   // The association has to be unique
   // So check in two directions
   // find the closest rec to a gen
   // and check if there is no other rec which is closer
   // Caveat: Close low energy/split jets may disturb this correlation
   //
   // Idea: search in two directions generated e.g (a--e) and rec (1--3)
   // Fill a matrix with Flags (1 for closest rec jet, 2 for closest rec jet
   // in the end we have something like this
   //   1   2   3
   // ------------
   //a| 3   2   0
   //b| 0   1   0
   //c| 0   0   3
   //d| 0   0   1
   //e| 0   0   1
   // Topology
   //  1     2
   //    a         b        
   //
   //  d      c
   //       3     e
   // Only entries with "3" match from both sides

        iGenIndex.Reset(-1);
        iRecIndex.Reset(-1);
	
	Int_t size_gen = genJets->size();
	Int_t size_rec = recJets->size();
        const int kMode = 3;
        const Int_t nGenJets = TMath::Min(size_gen,kGenJets);
        const Int_t nRecJets = TMath::Min(size_rec,kRecJets);
        if(nRecJets==0||nGenJets==0)return;

        static TArrayS iFlag(nGenJets*nRecJets);
        if(iFlag.GetSize()<(nGenJets*nRecJets)){
        iFlag.Set(nGenJets*nRecJets);
        }
        iFlag.Reset(0);

       // find the closest distance to the generated
            for(int ig = 0;ig<nGenJets;++ig){
       fastjet::PseudoJet &genJet = genJets->at(ig);
         // if(!genJet)continue;

      Float_t dist = maxDist;
      if(iDebug>1)Printf("Gen (%d) p_T %3.3f eta %3.3f ph %3.3f ",ig,genJet.pt(),genJet.eta(),genJet.phi());
      for(int ir = 0;ir<nRecJets;++ir){
        fastjet::PseudoJet const &recJet = recJets->at(ir);
   	// if(!recJet)continue;

      if(iDebug>1){
      printf("Rec (%d) ",ir);
      Printf("p_T %3.3f eta %3.3f ph %3.3f ",recJet.pt(),recJet.eta(),recJet.phi());
          }
        Double_t dR = genJet.delta_R(recJet);
        if(iDebug>1)Printf("Distance (%d)--(%d) %g ",ig,ir,dR);

        if(dR<dist){
    iRecIndex[ig] = ir;
    dist = dR;
        }
      }
      if(iRecIndex[ig]>=0)iFlag[ig*nRecJets+iRecIndex[ig]]+=1;
      // Resetting...
	 iRecIndex[ig] = -1;
    }
    // other way around
    for(int ir = 0;ir<nRecJets;++ir){
       fastjet::PseudoJet const &recJet = recJets->at(ir);
     //        if(!recJet)continue;
               Float_t dist = maxDist;
               for(int ig = 0;ig<nGenJets;++ig){
          fastjet::PseudoJet &genJet = genJets->at(ig);
//	  if(!genJet)continue;
     Double_t dR = genJet.delta_R(recJet);
         if(dR<dist){
           iGenIndex[ir] = ig;
           dist = dR;
               }
             }
             if(iGenIndex[ir]>=0)iFlag[iGenIndex[ir]*nRecJets+ir]+=2;
             // Resetting...
            iGenIndex[ir] = -1;
           }

if(iDebug>1)Printf(">>>>>> Matrix Size %d",iFlag.GetSize());

     for(int ig = 0;ig<nGenJets;++ig){
       for(int ir = 0;ir<nRecJets;++ir){
         // Print
         if(iDebug>1)printf("Flag2[%d][%d] %d ",ig,ir,iFlag[ig*nRecJets+ir]);
        	 if(kMode==3){
         // we have a unique correlation
                if(iFlag[ig*nRecJets+ir]==3){
             iGenIndex[ir] = ig;
             iRecIndex[ig] = ir;
         }
         }
         else{
     // we just take the correlation from on side                 
           if((iFlag[ig*nRecJets+ir]&2)==2){
             iGenIndex[ir] = ig;
           }
           if((iFlag[ig*nRecJets+ir]&1)==1){
         iRecIndex[ig] = ir;
     }
               }
             }
             if(iDebug>1)printf("\n");
         }
};

int main(int argc, char **argv) {
//~~~~~~~~~~~~~~~~~~~~~Analysis Level (0th Level)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 if (argc < 2) {
        cerr << "Need two arguments: infile outfile" << endl << "infile is HEPMC ascii format; outfile will be root format" << endl;
        return 1;
  }

// Input and Output name manipulation
char *inname = argv[1];
  // specify an input file
HepMC::IO_GenEvent ascii_in(inname,std::ios::in);

//Specify the name of the output file
 char *outbase = argv[2];
  TString outlabel;
  if (charged_jets == 1) {outlabel += "_chargedNEW";}
  if (do_bkg == 1) {outlabel += "_bkgsub";}
TString outname(Form("%s%s.root",outbase,outlabel.Data()));
  cout << "Input: " << inname << ", output " << outname << endl;

TFile fout(outname,"RECREATE");
 static const int nR = 4; // radii 0.2, 0.4, 0.6, 0.8

  TH1F *hJetPt[nR] = {0};
  TH1F *hJetEta[nR] = {0};
  TH1F *hJetPtMatched[nR] = {0};
  TH1F *hJetEtaMatched[nR] = {0};
  TH2F *hPtJetDeltaPt[nR] = {0};
  TH2F *hPtJetDeltaR[nR] = {0};
  TH3F *hPtJetDeltaRDeltaPt[nR] = {0};
  TH2F *hPtJetDeltaPtoverPt_low[nR] = {0};
  TH2F *hPtJetDeltaPtoverPt_high[nR] = {0};
  TH2F *hEtaJetDeltaR[nR] = {0};  
  TH3F *hDptPtDEta[nR] = {0};
  TH3F *hDptPtEta_low[nR] = {0};
  TH3F *hDptPtEta_high[nR] = {0};  
  TH3F* hDptPtMultiplicity[nR] = {0};

//Binning choices;
 const int NPtBins = 100;
 const int PtBin_min = 0;
 const int PtBin_max = 100;
 const int NEtaBins = 100;
 const float EtaBin_max = 1.0;
 const float maxEtapart = 2;
 const int NDptBins = 100;
 const int DPtBin_min = -20;
 const int DPtBin_max = 80;
 const int NDRBins = 100;
 const float DRBin_min = 0;
 const float DRBin_max = 0.4;

//General histograms
  TH1F *hNEvent = new TH1F("hNEvent","number of events; N",1,0,1);
  hNEvent->Sumw2();
  TH2D *hPtPartEta = new TH2D("hPtPartEta","Particle p_{t} and #eta distribution ;p_{t};#eta",100,0,100,NEtaBins,-maxEtapart,maxEtapart);
  hPtPartEta->Sumw2();
 
  static float Rstep = 0.1;
  TString histname;
  TString htitle;
/* Bkg subtracted histograms applicable to Non-pp case- Not implemented for now
   TH2F *hbkgdensity = new TH2F("hbkgdensity","main jet rho,matched jet rho;main;matched",100,0,30, 100,0,30);
  hbkgdensity->Sumw2();
*/


// Creation of R-dependent histograms
for (Int_t iR=0; iR<nR; iR++){
	Float_t Rjet = (iR+1)*Rstep;

	histname = TString::Format("hJetPt_R%02d",int(Rjet*10));
	htitle = TString::Format("Jet p_{t} spectrum for R=%.1f;p_{t} (GeV/c)",Rjet);
	hJetPt[iR] = new TH1F(histname,htitle,NPtBins,PtBin_min,PtBin_max);

        histname = TString::Format("hJetEta_R%02d",int(Rjet*10));
        htitle = TString::Format("Jet #eta distribution for R=%.1f;#eta",Rjet);
        hJetEta[iR] = new TH1F(histname,htitle,NEtaBins,-EtaBin_max,EtaBin_max);

        histname = TString::Format("hJetPtMatched_R%02d",int(Rjet*10));
        htitle = TString::Format("Matched jet p_{t} spectrum for R=%.1f;p_{t} (GeV/c)",Rjet);
        hJetPtMatched[iR] = new TH1F(histname,htitle,NPtBins,PtBin_min,PtBin_max);

        histname = TString::Format("hJetEtaMatched_R%02d",int(Rjet*10));
        htitle = TString::Format("Matched jet #eta distribution for R=%.1f;#eta",Rjet);
        hJetEtaMatched[iR] = new TH1F(histname,htitle,NEtaBins,-EtaBin_max,EtaBin_max);

	if(iR<nR-1){
	histname = TString::Format("hPtJetDeltaPt_R%02d",int(Rjet*10));
        htitle = TString::Format("#Deltap_{t} between R=%.1f and R=%.1f;p_{t,R=%.1f} (GeV/c);#Delta P_{t} (GeV/c)",Rjet,Rjet+Rstep,Rjet);
        hPtJetDeltaPt[iR] = new TH2F(histname,htitle,NPtBins,PtBin_min,PtBin_max,NDptBins,DPtBin_min,DPtBin_max);

        histname = TString::Format("hPtJetDeltaRDeltaPt_R%02d",int(Rjet*10));
        htitle = TString::Format("#Deltap_{t} vs #DeltaR between R=%.1f and R=%.1f;p_{t,R=%.1f} (GeV/c);#DeltaR;#Deltap_{t} (GeV/c)",Rjet,Rjet+Rstep,Rjet);
        hPtJetDeltaRDeltaPt[iR] = new TH3F(histname,htitle,NPtBins,PtBin_min,PtBin_max,NDRBins,DRBin_min,DRBin_max,NDptBins,DPtBin_min,DPtBin_max);

        histname = TString::Format("hPtJetDeltaPtoverPt_R%02d_low",int(Rjet*10));
        htitle = TString::Format("#Deltap_{t} between R=%.1f and R=%.1f over the p_{t} of R=%.1f;p_{t,R=%.1f} (GeV/c);#Deltap_{t}/p_{t}",Rjet,Rjet+Rstep,Rjet,Rjet);
        hPtJetDeltaPtoverPt_low[iR] = new TH2F(histname,htitle,NPtBins,PtBin_min,PtBin_max,100,0,2);

        histname = TString::Format("hPtJetDeltaPtoverPt_R%02d_high",int((Rjet+Rstep)*10));
        htitle = TString::Format("#Deltap_{t} between R=%.1f and R=%.1f over the p_{t} of R=%.1f;p_{t,R=%.1f} (GeV/c);#Deltap_{t}/p_{t}",Rjet,Rjet+Rstep,Rjet+Rstep,Rjet);
        hPtJetDeltaPtoverPt_high[iR] = new TH2F(histname,htitle,NPtBins,PtBin_min,PtBin_max,100,0,2);

        histname = TString::Format("hEtaJetDeltaR_R%02d",int(Rjet*10));
        htitle = TString::Format("#DeltaR between R=%.1f and R=%.1f;#eta{R=%.1f};#DeltaR",Rjet,Rjet+Rstep,Rjet);
        hEtaJetDeltaR[iR] = new TH2F(histname,htitle,NEtaBins,-EtaBin_max,EtaBin_max,NDRBins,DRBin_min,DRBin_max);

        histname = TString::Format("hPtJetDeltaR_R%02d",int(Rjet*10));
        htitle = TString::Format("#DeltaR between R=%.1f and R=%.1f;p_{t,R=%.1f} (GeV/c);#DeltaR",Rjet,Rjet+Rstep,Rjet);
        hPtJetDeltaR[iR] = new TH2F(histname,htitle,NPtBins,PtBin_min,PtBin_max,NDRBins,DRBin_min,DRBin_max);
	
	histname = TString::Format("hDptPtDEta_R%02d",int(Rjet*10));
        htitle = TString::Format("#DeltaP_{t} between R=%.1f and R=%.1fvs P_{t} vs #Delta#eta;#Deltap_{t,R=%.1f} (GeV/c);p_{t,R=%.1f} (GeV/c);#Delta #Eta",Rjet,Rjet+Rstep,Rjet);
        hDptPtDEta[iR] = new TH3F(histname,htitle,NDptBins,DPtBin_min,DPtBin_max,NPtBins,PtBin_min,PtBin_max,200,0,2);
	
	histname = TString::Format("hDptPtEta_lowR%02d",int(Rjet*10));
        htitle = TString::Format("#DeltaP_{t} between R=%.1f and R=%.1f vs P_{t} vs #eta;#Deltap_{t,R=%.1f} (GeV/c);p_{t,R=%.1f} (GeV/c); #Eta_{R=%.1f}",Rjet,Rjet+Rstep,Rjet,Rjet);
	hDptPtEta_low[iR] = new TH3F(histname,htitle,NDptBins,DPtBin_min,DPtBin_max,NPtBins,PtBin_min,PtBin_max,NEtaBins,-EtaBin_max,EtaBin_max);

        histname = TString::Format("hDptPtEta_highR%02d",int(Rjet*10));
        htitle = TString::Format("#DeltaP_{t} between R=%.1f and R=%.1fvs P_{t} vs #eta;#Deltap_{t,R=%.1f} (GeV/c);p_{t,R=%.1f} (GeV/c); #Eta_{R=%.1f}",Rjet,Rjet+Rstep,Rjet,Rjet+Rstep);
        hDptPtEta_high[iR] = new TH3F(histname,htitle,NDptBins,DPtBin_min,DPtBin_max,NPtBins,PtBin_min,PtBin_max,NEtaBins,-EtaBin_max,EtaBin_max);
	
	histname = TString::Format("hDptPtMultiplicity_R%02d",int(Rjet*10));
	htitle = TString::Format("#DeltaP_{t} between R=%.1f and R=%.1f vs P_{t} vs multiplicity;#Deltap_{t,R=%.1f} (GeV/c);p_{t,R=%.1f} (GeV/c); Multiplicity",Rjet,Rjet+Rstep,Rjet,Rjet+Rstep);
	hDptPtMultiplicity[iR]= new TH3F(histname,htitle,NDptBins,DPtBin_min,DPtBin_max,NPtBins,PtBin_min,PtBin_max,20,0,20);
	}
}

HepMC::GenEvent* evt = ascii_in.read_next_event();
  if (!evt) cerr << "Input file not found " << inname << endl;

//Event loop

 while(evt)
{

if (debug)cout << "Event " << endl;

    hNEvent->Fill(0.5,evt->weights()[0]); // count events

    float max_eta_track = 2; //2.8;
    float min_pt =10;
    float max_eta_jet = 0.9-Rstep*nR;
    int index = 0;
    std::vector <fastjet::PseudoJet> Tracks_in;
    for ( HepMC::GenEvent::particle_iterator pit = evt->particles_begin();
          pit != evt->particles_end(); ++pit )
      {
        const HepMC::GenParticle *p = *pit;
        if ((!p->end_vertex() && p->status()==1 && (!charged_jets || is_charged(p)))&& fabs(p->momentum().eta()) < max_eta_track) hPtPartEta->Fill(p->momentum().perp(), p->momentum().eta(),evt->weights()[0]);
        	if(p->momentum().perp()==0) continue;
                double mom = sqrt(p->momentum().x()*p->momentum().x() +
                                  p->momentum().y()*p->momentum().y() +
                                  p->momentum().z()*p->momentum().z());
                fastjet::PseudoJet jInp(p->momentum().x(),p->momentum().y(),p->momentum().z(),mom);
                jInp.set_user_index(index);
                Tracks_in.push_back(jInp);



                index++;
       }

//Jet Finding
	fastjet::GhostedAreaSpec ghostSpec(max_eta_track,1,0.01);
	fastjet::Strategy               strategy = fastjet::Best;
	fastjet::RecombinationScheme    recombScheme = fastjet::BIpt_scheme;
	fastjet::AreaType areaType =   fastjet::active_area;
	fastjet::AreaDefinition areaDef = fastjet::AreaDefinition(areaType,ghostSpec);

	vector <fastjet::PseudoJet> jets[nR];
	vector <fastjet::PseudoJet> AcceptedJets[nR];
	fastjet::ClusterSequenceArea *clustSeqCh[nR]={0};

	for (int iR =0; iR < nR; iR++) {
      	float jetR = 0.2+0.2*iR;
      	fastjet::JetDefinition jetDefCh(fastjet::antikt_algorithm, jetR,recombScheme, strategy);
      	clustSeqCh[iR]=new fastjet::ClusterSequenceArea(Tracks_in, jetDefCh,areaDef);
      	jets[iR] = clustSeqCh[iR]->inclusive_jets();
	// Collection of jets before acceptance cut
	for (auto j:jets[iR]) if(fabs(j.eta()) < max_eta_jet)AcceptedJets[iR].push_back(j);
	// Collection of jets after acceptance cut
	for (auto k:AcceptedJets[iR]){
		hJetPt[iR]->Fill(k.pt(),evt->weights()[0]);
		hJetEta[iR]->Fill(k.eta(),evt->weights()[0]);		
				     }  }
//Jet matching
	Double_t DeltaPt = 0.0;
	Double_t DeltaR = 0.0 ;
	Double_t DeltaEta = 0.0;
	Int_t Multiplicity = 0;
	Int_t Njets = 200;      //Just a high number so that the matching matrix created by the matcher task will always have the size of the input jet lists
	Int_t &kLowRJets = Njets;
	Int_t &kHighRJets = Njets;

//This array points to the low R jet that matches to each high R jet
	TArrayI iLowRIndex;

// This array points to the high R jet that matches to each low R jet
	TArrayI iHighRIndex;

	for (int iR =0; iR < nR-1; iR++) {
	
	Double_t Maxdist = 0.2;
	
        iLowRIndex.Set(AcceptedJets[iR+1].size());
        iHighRIndex.Set(AcceptedJets[iR].size());

        if(AcceptedJets[iR].size()==0||AcceptedJets[iR+1].size()==0) continue;
	
        JetMatcher(&AcceptedJets[iR],kLowRJets,&AcceptedJets[iR+1],kHighRJets, iLowRIndex,iHighRIndex,0,Maxdist);

	//Fill histograms
	for (Int_t j=0; j<iHighRIndex.GetSize()-1;j++)
                {
        if(iHighRIndex[j]>=0){
		Int_t match_index = iHighRIndex[j];
  //If a match exists then calculate the DeltaPt and DeltaEta quantities before filling the histograms
        if (iLowRIndex[match_index]==j){
		DeltaPt =AcceptedJets[iR+1].at(match_index).pt() -AcceptedJets[iR].at(j).pt();
	
//Multiplicity calculation
		Multiplicity =AcceptedJets[iR].at(j).constituents().size();

	if (AcceptedJets[iR].at(j).eta()*AcceptedJets[iR+1].at(match_index).eta()>=0)
		DeltaEta = fabs(AcceptedJets[iR+1].at(match_index).eta()-AcceptedJets[iR].at(j).eta());
	else{if(AcceptedJets[iR].at(j).eta()>0)
		DeltaEta = AcceptedJets[iR].at(j).eta() - AcceptedJets[iR+1].at(match_index).eta();
	     else DeltaEta =AcceptedJets[iR+1].at(match_index).eta() -AcceptedJets[iR].at(j).eta();
	    }

		hDptPtMultiplicity[iR]->Fill(DeltaPt,AcceptedJets[iR].at(j).pt(),Multiplicity,evt->weights()[0]);
                DeltaR = AcceptedJets[iR+1].at(match_index).delta_R(AcceptedJets[iR].at(j));
		hJetPtMatched[iR]->Fill(AcceptedJets[iR].at(j).pt(),evt->weights()[0]);
		hJetEtaMatched[iR]->Fill(AcceptedJets[iR].at(j).eta(),evt->weights()[0]);
		if(iR==nR-2){
			hJetPtMatched[iR+1]->Fill(AcceptedJets[iR+1].at(match_index).pt(),evt->weights()[0]);
			hJetEtaMatched[iR+1]->Fill(AcceptedJets[iR+1].at(match_index).eta(),evt->weights()[0]);	
				}
	
		hDptPtDEta[iR]->Fill(DeltaPt,AcceptedJets[iR].at(j).pt(),DeltaEta,evt->weights()[0]);
		hDptPtEta_low[iR]->Fill(DeltaPt,AcceptedJets[iR].at(j).pt(),AcceptedJets[iR].at(j).eta(),evt->weights()[0]);
		hDptPtEta_high[iR]->Fill(DeltaPt,AcceptedJets[iR].at(j).pt(),AcceptedJets[iR+1].at(match_index).eta(),evt->weights()[0]);

	
		hPtJetDeltaPt[iR]->Fill(AcceptedJets[iR].at(j).pt(),DeltaPt,evt->weights()[0]);
		hEtaJetDeltaR[iR]->Fill(AcceptedJets[iR].at(j).eta(),DeltaR,evt->weights()[0]);
		hPtJetDeltaRDeltaPt[iR]->Fill(AcceptedJets[iR].at(j).pt(),DeltaR,DeltaPt,evt->weights()[0]);
		hPtJetDeltaPtoverPt_low[iR]->Fill(AcceptedJets[iR].at(j).pt(),DeltaPt/AcceptedJets[iR].at(j).pt(),evt->weights()[0]);
		hPtJetDeltaPtoverPt_high[iR]->Fill(AcceptedJets[iR].at(j).pt(),DeltaPt/AcceptedJets[iR+1].at(match_index).pt(),evt->weights()[0]);
		hPtJetDeltaR[iR]->Fill(AcceptedJets[iR].at(j).pt(),DeltaR,evt->weights()[0]);
						}
					}
		}// End of loop over the lowR jets
		}// End of loop over the jet pairs

	//Clean up
	for (int iR =0; iR < nR; iR++) {
      	delete clustSeqCh[iR];
      	clustSeqCh[iR]=0;
   	 }
   
 // delete the created event from memory
     delete evt;
// read the next event
     ascii_in >> evt;
} //End of event loop

//----------------------------------------------------------
  fout.Write();

  fout.Close();

  return 0;
}
