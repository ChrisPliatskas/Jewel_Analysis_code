//--------------------------------Script description--------------------------------------------------------------------------------------------------------------------------------------------------------
// Study of the jet's pt dependence on jet radius using JEWEL jets embedded in a thermal background as input. An increase of pt is expected (due to increased jet area) but deviations from this behaviour could lead to identifying a new mechanism/ effect. The goal of this script is to provide the data for a Gaussian smearing study on AA data (with and without recoil) as well as pp.
// Notes: Trees were chosen as the preferred output format due to their plotting flexibility, howeveras a first approximation I chose to focus on the hardest jet of each event.Since one jet will be saved per event I relaxed the matching requirement of the jet for different radii but this could be also a possible mini result.

// Including libraries

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



using std::cout;
using std::cerr;
using std::endl;
using std::vector;

//Definition of supplemental functions
int is_charged(const HepMC::GenParticle *part) {
  int abs_kf = abs(part->pdg_id());

  if (abs_kf==211 || abs_kf==321 || abs_kf==2212 || abs_kf==11 || abs_kf==13)
    return 1;
  else if (abs_kf != 22 && abs_kf!=111 && abs_kf!=130 && abs_kf!=2112 && abs_kf!=311 && abs_kf!=12 && abs_kf !=14 && abs_kf!=16)
    cout << " Unexpected particle: kf=" << abs_kf << endl;
  return 0;
}

//--------------------------------------Main function---------------------------

//Flag definition

//using std::cout;
//using std::cerr;
//using std::endl;
//using std::vector;

static const int debug = 0;
static const int do_bkg = 1;
static const int charged_jets = 1;

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

 char *outbase = argv[2];
  TString outlabel;
  if (charged_jets == 1) {outlabel += "_charged";}
  if (do_bkg == 1) {outlabel += "_bkgsub";}
TString outname(Form("%s%s.root",outbase,outlabel.Data()));
  cout << "Input: " << inname << ", output " << outname << endl;

TFile fout(outname,"RECREATE");

//Analysis settings ( Rjet, Pt/eta cuts, number of thermal particles)
static const int nR = 4; // radii 0.2, 0.4, 0.6, 0.8
float max_eta_track = 2; //2.8;
float min_pt =10;
float jet_pt_max =80.0;
float jet_pt_min =40.0;
double Rjet=0.2;
float Rho=0.0;
float Rho_gen =0.0;
// Tree definition
typedef struct {Double_t pt,eta,phi;} JET;

JET R02, R04, R06, R08, jew_R02, jew_R04, jew_R06, jew_R08;
R02.pt = 0.0; R02.eta = -4.0; R02.phi = -4.0;
R04.pt = 0.0; R04.eta = -4.0; R04.phi = -4.0;
R06.pt = 0.0; R06.eta = -4.0; R06.phi = -4.0;
R08.pt = 0.0; R08.eta = -4.0; R08.phi = -4.0;

jew_R02.pt = 0.0; jew_R02.eta = -4.0; jew_R02.phi = -4.0;
jew_R04.pt = 0.0; jew_R04.eta = -4.0; jew_R04.phi = -4.0;
jew_R06.pt = 0.0; jew_R06.eta = -4.0; jew_R06.phi = -4.0;
jew_R08.pt = 0.0; jew_R08.eta = -4.0; jew_R08.phi = -4.0;

double Avg_pt_b[nR] = {0};
double Avg_pt_a[nR] = {0}; 


TTree* tree = new TTree("event_tree","event information");
tree->Branch("R02",&R02,"pt/D:eta:phi");
tree->Branch("R04",&R04,"pt/D:eta:phi");
tree->Branch("R06",&R06,"pt/D:eta:phi");
tree->Branch("R08",&R08,"pt/D:eta:phi");
tree->Branch("jew_R02",&jew_R02,"pt/D:eta:phi");
tree->Branch("jew_R04",&jew_R04,"pt/D:eta:phi");
tree->Branch("jew_R06",&jew_R06,"pt/D:eta:phi");
tree->Branch("jew_R08",&jew_R08,"pt/D:eta:phi");
tree->Branch("Rho", &Rho);
tree->Branch("Rho_gen", &Rho_gen);
tree->Branch("Avg_pt_b",&Avg_pt_b,"R02/D:R04:R06:R08");
tree->Branch("Avg_pt_a",&Avg_pt_a,"R02/D:R04:R06:R08");
//Thermal distribution
	int nThermalParticles = 4000;
        TF1* f_pT = new TF1("f_pT","x*exp(-x/0.3)", 0.0, 400.0);
        f_pT->SetNpx(40000);

        TF1* f_eta = new TF1("f_eta", "1", -2.0, 2.0);
        f_eta->SetNpx(200);

        TF1* f_phi = new TF1("f_phi", "1", (-1.0)*TMath::Pi(), TMath::Pi() );
        f_phi->SetNpx(700);

//~~~~~~~~~~~~~~~~~~~~Start of the event Loop (1st Level)~~~~~~~~~~~~~~~~~~~~~~ 

HepMC::GenEvent* evt = ascii_in.read_next_event();
if (!evt) cerr << "Input file not found " << inname << endl;
while (evt)
{

	R02.pt = 0.0; R02.eta = -4.0; R02.phi = -4.0;
	R04.pt = 0.0; R04.eta = -4.0; R04.phi = -4.0;
	R06.pt = 0.0; R06.eta = -4.0; R06.phi = -4.0;
        R08.pt = 0.0; R08.eta = -4.0; R08.phi = -4.0;
	jew_R02.pt = 0.0; jew_R02.eta = -4.0; jew_R02.phi = -4.0;
	jew_R04.pt = 0.0; jew_R04.eta = -4.0; jew_R04.phi = -4.0;
	jew_R06.pt = 0.0; jew_R06.eta = -4.0; jew_R06.phi = -4.0;
	jew_R08.pt = 0.0; jew_R08.eta = -4.0; jew_R08.phi = -4.0;

	float rho= 0.0; float rho_gen = 0.0;

	double fourvec[4];
	std::vector <fastjet::PseudoJet> fjInputs;
	std::vector <fastjet::PseudoJet> fjInputs2;
	bool foundjet = 0;
//Jewel particle definition
	for ( HepMC::GenEvent::particle_iterator pit = evt->particles_begin();pit != evt->particles_end(); ++pit )
	{ const HepMC::GenParticle *p = *pit;
	  if (!p->end_vertex() && p->status()==1 && (!charged_jets || is_charged(p)))
	  {if ( fabs(p->momentum().eta()) < max_eta_track )
	   { double mom = sqrt(p->momentum().x()*p->momentum().x() + p->momentum().y()*p->momentum().y() + p->momentum().z()*p->momentum().z());
             fastjet::PseudoJet Jewel_in(p->momentum().x(),p->momentum().y(),p->momentum().z(),mom);
             fjInputs.push_back(Jewel_in);
           }
	  }
	}

//cout<<"Checkpoint 1:Jewel particles"<<endl;

	
//Thermal particle definition

	for(int j = 0; j < nThermalParticles; j++)
                {
                        double pT = f_pT->GetRandom();

                        double Eta = f_eta->GetRandom();

                        double Phi = f_phi->GetRandom();

                       // if(pT < trackLowPtCut) continue;//pt cut

                        fourvec[0] = pT*TMath::Cos(Phi);
                        fourvec[1] = pT*TMath::Sin(Phi);
                        fourvec[2] = pT*TMath::SinH(Eta);
                        fourvec[3] = TMath::Sqrt(pT*pT + fourvec[2]*fourvec[2]);

                        fastjet::PseudoJet ThermalParticle(fourvec);

                        ThermalParticle.set_user_index(0);


            fjInputs2.push_back(ThermalParticle);
                }

	double therm_sum = 0.0;
        double Debug_sum = 0.0;
	for(unsigned int j=0;j<fjInputs2.size();j++){therm_sum +=fjInputs2[j].pt();} 

//	cout<< "Event start, number of particles: "<<fjInputs2.size()<<" with pt sum "<<therm_sum<<endl;



//cout<<"Checkpoint 2: Thermal particles"<<endl;
//Jet Finding algorithm definitions
        fastjet::GhostedAreaSpec ghostSpec(max_eta_track,1,0.01);
        fastjet::Strategy               strategy = fastjet::Best;
        fastjet::RecombinationScheme    recombScheme = fastjet::BIpt_scheme;
        fastjet::AreaType areaType =   fastjet::active_area;
        fastjet::AreaDefinition areaDef = fastjet::AreaDefinition(areaType,ghostSpec);

        vector <fastjet::PseudoJet> jets[nR];// = {{},{},{},{}};
	vector <fastjet::PseudoJet> BGJets[nR]; // = {{},{},{},{}};
	vector <fastjet::PseudoJet> NewJets[nR]; //= {{},{},{},{}};
       	fastjet::PseudoJet ProbeJet[nR];
	fastjet::ClusterSequenceArea *clustSeqCh[nR]={0};
	fastjet::ClusterSequenceArea *clustSeqBG[nR] = {0};
	double avg_pt_before[nR] = {0};
	double avg_pt_after[nR] = {0};

        double pt_sum=0;
        for(unsigned int j=0;j<fjInputs2.size();j++){
                pt_sum += fjInputs2[j].pt();    
        }
	double temp = pt_sum;
//	cout<<"Thermal particles pt sum is: "<<temp<<endl;
//~~~~~~~~~~~~~~~~~~~Start of the Rjet Loop (2nd Level)~~~~~~~~~~~~~~~~~~~~~~~~~

	for (int iR =0; ((iR < nR)&&(foundjet))||(iR==0); iR++) 
	{
	float jetR =Rjet+iR*Rjet;
	foundjet = 0;
//Jet finding
	fastjet::JetDefinition jetDefCh(fastjet::antikt_algorithm, jetR,recombScheme, strategy);
        clustSeqCh[iR] = new fastjet::ClusterSequenceArea(fjInputs, jetDefCh,areaDef);

                jets[iR] = sorted_by_pt(clustSeqCh[iR]->inclusive_jets(1.));
		if(jets[iR].size()==0) continue;
//Probe jet definition
		ProbeJet[iR]= jets[iR][0];
//Probe jet cuts
		if(ProbeJet[iR].pt()<jet_pt_min) continue;
        	if(iR==0){if(ProbeJet[iR].pt()>jet_pt_max) continue;}
		if(ProbeJet[iR].eta()>(max_eta_track-jetR)) continue;
		if(ProbeJet[iR].eta()<-(max_eta_track-jetR)) continue; 
		foundjet=1;
// Probe + Thermal particles sample
	std::vector<fastjet::PseudoJet> ProbeParticles = sorted_by_pt(ProbeJet[iR].constituents());
	for(unsigned int i=0;i<ProbeParticles.size();i++){
		ProbeParticles[i].set_user_index(1);
		fjInputs2.push_back(ProbeParticles[i]);
							}
	pt_sum = 0;						 
	for(unsigned int j=0;j<fjInputs2.size();j++){
		if(fjInputs2[j].user_index()==1) pt_sum += fjInputs2[j].pt();						
	}	
//	cout<<"Before sum:Probe: "<<pt_sum<<" and thermal:"<<temp<<endl;
	avg_pt_before[iR] = (pt_sum+temp)/(2*TMath::Pi()*2*max_eta_track);
//	cout<<"Before average: "<< avg_pt_before[iR]<<"for R ="<<jetR<< endl;

	Debug_sum = 0;
	for(unsigned int j=0;j<fjInputs2.size();j++){Debug_sum += fjInputs2[j].pt();}

  //      cout<<"Debug sum (Thermal+probe):"<<Debug_sum<<endl;
	

//cout<<"Checkpoint 3: Probe particles"<<endl;
//Bkg subtraction
	fastjet::JetMedianBackgroundEstimator bge;
	fastjet::Selector BGSelector = fastjet::SelectorAbsEtaMax(2.0);
        fastjet::JetDefinition jetDefBG(fastjet::kt_algorithm, jetR, recombScheme, strategy);
        fastjet::AreaDefinition fAreaDefBG(fastjet::active_area_explicit_ghosts,ghostSpec);
        clustSeqBG[iR] = new fastjet::ClusterSequenceArea(fjInputs2, jetDefBG,fAreaDefBG);
        BGJets[iR] = clustSeqBG[iR]->inclusive_jets();
	if(BGJets[iR].size()==0) cout<<"Error: No BGJets found at "<<jetR<<endl;
        bge.set_selector(BGSelector);
        bge.set_jets(BGJets[iR]);
        fastjet::contrib::ConstituentSubtractor subtractor(&bge);
 
       subtractor.set_common_bge_for_rho_and_rhom(true);
	rho_gen= bge.rho(ProbeJet[iR]);
	rho = bge.rho();
	subtractor.set_max_standardDeltaR(0.5*jetR);
	std::vector<fastjet::PseudoJet> fjInputs3 = subtractor.subtract_event(fjInputs2, max_eta_track);

	 pt_sum = 0;                                    
        for(unsigned int j=0;j<fjInputs3.size();j++){
                if(fjInputs3[j].user_index()==1) pt_sum += fjInputs3[j].pt();
        }
//	cout<<"After sum: Probe: "<<pt_sum<<"and thermal:"<<temp<<endl;
        avg_pt_after[iR] = (pt_sum+temp)/(2*TMath::Pi()*2*max_eta_track);
//	cout<<"After average: "<< avg_pt_after[iR]<<"for R ="<<jetR<< endl;



//	cout<<"Checkpoint 4: Bkg subtraction"<<endl;
// New jet finding
	fastjet::GhostedAreaSpec New_ghost_spec(max_eta_track,1,0.01 );//Ghosts to calculate the Jet Area

        fastjet::AreaDefinition New_fAreaDef(fastjet::passive_area,New_ghost_spec);//Area Definition

        fastjet::ClusterSequenceArea New_clustSeq_Sig(fjInputs3,jetDefCh , New_fAreaDef);//Cluster Sequence
	NewJets[iR] = sorted_by_pt(New_clustSeq_Sig.inclusive_jets(1.));//Vector with the Reconstructed Jets
	if(NewJets[iR].size()==0) cout<<"Error: No NewJets found at "<<jetR<<endl;

//Analysis cut on the reconstructed jets
//		if(NewJets[iR][0].eta()>(max_eta_track-jetR)) continue;
//                if(NewJets[iR][0].eta()<-(max_eta_track-jetR)) continue;
//                foundjet=1;


//	cout<< "Before cleanup, number of particles: "<<fjInputs2.size()<<endl;
	
	for(unsigned int k=0;k<fjInputs2.size();k++){
	if (fjInputs2[k].user_index()==1){fjInputs2.erase(fjInputs2.begin()+k);k--;}

							}
//	 therm_sum = 0.0;
//	for(unsigned int j=0;j<fjInputs2.size();j++){therm_sum +=fjInputs2[j].pt();} 

//cout<< "After cleanup, number of particles: "<<fjInputs2.size()<<" with pt sum "<<therm_sum<<endl;



	}
//~~~~~~~~~~~~~~~~End of the Rjet loop
//cout<<"Checkpoint 5: End of R loop"<<endl;
//Saving in tree
if(foundjet){
R02.pt = NewJets[0][0].pt(); R02.eta = NewJets[0][0].eta(); R02.phi = NewJets[0][0].phi();
R04.pt = NewJets[1][0].pt(); R04.eta = NewJets[1][0].eta(); R04.phi = NewJets[1][0].phi();
R06.pt = NewJets[2][0].pt(); R06.eta = NewJets[2][0].eta(); R06.phi = NewJets[2][0].phi();
R08.pt = NewJets[3][0].pt(); R08.eta = NewJets[3][0].eta(); R08.phi = NewJets[3][0].phi();

jew_R02.pt = ProbeJet[0].pt(); jew_R02.eta = ProbeJet[0].eta(); jew_R02.phi = ProbeJet[0].phi();
jew_R04.pt = ProbeJet[1].pt(); jew_R04.eta = ProbeJet[1].eta(); jew_R04.phi = ProbeJet[1].phi();
jew_R06.pt = ProbeJet[2].pt(); jew_R06.eta = ProbeJet[2].eta(); jew_R06.phi = ProbeJet[2].phi();
jew_R08.pt = ProbeJet[3].pt(); jew_R08.eta = ProbeJet[3].eta(); jew_R08.phi = ProbeJet[3].phi();

Avg_pt_b[0] =avg_pt_before[0];Avg_pt_b[1] =avg_pt_before[1];Avg_pt_b[2] =avg_pt_before[2];Avg_pt_b[3] =avg_pt_before[3];
Avg_pt_a[0] =avg_pt_after[0];Avg_pt_a[1] =avg_pt_after[1];Avg_pt_a[2] =avg_pt_after[2];Avg_pt_a[3] =avg_pt_after[3]; 
Rho = rho;
Rho_gen = rho_gen;
tree->Fill();
}
cout<<"Cleanup at event end"<<endl;
for (int iR =0; iR < nR; iR++) {
      delete clustSeqCh[iR];
      clustSeqCh[iR]=0;
      delete clustSeqBG[iR];
      clustSeqBG[iR] = 0;
	}

delete evt;
ascii_in >> evt;

}
//~~~~~~~~~~~~~~~~End of the event loop
//cout<<"Checkpoint 6: End of event loop"<<endl;

fout.cd();
tree->Write();
fout.Close();
return 0;}
