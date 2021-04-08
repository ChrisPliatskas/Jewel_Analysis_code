void EFPlotter_jewel(const char* fname, int Ptlow, int Pthigh){
	
	//Open the file and retrieve the histograms
	TFile*f1 = TFile::Open(fname);
	
	//--Dpt vs Pt histograms
	TH2F* hDptR01_jew;
	TH2F* hDptR02_jew;
	TH2F* hDptR03_jew;

	//--DR vs Pt histograms
	TH2F* hDeltaR01_jew;
        TH2F* hDeltaR02_jew;
        TH2F* hDeltaR03_jew;

	//-- Pt spectra
	TH1F* hJetPt_R01;
	TH1F* hJetPt_R02;
	TH1F* hJetPt_R03;
	TH1F* hJetPt_R04;
	
	//-- Matched jet Pt spectra
	TH1F* hJetPtMatched_R01;
	TH1F* hJetPtMatched_R02;
	TH1F* hJetPtMatched_R03;
	TH1F* hJetPtMatched_R04;

	//-- Eta distributions
	TH1F* hJetEta_R01;
	TH1F* hJetEta_R02;
	TH1F* hJetEta_R03;
	TH1F* hJetEta_R04;

	//--Matched jet eta distributions
	TH1F* hJetEtaMatched_R01;
        TH1F* hJetEtaMatched_R02;
        TH1F* hJetEtaMatched_R03;
        TH1F* hJetEtaMatched_R04;

	//-- Pt vs DR vs Dpt
	TH3F*hPtJetDeltaRDeltaPt_R01;
	TH3F*hPtJetDeltaRDeltaPt_R02;
	TH3F*hPtJetDeltaRDeltaPt_R03;

	//-- Eta vs DR
	TH2F*hEtaJetDeltaR_R01;
	TH2F*hEtaJetDeltaR_R02;
	TH2F*hEtaJetDeltaR_R03;
	
	//-- Dpt vs Pt vs DEta
	TH3F*hDptPtDEta_R01;
	TH3F*hDptPtDEta_R02;
	TH3F*hDptPtDEta_R03;

	//--Dpt vs Pt vs Multiplicity
	TH3F*hDptPtMultiplicity_R01;
	TH3F*hDptPtMultiplicity_R02;
	TH3F*hDptPtMultiplicity_R03;

	//--Acquiring the pt and eta spectra of the jets
	f1->GetObject("hJetPt_R01",hJetPt_R01);
	f1->GetObject("hJetPt_R02",hJetPt_R02);
	f1->GetObject("hJetPt_R03",hJetPt_R03);
	f1->GetObject("hJetPt_R04",hJetPt_R04);

	f1->GetObject("hJetEta_R01",hJetEta_R01);
        f1->GetObject("hJetEta_R02",hJetEta_R02);
        f1->GetObject("hJetEta_R03",hJetEta_R03);
        f1->GetObject("hJetEta_R04",hJetEta_R04);
	
	f1->GetObject("hJetPtMatched_R01",hJetPtMatched_R01);
        f1->GetObject("hJetPtMatched_R02",hJetPtMatched_R02);
        f1->GetObject("hJetPtMatched_R03",hJetPtMatched_R03);
        f1->GetObject("hJetPtMatched_R04",hJetPtMatched_R04);

	f1->GetObject("hJetEtaMatched_R01",hJetEtaMatched_R01);
        f1->GetObject("hJetEtaMatched_R02",hJetEtaMatched_R02);
        f1->GetObject("hJetEtaMatched_R03",hJetEtaMatched_R03);
        f1->GetObject("hJetEtaMatched_R04",hJetEtaMatched_R04);

	//Acquiring the Dpt and DR histograms from the file
	f1->GetObject("hPtJetDeltaPt_R01",hDptR01_jew);
	f1->GetObject("hPtJetDeltaPt_R02",hDptR02_jew);
	f1->GetObject("hPtJetDeltaPt_R03",hDptR03_jew);

        f1->GetObject("hPtJetDeltaR_R01",hDeltaR01_jew);
        f1->GetObject("hPtJetDeltaR_R02",hDeltaR02_jew);
        f1->GetObject("hPtJetDeltaR_R03",hDeltaR03_jew);

	f1->GetObject("hPtJetDeltaRDeltaPt_R01",hPtJetDeltaRDeltaPt_R01);
        f1->GetObject("hPtJetDeltaRDeltaPt_R02",hPtJetDeltaRDeltaPt_R02);
        f1->GetObject("hPtJetDeltaRDeltaPt_R03",hPtJetDeltaRDeltaPt_R03);

        f1->GetObject("hDptPtDEta_R01",hDptPtDEta_R01);
        f1->GetObject("hDptPtDEta_R02",hDptPtDEta_R02);
        f1->GetObject("hDptPtDEta_R03",hDptPtDEta_R03);

	f1->GetObject("hEtaJetDeltaR_R01",hEtaJetDeltaR_R01);
	f1->GetObject("hEtaJetDeltaR_R02",hEtaJetDeltaR_R02);
	f1->GetObject("hEtaJetDeltaR_R03",hEtaJetDeltaR_R03);

	//Acquiring the multiplicity histograms from the file
	f1->GetObject("hDptPtMultiplicity_R01",hDptPtMultiplicity_R01);
	f1->GetObject("hDptPtMultiplicity_R02",hDptPtMultiplicity_R02);
	f1->GetObject("hDptPtMultiplicity_R03",hDptPtMultiplicity_R03);



	//Calculating the Dpt projections in the pt interval of our interest
	TH1D*pro_R01_jew; TH1D*pro_R02_jew; TH1D*pro_R03_jew;
		
	
        pro_R01_jew = hDptR01_jew->ProjectionY("pro_R01_jew",Ptlow,Pthigh);
        pro_R02_jew = hDptR02_jew->ProjectionY("pro_R02_jew",Ptlow,Pthigh);
        pro_R03_jew = hDptR03_jew->ProjectionY("pro_R03_jew",Ptlow,Pthigh);


	TH1D*pro_R01_jew_dR; TH1D*pro_R02_jew_dR; TH1D*pro_R03_jew_dR;


        pro_R01_jew_dR = hDeltaR01_jew->ProjectionY("pro_R01_jew_dR",Ptlow,Pthigh);
        pro_R02_jew_dR = hDeltaR02_jew->ProjectionY("pro_R02_jew_dR",Ptlow,Pthigh);
        pro_R03_jew_dR = hDeltaR03_jew->ProjectionY("pro_R03_jew_dR",Ptlow,Pthigh);
        

	hPtJetDeltaRDeltaPt_R01->GetXaxis()->SetRange(Ptlow,Pthigh);
	TH1* pro_R01_jew_deltas=hPtJetDeltaRDeltaPt_R01->Project3D("zy");
        hPtJetDeltaRDeltaPt_R02->GetXaxis()->SetRange(Ptlow,Pthigh);
        TH1* pro_R02_jew_deltas=hPtJetDeltaRDeltaPt_R02->Project3D("zy");
        hPtJetDeltaRDeltaPt_R03->GetXaxis()->SetRange(Ptlow,Pthigh);
        TH1* pro_R03_jew_deltas=hPtJetDeltaRDeltaPt_R03->Project3D("zy");

	hDptPtDEta_R01->GetYaxis()->SetRange(Ptlow,Pthigh);
	TH1* pro_R01_jew_Deta=hDptPtDEta_R01->Project3D("zx");
        hDptPtDEta_R02->GetYaxis()->SetRange(Ptlow,Pthigh);
        TH1* pro_R02_jew_Deta=hDptPtDEta_R02->Project3D("zx");
        hDptPtDEta_R03->GetYaxis()->SetRange(Ptlow,Pthigh);
        TH1* pro_R03_jew_Deta=hDptPtDEta_R03->Project3D("zx");

	//Calculation the Dpt vs Multiplicity projections and plotting them
	hDptPtMultiplicity_R01->GetYaxis()->SetRange(Ptlow,Pthigh);
	TH1* pro_R01_jew_mult = hDptPtMultiplicity_R01->Project3D("zx");
	hDptPtMultiplicity_R02->GetYaxis()->SetRange(Ptlow,Pthigh);
        TH1* pro_R02_jew_mult = hDptPtMultiplicity_R02->Project3D("zx");
	hDptPtMultiplicity_R03->GetYaxis()->SetRange(Ptlow,Pthigh);
        TH1* pro_R03_jew_mult = hDptPtMultiplicity_R03->Project3D("zx");

	TCanvas *mult1 = new TCanvas("mult1","Dpt vs Multiplicity R01",1000,1000);
	pro_R01_jew_mult->SetTitle(TString::Format("Multiplicity vs #DeltaP_{t} projection between R=0.1 and R=0.2 [%d,%d]GeV/c (JEWEL)",Ptlow,Pthigh));
	pro_R01_jew_mult->Draw("colz");

	TCanvas *mult2 = new TCanvas("mult2","Dpt vs Multiplicity R02",1000,1000);
        pro_R02_jew_mult->SetTitle(TString::Format("Multiplicity vs #DeltaP_{t} between R=0.2 and R=0.3 [%d,%d]GeV/c (JEWEL)",Ptlow,Pthigh));
        pro_R02_jew_mult->Draw("colz");

	TCanvas *mult3 = new TCanvas("mult3","Dpt vs Multiplicity R03",1000,1000);
        pro_R03_jew_mult->SetTitle(TString::Format("Multiplicity vs #DeltaP_{t} between R=0.3 and R=0.4 [%d,%d]GeV/c (JEWEL)",Ptlow,Pthigh));
        pro_R03_jew_mult->Draw("colz");

//------------------------Histogram manipulation will happen here---------------------

	TCanvas *Delta1 = new TCanvas("Delta1","Dpt vs DR R01",1000,1000);
        pro_R01_jew_deltas->SetTitle(TString::Format("#DeltaP_{t} vs #DeltaR between R=0.1 and R=0.2 [%d,%d]GeV/c (JEWEL)",Ptlow,Pthigh));
        pro_R01_jew_deltas->Draw("colz");

	TCanvas *Delta2 = new TCanvas("Delta2","Dpt vs DR R02",1000,1000);
        pro_R02_jew_deltas->SetTitle(TString::Format("#DeltaP_{t} vs #DeltaR between R=0.2 and R=0.3 [%d,%d]GeV/c (JEWEL)",Ptlow,Pthigh));
        pro_R02_jew_deltas->Draw("colz");

	TCanvas *Delta3 = new TCanvas("Delta3","Dpt vs DR R03",1000,1000);
        pro_R03_jew_deltas->SetTitle(TString::Format("#DeltaP_{t} vs #DeltaR between R=0.3 and R=0.4 [%d,%d]GeV/c (JEWEL)",Ptlow,Pthigh));
        pro_R03_jew_deltas->Draw("colz");


	TCanvas *DReta1 = new TCanvas("DReta1","DR vs eta R01",1000,1000);
	hEtaJetDeltaR_R01->Draw("colz");

	TCanvas *DReta2 = new TCanvas("DReta2","DR vs eta R02",1000,1000);
        hEtaJetDeltaR_R02->Draw("colz");

	TCanvas *DReta3 = new TCanvas("DReta3","DR vs eta R03",1000,1000);
        hEtaJetDeltaR_R03->Draw("colz");

	TCanvas *DEta1 = new TCanvas("DEta1","Dpt vs DEta R01",1000,1000);
        pro_R01_jew_Deta->SetTitle(TString::Format("#DeltaP_{t} vs #Delta#eta between R=0.1 and R=0.2 [%d,%d]GeV/c (JEWEL)",Ptlow,Pthigh));
        pro_R01_jew_Deta->Draw("colz");

        TCanvas *DEta2 = new TCanvas("DEta2","Dpt vs DEta R02",1000,1000);
        pro_R02_jew_Deta->SetTitle(TString::Format("#DeltaP_{t} vs #Delta#eta between R=0.2 and R=0.3 [%d,%d]GeV/c (JEWEL)",Ptlow,Pthigh));
        pro_R02_jew_Deta->Draw("colz");

        TCanvas *DEta3 = new TCanvas("DEta3","Dpt vs DEta R03",1000,1000);
        pro_R03_jew_Deta->SetTitle(TString::Format("#DeltaP_{t} vs #Delta#eta between R=0.3 and R=0.4 [%d,%d]GeV/c (JEWEL)",Ptlow,Pthigh));
        pro_R03_jew_Deta->Draw("colz");


	//----------------------Dpt Projections---------------------------------------
	pro_R01_jew->SetTitle(TString::Format("#DeltaP_{t} projection between R=0.1 and R=0.2 [%d,%d]GeV/c (JEWEL)",Ptlow,Pthigh));
	pro_R02_jew->SetTitle(TString::Format("#DeltaP_{t} projection between R=0.2 and R=0.3 [%d,%d]GeV/c (JEWEL)",Ptlow,Pthigh));
	pro_R03_jew->SetTitle(TString::Format("#DeltaP_{t} projection between R=0.3 and R=0.4 [%d,%d]GeV/c (JEWEL)",Ptlow,Pthigh));

	//----------------------DR Projections-----------------------------------------
        pro_R01_jew_dR->SetTitle(TString::Format("#DeltaR projection between R=0.1 and R=0.2 [%d,%d]GeV/c (JEWEL)",Ptlow,Pthigh));
        pro_R02_jew_dR->SetTitle(TString::Format("#DeltaR projection between R=0.2 and R=0.3 [%d,%d]GeV/c (JEWEL)",Ptlow,Pthigh));
        pro_R03_jew_dR->SetTitle(TString::Format("#DeltaR projection between R=0.3 and R=0.4 [%d,%d]GeV/c (JEWEL)",Ptlow,Pthigh));	

	//General plot style settings
	gStyle->SetOptStat(0);
        gStyle->SetTitleFont(42);
        gStyle->SetTitleSize(0.05);
        gStyle->SetTitleSize(0.05,"Y");
        gStyle->SetTitleOffset(0.9,"xy");

	//Adjusting colors and marker style/sizes	
	pro_R01_jew->SetLineColor(kBlue+2);
	pro_R01_jew->SetLineWidth(2);
	pro_R01_jew_dR->SetLineColor(kBlue+2);
	pro_R01_jew_dR->SetLineWidth(2);	

	pro_R02_jew->SetLineColor(2);
        pro_R02_jew->SetLineWidth(2);
	pro_R02_jew_dR->SetLineColor(2);
        pro_R02_jew_dR->SetLineWidth(2);
		
	pro_R03_jew->SetLineColor(kGreen+2);
        pro_R03_jew->SetLineWidth(2);
	pro_R03_jew_dR->SetLineColor(kGreen+2);
        pro_R03_jew_dR->SetLineWidth(2);	


	pro_R01_jew->SetMarkerSize(3);
	pro_R01_jew->SetMarkerStyle(87);
	pro_R01_jew->SetMarkerColor(kBlue+2);
	
	pro_R01_jew_dR->SetMarkerSize(3);
        pro_R01_jew_dR->SetMarkerStyle(87);
        pro_R01_jew_dR->SetMarkerColor(kBlue+2);

	pro_R02_jew->SetMarkerSize(3);
	pro_R02_jew->SetMarkerStyle(87);
        pro_R02_jew->SetMarkerColor(2);
	
	pro_R02_jew_dR->SetMarkerSize(3);
        pro_R02_jew_dR->SetMarkerStyle(87);
        pro_R02_jew_dR->SetMarkerColor(2);

	pro_R03_jew->SetMarkerSize(3);
        pro_R03_jew->SetMarkerStyle(87);
	pro_R03_jew->SetMarkerColor(kGreen+2);

        pro_R03_jew_dR->SetMarkerSize(3);
        pro_R03_jew_dR->SetMarkerStyle(87);
        pro_R03_jew_dR->SetMarkerColor(kGreen+2);	

	
	//Jet Pt spectrum
	hJetPt_R01->SetLineColor(kBlue+2);
        hJetPt_R01->SetLineWidth(2);
	hJetPt_R01->SetMarkerSize(2);
        hJetPt_R01->SetMarkerStyle(87);
        hJetPt_R01->SetMarkerColor(kBlue+2);
	
	hJetPt_R02->SetLineColor(2);
        hJetPt_R02->SetLineWidth(2);
	hJetPt_R02->SetMarkerSize(2);
        hJetPt_R02->SetMarkerStyle(87);
        hJetPt_R02->SetMarkerColor(2);
	
	hJetPt_R03->SetLineColor(kGreen+2);
        hJetPt_R03->SetLineWidth(2);
	hJetPt_R03->SetMarkerSize(2);
        hJetPt_R03->SetMarkerStyle(87);
        hJetPt_R03->SetMarkerColor(kGreen+2);

	hJetPt_R04->SetLineColor(kMagenta+2);
        hJetPt_R04->SetLineWidth(2);
        hJetPt_R04->SetMarkerSize(2);
        hJetPt_R04->SetMarkerStyle(87);
        hJetPt_R04->SetMarkerColor(kMagenta+2);

	//Jet Eta distribution
	hJetEta_R01->SetLineColor(kBlue+2);
        hJetEta_R01->SetLineWidth(2);
        hJetEta_R01->SetMarkerSize(3);
        hJetEta_R01->SetMarkerStyle(87);
        hJetEta_R01->SetMarkerColor(kBlue+2);

        hJetEta_R02->SetLineColor(2);
        hJetEta_R02->SetLineWidth(2);
        hJetEta_R02->SetMarkerSize(3);
        hJetEta_R02->SetMarkerStyle(87);
        hJetEta_R02->SetMarkerColor(2);

        hJetEta_R03->SetLineColor(kGreen+2);
        hJetEta_R03->SetLineWidth(2);
        hJetEta_R03->SetMarkerSize(3);
        hJetEta_R03->SetMarkerStyle(87);
        hJetEta_R03->SetMarkerColor(kGreen+2);

        hJetEta_R04->SetLineColor(kMagenta+2);
        hJetEta_R04->SetLineWidth(2);
        hJetEta_R04->SetMarkerSize(3);
        hJetEta_R04->SetMarkerStyle(87);
        hJetEta_R04->SetMarkerColor(kMagenta+2);

	//Matched Jet Pt spectrum
	hJetPtMatched_R01->SetLineColor(kBlue+2);
        hJetPtMatched_R01->SetLineWidth(2);
        hJetPtMatched_R01->SetMarkerSize(2);
        hJetPtMatched_R01->SetMarkerStyle(87);
        hJetPtMatched_R01->SetMarkerColor(kBlue+2);

        hJetPtMatched_R02->SetLineColor(2);
        hJetPtMatched_R02->SetLineWidth(2);
        hJetPtMatched_R02->SetMarkerSize(2);
        hJetPtMatched_R02->SetMarkerStyle(87);
        hJetPtMatched_R02->SetMarkerColor(2);

        hJetPtMatched_R03->SetLineColor(kGreen+2);
        hJetPtMatched_R03->SetLineWidth(2);
        hJetPtMatched_R03->SetMarkerSize(2);
        hJetPtMatched_R03->SetMarkerStyle(87);
        hJetPtMatched_R03->SetMarkerColor(kGreen+2);

        hJetPtMatched_R04->SetLineColor(kMagenta+2);
        hJetPtMatched_R04->SetLineWidth(2);
        hJetPtMatched_R04->SetMarkerSize(2);
        hJetPtMatched_R04->SetMarkerStyle(87);
        hJetPtMatched_R04->SetMarkerColor(kMagenta+2);
	
        hJetPtMatched_R01->SetLineColor(kBlue+2);
        hJetPtMatched_R01->SetLineWidth(2);
        hJetPtMatched_R01->SetMarkerSize(2);
        hJetPtMatched_R01->SetMarkerStyle(87);
        hJetPtMatched_R01->SetMarkerColor(kBlue+2);

        hJetPtMatched_R02->SetLineColor(2);
        hJetPtMatched_R02->SetLineWidth(2);
        hJetPtMatched_R02->SetMarkerSize(2);
        hJetPtMatched_R02->SetMarkerStyle(87);
        hJetPtMatched_R02->SetMarkerColor(2);

        hJetPtMatched_R03->SetLineColor(kGreen+2);
        hJetPtMatched_R03->SetLineWidth(2);
        hJetPtMatched_R03->SetMarkerSize(2);
        hJetPtMatched_R03->SetMarkerStyle(87);
        hJetPtMatched_R03->SetMarkerColor(kGreen+2);

        hJetPtMatched_R04->SetLineColor(kMagenta+2);
        hJetPtMatched_R04->SetLineWidth(2);
        hJetPtMatched_R04->SetMarkerSize(2);
        hJetPtMatched_R04->SetMarkerStyle(87);
        hJetPtMatched_R04->SetMarkerColor(kMagenta+2);

	//Matched Jet Eta spectrum
	hJetEtaMatched_R01->SetLineColor(kBlue+2);
        hJetEtaMatched_R01->SetLineWidth(2);
        hJetEtaMatched_R01->SetMarkerSize(3);
        hJetEtaMatched_R01->SetMarkerStyle(87);
        hJetEtaMatched_R01->SetMarkerColor(kBlue+2);

        hJetEtaMatched_R02->SetLineColor(2);
        hJetEtaMatched_R02->SetLineWidth(2);
        hJetEtaMatched_R02->SetMarkerSize(3);
        hJetEtaMatched_R02->SetMarkerStyle(87);
        hJetEtaMatched_R02->SetMarkerColor(2);

        hJetEtaMatched_R03->SetLineColor(kGreen+2);
        hJetEtaMatched_R03->SetLineWidth(2);
        hJetEtaMatched_R03->SetMarkerSize(3);
        hJetEtaMatched_R03->SetMarkerStyle(87);
        hJetEtaMatched_R03->SetMarkerColor(kGreen+2);

        hJetEtaMatched_R04->SetLineColor(kMagenta+2);
        hJetEtaMatched_R04->SetLineWidth(2);
        hJetEtaMatched_R04->SetMarkerSize(3);
        hJetEtaMatched_R04->SetMarkerStyle(87);
        hJetEtaMatched_R04->SetMarkerColor(kMagenta+2);

	
	//Histogram normalization
	
	pro_R01_jew->Scale(1/pro_R01_jew->Integral());
	pro_R02_jew->Scale(1/pro_R02_jew->Integral());
	pro_R03_jew->Scale(1/pro_R03_jew->Integral());

	pro_R01_jew_dR->Scale(1/pro_R01_jew_dR->Integral());
        pro_R02_jew_dR->Scale(1/pro_R02_jew_dR->Integral());
        pro_R03_jew_dR->Scale(1/pro_R03_jew_dR->Integral());




	//Mean and RMS calculation
	Double_t R[3] = {.1,.2,.3};
        Double_t mean[3] ={pro_R01_jew->GetMean(),pro_R02_jew->GetMean(),pro_R03_jew->GetMean()};
        Double_t RMS[3] = {pro_R01_jew->GetRMS(),pro_R02_jew->GetRMS(),pro_R03_jew->GetRMS()};
        Double_t mean_err[3] = {pro_R01_jew->GetMeanError(),pro_R02_jew->GetMeanError(),pro_R03_jew->GetMeanError()};
        Double_t RMS_err[3] = {pro_R01_jew->GetRMSError(),pro_R02_jew->GetRMSError(),pro_R03_jew->GetRMSError()};
        Double_t R_err[3] = {0,0,0};

        auto Mean_gr = new TGraphErrors(3,R,mean,R_err,mean_err);
        Mean_gr->SetMarkerSize(3);
        Mean_gr->SetMarkerStyle(33);
	Mean_gr->SetTitle(TString::Format("Mean of the #DeltaP_{t} projection distributions [%d,%d](JEWEL);R_{low}; Mean value",Ptlow,Pthigh));
//      Mean_gr->SetAxisRange(0,3,"Y");

        auto RMS_gr = new TGraphErrors(3,R,RMS,R_err,RMS_err);
        RMS_gr->SetMarkerSize(3);
        RMS_gr->SetMarkerStyle(33);
        RMS_gr->SetTitle(TString::Format("RMS of the #DeltaP_{t} projection distributions [%d,%d](JEWEL);R_{low}; RMS value",Ptlow,Pthigh));


	/* Old implementation for the Mean and RMS plot as a histogram
	TH2D*h_mean_pp_jew = new TH2D("h_mean_pp_jew",TString::Format("Mean of the #DeltaP_{t} projection distributions [%d,%d]GeV/c (JEWEL);R_{low}; Mean value",Ptlow,Pthigh),100,0,.8,100,0,10);
	
	TH2D*h_rms_pp_jew = new TH2D("h_rms_pp_jew",TString::Format("RMS of the #DeltaP_{t} projection distributions [%d,%d]GeV/c (JEWEL);R_{low}; RMS value",Ptlow,Pthigh),100,0,.8,100,0,10);

	h_mean_pp_jew->Fill(.2,pro_R02_jew->GetMean());
	h_mean_pp_jew->Fill(.4,pro_R04_jew->GetMean());
	h_mean_pp_jew->Fill(.6,pro_R06_jew->GetMean());

	h_mean_pp_jew->SetMarkerSize(2);
	h_mean_pp_jew->SetMarkerStyle(33);

	h_rms_pp_jew->Fill(.2,pro_R02_jew->GetRMS());
        h_rms_pp_jew->Fill(.4,pro_R04_jew->GetRMS());
        h_rms_pp_jew->Fill(.6,pro_R06_jew->GetRMS());

        h_rms_pp_jew->SetMarkerSize(2);
        h_rms_pp_jew->SetMarkerStyle(33);
*/
//Plotting section
	TCanvas *c1 = new TCanvas("c1","Jet Pt spectra",1000,1000);
        c1->SetLogy(1);
	hJetPt_R01->SetMaximum(1.2);
	hJetPt_R01->Draw();
	hJetPt_R02->Draw("same");
	hJetPt_R03->Draw("same");
	hJetPt_R04->Draw("same");

	auto legend = new TLegend(.7,.8,.9,.9);
	legend->AddEntry(hJetPt_R01,"R = 0.1","l");
        legend->AddEntry(hJetPt_R02,"R = 0.2","l");
        legend->AddEntry(hJetPt_R03,"R = 0.3","l");
	legend->AddEntry(hJetPt_R04,"R = 0.4","l");
        legend->Draw();

        TCanvas *c1a = new TCanvas("c1a","Jet Eta distribution",1000,1000);
        hJetEta_R01->Draw();
        hJetEta_R02->Draw("same");
        hJetEta_R03->Draw("same");
        hJetEta_R04->Draw("same");
	legend->Draw();

        TCanvas *c1b = new TCanvas("c1b","Matched jet Pt spectra",1000,1000);
        c1b->SetLogy(1);
        hJetPtMatched_R01->SetMaximum(1.2);
        hJetPtMatched_R01->Draw();
        hJetPtMatched_R02->Draw("same");
        hJetPtMatched_R03->Draw("same");
        hJetPtMatched_R04->Draw("same");
	legend->Draw();

	TCanvas *c1c = new TCanvas("c1c","Matched jet Eta distribution",1000,1000);
        hJetEtaMatched_R01->Draw();
        hJetEtaMatched_R02->Draw("same");
        hJetEtaMatched_R03->Draw("same");
        hJetEtaMatched_R04->Draw("same");
        legend->Draw();

	TCanvas*c2 = new TCanvas("c2","DeltaPt projections overview",1000,1000);
        c2->SetLogy(1);
        pro_R01_jew->SetMaximum(.8);
        pro_R01_jew->SetAxisRange(-20,30,"X");
        pro_R01_jew->Draw();
        pro_R02_jew->Draw("same");
        pro_R03_jew->Draw("same");

        auto legend1 = new TLegend(.7,.8,.9,.9);
        legend1->AddEntry(pro_R01_jew,"R = 0.1 to R = 0.2","l");
        legend1->AddEntry(pro_R02_jew,"R = 0.2 to R = 0.3","l");
        legend1->AddEntry(pro_R03_jew,"R = 0.3 to R = 0.4","l");
        legend1->Draw();

        TCanvas*c3 = new TCanvas("c3","DeltaR overview",1000,1000);
        pro_R01_jew_dR->SetMaximum(.5);
	pro_R01_jew_dR->Draw();
	pro_R02_jew_dR->Draw("same");
	pro_R03_jew_dR->Draw("same");
	legend1->Draw();

	//Mean and RMS plots
	TCanvas*c4 = new TCanvas("c4","Dpt projection Mean",1000,1000);
	c4->DrawFrame(0,0,.5,12);
        Mean_gr->Draw("ALP");

        TCanvas*c4a = new TCanvas("c4a","Dpt projection RMS",1000,1000);
        c4a->DrawFrame(0,0,.5,10);
	RMS_gr->Draw("ALP");



//------------------------Writing the output to a new file
	TFile*sum = new TFile(TString::Format("EFFluctuations_summary_%d_%djewel_pp.root",Ptlow,Pthigh),"RECREATE",TString::Format("EFFluctuations_summary_%d_%djewel_pp.root",Ptlow,Pthigh));
	sum->cd();

//	h_rms_pp_jew->Write();
//	h_mean_pp_jew->Write();

	RMS_gr->Write(TString::Format("RMS_graph_%d_%d",Ptlow,Pthigh));
        Mean_gr->Write(TString::Format("Mean_graph_%d_%d",Ptlow,Pthigh));
	
	pro_R01_jew_dR->Write();
        pro_R02_jew_dR->Write();
        pro_R03_jew_dR->Write();
	
	pro_R01_jew->Write(); pro_R02_jew->Write(); pro_R03_jew->Write();

	pro_R01_jew_deltas->Write(); pro_R02_jew_deltas->Write(); pro_R03_jew_deltas->Write();

        pro_R01_jew_mult->Write(); pro_R02_jew_mult->Write(); pro_R03_jew_mult->Write();

	pro_R01_jew_Deta->Write(); pro_R02_jew_Deta->Write(); pro_R03_jew_Deta->Write();

	hJetPt_R01->Write(); hJetPt_R02->Write(); hJetPt_R03->Write(); hJetPt_R04->Write();
	hJetEta_R01->Write(); hJetEta_R02->Write(); hJetEta_R03->Write(); hJetEta_R04->Write();
	hJetPtMatched_R01->Write(); hJetPtMatched_R02->Write(); 
	hJetPtMatched_R03->Write(); hJetPtMatched_R04->Write();
	hJetEtaMatched_R01->Write(); hJetEtaMatched_R02->Write();
        hJetEtaMatched_R03->Write(); hJetEtaMatched_R04->Write();

	hEtaJetDeltaR_R01->Write();hEtaJetDeltaR_R02->Write();hEtaJetDeltaR_R03->Write();

	sum->Write();
		//----Adjustment of the plot titles
        hJetPt_R01->SetTitle("Jet p_{t} spectra overview (JEWEL)");
	hJetEta_R01->SetTitle("Jet #eta distribution overview (JEWEL)");
	hJetPtMatched_R01->SetTitle("Matched jet p_{t} spectra overview (JEWEL)");
	hJetEtaMatched_R01->SetTitle("Matched jet #eta distribution overview (JEWEL)");
	
	
	
	pro_R01_jew->SetTitle(TString::Format("#DeltaP_{t} projection overview [%d,%d](GeV/c) (JEWEL)",Ptlow,Pthigh));
        pro_R01_jew_dR->SetTitle(TString::Format("#DeltaR projection overview [%d,%d](GeV/c) (JEWEL)",Ptlow,Pthigh));
	sum->Close();
}
