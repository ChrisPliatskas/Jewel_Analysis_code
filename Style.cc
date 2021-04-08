//Plotting style macro by Marco

void Style(){

//_______________________________________________________________General options______________________
gStyle->SetOptStat(0);
gStyle->SetLineWidth(2);
gStyle->SetFrameLineWidth(2);
gStyle->SetOptTitle(1);
gStyle->SetMarkerSize(1.5);


//_______________________________________________________________Title options________________________
gStyle->SetTitleFont(42);
gStyle->SetTitleSize(0.05);
gStyle->SetTitleSize(0.05,"Y");
gStyle->SetTitleOffset(0.9);

//_______________________________________________________________Pad options__________________________
gStyle->SetNdivisions(505);
gStyle->SetNdivisions(505,"Y");
gStyle->SetPadRightMargin(0.05);
gStyle->SetPadLeftMargin(0.12);
gStyle->SetPadTickX(1);
gStyle->SetPadTickY(1);

//_______________________________________________________________Text Options_________________________
gStyle->SetLabelSize(0.04);
gStyle->SetLabelSize(0.04,"Y");
gStyle->SetTextFont(42);
}
