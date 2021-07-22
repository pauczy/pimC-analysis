void draw_histos_2()
{
  gStyle->SetOptStat(1111);
  gStyle->SetPalette(1,0);
  gStyle->SetOptTitle(0);
  //TPaveLabel pl, pl1;
  //TPaveText *pt;

  TFile *file1 = new TFile("pimC_incl_dst.root");
  TFile *file2 = new TFile("pimC_incl_dst_M2.root");

  TH1D *hist[4];
  TH1D *hist_M2[4];

  file1->cd();
  hist[0]=(TH1D*)file1->Get("p_T");
  hist[1]=(TH1D*)file1->Get("d_T");
  hist[2]=(TH1D*)file1->Get("pip_T");
  hist[3]=(TH1D*)file1->Get("pim_T");

  file2->cd();
  hist_M2[0]=(TH1D*)file2->Get("p_T");
  hist_M2[1]=(TH1D*)file2->Get("d_T");
  hist_M2[2]=(TH1D*)file2->Get("pip_T");
  hist_M2[3]=(TH1D*)file2->Get("pim_T");


  TCanvas *canv[4];
  TLegend *leg[4];

  leg[0] = new TLegend(0.5687,0.5831,0.8782,0.7536);
  leg[1] = new TLegend(0.5687,0.5831,0.8782,0.7536);
  leg[2] = new TLegend(0.5687,0.5831,0.8782,0.7536);
  leg[3] = new TLegend(0.5687,0.5831,0.8782,0.7536);
  leg[0]->SetHeader("proton: T with M2 trig.");
  leg[1]->SetHeader("deuteron: T with M2 trig.");
  leg[2]->SetHeader("pip: T with M2 trig.");
  leg[3]->SetHeader("pim: T with M2 trig.");

  for(int i = 0; i < 4; i++)
  {
    canv[i] = new TCanvas();
    hist[i]->SetLineColor(kRed);
    hist[i]->Draw();
    hist_M2[i]->Draw("SAME");

    leg[i]->AddEntry(hist_M2[i],"M2 trig.","l");
    leg[i]->AddEntry(hist[i],"without trig.","l");
    leg[i]->Draw("same");

  }

  for(int i = 0; i < 4; i++)
  {
    cout<<i<<": "<<hist_M2[i]->Integral()/hist[i]->Integral()<<endl;
  }

}
