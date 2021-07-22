void draw_histos()
{
  gStyle->SetOptStat(1111);
  gStyle->SetPalette(1,0);
  gStyle->SetOptTitle(0);
  TPaveLabel pl, pl1;
  TPaveText *pt;

  TFile *file1 = new TFile("pimC_incl_dst.root");

  TH1D *hist[4];
  TH1D *hist_4pi[4];

  file1->cd();
  hist[0]=(TH1D*)file1->Get("p_T");
  hist_4pi[0]=(TH1D*)file1->Get("p_T_4pi");
  hist[1]=(TH1D*)file1->Get("d_T");
  hist_4pi[1]=(TH1D*)file1->Get("d_T_4pi");
  hist[2]=(TH1D*)file1->Get("pip_T");
  hist_4pi[2]=(TH1D*)file1->Get("pip_T_4pi");
  hist[3]=(TH1D*)file1->Get("pim_T");
  hist_4pi[3]=(TH1D*)file1->Get("pim_T_4pi");

  TCanvas *canv[4];
  TLegend *leg[4];

  leg[0] = new TLegend(0.5687,0.5831,0.8782,0.7536);
  leg[1] = new TLegend(0.5687,0.5831,0.8782,0.7536);
  leg[2] = new TLegend(0.5687,0.5831,0.8782,0.7536);
  leg[3] = new TLegend(0.5687,0.5831,0.8782,0.7536);
  leg[0]->SetHeader("proton: T acc vs 4pi");
  leg[1]->SetHeader("deuteron: T acc vs 4pi");
  leg[2]->SetHeader("pip: T acc vs 4pi");
  leg[3]->SetHeader("pim: T acc vs 4pi");

  for(int i = 0; i < 4; i++)
  {
    canv[i] = new TCanvas();
    hist_4pi[i]->Draw();
    hist[i]->SetLineColor(kRed);
    hist[i]->Draw("SAME");

    leg[i]->AddEntry(hist_4pi[i],"4pi","l");
    leg[i]->AddEntry(hist[i],"acceptance","l");
    leg[i]->Draw("same");

  }

}
