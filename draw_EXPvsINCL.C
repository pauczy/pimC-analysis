void draw_EXPvsINCL()
{
  gStyle->SetOptStat(1111);
  gStyle->SetPalette(1,0);
  gStyle->SetOptTitle(0);

  TFile *file1 = new TFile("pimC_incl_full3.root");
  TFile *file2 = new TFile("sum_pimC_690.root");

  //histos in acceptance
  TH1F *hist_incl[4];
  TH1F *hist_exp[4];

  //h*_enth_40, *_pt, *_y, *_thLab (theta), *_momLab (p)
  file1->cd();
  hist_incl[0]=(TH1F*)file1->Get("p_theta");
  hist_incl[1]=(TH1F*)file1->Get("d_theta");
  hist_incl[2]=(TH1F*)file1->Get("pip_theta");
  hist_incl[3]=(TH1F*)file1->Get("pim_theta");
  //całkowity przekroj czynny [mb] podzielony przez ilosc eventow w symulacji:1462.32 mb/100 000 000
  hist_incl[0]->Scale(1.462*TMath::Power(10, -5));
  hist_incl[1]->Scale(1.462*TMath::Power(10, -5));
  hist_incl[2]->Scale(1.462*TMath::Power(10, -5));
  hist_incl[3]->Scale(1.462*TMath::Power(10, -5));
  // hist_incl[0]->Scale(1./hist_incl[0]->Integral());
  // hist_incl[1]->Scale(1./hist_incl[1]->Integral());
  // hist_incl[2]->Scale(1./hist_incl[2]->Integral());
  // hist_incl[3]->Scale(1./hist_incl[3]->Integral());
  hist_incl[0]->Sumw2(kFALSE);
  hist_incl[1]->Sumw2(kFALSE);
  hist_incl[2]->Sumw2(kFALSE);
  hist_incl[3]->Sumw2(kFALSE);

  file2->cd();
  hist_exp[0]=(TH1F*)file2->Get("p_thLab");
  hist_exp[1]=(TH1F*)file2->Get("d_thLab");
  hist_exp[2]=(TH1F*)file2->Get("pip_thLab");
  hist_exp[3]=(TH1F*)file2->Get("pim_thLab");
  //czynnik wyznaczony z elastycznego rozpraszania pim-p;
  hist_exp[0]->Scale(10.165*TMath::Power(10, -7));
  hist_exp[1]->Scale(10.165*TMath::Power(10, -7));
  hist_exp[2]->Scale(10.165*TMath::Power(10, -7));
  hist_exp[3]->Scale(10.165*TMath::Power(10, -7));
  // hist_exp[0]->Scale(1./hist_exp[0]->Integral());
  // hist_exp[1]->Scale(1./hist_exp[1]->Integral());
  // hist_exp[2]->Scale(1./hist_exp[2]->Integral());
  // hist_exp[3]->Scale(1./hist_exp[3]->Integral());
  hist_exp[0]->Sumw2(kFALSE);
  hist_exp[1]->Sumw2(kFALSE);
  hist_exp[2]->Sumw2(kFALSE);
  hist_exp[3]->Sumw2(kFALSE);

  TCanvas *canv[4];
  TLegend *leg[4];

  leg[0] = new TLegend(0.6691,0.6126,0.8682,0.7389);
  leg[0]->SetBorderSize(0);
  leg[1] = new TLegend(0.6691,0.6126,0.8682,0.7389);
  leg[1]->SetBorderSize(0);
  leg[2] = new TLegend(0.6691,0.6126,0.8682,0.7389);
  leg[2]->SetBorderSize(0);
  leg[3] = new TLegend(0.6691,0.6126,0.8682,0.7389);
  leg[3]->SetBorderSize(0);
  // leg[0]->SetHeader("p: pt experiment vs INCL (acceptance)");
  // leg[1]->SetHeader("d: pt experiment vs INCL (acceptance)");
  // leg[2]->SetHeader("pip: pt experiment vs INCL (acceptance)");
  // leg[3]->SetHeader("pim: pt experiment vs INCL (acceptance)");

  for(int i = 0; i < 4; i++)
  {
    canv[i] = new TCanvas();
    double maxy = max(hist_exp[i]->GetMaximum(), hist_incl[i]->GetMaximum());
    hist_incl[i]->GetYaxis()->SetRangeUser(0., maxy*1.1);
    hist_incl[i]->SetLineColor(kRed);
    hist_incl[i]->Draw();
    hist_exp[i]->Draw("SAME");

    leg[i]->AddEntry(hist_incl[i],"INCL","l");
    leg[i]->AddEntry(hist_exp[i],"data","l");
    leg[i]->Draw("same");
  }

}
