void draw_histos_2()
{
  gStyle->SetOptStat(1111);
  gStyle->SetPalette(1,0);
  gStyle->SetOptTitle(0);
  //TPaveLabel pl, pl1;
  //TPaveText *pt;

  TFile *file1 = new TFile("pimC_incl_full2.root");
  //TFile *file2 = new TFile("pimC_incl_dst_M2.root");

  TH1D *hist[4];
  TH1D *hist_4pi[4];
  TH1D *hist_ratio[4];

  file1->cd();
  hist[0]=(TH1D*)file1->Get("p_T");
  hist[1]=(TH1D*)file1->Get("d_T");
  hist[2]=(TH1D*)file1->Get("pip_T");
  hist[3]=(TH1D*)file1->Get("pim_T");

  hist_4pi[0]=(TH1D*)file1->Get("p_T_4pi");
  hist_4pi[1]=(TH1D*)file1->Get("d_T_4pi");
  hist_4pi[2]=(TH1D*)file1->Get("pip_T_4pi");
  hist_4pi[3]=(TH1D*)file1->Get("pim_T_4pi");


  for(int i =0; i < 4; i++)
  {
    hist_ratio[i]=(TH1D*)hist[i]->Clone();
    hist_ratio[i]->Divide(hist_4pi[i]);

    cout<<i<<" "<<hist[i]->GetNbinsX()<<" "<<hist_4pi[i]->GetNbinsX()<<" "<<hist_ratio[i]->GetNbinsX()<<endl;
    cout<<i<<" acc: "<<hist[i]->GetBinContent(35)<<" 4pi: "<<hist_4pi[i]->GetBinContent(35)<<endl;
    cout<<"div: "<<hist[i]->GetBinContent(35)/hist_4pi[i]->GetBinContent(35)<<" ratio: "<<hist_ratio[i]->GetBinContent(35)<<endl;
  }


  TCanvas *canv[4];
  TLegend *leg[4];

  leg[0] = new TLegend(0.6848,0.6611,0.8797,0.7516);
  leg[1] = new TLegend(0.6848,0.6611,0.8797,0.7516);
  leg[2] = new TLegend(0.6848,0.6611,0.8797,0.7516);
  leg[3] = new TLegend(0.6848,0.6611,0.8797,0.7516);
  leg[0]->SetHeader("proton: T acc vs 4pi");
  leg[1]->SetHeader("deuteron: T acc vs 4pi");
  leg[2]->SetHeader("pip: T acc vs 4pi");
  leg[3]->SetHeader("pim: T acc vs 4pi");

  for(int i = 0; i < 4; i++)
  {
    canv[i] = new TCanvas();
    hist_ratio[i]->Draw();
    leg[i]->AddEntry(hist_ratio[i],"acceptance/4pi","l");
    leg[i]->Draw("same");

  }

  // for(int i = 0; i < 4; i++)
  // {
  //   cout<<i<<": "<<hist[i]->Integral()/hist_4pi[i]->Integral()<<endl;
  // }






}
