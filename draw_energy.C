void draw_energy()
{
  gStyle->SetOptStat(1111);
  gStyle->SetPalette(1,0);
  gStyle->SetOptTitle(0);

  TFile *file1 = new TFile("pimC_incl_full3.root");

  TH1F* hpim_enth[13];
  TH1F* hpip_enth[13];
  TH1F* hp_enth[13];
  TH1F* hd_enth[13];
  char name[200];


  file1->cd();
  for (int i=0;i<13;i++){

    int k=20+i*5;
    sprintf(name,"hpim_enth_%d",k);
    hpim_enth[i]= (TH1F*)file1->Get(name);

    sprintf(name,"hpip_enth_%d",k);
    hpip_enth[i]= (TH1F*)file1->Get(name);

    sprintf(name,"hp_enth_%d",k);
    hp_enth[i]= (TH1F*)file1->Get(name);

    sprintf(name,"hd_enth_%d",k);
    hd_enth[i]= (TH1F*)file1->Get(name);

  }


    TCanvas *canv[13];
    for(int i = 0; i < 13; i++)
    {
      canv[i] = new TCanvas();
      hd_enth[i]->Draw();//tu zmienić na histogram dla odp. cząstki
      int k=20+i*5;
      sprintf(name,"hd_enth_%d.pdf",k);//i tu też
      canv[i]->SaveAs(name);
    }



}
