#include "fwdet_res.h"

//#include "hgeantfwdet.h"
//#include "fwdetdef.h"
//#include "hfwdetstrawcalsim.h"
//#include "hfwdetcand.h"
//#include "hfwdetcandsim.h"

#include "hparticlecandsim.h"
#include "hparticlecand.h"
#include "hparticletool.h"
//#include "hparticlegeantdecay.h"

#include "TIterator.h"
#include "heventheader.h"
//#include "hpidtrackcand.h"
#include "hmdcdef.h"
#include "hmdcsetup.h"
#include "hmdccal1.h"
#include "hmdccal2.h"

#include "tofdef.h"
#include "htofhit.h"


#include "rpcdef.h"
#include "hrpcraw.h"
#include "hrpchit.h"
#include "hrpccluster.h"
#include "hgeantkine.h"


#include "hloop.h"
#include "hcategory.h"

#include <TCanvas.h>
#include <TStyle.h>
#include <sstream>
#include <TGraph.h>
#include <TGraphErrors.h>

#define PR(x) std::cout << "++DEBUG: " << #x << " = |" << x << "| (" << __FILE__ << ", " << __LINE__ << ")\n";

using namespace std;

Int_t   nentries;

const double D2R = 1.74532925199432955e-02;
const double R2D = 57.2957795130823229;


    double trackDistance(HParticleCand* track1, HParticleCand*  track2)
    {
      double dist;
      HGeomVector base_1, base_2, dir_1, dir_2;
      HParticleTool p_tool;

      p_tool.calcSegVector(track1->getZ(),track1->getR(),TMath::DegToRad()*track1->getPhi(),TMath::DegToRad()*track1->getTheta(),base_1,dir_1);
      p_tool.calcSegVector(track2->getZ(),track2->getR(),TMath::DegToRad()*track2->getPhi(),TMath::DegToRad()*track2->getTheta(),base_2,dir_2);
      dist=p_tool.calculateMinimumDistance(base_1,dir_1,base_2,dir_2);
      return dist;
    }

    HGeomVector trackVertex(HParticleCand* track1, HParticleCand*  track2)
    {
      HGeomVector ver;
      HGeomVector base_1, base_2, dir_1, dir_2;
      HParticleTool p_tool;

      p_tool.calcSegVector(track1->getZ(),track1->getR(),TMath::DegToRad()*track1->getPhi(),TMath::DegToRad()*track1->getTheta(),base_1,dir_1);
      p_tool.calcSegVector(track2->getZ(),track2->getR(),TMath::DegToRad()*track2->getPhi(),TMath::DegToRad()*track2->getTheta(),base_2,dir_2);
      ver=p_tool.calcVertexAnalytical(base_1,dir_1,base_2,dir_2);
      return ver;
    }

    Int_t getMotherIndex(HGeantKine* particle)
    {
      Int_t trackID=particle->getTrack();
      HGeantKine* particleParent=particle->getParent(trackID);
      Int_t parentID=0;
      if(particleParent!=0)
	parentID=particleParent->getID();
      return parentID;

      //return particle->getGeneratorInfo1();
    }


    bool isLepton(HParticleCand* particle)
    {
      double delta=0.05;
      double mquality=particle->getRichMatchingQuality();

      double dphi=particle->getDeltaPhi();
      double dtheta=particle->getDeltaTheta();


      if(particle->isFlagBit(kIsUsed))
	{
	  if(mquality==-1 || mquality>5)
	    return false;
	  if(particle->getBeta()<(1-delta) || particle->getBeta()>(1+delta))
	    return false;
	  // if(dtheta<-0.4 || dtheta>0.4)
	  // return false;
	  // if(dphi*TMath::Sin(particle->getTheta())>0.4 || dphi*TMath::Sin(particle->getTheta())<-0.4)
	  // return false;
	}
      else
	return false;

      return true;
    }


    int th_bin(float theta){
      int th=-1;

      if(theta>17.5 && theta<22.5) th=0;//theta=20+/-2.5 deg, Dtheta= 5 deg
      else if(theta>=22.5 && theta<27.5) th=1;
      else if(theta>=27.5 && theta<32.5) th=2;
      else if(theta>=32.5 && theta<37.5) th=3;
      else if(theta>=37.5 && theta<42.5) th=4;
      else if(theta>=42.5 && theta<47.5) th=5;
      else if(theta>=47.5 && theta<52.5) th=6;
      else if(theta>=52.5 && theta<57.5) th=7;
      else if(theta>=57.5 && theta<62.5) th=8;
      else if(theta>=62.5 && theta<67.5) th=9;
      else if(theta>=67.5 && theta<72.5) th=10;
      else if(theta>=72.5 && theta<77.5) th=11;
      else if(theta>=77.5 && theta<82.5) th=12;

      return th;
    }


//--------------------------------------------------------------------

Int_t fwdet_tests(HLoop * loop, const AnaParameters & anapars)
{


	    if (!loop->setInput(""))
	      {                                                    // reading file structure
		std::cerr << "READBACK: ERROR : cannot read input !" << std::endl;
		std::exit(EXIT_FAILURE);
	      }
	    gStyle->SetOptStat(1);
	    gStyle->SetOptFit(1);

	    TStopwatch timer;
	    timer.Reset();
	    timer.Start();



	    loop->printCategories();

	    HCategory* pGeantCat = (HCategory*)gHades->getCurrentEvent()->getCategory(catGeantKine);

	    HIterator* pitGeant = (HIterator *)pGeantCat->MakeIterator();


	    //GeantKine idealne czastki z modelu
    HCategory * fCatGeantKine = nullptr;
    fCatGeantKine = HCategoryManager::getCategory(catGeantKine, kTRUE, "catGeantKine");
    if (!fCatGeantKine)
    {
        cout << "No catGeantKine!" << endl;
        //exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }



    //particleCatSim kandydaci, po przefiltrowaniu przez detektor
    HCategory * particleCatSim= nullptr;
    particleCatSim = HCategoryManager::getCategory(catParticleCand, kTRUE, "catParticleCand");
    if(!particleCatSim)
      {
	cout<< "No catParticleCandSim!"<<endl;
      }


    //*****************
    HCategory * tofHit= nullptr;
    tofHit = HCategoryManager::getCategory(catTofHit, kTRUE, "catTofHit");
    if(!tofHit)
      {
	cout<< "No catTofHit!"<<endl;
      }

    HCategory * rpcHit= nullptr;
    rpcHit = HCategoryManager::getCategory(catRpcHit, kTRUE, "catRpcHit");
    if(!rpcHit)
      {
	cout<< "No catRpcHit!"<<endl;
      }

    //*******************************************************
    HEnergyLossCorrPar enLossCorr;
    enLossCorr.setDefaultPar("aug14_C7"); // "jul14_W" - Wolfram target, "jul14_C3" - Carbon 3 segments, "jul14_PE" & "aug14_PE" - PE target, "au\g14_C7" - Carbon 7 segments  
    //******************************************************
    TFile *cutFile = new TFile("./cuts_pimC_690.root");

    //cuts mom in MeV
    TCutG* cutg_pim_mass2Mom= (TCutG*)cutFile->Get("cutg_pim_mass2Mom");
    TCutG* cutg_pip_mass2Mom= (TCutG*)cutFile->Get("cutg_pip_mass2Mom");
    TCutG* cutg_p_mass2Mom= (TCutG*)cutFile->Get("cutg_p_mass2Mom");
    TCutG* cutg_d_mass2Mom= (TCutG*)cutFile->Get("cutg_d_mass2Mom");

    TCutG* cutg_pim_betaMom= (TCutG*)cutFile->Get("cutg_pim_betaMom");
    TCutG* cutg_pip_betaMom= (TCutG*)cutFile->Get("cutg_pip_betaMom");
    TCutG* cutg_p_betaMom= (TCutG*)cutFile->Get("cutg_p_betaMom");
    TCutG* cutg_d_betaMom= (TCutG*)cutFile->Get("cutg_d_betaMom");

    //same cuts but mom in GeV
    TCutG* cutg_pim_mass2Mom1= (TCutG*)cutFile->Get("cutg_pim_mass2Mom1");
    TCutG* cutg_pip_mass2Mom1= (TCutG*)cutFile->Get("cutg_pip_mass2Mom1");
    TCutG* cutg_p_mass2Mom1= (TCutG*)cutFile->Get("cutg_p_mass2Mom1");
    TCutG* cutg_d_mass2Mom1= (TCutG*)cutFile->Get("cutg_d_mass2Mom1");

    TCutG* cutg_pim_betaMom1= (TCutG*)cutFile->Get("cutg_pim_betaMom1");
    TCutG* cutg_pip_betaMom1= (TCutG*)cutFile->Get("cutg_pip_betaMom1");
    TCutG* cutg_p_betaMom1= (TCutG*)cutFile->Get("cutg_p_betaMom1");
    TCutG* cutg_d_betaMom1= (TCutG*)cutFile->Get("cutg_d_betaMom1");


    //************************** HISTOS ********************

    TH2F *hBetavsMom=new TH2F("betavsMom","betavsMom; p*q [GeV]; beta",500,-1.4,1.4,500,0,1.5);
    TH2F *hmomvsMass=new TH2F("momvsMass","momvsMass; p [GeV]; mass2",500,0,1.4,1000,-1000000,10000000);
    TH2F *hmomvsMass1=new TH2F("momvsMass1","momvsMass1; p*q [GeV]; mass2",1000,-1.4,1.4,1000,-1000000,10000000);
    TH1F* hmass= new TH1F("mass","mass",1000,-1000000,10000000);
    TH1F* htof = new TH1F("tof","tof",100,0,100);
    TH1F* hvertReco_z = new TH1F("hvertReco_z","hvertReco_z",10000,-1200, 500);
    TH1F* hvert_z = new TH1F("hvert_z","hvert_z",10000,-1200,500);
    TH1F* hvertReco_x = new TH1F("hvertReco_x", "hvertReco_x", 10000, -1200, 500);
    TH1F* hvert_x = new TH1F("hvert_x", "hvert_x", 10000, -1200, 500);
    TH1F* hvertReco_y = new TH1F("hvertReco_y", "hvertReco_y", 10000, -1200, 500);
    TH1F* hvert_y = new TH1F("hvert_y", "hvert_y", 10000, -1200, 500);
    TH2F* hvertReco_xy = new TH2F("hvertReco_xy", "hvertReco_xy;x;y", 1000, -50, 50, 1000, -50, 50);
    TH2F* hvertReco_xz = new TH2F("hvertReco_xz", "hvertReco_xz;x;z", 1000, -50, 50, 1500, -100, 50);
    TH2F* hvertReco_yz = new TH2F("hvertReco_yz", "hvertReco_yz;y;z", 1000, -50, 50, 1500, -100, 50);
    TH1F* htrMult = new TH1F("htrMult","htrMult",15,0,15);


    TH2F *hp_BetavsMom=new TH2F("p_betavsMom","p_betavsMom; p [GeV]; beta",200,0,1.4,200,0,1.2);
    TH2F *hd_BetavsMom=new TH2F("d_betavsMom","d_betavsMom; p [GeV]; beta",200,0,1.4,200,0,1.2);
    TH2F *hpip_BetavsMom=new TH2F("pip_betavsMom","pip_betavsMom; p [GeV]; beta",200,0,1.4,200,0,1.2);
    TH2F *hpim_BetavsMom=new TH2F("pim_betavsMom","pim_betavsMom; p [GeV]; beta",200,0,1.4,200,0,1.2);


    TH2F *hp_ThetavsMom=new TH2F("p_thetavsMom","p_thetavsMom;p [GeV];theta",100,0,1.4,100,0,100);
    TH2F *hd_ThetavsMom=new TH2F("d_thetavsMom","d_thetavsMom;p [GeV];theta",100,0,1.6,100,0,100);
    TH2F *hpip_ThetavsMom=new TH2F("pip_thetavsMom","pip_thetavsMom;p [GeV];theta",100,0,1.4,100,0,100);
    TH2F *hpim_ThetavsMom=new TH2F("pim_thetavsMom","pim_thetavsMom;p [GeV];theta",100,0,1.4,100,0,100);

    TH1F* hp_tof = new TH1F("p_tof","p_tof",100,0,100);
    TH1F* hd_tof = new TH1F("d_tof","d_tof",100,0,100);
    TH1F* hpip_tof = new TH1F("pip_tof","pip_tof",100,0,100);
    TH1F* hpim_tof = new TH1F("pim_tof","pim_tof",100,0,100);


    //_4pi dla idealnych czastek, przed detektorem
    TH2F *hp_ThetavsMom_4pi=new TH2F("p_thetavsMom_4pi","p_thetavsMom_4pi;p [GeV];theta",200,0,1.4,200,0,200);
    TH2F *hd_ThetavsMom_4pi=new TH2F("d_thetavsMom_4pi","d_thetavsMom_4pi;p [GeV];theta",200,0,1.6,200,0,200);
    TH2F *hpip_ThetavsMom_4pi=new TH2F("pip_thetavsMom_4pi","pip_thetavsMom_4pi;p [GeV];theta",200,0,1.4,200,0,200);
    TH2F *hpim_ThetavsMom_4pi=new TH2F("pim_thetavsMom_4pi","pim_thetavsMom_4pi;p [GeV];theta",200,0,1.4,200,0,200);

    //************************************************
    //************************************************

    TH1F* hp_theta = new TH1F("p_theta","p_theta;theta;",90,0,180);
    TH1F* hp_theta_4pi = new TH1F("p_theta_4pi","p_theta_4pi;theta;",90,0,180);
    TH1F* hd_theta = new TH1F("d_theta", "d_theta;theta;", 90, 0, 180);
    TH1F* hd_theta_4pi = new TH1F("d_theta_4pi", "d_theta_4pi;theta;", 90, 0, 180);
    TH1F* hpip_theta = new TH1F("pip_theta", "pip_theta;theta;", 90, 0, 180);
    TH1F* hpip_theta_4pi = new TH1F("pip_theta_4pi", "pip_theta_4pi;theta;", 90, 0 ,180);
    TH1F* hpim_theta = new TH1F("pim_theta", "pim_theta;theta;", 90, 0, 180);
    TH1F* hpim_theta_4pi = new TH1F("pim_theta_4pi", "pim_theta_4pi;theta;", 90, 0, 180);

    TH1F* hp_y = new TH1F("p_y", "p_y;y;", 100, 0, 1.);
    TH1F* hp_y_4pi = new TH1F("p_y_4pi", "p_y_4pi;y;", 100, -1., 2.);
    TH1F* hp_y_4pi_crop = new TH1F("p_y_4pi_crop", "p_y_4pi_crop;y;", 100, 0, 1.);
    TH1F* hd_y = new TH1F("d_y", "d_y;y", 100, 0., 1.);
    TH1F* hd_y_4pi = new TH1F("d_y_4pi", "d_y_4pi;y;", 100, -1., 2.);
    TH1F* hd_y_4pi_crop = new TH1F("d_y_4pi_crop", "d_y_4pi_crop;y;", 100, 0., 1.);
    TH1F* hpip_y = new TH1F("pip_y", "pip_y;y;", 100, 0., 2.);
    TH1F* hpip_y_4pi = new TH1F("pip_y_4pi", "pip_y_4pi;y;", 100, -2., 3.);
    TH1F* hpip_y_4pi_crop = new TH1F("pip_y_4pi_crop", "pip_y_4pi_crop;y;", 100, 0., 2.);
    TH1F* hpim_y = new TH1F("pim_y", "pim_y;y;", 100, 0., 2.);
    TH1F* hpim_y_4pi = new TH1F("pim_y_4pi", "pim_y_4pi;y;", 100, -2., 3.);
    TH1F* hpim_y_4pi_crop = new TH1F("pim_y_4pi_crop", "pim_y_4pi_crop;y;", 100, 0., 2.);

    TH1F* hp_pt = new TH1F("p_pt", "p_pt;pt [GeV/c];", 75, 0., 1.5);
    TH1F* hp_pt_4pi = new TH1F("p_pt_4pi", "p_pt_4pi;pt [GeV/c];", 75, 0, 1.5);
    TH1F* hd_pt = new TH1F("d_pt", "d_pt;pt [GeV/c];", 75, 0., 1.5);
    TH1F* hd_pt_4pi = new TH1F("d_pt_4pi", "d_pt_4pi;pt [GeV/c];", 75, 0., 1.5);
    TH1F* hpip_pt = new TH1F("pip_pt", "pip_pt;pt [GeV/c];", 75, 0., 1.5);
    TH1F* hpip_pt_4pi = new TH1F("pip_pt_4pi", "pip_pt_4pi;pt [GeV/c];", 75, 0., 1.5);
    TH1F* hpim_pt = new TH1F("pim_pt", "pim_pt; pt [GeV/c];", 75, 0., 1.5);
    TH1F* hpim_pt_4pi = new TH1F("pim_pt_4pi", "pim_pt_4pi; pt [GeV/c];", 75, 0., 1.5);

    TH1F* hp_T = new TH1F("p_T", "p_T;T [GeV];", 80, 0., 0.8);
    TH1F* hp_T_4pi = new TH1F("p_T_4pi", "p_T_4pi;T [GeV];", 80, 0., 0.8);
    TH1F* hd_T = new TH1F("d_T", "d_T;T [GeV];", 80, 0., 0.8);
    TH1F* hd_T_4pi = new TH1F("d_T_4pi", "d_T_4pi;T [GeV];", 80, 0., 0.8);
    TH1F* hpip_T = new TH1F("pip_T", "pip_T;T [GeV];", 70, 0, 0.7);
    TH1F* hpip_T_4pi = new TH1F("pip_T_4pi", "pip_T_4pi;T [GeV];", 70, 0, 0.7);
    TH1F* hpim_T = new TH1F("pim_T", "pim_T;T [GeV];", 100, 0, 1.);
    TH1F* hpim_T_4pi = new TH1F("pim_T_4pi", "pim_T_4pi; T[GeV];", 100, 0, 1.);

    TH1F* hp_p = new TH1F("p_p", "p_p; p [GeV];", 140, 0., 1.4);
    TH1F* hp_p_4pi = new TH1F("p_p_4pi", "p_p_4pi; p [GeV];", 140, 0., 1.4);
    TH1F* hd_p = new TH1F("d_p", "d_p; p [GeV];", 160, 0., 1.6);
    TH1F* hd_p_4pi = new TH1F("d_p_4pi", "d_p_4pi; p[GeV];", 160, 0., 1.6);
    TH1F* hpip_p = new TH1F("pip_p", "pip_p; p [GeV];", 70, 0., 0.7);
    TH1F* hpip_p_4pi = new TH1F("pip_p_4pi", "pip_p_4pi; p [GeV];", 70, 0., 0.7);
    TH1F* hpim_p = new TH1F("pim_p", "pim_p; p [GeV];", 130, 0., 1.3);
    TH1F* hpim_p_4pi = new TH1F("pim_p_4pi", "pim_p_4pi; p [GeV];", 130, 0., 1.3);


    TH1F* hpim_enth[13];
    TH1F* hpip_enth[13];
    TH1F* hp_enth[13];
    TH1F* hd_enth[13];
    char name[200];
    char name2[200];

    for (int i=0;i<13;i++){

      int k=20+i*5;
      sprintf(name,"hpim_enth_%d",k);
      sprintf(name2, "hpim_enth_%d;T [GeV];", k);
      hpim_enth[i]= new TH1F(name, name2,160,0,0.8);;

      sprintf(name,"hpip_enth_%d",k);
      sprintf(name2, "hpip_enth_%d;T [GeV];", k);
      hpip_enth[i]= new TH1F(name,name2,160,0,0.8);;

      sprintf(name,"hp_enth_%d",k);
      sprintf(name2, "hp_enth_%d;T [GeV];", k);
      hp_enth[i]= new TH1F(name,name2,160,0,0.8);;

      sprintf(name,"hd_enth_%d",k);
      sprintf(name2, "hd_enth_%d;T [GeV];", k);
      hd_enth[i]= new TH1F(name,name2,160,0,0.8);;

    }

    //***************************************************************
    //********************************************

    double pion_momentum = 685.;//MeV/c
    double pion_energy = sqrt(pion_momentum*pion_momentum + 139.56995*139.56995);

   TLorentzVector *proj;
   proj = new TLorentzVector(0,0,pion_momentum, pion_energy);
   //TLorentzVector proj(0,0,pion_momentum, pion_energy);


    //-------TARGET
   TLorentzVector *targ;
   targ = new TLorentzVector(0,0,0, 938.27231);
    //TLorentzVector targ(0,0,0, 938.27231);


    //--------------   BEAM+TARGET
    TLorentzVector *beam;
    beam = new TLorentzVector(0,0,0,0);
    //TLorentzVector beam(0,0,0,0);
    *beam = *proj + *targ;
    //beam = proj + targ;

     Int_t entries = loop->getEntries();
    //     //setting numbers of events regarding the input number of events by the user
    if (anapars.events < entries and anapars.events >= 0 ) entries = anapars.events;



    TFile * output_file = TFile::Open(anapars.outfile, "RECREATE");
    output_file->cd();

    cout << "NEW ROOT TREE , vertex_ana" << endl;
    Int_t evnb;

 //  if (fChain == 0) return;

    HTofHit *tof = nullptr;
    HRpcHit *rpc = nullptr;


    vector<HParticleCandSim*> prot, pim,pip;
    vector<HGeantKine*> protKine, pimKine, pipKine;

    HGeomVector epVertex, emVertex, vertex;
    HGeomVector xVertex, vertex_kine;
    Int_t aTrack, aID;
    Float_t xMom=0.,yMom=0.,zMom=0.;
    Int_t aPar, aMed, aMech;
    Float_t genInfo=0., genWeight=0.;
    Float_t genInfo1=0., genWeight1=0.;
    Float_t genInfo2=0., genWeight2=0.;
    Float_t genWeight3=0.;
    TVector3    vec1, vec2, vec3, vec4;
    TLorentzVector*  pim_LAB = new TLorentzVector(0,0,0,0);
    TLorentzVector*  pip_LAB = new TLorentzVector(0,0,0,0);
    TLorentzVector*  p_LAB = new TLorentzVector(0,0,0,0);
    TLorentzVector*  d_LAB = new TLorentzVector(0, 0, 0, 0);
    int charTr=0;


    for (Long_t event=0; event<entries; event++)
    {
      charTr=0;
      loop->nextEvent(event);
      if(!(event%10000))   cout<<"event no. "<<event<<endl;


      evnb=event;
      //fChain->GetEntry(event);
      //cout<< eventNum<<endl;

      int tof_ev=tofHit->getEntries();
      int rpc_ev=rpcHit->getEntries();

      prot.clear();
      pim.clear();
      pip.clear();
      protKine.clear();
      pimKine.clear();
      pipKine.clear();



	int hnum=particleCatSim->getEntries();
	int knum=fCatGeantKine->getEntries();

	HParticleCandSim* partH=nullptr;
	HGeantKine* kine=nullptr;


	HParticleTool tool;

	HGeomVector vertexL;
	HGeomVector vertexDL;
	HGeomVector vertexL1520;
	HGeomVector dirL;
	HGeomVector dirDL;
	HGeomVector dirL1520,dirL1520_1;


	HGeomVector base_Tg, dir_Tg;

	base_Tg.setX(0);
	base_Tg.setY(0);
	base_Tg.setZ(-25.);
	dir_Tg.setX(0);
	dir_Tg.setY(0);
	dir_Tg.setZ(1);

	HGeomVector ver_L1520Tg, ver_L1520TgFT, baseL1520, ver_LTg, ver_LTgFT;

	int tofhit=0;
	int rpchit=0;

	//************************
	if(tof_ev){
	  for (int i=0; i<tof_ev; i++){
	    tof = HCategoryManager::getObject(tof, catTofHit, i);

	    int sec = tof->getSector();
	    //int mod = mdc->getModule();

	    //tof_sec[sec]++;
	    tofhit++;
	    //cout<<event<<" "<<i<<" "<<sec<<endl;
	  }
	}//tof_ev


	if(rpc_ev){
	  for(int i=0; i<rpc_ev; i++){
	    rpc = HCategoryManager::getObject(rpc, catRpcHit, i);

	    int sec = rpc->getSector();
	    //int mod = mdc->getModule();

	    //rpc_sec[sec]++;
	    rpchit++;
	    //cout<<event<<" "<<i<<" "<<sec<<endl;
	  }
	}//rpc_ev

	//****************************************
	//cout<<rpchit<<" "<<tofhit<<endl;
	int trigM2=rpchit+tofhit;

	//cout<<hnum<<" "<<trigM2<<endl;
	//if(hnum && trigM2>=2) cout<<"ok"<<endl;


	//****************************************
	if (hnum){//without M2 trigger
	//if (hnum && trigM2>=2){ //with M2

          //HADES
	  //##############################
	  for (int i=0;i<hnum;i++){


	  partH=HCategoryManager::getObject(partH, particleCatSim,i);


	  //cout<<i<<"  "<<partH->getSystem()<<" "<<partH->getChi2()<<" "<<partH->getInnerSegmentChi2()<<endl;
	  float nerbyFit=0;
	  float nerbyUnFit=0;

	  //nerbyFit=fabs(partH->getAngleToNearbyFittedInner())*TMath::RadToDeg();
	  //nerbyUnFit=fabs(partH->getAngleToNearbyUnfittedInner())*TMath::RadToDeg();
	  nerbyFit=fabs(partH->getAngleToNearbyFittedInner());
	  nerbyUnFit=fabs(partH->getAngleToNearbyUnfittedInner());


	  //if(partH->isFlagBit(kIsUsed) && partH->getGeantParentTrackNum()==0 && partH->getGeantGrandParentPID()==-1){
	  //if(partH->isFlagBit(kIsUsed) && partH->getGeantParentTrackNum()==0 && partH->getGeantParentPID()==-1){

	  //cout<<partH->getGeantPID()<<endl;
	  if(partH->isFlagBit(kIsUsed) ){//wybiera dobrze zrekonstruowane trajektorie, przechodza tylko M2trig >= 1


	    double mom=partH->getMomentum();
	    double beta=partH->getBeta();
	    float q=partH->getCharge();
	    double mass = mom*mom * (  1. / (beta*beta)  - 1. ) ;
	    double tof=partH->getTof();
	    double theta = partH->getTheta();
	    double phi = partH->getPhi();
	    double y=partH->Rapidity();
	    double pt=partH->Pt();


	    double tracklength = partH->getDistanceToMetaHit();
	    HVertex& vertexreco = gHades->getCurrentEvent()->getHeader()->getVertexReco();
	    double eVertReco_x=vertexreco.getX();
	    double eVertReco_y=vertexreco.getY();
	    double eVertReco_z=vertexreco.getZ();

	    HVertex& vertex = gHades->getCurrentEvent()->getHeader()->getVertex();
	    double eVert_x=vertex.getX();
	    double eVert_y=vertex.getY();
	    double eVert_z=vertex.getZ();


	    hvertReco_z->Fill(eVertReco_z);
	    hvert_z->Fill(eVert_z);
	    hvertReco_x->Fill(eVertReco_x);
	    hvert_x->Fill(eVert_x);
	    hvertReco_y->Fill(eVertReco_y);
	    hvert_y->Fill(eVert_y);
	    hvertReco_xy->Fill(eVertReco_x, eVertReco_y);
	    hvertReco_xz->Fill(eVertReco_x, eVertReco_z);
	    hvertReco_yz->Fill(eVertReco_y, eVertReco_z);

	    hBetavsMom->Fill(mom*q/1000,beta);
	    hmass->Fill(mass);
	    htof->Fill(tof);
	    hmomvsMass->Fill(mom/1000,mass);
	    hmomvsMass1->Fill(mom*q/1000,mass);


	    if(eVertReco_z>-500){ //dobrze zrekonstruowane trajektorie (punkt reakcji musi byc w obrebie tarczy) przechodza tylko hnum>=2 M2trig>=2



	    if(partH->getCharge()!=0)charTr++; ;

	    //cout<<i<<"  "<<partH->getSystem()<<" "<<partH->getChi2()<<" "<<partH->getInnerSegmentChi2()<<endl;

	    //14 = proton
	    //if(partH->getGeantPID()==14){
	    //if(mass > 200000. && mass < 2200000. && mom*q > 0.){
	    //if(partH->getGeantPID()==14 && mass > 100000. && mass < 2200000. && mom*q >0.){
	    //if(cutg_p_mass2Mom->IsInside(mom*q,mass)){
	      if(cutg_p_betaMom->IsInside(mom*q,beta)){

	      float mom_corr=partH->getCorrectedMomentumPID(14);
	      partH->calc4vectorProperties(HPhysicsConstants::mass(14));
	      hp_ThetavsMom->Fill(mom_corr/1000,partH->getTheta());
	      hp_BetavsMom->Fill(mom_corr/1000,partH->getBeta());
	      hp_theta->Fill(partH->getTheta());
	      //hp_ThetavsMom->Fill(partH->getMomentum()/1000,partH->getTheta());
	      //hp_BetavsMom->Fill(partH->getMomentum()/1000,partH->getBeta());
	      //hp_T->Fill((partH->T()-partH->M())/1000);
	      //hp_y->Fill(partH->Rapidity());
	      //hp_pt->Fill(partH->Pt()/1000);
	      //hp_p->Fill(partH->getMomentum()/1000);
	      hp_p->Fill(mom_corr/1000);

	      TLorentzVector p_corr1;
	      p_corr1.SetXYZM(mom_corr*sin(D2R*theta)*cos(D2R*phi),mom_corr*sin(D2R*theta)*sin(D2R*phi),mom_corr*cos(D2R*theta),partH->M());
	      hp_pt->Fill(p_corr1.Pt()/1000);
	      hp_y->Fill(p_corr1.Rapidity());
	      
	      float ekin=sqrt(mom_corr*mom_corr+partH->M()*partH->M())-partH->M();
	      hp_T->Fill(ekin/1000);
	      
	      int thbin=th_bin(theta);
	      if(thbin>-1)hp_enth[thbin]->Fill(ekin/1000.);


	    }

	    //45 = deuteron
	    //if(partH->getGeantPID()==45){
	    //if(mass > 2200000. && mass < 6100000. && mom*q > 0.){
	    //if(partH->getGeantPID()==45 && mass > 2200000. && mass < 6100000. && mom*q > 0.){
	      //if(cutg_d_mass2Mom->IsInside(mom*q, mass)){
	      if(cutg_d_betaMom->IsInside(mom*q, beta)){
	      
	      float mom_corr=partH->getCorrectedMomentumPID(45);

	      partH->calc4vectorProperties(HPhysicsConstants::mass(45));
	      //hd_ThetavsMom->Fill(partH->getMomentum()/1000,partH->getTheta());
	      //hd_BetavsMom->Fill(partH->getMomentum()/1000,partH->getBeta());
	      hd_ThetavsMom->Fill(mom_corr/1000,partH->getTheta());
	      hd_BetavsMom->Fill(mom_corr/1000,partH->getBeta());
	      hd_theta->Fill(partH->getTheta());
	      //hd_y->Fill(partH->Rapidity());
	      //hd_pt->Fill(partH->Pt()/1000);
	      //hd_T->Fill((partH->T()-partH->M())/1000);
	      //hd_p->Fill(partH->getMomentum()/1000);
	      hd_p->Fill(mom_corr/1000);

	      TLorentzVector p_corr1;
	      p_corr1.SetXYZM(mom_corr*sin(D2R*theta)*cos(D2R*phi),mom_corr*sin(D2R*theta)*sin(D2R*phi),mom_corr*cos(D2R*theta),partH->M());

	      hd_pt->Fill(p_corr1.Pt()/1000);
	      hd_y->Fill(p_corr1.Rapidity());
	      
	      float ekin=sqrt(mom_corr*mom_corr+partH->M()*partH->M())-partH->M();
	      hd_T->Fill(ekin/1000);
	      int thbin=th_bin(theta);
	      if(thbin>-1)hd_enth[thbin]->Fill(ekin/1000.);
	    }

	    //9 = pion-
	    //if(partH->getGeantPID()==9){
	    //if(mass > -20000. && mass < 50000. && mom*q < 0.){
	    //if(partH->getGeantPID()==9 && mass > -20000. && mass < 50000. && mom*q < 0.){
	      //if(cutg_pim_mass2Mom->IsInside(mom*q, mass)){
	    if(cutg_pim_betaMom->IsInside(mom*q, beta)){
	      
	      float mom_corr=partH->getCorrectedMomentumPID(9);

	      partH->calc4vectorProperties(HPhysicsConstants::mass(9));
	      //hpim_ThetavsMom->Fill(partH->getMomentum()/1000,partH->getTheta());
	      //hpim_BetavsMom->Fill(partH->getMomentum()/1000,partH->getBeta());
	      hpim_ThetavsMom->Fill(mom_corr/1000,partH->getTheta());
	      hpim_BetavsMom->Fill(mom_corr/1000,partH->getBeta());
	      hpim_theta->Fill(partH->getTheta());
	      //hpim_y->Fill(partH->Rapidity());
	      //hpim_pt->Fill(partH->Pt()/1000);
	      //hpim_T->Fill((partH->T()-partH->M())/1000);
	      //hpim_p->Fill(partH->getMomentum()/1000);
	      hpim_p->Fill(mom_corr/1000);

	      TLorentzVector p_corr1;
	      p_corr1.SetXYZM(mom_corr*sin(D2R*theta)*cos(D2R*phi),mom_corr*sin(D2R*theta)*sin(D2R*phi),mom_corr*cos(D2R*theta),partH->M());
	      hpim_pt->Fill(p_corr1.Pt()/1000);
	      hpim_y->Fill(p_corr1.Rapidity());
	      
	      float ekin=sqrt(mom_corr*mom_corr+partH->M()*partH->M())-partH->M();
	      hpim_T->Fill(ekin/1000);
	      int thbin=th_bin(theta);
	      if(thbin>-1)hpim_enth[thbin]->Fill(ekin/1000.);

	    }

	    //8 = pion+
	    //if(partH->getGeantPID()==8){
	    //if(mass > 10000. && mass < 70000. && mom*q > 0.){
	    //if(partH->getGeantPID()==8 && mass > 10000. && mass < 70000. && mom*q > 0.){
	    //if(cutg_pip_mass2Mom->IsInside(mom*q, mass)){
	    if(cutg_pip_betaMom->IsInside(mom*q, beta)){

	      
	      float mom_corr=partH->getCorrectedMomentumPID(8);
	      partH->calc4vectorProperties(HPhysicsConstants::mass(8));
	      //hpip_ThetavsMom->Fill(partH->getMomentum()/1000,partH->getTheta());
	      //hpip_BetavsMom->Fill(partH->getMomentum()/1000,partH->getBeta());
	      //hpip_y->Fill(partH->Rapidity());
	      //hpip_pt->Fill(partH->Pt()/1000);
	      //hpip_T->Fill((partH->T()-partH->M())/1000);

	      hpip_theta->Fill(partH->getTheta());
	      hpip_ThetavsMom->Fill(mom_corr/1000,partH->getTheta());
	      hpip_BetavsMom->Fill(mom_corr/1000,partH->getBeta());
	      //hpip_p->Fill(partH->getMomentum()/1000);
	      hpip_p->Fill(mom_corr/1000);

	      TLorentzVector p_corr1;
	      p_corr1.SetXYZM(mom_corr*sin(D2R*theta)*cos(D2R*phi),mom_corr*sin(D2R*theta)*sin(D2R*phi),mom_corr*cos(D2R*theta),partH->M());
	      hpip_pt->Fill(p_corr1.Pt()/1000);
	      hpip_y->Fill(p_corr1.Rapidity());
	      
	      float ekin=sqrt(mom_corr*mom_corr+partH->M()*partH->M())-partH->M();
	      hpip_T->Fill(ekin/1000);
	      int thbin=th_bin(theta);
	      if(thbin>-1)hpip_enth[thbin]->Fill(ekin/1000);

	    }


	  }//eVertReco_z
	  }//kIsUsed

	  }
	}
	htrMult->Fill(charTr);


	//**********************************************************************
	//*************** HGEANT KINE analysis*****************************************
     pitGeant->Reset();
     while((kine = (HGeantKine*)pitGeant->Next()) != 0){

       //for(int p=0;p<knum;p++) {
       //kine=HCategoryManager::getObject(kine, fCatGeantKine, p);

       kine->getParticle( aTrack, aID );
       kine->getCreator ( aPar, aMed, aMech);
       kine->getGenerator( genInfo, genInfo1, genInfo2);
       kine->getGenerator( genInfo, genWeight );
       kine->getMomentum( xMom, yMom, zMom );

       //int kineID=kine->getID();
       //int mech=kine->getMechanism();
       //int kineparentID=getMotherIndex(kine);
       //HGeomVector lambdaVertex;
       //      cout<<aID<<" "<<aMech<<endl;

       //aMech == 0 zostawiamy tylko pierwotne czastki, bez reakcji wtornych po drodze
       if(aID==14 && aMech==0){//proton from primary vertex

	 vec1.SetXYZ( xMom, yMom, zMom );
	 p_LAB->SetVectM(vec1, 938.272);

	 hp_ThetavsMom_4pi->Fill(p_LAB->P()/1000,p_LAB->Theta()*TMath::RadToDeg());
	 hp_theta_4pi->Fill(p_LAB->Theta()*TMath::RadToDeg());
	 hp_y_4pi->Fill(p_LAB->Rapidity());
	 hp_y_4pi_crop->Fill(p_LAB->Rapidity());
	 hp_pt_4pi->Fill(p_LAB->Pt()/1000);
	 hp_T_4pi->Fill((p_LAB->T()-p_LAB->M())/1000);
	 hp_p_4pi->Fill(p_LAB->P()/1000);

       }

       if(aID==9 && aMech==0)//pion-
	 {
	   vec2.SetXYZ( xMom, yMom, zMom );
	   pim_LAB->SetVectM(vec2, 139.56995);

	   hpim_ThetavsMom_4pi->Fill(pim_LAB->P()/1000,pim_LAB->Theta()*TMath::RadToDeg());
	   hpim_theta_4pi->Fill(pim_LAB->Theta()*TMath::RadToDeg());
	   hpim_y_4pi->Fill(pim_LAB->Rapidity());
	   hpim_y_4pi_crop->Fill(pim_LAB->Rapidity());
	   hpim_pt_4pi->Fill(pim_LAB->Pt()/1000);
	   hpim_T_4pi->Fill((pim_LAB->T()-pim_LAB->M())/1000);
	   hpim_p_4pi->Fill(pim_LAB->P()/1000);

	   //kine->getVertex(lambdaVertex);
	   //h2IIpions->Fill(kine->getTotalMomentum(),kine->getThetaDeg());
		//hEIZpionSim->Fill(lambdaVertex.getZ());
	 }

       if(aID==8 && aMech==0){//pion+ from primary vertex
	 vec3.SetXYZ( xMom, yMom, zMom );
	 pip_LAB->SetVectM(vec3, 139.56995);

	 hpip_ThetavsMom_4pi->Fill(pip_LAB->P()/1000,pip_LAB->Theta()*TMath::RadToDeg());
	 hpip_theta_4pi->Fill(pip_LAB->Theta()*TMath::RadToDeg());
	 hpip_y_4pi->Fill(pip_LAB->Rapidity());
	 hpip_y_4pi_crop->Fill(pip_LAB->Rapidity());
	 hpip_pt_4pi->Fill(pip_LAB->Pt()/1000);
	 hpip_T_4pi->Fill((pip_LAB->T()-pip_LAB->M())/1000);
	 hpip_p_4pi->Fill(pip_LAB->P()/1000);

       }

       if(aID==45 && aMech==0)//deuteron
	 {
	   vec4.SetXYZ(xMom, yMom, zMom);
	   d_LAB->SetVectM(vec4, 1875.612943);

	   hd_ThetavsMom_4pi->Fill(d_LAB->P()/1000, d_LAB->Theta()*TMath::RadToDeg());
	   hd_theta_4pi->Fill(d_LAB->Theta()*TMath::RadToDeg());
	   hd_y_4pi->Fill(d_LAB->Rapidity());
	   hd_y_4pi_crop->Fill(d_LAB->Rapidity());
	   hd_pt_4pi->Fill(d_LAB->Pt()/1000);
	   hd_T_4pi->Fill((d_LAB->T()-d_LAB->M())/1000);
	   hd_p_4pi->Fill(d_LAB->P()/1000);
	 }

	    //if(kineID==14 && kine->getThetaDeg()<6.5)h2FDsimProtons->Fill(kine->getTotalMomentum(),kine->getThetaDeg());
    }




	//******************************************************end kine

    }

    htrMult->Write();
    hBetavsMom->Write();
    hmass->Write();
    hmomvsMass->Write();
    hmomvsMass1->Write();
    hvertReco_z->Write();
    hvert_z->Write();
    hvertReco_x->Write();
    hvert_x->Write();
    hvertReco_y->Write();
    hvert_y->Write();
    hvertReco_xy->Write();
    hvertReco_xz->Write();
    hvertReco_yz->Write();



    hp_ThetavsMom->Write();
    hd_ThetavsMom->Write();
    hpip_ThetavsMom->Write();
    hpim_ThetavsMom->Write();

    hp_ThetavsMom_4pi->Write();
    hd_ThetavsMom_4pi->Write();
    hpip_ThetavsMom_4pi->Write();
    hpim_ThetavsMom_4pi->Write();

    hp_BetavsMom->Write();
    hd_BetavsMom->Write();
    hpip_BetavsMom->Write();
    hpim_BetavsMom->Write();

    htof->Write();
    hp_tof->Write();
    hd_tof->Write();
    hpip_tof->Write();
    hpim_tof->Write();

    hp_theta->Write();
    hd_theta->Write();
    hpip_theta->Write();
    hpim_theta->Write();

    hp_theta_4pi->Write();
    hd_theta_4pi->Write();
    hpip_theta_4pi->Write();
    hpim_theta_4pi->Write();

    hp_y->Write();
    hd_y->Write();
    hpip_y->Write();
    hpim_y->Write();

    hp_y_4pi->Write();
    hd_y_4pi->Write();
    hpip_y_4pi->Write();
    hpim_y_4pi->Write();

    hp_y_4pi_crop->Write();
    hd_y_4pi_crop->Write();
    hpip_y_4pi_crop->Write();
    hpim_y_4pi_crop->Write();

    hp_pt->Write();
    hd_pt->Write();
    hpip_pt->Write();
    hpim_pt->Write();

    hp_pt_4pi->Write();
    hd_pt_4pi->Write();
    hpip_pt_4pi->Write();
    hpim_pt_4pi->Write();

    hp_T->Write();
    hd_T->Write();
    hpip_T->Write();
    hpim_T->Write();

    hp_T_4pi->Write();
    hd_T_4pi->Write();
    hpip_T_4pi->Write();
    hpim_T_4pi->Write();

    hp_p->Write();
    hd_p->Write();
    hpip_p->Write();
    hpim_p->Write();

    hp_p_4pi->Write();
    hd_p_4pi->Write();
    hpip_p_4pi->Write();
    hpim_p_4pi->Write();

    for (int i=0;i<13;i++){

      hpim_enth[i]->Write();
      hpip_enth[i]->Write();
      hp_enth[i]->Write();
      hd_enth[i]->Write();

   }


}
