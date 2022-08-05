
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include <vector>
#include "TVector3.h"


using namespace std;

void analyzeBhadronTree()
{   
    //  read in file 

  TString inFileName = "./pythia";                                // . one dot means current working directory, and two dots means parent directory
  inFileName += "8Jets_example.root";
  TFile *inFile = TFile::Open(inFileName);                       //"Open" is an object and the "TFile" class function

    //Histograms
  TH1F* zH  = new TH1F("zH",  ";z;",            200,   0., 2.);            
  TH1F* jtH  = new TH1F("jtH",  ";jt;",              200,   0., 10.);
  TH1F* rH  = new TH1F("rH",  "; R;",              200,   0., 10.);
	//TH1F* pzH  = new TH1F("pzH",  "pz",              100,   0., 5.); 
	//TH1F* eH  = new TH1F("eH",  "e",              100,   0., 50.);

  double z, Jt;
  
	// create and read in tree
  double jPx, jPy, jPz, jE, bPx, bPy, bPz, bE ;
  int  nConstituents, eventNum, beventNum ;

  double JPsiPx, JPsiPy, JPsiPz, JPsiE; 
  double KaonPx, KaonPy, KaonPz, KaonE; 
  double Mu1Px, Mu1Py, Mu1Pz, Mu1E;
  double Mu2Px, Mu2Py, Mu2Pz, Mu2E;

  TTree *JetsTree = 	(TTree*) inFile->Get("jets");   
    
  JetsTree->SetBranchAddress("eventNum", &eventNum);
  JetsTree->SetBranchAddress("jPx", &jPx);
  JetsTree->SetBranchAddress("jPy", &jPy);
  JetsTree->SetBranchAddress("jPz", &jPz);
  JetsTree->SetBranchAddress("jE", &jE);
  JetsTree->SetBranchAddress("nConstituents", &nConstituents);

    
  TTree *BhadronTree = 	(TTree*) inFile->Get("B-hadron");
    
  BhadronTree->SetBranchAddress("eventNum", &beventNum);
  BhadronTree->SetBranchAddress("bPx", &bPx);
  BhadronTree->SetBranchAddress("bPy", &bPy);
  BhadronTree->SetBranchAddress("bPz", &bPz);
  BhadronTree->SetBranchAddress("bE", &bE);

  BhadronTree->SetBranchAddress("JPsiPx", &JPsiPx);
  BhadronTree->SetBranchAddress("JPsiPy", &JPsiPy);
  BhadronTree->SetBranchAddress("JPsiPz", &JPsiPz);
  BhadronTree->SetBranchAddress("JPsiE", &JPsiE);

  
  BhadronTree->SetBranchAddress("KaonPx", &KaonPx);
  BhadronTree->SetBranchAddress("KaonPy", &KaonPy);
  BhadronTree->SetBranchAddress("KaonPz", &KaonPz);
  BhadronTree->SetBranchAddress("KaonE", &KaonE);


  BhadronTree->SetBranchAddress("Mu1Px", &Mu1Px);
  BhadronTree->SetBranchAddress("Mu1Py", &Mu1Py);
  BhadronTree->SetBranchAddress("Mu1Pz", &Mu1Pz);
  BhadronTree->SetBranchAddress("Mu1E", &Mu1E);


  BhadronTree->SetBranchAddress("Mu2Px", &Mu2Px);
  BhadronTree->SetBranchAddress("Mu2Py", &Mu2Py);
  BhadronTree->SetBranchAddress("Mu2Pz", &Mu2Pz);
  BhadronTree->SetBranchAddress("Mu2E", &Mu2E);
    
  // loop over the entries of trees 
                         
  for (int i=0; i<BhadronTree->GetEntries(); ++i) 
  { 
    BhadronTree->GetEntry(i); 
    TLorentzVector b4vec = TLorentzVector( bPx, bPy, bPz, bE);
	TVector3 b3vec = TVector3( bPx, bPy, bPz);    // same as <= b4vec.Vect()
    
    TLorentzVector Kaon4vec = TLorentzVector( KaonPx, KaonPy, KaonPz,KaonE);
    TLorentzVector Mu14vec = TLorentzVector( Mu1Px, Mu1Py, Mu1Pz, Mu1E);
    TLorentzVector Mu24vec = TLorentzVector( Mu2Px, Mu2Py, Mu2Pz, Mu2E);
   

    // jet loop 
    for (int j=0; j<JetsTree->GetEntries(); ++j)
	{
	  JetsTree->GetEntry(j);
      if (eventNum != beventNum) 
      {
        continue;
      }
    
      TLorentzVector j4vec = TLorentzVector( jPx, jPy, jPz, jE); 
      TVector3 j3vec = TVector3( jPx, jPy, jPz);
      double br = sqrt( pow( j4vec.Eta() - b4vec.Eta() , 2) + pow( j4vec.Phi() - b4vec.Phi() , 2) );       //fix bump at phi 
      double Kaonr = sqrt( pow( j4vec.Eta() - Kaon4vec.Eta() , 2) + pow( j4vec.Phi() - Kaon4vec.Phi() , 2) );
      double Mu1r = sqrt( pow( j4vec.Eta() - Mu14vec.Eta() , 2) + pow( j4vec.Phi() - Mu14vec.Phi() , 2) );
      double Mu2r = sqrt( pow( j4vec.Eta() - Mu24vec.Eta() , 2) + pow( j4vec.Phi() - Mu24vec.Phi() , 2) );

      
      rH ->Fill(br);

      bool Kinjet = Kaonr <= 0.5;
      bool Mu1injet = Mu1r <= 0.5;
      bool Mu2injet = Mu2r <= 0.5;


      if (br <= 0.5 && Kaonr <= 0.5 && Mu1r <= 0.5 && Mu2r <= 0.5)  // Should we check  b isninside the jet again after adding particle in there?
      { 
        z = ( j3vec * b3vec ) / ( j3vec.Mag2() );
	    Jt = (j3vec.Cross(b3vec).Mag())/(j3vec.Mag());
        
        zH->Fill(z);
        jtH->Fill(Jt); 
      }
      else
      {  
      	if (Kinjet && Mu1injet && !Mu2injet) 
        {  
          j4vec += Mu24vec;
          j3vec = j4vec.Vect();
          z = ( j3vec * b3vec ) / ( j3vec.Mag2() );
	      Jt = (j3vec.Cross(b3vec).Mag())/(j3vec.Mag());
        
          zH->Fill(z);
          jtH->Fill(Jt);

        }
		
		else if (Kinjet && !Mu1injet && Mu2injet) 
        {
          
          j4vec += Mu14vec;
          j3vec = j4vec.Vect();
          z = ( j3vec * b3vec ) / ( j3vec.Mag2() );
	      Jt = (j3vec.Cross(b3vec).Mag())/(j3vec.Mag());
        
          zH->Fill(z);
          jtH->Fill(Jt);

		   
        }
      
        else if (!Kinjet && Mu1injet && Mu2injet) 
        {
          j4vec += Kaon4vec;
          j3vec = j4vec.Vect();
          z = ( j3vec * b3vec ) / ( j3vec.Mag2() );
	      Jt = (j3vec.Cross(b3vec).Mag())/(j3vec.Mag());
        
          zH->Fill(z);
          jtH->Fill(Jt);
        }
        
        else 
        { 
         cout<< "Problem encountered" << endl;
        }
	  } 
     
    } 
        
  }
	

	auto c=new TCanvas("Canvas","Canvas",800,800);
	c->Divide(1,3);
	c->cd(1); zH->Draw();
	c->cd(2); jtH->Draw();
	c->cd(3); rH->Draw();


 
}




