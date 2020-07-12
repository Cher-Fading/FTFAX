#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>

#include "../atlasstyle-00-04-02/AtlasUtils.h"
#include "../atlasstyle-00-04-02/AtlasStyle.h"
#include "../atlasstyle-00-04-02/AtlasLabels.h"
#include "../atlasstyle-00-04-02/AtlasStyle.C"

#ifdef __CLING__
// these are not headers - do not treat them as such - needed for ROOT6
#include "../atlasstyle-00-04-02/AtlasLabels.C"
#include "../atlasstyle-00-04-02/AtlasUtils.C"
#endif

#ifdef __CINT__
gROOT->LoadMacro("../atlasstyle-00-04-02/AtlasLabels.C");
gROOT->LoadMacro("../atlasstyle-00-04-02/AtlasUtils.C");
#endif

const Float_t FCal_range[] = {0, 0.063719, 0.14414, 0.289595, 0.525092, 0.87541, 1.36875, 2.04651, 2.98931, 5}; // fcal_cuts options

const float Weight[] = {6.7890E+07, 6.789E+07, 6.3996E+05, 4.7195E+03, 2.6602E+01, 2.2476E-01};
const float Filter[] = {9.9713E-01, 2.8748E-03, 4.2952E-03, 5.2994E-03, 4.5901E-03, 2.1846E-03};
const int grid_size = sizeof(Weight) / sizeof(float);
const int myColor[] = {kBlue, kViolet, kMagenta, kPink, kOrange, kYellow, kSpring, kTeal, kCyan, kAzure, kGray, kGray + 1, kGray + 3};
const int cet[] = {0, 2, 2, 5, 5, 8}; //selected centrality sections
const int cet_N = (sizeof(cet) / sizeof(int)) / 2;
const bool PbPb = true;
char Type[2][10] = {"pp", "PbPb"};
const char *dataType = "WorkingDefaultpp";
Float_t Eta_range[] = {-2.1, -1.5, -0.9, -0.3, 0.3, 0.9, 1.5, 2.1};
const int Eta_N = sizeof(Eta_range) / sizeof(float) - 1;

const int count_max = 1500;
const int pt_min = 50;
const int pt_max = 600;
const int pt_jf_min = 50;
const int pt_max_jf = 2000;
const int pt_bin = 16;
const int min_dist = 0;
const float min_dist0 = 1;
const int max_dist = 40;
const int dist_bin = 50;

const float min_distx = -0.5;
const float max_distx = -0.25;
const float min_disty = -0.95;
const float max_disty = -0.75;
const float min_distz = -200;
const float max_distz = 200;

const float min_distxs = -20;
const float max_distxs = 20;
const float min_distys = -20;
const float max_distys = 20;
const float min_distzs = -200;
const float max_distzs = 200;

const int dist_bin3 = 100;

const int bHadCut = 5;
const int cHadCut = 5;
const Float_t eta_selection = 2.1;

const int kTruth = kBlue;
const int kReco = kRed;
const char *track_selection = "No Selection";
const bool unique_B = false;

const int dist_bin2 = 20;
const float low_log = 1e-4;
const float high_log = 0.9;
const float distsize = 2;
//const float min_dist2 = 1e-5;
//const float max_dist2 = 0.3;

const bool jetTruth = true;
const bool TRK = false;

void initBranches(TChain *myChain)
{

   myChain->SetBranchStatus("*", 1);

   myChain->SetBranchStatus("jet_pt", 1);
   myChain->SetBranchStatus("jet_eta", 1);
   myChain->SetBranchStatus("jet_truthMatch", 1);
   myChain->SetBranchStatus("njets", 1);
   myChain->SetBranchStatus("eventnb", 1);
   myChain->SetBranchStatus("Fcal", 1);
   myChain->SetBranchStatus("jet_jf_nvtx", 1);
   myChain->SetBranchStatus("jet_jf_vtx_L3D", 1);
   myChain->SetBranchStatus("runnb", 1);
   myChain->SetBranchStatus("PVx", 1);
   myChain->SetBranchStatus("PVy", 1);
   myChain->SetBranchStatus("PVz", 1);
   myChain->SetBranchStatus("truth_PVx", 1);
   myChain->SetBranchStatus("truth_PVy", 1);
   myChain->SetBranchStatus("truth_PVz", 1);
   myChain->SetBranchStatus("jet_bH_pdgId", 1);
   myChain->SetBranchStatus("jet_cH_pdgId", 1);
   myChain->SetBranchStatus("jet_bH_Lxy", 1);
   myChain->SetBranchStatus("jet_cH_Lxy", 1);
   myChain->SetBranchStatus("jet_bH_x", 1);
   myChain->SetBranchStatus("jet_bH_y", 1);
   myChain->SetBranchStatus("jet_bH_z", 1);
   myChain->SetBranchStatus("jet_m", 1);
   myChain->SetBranchStatus("jet_truthMatch", 1);
   myChain->SetBranchStatus("jet_truthPt", 1);
   myChain->SetBranchStatus("jet_truthEta", 1);
   myChain->SetBranchStatus("jet_dRminToB", 1);
   myChain->SetBranchStatus("jet_dRminToC", 1);
   myChain->SetBranchStatus("jet_dRminToT", 1);
   myChain->SetBranchStatus("jet_bH_prod_x", 1);
   myChain->SetBranchStatus("jet_bH_prod_y", 1);
   myChain->SetBranchStatus("jet_bH_prod_z", 1);
   myChain->SetBranchStatus("jet_jf_vtx_x", 1);
   myChain->SetBranchStatus("jet_jf_vtx_y", 1);
   myChain->SetBranchStatus("jet_jf_vtx_z", 1);
   myChain->SetBranchStatus("closestVtx_L3D", 1);
   myChain->SetBranchStatus("jet_btag_ntrk", 1);
   myChain->SetBranchStatus("jet_trk_orig", 1);
   myChain->SetBranchStatus("jet_jf_llr", 1);
}

void bTagJF_grid()
{
   Float_t ptbins[pt_bin + 1];
   float initial = log(pt_min);
   float incre = log(pt_max / pt_min) / pt_bin;
   for (int i = 0; i < (pt_bin + 1); i++)
   {
      ptbins[i] = TMath::Power(TMath::E(), initial);
      initial = initial + incre;
   }

   Float_t distbins[pt_bin + 1];
   float initialdist = log(min_dist0);
   float incredist = log(max_dist / min_dist0) / pt_bin;
   for (int i = 0; i < (pt_bin + 1); i++)
   {
      distbins[i] = TMath::Power(TMath::E(), initialdist);
      initialdist = initialdist + incredist;
   }
   std::string chain_name = "bTag_AntiKt4HIJets";
   TChain *myChain = new TChain(chain_name.c_str());

   std::ifstream file(Form("/usatlas/u/cher97/GetStuff/%s.txt", dataType));
   std::string line;
   while (std::getline(file, line))
   {
      std::stringstream linestream(line);
      std::string item;
      int linePos = 0;
      std::string fileName;
      while (std::getline(linestream, item, '|'))
      {

         //	  std::cout <<  item << " linePos " << linePos << endl;
         if (linePos == 5)
         {
            if (item.find("REPLICA") != std::string::npos)
               continue;
            std::string::iterator end_pos = std::remove(item.begin(), item.end(), ' ');
            item.erase(end_pos, item.end());
            //cout << end_pos << endl;
            //cout << item << endl;
            int start_pos = item.find_last_of("=");
            //cout << start_pos << endl;
            fileName = item.substr(start_pos + 1, item.length() - start_pos - 1);
            //cout << fileName << endl;

            if (fileName.find(".root") == std::string::npos)
            {
               cout << fileName << "file missing" << endl;
               return;
            }
            if (fileName.find("user.xiaoning") == std::string::npos)
            {
               cout << fileName << "file missing" << endl;
               break;
            }
            //cout << fileName << endl;
            myChain->Add(fileName.c_str());
         }
         ++linePos;
      }
   }

   file.close();

   int JZ_ID[grid_size];

   std::ifstream filej("/usatlas/u/cher97/GetStuff/JZ_ID.txt");
   std::string linej;
   while (std::getline(filej, linej))
   {
      std::stringstream linestreamj(linej);
      std::string itemj;
      int linePosj = 0;
      std::string id;
      cout << linej << endl;
      while (std::getline(linestreamj, itemj, ' '))
      {
         if (itemj == "")
            continue;

         if (linePosj == 0)
         {
            id = itemj;
            cout << "id: " << id << endl;
         }
         if (linePosj == 4)
         {
            if (itemj.find(dataType) == std::string::npos)
            {
               cout << "Wrong file name" << itemj << endl;
               goto here;
            }
            int k = itemj.find("JZ");
            cout << itemj << endl;
            if (k == std::string::npos)
            {
               cout << "Wrong name" << itemj << endl;
               break;
            }
            JZ_ID[itemj[k + 2] - 48] = std::stoi(id);
            cout << itemj[k + 2] - 48 << ": " << id << endl;
         }
         ++linePosj;
      }
   }

here:
   std::cout << "Chain Entries:" << myChain->GetEntries() << std::endl;
   initBranches(myChain);

   std::vector<float> *jet_pt = 0;
   std::vector<float> *jet_eta = 0;
   std::vector<int> *jet_truthMatch = 0;
   Int_t njets = 0;
   Int_t eventnb = 0;
   Int_t runnb = 0;
   Float_t Fcal = 0.;
   Float_t truth_PVx;
   Float_t truth_PVy;
   Float_t truth_PVz;
   Float_t PVx;
   Float_t PVy;
   Float_t PVz;
   std::vector<float> *jet_truthPt = 0;
   std::vector<float> *jet_truthEta = 0;
   std::vector<float> *jet_m = 0;
   std::vector<float> *jet_dRminToT = 0;
   std::vector<float> *jet_dRminToC = 0;
   std::vector<float> *jet_dRminToB = 0;
   std::vector<std::vector<float>> *jet_bH_prod_x = 0;
   std::vector<std::vector<float>> *jet_bH_prod_y = 0;
   std::vector<std::vector<float>> *jet_bH_prod_z = 0;
   std::vector<std::vector<float>> *jet_bH_Lxy = 0;
   std::vector<std::vector<float>> *jet_cH_Lxy = 0;
   std::vector<std::vector<int>> *jet_bH_pdgId = 0;
   std::vector<std::vector<int>> *jet_cH_pdgId = 0;
   std::vector<float> *jet_jf_m = 0;
   std::vector<std::vector<float>> *jet_bH_x = 0;
   std::vector<std::vector<float>> *jet_bH_y = 0;
   std::vector<std::vector<float>> *jet_bH_z = 0;
   std::vector<float> *jet_jf_nvtx = 0;
   std::vector<std::vector<float>> *jet_jf_vtx_L3D = 0;
   std::vector<float> *closestVtx_L3D = 0;
   std::vector<std::vector<float>> *jet_jf_vtx_x = 0;
   std::vector<std::vector<float>> *jet_jf_vtx_y = 0;
   std::vector<std::vector<float>> *jet_jf_vtx_z = 0;
   std::vector<int> *jet_btag_ntrk = 0;
   std::vector<std::vector<int>> *jet_trk_orig = 0;
   std::vector<float> *jet_jf_llr = 0;

   TBranch *b_jet_pt, *b_jet_eta, *b_runnb, *b_PVx, *b_PVy, *b_PVz, *b_jet_bH_pdgId, *b_jet_cH_pdgId, *b_jet_bH_Lxy, *b_jet_cH_Lxy, *b_jet_m, *b_jet_truthMatch, *b_njets, *b_eventnb, *b_Fcal, *b_truth_PVx, *b_truth_PVy, *b_truth_PVz, *b_jet_truthEta, *b_jet_truthPt, *b_jet_dRminToB, *b_jet_dRminToC, *b_jet_dRminToT, *b_jet_bH_prod_x, *b_jet_bH_prod_y, *b_jet_bH_prod_z, *b_jet_jf_m, *b_jet_bH_x, *b_jet_bH_y, *b_jet_bH_z, *b_jet_jf_nvtx, *b_jet_jf_vtx_L3D, *b_closestVtx_L3D, *b_jet_jf_vtx_x, *b_jet_jf_vtx_y, *b_jet_jf_vtx_z, *b_jet_btag_ntrk, *b_jet_trk_orig, *b_jet_jf_llr;

   /*TBranch* b_mcwg;
	Float_t mcwg = 0;
	myChain->SetBranchAddress("mcwg", &mcwg, &b_mcwg);*/

   myChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   myChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   myChain->SetBranchAddress("jet_jf_m", &jet_jf_m, &b_jet_jf_m);
   myChain->SetBranchAddress("jet_truthMatch", &jet_truthMatch, &b_jet_truthMatch);
   myChain->SetBranchAddress("njets", &njets, &b_njets);
   myChain->SetBranchAddress("eventnb", &eventnb, &b_eventnb);
   myChain->SetBranchAddress("Fcal", &Fcal, &b_Fcal);
   myChain->SetBranchAddress("jet_bH_pdgId", &jet_bH_pdgId, &b_jet_bH_pdgId);
   myChain->SetBranchAddress("jet_cH_pdgId", &jet_cH_pdgId, &b_jet_cH_pdgId);
   myChain->SetBranchAddress("jet_bH_Lxy", &jet_bH_Lxy, &b_jet_bH_Lxy);
   myChain->SetBranchAddress("jet_cH_Lxy", &jet_cH_Lxy, &b_jet_cH_Lxy);
   myChain->SetBranchAddress("jet_bH_prod_x", &jet_bH_prod_x, &b_jet_bH_prod_x);
   myChain->SetBranchAddress("jet_bH_prod_y", &jet_bH_prod_y, &b_jet_bH_prod_y);
   myChain->SetBranchAddress("jet_bH_prod_z", &jet_bH_prod_z, &b_jet_bH_prod_z);
   myChain->SetBranchAddress("jet_jf_nvtx", &jet_jf_nvtx, &b_jet_jf_nvtx);
   myChain->SetBranchAddress("jet_jf_vtx_L3D", &jet_jf_vtx_L3D, &b_jet_jf_vtx_L3D);
   myChain->SetBranchAddress("runnb", &runnb, &b_runnb);
   myChain->SetBranchAddress("jet_jf_llr", &jet_jf_llr, &b_jet_jf_llr);
   myChain->SetBranchAddress("PVx", &PVx, &b_PVx);
   myChain->SetBranchAddress("PVy", &PVy, &b_PVy);
   myChain->SetBranchAddress("PVz", &PVz, &b_PVz);
   myChain->SetBranchAddress("truth_PVx", &truth_PVx, &b_truth_PVx);
   myChain->SetBranchAddress("truth_PVy", &truth_PVy, &b_truth_PVy);
   myChain->SetBranchAddress("truth_PVz", &truth_PVz, &b_truth_PVz);
   myChain->SetBranchAddress("jet_bH_x", &jet_bH_x, &b_jet_bH_x);
   myChain->SetBranchAddress("jet_bH_y", &jet_bH_y, &b_jet_bH_y);
   myChain->SetBranchAddress("jet_bH_z", &jet_bH_z, &b_jet_bH_z);
   myChain->SetBranchAddress("jet_m", &jet_m, &b_jet_m);
   myChain->SetBranchAddress("jet_truthPt", &jet_truthPt, &b_jet_truthPt);
   myChain->SetBranchAddress("jet_truthEta", &jet_truthEta, &b_jet_truthEta);
   myChain->SetBranchAddress("jet_dRminToB", &jet_dRminToB, &b_jet_dRminToB);
   myChain->SetBranchAddress("jet_dRminToC", &jet_dRminToC, &b_jet_dRminToC);
   myChain->SetBranchAddress("jet_dRminToT", &jet_dRminToT, &b_jet_dRminToT);
   myChain->SetBranchAddress("jet_jf_vtx_x", &jet_jf_vtx_x, &b_jet_jf_vtx_x);
   myChain->SetBranchAddress("jet_jf_vtx_y", &jet_jf_vtx_y, &b_jet_jf_vtx_y);
   myChain->SetBranchAddress("jet_jf_vtx_z", &jet_jf_vtx_z, &b_jet_jf_vtx_z);
   myChain->SetBranchAddress("closestVtx_L3D", &closestVtx_L3D, &b_closestVtx_L3D);
   myChain->SetBranchAddress("jet_btag_ntrk", &jet_btag_ntrk, &b_jet_btag_ntrk);
   myChain->SetBranchAddress("jet_trk_orig", &jet_trk_orig, &b_jet_trk_orig);

   if (myChain == 0)
      return;

   const char *TorR = jetTruth ? "truth" : "reco";
   const int cent_N = PbPb ? cet_N : 1;

   Long64_t nentries = myChain->GetEntries();
   std::string fname;
   float weight;

   TFile *out = TFile::Open(Form("%s%srapidityJF%.1f%d.root", Type[PbPb], dataType, eta_selection, pt_min), "RECREATE");

   Long64_t nbytes = 0, nb = 0;

   //JFV resolution
   TH1F *SV_resolution_b[cent_N];
   TH1F *JFV_x[cent_N];
   TH1F *JFV_y[cent_N];
   TH1F *JFV_z[cent_N];

   TH1F *JFV_truth_x[cent_N];
   TH1F *JFV_truth_y[cent_N];
   TH1F *JFV_truth_z[cent_N];

   //JFV
   TH1F *reco_jf_b[cent_N];
   TH1F *all_jf_b[cent_N];
   TH1F *reco_jf_c[cent_N];
   TH1F *all_jf_c[cent_N];
   TH1F *reco_jf_l[cent_N];
   TH1F *all_jf_l[cent_N];

   //JF intergated efficiency:
   Float_t recoJFb[cent_N];
   Float_t allJFb[cent_N];
   Float_t recoJFc[cent_N];
   Float_t allJFc[cent_N];
   Float_t recoJFl[cent_N];
   Float_t allJFl[cent_N];

   Float_t secJFt[cent_N];
   Float_t recoJFt[cent_N];
   Float_t truthJFt[cent_N];

   Float_t recoJFt_cate[cent_N][5];
   Float_t secJFt_cate[cent_N][5];

   //JF track purity
   TH1F *reco_trk_b[cent_N];
   TH1F *reco_trk_d[cent_N];
   TH1F *reco_trk_f[cent_N];
   /*TH1F* truth_trk_b[cent_N];
   TH1F* truth_trk_d[cent_N];
   TH1F* truth_trk_f[cent_N];*/
   TH1F *reco_b_truth_b[cent_N];
   TH1F *reco_b_truth_d[cent_N];
   TH1F *reco_b_truth_f[cent_N];
   TH1F *reco_d_truth_b[cent_N];
   TH1F *reco_d_truth_d[cent_N];
   TH1F *reco_d_truth_f[cent_N];

   for (int i = 0; i < cent_N; i++)
   {

      SV_resolution_b[i] = hotTH1F(Form("SV_resolution_b_cent_%d", i), "Distance between Truth and JF Reco Secondary Vertices", dist_bin2, 0, 80, "", "", myColor[i * 2], 0.3, 21, 1, true);
      JFV_x[i] = hotTH1F(Form("JFV_x_cent_%d", i), "JF Secondary Vertex x", dist_bin2, min_distxs, max_distxs, "", "", kReco, 0.3, 21, 1, true);
      JFV_y[i] = hotTH1F(Form("JFV_y_cent_%d", i), "JF Secondary Vertex y", dist_bin2, min_distys, max_distys, "", "", kReco, 0.3, 21, 1, true);
      JFV_z[i] = hotTH1F(Form("JFV_z_cent_%d", i), "JF Secondary Vertex z", dist_bin2, min_distzs, max_distzs, "", "", kReco, 0.3, 21, 1, true);

      JFV_truth_x[i] = hotTH1F(Form("JFV_truth_x_cent_%d", i), "JF Truth Secondary Vertex x", dist_bin2, min_distxs, max_distxs, "", "", kTruth, 0.3, 21, 1, true);
      JFV_truth_y[i] = hotTH1F(Form("JFV_truth_y_cent_%d", i), "JF Truth Secondary Vertex y", dist_bin2, min_distys, max_distys, "", "", kTruth, 0.3, 21, 1, true);
      JFV_truth_z[i] = hotTH1F(Form("JFV_truth_z_cent_%d", i), "JF Truth Secondary Vertex z", dist_bin2, min_distzs, max_distzs, "", "", kTruth, 0.3, 21, 1, true);

      //dL3d_b[i][j] = hotTH1F(Form("dL3d_b_cent_%d_pT_%d",i,j), "Distribution of L3d difference between Truth and Reco in b-jet", dist_bin, min_dist, max_dist,"","",myColor[j],0.3,21,1,true);

      //JetFitter secondary vertices
      reco_jf_b[i] = hotTH1F(Form("reco_jf_b_cent_%d", i), "Number of jets with JFV reconstructed for b jet with Truth Match versus Truth Jet pt", pt_bin, ptbins, "", "", kBlack, 0.3, 21, 1, true);
      all_jf_b[i] = hotTH1F(Form("all_jf_b_cent_%d", i), "Number of b jet with with Truth Match versus Truth Jet pt", pt_bin, ptbins, "", "", kBlack, 0.3, 21, 1, true);

      reco_jf_c[i] = hotTH1F(Form("reco_jf_c_cent_%d", i), "Number of jets with JFV reconstructed for c jet with Truth Match versus Truth Jet pt", pt_bin, ptbins, "", "", kBlack, 0.3, 21, 1, true);
      all_jf_c[i] = hotTH1F(Form("all_jf_c_cent_%d", i), "Number of c jet with with Truth Match versus Truth Jet pt", pt_bin, ptbins, "", "", kBlack, 0.3, 21, 1, true);

      reco_jf_l[i] = hotTH1F(Form("reco_jf_l_cent_%d", i), "Number of jets with JFV reconstructed for light jet (fake SV) with Truth Match versus Truth Jet pt", pt_bin, ptbins, "", "", kBlack, 0.3, 21, 1, true);
      all_jf_l[i] = hotTH1F(Form("all_jf_l_cent_%d", i), "Number of light jet with with Truth Match versus Truth Jet pt", pt_bin, ptbins, "", "", kBlack, 0.3, 21, 1, true);
   }

   int multiB = 0;
   int nJets = 0;
   int NJets = 0;

   std::vector<float> wgsum;

   for (int jz = 0; jz < grid_size; jz++)
   {
      wgsum.push_back(0);
   }
   //std::vector<int> order;
   for (Long64_t jentry0 = 0; jentry0 < nentries; jentry0++)
   {
      Long64_t ientry0 = myChain->LoadTree(jentry0);
      //b_mcwg->GetEntry(ientry0);
      if (ientry0 < 0)
      {
         cout << "invalid index" << endl;
         break;
      }

      //if (jentry0 == 100000)
      //   break;

      fname = myChain->GetCurrentFile()->GetName();
      //cout << fname << endl;
      int k = fname.find("Akt4HIJets");
      if (k == std::string::npos)
      {
         cout << "Wrong name" << fname << endl;
         break;
      }
      bool found = false;
      //cout << stoi(fname.substr(k - 9, 8)) << endl;
      for (int j = 0; j < grid_size; j++)
      {
         if (stoi(fname.substr(k - 9, 8)) == JZ_ID[j])
         {
            found = true;
            wgsum[j] = wgsum[j] + 1;
         }
      }
      if (!found)
      {
         cout << fname << endl;
         cout << "not found" << endl;
         return;
      }
   }
   cout << "number of slices:" << wgsum.size() << endl;
   std::ofstream outfile;
   outfile.open(Form("%s.txt",dataType), std::ios_base::trunc);
   for (int i = 0; i < wgsum.size(); i++)
   {
      cout << "JZ Slice " << i << ": " <<  std::setprecision(9) << wgsum[i] << endl;
      outfile << i << ": " << std::setprecision(9) << wgsum[i] << endl;
   }
return;
   //loop over for weight
   //float min_dist3d = 10;
   //TH1F* l3d_truth = new TH1F("l3d_truth","l3d_truth");
   for (Long64_t jentry = 0; jentry < nentries; jentry++)
   {
      Long64_t ientry = myChain->LoadTree(jentry);
      if (ientry < 0)
         break;

      if (ientry < 0)
         break;
      if (jentry % 10000 == 0)
         std::cout << jentry << std::endl;

      //if (jentry == 10000)
      //   break;

      b_jet_pt->GetEntry(ientry);
      b_jet_eta->GetEntry(ientry);
      b_jet_truthMatch->GetEntry(ientry);
      b_njets->GetEntry(ientry);
      b_runnb->GetEntry(ientry);
      b_eventnb->GetEntry(ientry);
      b_Fcal->GetEntry(ientry);
      b_truth_PVx->GetEntry(ientry);
      b_truth_PVy->GetEntry(ientry);
      b_truth_PVz->GetEntry(ientry);
      b_PVx->GetEntry(ientry);
      b_PVy->GetEntry(ientry);
      b_PVz->GetEntry(ientry);
      b_jet_bH_pdgId->GetEntry(ientry);
      b_jet_cH_pdgId->GetEntry(ientry);
      b_jet_bH_Lxy->GetEntry(ientry);
      b_jet_cH_Lxy->GetEntry(ientry);
      b_jet_jf_m->GetEntry(ientry);
      b_jet_truthEta->GetEntry(ientry);
      b_jet_truthPt->GetEntry(ientry);
      b_jet_dRminToB->GetEntry(ientry);
      b_jet_dRminToC->GetEntry(ientry);
      b_jet_dRminToT->GetEntry(ientry);
      b_jet_bH_prod_x->GetEntry(ientry);
      b_jet_bH_prod_y->GetEntry(ientry);
      b_jet_bH_prod_z->GetEntry(ientry);
      b_jet_m->GetEntry(ientry);
      b_jet_bH_x->GetEntry(ientry);
      b_jet_bH_y->GetEntry(ientry);
      b_jet_bH_z->GetEntry(ientry);
      b_jet_jf_nvtx->GetEntry(ientry);
      b_jet_jf_vtx_L3D->GetEntry(ientry);
      b_closestVtx_L3D->GetEntry(ientry);
      b_jet_jf_vtx_x->GetEntry(ientry);
      b_jet_jf_vtx_y->GetEntry(ientry);
      b_jet_jf_vtx_z->GetEntry(ientry);
      b_jet_btag_ntrk->GetEntry(ientry);
      b_jet_trk_orig->GetEntry(ientry);
      b_jet_jf_llr->GetEntry(ientry);
      //
      //cout << "line 564" << endl;
      int pos;

      if (ientry == 0){
      fname = myChain->GetCurrentFile()->GetName();

      int k = fname.find("Akt4HIJets");
      if (k == std::string::npos)
      {
         cout << "Wrong name" << fname << endl;
         break;
      }
      bool found = false;
      for (int j = 0; j < grid_size; j++)
      {
         if (stoi(fname.substr(k - 9, 8)) == JZ_ID[j])
         {
            found = true;
            weight = Weight[j] * Filter[j] / wgsum[j];
         }
      }
      if (!found)
      {
         cout << fname << endl;
         cout << "not found" << endl;
         break;
      }
      }
      if (weight == 0)
      {
         cout << "wrong weight" << endl;
         break;
      }

      float FCal_et = Fcal;
      //cout << "FCal_et: " << FCal_et << endl;
      int central = 0;
      if (PbPb)
      {
         for (int f = 0; f < cet_N; f++)
         {
            if (FCal_et > FCal_range[9 - cet[f * 2 + 1]] && FCal_et < FCal_range[9 - cet[f * 2]])
               central = f;
         }
      }
      //cout << "fcal:" << central << Fcal << endl;

      for (int j = 0; j < jet_truthMatch->size(); j++)
      {
         if (jetTruth && jet_truthMatch->at(j) != 1)
            continue;
         float jetPt = jetTruth ? (jet_truthPt->at(j) * 1e-3) : (jet_pt->at(j) * 1e-3);
         //cout << "jetpt " << jetPt << endl;
         if (jetPt < pt_min && jetPt > pt_max)
            continue;
         float jetEta = jetTruth ? (jet_truthEta->at(j)) : (jet_eta->at(j));
         int jetType = 0; //1 for B jet, 2 for C jet, 3 for T jet, 0 for light jet

         if (fabs(jetEta) > eta_selection)
            continue;

         if (jet_dRminToT->at(j) < 0.3)
            jetType = 3; //b hadron particle pt cut 3GeV
         if (jet_dRminToC->at(j) < 0.3)
            jetType = 2;
         if (jet_dRminToB->at(j) < 0.3)
            jetType = 1;

         //if (jet_cH_pdgId->at(j)[0]!=-99) jetType = 2;
         //if (jet_bH_pdgId->at(j)[0]!=-99) jetType = 1;

         if (jetType == 1)
         {

            NJets++;
            if (unique_B && jet_bH_pdgId->at(j).size() > 1)
            {
               //cout << "here" << endl;
               multiB++;
               continue;
            }
            nJets++;

            //JF
            //cout << jetPt << weight << endl;
            all_jf_b[central]->Fill(jetPt, weight);
            if (jet_jf_llr->at(j) != -99)
               reco_jf_b[central]->Fill(jetPt, weight);
            //return;
            /*if (jetPt >= pt_min && jetPt <= pt_max)
            {
               if (jetPt >= pt_jf_min)
                  allJFb[central] = allJFb[central] + weight;
               if (jet_jf_llr->at(j) != -99)
               {
                  reco_jf_b[central]->Fill(jetPt, weight);
                  if (jetPt >= pt_jf_min)
                     recoJFb[central] = recoJFb[central] + weight;
               }
            }*/
         }

         if (jetType == 2)
         {
            //JF
            all_jf_c[central]->Fill(jetPt, weight);
            if (jet_jf_llr->at(j) != -99)
               reco_jf_c[central]->Fill(jetPt, weight);
            /*if (jetPt >= pt_min && jetPt <= pt_max)
               allJFc[central] = allJFc[central] + weight;
            if (jet_jf_llr->at(j) != -99)
            {
               reco_jf_c[central]->Fill(jetPt, weight);
               if (jetPt >= pt_min && jetPt <= pt_max)
                  recoJFc[central] = recoJFc[central] + weight;
            }*/
         }

         if (jetType == 0)
         {
            //JF
            all_jf_l[central]->Fill(jetPt, weight);
            if (jet_jf_llr->at(j) != -99)
               reco_jf_l[central]->Fill(jetPt, weight);
            /*if (jetPt >= 20)
               allJFl[central] = allJFl[central] + weight;
            if (jet_jf_llr->at(j) != -99)
            {
               reco_jf_l[central]->Fill(jetPt, weight);
               if (jetPt >= 20)
                  recoJFl[central] = recoJFl[central] + weight;
            }*/
         }

         jetType = 0;
         //if (jetPt<=pt_min || jetPt>=pt_max) continue;
         if (jet_dRminToT->at(j) < 0.3)
            jetType = 3;
         //if (jet_nGhostTau) jetType = 3;
         if (jet_cH_pdgId->at(j)[0] != -99)
            jetType = 2;
         if (jet_bH_pdgId->at(j)[0] != -99)
            jetType = 1;
         //if (jet_dRminToB->at(j) < 0.3) jetType = 1;
         //int nbhad_t = 0;//number of truth hadrons matched to each jet
         //int nchad_t = 0;
         //use the closest truth hadron matched to the jet

         if (jetType == 1)
         {
            //if (jet_bH_pt->at(j)[0]*1e-3 < bHadCut) continue;
            //if ((jet_pt->at(j)/jet_truthPt->at(j))>2 || (jet_truthPt->at(j)/jet_pt->at(j))>2) continue;
            //int pT = (int)(TMath::Log(jetPt/pt_min)/incre);

            if (jet_bH_Lxy->at(j)[0] <= -99)
            {
               cout << "出大问题 (big problem) " << endl;
               cout << ientry << endl;
               cout << fname << endl;
               cout << j << endl;
               cout << jet_bH_pdgId->at(j)[0] << endl;
               return;
            }
            if (!(jet_bH_Lxy->at(j)[0] >= 0))
            {
               cout << "Has no decay vertex" << endl;
               continue;
            }
            //exclude b hadrons with production vertex not at truth PV.
            bool truthB = (jet_bH_prod_x->at(j)[0] == truth_PVx && jet_bH_prod_y->at(j)[0] == truth_PVy && jet_bH_prod_z->at(j)[0] == truth_PVz);
            bool recoB = (jet_jf_m->at(j) > 0);

            float truth_L3d;
            float Diff_3d;
            float Diff;

            if (truthB)
            {
               JFV_truth_x[central]->Fill(jet_bH_x->at(j)[0], weight);
               JFV_truth_y[central]->Fill(jet_bH_y->at(j)[0], weight);
               JFV_truth_z[central]->Fill(jet_bH_z->at(j)[0], weight);
            }

            int jfsv = -1;

            if (recoB)
            {
               for (int v = 0; v < jet_jf_nvtx->at(j); v++)
               {
                  if (jet_jf_vtx_L3D->at(j)[v] == closestVtx_L3D->at(j))
                  {
                     jfsv = v;
                     JFV_x[central]->Fill(jet_jf_vtx_x->at(j)[v], weight);
                     JFV_y[central]->Fill(jet_jf_vtx_y->at(j)[v], weight);
                     JFV_z[central]->Fill(jet_jf_vtx_z->at(j)[v], weight);
                  }
               }
            }

            if (recoB && truthB)
            {
               truth_L3d = jet_jf_vtx_L3D->at(j)[jfsv];
               Diff_3d = sqrt(pow(jet_jf_vtx_x->at(j)[jfsv] - jet_bH_x->at(j)[0], 2) + pow(jet_jf_vtx_y->at(j)[jfsv] - jet_bH_y->at(j)[0], 2) + pow(jet_jf_vtx_z->at(j)[jfsv] - jet_bH_z->at(j)[0], 2));
               SV_resolution_b[central]->Fill(Diff_3d, weight);
            }
         }

         if (jetType == 2)
         {
            if (jet_cH_Lxy->at(j)[0] <= -99)
            {
               cout << "出大问题 (big problem) " << endl;
               continue;
            }
            if (!(jet_cH_Lxy->at(j)[0] >= 0))
            {
               cout << "Has no decay vertex" << endl;
               continue;
            }
         }
      }
   }

   //JFV resolution

   TCanvas *cjf3 = new TCanvas("cjf3", "cjf3", 500, 500);
   TPad *p0 = (TPad *)cjf3->cd();
   TH1F *h0 = (TH1F *)p0->DrawFrame(0, 1e-4, 80, 1);
   h0->GetXaxis()->SetTitle("Distance (Truth-Reco) (mm)");
   h0->GetXaxis()->SetTitle(Form("Normalized Fraction/%.2f mm", 40. / dist_bin2));
   h0->SetTitle("Distance between Truth and JF Reco Secondary Vertices");

   for (int i = 0; i < cent_N; i++)
   {
      const char *Centrality = PbPb ? Form(" %d %% %d %%", 10 * cet[2 * i], 10 * cet[2 * i + 1]) : "";
      SV_resolution_b[i]->Scale(1. / SV_resolution_b[i]->GetSumOfWeights());
      SV_resolution_b[i]->Draw("SAME");
      SV_resolution_b[i]->Write();
      myBoxText(0.4, 0.85 - i * 0.05, 0.05, myColor[2 * i], 0, Form("%s%s", Type[PbPb], Centrality), 1, myColor[2 * i], 21, true);
   }
   gPad->SetGrid(1);
   gPad->SetTicks(1);
   cjf3->SetLogy();
   cjf3->SaveAs(Form("JFV reco truth distance %s %s rapidity %.1f.pdf", Type[PbPb], dataType, eta_selection));

   TCanvas *cjf4 = new TCanvas("cjf4", "cjf4", 1500, 500);
   cjf4->Divide(3, 1);
   TPad *px;
   TPad *py;
   TPad *pz;

   //cout << min_distxs << " " << max_distxs << endl;
   //cout << min_distys << " " << max_distys << endl;
   //cout << min_distzs << " " << max_distzs << endl;

   TH1F *hx;
   TH1F *hy;
   TH1F *hz;

   for (int i = 0; i < cent_N; i++)
   {
      //gPad->SetGrid(1);
      //gPad->SetTicks(1);
      const char *Centrality = PbPb ? Form(" %d %% %d %%", 10 * cet[2 * i], 10 * cet[2 * i + 1]) : "";

      px = (TPad *)cjf4->cd(1);
      hx = (TH1F *)px->DrawFrame(min_distxs, 1e-4, max_distxs, 0.7);
      hx->GetXaxis()->SetTitle("SV x coord (mm)");
      hx->GetYaxis()->SetTitle(Form("Nornalized Fraction/%.2f mm", (max_distxs - min_distxs) / dist_bin2));
      hx->SetTitle("x Coordinate of JF Secondary Vertex");

      cjf4->cd(1);
      hx->Draw();
      JFV_x[i]->Scale(1. / JFV_x[i]->GetSumOfWeights());
      JFV_x[i]->Draw("SAME");
      JFV_x[i]->Write();
      myBoxText(0.4, 0.85, 0.05, kReco, 0, Form("JF SV: %s%s", Type[PbPb], Centrality), 1, kReco, 21, true);
      JFV_truth_x[i]->Scale(1. / JFV_truth_x[i]->GetSumOfWeights());
      JFV_truth_x[i]->Draw("SAME");
      JFV_truth_x[i]->Write();
      myBoxText(0.4, 0.80, 0.05, kTruth, 0, Form("Truth SV: %s%s", Type[PbPb], Centrality), 1, kTruth, 21, true);

      py = (TPad *)cjf4->cd(2);
      hy = (TH1F *)py->DrawFrame(min_distys, 1e-4, max_distys, 0.3);
      hy->GetXaxis()->SetTitle("SV y coord (mm)");
      hy->GetYaxis()->SetTitle(Form("Nornalized Fraction/%.2f mm", (max_distxs - min_distxs) / dist_bin2));
      hy->SetTitle("y Coordinate of JF Secondary Vertex");

      cjf4->cd(2);
      hy->Draw();
      JFV_y[i]->Scale(1. / JFV_y[i]->GetSumOfWeights());
      JFV_y[i]->Draw("SAME");
      JFV_y[i]->Write();
      myBoxText(0.4, 0.85, 0.05, kReco, 0, Form("JF SV: %s%s", Type[PbPb], Centrality), 1, kReco, 21, true);
      JFV_truth_y[i]->Scale(1. / JFV_truth_y[i]->GetSumOfWeights());
      JFV_truth_y[i]->Draw("SAME");
      JFV_truth_y[i]->Write();
      myBoxText(0.4, 0.80, 0.05, kTruth, 0, Form("Truth SV: %s%s", Type[PbPb], Centrality), 1, kTruth, 21, true);

      pz = (TPad *)cjf4->cd(3);
      hz = (TH1F *)pz->DrawFrame(min_distzs, 1e-4, max_distzs, 0.3);
      hz->GetXaxis()->SetTitle("SV z coord (mm)");
      hz->GetYaxis()->SetTitle(Form("Nornalized Fraction/%.2f mm", (max_distxs - min_distxs) / dist_bin2));
      hz->SetTitle("z Coordinate of JF Secondary Vertex");

      cjf4->cd(3);
      hz->Draw();
      JFV_z[i]->Scale(1. / JFV_z[i]->GetSumOfWeights());
      JFV_z[i]->Draw("SAME");
      JFV_z[i]->Write();
      myBoxText(0.4, 0.85, 0.05, kReco, 0, Form("JF SV: %s%s", Type[PbPb], Centrality), 1, kReco, 21, true);
      JFV_truth_z[i]->Scale(1. / JFV_truth_z[i]->GetSumOfWeights());
      JFV_truth_z[i]->Draw("SAME");
      JFV_truth_z[i]->Write();
      myBoxText(0.4, 0.80, 0.05, kTruth, 0, Form("Truth SV: %s%s", Type[PbPb], Centrality), 1, kTruth, 21, true);

      /*px->SetLogy();
      py->SetLogy();
      pz->SetLogy();*/

      //gPad->SetGrid(1);
      //gPad->SetTicks(1);

      cjf4->SaveAs(Form("JFV reco and truth resolution %s%s %s rapidity %.1f_morehist.pdf", Type[PbPb], Centrality, dataType, eta_selection));
   }

   //JF Efficiency vs pt
   TCanvas *cjf0 = new TCanvas("cjf0", "cjf0", 500, 500);
   TGraphAsymmErrors *eff_jf_b[cent_N];
   TGraphAsymmErrors *eff_jf_c[cent_N];
   TGraphAsymmErrors *eff_jf_l[cent_N];
   TMultiGraph *hjf0[cent_N];
   gPad->SetTicks(1);

   for (int i = 0; i < cent_N; i++)
   {
      eff_jf_b[i] = new TGraphAsymmErrors(reco_jf_b[i], all_jf_b[i], "cl=0.683 b(1,1) mode");
      eff_jf_b[i]->SetMarkerSize(1);
      eff_jf_b[i]->SetMarkerColor(kBlue);
      eff_jf_b[i]->SetMarkerStyle(21);
      eff_jf_b[i]->SetLineColor(kBlue);
      eff_jf_b[i]->SetLineWidth(1);
      eff_jf_c[i] = new TGraphAsymmErrors(reco_jf_c[i], all_jf_c[i], "cl=0.683 b(1,1) mode");
      eff_jf_c[i]->SetMarkerSize(1);
      eff_jf_c[i]->SetMarkerColor(kRed);
      eff_jf_c[i]->SetMarkerStyle(21);
      eff_jf_c[i]->SetLineColor(kRed);
      eff_jf_c[i]->SetLineWidth(1);
      eff_jf_c[i]->SetLineStyle(10);
      eff_jf_l[i] = new TGraphAsymmErrors(reco_jf_l[i], all_jf_l[i], "cl=0.683 b(1,1) mode");
      eff_jf_l[i]->SetMarkerSize(1);
      eff_jf_l[i]->SetMarkerColor(418);
      eff_jf_l[i]->SetMarkerStyle(21);
      eff_jf_l[i]->SetLineColor(418);
      eff_jf_l[i]->SetLineWidth(1);
      eff_jf_l[i]->SetLineStyle(9);
      hjf0[i] = new TMultiGraph();
      //cout << a << endl;

      hjf0[i]->Add(eff_jf_b[i]);
      hjf0[i]->Add(eff_jf_c[i]);
      hjf0[i]->Add(eff_jf_l[i]);
      hjf0[i]->Draw("APE");

      hjf0[i]->GetXaxis()->SetLimits(pt_min, pt_max);
      hjf0[i]->GetYaxis()->SetLimits(0, 1.4);
      hjf0[i]->GetXaxis()->SetRangeUser(pt_min, pt_max);
      hjf0[i]->GetYaxis()->SetRangeUser(0, 1.4);
      hjf0[i]->GetXaxis()->SetTitle(Form("p_{T}^{%s} (GeV)", TorR));
      hjf0[i]->GetYaxis()->SetTitle("JFV Reco Efficiency");
      hjf0[i]->SetTitle(Form("JFV Reconstruction Efficiency in %s MC", Type[PbPb]));
      reco_jf_b[i]->Write();
      all_jf_b[i]->Write();
      eff_jf_b[i]->Write();

      reco_jf_c[i]->Write();
      all_jf_c[i]->Write();
      eff_jf_c[i]->Write();

      reco_jf_l[i]->Write();
      all_jf_l[i]->Write();
      eff_jf_l[i]->Write();

      myBoxText(0.4, 0.85, 0.05, kBlue, 0, "b-jet", 1, kBlue, 21, true);
      myBoxText(0.4, 0.8, 0.05, kRed, 0, "c-jet", 10, kRed, 21, true);
      myBoxText(0.4, 0.75, 0.05, 418, 0, "light jet", 9, 418, 21, true);

      const char *Centrality = PbPb ? Form(" %d %% %d %%", 10 * cet[2 * i], 10 * cet[2 * i + 1]) : "";

      //c0->SetLogy();
      myText(0.6, 0.85, kBlack, "|#eta|<2.1");
      myText(0.6, 0.8, kBlack, "#sqrt{s} = 5 TeV");
      myText(0.6, 0.75, kBlack, Form("%s b#bar{b} filtered MC%s", Type[PbPb], Centrality));

      //cout << Form("Integrated Purity for %s%s: %.2f",Type[PbPb],Centrality,reco_jf_b->GetSumOfWeights()/all_jf_b->GetSumOfWeights()) << endl;
      cout << Form("Integrated Fraction for b Jet in %s%s: %.3f", Type[PbPb], Centrality, recoJFb[i] / allJFb[i]) << endl;
      cout << Form("Integrated Fraction for c Jet in %s%s: %.3f", Type[PbPb], Centrality, recoJFc[i] / allJFc[i]) << endl;
      cout << Form("Integrated Fraction for light Jet (fake) in %s%s: %.3f", Type[PbPb], Centrality, recoJFl[i] / allJFl[i]) << endl;

      cjf0->SaveAs(Form("JFV reco efficiency %s%s %s rapidity %.1f_morehist.pdf", Type[PbPb], Centrality, dataType, eta_selection));
   }

   cjf0->Close();
   cjf3->Close();
   cjf4->Close();

   out->Close();
   cout << "Number of Jets Used: " << nJets << endl;
   cout << "Number of Jets not used for having multiple B: " << multiB << endl;
   cout << "Total number of Jets looped over: " << NJets << endl;
   //cout << "Minimum distance between truth and reco SV: " << min_dist3d << endl;
}
