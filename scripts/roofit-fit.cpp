#include <ROOT/RDataFrame.hxx>
#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include <cmath>
#include <iostream>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TF1.h>
#include <TSystem.h>
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooUniform.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistFunc.h"
#include "RooRealSumPdf.h"
#include "RooParamHistFunc.h"
#include "RooHistConstraint.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h" 
#include "TCanvas.h"
#include "TPaveText.h"
#include <memory>
#include <map>
#include <algorithm>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <TGaxis.h>
#include <RooProduct.h>
#include <RooPolynomial.h>
#include <RooCrystalBall.h>
#include "RooStats/SPlot.h"
#include <RooBernstein.h>
#include <RooGamma.h>
#include <RooConstVar.h>
#include <THStack.h>
#include <TMatrixD.h>
#include <RooExponential.h>

using namespace ROOT;
using namespace ROOT::Math;
using namespace RooFit;
using namespace RooStats;


/*
--------------------------------------------
    CUTS - PRESELECTION CODE
--------------------------------------------
*/

// Helper to merge root files from CERNbox
std::vector<std::string> GetFiles(const std::string &pattern)
{
    void *dirp = gSystem->OpenDirectory(".");
    const char *entry;
    std::vector<std::string> files;
    while ((entry = gSystem->GetDirEntry(dirp))) {
        std::string fname(entry);
        if (fname.find(pattern) == 0 && fname.find("group_3.root") != std::string::npos)
            files.push_back(fname);
    }
    std::sort(files.begin(), files.end());
    return files;
}

void data_merge_Preselection() {
    std::vector<std::string> files = GetFiles("file_13703527_");
    ROOT::RDataFrame df("Tuple_SpruceSLB_B0ToDpMuNu_DpToKPiPi_WS", files);
    
    auto df_cut = df
        // Flight direction
        .Filter("cand_Bz_OWNPV_FD_Z > 2", "Flight direction Z")

        // B0 and D+ vertex quality
        .Filter("cand_Bz_ENDVERTEX_CHI2DOF < 4", "B0 vertex chi2/ndof < 4")
        .Filter("cand_Bz_Dp_ENDVERTEX_CHI2DOF < 4", "D vertex chi2/ndof < 4")

        // K+ track quality and PID
        .Filter("cand_Dp_Kp_TRACKCHI2DOF < 2.5", "K+ track chi2/ndof < 2.5")
        .Filter("cand_Dp_Kp_OWNPV_IPCHI2 > 10", "K+ IP chi2 > 10")
        .Filter("cand_Dp_Kp_TRACKGHOSTPROB < 0.3", "K+ ghost prob < 0.3")
        .Filter("cand_Dp_Kp_PROBNN_K > 0.3", "K+ probNN(Kp) > 0.3")

        // pi1
        .Filter("cand_Dp_pim1_TRACKCHI2DOF < 2.5", "pi1 track chi2/ndof < 2.5")
        .Filter("cand_Dp_pim1_OWNPV_IPCHI2 > 10", "pi1 IP chi2 > 10")
        .Filter("cand_Dp_pim1_TRACKGHOSTPROB < 0.3", "pi1 ghost prob < 0.3")
        .Filter("cand_Dp_pim1_PROBNN_PI > 0.3", "pi1 probNN(pi1) > 0.3")

        // pi2
        .Filter("cand_Dp_pim2_TRACKCHI2DOF < 2.5", "pi2 track chi2/ndof < 2.5")
        .Filter("cand_Dp_pim2_OWNPV_IPCHI2 > 10", "pi2 IP chi2 > 10")
        .Filter("cand_Dp_pim2_TRACKGHOSTPROB < 0.3", "pi2 ghost prob < 0.3")
        .Filter("cand_Dp_pim2_PROBNN_PI > 0.3", "pi2 probNN(pi2) > 0.3")

        // Muon
        .Filter("cand_Bz_mum_PT > 1000", "mu pT > 1000 MeV")
        .Filter("cand_Bz_mum_P > 6000", "mu p > 6000 MeV")
        .Filter("cand_Bz_mum_TRACKCHI2DOF < 2.5", "mu track chi2/ndof < 2.5")
        .Filter("cand_Bz_mum_OWNPV_IPCHI2 > 9", "mu IP chi2 > 9")
        .Filter("cand_Bz_mum_PROBNN_MU > 0.3", "mu probNN(mu) > 0.3")
        
        // PV
        .Filter("cand_nPVs > 1", "nPVs > 1")

        // Triggers
        .Filter("cand_Bz_Hlt2SLB_B0ToDpMuNu_DpToKPiPiDecision_TOS == 1", "Hlt2 B0->D+μν TOS")
        .Filter("cand_Bz_Hlt1TwoTrackMVADecision_TOS == 1", "Hlt1 TwoTrack MVA TOS");

        // Dp Mass window
        //.Filter("cand_Bz_Dp_M > 1828 && cand_Bz_Dp_M < 1913");

    auto nTotal = *df.Count();
    auto nSel   = *df_cut.Count();
    std::cout << "Entries before cuts: " << nTotal << std::endl;
    std::cout << "Entries after cuts:  " << nSel   << std::endl;
    std::cout << "Efficiency: " << 100.0 * nSel / nTotal << " %" << std::endl;

    df_cut.Snapshot("DecayTree", "Data_1_WS_cuts.root");
}

void merge_incremental() 
{
    std::vector<std::string> newFiles = {
        "Data_ALL.root",
        "Data_3_cuts.root"
    };

    const std::string mergedFile = "Data_ALL_1.root";

    TFileMerger merger;
    merger.SetFastMethod(true);

    for (const auto &f : newFiles) {
        merger.AddFile(f.c_str());
    }

    merger.OutputFile(mergedFile.c_str());

    merger.Merge();

    std::cout << "Merged into " << mergedFile << std::endl;
}

// Helper for opening angle calculation
double OpeningAngle(const ROOT::Math::XYZVector& p1,
                    const ROOT::Math::XYZVector& p2)
{
    const double mag = p1.R() * p2.R();
    if (mag == 0.0) return -1.0;

    double c = p1.Dot(p2) / mag;
    c = std::clamp(c, -1.0, 1.0);

    return std::acos(c);
}


void new_cuts(){
    ROOT::RDataFrame df("DecayTree", "Data_1_cuts.root"); // For WS you need to manually change _mup->_mum (RS is _mup)
    
    auto df_helper = df
    .Define("p_Dp",
        [](float px, float py, float pz) {
            return ROOT::Math::XYZVector(px, py, pz);
        },
        {"cand_Bz_Dp_PX","cand_Bz_Dp_PY","cand_Bz_Dp_PZ"})
    .Define("p_Kp",
            [](float px, float py, float pz) {
                return ROOT::Math::XYZVector(px, py, pz);
            },
            {"cand_Dp_Kp_PX","cand_Dp_Kp_PY","cand_Dp_Kp_PZ"})
    .Define("p_pi1",
            [](float px, float py, float pz) {
                return ROOT::Math::XYZVector(px, py, pz);
            },
            {"cand_Dp_pim1_PX","cand_Dp_pim1_PY","cand_Dp_pim1_PZ"})
    .Define("p_pi2",
            [](float px, float py, float pz) {
                return ROOT::Math::XYZVector(px, py, pz);
            },
            {"cand_Dp_pim2_PX","cand_Dp_pim2_PY","cand_Dp_pim2_PZ"})
    .Define("p_mup",
            [](float px, float py, float pz) {
                return ROOT::Math::XYZVector(px, py, pz);
            },
            {"cand_Bz_mup_PX","cand_Bz_mup_PY","cand_Bz_mup_PZ"});

    // Composite Hadrons                
    auto df_cuts1 = df_helper
    .Filter("cand_Bz_ENDVERTEX_POS_Z-cand_Bz_Dp_ENDVERTEX_POS_Z < 0")
    .Filter("cand_Bz_ETA < 5 && cand_Bz_ETA > 2")
    .Filter("cand_Bz_OWNPV_CORRM < 10000")
    .Filter("log(cand_Bz_Dp_OWNPV_IP) > -3 && log(cand_Bz_Dp_OWNPV_IP) < 5")
    .Filter("cand_nPVs < 11")

    .Define("Delta_OpeningAngle_Dp_mup",OpeningAngle,{"p_Dp","p_mup"})
    .Define("Delta_OpeningAngle_Kp_mup",OpeningAngle,{"p_Kp","p_mup"})
    .Define("Delta_OpeningAngle_pi1_mup",OpeningAngle,{"p_pi1","p_mup"})
    .Define("Delta_OpeningAngle_pi2_mup",OpeningAngle,{"p_pi2","p_mup"})

    .Define("Delta_OpeningAngle_Dp_pi1",OpeningAngle,{"p_Dp","p_pi1"})
    .Define("Delta_OpeningAngle_Dp_pi2",OpeningAngle,{"p_Dp","p_pi2"})
    .Define("Delta_OpeningAngle_Dp_Kp",OpeningAngle,{"p_Dp","p_Kp"})
    .Define("Delta_OpeningAngle_Kp_pi1",OpeningAngle,{"p_Kp","p_pi1"})
    .Define("Delta_OpeningAngle_Kp_pi2",OpeningAngle,{"p_Kp","p_pi2"})
    .Define("Delta_OpeningAngle_pi1_pi2",OpeningAngle,{"p_pi1","p_pi2"})
        
    //.Filter("cand_Bz_mum_PT > 1500")

    .Define("cand_Bz_FD_rho", [](double svx, double svy, double svz) {
            XYZVector dir(svx, svy, svz);
            return dir.Rho();
        }, {"cand_Bz_OWNPV_FD_X", "cand_Bz_OWNPV_FD_Y", "cand_Bz_OWNPV_FD_Z"}) // radial FD
    
    .Define("cand_Bz_Dp_MASSWITHHYPOTHESIS_kkpi2",
            [](float K1_px, float K1_py, float K1_pz,
               float K2_px, float K2_py, float K2_pz,
               float pi_px, float pi_py, float pi_pz)
            {
                const double mK  = 493.677;
                const double mpi = 139.570;

                ROOT::Math::PxPyPzMVector pK1(K1_px, K1_py, K1_pz, mK);
                ROOT::Math::PxPyPzMVector pK2(K2_px, K2_py, K2_pz, mK);
                ROOT::Math::PxPyPzMVector pPi(pi_px, pi_py, pi_pz, mpi);

                auto pD = pK1 + pK2 + pPi;
                return pD.M();
            },
            {
                "cand_Dp_Kp_PX","cand_Dp_Kp_PY","cand_Dp_Kp_PZ",
                "cand_Dp_pim1_PX","cand_Dp_pim1_PY","cand_Dp_pim1_PZ",
                "cand_Dp_pim2_PX","cand_Dp_pim2_PY","cand_Dp_pim2_PZ"
            }
            )
    .Define("cand_Bz_Dp_MASSWITHHYPOTHESIS_kkpi1",
            [](float K1_px, float K1_py, float K1_pz,
               float K2_px, float K2_py, float K2_pz,
               float pi_px, float pi_py, float pi_pz)
            {
                const double mK  = 493.677;
                const double mpi = 139.570;

                ROOT::Math::PxPyPzMVector pK1(K1_px, K1_py, K1_pz, mK);
                ROOT::Math::PxPyPzMVector pK2(K2_px, K2_py, K2_pz, mK);
                ROOT::Math::PxPyPzMVector pPi(pi_px, pi_py, pi_pz, mpi);

                auto pD = pK1 + pK2 + pPi;
                return pD.M();
            },
            {
                "cand_Dp_Kp_PX","cand_Dp_Kp_PY","cand_Dp_Kp_PZ",
                "cand_Dp_pim2_PX","cand_Dp_pim2_PY","cand_Dp_pim2_PZ",
                "cand_Dp_pim1_PX","cand_Dp_pim1_PY","cand_Dp_pim1_PZ"
            }
            );

    auto df_cuts2 = df_cuts1
    .Filter("Delta_OpeningAngle_Kp_mup > 0.0005")
    .Filter("Delta_OpeningAngle_Dp_mup > 0.0005")
    .Filter("Delta_OpeningAngle_pi1_mup > 0.0005")
    .Filter("Delta_OpeningAngle_pi2_mup > 0.0005")

    .Filter("Delta_OpeningAngle_Dp_pi1 > 0.0005")
    .Filter("Delta_OpeningAngle_Dp_pi2 > 0.0005")
    .Filter("Delta_OpeningAngle_Dp_Kp > 0.0005")
    .Filter("Delta_OpeningAngle_Kp_pi1 > 0.0005")
    .Filter("Delta_OpeningAngle_Kp_pi2 > 0.0005")
    .Filter("Delta_OpeningAngle_pi1_pi2 > 0.0005")
    .Filter("cand_Bz_FD_rho < 4.8");

    auto nTotal = *df.Count();
    auto nSel   = *df_cuts1.Count();
    auto nSel2   = *df_cuts2.Count();
    std::cout << "Entries before cuts: " << nTotal << std::endl;
    std::cout << "Entries after cuts:  " << nSel   << std::endl;
    std::cout << "Entries after cuts 2:  " << nSel2   << std::endl;
    std::cout << "Efficiency: " << 100.0 * nSel2 / nTotal << " %" << std::endl;

    df_cuts2.Snapshot("DecayTree", "Data_1_new.root");
}

/*
--------------------------------------------
    VARIABLES
--------------------------------------------
*/

int nBins = 100;
int mDp_low = 1825;
int mDp_max = 1915;

/*
--------------------------------------------
    FILE NAMES
--------------------------------------------
*/

std::string data_1 = "/Volumes/Rome Drive/Latest_INFN/INFN/root_files/data/Data_ALL_sw.root";
std::string data_folder = "/Volumes/Rome Drive/Latest_INFN/INFN/root_files/data/";
std::string plots = "/Volumes/Rome Drive/Latest_INFN/INFN/plots/";
std::string WS = "/Volumes/Rome Drive/Latest_INFN/INFN/root_files/data/Data_ALL_WS_sw.root";
std::string MC = "/Volumes/Rome Drive/Latest_INFN/INFN/root_files/MC/MC_2_eff_with_mCorr.root";
std::string MC_templates_50 = "/Volumes/Rome Drive/Latest_INFN/INFN/root_files/MC/templates_MC/templates_mcorr_50.root";
std::string MC_templates_100 = "/Volumes/Rome Drive/Latest_INFN/INFN/root_files/MC/templates_MC/templates_mcorr_100.root";

/*
--------------------------------------------
    ADDITION p_PERP TO DECAYTREE
--------------------------------------------
*/

void TreeAdd_pPerp() {
    ROOT::RDataFrame df("DecayTree", data_1.c_str());

    //auto df_filtered = df.Filter("cand_Bz_OWNPV_CORRM > 4400 && (cand_Bz_Dp_M < 1845 || cand_Bz_Dp_M > 1895)");

    std::cout << "Entries before filter: " << *df.Count() << std::endl;
    //std::cout << "Entries after filter: " << *df_filtered.Count() << std::endl;

    auto df2 = df
        .Define("flight_dir", [](double svx, double svy, double svz,
                                 float pvx, float pvy, float pvz) {
            XYZVector dir(svx - pvx, svy - pvy, svz - pvz);
            return dir.Unit();
        }, {"cand_Bz_ENDVERTEX_POS_X", "cand_Bz_ENDVERTEX_POS_Y", "cand_Bz_ENDVERTEX_POS_Z",
            "cand_Bz_OWNPV_X", "cand_Bz_OWNPV_Y", "cand_Bz_OWNPV_Z"})

        .Define("p_perp_vis", [](float Dpx, float Dpy, float Dpz, double D_M,
                                float MU_px, float MU_py, float MU_pz, double MU_M, 
                                const XYZVector &flight_dir) {
                PxPyPzMVector pD_m(Dpx, Dpy, Dpz, D_M);
                PxPyPzEVector pD(pD_m);
                PxPyPzMVector pMu_m(MU_px, MU_py, MU_pz, MU_M);
                PxPyPzEVector pMu(pMu_m);
                PxPyPzEVector pVis = pD + pMu;
                double p_par = pVis.Vect().Dot(flight_dir);
                XYZVector p_perp_vec = pVis.Vect() - p_par * flight_dir;
                return p_perp_vec.R();    

        }, {"cand_Bz_Dp_PX", "cand_Bz_Dp_PY", "cand_Bz_Dp_PZ", "cand_Bz_Dp_M",
            "cand_Bz_mup_PX", "cand_Bz_mup_PY", "cand_Bz_mup_PZ", "cand_Bz_mup_M",
            "flight_dir"})
            
        .Define("p_perp_D", [](float Dpx, float Dpy, float Dpz, double D_M, const XYZVector &flight_dir) {
                PxPyPzMVector pD_m(Dpx, Dpy, Dpz, D_M);
                PxPyPzEVector pD(pD_m);
                double p_par = pD.Vect().Dot(flight_dir);
                XYZVector p_perp_D_vec = pD.Vect() - p_par * flight_dir;
                return p_perp_D_vec.R();
        }, {"cand_Bz_Dp_PX", "cand_Bz_Dp_PY", "cand_Bz_Dp_PZ", "cand_Bz_Dp_M", "flight_dir"})

        .Define("mCorr", [](double Bz_M, double p_perp) {
            return std::sqrt(Bz_M * Bz_M + p_perp * p_perp) + p_perp;
        }, {"cand_Bz_M", "p_perp_vis"});
    

    df2.Snapshot("DecayTree", data_folder + "Data_ALL_pPERP.root");

}

void plot_pPerp() {
    using namespace ROOT;

    RDataFrame df0("DecayTree", data_folder + "Data_ALL_pPERP.root");

    // Gev conversion
    auto df = df0
        .Define("p_perp_D_GeV", "p_perp_D / 1000.0")
        .Define("mCorr_GeV", "cand_Bz_OWNPV_CORRM / 1000.0");

    const int nx = 100;  const double xlo = 3500, xhi = 5500;   // m_corr [meV]
    const int ny = 100;  const double ylo = 0.0, yhi = 3000;   // p_perp [meV]

    auto h_all = df.Histo2D(
        {"h_all",
         "All candidates; m_{corr} [MeV/c^{2}]; p_{perp}(D^{+}) [MeV/c]",
         nx, xlo, xhi, ny, ylo, yhi},
         "cand_Bz_OWNPV_CORRM", "p_perp_D"
    );

    // Normalisation to 1
    double max = h_all->GetMaximum();
    if (max > 0) h_all->Scale(1.0 / max);
    h_all->SetMinimum(0.0);
    h_all->SetMaximum(1.0);

    gStyle->SetOptStat(0);

    TCanvas* c = new TCanvas("c_pperp_mcorr_data", "p_perp vs mCorr (data)", 800, 700);
    c->SetLeftMargin(0.13);
    c->SetRightMargin(0.15);
    c->SetBottomMargin(0.13);
    c->SetTopMargin(0.08);

    h_all->Draw("COLZ");

    c->SaveAs((plots+"pPerp_vs_mCorr.png").c_str());

    delete c;
}


/*
-------------------------------------------
    KINEMATICS - BACKGROUND CHECKS
-------------------------------------------
*/


void plot_Dp_M(){

    TFile *fInput1 = TFile::Open(data_1.c_str());
    TTree *tree1 = (TTree*)fInput1->Get("DecayTree");

    TFile *fInput2 = TFile::Open(WS.c_str());
    TTree *tree2 = (TTree*)fInput2->Get("DecayTree");

    TCanvas *c = new TCanvas("c","Dp_M",1600,800);
    c->Divide(2,1);

    c->cd(1);
    tree1->Draw("cand_Bz_Dp_M>>hRS(100,1820,1920)");
    TH1 *hRS = (TH1*)gDirectory->Get("hRS");
    hRS->GetXaxis()->SetTitle("Dp_M [MeV]");
    hRS->SetLineColor(kRed);
    hRS->Draw();

    c->cd(2);
    tree2->Draw("cand_Bz_Dp_M>>hWS(100,1820,1920)");
    TH1 *hWS = (TH1*)gDirectory->Get("hWS");
    hWS->GetXaxis()->SetTitle("Dp_M [MeV]");
    hWS->SetLineColor(kBlue);
    hWS->Draw();

    c->SaveAs("Dp_M_comparison.png");
}

void SB_RS(){

    TFile *fInput = TFile::Open(data_1.c_str());
    TTree *tree = (TTree*)fInput->Get("DecayTree");

    TCanvas *c = new TCanvas("c","SB removal",800,600);

    // Data
    tree->Draw("cand_Bz_OWNPV_CORRM>>hData(100,2000,8000)");

    // Sidebands
    tree->Draw("cand_Bz_OWNPV_CORRM>>hSB(100,2000,8000)",
               "cand_Bz_Dp_M < 1840 || cand_Bz_Dp_M > 1900",
               "SAME");

    TH1 *hData = (TH1*)gDirectory->Get("hData");
    TH1 *hSB  = (TH1*)gDirectory->Get("hSB");

    double nSB = hSB->Integral();
    double nSigWinExpected = nSB * ((1900-1840) / ((1840-1826)+(1916-1900))); // This window can be redefined if needed
    hSB->Scale(nSigWinExpected / nSB);

    hData->SetLineColor(kRed);
    hSB->SetLineColor(kBlue);

    TLegend *leg = new TLegend(0.65,0.75,0.88,0.88);
    leg->AddEntry(hData,"Data all","l");
    leg->AddEntry(hSB,"Scaled sidebands","l");
    leg->Draw();

    c->SaveAs("SB_RS.png");
}

void SB_WS(){

    TFile *fInput = TFile::Open(WS.c_str());
    TTree *tree = (TTree*)fInput->Get("DecayTree");

    TCanvas *c = new TCanvas("c","SB removal",800,600);

    // Data
    tree->Draw("cand_Bz_OWNPV_CORRM>>hData(100,2000,8000)");

    // Sidebands
    tree->Draw("cand_Bz_OWNPV_CORRM>>hSB(100,2000,8000)",
               "cand_Bz_Dp_M < 1840 || cand_Bz_Dp_M > 1900",
               "SAME");

    TH1 *hData = (TH1*)gDirectory->Get("hData");
    TH1 *hSB  = (TH1*)gDirectory->Get("hSB");

    double nSB = hSB->Integral();
    double nSigWinExpected = nSB * ((1895-1845) / ((1840-1826)+(1914-1900)));
    hSB->Scale(nSigWinExpected / nSB);

    hData->SetLineColor(kRed);
    hSB->SetLineColor(kBlue);

    TLegend *leg = new TLegend(0.65,0.75,0.88,0.88);
    leg->AddEntry(hData,"Data all","l");
    leg->AddEntry(hSB,"Scaled sidebands","l");
    leg->Draw();

    c->SaveAs("SB_WS.png");
}

void RS_WS_SB_comparison(){

    TFile *fRS = TFile::Open(data_1.c_str());
    TFile *fWS = TFile::Open(WS.c_str());

    TTree *tRS = (TTree*)fRS->Get("DecayTree");
    TTree *tWS = (TTree*)fWS->Get("DecayTree");

    TCanvas *c = new TCanvas("c","RS vs WS",1000,800);

    tRS->Draw("cand_Bz_OWNPV_CORRM>>hRS(100,2000,8000)", "cand_Bz_Dp_M < 1840 || cand_Bz_Dp_M > 1900");
    tWS->Draw("cand_Bz_OWNPV_CORRM>>hWS(100,2000,8000)","cand_Bz_Dp_M < 1840 || cand_Bz_Dp_M > 1900", "SAME");

    TH1 *hRS = (TH1*)gDirectory->Get("hRS");
    TH1 *hWS = (TH1*)gDirectory->Get("hWS");

    hRS->Scale(1.0 / hRS->GetMaximum());
    hWS->Scale(1.0 / hWS->GetMaximum());

    hRS->SetLineColor(kRed);
    hWS->SetLineColor(kBlue);

    hRS->SetTitle("RS vs WS;Corrected mass;Normalized to 1");

    hRS->Draw("HIST");
    hWS->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.65,0.75,0.88,0.88);
    leg->AddEntry(hRS,"RS","l");
    leg->AddEntry(hWS,"WS","l");
    leg->Draw();

    c->SaveAs("RS_WS_SB_comparison.png");
}

void correlation_check() {
    
    TFile f(data_1.c_str(),"READ");
    TTree* t = (TTree*)f.Get("DecayTree");

    TCanvas *c1 = new TCanvas("c1","Correlation",800,600);

    TH2F *h2 = new TH2F(
        "h2",
        "Correlation;Dp mass;Corrected mass",
        50,1820,1920,
        50,2000,8000
    );

    t->Draw("cand_Bz_OWNPV_CORRM:cand_Bz_Dp_M >> h2", "", "COLZ");

    c1->Update();

    std::cout << "Entries: " << h2->GetEntries() << std::endl;
    std::cout << "Correlation factor: "
              << h2->GetCorrelationFactor() << std::endl;

    c1->SaveAs("Correlation_RS.png");
}

void plot_logIP() {
    using namespace ROOT;

    RDataFrame df0("DecayTree", data_1.c_str());

    auto df = df0
        .Define("logIP", "log(cand_Bz_Dp_OWNPV_IP)");

    const int nx = 100;  const double xlo = 2000, xhi = 8000;   // m_corr [GeV]
    const int ny = 100;  const double ylo = -7, yhi = 5;   // p_perp [GeV]

    auto h_all = df.Histo2D(
        {"h_all",
         "LogIP-OWNPV_CORRM: WS after cuts; m_{corr} [MeV]; logIP",
         nx, xlo, xhi, ny, ylo, yhi},
         "cand_Bz_OWNPV_CORRM", "logIP"
    );

    // Normalisation
    double max = h_all->GetMaximum();
    if (max > 0) h_all->Scale(1.0 / max);
    h_all->SetMinimum(0.0);
    h_all->SetMaximum(1.0);

    gStyle->SetOptStat(0);

    TCanvas* c = new TCanvas("c_pperp_mcorr_data", "p_perp vs mCorr (data)", 800, 700);
    c->SetLeftMargin(0.13);
    c->SetRightMargin(0.15);
    c->SetBottomMargin(0.13);
    c->SetTopMargin(0.08);

    h_all->Draw("COLZ");

    c->SaveAs("mcorr-logip_RS_before.png");

    delete c;
}


/*
-------------------------------------------
    Sidebands - Binned Fit - MC treated as PDF
-------------------------------------------
*/


void fit_SB_Binned()
{
    const double mCorr_min = 2000;
    const double mCorr_max = 8000;
    const int nBins = 100;

    const double mDp_sig_lo = 1840;
    const double mDp_sig_hi = 1900;

    const double mDp_ws_sb_lo = 1840;
    const double mDp_ws_sb_hi = 1900;

    TFile fmc(MC_templates_50.c_str(),"READ"); // File acquired from main_MC.cpp. It contains the histograms of the signal components.

    TH1D* h_Dp = (TH1D*)fmc.Get("h_Dp_w")->Clone("h_Dp");
    TH1D* h_Dst = (TH1D*)fmc.Get("h_Dst_w")->Clone("h_Dst");
    TH1D* h_Db = (TH1D*)fmc.Get("h_Dstst_broad_w")->Clone("h_Db"); 
    TH1D* h_Dn = (TH1D*)fmc.Get("h_Dstst_narrow_w")->Clone("h_Dn");

    h_Dp->SetDirectory(nullptr);
    h_Dst->SetDirectory(nullptr);
    h_Db->SetDirectory(nullptr);
    h_Dn->SetDirectory(nullptr);

    fmc.Close();

    ROOT::RDataFrame df("DecayTree", data_1.c_str());

    ROOT::RDF::TH1DModel modelData("hData","",nBins,mCorr_min,mCorr_max);
    TH1D* hData = (TH1D*)df.Histo1D(modelData,"cand_Bz_OWNPV_CORRM").GetPtr()->Clone("hData");
    hData->SetDirectory(nullptr);

    ROOT::RDF::TH1DModel modelSB("hSB_raw","",nBins,mCorr_min,mCorr_max);
    TH1D* hSB_raw = (TH1D*)df
        .Filter(Form("cand_Bz_Dp_M < %f || cand_Bz_Dp_M > %f", mDp_sig_lo, mDp_sig_hi))
        .Histo1D(modelSB,"cand_Bz_OWNPV_CORRM")
        .GetPtr()->Clone("hSB_raw");

    hSB_raw->SetDirectory(nullptr);

    double nSB = hSB_raw->Integral();
    double scale_factor = (1900.-1840) / ((1840-1826)+(1916-1900));
    double N_SB_expected = hSB_raw->Integral() * scale_factor;

    hSB_raw->Scale(1.0 / hSB_raw->Integral());

    ROOT::RDataFrame dfWS("DecayTree", WS.c_str());

    ROOT::RDF::TH1DModel modelWS("hWS_SB_raw","",nBins,mCorr_min,mCorr_max);
    TH1D* hWS_SB = (TH1D*)dfWS
        .Filter(Form("cand_Bz_Dp_M > %f && cand_Bz_Dp_M < %f", mDp_ws_sb_lo, mDp_ws_sb_hi))
        .Histo1D(modelWS,"cand_Bz_OWNPV_CORRM")
        .GetPtr()->Clone("hWS_SB");

    hWS_SB->SetDirectory(nullptr);

    hWS_SB->Scale(1.0 / hWS_SB->Integral());

    RooRealVar x("x","cand_Bz_OWNPV_CORRM [MeV]",mCorr_min,mCorr_max);
    x.setBins(nBins);

    RooDataHist dh_Dp ("dh_Dp","",RooArgList(x),h_Dp);
    RooDataHist dh_Dst("dh_Dst","",RooArgList(x),h_Dst);
    RooDataHist dh_Db ("dh_Db","",RooArgList(x),h_Db);
    RooDataHist dh_Dn ("dh_Dn","",RooArgList(x),h_Dn);
    RooDataHist dh_SB ("dh_SB","",RooArgList(x),hSB_raw);
    RooDataHist dh_WS ("dh_WS","",RooArgList(x),hWS_SB);

    RooHistPdf pdf_Dp ("pdf_Dp","",x,dh_Dp);
    RooHistPdf pdf_Dst("pdf_Dst","",x,dh_Dst);
    RooHistPdf pdf_Db ("pdf_Db","",x,dh_Db);
    RooHistPdf pdf_Dn ("pdf_Dn","",x,dh_Dn);
    RooHistPdf pdf_SB ("pdf_SB","",x,dh_SB);
    RooHistPdf pdf_WS ("pdf_WS","",x,dh_WS);

    RooRealVar N_Dp ("N_Dp","",0,0,1e7);
    RooRealVar N_Dst("N_Dst","",0,0,1e7);

    RooRealVar N_Db ("N_Db","",0,0,1e7);
    RooRealVar N_Dn ("N_Dn","",0,0,1e7);

    RooRealVar N_SB("N_SB","",N_SB_expected,0,1e7);
    N_SB.setConstant(true);

    RooRealVar N_WS("N_WS","", 0, 0, 1e7);

    RooAddPdf model(
        "model","",
        RooArgList(pdf_Dp,pdf_Dst,pdf_Db,pdf_Dn,pdf_SB,pdf_WS),
        RooArgList(N_Dp,N_Dst,N_Db,N_Dn,N_SB,N_WS) 
    );

    RooDataHist dh_data("dh_data","",RooArgList(x),hData);
    std::unique_ptr<RooFitResult> fr(    
        model.fitTo(dh_data,Extended(true), Save(true))
    );

    TCanvas c("c","",900,700);
    RooPlot* frame = x.frame();

    dh_data.plotOn(frame, DataError(RooAbsData::SumW2));

    model.plotOn(frame, LineColor(kBlack), LineWidth(3));

    model.plotOn(frame, Components("pdf_Dp"), LineColor(kRed), LineWidth(2));
    model.plotOn(frame, Components("pdf_Dst"), LineColor(kBlue), LineWidth(2));
    model.plotOn(frame, Components("pdf_Db"), LineColor(kMagenta), LineWidth(2));
    model.plotOn(frame, Components("pdf_Dn"), LineColor(kOrange+7), LineWidth(2));
    model.plotOn(frame, Components("pdf_SB"), LineColor(kYellow+2), LineWidth(2));
    model.plotOn(frame, Components("pdf_WS"), LineColor(kGreen+2), LineWidth(2));


    frame->SetTitle("Binned fit of cand_Bz_OWNPV_CORRM with sideband background");
    frame->Draw();
    gPad->cd();

    TLegend leg(0.60,0.60,0.88,0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);

    leg.AddEntry(frame->getObject(0), "Data", "pe");
    leg.AddEntry(frame->getObject(1), "Total fit", "l");
    leg.AddEntry(frame->getObject(2), "D+", "l");
    leg.AddEntry(frame->getObject(3), "D*+", "l");
    leg.AddEntry(frame->getObject(4), "D** broad", "l");
    leg.AddEntry(frame->getObject(5), "D** narrow", "l");
    leg.AddEntry(frame->getObject(6), "RS sideband bkg", "l");
    leg.AddEntry(frame->getObject(7), "WS sideband bkg", "l");


    leg.Draw();

    c.SaveAs("Fit_SB_Binned.png");
    fr->Print("v");
    fr->correlationMatrix().Print();

}


/*
-------------------------------------------
    Sidebands - Unbinned Fit: RS BKG (from sidebands) + WS BKG (from signal region)
-------------------------------------------
*/


void fit_SB_Unbinned()
{
    TFile *fData = TFile::Open(data_1.c_str());
    TTree *tData = (TTree*)fData->Get("DecayTree");

    TFile *fWS = TFile::Open(WS.c_str());
    TTree *tWS = (TTree*)fWS->Get("DecayTree");

    RooRealVar x("cand_Bz_OWNPV_CORRM","Corrected mass [MeV]",2250.,8000.);
    RooRealVar m("cand_Bz_Dp_M","m(D^{+}) [MeV]",1830.,1910.);
    x.setBins(100);

    const double sigLo = 1840.;
    const double sigHi = 1900.;
    const double SBscale_RS = (1900-1840) / ((1840-1826)+(1916-1900));

    RooDataSet data_sig(
        "data_sig","",
        tData,
        RooArgSet(x,m),
        Form("cand_Bz_Dp_M>%f && cand_Bz_Dp_M<%f",sigLo,sigHi)
    );

    RooDataSet data_sb(
        "data_sb","",
        tData,
        RooArgSet(x,m),
        Form("cand_Bz_Dp_M<%f || cand_Bz_Dp_M>%f",sigLo,sigHi)
    );

    RooDataSet data_ws(
        "data_ws_sb","",
        tWS,
        RooArgSet(x,m),
        Form("cand_Bz_Dp_M>%f && cand_Bz_Dp_M<%f",sigLo,sigHi)
    );

    const double Nsig_obs = data_sig.numEntries();
    const double NSB_obs  = data_sb.numEntries();

    RooDataHist sbHist("sbHist","",RooArgSet(x),data_sb);
    RooHistPdf pdf_bkg("pdf_bkg","",x,sbHist);

    RooDataHist wsHist("wsHist","",RooArgSet(x),data_ws);
    RooHistPdf pdf_ws("pdf_ws","",x,wsHist);

    auto makeDCB = [&](const char* name,
                       double mean,double sigma,
                       double aL,double nL,
                       double aR,double nR)
    {
        RooConstVar* m  = new RooConstVar(Form("mean_%s",name),"",mean);
        RooConstVar* s  = new RooConstVar(Form("sigma_%s",name),"",sigma);
        RooConstVar* al = new RooConstVar(Form("aL_%s",name),"",aL);
        RooConstVar* nl = new RooConstVar(Form("nL_%s",name),"",nL);
        RooConstVar* ar = new RooConstVar(Form("aR_%s",name),"",aR);
        RooConstVar* nr = new RooConstVar(Form("nR_%s",name),"",nR);

        return new RooCrystalBall(
            Form("pdf_%s",name),"",
            x,*m,*s,*al,*nl,*ar,*nr
        );
    };

    auto makeGauss = [&](const char* name,double mean,double sigma)
    {
        RooConstVar* m = new RooConstVar(Form("mean_%s",name),"",mean);
        RooConstVar* s = new RooConstVar(Form("sigma_%s",name),"",sigma);

        return new RooGaussian(Form("pdf_%s",name),"",x,*m,*s);
    };

    // Values manually copied from the MC fit => no propagation of error
    auto pdf_Dp = makeDCB("Dp", 5161.4,209.405,0.456019,115,1.10143,4.05265);
    auto pdf_Dpst  = makeDCB("Dpst", 4944.26,250.588,0.618174,114.995,1.42973,4.07023);
    auto pdf_Dst_narrow = makeGauss("Dst_narrow",4446.08,574.064);
    auto pdf_Dst_broad = makeGauss("Dst_broad",4407.21,517.073);

    RooConstVar n_bkg("n_bkg","",SBscale_RS*NSB_obs);

    RooRealVar n_ws("n_ws","",0.1*Nsig_obs,0.,0.5*Nsig_obs);
    RooRealVar n_Dp("n_Dp","",0.4*Nsig_obs,0.,Nsig_obs);
    RooRealVar n_Dpst("n_Dpst","",0.3*Nsig_obs,0.,Nsig_obs);
    RooRealVar n_Dst_narrow("n_Dst_narrow","",0.05*Nsig_obs,0.,0.3*Nsig_obs);
    RooRealVar n_Dst_broad("n_Dst_broad","",0.02*Nsig_obs,0.,0.2*Nsig_obs);

    RooAddPdf model(
        "model","",
        RooArgList(*pdf_Dp,
                   *pdf_Dpst,
                   *pdf_Dst_narrow,
                   *pdf_Dst_broad,
                   pdf_bkg,
                   pdf_ws),
        RooArgList(n_Dp,
                   n_Dpst,
                   n_Dst_narrow,
                   n_Dst_broad,
                   n_bkg,
                   n_ws)
    );

    auto fitres = model.fitTo(data_sig,Save(true),Extended(true));
    fitres->Print("v");

    TCanvas c("c","",900,700);
    RooPlot* frame = x.frame();
    data_sig.plotOn(frame);
    model.plotOn(frame,LineColor(kBlack),LineWidth(3));
    model.plotOn(frame,Components(pdf_bkg),LineColor(kGreen+2));
    model.plotOn(frame,Components(pdf_ws),LineColor(kCyan+2),LineStyle(kDashed));
    model.plotOn(frame,Components("pdf_Dp"),LineColor(kRed));
    model.plotOn(frame,Components("pdf_Dpst"),LineColor(kBlue));
    model.plotOn(frame,Components("pdf_Dst_narrow"),LineColor(kMagenta));
    model.plotOn(frame,Components("pdf_Dst_broad"),LineColor(kOrange+7));
    frame->Draw();

    TLegend leg(0.60,0.60,0.88,0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(frame->getObject(0),"Data","pe");
    leg.AddEntry(frame->getObject(1),"Total fit","l");
    leg.AddEntry(frame->getObject(2),"RS background","l");
    leg.AddEntry(frame->getObject(3),"WS background","l");
    leg.AddEntry(frame->getObject(4),"D^{+}","l");
    leg.AddEntry(frame->getObject(5),"D*","l");
    leg.AddEntry(frame->getObject(6),"D** broad","l");
    leg.AddEntry(frame->getObject(7),"D** narrow","l");
    leg.Draw();

    c.SaveAs("Unbinned_fit_SB_WS.png");

    std::cout << "Observed SB events = " << NSB_obs << std::endl;
    std::cout << "Predicted RS bkg   = " << n_bkg.getVal() << std::endl;
    std::cout << "Observed signal   = " << Nsig_obs << std::endl;
    std::cout << "Fitted WS yield   = " << n_ws.getVal()
              << " +- " << n_ws.getError() << std::endl;
    fitres->correlationMatrix().Print();
}


/*
----------------------------------------
    BEESTOW-BARLOW implementation --- NOT WORKING due to MC statistics
----------------------------------------
*/

void fit_BB()
{
    const double mCorr_min = 2000;
    const double mCorr_max = 8000;
    const int nBins = 50;

    const double mDp_sig_lo = 1840;
    const double mDp_sig_hi = 1900;

    TFile *fmc = TFile::Open(MC_templates_50.c_str(), "READ");

    TH1D* h_Dp  = (TH1D*)fmc->Get("h_Dp_w")->Clone("h_Dp_c");
    TH1D* h_Dst = (TH1D*)fmc->Get("h_Dst_w")->Clone("h_Dst_c");
    TH1D* h_Db  = (TH1D*)fmc->Get("h_Dstst_broad_w")->Clone("h_Db_c");
    TH1D* h_Dn  = (TH1D*)fmc->Get("h_Dstst_narrow_w")->Clone("h_Dn_c");

    h_Dp->SetDirectory(nullptr);
    h_Dst->SetDirectory(nullptr);
    h_Db->SetDirectory(nullptr);
    h_Dn->SetDirectory(nullptr);

    fmc->Close();

    ROOT::RDataFrame dfData("DecayTree",data_1.c_str());

    ROOT::RDF::TH1DModel modelData("hData",";m_{corr};Entries",nBins,mCorr_min,mCorr_max);
    TH1D* hData = (TH1D*)dfData
        .Histo1D(modelData,"cand_Bz_OWNPV_CORRM")
        .GetPtr()->Clone("hData_c");
    hData->SetDirectory(nullptr);

    RooRealVar x("cand_Bz_OWNPV_CORRM", "m_{corr}", mCorr_min, mCorr_max);
    x.setBins(nBins);

    RooDataHist data("data","binned data",x,hData);

    ROOT::RDF::TH1DModel modelSB("hSB","hSB",nBins,mCorr_min,mCorr_max);
    TH1D* hSB = (TH1D*)dfData
        .Filter(Form("cand_Bz_Dp_M<%f || cand_Bz_Dp_M>%f",mDp_sig_lo, mDp_sig_hi))
        .Histo1D(modelSB,"cand_Bz_OWNPV_CORRM")
        .GetPtr()->Clone("hSB_c");
    hSB->SetDirectory(nullptr);

    ROOT::RDF::TH1DModel modelSignal("hSig","hSig",nBins,mCorr_min,mCorr_max);
    TH1D* hSig = (TH1D*)dfData
        .Filter(Form("cand_Bz_Dp_M>%f && cand_Bz_Dp_M<%f",mDp_sig_lo, mDp_sig_hi))
        .Histo1D(modelSB,"cand_Bz_OWNPV_CORRM")
        .GetPtr()->Clone("hSig_c");
    hSig->SetDirectory(nullptr);

    double nSB = hSB->Integral();
    double scale_factor = (1900.-1840) / ((1840-1821)+(1916-1900)); // 1.71
    double N_SB_expected = nSB * scale_factor;  

    ROOT::RDataFrame dfWS("DecayTree",WS.c_str());
    ROOT::RDF::TH1DModel modelWS("hWS","hWS",nBins,mCorr_min,mCorr_max);
    TH1D* hWS = (TH1D*)dfWS
        .Filter(Form("cand_Bz_Dp_M>%f && cand_Bz_Dp_M<%f",mDp_sig_lo,mDp_sig_hi))
        .Histo1D(modelWS,"cand_Bz_OWNPV_CORRM")
        .GetPtr()->Clone("hWS");
    hWS->SetDirectory(nullptr);

    RooDataHist dh_Dp("dh_Dp","dh_Dp",x,h_Dp);
    RooDataHist dh_Dst("dh_Dst","dh_Dst",x,h_Dst);
    RooDataHist dh_Db("dh_Db","dh_Db",x,h_Db);
    RooDataHist dh_Dn("dh_Dn","dh_Dn",x,h_Dn);
    RooDataHist dh_SB("dh_SB","dh_SB",x,hSB);
    RooDataHist dh_WS("dh_WS","dh_WS",x,hWS);

    // ParamHistFunc -> for MC only as it takes into account MC fluctuations
    RooParamHistFunc p_Dp ("p_Dp","",dh_Dp,x);
    RooParamHistFunc p_Dst("p_Dst","",dh_Dst,x, &p_Dp, true);
    RooParamHistFunc p_Db ("p_Db","",dh_Db,x);
    RooParamHistFunc p_Dn ("p_Dn","",dh_Dn ,x, &p_Db, true);
    RooHistFunc hf_SB ("hf_SB","",x,dh_SB);
    RooHistFunc hf_WS ("hf_WS","",x,dh_WS);

    RooRealVar n_Dp   ("n_Dp"  ,"n_Dp"  ,10000,0.,100000); 
    RooRealVar n_Dpst ("n_Dpst","n_Dpst",240,0.,1000000);
    RooRealVar n_Db   ("n_Db"  ,"n_Db"  ,200,0.,100000);
    RooRealVar n_Dn   ("n_Dn"  ,"n_Dn"  ,182,0.,100000);

    RooRealVar n_SB("n_SB", "",0.014285714,0,10.); //3.3*scale_factor -> 3.27 stabilises the fit
    RooRealVar n_WS("n_WS","n_WS",1.5,0,100.);

    RooRealSumPdf totalModel(
        "totalModel","",
        RooArgList(p_Dp,p_Dst,p_Db,p_Dn,hf_SB,hf_WS),
        RooArgList(n_Dp,n_Dpst,n_Db,n_Dn,n_SB,n_WS),
        true
    );

    RooHistConstraint hc_Dp("hc_Dp","hc_Dp",p_Dp);
    RooHistConstraint hc_Dst("hc_Dst","hc_Dst",p_Dst);
    RooHistConstraint hc_Db("hc_Db","hc_Db",p_Db);
    RooHistConstraint hc_Dn("hc_Dn","hc_Dn",p_Dn);

    RooHistConstraint hc_1("hc_1","hc_1",RooArgSet(p_Dp,p_Dst));
    RooHistConstraint hc_2("hc_2","hc_2",RooArgSet(p_Db,p_Dn));

    RooProdPdf model(
        "model","", 
        RooArgSet(hc_1,hc_2),
        Conditional(totalModel, x)
    );

    std::unique_ptr<RooFitResult> fr(
        model.fitTo(
            data,
            Save(true),
            Strategy(2),
            MaxCalls(2147483647),
            Hesse(true)
        )
    );

    TCanvas c("c","Fit",900,700);
    RooPlot* frame = x.frame();
    frame->SetTitle("Barlow-Beeston light");

    data.plotOn(frame,Name("data"));
    model.plotOn(frame,LineColor(kRed),LineWidth(3),Name("total"));
    model.plotOn(frame,Components(p_Dp),LineColor(kBlue),Name("Dp"));
    model.plotOn(frame,Components(p_Dst),LineColor(kGreen+2),Name("Dst"));
    model.plotOn(frame,Components(p_Db),LineColor(kMagenta),Name("Db"));
    model.plotOn(frame,Components(p_Dn),LineColor(kOrange+7),Name("Dn"));
    model.plotOn(frame,Components(hf_SB),LineColor(kYellow+3),Name("SB")); 
    model.plotOn(frame,Components(hf_WS),LineColor(kRed+2),Name("WS")); 

    frame->Draw();

    TLegend leg(0.60,0.50,0.88,0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(frame->findObject("data"),"Data","lep");
    leg.AddEntry(frame->findObject("total"),"Total fit","l");
    leg.AddEntry(frame->findObject("Dp"),"MC Dp","l");
    leg.AddEntry(frame->findObject("Dst"),"MC D*","l");
    leg.AddEntry(frame->findObject("Db"),"MC D** broad","l");
    leg.AddEntry(frame->findObject("Dn"),"MC D** narrow","l");
    leg.AddEntry(frame->findObject("SB"),"RS sidebands","l");
    leg.AddEntry(frame->findObject("WS"),"WS","l");
    leg.Draw();

    c.SaveAs("Fit_BB_DATASET.png");
    fr->Print("v");
    

}

    /* Stacked plot form for fit_BB()

    TH1* hDp_fit  = p_Dp .createHistogram("hDp_fit",  x, Binning(nBins));
    TH1* hDst_fit = p_Dst.createHistogram("hDst_fit", x, Binning(nBins));
    TH1* hDb_fit  = p_Db .createHistogram("hDb_fit",  x, Binning(nBins));
    TH1* hDn_fit  = p_Dn .createHistogram("hDn_fit",  x, Binning(nBins));
    TH1* hSB_fit  = hf_SB.createHistogram("hSB_fit",  x, Binning(nBins));
    TH1* hWS_fit  = hf_WS.createHistogram("hWS_fit",  x, Binning(nBins));

    hDp_fit->Scale(n_Dp.getVal());
    hDst_fit->Scale(n_Dpst.getVal());
    hDb_fit->Scale(n_Db.getVal());
    hDn_fit->Scale(n_Dn.getVal());
    hSB_fit->Scale(n_SB.getVal());
    hWS_fit->Scale(n_WS.getVal());

    hDp_fit ->SetFillColor(kBlue-7);
    hDst_fit->SetFillColor(kGreen-6);
    hDb_fit ->SetFillColor(kMagenta-6);
    hDn_fit ->SetFillColor(kOrange-3);
    hSB_fit ->SetFillColor(kYellow-7);
    hWS_fit ->SetFillColor(kRed-7);

    THStack hs("hs",";m_{corr} [MeV];Events");
    hs.Add(hWS_fit);
    hs.Add(hSB_fit);
    hs.Add(hDn_fit);
    hs.Add(hDb_fit);
    hs.Add(hDst_fit);
    hs.Add(hDp_fit);

    hs.Draw("HIST");

    hData->SetMarkerStyle(kFullCircle);
    hData->Draw("E1 SAME");

    TLegend leg(0.60,0.50,0.88,0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);

    leg.AddEntry(hData   ,"Data","lep");
    leg.AddEntry(hDp_fit ,"MC D^{+}","f");
    leg.AddEntry(hDst_fit,"MC D^{*+}","f");
    leg.AddEntry(hDb_fit ,"MC D^{**} broad","f");
    leg.AddEntry(hDn_fit ,"MC D^{**} narrow","f");
    leg.AddEntry(hSB_fit ,"RS sidebands","f");
    leg.AddEntry(hWS_fit ,"WS","f");

    leg.Draw();
*/


/*
----------------------------------------
    sWeights - Binned 
----------------------------------------

    NOTE!!!
    The range of Dp_M is set to 1830-1910MeV.

    The sWeights are extracted from the MC shape of the Dp_M and then the same parameters
    are used for the fitting on data.

*/


void sweights_extraction_MC()
{
    ROOT::RDataFrame df_MC("DecayTree", MC.c_str());
    ROOT::RDF::TH1DModel h_mD("h_mD", ";m_{D+} [MeV];Entries", nBins, mDp_low, mDp_max);
    auto hMC_temp = df_MC.Histo1D(h_mD, "Bz_Dp_M");

    auto& hMC_ref = hMC_temp.GetValue();
    TH1D* hMC = (TH1D*)hMC_ref.Clone("h_mD");
    hMC->SetDirectory(nullptr);

    RooRealVar x("cand_Bz_Dp_M", "m_{D+}", 1820, 1920);
    x.setRange("fitWindow",1830,1910);
    x.setBins(nBins);

    RooDataHist dh_MC("dh_MC", "dh_MC", RooArgList(x), hMC);

    RooRealVar mean("mean","mean",1869,1860,1880);
    RooRealVar sigma("sigma","sigma",7,1,20);

    RooRealVar alphaL("alphaL","alphaL",1.5,0.1,5);
    RooRealVar nL("nL","nL",5,0.1,100);

    RooRealVar alphaR("alphaR","alphaR",2.0,0.1,5);
    RooRealVar nR("nR","nR",10,0.1,200);

    RooCrystalBall crystal_ball("crystal_ball","crystal_ball",
                       x,mean,sigma,
                       alphaL,nL,
                       alphaR,nR);

    auto fitres = crystal_ball.fitTo(dh_MC, Save(true), Range("fitWindow"), Strategy(2));

    alphaL.setConstant(true);
    nL.setConstant(true);
    alphaR.setConstant(true);
    nR.setConstant(true);

    TCanvas c("c", "Fit", 900, 700);
    RooPlot* frame = x.frame();

    dh_MC.plotOn(frame, DataError(RooAbsData::SumW2));
    crystal_ball.plotOn(frame, Range("fitWindow"), NormRange("fitWindow"), LineColor(kRed), LineWidth(3)); // total fit
    crystal_ball.plotOn(frame, Range("fitWindow"), NormRange("fitWindow"), LineColor(kGreen+1), LineStyle(kDotted), Components(crystal_ball)); // signal

    frame->SetTitle("Fit of m_{D+} of MC with a Crystal Ball");
    frame->SetXTitle("m_{D+} [MeV]");
    frame->SetYTitle("Entries");
    frame->Draw();

    TLegend legend1(0.65, 0.65, 0.88, 0.88);
    legend1.SetBorderSize(0);
    legend1.SetFillStyle(0);
    legend1.AddEntry(frame->getObject(0), "Data", "pe"); 
    legend1.AddEntry(frame->getObject(1), "Total fit", "l");
    legend1.AddEntry(frame->getObject(2), "CB", "l");
    legend1.Draw();

    c.SaveAs((plots + "Dp_M_fit_MC.png").c_str());
    c.Close();

    // Fit on Data
    ROOT::RDataFrame df_data("DecayTree",data_1);
    ROOT::RDF::TH1DModel model("hData", ";m_{D+} [MeV];Entries", nBins, mDp_low, mDp_max);
    auto hData_temp = df_data.Histo1D(model, "cand_Bz_Dp_M");

    auto& hData_ref = hData_temp.GetValue();
    TH1D* hData = (TH1D*)hData_ref.Clone("hData");
    hData->SetDirectory(nullptr);

    RooRealVar tau("tau", "background slope", -0.002, -0.02, 0.0);
    RooExponential exp_bkg("exp_bkg", "Background", x, tau);
    
    RooRealVar N_sig("N_sig", "MC signal", hData->Integral()/2, 0, 1e9);
    RooRealVar N_bkg("N_bkg", "Background yield", hData->Integral()/5, 0, 1e9);

    RooAddPdf model0("model0", "MC template + exp bkg",
                    RooArgList(crystal_ball, exp_bkg),
                    RooArgList(N_sig, N_bkg));

    RooDataHist dh_data("dh_data", "Data", RooArgList(x), hData);

    auto fitres_data = model0.fitTo(dh_data, Save(true), Extended(true), Range("fitWindow"), Strategy(2));

    TCanvas c2("c2", "Fit", 900, 700);
    RooPlot* frame2 = x.frame();

    dh_data.plotOn(frame2, DataError(RooAbsData::SumW2));
    model0.plotOn(frame2, Range("fitWindow"), NormRange("fitWindow"), LineColor(kRed), LineWidth(3)); // total fit
    model0.plotOn(frame2, Range("fitWindow"), NormRange("fitWindow"), LineColor(kBlue), LineStyle(kDotted), Components(exp_bkg)); // background
    model0.plotOn(frame2, Range("fitWindow"), NormRange("fitWindow"), LineColor(kGreen+1), LineStyle(kDotted), Components(crystal_ball)); // signal

    frame2->SetTitle("Fit of m_{D+} using MC template DCB + exponential");
    frame2->SetXTitle("m_{D+} [MeV]");
    frame2->SetYTitle("Entries");
    frame2->Draw();

    TLegend legend2(0.65, 0.65, 0.88, 0.88);
    legend2.SetBorderSize(0);
    legend2.SetFillStyle(0);
    legend2.AddEntry(frame2->getObject(0), "Data", "pe"); 
    legend2.AddEntry(frame2->getObject(1), "Total fit", "l");
    legend2.AddEntry(frame2->getObject(2), "Background (exponential)", "l");
    legend2.AddEntry(frame2->getObject(3), "MC Signal (DCB)", "l");


    legend2.Draw();

    c2.SaveAs((plots + "Dp_M_fit_data.png").c_str());
    c2.Close();

    TFile *fInput = TFile::Open(data_1.c_str());
    TTree *tree = (TTree*)fInput->Get("DecayTree");

    RooRealVar Dp_M("cand_Bz_Dp_M", "m_{D+}", 1820, 1920);
    RooDataSet data("data", "Unbinned Data", RooArgSet(Dp_M), Import(*tree));

    RooStats::SPlot* sData = new RooStats::SPlot("sData", "SPlot", data, &model0, RooArgList(N_sig, N_bkg));
    RooDataSet* ds = sData->GetSDataSet();

    TTree* t_orig = (TTree*)fInput->Get("DecayTree");

    TFile fout(data_1.c_str(),"RECREATE");
    TTree* t_out = t_orig->CloneTree(0);
    t_out->SetName("DecayTree");

    double N_sig_sw_val, N_bkg_sw_val;
    t_out->Branch("N_sig_sw",&N_sig_sw_val,"N_sig_sw/D");
    t_out->Branch("N_bkg_sw",&N_bkg_sw_val,"N_bkg_sw/D");

    for (Long64_t i=0;i<t_orig->GetEntries();i++){
        t_orig->GetEntry(i);
        N_sig_sw_val = ds->get(i)->getRealValue("N_sig_sw");
        N_bkg_sw_val = ds->get(i)->getRealValue("N_bkg_sw");
        t_out->Fill();
    }

    t_out->Write("", TObject::kOverwrite);
    fout.Close();
}

void fit_sw_Binned()
{
    const double mCorr_min = 2200;
    const double mCorr_max = 7200;
    const int nBins = 50;

    ROOT::RDataFrame df_mc("DecayTree",MC.c_str());

    ROOT::RDF::TH1DModel modelDp("h_Dp_w","",nBins,mCorr_min,mCorr_max);
    TH1D* h_Dp = (TH1D*)df_mc
        .Filter("SL_category == 1")
        .Histo1D(modelDp,"Bz_OWNPV_CORRM","eff_Bz")
        .GetPtr()->Clone("h_Dp");
    h_Dp->SetDirectory(nullptr);

    ROOT::RDF::TH1DModel modelDst("h_Dst_w","",nBins,mCorr_min,mCorr_max);
    TH1D* h_Dst = (TH1D*)df_mc
        .Filter("SL_category == 2")
        .Histo1D(modelDst,"Bz_OWNPV_CORRM","eff_Bz")
        .GetPtr()->Clone("h_Dst");
    h_Dst->SetDirectory(nullptr);

    ROOT::RDF::TH1DModel modelDb("h_Dstst_broad_w","",nBins,mCorr_min,mCorr_max);
    TH1D* h_Db = (TH1D*)df_mc
        .Filter("SL_category == 31")
        .Histo1D(modelDb,"Bz_OWNPV_CORRM","eff_Bz")
        .GetPtr()->Clone("h_Db");
    h_Db->SetDirectory(nullptr);

    ROOT::RDF::TH1DModel modelDn("h_Dstst_narrow_w","",nBins,mCorr_min,mCorr_max);
    TH1D* h_Dn = (TH1D*)df_mc
        .Filter("SL_category == 32")
        .Histo1D(modelDn,"Bz_OWNPV_CORRM","eff_Bz")
        .GetPtr()->Clone("h_Dn");
    h_Dn->SetDirectory(nullptr);

    ROOT::RDataFrame df("DecayTree", data_1.c_str());

    ROOT::RDF::TH1DModel modelData("hData","",nBins,mCorr_min,mCorr_max);
    TH1D* hData = (TH1D*)df
        .Histo1D(modelData,"cand_Bz_OWNPV_CORRM", "N_sig_sw")
        .GetPtr()->Clone("hData");
    hData->SetDirectory(nullptr);

    ROOT::RDataFrame dfWS("DecayTree", WS.c_str());

    ROOT::RDF::TH1DModel modelWS("hWS_SB_raw","",nBins,mCorr_min,mCorr_max);
    TH1D* hWS = (TH1D*)dfWS
        .Histo1D(modelWS,"cand_Bz_OWNPV_CORRM", "N_sig_sw")
        .GetPtr()->Clone("hWS_SB");
    hWS->SetDirectory(nullptr);

    for (int i = 1; i <= hWS->GetNbinsX(); ++i) {
        double content = hWS->GetBinContent(i);
        if (content <= 0) {
            std::cout << "Bin " << i
                    << " | content = " << content
                    << std::endl;
        }
    }

    RooRealVar x("cand_Bz_OWNPV_CORRM","B_{0} CORRM [MeV]",mCorr_min,mCorr_max);
    x.setBins(nBins);

    RooDataHist dh_Dp ("dh_Dp","",RooArgList(x),h_Dp);
    RooDataHist dh_Dst("dh_Dst","",RooArgList(x),h_Dst);
    RooDataHist dh_Db ("dh_Db","",RooArgList(x),h_Db);
    RooDataHist dh_Dn ("dh_Dn","",RooArgList(x),h_Dn);
    RooDataHist dh_WS ("dh_WS","",RooArgList(x),hWS);

    RooHistPdf pdf_Dp ("pdf_Dp","",x,dh_Dp);
    RooHistPdf pdf_Dst("pdf_Dst","",x,dh_Dst);
    RooHistPdf pdf_Db ("pdf_Db","",x,dh_Db);
    RooHistPdf pdf_Dn ("pdf_Dn","",x,dh_Dn);
    RooHistPdf pdf_WS ("pdf_WS","",x,dh_WS);

    RooRealVar N_Dp ("N_Dp","",20000,0,1e7);
    RooRealVar N_Dst("N_Dst","",10000,0,1e7);

    RooRealVar N_Db ("N_Db","",1000,0,1e7);
    RooRealVar N_Dn ("N_Dn","",500,0,1e7);

    RooRealVar N_WS("N_WS","",0, 0, 1e7);

    RooAddPdf model(
        "model","",
        RooArgList(pdf_Dp,pdf_Dst,pdf_Db,pdf_Dn,pdf_WS),
        RooArgList(N_Dp,N_Dst,N_Db,N_Dn,N_WS) 
    );

    RooDataHist dh_data("dh_data","",RooArgList(x),hData);
    std::unique_ptr<RooFitResult> fr(    
        model.fitTo(dh_data,Extended(true), Save(true), SumW2Error(true))
    );

    TCanvas c("c","",900,700);
    RooPlot* frame = x.frame();

    dh_data.plotOn(frame, DataError(RooAbsData::SumW2));

    model.plotOn(frame, LineColor(kBlack), LineWidth(3));

    model.plotOn(frame, Components("pdf_WS"), LineColor(kGreen+2), LineWidth(2));
    model.plotOn(frame, Components("pdf_Dp"), LineColor(kRed), LineWidth(2));
    model.plotOn(frame, Components("pdf_Dst"), LineColor(kBlue), LineWidth(2));
    model.plotOn(frame, Components("pdf_Db"), LineColor(kMagenta), LineWidth(2));
    model.plotOn(frame, Components("pdf_Dn"), LineColor(kOrange+7), LineWidth(2));

    frame->SetTitle("Binned fit of cand_Bz_OWNPV_CORRM with sWeighted data and WS");
    frame->Draw();
    gPad->cd();

    TLegend leg(0.60,0.60,0.88,0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);

    leg.AddEntry(frame->getObject(0), "Data - sweighted", "pe");
    leg.AddEntry(frame->getObject(1), "Total fit", "l");
    leg.AddEntry(frame->getObject(2), "WS - sweighted", "l");
    leg.AddEntry(frame->getObject(3), "D+", "l");
    leg.AddEntry(frame->getObject(4), "D*+", "l");
    leg.AddEntry(frame->getObject(5), "D** narrow", "l");
    leg.AddEntry(frame->getObject(6), "D** broad", "l");

    leg.Draw();

    c.SaveAs("Fit_sw_Binned.png");
    fr->Print("v");
    fr->correlationMatrix().Print();
}

void plot_RS_sw() {

    RooRealVar x("cand_Bz_OWNPV_CORRM",
                 "Corrected mass [MeV]",
                 2000.,8000.);

    TFile *f = TFile::Open(WS.c_str());
    TTree *t = (TTree*)f->Get("DecayTree");

    TH1D *h_data = new TH1D("h_data",
                            "WS data;Corrected mass [MeV];Candidates",
                            nBins,2000,8000);
    h_data->Sumw2();
    t->Draw("cand_Bz_OWNPV_CORRM>>h_data","","goff");

    TH1D *h_sig = new TH1D("h_sig",
                           "Signal (sWeighted);Corrected mass [MeV];Candidates",
                           nBins,2000,8000);
    h_sig->Sumw2();
    t->Draw("cand_Bz_OWNPV_CORRM>>h_sig",
            "N_sig_sw",
            "goff");

    TH1D *h_bkg = new TH1D("h_bkg",
                           "Background (sWeighted);Corrected mass [MeV];Candidates",
                           nBins,2000,8000);
    h_bkg->Sumw2();
    t->Draw("cand_Bz_OWNPV_CORRM>>h_bkg",
            "N_bkg_sw",
            "goff");

    TCanvas *c = new TCanvas("c","sWeights",1000,800);

    h_data->SetLineColor(kBlack);
    h_data->SetMarkerStyle(20);

    h_sig->SetLineColor(kBlue);
    h_sig->SetMarkerColor(kBlue);

    h_bkg->SetLineColor(kRed);
    h_bkg->SetMarkerColor(kRed);

    h_data->Draw("E");
    h_sig->Draw("E SAME");
    h_bkg->Draw("E SAME");

    TLegend *leg = new TLegend(0.60,0.70,0.88,0.88);
    leg->AddEntry(h_data,"WS data (raw)","lep");
    leg->AddEntry(h_sig,"Combinatorial signal after sWeights","lep");
    leg->AddEntry(h_bkg,"Background after sWeights","lep");
    leg->SetBorderSize(0);
    leg->Draw();

    c->SaveAs("sweights_ws.png");

    std::cout << "Raw data integral: " << h_data->Integral() << std::endl;
    std::cout << "Signal sWeighted yield: " << h_sig->Integral() << std::endl;
    std::cout << "Background sWeighted yield: " << h_bkg->Integral() << std::endl;
}

void extra_tracks() {
    ROOT::RDataFrame df("DecayTree", data_1.c_str());

    auto df_tracks = df
        .Define("pass_mask",
            "EXTRA_PARTS_PROBNN_PI > 0.3 && EXTRA_PARTS_PT > 1000 && EXTRA_PARTS_ETA > 2 && EXTRA_PARTS_ETA");

    auto df_clean = df_tracks.Define("overlap_mask",
        [](ROOT::RVec<float> ep_px, ROOT::RVec<float> ep_py, ROOT::RVec<float> ep_pz,
           float d_K_px, float d_K_py, float d_K_pz,
           float d_pi1_px, float d_pi1_py, float d_pi1_pz, 
           float d_pi2_px, float d_pi2_py, float d_pi2_pz) {
            
            ROOT::RVec<int> is_clean(ep_px.size(), 1);
            float tol = 1;
            
            for(size_t i = 0; i < ep_px.size(); ++i) {
                if (std::abs(ep_px[i] - d_K_px) < tol && std::abs(ep_py[i] - d_K_py) < tol && std::abs(ep_pz[i] - d_K_pz) < tol) is_clean[i] = 0;
                if (std::abs(ep_px[i] - d_pi1_px) < tol && std::abs(ep_py[i] - d_pi1_py) < tol && std::abs(ep_pz[i] - d_pi1_pz) < tol) is_clean[i] = 0;
                if (std::abs(ep_px[i] - d_pi2_px) < tol && std::abs(ep_py[i] - d_pi2_py) < tol && std::abs(ep_pz[i] - d_pi2_pz) < tol) is_clean[i] = 0;
            }
            return is_clean;
            
        }, {"EXTRA_PARTS_PX", "EXTRA_PARTS_PY", "EXTRA_PARTS_PZ", 
            "cand_Dp_Kp_PX", "cand_Dp_Kp_PY", "cand_Dp_Kp_PZ",
            "cand_Dp_pim1_PX", "cand_Dp_pim1_PY", "cand_Dp_pim1_PZ", 
            "cand_Dp_pim2_PX", "cand_Dp_pim2_PY", "cand_Dp_pim2_PZ"});

    auto df_final_mask = df_clean.Define("final_mask", "pass_mask && overlap_mask");
    auto df_pass = df_final_mask.Filter("Sum(final_mask) > 0", "No overlap + passing extra part");

    auto df_filtered = df_pass
        .Define("FILT_EP_PX", "EXTRA_PARTS_PX[pass_mask]")
        .Define("FILT_EP_PY", "EXTRA_PARTS_PY[pass_mask]")
        .Define("FILT_EP_PZ", "EXTRA_PARTS_PZ[pass_mask]");

    float mass_pi = 139.57; 
    float mass_D_PDG = 1869.66; 
    
    ROOT::RDF::TH1DModel model_m_DEP("h_m_DEP", "EXTRA_PARTS peak; m(Dpi) - m(D) + m(Dpdg) [MeV]; Events", 100, 2100, 2600);
        
    auto h_m_DEP = df_filtered
        .Define("m_DEP", [mass_pi, mass_D_PDG](float Dpx, float Dpy, float Dpz, double D_M,
                                   ROOT::RVec<float> EP_px, ROOT::RVec<float> EP_py, ROOT::RVec<float> EP_pz) {
                
                RVec<float> masses;
                
                PxPyPzMVector pD_m(Dpx, Dpy, Dpz, D_M);
                PxPyPzEVector pD(pD_m);
                
                for (size_t i = 0; i < EP_px.size(); ++i) {
                    PxPyPzMVector pEP_m(EP_px[i], EP_py[i], EP_pz[i], mass_pi);
                    PxPyPzEVector pEP(pEP_m);
                    
                    PxPyPzEVector pDEP = pD + pEP;
                    
                    // m(Dpi) - m(D) + m(D_PDG)
                    float improved_mass = pDEP.M() - pD.M() + mass_D_PDG;
                    
                    masses.push_back(improved_mass);    
                }
                return masses;

        }, {"cand_Bz_Dp_PX", "cand_Bz_Dp_PY", "cand_Bz_Dp_PZ", "cand_Bz_Dp_M",
            "FILT_EP_PX", "FILT_EP_PY", "FILT_EP_PZ"})
        .Histo1D(model_m_DEP, "m_DEP");
    
    TCanvas* c1 = new TCanvas("c1", "D + Extra Track Mass", 800, 600);
    
    h_m_DEP->SetLineColor(kBlue+2);
    h_m_DEP->SetLineWidth(2);
    h_m_DEP->SetFillColorAlpha(kBlue, 0.1); 
    
    h_m_DEP->Draw("HIST");
    c1->Update(); 
    
    double y_min = h_m_DEP->GetMinimum() > 0 ? h_m_DEP->GetMinimum() : 0.1;
    double y_max = h_m_DEP->GetMinimum() * 2.2; 
    
    // Line for x = 2460
    TLine *line_2460 = new TLine(2460, y_min, 2460, y_max);
    line_2460->SetLineColor(kRed);
    line_2460->SetLineStyle(2); 
    line_2460->SetLineWidth(2);
    line_2460->Draw(); 

    // Line for x = 2280
    TLine *line_2280 = new TLine(2280, y_min, 2280, y_max);
    line_2280->SetLineColor(kGreen);
    line_2280->SetLineStyle(2); 
    line_2280->SetLineWidth(2);
    line_2280->Draw(); 
    
    c1->SaveAs("m_DEP_plot_RS.png");
}


void variable_difference_scan() {
    ROOT::EnableImplicitMT();

    std::string tree_name = "DecayTree";
    std::string sideband_cut = "cand_Bz_Dp_M < 1845 || cand_Bz_Dp_M > 1895";

    TFile* f_rs = TFile::Open(data_1.c_str(), "READ");
    TTree* tree_rs = (TTree*)f_rs->Get(tree_name.c_str());
    
    // Open WS file purely to validate branch existence
    TFile* f_ws = TFile::Open(WS.c_str(), "READ");
    TTree* tree_ws = (TTree*)f_ws->Get(tree_name.c_str());
    
    std::vector<std::string> variables;
    struct VarLimits { double min; double max; };
    std::map<std::string, VarLimits> limits;

    std::cout << "Validating branches and extracting limits safely..." << std::endl;

    TObjArray* branches = tree_rs->GetListOfBranches();
    long n_scan_events = std::min(tree_rs->GetEntries(), 50000LL);
    tree_rs->SetEstimate(n_scan_events);

    for (int i = 0; i < branches->GetEntries(); ++i) {
        TBranch* b_rs = (TBranch*)branches->At(i);
        std::string varName = b_rs->GetName();
        
        if (varName.find("cand_") != 0) continue; // no extra parts and other variables
        if (varName == "cand_Bz_Dp_M") continue;
        if (!tree_ws->GetBranch(varName.c_str())) continue;
        
        // Scalar check
        TLeaf* leaf = b_rs->GetLeaf(varName.c_str());
        if (!leaf || leaf->GetLen() != 1) continue; 

        tree_rs->Draw(varName.c_str(), "", "goff", n_scan_events);
        long rows = tree_rs->GetSelectedRows();
        
        if (rows > 0) {
            double* v1 = tree_rs->GetV1();
            double vmin = TMath::MinElement(rows, v1);
            double vmax = TMath::MaxElement(rows, v1);
            
            // Only add if we have valid limits
            if (!std::isnan(vmin) && !std::isnan(vmax)) {
                variables.push_back(varName);
                limits[varName] = {vmin, vmax};
            }
        }
    }
    f_rs->Close();
    f_ws->Close();

    std::cout << "Filtered down to " << variables.size() << " safe scalar variables.\n" << std::endl;

    ROOT::RDataFrame df_RS(tree_name, data_1);
    ROOT::RDataFrame df_WS(tree_name, WS);

    auto df0_RS = df_RS.Filter(sideband_cut);
    auto df0_WS = df_WS.Filter(sideband_cut);

    std::vector<ROOT::RDF::RResultPtr<::TH1D>> histos_RS;
    std::vector<ROOT::RDF::RResultPtr<::TH1D>> histos_WS;

    const int nx = 100;
    for (const auto& var : variables) {
        double vmin = limits[var].min;
        double vmax = limits[var].max;
        
        if (vmin == vmax) { vmin -= 1; vmax += 1; }

        histos_RS.push_back(df0_RS.Histo1D({("RS_"+var).c_str(), var.c_str(), nx, vmin, vmax}, var));
        histos_WS.push_back(df0_WS.Histo1D({("WS_"+var).c_str(), var.c_str(), nx, vmin, vmax}, var));
    }

    std::cout << "Executing RDataFrame event loops...\n" << std::endl;

    std::vector<std::pair<double, std::string>> diff_results;

    // Bin content difference of RS-WS cand_ variables. Max difference = 1, Min difference = 0.
    for (size_t i = 0; i < variables.size(); ++i) {
        TH1D* hRS = histos_RS[i].GetPtr();
        TH1D* hWS = histos_WS[i].GetPtr();
        
        if (hRS->GetEntries() == 0 || hWS->GetEntries() == 0) continue;

        hRS->Scale(1.0 / hRS->Integral());
        hWS->Scale(1.0 / hWS->Integral());

        double current_diff = 0.0;
        
        for (int bin = 1; bin <= hRS->GetNbinsX(); ++bin) {
            current_diff += std::abs(hRS->GetBinContent(bin) - hWS->GetBinContent(bin));
        }

        diff_results.push_back({current_diff, variables[i]});
    }

    std::sort(diff_results.rbegin(), diff_results.rend());

    std::cout << "==========================================\n";
    std::cout << "Top 10 Variables with Highest Shape Difference\n";
    std::cout << "==========================================\n";
    
    int print_limit = std::min(static_cast<int>(diff_results.size()), 20);
    
    for (int i = 0; i < print_limit; ++i) {
        std::cout << i + 1 << ". " << diff_results[i].second 
                  << " (Score: " << diff_results[i].first/2 << ")\n";
    }
    std::cout << "==========================================\n";
}