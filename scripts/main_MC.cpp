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
#include <RooRealVar.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooProduct.h>
#include <RooFitResult.h>
#include <RooCrystalBall.h>
#include "RooGaussian.h"
#include <RooPlot.h>
#include <RooJohnson.h>
#include <Math/GoFTest.h>
#include <RooPolynomial.h>
using namespace RooFit;
using namespace ROOT;
using namespace ROOT::Math;

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


/*
--------------------------------------------
    FILE NAMES
--------------------------------------------
*/

std::string MC_eff = "/Volumes/Rome Drive/Latest_INFN/INFN/root_files/MC/MC_2_eff_new.root"; // MC with efficiencies without any cuts
std::string MC = "/Volumes/Rome Drive/Latest_INFN/INFN/root_files/MC/MC_2_eff_with_mCorr.root"; // MC with efficiencies and cuts
std::string MC_templates_50 = "/Volumes/Rome Drive/Latest_INFN/INFN/root_files/MC/templates_MC/templates_mcorr_50.root";
std::string MC_templates_100 = "/Volumes/Rome Drive/Latest_INFN/INFN/root_files/MC/templates_MC/templates_mcorr_100.root";


/*
--------------------------------------------
    Cuts - Preselection - Plots
--------------------------------------------
*/

void Cuts_TreeAdd() {
    ROOT::RDataFrame df("DecayTree", MC_eff.c_str());

    auto df_filtered = df
        // Truth-matching
        .Filter([](int Bz_TRUEID, int MU_TRUEID, int MU_mother, int Dp_TRUEID) {
        return (std::abs(Bz_TRUEID) == 511 && std::abs(MU_mother) == 511 && std::abs(Dp_TRUEID) == 411 && std::abs(MU_TRUEID) == 13);
            }, {"Bz_TRUEID", "Bz_mup_TRUEID", "Bz_mup_MC_MOTHER_ID", "Bz_Dp_TRUEID"})
        
        // Flight direction cut
        .Filter("Bz_OWNPV_FD_Z > 2", "Flight direction Z")

        // B0 and D+ vertex quality
        .Filter("Bz_ENDVERTEX_CHI2DOF < 4", "B0 vertex chi2/ndof < 4")
        .Filter("Bz_Dp_ENDVERTEX_CHI2DOF < 4", "D vertex chi2/ndof < 4")

        // K+ track quality and PID
        .Filter("Dp_Kp_TRACKCHI2DOF < 2.5", "K+ track chi2/ndof < 2.5")
        .Filter("Dp_Kp_OWNPV_IPCHI2 > 10", "K+ IP chi2 > 10")
        .Filter("Dp_Kp_TRACKGHOSTPROB < 0.3", "K+ ghost prob < 0.3")
        .Filter("Dp_Kp_PROBNN_K > 0.3", "K+ probNN(Kp) > 0.3")

        // pi1
        .Filter("Dp_pim1_TRACKCHI2DOF < 2.5", "pi1 track chi2/ndof < 2.5")
        .Filter("Dp_pim1_OWNPV_IPCHI2 > 10", "pi1 IP chi2 > 10")
        .Filter("Dp_pim1_TRACKGHOSTPROB < 0.3", "pi1 ghost prob < 0.3")
        .Filter("Dp_pim1_PROBNN_PI > 0.3", "pi1 probNN(pi1) > 0.3")

        // pi2
        .Filter("Dp_pim2_TRACKCHI2DOF < 2.5", "pi2 track chi2/ndof < 2.5")
        .Filter("Dp_pim2_OWNPV_IPCHI2 > 10", "pi2 IP chi2 > 10")
        .Filter("Dp_pim2_TRACKGHOSTPROB < 0.3", "pi2 ghost prob < 0.3")
        .Filter("Dp_pim2_PROBNN_PI > 0.3", "pi2 probNN(pi2) > 0.3")

        // Muon
        .Filter("Bz_mup_PT > 1500", "mu pT > 1000 MeV")
        .Filter("Bz_mup_P > 6000", "mu p > 6000 MeV")
        .Filter("Bz_mup_TRACKCHI2DOF < 2.5", "mu track chi2/ndof < 2.5")
        .Filter("Bz_mup_OWNPV_IPCHI2 > 9", "mu IP chi2 > 9")
        .Filter("Bz_mup_PROBNN_MU > 0.3", "mu probNN(mu) > 0.3")
        
        // PV - composite hadrons
        .Filter("nPVs > 1", "nPVs > 1")
        .Filter("Bz_ENDVERTEX_POS_Z-Bz_Dp_ENDVERTEX_POS_Z < 0")
        .Filter("Bz_ETA < 5 && Bz_ETA > 2")
        .Filter("Bz_OWNPV_CORRM < 10000")
        .Filter("log(Bz_Dp_OWNPV_IP) > -3 && log(Bz_Dp_OWNPV_IP) < 5")
        .Filter("nPVs < 11")

        // Triggers
        .Filter("Bz_Hlt2SLB_B0ToDpMuNu_DpToKPiPiDecision_TOS == 1", "Hlt2 B0->D+μν TOS")
        .Filter("Bz_Hlt1TwoTrackMVADecision_TOS == 1", "Hlt1 TwoTrack MVA TOS")

        .Define("p_Dp",
                [](float px, float py, float pz) {
                    return ROOT::Math::XYZVector(px, py, pz);
                },
                {"Bz_Dp_PX","Bz_Dp_PY","Bz_Dp_PZ"})

        .Define("p_Kp",
                [](float px, float py, float pz) {
                    return ROOT::Math::XYZVector(px, py, pz);
                },
                {"Dp_Kp_PX","Dp_Kp_PY","Dp_Kp_PZ"})
        .Define("p_pi1",
                [](float px, float py, float pz) {
                    return ROOT::Math::XYZVector(px, py, pz);
                },
                {"Dp_pim1_PX","Dp_pim1_PY","Dp_pim1_PZ"})
        .Define("p_pi2",
                [](float px, float py, float pz) {
                    return ROOT::Math::XYZVector(px, py, pz);
                },
                {"Dp_pim2_PX","Dp_pim2_PY","Dp_pim2_PZ"})
        .Define("p_mup",
                [](float px, float py, float pz) {
                    return ROOT::Math::XYZVector(px, py, pz);
                },
                {"Bz_mup_PX","Bz_mup_PY","Bz_mup_PZ"})

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

        .Define("Bz_FD_rho", [](double svx, double svy, double svz) {
                XYZVector dir(svx, svy, svz);
                return dir.Rho();
            }, {"Bz_OWNPV_FD_X", "Bz_OWNPV_FD_Y", "Bz_OWNPV_FD_Z"}); // radial FD


    std::cout << "Entries before filter: " << *df.Count() << std::endl;
    std::cout << "Entries after filter: " << *df_filtered.Count() << std::endl;

    auto df2 = df_filtered
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
        .Filter("Bz_FD_rho < 4.8")  
        
        .Define("flight_dir", [](double svx, double svy, double svz,
                                    float pvx, float pvy, float pvz) {
            XYZVector dir(svx - pvx, svy - pvy, svz - pvz);
            return dir.Unit();
        }, {"Bz_ENDVERTEX_POS_X", "Bz_ENDVERTEX_POS_Y", "Bz_ENDVERTEX_POS_Z",
            "Bz_OWNPV_X", "Bz_OWNPV_Y", "Bz_OWNPV_Z"})

        // Perpendicular to FD for the visible mass D+mu
        .Define("p_perp_vis", [](float Dpx, float Dpy, float Dpz, float D_E, float MU_px, float MU_py, float MU_pz, 
                             float MU_E, const XYZVector &flight_dir) {
                PxPyPzEVector pD(Dpx, Dpy, Dpz, D_E);
                PxPyPzEVector pMu(MU_px, MU_py, MU_pz, MU_E);
                PxPyPzEVector pVis = pD + pMu;
                double p_par = pVis.Vect().Dot(flight_dir); // Momentum parallel to FD
                XYZVector p_perp_vec = pVis.Vect() - p_par * flight_dir; // Momentum only perpendicular to FD
                return p_perp_vec.R();              
        }, {"Bz_Dp_PX", "Bz_Dp_PY", "Bz_Dp_PZ", "Bz_Dp_ENERGY",
            "Bz_mup_PX", "Bz_mup_PY", "Bz_mup_PZ", "Bz_mup_ENERGY",
            "flight_dir"})
            
        // Perpendicular to FD only for the D mass
        .Define("p_perp_D", [](float Dpx, float Dpy, float Dpz, float D_E, const XYZVector &flight_dir) {
                PxPyPzEVector pD(Dpx, Dpy, Dpz, D_E);
                double p_par = pD.Vect().Dot(flight_dir);
                XYZVector p_perp_D_vec = pD.Vect() - p_par * flight_dir;
                return p_perp_D_vec.R();
        }, {"Bz_Dp_PX", "Bz_Dp_PY", "Bz_Dp_PZ", "Bz_Dp_ENERGY", "flight_dir"})

        .Define("mCorr", [](double Bz_M, double p_perp) {
            return std::sqrt(Bz_M * Bz_M + p_perp * p_perp) + p_perp;
        }, {"Bz_M", "p_perp_vis"})

        // Bbar or B - not used in the final code
        .Define("Bz_barness", [](int B_id) {
            return (B_id > 0 ? +1 : -1);
        }, {"Bz_TRUEID"})

        // Categorisation of decay products [31: BROAD D**] [32: NARROW D**]
        .Define("SL_category", [](int B_id, int D_id, int mu_id, int D_mother, int D_gdmother, int D_gdgdmother) {

                bool corr = ( (B_id == 511 && D_id == -411 && mu_id == -13) || // B0->D-mu+
                            (B_id == -511 && D_id == +411 && mu_id == 13) ); // B0bar->D+ mu-
                if (!corr) return 0;

                int m = std::abs(D_mother);
                int gm = std::abs(D_gdmother);
                int gg = std::abs(D_gdgdmother);

                if (m == 511) return 1; // B-> D
                if (m == 413 && gm == 511) return 2; // B-> D*

                // Narrow states
                if (m == 413 && (gm == 20413 || gm == 415) && gg == 511) return 31;  // [B0 -> D*_1+ -> D*+ (pi0) -> D+] OR [B0 -> (D*_2+) -> (D*+) + (pi0) -> D+]
                if (m == 415 && gm == 511) return 31; // [B0 -> (D*_2+) -> (D+) + (pi0) -> D+]

                // Broad states
                if (m == 413 && gm == 10413 && gg == 511) return 32; // [B0 -> D*_1+ (broad) -> D*+ (pi0) -> D+]
                if (m == 10411 && gm == 511) return 32; // [B0 -> (D*_0+) -> (D+) + (pi0) -> D+]

                return 0;
            },
            {"Bz_TRUEID","Bz_Dp_TRUEID","Bz_mup_TRUEID",
            "Bz_Dp_MC_MOTHER_ID","Bz_Dp_MC_GD_MOTHER_ID","Bz_Dp_MC_GD_GD_MOTHER_ID"});

    df2.Snapshot("DecayTree", MC.c_str());
}

void plot_CORRM() {
    ROOT::RDataFrame df_same("DecayTree", MC.c_str());

    // Unweighted
    auto h_Dp = df_same.Filter("SL_category == 1")
        .Histo1D({"h_Dp", "Bz_OWNPV_CORRM (unweighted);m_{corr} [MeV];Events",50,2000,8000}, "Bz_OWNPV_CORRM");
    auto h_Dst  = df_same.Filter("SL_category == 2")
        .Histo1D({"h_Dst", "Bz_OWNPV_CORRM (unweighted);m_{corr} [MeV];Events",50,2000,8000}, "Bz_OWNPV_CORRM");
    auto h_Dstst_narrow = df_same.Filter("SL_category == 31")
        .Histo1D({"h_Dstst_narrow","Bz_OWNPV_CORRM (unweighted);m_{corr} [MeV];Events",50,2000,8000}, "Bz_OWNPV_CORRM");
    auto h_Dstst_broad = df_same.Filter("SL_category == 32")
        .Histo1D({"h_Dstst_broad","Bz_OWNPV_CORRM (unweighted);m_{corr} [MeV];Events",50,2000,8000}, "Bz_OWNPV_CORRM");

    TCanvas *c1 = new TCanvas("c1", "Unweighted mCorr", 900, 700);

    h_Dp->SetLineColor(kBlue);
    h_Dst->SetLineColor(kRed);
    h_Dstst_broad->SetLineColor(kGreen+2);
    h_Dstst_narrow->SetLineColor(kMagenta+1);

    h_Dp->SetLineWidth(2);
    h_Dst->SetLineWidth(2);
    h_Dstst_broad->SetLineWidth(2);
    h_Dstst_narrow->SetLineWidth(2);

    h_Dp->Draw("HIST");
    h_Dst->Draw("HIST SAME");
    h_Dstst_broad->Draw("HIST SAME");
    h_Dstst_narrow->Draw("HIST SAME");

    auto leg1 = new TLegend(0.65, 0.60, 0.88, 0.80);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->SetMargin(0.20);
    leg1->AddEntry(h_Dp.GetPtr(), "D+", "l");
    leg1->AddEntry(h_Dst.GetPtr(), "D*", "l");
    leg1->AddEntry(h_Dstst_broad.GetPtr(), "D** broad", "l");
    leg1->AddEntry(h_Dstst_narrow.GetPtr(), "D** narrow", "l");
    leg1->Draw();

    c1->SaveAs("MC_truthmatched_unweighted.png");

    // Weighted
    auto h_Dp_w = df_same.Filter("SL_category == 1")
        .Histo1D({"h_Dp_w","Bz_OWNPV_CORRM (weighted);m_{corr} [MeV];Weighted Events",50,2000,8000}, "Bz_OWNPV_CORRM", "eff_Bz");
    auto h_Dst_w = df_same.Filter("SL_category == 2")
        .Histo1D({"h_Dst_w","Bz_OWNPV_CORRM (weighted);m_{corr} [MeV];Weighted Events",50,2000,8000}, "Bz_OWNPV_CORRM", "eff_Bz");
    auto h_Dstst_narrow_w = df_same.Filter("SL_category == 31")
        .Histo1D({"h_Dstst_narrow_w","Bz_OWNPV_CORRM (weighted);m_{corr} [MeV];Weighted Events",50,2000,8000}, "Bz_OWNPV_CORRM", "eff_Bz");
    auto h_Dstst_broad_w = df_same.Filter("SL_category == 32")
        .Histo1D({"h_Dstst_broad_w","Bz_OWNPV_CORRM (weighted);m_{corr} [MeV];Weighted Events",50,2000,8000}, "Bz_OWNPV_CORRM", "eff_Bz");

    TCanvas *c2 = new TCanvas("c2", "Weighted mCorr", 900, 700);

    h_Dp_w->SetLineColor(kBlue);
    h_Dst_w->SetLineColor(kRed);
    h_Dstst_narrow_w->SetLineColor(kMagenta+1);
    h_Dstst_broad_w->SetLineColor(kGreen+2);


    h_Dp_w->SetLineWidth(2);
    h_Dst_w->SetLineWidth(2);
    h_Dstst_narrow_w->SetLineWidth(2);
    h_Dstst_broad_w->SetLineWidth(2);

    h_Dp_w->Draw("HIST");
    h_Dst_w->Draw("HIST SAME");
    h_Dstst_narrow_w->Draw("HIST SAME");
    h_Dstst_broad_w->Draw("HIST SAME");

    auto leg2 = new TLegend(0.65, 0.60, 0.88, 0.80);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->SetMargin(0.20);
    leg2->AddEntry(h_Dp_w.GetPtr(), "D+", "l");
    leg2->AddEntry(h_Dst_w.GetPtr(), "D*", "l");
    leg2->AddEntry(h_Dstst_narrow_w.GetPtr(), "D** narrow", "l");
    leg2->AddEntry(h_Dstst_broad_w.GetPtr(), "D** broad", "l");
    leg2->Draw();

    c2->SaveAs("MC_truthmatched_weighted.png");

    // Unweighted D**
    TCanvas *c3 = new TCanvas("c3", "D** only (unweighted)", 900, 700);

    h_Dstst_broad->SetTitle("Bz_OWNPV_CORRM for D**p only (unweighted);m_{corr} [MeV];Events");
    h_Dstst_broad->GetYaxis()->SetRangeUser(0, 50);
    h_Dstst_broad->GetYaxis()->SetTitleOffset(1.3);
    h_Dstst_broad->GetXaxis()->SetTitleOffset(1.2);

    h_Dstst_broad->Draw("HIST");
    h_Dstst_narrow->Draw("HIST SAME");

    auto leg3 = new TLegend(0.60, 0.75, 0.85, 0.90);
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
    leg3->AddEntry(h_Dstst_broad.GetPtr(),  "D** broad", "l");
    leg3->AddEntry(h_Dstst_narrow.GetPtr(), "D** narrow", "l");
    leg3->Draw();

    c3->SaveAs("mCorr_unweighted_DssOnly.png");

    // Weighted D**
    TCanvas *c4 = new TCanvas("c4", "D** only (weighted)", 900, 700);

    h_Dstst_broad_w->SetTitle("Bz_OWNPV_CORRM for D**p only (weighted);m_{corr} [MeV];Weighted Events");
    h_Dstst_broad_w->GetYaxis()->SetRangeUser(0, 25);
    h_Dstst_broad_w->GetYaxis()->SetTitleOffset(1.3);
    h_Dstst_broad_w->GetXaxis()->SetTitleOffset(1.2);

    h_Dstst_broad_w->Draw("HIST");
    h_Dstst_narrow_w->Draw("HIST SAME");

    auto leg4 = new TLegend(0.60, 0.75, 0.85, 0.90);
    leg4->SetBorderSize(0);
    leg4->SetFillStyle(0);
    leg4->AddEntry(h_Dstst_broad_w.GetPtr(),  "D** broad", "l");
    leg4->AddEntry(h_Dstst_narrow_w.GetPtr(), "D** narrow", "l");
    leg4->Draw();

    c4->SaveAs("mCorr_weighted_DssOnly.png");

    // Saving template histograms to be used later if needed
    TFile fout("templates_mcorr.root", "RECREATE");
    h_Dp_w->Write();
    h_Dst_w->Write();
    h_Dstst_broad_w->Write();
    h_Dstst_narrow_w->Write();
    fout.Close();
}

void plot_pPerp() {
    ROOT::RDataFrame df0("DecayTree", MC.c_str());
    auto df = df0
                .Define("p_perp_D_GeV", "p_perp_D/1000.0")
                .Define("mCorr_GeV", "mCorr/1000.0");

    auto d_D = df.Filter("SL_category == 1"); // B0 -> Dp
    auto d_Dst = df.Filter("SL_category == 2"); // B0 -> Dp*
    auto d_Dss = df.Filter("SL_category == 31 || SL_category == 32"); // B0 -> Dp**

    const int nx = 25; const double xlo = 3.8, xhi = 5.9; // [GeV]
    const int ny = 25; const double ylo = 0, yhi = 2.5;

    auto h_D   = d_D.Histo2D({"h_D", "B^{0} -> D mu; m_{corr} [GeV]; p_{#perp}(D) [GeV/c]",
                                nx,xlo,xhi, ny,ylo,yhi}, "mCorr_GeV", "p_perp_D_GeV");
    auto h_Dst = d_Dst.Histo2D({"h_Dst", "B^{0} -> D* mu; m_{corr} [GeV]; p_{#perp}(D) [GeV/c]",
                                nx,xlo,xhi, ny,ylo,yhi}, "mCorr_GeV", "p_perp_D_GeV");
    auto h_Dss = d_Dss.Histo2D({"h_Dss", "B^{0} -> D** mu; m_{corr} [GeV]; p_{#perp}(D) [GeV/c]",
                                nx,xlo,xhi, ny,ylo,yhi}, "mCorr_GeV", "p_perp_D_GeV");

    auto normalize = [](TH2D* h){
        double m = h->GetMaximum();
        if (m > 0) h->Scale(1.0/m);
        h->SetMinimum(0.0);
        h->SetMaximum(1.0);
    };

    // D only
    {
        TCanvas* cD = new TCanvas("c_pperp_mcorr_D", "D only", 800, 700);
        cD->SetLeftMargin(0.13);
        cD->SetRightMargin(0.15);
        cD->SetBottomMargin(0.13);
        cD->SetTopMargin(0.08);

        normalize(h_D.GetPtr());
        h_D->Draw("COLZ");

        cD->SaveAs("pPerp_vs_mCorr_D_only.png");
        delete cD;
    }

    // D* D**
    {
        TCanvas* cDS = new TCanvas("c_pperp_mcorr_DstDss", "D* and D**", 1400, 700);
        cDS->Divide(2,1,0.001,0.001);

        auto draw_pad = [&](int pad, TH2D* h){
            cDS->cd(pad);
            gPad->SetLeftMargin(0.13);
            gPad->SetRightMargin(0.15);
            gPad->SetBottomMargin(0.13);
            gPad->SetTopMargin(0.08);
            normalize(h);
            h->Draw("COLZ");
        };

        draw_pad(1, h_Dst.GetPtr());  // left D*
        draw_pad(2, h_Dss.GetPtr());  // right D**

        cDS->SaveAs("pPerp_vs_mCorr_Ds.png");
        delete cDS;
    }
}

// For Siracuse comparison
void plot_CORRM_bins()
{
    TFile *f = TFile::Open(MC.c_str());
    TTree *t = (TTree*)f->Get("DecayTree");

    float Bz_TRUEPT, Bz_OWNPV_CORRM, Bz_TRUEETA;
    int SL_category;

    t->SetBranchAddress("Bz_TRUEPT", &Bz_TRUEPT);
    t->SetBranchAddress("Bz_TRUEETA", &Bz_TRUEETA);
    t->SetBranchAddress("Bz_OWNPV_CORRM", &Bz_OWNPV_CORRM);
    t->SetBranchAddress("SL_category", &SL_category);

    const int NBINS = 50;
    const double XMIN = 2000, XMAX = 8000;

    TString titles_pT[3] = {
        "p_{T} < 7.5 GeV",
        "7.5 < p_{T} < 15 GeV",
        "p_{T} > 15 GeV"
    };

    TString titles_eta[3] = {
        "2< eta < 3",
        "3<= eta < 4",
        "4<= eta < 5"
    };

    TCanvas *c = new TCanvas("c","",1800,600);
    c->Divide(3,1);

    TH1D *h[3][2];
    TPad *pads[3];

    for(int i = 0; i < 3; i++){
        pads[i] = (TPad*)c->cd(i+1);

        h[i][0] = new TH1D(Form("hDp_%d", i), "", NBINS, XMIN, XMAX);
        h[i][1] = new TH1D(Form("hDst_%d", i), "", NBINS, XMIN, XMAX);

        h[i][0]->SetLineColor(kBlue);
        h[i][1]->SetLineColor(kRed);
        h[i][0]->SetLineWidth(2);
        h[i][1]->SetLineWidth(2);

        h[i][0]->SetDirectory(nullptr);
        h[i][1]->SetDirectory(nullptr);

        h[i][0]->SetTitle(titles_eta[i]);
        h[i][1]->SetTitle(titles_eta[i]);

        h[i][0]->GetXaxis()->SetTitle("m_{corr} [MeV]");
        h[i][0]->GetYaxis()->SetTitle("Events");

        h[i][0]->Draw("HIST");
        h[i][1]->Draw("HIST SAME");

        auto *leg = new TLegend(0.7,0.75,0.9,0.9);
        leg->AddEntry(h[i][0],"Dp","l");
        leg->AddEntry(h[i][1],"Dst","l");
        leg->Draw();
    }

    c->Update();

    Long64_t n = t->GetEntries();
    for(Long64_t i = 0; i < n; i++){
        t->GetEntry(i);

        int ieta = -1;
        if (Bz_TRUEETA > 2 && Bz_TRUEETA < 3)      ieta = 0;
        else if (Bz_TRUEETA >= 3 && Bz_TRUEETA < 4) ieta = 1;
        else if (Bz_TRUEETA >= 4 && Bz_TRUEETA < 5) ieta = 2;
        else continue;
 
        if(SL_category == 1) h[ieta][0]->Fill(Bz_OWNPV_CORRM);
        if(SL_category == 2) h[ieta][1]->Fill(Bz_OWNPV_CORRM);
    }

    c->Modified();
    c->Update();
    f->Close();
}

void KolmogorovSmirnov() {
    TFile* f = TFile::Open(MC.c_str());
    TTree* t = (TTree*)f->Get("DecayTree");

    int SL_category;
    float Bz_TRUEPT, Bz_OWNPV_CORRM, Bz_TRUEETA;

    t->SetBranchAddress("Bz_TRUEPT", &Bz_TRUEPT);
    t->SetBranchAddress("Bz_OWNPV_CORRM", &Bz_OWNPV_CORRM);
    t->SetBranchAddress("Bz_TRUEETA", &Bz_TRUEETA);
    t->SetBranchAddress("SL_category", &SL_category);

    std::vector<double> etaBin[3];

    Long64_t n = t->GetEntries();
    for (Long64_t i = 0; i < n; i++) {
        t->GetEntry(i);

        if (SL_category != 2) continue;

        int ieta = -1;

        // Three eta bins: 2-3, 3-4, 4-5
        if (Bz_TRUEPT < 7500) ieta = 0;
        else if (Bz_TRUEPT >= 7500 && Bz_TRUEPT < 15000) ieta = 1;
        else if (Bz_TRUEPT > 15000) ieta = 2;
        else continue;

        etaBin[ieta].push_back((double)Bz_OWNPV_CORRM);
    }

    // Kolmogorov Smirnov tests between the two consecutive bins
    for (int i = 0; i < 2; i++) {
        const auto& sampleA = etaBin[i];
        const auto& sampleB = etaBin[i+1];

        ROOT::Math::GoFTest gof(sampleA.size(), sampleA.data(), sampleB.size(), sampleB.data());

        double pvalue = gof.KolmogorovSmirnov2SamplesTest("p");
        double Dstat  = gof.KolmogorovSmirnov2SamplesTest("t");

        std::cout << "Kolmogorov-Smirnov test pT bin " << i << " vs " << i+1
                  << " | p-value = " << pvalue
                  << ", D = " << Dstat << std::endl;
    }
        const auto& sampleA = etaBin[0];
        const auto& sampleB = etaBin[2];

        ROOT::Math::GoFTest gof(sampleA.size(), sampleA.data(), sampleB.size(), sampleB.data());

        double pvalue = gof.KolmogorovSmirnov2SamplesTest("p");
        double Dstat  = gof.KolmogorovSmirnov2SamplesTest("t");

        std::cout << "Kolmogorov-Smirnov test pT bin " << 0 << " vs " << 2
                  << " | p-value = " << pvalue
                  << ", D = " << Dstat << std::endl;

    // Graph of CDFs for bins 1,2
    std::sort(etaBin[1].begin(), etaBin[1].end());
    std::sort(etaBin[2].begin(), etaBin[2].end());

    TGraph* g1 = new TGraph(etaBin[1].size());
    TGraph* g2 = new TGraph(etaBin[2].size());
    for (size_t i=0;i<etaBin[1].size();i++) g1->SetPoint(i, etaBin[1][i], double(i+1)/etaBin[1].size());
    for (size_t i=0;i<etaBin[2].size();i++) g2->SetPoint(i, etaBin[2][i], double(i+1)/etaBin[2].size());

    TCanvas* c = new TCanvas("c", "CDFs", 800, 600);
    g1->SetLineColor(kBlue);
    g1->SetLineWidth(2);
    g1->SetTitle("Comparison of CDFs for Kolmogorov-Smirnov");
    g1->GetXaxis()->SetTitle("Bz_OWNPV_CORRM");
    g1->GetYaxis()->SetTitle("CDF");
    g1->Draw("AL");

    g2->SetLineColor(kRed);
    g2->SetLineWidth(2);
    g2->Draw("L SAME");

    auto legend = new TLegend(0.65, 0.15, 0.88, 0.35);
    legend->AddEntry(g1, "Bin1: 7.5GeV <= pT < 15GeV", "l");
    legend->AddEntry(g2, "Bin2: pT > 15GeV", "l");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.03);
    legend->Draw();

    f->Close();
}



/*
--------------------------------------------
    Sidebands Dp_M checks
--------------------------------------------
*/

void SB_check_MC() {
    ROOT::RDataFrame dfRS("DecayTree", MC.c_str());
    auto hRS_tmp = dfRS.Filter("SL_category != 0")
                       .Filter("Bz_Dp_M < 1840 || Bz_Dp_M > 1900")
                       .Histo1D({"hRS",";m_{corr} [MeV];A.U.",200,2000,8000}, "Bz_OWNPV_CORRM");
    TH1D* hRS = (TH1D*)hRS_tmp.GetPtr()->Clone("hRS");
    hRS->Scale(1.0 / hRS->GetMaximum());
    hRS->SetDirectory(nullptr);

    TCanvas* c1 = new TCanvas("c", "mCorr MC", 900, 700);

    hRS->SetTitle("mCorr for Dp sidebands");
    hRS->SetLineColor(kBlue); hRS->SetLineWidth(2);

    hRS->Draw("HIST");
    hRS->SetMaximum(1.2);

    c1->SaveAs("Dp_sidebands_mCorr.png");
}


void fit_Dp_M_MC() {
    ROOT::RDataFrame df("DecayTree", MC.c_str());
    auto h_Dp_M_tmp = df.Filter("SL_category != 0") // Selecting signal
                        .Histo1D({"h_Dp_M", "D+ mass;Mass [MeV];Entries", 50, 1820, 1920}, "Bz_Dp_M");
    TH1D* h_Dp_M = (TH1D*)h_Dp_M_tmp.GetPtr()->Clone("h_Dp_M");
    h_Dp_M->SetDirectory(nullptr);

    for (int i = 1; i <= h_Dp_M->GetNbinsX(); ++i) {
        if (h_Dp_M->GetBinContent(i) <= 0) h_Dp_M->SetBinContent(i, 1e-6);
    }

    RooRealVar x("x", "m", 1820, 1920);
    x.setBins(h_Dp_M->GetNbinsX());

    RooDataHist dh_Dp_M("dh_Dp_M", "dh_Dp_M", RooArgList(x), h_Dp_M);

    RooRealVar mean("mean", "DCB mean", 1870, 1860, 1880);
    RooRealVar sigma("sigma", "DCB width", 8, 2, 15);
    RooRealVar alphaL("alphaL", "Left alpha", 1.7, 0.1, 4);
    RooRealVar nL("nL", "Left n", 2, 0.1, 20);
    RooRealVar alphaR("alphaR", "Right alpha", 3, 0.1, 4);
    RooRealVar nR("nR", "Right n", 3, 0.1, 50);

    RooCrystalBall pdf_DCB("pdf_DCB", "Double Crystal Ball",
                           x,
                           mean, sigma,
                           alphaL, nL,
                           alphaR, nR);

    std::unique_ptr<RooFitResult> fr(pdf_DCB.fitTo(dh_Dp_M, Save(true), SumW2Error(true)));

    TCanvas cMC("cMC", "Double Crystal Ball Fit", 900, 700);
    RooPlot* frameMC = x.frame();

    dh_Dp_M.plotOn(frameMC, DataError(RooAbsData::SumW2));
    pdf_DCB.plotOn(frameMC, LineColor(kRed), LineWidth(3));

    frameMC->SetTitle("MC Mass of D+ DCB fit");
    frameMC->SetXTitle("m_{D+} [MeV]");
    frameMC->SetYTitle("Entries");
    frameMC->SetMinimum(0.1);
    gPad->SetLogy(); // for better visualisation
    frameMC->Draw();

    cMC.SaveAs("MC_DCB_fit.png");

    fr->Print("v");
}


/*
--------------------------------------------
    Unbinned Fit - MC Components
--------------------------------------------
*/


void fit_DCB() {

    ROOT::RDataFrame df("DecayTree", MC.c_str());

    // These have to be fitted one by one, uncomment whichever component you need
    std::vector<std::pair<std::string, std::string>> comps = {
        {"Dp",      "SL_category == 1"} 
        //{"Dst",     "SL_category == 2"}
        //{"Dstst31", "SL_category == 31"}
        //{"Dstst32", "SL_category == 32"}
    };

    RooRealVar mCorr("mCorr", "Corrected mass", 2250, 8000);
    RooRealVar eff_Bz("eff_Bz", "weight", 0, 1);

    for (auto &c : comps) {

        TString name = c.first;
        TString cut  = c.second;

        std::cout << "\n==========================\n";
        std::cout << "  Fitting category: " << name << "\n";
        std::cout << "==========================\n";

        auto df_cut = df.Filter(cut.Data()); // Selecting the component we study

        TString tmpfile = Form("tmp_%s.root", name.Data());
        TString tmptree = "DecayTree";

        df_cut.Snapshot(tmptree.Data(), tmpfile.Data()); // Create a temp file

        TFile f(tmpfile);
        TTree *t = (TTree*)f.Get(tmptree);

        RooDataSet data(
            "data", "Weighted MC sample",
            RooArgSet(mCorr, eff_Bz),
            RooFit::Import(*t),
            RooFit::WeightVar("eff_Bz")
        );

        std::cout << "Loaded entries (weighted): " << data.sumEntries() << "\n";

        RooRealVar mean ("mean",  "Mean", 4700, 4500, 5500);
        RooRealVar sigma("sigma", "Sigma", 150, 50, 300);

        RooRealVar alphaL("alphaL", "Left alpha", 1.5, 0, 5);
        RooRealVar nL    ("nL",     "Left n",     5.0, 0.5, 115);

        RooRealVar alphaR("alphaR", "Right alpha", 1.5, 0.5, 5);
        RooRealVar nR    ("nR",     "Right n",     5.0, 0.5, 50);

        RooCrystalBall dcb(
            "dcb", "Double Crystal Ball",
            mCorr,
            mean, sigma,
            alphaL, nL,
            alphaR, nR
        );

        dcb.fitTo(data,
                  RooFit::Save(true),
                  RooFit::SumW2Error(true));

        TCanvas *c1 = new TCanvas(Form("c_%s", name.Data()), "DCB Fit", 900, 700);
        RooPlot *frame = mCorr.frame(Title(Form("DCB Fit for %s", name.Data())));

        data.plotOn(frame, RooFit::Name("data"));
        dcb.plotOn(frame, RooFit::Name("dcb"));

        frame->Draw();
        c1->SaveAs(Form("Fit_%s.png", name.Data()));

        f.Close();
    }
}


void fit_Gauss() {
    ROOT::RDataFrame df("DecayTree", MC.c_str());

    // Same philosophy as the fit_DCB() here selecting for a combined D**; use whichever you like to study
    std::vector<std::pair<std::string, std::string>> comps = {
        {"Dstst_broad", "SL_category == 32"}
    };

    RooRealVar mCorr("Bz_OWNPV_CORRM", "Corrected mass", 2250, 8000);
    RooRealVar eff_Bz("eff_Bz", "weight", 0, 1);

    for (auto &c : comps) {

        TString name = c.first;
        TString cut  = c.second;

        std::cout << "\n==========================\n";
        std::cout << "  Fitting category: " << name << "\n";
        std::cout << "==========================\n";

        auto df_cut = df.Filter(cut.Data());

        TString tmpfile = Form("tmp_%s.root", name.Data());
        TString tmptree = "DecayTree";
        df_cut.Snapshot(tmptree.Data(), tmpfile.Data());

        TFile f(tmpfile);
        TTree *t = (TTree*)f.Get(tmptree);

        RooDataSet data(
            "data", "Weighted MC sample",
            RooArgSet(mCorr, eff_Bz),
            RooFit::Import(*t),
            RooFit::WeightVar("eff_Bz")
        );

        std::cout << "Loaded entries (weighted): " << data.sumEntries() << "\n";

        // Gaussian
        RooRealVar mean ("mean",  "Mean", 4100, 3500, 5000);
        RooRealVar sigma("sigma", "Sigma", 400, 10, 10000);
        //sigma.setConstant();

        RooGaussian gauss(
            "gauss",
            "Gaussian PDF",
            mCorr,
            mean,
            sigma
        );

        gauss.fitTo(data,
                    RooFit::Save(true),
                    RooFit::SumW2Error(true));

        TCanvas *c1 = new TCanvas(Form("c_%s", name.Data()), "Gaussian Fit", 900, 700);
        RooPlot *frame = mCorr.frame(Title(Form("Gaussian Fit for %s", name.Data())));

        data.plotOn(frame, RooFit::Name("data"));
        gauss.plotOn(frame, RooFit::Name("gauss"));

        frame->Draw();
        c1->SaveAs(Form("Fit_%s.png", name.Data()));

        f.Close();
    }
}


void fit_Johnson() {
    ROOT::RDataFrame df("DecayTree", MC.c_str());

    std::vector<std::pair<std::string, std::string>> comps = {
        //{"Dp",      "SL_category == 1"},
        {"Dst",     "SL_category == 2"}
        //{"Dstst31", "SL_category == 31"}
        //{"Dstst32", "SL_category == 32"}
    };

    RooRealVar mCorr("mCorr", "Corrected mass", 2000, 9000);
    RooRealVar eff_Bz("eff_Bz", "weight", 0, 1);

    for (auto &c : comps) {

        TString name = c.first;
        TString cut  = c.second;

        std::cout << "\n==========================\n";
        std::cout << "  Fitting category: " << name << "\n";
        std::cout << "==========================\n";

        auto df_cut = df.Filter(cut.Data());

        TString tmpfile = Form("tmp_%s.root", name.Data());
        TString tmptree = "DecayTree";

        df_cut.Snapshot(tmptree.Data(), tmpfile.Data());

        TFile f(tmpfile);
        TTree *t = (TTree*)f.Get(tmptree);

        RooDataSet data(
            "data", "Weighted MC sample",
            RooArgSet(mCorr, eff_Bz),
            RooFit::Import(*t),
            RooFit::WeightVar("eff_Bz")
        );

        std::cout << "Loaded entries (weighted): " << data.sumEntries() << "\n";

        RooRealVar mu    ("mu",    "Location", 5339, 4500, 5500);
        RooRealVar lambda("lambda","Scale",    200,  50,   800);
        RooRealVar gamma ("gamma", "Gamma",    0.0, -5.0,  5.0);
        RooRealVar delta ("delta", "Delta",    1.0,  0.2,  5.0);

        RooJohnson johnson(
            "johnson", "Johnson",
            mCorr,
            mu,
            lambda,
            gamma,
            delta
        );

        johnson.fitTo(
            data,
            RooFit::Save(true),
            RooFit::SumW2Error(true)
        );

        TCanvas *c1 = new TCanvas(Form("c_%s", name.Data()),
                                 "Johnson Fit", 900, 700);

        RooPlot *frame = mCorr.frame(
            Title(Form("Johnson SU Fit for %s", name.Data()))
        );

        data.plotOn(frame, RooFit::Name("data"));
        johnson.plotOn(frame, RooFit::Name("johnson"));

        frame->Draw();
        c1->SaveAs(Form("Fit_%s_johnson.pdf", name.Data()));

        f.Close();
    }
}


/*
--------------------------------------------
    sWeights - Extraction
--------------------------------------------
*/


void sweights_MC()
{
    TFile fmc("Dp_M_template.root","READ");
    TH1D* h_Dp_M = (TH1D*)fmc.Get("h_Dp_M_temp__x");
    h_Dp_M->SetDirectory(nullptr);
    fmc.Close();

    ROOT::RDataFrame df("DecayTree", "Data_1_cuts.root");
    ROOT::RDF::TH1DModel model("hData", ";m_{D+} [MeV];Entries", 100, 1820, 1920);
    auto hDataR = df.Histo1D(model, "cand_Bz_Dp_M");

    TH1D* hData = (TH1D*) hDataR.GetPtr()->Clone("hData");
    hData->SetDirectory(nullptr);

    RooRealVar x("x", "m_{D+}", 1820, 1920);
    x.setBins(100);

    RooDataHist dh_MC("dh_MC", "dh_MC", RooArgList(x), h_Dp_M);
    
    RooHistPdf pdf_MC("pdf_MC", "MC template PDF", x, dh_MC);
    RooRealVar N_MC("N_MC", "Signal yield", hData->Integral()/2, 0, 1e9);

    RooPolynomial pdf_bkg("pdf_bkg", "Background PDF", x, RooArgList());
    RooRealVar N_bkg("N_bkg", "Background yield", hData->Integral()/5, 0, 1e9);

    RooAddPdf model0("model0", "MC template + constant polynomial", 
                    RooArgList(pdf_MC, pdf_bkg), 
                    RooArgList(N_MC, N_bkg));

    RooDataHist dh_data("dh_data", "Data", RooArgList(x), hData);

    auto fitres = model0.fitTo(dh_data, Save(true), Extended(true));

    TCanvas c("c", "Fit", 900, 700);
    RooPlot* frame = x.frame();

    dh_data.plotOn(frame, DataError(RooAbsData::SumW2));
    model0.plotOn(frame, LineColor(kRed), LineWidth(3)); // total fit
    model0.plotOn(frame, LineColor(kBlue), LineStyle(kDotted), Components(pdf_bkg)); // background
    model0.plotOn(frame, LineColor(kGreen+1), LineStyle(kDotted), Components(pdf_MC)); // signal

    frame->SetTitle("Fit of m_{D+} using MC template + const poly");
    frame->SetXTitle("m_{D+} [MeV]");
    frame->SetYTitle("Entries");
    frame->Draw();

    TLegend legend(0.65, 0.65, 0.88, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.AddEntry(frame->getObject(0), "Data", "pe"); 
    legend.AddEntry(frame->getObject(1), "Total fit", "l");
    legend.AddEntry(frame->getObject(2), "Background", "l");
    legend.AddEntry(frame->getObject(3), "MC Signal", "l");
    legend.Draw();

    c.SaveAs("Dp_M_fit_poly.png");

}

/*
--------------------------------------------
    Extra Tracks 
--------------------------------------------
*/

void extra_tracks() {
    ROOT::RDataFrame df("DecayTree", MC.c_str());

    auto df_tracks = df
        .Define("pass_mask",
            "EXTRA_PARTS_PROBNN_PI > 0.3 && EXTRA_PARTS_PT > 1000 && EXTRA_PARTS_ETA > 2 && EXTRA_PARTS_ETA < 5"); // && (EXTRA_PARTS_CHARGE * Bz_Dp_CHARGE ==+1)"); // WS=+1, RS=-1

    auto df_clean = df_tracks.Define("overlap_mask",
        [](ROOT::RVec<float> ep_px, ROOT::RVec<float> ep_py, ROOT::RVec<float> ep_pz,
           float d_K_px, float d_K_py, float d_K_pz,
           float d_pi1_px, float d_pi1_py, float d_pi1_pz, 
           float d_pi2_px, float d_pi2_py, float d_pi2_pz) {
            
            ROOT::RVec<int> is_clean(ep_px.size(), 1);
            float tol = 0.0001;
            
            for(size_t i = 0; i < ep_px.size(); ++i) {
                if (std::abs(ep_px[i] - d_K_px) < tol && std::abs(ep_py[i] - d_K_py) < tol && std::abs(ep_pz[i] - d_K_pz) < tol) is_clean[i] = 0;
                if (std::abs(ep_px[i] - d_pi1_px) < tol && std::abs(ep_py[i] - d_pi1_py) < tol && std::abs(ep_pz[i] - d_pi1_pz) < tol) is_clean[i] = 0;
                if (std::abs(ep_px[i] - d_pi2_px) < tol && std::abs(ep_py[i] - d_pi2_py) < tol && std::abs(ep_pz[i] - d_pi2_pz) < tol) is_clean[i] = 0;
            }
            return is_clean;
            
        }, {"EXTRA_PARTS_TRUEPX", "EXTRA_PARTS_TRUEPY", "EXTRA_PARTS_TRUEPZ", 
            "Dp_Kp_TRUEPX", "Dp_Kp_TRUEPY", "Dp_Kp_TRUEPZ",
            "Dp_pim1_TRUEPX", "Dp_pim1_TRUEPY", "Dp_pim1_TRUEPZ", 
            "Dp_pim2_TRUEPX", "Dp_pim2_TRUEPY", "Dp_pim2_TRUEPZ"});

    auto df_final_mask = df_clean.Define("final_mask", "pass_mask && overlap_mask");
    auto df_pass = df_final_mask.Filter("Sum(final_mask) > 0", "No overlap + passing extra part");

    auto df_filtered = df_pass
        .Define("FILT_EP_PX", "EXTRA_PARTS_PX[final_mask]")
        .Define("FILT_EP_PY", "EXTRA_PARTS_PY[final_mask]")
        .Define("FILT_EP_PZ", "EXTRA_PARTS_PZ[final_mask]");

    auto df_cat = df_filtered
        .Define("category",
            [](int B_id, int D_id, int mu_id,
               int D_mother, int D_gdmother, int D_gdgdmother) {

                int m  = std::abs(D_mother);
                int gm = std::abs(D_gdmother);
                int gg = std::abs(D_gdgdmother);
                
                // Category 1: D_2*(2460) as mother
                if (m == 415) return 1;
                
                // Category 2: D_1 states (2420)
                if (gm == 10413 || gm == 20413 || m == 10413 || m == 20413) return 2;                
                
                // Category 3: D_2*(2460) as grandmother
                if (gm == 415) return 3;

                // Category 4: Direct B -> D pi
                if (m == 511) return 4;

                // Category 5: Direct B -> D* pi
                if (m == 413 && gm == 511) return 5;

                return 0;
            },
            {"Bz_TRUEID", "Bz_Dp_TRUEID", "Bz_mup_TRUEID",
             "Bz_Dp_MC_MOTHER_ID",
             "Bz_Dp_MC_GD_MOTHER_ID",
             "Bz_Dp_MC_GD_GD_MOTHER_ID"});

    float mass_pi = 139.57;
    float mass_D_PDG = 1869.66;

    ROOT::RDF::TH1DModel model_m_DEP("h_m_DEP", "EXTRA_PARTS peak; m(Dpi) - m(D) + m(Dpdg) [MeV]; Events", 100, 2100, 2600);

    auto df_mdep = df_cat
        .Define("m_DEP", [mass_pi, mass_D_PDG](float Dpx, float Dpy, float Dpz, double D_M,
                                   ROOT::RVec<float> EP_px, ROOT::RVec<float> EP_py, ROOT::RVec<float> EP_pz) {

                RVec<float> masses;

                PxPyPzMVector pD_m(Dpx, Dpy, Dpz, D_M);
                PxPyPzEVector pD(pD_m);

                for (size_t i = 0; i < EP_px.size(); ++i) {
                    PxPyPzMVector pEP_m(EP_px[i], EP_py[i], EP_pz[i], mass_pi);
                    PxPyPzEVector pEP(pEP_m);

                    PxPyPzEVector pDEP = pD + pEP;

                    float improved_mass = pDEP.M() - pD.M() + mass_D_PDG;
                    masses.push_back(improved_mass);
                }

                std::sort(masses.begin(), masses.end());

                std::vector<float> unique_masses;
                for (float m : masses) {
                    if (unique_masses.empty() || std::abs(m - unique_masses.back()) >= 0.0001) {
                        unique_masses.push_back(m);
                    }
                }
                return unique_masses;

        }, {"Bz_Dp_PX", "Bz_Dp_PY", "Bz_Dp_PZ", "Bz_Dp_M",
            "FILT_EP_PX", "FILT_EP_PY", "FILT_EP_PZ"});

    auto h_cat0 = df_mdep.Filter("category == 0").Histo1D(model_m_DEP, "m_DEP");
    auto h_cat1 = df_mdep.Filter("category == 1").Histo1D(model_m_DEP, "m_DEP");
    auto h_cat2 = df_mdep.Filter("category == 2").Histo1D(model_m_DEP, "m_DEP");
    auto h_cat3 = df_mdep.Filter("category == 3").Histo1D(model_m_DEP, "m_DEP");
    auto h_cat4 = df_mdep.Filter("category == 4").Histo1D(model_m_DEP, "m_DEP");
    auto h_cat5 = df_mdep.Filter("category == 5").Histo1D(model_m_DEP, "m_DEP");
    auto h_cat6 = df_mdep.Histo1D(model_m_DEP, "m_DEP");

    df_mdep.Filter("category == 0").Display({"m_DEP"}, 20)->Print(); // for debug of clones

    TCanvas* c = new TCanvas("c_3x2", "Categories 3x2", 1500, 1000);
    c->Divide(3, 2);

    auto draw_pad = [](TH1D* h, int color, const char* labelText) {
        h->SetLineColor(color);
        h->SetLineWidth(2);
        h->Draw("HIST");

        double y_min = h->GetMinimum() > 0 ? h->GetMinimum() : 0.1;
        double y_max = h->GetMaximum() * 1.2;

        TLine *line_2460 = new TLine(2460, y_min, 2460, y_max);
        line_2460->SetLineColor(kRed);
        line_2460->SetLineStyle(2);
        line_2460->Draw();

        TLine *line_2280 = new TLine(2280, y_min, 2280, y_max);
        line_2280->SetLineColor(kGreen+2);
        line_2280->SetLineStyle(2);
        line_2280->Draw();

        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(42);
        latex.SetTextSize(0.06);
        latex.DrawLatex(0.15, 0.82, labelText); 
    };

    c->cd(1);
    draw_pad(h_cat6.GetPtr(), kBlack, "Other ancestry");

    c->cd(2);
    draw_pad(h_cat1.GetPtr(), kRed, "D_{2}*+ (2460)");

    c->cd(3);
    draw_pad(h_cat2.GetPtr(), kBlue, "D_{1} (2420)");

    c->cd(4);
    draw_pad(h_cat3.GetPtr(), kGreen+2, "D_{2}*+ as gm");

    c->cd(5);
    draw_pad(h_cat4.GetPtr(), kYellow+2, "B0#rightarrowD+");

    c->cd(6);
    draw_pad(h_cat5.GetPtr(), kCyan-1, "B0#rightarrowD*#pi^{0}#rightarrowD+");

    c->SaveAs("m_DEP_WS_FINAL.png");
}


void extra_tracks_bp() {
    ROOT::RDataFrame df("DecayTree", "MC_Bp_new_eff.root");
    
    auto df_tracks = df
        .Define("pass_mask",
                "EXTRA_PARTS_PROBNN_PI > 0.3 && EXTRA_PARTS_PT > 1000 && EXTRA_PARTS_ETA > 2 && EXTRA_PARTS_ETA < 5 && (abs(EXTRA_PARTS_MC_MOTHER_ID) == 511 || abs(EXTRA_PARTS_MC_GD_MOTHER_ID) == 511 || abs(EXTRA_PARTS_MC_GD_GD_MOTHER_ID) == 511)"); // && (EXTRA_PARTS_CHARGE * Bz_Dp_CHARGE ==+1)"); // WS=+1, RS=-1
    
    auto df_clean = df_tracks.Define("overlap_mask",
                                     [](ROOT::RVec<float> ep_px, ROOT::RVec<float> ep_py, ROOT::RVec<float> ep_pz,
                                        float d_K_px, float d_K_py, float d_K_pz,
                                        float d_pi1_px, float d_pi1_py, float d_pi1_pz,
                                        float d_pi2_px, float d_pi2_py, float d_pi2_pz) {
        
        ROOT::RVec<int> is_clean(ep_px.size(), 1);
        float tol = 0.0001;
        
        for(size_t i = 0; i < ep_px.size(); ++i) {
            if (std::abs(ep_px[i] - d_K_px) < tol && std::abs(ep_py[i] - d_K_py) < tol && std::abs(ep_pz[i] - d_K_pz) < tol) is_clean[i] = 0;
            if (std::abs(ep_px[i] - d_pi1_px) < tol && std::abs(ep_py[i] - d_pi1_py) < tol && std::abs(ep_pz[i] - d_pi1_pz) < tol) is_clean[i] = 0;
            if (std::abs(ep_px[i] - d_pi2_px) < tol && std::abs(ep_py[i] - d_pi2_py) < tol && std::abs(ep_pz[i] - d_pi2_pz) < tol) is_clean[i] = 0;
        }
        return is_clean;
        
    }, {"EXTRA_PARTS_TRUEPX", "EXTRA_PARTS_TRUEPY", "EXTRA_PARTS_TRUEPZ",
        "Dp_Kp_TRUEPX", "Dp_Kp_TRUEPY", "Dp_Kp_TRUEPZ",
        "Dp_pim1_TRUEPX", "Dp_pim1_TRUEPY", "Dp_pim1_TRUEPZ",
        "Dp_pim2_TRUEPX", "Dp_pim2_TRUEPY", "Dp_pim2_TRUEPZ"});
    
    auto df_final_mask = df_clean.Define("final_mask", "pass_mask && overlap_mask");
    auto df_pass = df_final_mask.Filter("Sum(final_mask) > 0", "No overlap + passing extra part");
    
    auto df_filtered = df_pass
        .Define("FILT_EP_PX", "EXTRA_PARTS_PX[final_mask]")
        .Define("FILT_EP_PY", "EXTRA_PARTS_PY[final_mask]")
        .Define("FILT_EP_PZ", "EXTRA_PARTS_PZ[final_mask]")
        .Define("FILT_EP_MOTHER", "EXTRA_PARTS_MC_MOTHER_ID[final_mask]")
        .Define("FILT_EP_GM", "EXTRA_PARTS_MC_GD_MOTHER_ID[final_mask]")
        .Define("FILT_EP_GM_GM", "EXTRA_PARTS_MC_GD_GD_MOTHER_ID[final_mask]");

    
    auto df_cat = df_filtered.Define("category",
        [](int D_m, int D_gm, int D_ggm,
           ROOT::RVec<float> ep_m, ROOT::RVec<float> ep_gm, ROOT::RVec<float> ep_gmgm) {

            int dm = std::abs(D_m);
            int dgm = std::abs(D_gm);
            int dggm = std::abs(D_ggm);

            auto is_from_same_tree = [&](int m1, int m2, int m3) {
                for (size_t i = 0; i < ep_m.size(); ++i) {
                    int em = std::abs((int)ep_m[i]);
                    int egm = std::abs((int)ep_gm[i]);
                    int eggm = std::abs((int)ep_gmgm[i]);

                    if (em == m1 || em == m2 || em == m3 ||
                        egm == m1 || egm == m2 || egm == m3 ||
                        eggm == m1 || eggm == m2 || eggm == m3) {
                        return true;
                    }
                }
                return false;
            };

            if (!is_from_same_tree(dm, dgm, dggm)) return 3; // combinatorial 

            if (dm == 425 || dgm == 425) return 1; // D_2*0 family
            if (dm == 10423 || dgm == 10423 || dm == 20423 || dgm == 20423) return 2; // D_1 family
            
            bool is_B_m = (dm == 511 || dm == 521);
            bool is_B_gm = (dgm == 511 || dgm == 521 || dggm == 511 || dggm == 521);
            if (is_B_m) return 4;
            if (is_B_gm) return 5;

            return 0;
        },
        {"Bz_Dp_MC_MOTHER_ID", "Bz_Dp_MC_GD_MOTHER_ID", "Bz_Dp_MC_GD_GD_MOTHER_ID",
         "FILT_EP_MOTHER", "FILT_EP_GM", "FILT_EP_GM_GM"});

    float mass_pi = 139.57;
    float mass_D_PDG = 1869.66;
    
    ROOT::RDF::TH1DModel model_m_DEP("h_m_DEP", "EXTRA_PARTS peak; m(Dpi) - m(D) + m(Dpdg) [MeV]; Events", 100, 2100, 2600);
    
    auto df_mdep = df_cat
        .Define("m_DEP", [mass_pi, mass_D_PDG](float Dpx, float Dpy, float Dpz, double D_M,
                                               ROOT::RVec<float> EP_px, ROOT::RVec<float> EP_py, ROOT::RVec<float> EP_pz) {
            
            RVec<float> masses;
            
            PxPyPzMVector pD_m(Dpx, Dpy, Dpz, D_M);
            PxPyPzEVector pD(pD_m);
            
            for (size_t i = 0; i < EP_px.size(); ++i) {
                PxPyPzMVector pEP_m(EP_px[i], EP_py[i], EP_pz[i], mass_pi);
                PxPyPzEVector pEP(pEP_m);
                
                PxPyPzEVector pDEP = pD + pEP;
                
                float improved_mass = pDEP.M() - pD.M() + mass_D_PDG;
                masses.push_back(improved_mass);
            }
            
            std::sort(masses.begin(), masses.end()); // check for unique masses within a vector of extra tracks. not sure if this should be kept.
            
            std::vector<float> unique_masses;
            for (float m : masses) {
                if (unique_masses.empty() || std::abs(m - unique_masses.back()) >= 0.0001) {
                    unique_masses.push_back(m);
                }
            }
            return unique_masses;
            
        }, {"Bz_Dp_PX", "Bz_Dp_PY", "Bz_Dp_PZ", "Bz_Dp_M",
            "FILT_EP_PX", "FILT_EP_PY", "FILT_EP_PZ"});
    
    auto h_cat0 = df_mdep.Filter("category == 0").Histo1D(model_m_DEP, "m_DEP");
    auto h_cat1 = df_mdep.Filter("category == 1").Histo1D(model_m_DEP, "m_DEP");
    auto h_cat2 = df_mdep.Filter("category == 2").Histo1D(model_m_DEP, "m_DEP");
    auto h_cat3 = df_mdep.Filter("category == 3").Histo1D(model_m_DEP, "m_DEP");
    auto h_cat4 = df_mdep.Filter("category == 4").Histo1D(model_m_DEP, "m_DEP");
    auto h_cat5 = df_mdep.Filter("category == 5").Histo1D(model_m_DEP, "m_DEP");
    auto h_cat6 = df_mdep.Filter("category != -1").Histo1D(model_m_DEP, "m_DEP");
    
    df_mdep.Filter("category == 0").Display({"EXTRA_PARTS_MC_MOTHER_ID"}, 50)->Print(); // for debug of clones

    TCanvas* c = new TCanvas("c_3x2", "Categories 3x2", 1500, 1000);
    c->Divide(3, 2);

    auto draw_pad = [](TH1D* h, int color, const char* labelText) {
        h->SetLineColor(color);
        h->SetLineWidth(2);
        h->Draw("HIST");

        double y_min = h->GetMinimum() > 0 ? h->GetMinimum() : 0.1;
        double y_max = h->GetMaximum() * 1.2;

        TLine *line_2460 = new TLine(2460, y_min, 2460, y_max);
        line_2460->SetLineColor(kRed);
        line_2460->SetLineStyle(2);
        line_2460->Draw();

        TLine *line_2280 = new TLine(2280, y_min, 2280, y_max);
        line_2280->SetLineColor(kGreen+2);
        line_2280->SetLineStyle(2);
        line_2280->Draw();

        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(42);
        latex.SetTextSize(0.06);
        latex.DrawLatex(0.15, 0.82, labelText); 
    };

    c->cd(1);
    draw_pad(h_cat0.GetPtr(), kBlack, "Other with same ancestry of EP / Dp");

    c->cd(2);
    draw_pad(h_cat1.GetPtr(), kRed, "D_{2}*0 (2460) mother / gm");

    c->cd(3);
    draw_pad(h_cat2.GetPtr(), kBlue, "D_{1} (2420)");

    c->cd(4);
    draw_pad(h_cat3.GetPtr(), kGreen+2, "Other in which no same ancestry");

    c->cd(5);
    draw_pad(h_cat4.GetPtr(), kYellow+2, "B meson mother");

    c->cd(6);
    draw_pad(h_cat5.GetPtr(), kCyan-1, "B meson grandmother");

    c->SaveAs("m_DEP_WS_FINAL.png");
}
