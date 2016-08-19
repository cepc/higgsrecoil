#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
using namespace std;

draw_mrecoil()
{
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  //gStyle->SetTitleColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetTitleFillColor(0);

  // set the paper & margin sizes
  gStyle->SetPaperSize(20,26);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.17);
  gStyle->SetPadLeftMargin(0.17);

  // use large Times-Roman fonts
  gStyle->SetTitleFont(22,"xyz");  // set the all 3 axes title font
  gStyle->SetTitleFont(22," ");    // set the pad title font
  gStyle->SetTitleSize(0.06,"xyz"); // set the 3 axes title size
  gStyle->SetLabelFont(22,"xyz");
  gStyle->SetLabelSize(0.06,"xyz");
  gStyle->SetTextFont(22);
  gStyle->SetTextSize(0.08);
  gStyle->SetStatFont(22);

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetLineWidth(2);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  //..Get rid of X error bars
  gStyle->SetErrorX(0.001);

    // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  
  gROOT->Reset();

  TH1D *hsig = new TH1D("hsig","hsig",80,80,160);
  TH1D *hbkg = new TH1D("hbkg","hbkg",80,80,160);

  TChain *fChain1 = new TChain("MCPart");
  fChain1->Add("result/signal.root");
  TChain *fChain2 = new TChain("MCPart");
  fChain2->Add("result/bkg.root");

  TCanvas *myCanvas = new TCanvas();
  myCanvas->Divide(1,1);
  myCanvas->SetFrameLineWidth(3);
  myCanvas->SetFillColor(10);
  myCanvas->cd(1);
  gPad->SetTopMargin(0.015);
  gPad->SetRightMargin(0.08);
  gPad->SetLeftMargin(0.18);
  gPad->SetBottomMargin(0.16);
  gPad->SetFrameLineWidth(2);
  gPad->SetFillColor(0);
  gPad->SetFrameFillColor(0);
//  gPad->SetLogy();

  fChain1->Draw("RecoilMass[4]>>hsig","");
  fChain2->Draw("RecoilMass[4]>>hbkg", "");

  hsig->GetXaxis()->SetTitle("M_{recoil}^{#mu^{+}#mu^{-}}");
  hsig->GetYaxis()->SetTitle("Events");
  hsig->GetXaxis()->SetLabelFont(22);
  hsig->GetYaxis()->SetLabelFont(22);
  hsig->GetXaxis()->CenterTitle();
  hsig->GetYaxis()->CenterTitle();
  hsig->GetXaxis()->SetNdivisions(404);
  hsig->GetYaxis()->SetNdivisions(404);
  
  hsig->GetYaxis()->SetTitleOffset(1.10);
  hsig->GetYaxis()->SetLabelSize(0.05);
  hsig->GetYaxis()->SetLabelOffset(0.01);
  hsig->GetYaxis()->SetTitleSize(0.06);
  hsig->GetYaxis()->SetTitleFont(22);
  hsig->GetXaxis()->SetTitleOffset(1.10);
  hsig->GetXaxis()->SetLabelSize(0.05);
  hsig->GetXaxis()->SetLabelOffset(0.015);
  hsig->GetXaxis()->SetTitleSize(0.06);
  hsig->GetXaxis()->SetTitleFont(22);

  hsig->Draw("");
  hsig->SetMarkerStyle(4);
  hsig->SetMarkerSize(0.8);
  hsig->SetMarkerColor(2);
  hsig->SetLineColor(2);
  hbkg->Draw("same");
  hbkg->SetLineWidth(2);
  hbkg->SetLineColor(kBlue); 

  TLatex * CEPC = new TLatex(0.89,0.90, "CEPC");
  CEPC->SetNDC();
  CEPC->SetTextFont(22);
  CEPC->SetTextSize(0.08);
  CEPC->SetTextAlign(33);
  CEPC->Draw();
  
  TLatex * prelim = new TLatex(0.92,0.82, "Preliminary");
  prelim->SetNDC();
  prelim->SetTextFont(22);
  prelim->SetTextSize(0.055);
  prelim->SetTextAlign(33);
  prelim->Draw();
 
  lg =new TLegend(0.2, 0.75, 0.5, 0.9);//   (0.75,0.75,0.965555,0.95);
  lg->SetBorderSize(0);
  lg->SetFillColor(10);
  lg->AddEntry(hsig," Signal",  "L");
  lg->AddEntry(hbkg," Background",  "L");
  lg->SetTextFont(22);
  lg->SetTextSize(0.045);
  lg->Draw("same");

  myCanvas->SaveAs("plot/signal_bak.pdf"); 

}
                                   
