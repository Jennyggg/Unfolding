import numpy as np
import json
from argparse import ArgumentParser
import os,ast
import time
from math import *
import sys
import ROOT as rt
import  tdrstyle
import CMS_lumi
from Plotting_cfg import *
import itertools
import h5py
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from sklearn import metrics
rt.PyConfig.IgnoreCommandLineOptions = True
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptStat(0)
tdrstyle.setTDRStyle()
rt.TH1.SetDefaultSumw2()
rt.TH2.SetDefaultSumw2()
colors = [rt.kRed,rt.kBlue,rt.kOrange,rt.kViolet]

def PlotUnfoldMC(HistUnfold,HistListMC,LegendListMC,Title,isLogY,isRatio,SavePath):
  c = rt.TCanvas("c","c",5,30,W_long,H_long)
#  c.SetRightMargin(0.09)
  pad1 = rt.TPad("pad1", "pad1", 0.0, 0.2 if isRatio else 0.0,1, 1)
  pad1.SetRightMargin(0.03)
  pad1.Draw()
  pad1.cd()
  pad1.SetLogy(isLogY)
  pad1.Update()
  total = 0.
  for ibin in range(HistUnfold.GetNbinsX()):
    total += HistUnfold.GetBinContent(ibin+1)
  for ibin in range(HistUnfold.GetNbinsX()):
    HistUnfold.SetBinError(ibin+1,HistUnfold.GetBinError(ibin+1)/total)
    HistUnfold.SetBinContent(ibin+1,HistUnfold.GetBinContent(ibin+1)/total)
  HistUnfold.SetFillStyle(0)
  HistUnfold.SetLineColor(rt.kBlack)
  HistUnfold.SetLineWidth(2)
  HistUnfold.SetMarkerStyle(20)
  HistUnfold.SetMarkerColor(rt.kBlack) 
  for HistMC in HistListMC:
    total_MC=0.
    for ibin in range(HistMC.GetNbinsX()):
      total_MC += HistMC.GetBinContent(ibin+1)
    for ibin in range(HistMC.GetNbinsX()):
      HistMC.SetBinError(ibin+1,HistMC.GetBinError(ibin+1)/total_MC)
      HistMC.SetBinContent(ibin+1,HistMC.GetBinContent(ibin+1)/total_MC)

  xAxis = HistUnfold.GetXaxis()
  xAxis.SetTitle(Title)
  xAxis.SetTitleSize(FTS)
  xAxis.SetTitleOffset(1.05)
  if isRatio:
    xAxis.SetLabelSize(0)
  else:
    xAxis.SetLabelSize(FLS)

  yAxis = HistUnfold.GetYaxis()
  yAxis.SetNdivisions(6,5,0)
  yAxis.SetLabelSize(FLS)
  yAxis.SetTitleSize(FTS)
  yAxis.SetTitleOffset(FTO)
  yAxis.SetMaxDigits(3)
  yAxis.SetTitle("Events normalized to 1")
  yAxis.SetTitleOffset(1.1)
  ymax = max([HistUnfold.GetBinContent(ibin+1) for ibin in range(HistUnfold.GetNbinsX())])
  ymin = min([HistUnfold.GetBinContent(ibin+1) for ibin in range(HistUnfold.GetNbinsX())])
  HistUnfold.SetAxisRange(ymin*0.1 if isLogY else 0, ymax*50 if isLogY else ymax*1.5,"Y")
  HistUnfold.Draw("pe")
  HistListMC_step = []
  for index,HistMC in enumerate(HistListMC):
    HistMC.SetLineColor(colors[index])
    HistMC.SetLineWidth(2)
#    HistMC.SetFillStyle(3354)
    HistMC.Draw("esame")
    HistListMC_step.append(HistMC.Clone())
    HistListMC_step[index].SetFillStyle(0)
    HistListMC_step[index].SetLineWidth(2)
    HistListMC_step[index].Draw("histsame")
  Legend = rt.TLegend(0.6,0.72,0.95,0.92)
  Legend.SetBorderSize(0)
  Legend.SetFillColor(0)
  Legend.SetFillStyle(0)
  Legend.SetTextFont(42)
  Legend.SetTextSize(0.035)
  Legend.AddEntry(HistUnfold,"Data","lep")
  for index,HistMC in enumerate(HistListMC):
    Legend.AddEntry(HistMC,LegendListMC[index],"lep")
  Legend.Draw()
  pad1.Update()
  box = create_paves(1, "DataPAS", CMSposX=0.17, CMSposY=0.85,
                   prelimPosX=0.17, prelimPosY=0.80,
                   lumiPosX=0.975, lumiPosY=0.91, alignRight=False)

  box["CMS"].Draw("same")
  box["label"].Draw()
  pad1.Update()
  c.Update()
  if isRatio:
    rt.gStyle.SetHatchesLineWidth(1)
    c.cd()
    p1r = rt.TPad("p4","",0,0,1,0.27)
    p1r.SetRightMargin(0.03)
    p1r.SetTopMargin(P2TM)
    p1r.SetBottomMargin(P2BM)
    p1r.SetTicks()
    p1r.Draw()
    p1r.cd()
    xmin = float(HistUnfold.GetXaxis().GetXmin())
    xmax = float(HistUnfold.GetXaxis().GetXmax())
    one = rt.TF1("one","1",xmin,xmax)
    one.SetLineColor(1)
    one.SetLineStyle(2)
    one.SetLineWidth(1)
    HistUnfoldRatio = HistUnfold.Clone()
    HistListMCRatio = []
    for HistMC in HistListMC:
      HistListMCRatio.append(HistMC.Clone())
    for ibin in range(HistUnfold.GetNbinsX()):
      ndata = HistUnfold.GetBinContent(ibin+1)
      edata = HistUnfold.GetBinError(ibin+1)
      HistUnfoldRatio.SetBinContent(ibin+1,1)
      HistUnfoldRatio.SetBinError(ibin+1,edata/ndata if ndata>0 else 0)
      for index,HistMC in enumerate(HistListMC):
        HistListMCRatio[index].SetBinContent(ibin+1,HistMC.GetBinContent(ibin+1)/ndata if ndata>0 else 0)
        HistListMCRatio[index].SetBinError(ibin+1,HistMC.GetBinError(ibin+1)/ndata if ndata>0 else 0)

    HistUnfoldRatio.SetAxisRange(0,2,"Y")
    HistUnfoldRatio.SetTitle("")
    HistUnfoldRatio.GetXaxis().SetTitle(Title)
    HistUnfoldRatio.GetXaxis().SetTitleSize(RTSX)
    HistUnfoldRatio.GetXaxis().SetTitleOffset(RTOX)
    HistUnfoldRatio.GetXaxis().SetLabelSize(RLSX)
    HistUnfoldRatio.GetXaxis().SetLabelOffset(0.02)
    HistUnfoldRatio.GetYaxis().SetTitleSize(RTSY)
    HistUnfoldRatio.GetYaxis().SetLabelSize(RLSY)
    HistUnfoldRatio.GetYaxis().SetTitleOffset(RTOY)
    HistUnfoldRatio.GetYaxis().SetTitle("            MC / Data")
    HistUnfoldRatio.GetYaxis().SetDecimals(1)
    HistUnfoldRatio.GetYaxis().SetMaxDigits(3)
    HistUnfoldRatio.GetYaxis().SetNdivisions(4,2,0)
    HistUnfoldRatio.Draw("pe")
    HistListMCRatio_step = []
    for index,HistMCRatio in enumerate(HistListMCRatio):
      HistMCRatio.Draw("esame")
      HistListMCRatio_step.append(HistMCRatio.Clone())
      HistListMCRatio_step[index].SetFillStyle(0)
      HistListMCRatio_step[index].Draw("histsame")
    one.Draw("same")
    p1r.Update()
  c.SaveAs(SavePath+".pdf")
  c.SaveAs(SavePath+".png")


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('--varunfold', default='VarUnfold.json', help="json file storing the binning of the variables to be unfolded")
    parser.add_argument('--var', default='nparticle', help="The variable to be unfolded")
    parser.add_argument('--outputdir', default='plots', help="directory to store the output plots")
    parser.add_argument('--inputdir', default='results', help="directory of the unfolding results")

    args = parser.parse_args()
    with open(args.varunfold, 'r') as fjson:
        info_var = json.load(fjson)
    if not os.path.exists(args.outputdir):
      os.makedirs(args.outputdir)
    if not os.path.exists(args.inputdir):
      sys.exit("The directory "+args.inputdir+" does not exist!")
    f_norminal = rt.TFile(args.inputdir+"/unfold_"+args.var+"_norminal.root","read")
    f_tuneup = rt.TFile(args.inputdir+"/unfold_"+args.var+"_tuneup.root","read")
    f_tunedown = rt.TFile(args.inputdir+"/unfold_"+args.var+"_tunedown.root","read")
    PlotUnfoldMC(f_norminal.Get("HistGenUnfold"),[f_norminal.Get("HistGenReweight"),f_tuneup.Get("HistGenReweight"),f_tunedown.Get("HistGenReweight")],["Reweighted norminal MC","Reweighted tuneup MC","Reweighted tunedown MC"],info_var[args.var]["gen_name"],1,1,args.outputdir+"/unfold_reweight_gen_"+args.var+"_logy")
    PlotUnfoldMC(f_norminal.Get("HistGenUnfold"),[f_norminal.Get("HistGenReweight"),f_tuneup.Get("HistGenReweight"),f_tunedown.Get("HistGenReweight")],["Reweighted norminal MC","Reweighted tuneup MC","Reweighted tunedown MC"],info_var[args.var]["gen_name"],0,1,args.outputdir+"/unfold_reweight_gen_"+args.var)

    PlotUnfoldMC(f_norminal.Get("HistGenUnfold"),[f_norminal.Get("HistGen"),f_tuneup.Get("HistGen"),f_tunedown.Get("HistGen")],["Norminal MC","Tuneup MC","Tunedown MC"],info_var[args.var]["gen_name"],1,1,args.outputdir+"/unfold_MC_gen_"+args.var+"_logy")
    PlotUnfoldMC(f_norminal.Get("HistGenUnfold"),[f_norminal.Get("HistGen"),f_tuneup.Get("HistGen"),f_tunedown.Get("HistGen")],["Norminal MC","Tuneup MC","Tunedown MC"],info_var[args.var]["gen_name"],0,1,args.outputdir+"/unfold_MC_gen_"+args.var)

    PlotUnfoldMC(f_norminal.Get("HistReco_data"),[f_norminal.Get("HistReco_MC_norm2data"),f_tuneup.Get("HistReco_MC_norm2data"),f_tunedown.Get("HistReco_MC_norm2data")],["Norminal MC","Tuneup MC","Tunedown MC"],info_var[args.var]["reco_name"],1,1,args.outputdir+"/data_MC_reco_"+args.var+"_logy")
    PlotUnfoldMC(f_norminal.Get("HistReco_data"),[f_norminal.Get("HistReco_MC_norm2data"),f_tuneup.Get("HistReco_MC_norm2data"),f_tunedown.Get("HistReco_MC_norm2data")],["Norminal MC","Tuneup MC","Tunedown MC"],info_var[args.var]["reco_name"],0,1,args.outputdir+"/data_MC_reco_"+args.var)



