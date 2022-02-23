import numpy as np
import json
from argparse import ArgumentParser
import os,ast
import time
import math
import sys,stat
import ROOT
import itertools
import h5py
import pandas as pd
from math import sqrt
def ReweightGenHist(Tree,Tree_data,Cut,RecoObs, GenVar, Weight, NbinsReco, MinReco, MaxReco, NbinsGen, MinGen, MaxGen):
  hist2d = ROOT.TH2F("Hist2D","Hist2D",NbinsReco, MinReco, MaxReco, NbinsGen, MinGen, MaxGen)
  histreco_MC = ROOT.TH1F("HistReco_MC","HistReco_MC",NbinsReco, MinReco, MaxReco)
  histreco_data = ROOT.TH1F("HistReco_data","HistReco_data",NbinsReco, MinReco, MaxReco)
  histreco_MC_total = ROOT.TH1F("HistReco_MC_total","HistReco_MC_total",1, MinReco, MaxReco)
  histreco_data_total = ROOT.TH1F("HistReco_data_total","HistReco_data_total",1, MinReco, MaxReco)
  Tree.Draw(GenVar+":"+RecoObs+">>Hist2D",Weight+"*"+Cut,"colzgoff")
  Tree.Draw(RecoObs+">>HistReco_MC",Weight+"*"+Cut,"goff")
  Tree.Draw(RecoObs+">>HistReco_MC_total",Weight+"*"+Cut,"goff")
  if Tree_data.GetBranchStatus(Weight):
    Tree_data.Draw(RecoObs+">>HistReco_data",Weight+"*"+Cut,"goff")
    Tree_data.Draw(RecoObs+">>HistReco_data_total",Weight+"*"+Cut,"goff")
  else:
    Tree_data.Draw(RecoObs+">>HistReco_data",Cut,"goff")
    Tree_data.Draw(RecoObs+">>HistReco_data_total",Cut,"goff")
  ratio = histreco_data_total.GetBinContent(1)/histreco_MC_total.GetBinContent(1)
  for ibinreco in range(NbinsReco):
    if histreco_MC.GetBinContent(ibinreco+1)>0:
      reco_data=histreco_data.GetBinContent(ibinreco+1)
      reco_MC=histreco_MC.GetBinContent(ibinreco+1)
      reco_data_e=histreco_data.GetBinError(ibinreco+1)
      reco_MC_e=histreco_MC.GetBinError(ibinreco+1)
      data_to_MC = reco_data/reco_MC
      for ibingen in range(NbinsGen):
        if hist2d.GetBinContent(ibinreco+1,ibingen+1)>0:
          hist2d.SetBinError(ibinreco+1,ibingen+1,hist2d.GetBinContent(ibinreco+1,ibingen+1)*data_to_MC*sqrt(reco_data_e**2/reco_data**2+reco_MC_e**2/reco_MC**2+hist2d.GetBinError(ibinreco+1,ibingen+1)**2/hist2d.GetBinContent(ibinreco+1,ibingen+1)**2))
          hist2d.SetBinContent(ibinreco+1,ibingen+1,hist2d.GetBinContent(ibinreco+1,ibingen+1)*data_to_MC)
  histgen_MC_reweight = hist2d.ProjectionY("HistGenReweight",1,NbinsReco)
  return histgen_MC_reweight,ratio,histreco_MC,histreco_data

def UnfoldGenHist(FitResultFile,Tree,Cut,Weight,GenVar,NbinsGen, MinGen, MaxGen):
  histgenunfold = ROOT.TH1F("HistGenUnfold","HistGenUnfold",NbinsGen, MinGen, MaxGen)
  histgen = ROOT.TH1F("HistGen","HistGen",NbinsGen, MinGen, MaxGen)
  Tree.Draw(GenVar+">>HistGen",Weight+"*"+Cut,"goff")
  Tree.Draw(GenVar+">>HistGenUnfold",Weight+"*"+Cut,"goff")
  f_fit = ROOT.TFile.Open(FitResultFile,'read')
  fit_mdf = f_fit.Get("fit_mdf")
  for ibingen in range(NbinsGen):
    yield_genbin_var = fit_mdf.floatParsFinal().find("yield_genbin"+str(ibingen+1))
    yield_genbin = yield_genbin_var.getValV()
    yield_genbin_e = yield_genbin_var.getError()
    if histgenunfold.GetBinContent(ibingen+1)>0 and yield_genbin>0 :
      histgenunfold.SetBinError(ibingen+1,histgenunfold.GetBinContent(ibingen+1)*sqrt(histgenunfold.GetBinError(ibingen+1)**2/histgenunfold.GetBinContent(ibingen+1)**2+yield_genbin_e**2/yield_genbin**2))
      histgenunfold.SetBinContent(ibingen+1,histgenunfold.GetBinContent(ibingen+1)*yield_genbin)
    else:
      histgenunfold.SetBinError(ibingen+1,0)
      histgenunfold.SetBinContent(ibingen+1,0)
  return histgenunfold,histgen

if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('--varunfold', default='VarUnfold.json', help="json file storing the binning of the variables to be unfolded")
    parser.add_argument('--var', default='nparticle', help="The variable to be unfolded")
    parser.add_argument('--outputdir', default='results', help="directory to store the unfolded histograms and plots")
    parser.add_argument('--inputfilesim', default='/work/jinw/CMSSW_11_1_0/src/EXOVVNtuplizerRunII/Ntuplizer/flatTuple_MB_trk_noPU_new.root', help="The input root file containing the MC events for unfolding")
    parser.add_argument('--inputfiledata', default='/work/jinw/CMSSW_11_1_0/src/EXOVVNtuplizerRunII/Ntuplizer/flatTuple_ZeroBias_2018lowPU_new.root', help="The input root file containing the data events for unfolding")
    parser.add_argument('--MCtag', default='norminal', help="The tag for MC")
    parser.add_argument('--cut', default='(PV_N_good==1&&PV_isgood&&(Instanton_N_Trk_highPurity_pt05>2))', help="The event selection")
    parser.add_argument('--fitdir',default = "rootfiles",help="The directory storing the fit results")
    parser.add_argument('--nonnormalize',action="store_true",default=False,help="unset the normalization")
    args = parser.parse_args()

    with open(args.varunfold, 'r') as fjson:
        info_var = json.load(fjson)
    fin = ROOT.TFile(args.inputfilesim,"READ")
    tree = fin.Get("ntuplizer/tree")
    fin_data = ROOT.TFile(args.inputfiledata,"READ")
    tree_data = fin_data.Get("ntuplizer/tree")
    if not os.path.exists(args.outputdir):
      os.makedirs(args.outputdir)
    if not os.path.exists(args.fitdir):
      sys.exit("The directory "+args.fitdir+" does not exist!")

    hist_gen_reweight,data_to_MC_ratio,hist_reco_MC,hist_reco_data = ReweightGenHist(tree,tree_data,args.cut,info_var[args.var]["reco"],info_var[args.var]["gen"], "genWeight", info_var[args.var]["nbinsreco"],info_var[args.var]["minreco"],info_var[args.var]["maxreco"],info_var[args.var]["nbinsgen"],info_var[args.var]["mingen"],info_var[args.var]["maxgen"])
    hist_gen_unfold,hist_gen = UnfoldGenHist(args.fitdir+"/multidimfit"+args.var+"_"+args.MCtag+"_data.root",tree,args.cut,"genWeight",info_var[args.var]["gen"],info_var[args.var]["nbinsgen"],info_var[args.var]["mingen"],info_var[args.var]["maxgen"])
    fout = ROOT.TFile(args.outputdir+"/unfold_"+args.var+"_"+args.MCtag+".root","recreate")
    hist_gen_reweight.Write()
    if args.nonnormalize:
      hist_gen_unfold.Write()
    else:
      for ibingen in range(info_var[args.var]["nbinsgen"]):
        hist_gen_unfold.SetBinContent(ibingen+1,hist_gen_unfold.GetBinContent(ibingen+1)*data_to_MC_ratio)
        hist_gen_unfold.SetBinError(ibingen+1,hist_gen_unfold.GetBinError(ibingen+1)*data_to_MC_ratio)
        hist_gen.SetBinContent(ibingen+1,hist_gen.GetBinContent(ibingen+1)*data_to_MC_ratio)
        hist_gen.SetBinError(ibingen+1,hist_gen.GetBinError(ibingen+1)*data_to_MC_ratio)
      hist_gen_unfold.Write()
      hist_gen.Write()
    hist_reco_MC_norm2data = ROOT.TH1F("HistReco_MC_norm2data","HistReco_MC_norm2data",info_var[args.var]["nbinsreco"],info_var[args.var]["minreco"],info_var[args.var]["maxreco"])
    for ibinreco in range(info_var[args.var]["nbinsreco"]):
      hist_reco_MC_norm2data.SetBinContent(ibinreco+1,hist_reco_MC.GetBinContent(ibinreco+1)*data_to_MC_ratio)
      hist_reco_MC_norm2data.SetBinError(ibinreco+1,hist_reco_MC.GetBinError(ibinreco+1)*data_to_MC_ratio)
    hist_reco_MC.Write()
    hist_reco_data.Write()
    hist_reco_MC_norm2data.Write()
    fout.Close()
