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
import CombineHarvester.CombineTools.ch as ch


def WriteHist(Tree,Tree_data, FileName, Cut, RecoObs, GenVar, Weight, NbinsReco, MinReco, MaxReco, NbinsGen, MinGen, MaxGen,nonnormalize):
  hist2d = ROOT.TH2F("Hist2D","Hist2D",NbinsReco, MinReco, MaxReco, NbinsGen, MinGen, MaxGen)
  histreco = ROOT.TH1F("HistReco","HistReco",NbinsReco, MinReco, MaxReco)
  hist2d_total = ROOT.TH2F("Hist2D_total","Hist2D_total",1, MinReco, MaxReco, 1, MinGen, MaxGen)
  histreco_total = ROOT.TH1F("HistReco_total","HistReco_total",1, MinReco, MaxReco)
  Tree.Draw(GenVar+":"+RecoObs+">>Hist2D",Weight+"*"+Cut,"colzgoff")
  Tree.Draw(GenVar+":"+RecoObs+">>Hist2D_total",Weight+"*"+Cut,"colzgoff")
  if Tree_data.GetBranchStatus(Weight):
    Tree_data.Draw(RecoObs+">>HistReco",Weight+"*"+Cut,"goff")
    Tree_data.Draw(RecoObs+">>HistReco_total",Weight+"*"+Cut,"goff")
  else:
    Tree_data.Draw(RecoObs+">>HistReco",Cut,"goff")
    Tree_data.Draw(RecoObs+">>HistReco_total",Cut,"goff")
  ratio = histreco_total.GetBinContent(1)/hist2d_total.GetBinContent(1,1)
  print("ratio",ratio)
  f = ROOT.TFile(FileName,"RECREATE")
  N_data=0
  N_MC=0
  for ibinreco in range(NbinsReco):
    hist_obs = ROOT.TH1F("recobin"+str(ibinreco+1)+"_data_obs","observed events in recobin"+str(ibinreco+1),1,0,1)
    hist_obs.SetBinContent(1,histreco.GetBinContent(ibinreco+1))
    N_data+=hist_obs.GetBinContent(1)
    hist_obs.SetBinError(1,histreco.GetBinError(ibinreco+1))
    hist_obs.Write()
    for ibingen in range(NbinsGen):
      hist = ROOT.TH1F("recobin"+str(ibinreco+1)+"_"+"genbin"+str(ibingen+1), "events in recobin"+str(ibinreco+1)+" genbin"+str(ibingen+1), 1,0,1)
      if nonnormalize:
        hist.SetBinContent(1,hist2d.GetBinContent(ibinreco+1,ibingen+1))
        hist.SetBinError(1,hist2d.GetBinError(ibinreco+1,ibingen+1))
      else:
        hist.SetBinContent(1,hist2d.GetBinContent(ibinreco+1,ibingen+1)*ratio)
        hist.SetBinError(1,hist2d.GetBinError(ibinreco+1,ibingen+1)*ratio)
      N_MC+=hist.GetBinContent(1)
      hist.Write()
  print("data",N_data,"MC",N_MC)
  print("Writing histograms to root file: ",FileName)
  f.Close()

def WriteDatacard(analysis, era, NbinsReco, NbinsGen, rootfile, output_shape, output_datacard):
  cb = ch.CombineHarvester()
  cat_reco = [(i+1,"recobin"+str(i+1)) for i in range(NbinsReco) ]
  cb.AddObservations(['*'], [analysis], [era], ["reco"], cat_reco)
  cb.AddProcesses(['*'], [analysis], [era], ["reco"],["genbin"+str(i+1) for i in range(NbinsGen)], cat_reco, True)
#  cb.cp().signals().AddSyst(cb,"lumi","lnN",ch.SystMap()(1.10))
  print("Extracting shapes from ",rootfile)
  cb.cp().signals().ExtractShapes(
     rootfile,
     "$BIN_$PROCESS",
     "$BIN_$PROCESS_$SYSTEMATIC"
  )
  extrastring=""
  for i in range(NbinsGen):
    extrastring += "yield_genbin"+str(i+1)+" rateParam * genbin"+str(i+1)+" 1.000000 [0,10.000000]\n"
  cb.AddDatacardLineAtEnd(extrastring)
  cb.AddDatacardLineAtEnd('* autoMCStats 0')
  cb.WriteDatacard(output_datacard,ROOT.TFile.Open(output_shape,'recreate'))
  print("Writing datacard: ",output_datacard)
  print("Writing shapes: ",output_shape)

def WriteScript(datacard_name,output_script,NbinsGen,var,MCtag):
  GenParaInit = "--PO verbose "
  GenParaFit = ""
  GenParaFit0 = ""
  for i in range(NbinsGen):
    GenParaInit += "--PO \'map=.*/genbin{0}:yield_genbin{0}[1,0,10]\' ".format(i+1)
    if i != NbinsGen-1:
      GenParaFit += "yield_genbin{}=5,".format(i+1)
      GenParaFit0 += "yield_genbin{}=0,".format(i+1)
    else:
      GenParaFit += "yield_genbin{}=5".format(i+1)
      GenParaFit0 += "yield_genbin{}=0".format(i+1)
  InitWorkSpace = """
#!/bin/bash
text2workspace.py {0} -o workspace.root --X-allow-no-background -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel {1}
""".format(datacard_name,GenParaInit)
#  FitScript = """
#combine -M MultiDimFit -d workspace.root --setParameters={} -t -1
#""".format(GenParaFit)
  FitScript = """
combine -M MultiDimFit -d workspace.root -t -1 --setParameters={0} -n {2}_{3}_asimovs --robustFit 1 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=999999999  --cminDefaultMinimizerStrategy 0 --verbose 3 --algo singles --saveFitResult | tee MultiDimFit_{2}_{3}_asimovs.log
combine -M MultiDimFit -d workspace.root -t -1 --setParameters={1} -n {2}_{3}_asimovb --robustFit 1 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=999999999  --cminDefaultMinimizerStrategy 0 --verbose 3 --algo singles --saveFitResult | tee MultiDimFit_{2}_{3}_asimovb.log
combine -M MultiDimFit -d workspace.root  -n {2}_{3}_data --robustFit 1 --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=999999999  --cminDefaultMinimizerStrategy 0 --verbose 3 --algo singles --saveFitResult | tee MultiDimFit_{2}_{3}_data.log
#combine -M FitDiagnostics -d workspace.root -t -1 --setParameters={0} --saveShapes --saveNormalizations --saveWithUncertainties --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=999999999 --robustFit=1  -n {2}_asimovs --saveOverallShapes --cminDefaultMinimizerStrategy 0 2
#combine -M FitDiagnostics -d workspace.root -t -1 --setParameters={1} --saveShapes --saveNormalizations --saveWithUncertainties --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=999999999 --robustFit=1  -n {2}_asimovb --saveOverallShapes --cminDefaultMinimizerStrategy 0 2
#combine -M FitDiagnostics -d workspace.root --saveShapes --saveNormalizations --saveWithUncertainties --X-rtd MINIMIZER_analytic --X-rtd MINIMIZER_MaxCalls=999999999 --robustFit=1  -n {2}_data --saveOverallShapes --cminDefaultMinimizerStrategy 0 2
""".format(GenParaFit,GenParaFit0,var,MCtag)
  with open(output_script, 'w') as f:
    f.write(InitWorkSpace)
    f.write(FitScript)
  st = os.stat(output_script)
  os.chmod(output_script, st.st_mode | stat.S_IEXEC)
  print("Writing to script: ",output_script)


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('--varunfold', default='VarUnfold.json', help="json file storing the binning of the variables to be unfolded")
    parser.add_argument('--var', default='nparticle', help="The variable to be unfolded")
    parser.add_argument('--outputdir', default='rootfiles', help="directory to store the root files")
    parser.add_argument('--inputfilesim', default='/work/jinw/CMSSW_11_1_0/src/EXOVVNtuplizerRunII/Ntuplizer/flatTuple_MB_trk_noPU_new.root', help="The input root file containing the MC events for unfolding")
    parser.add_argument('--inputfiledata', default='/work/jinw/CMSSW_11_1_0/src/EXOVVNtuplizerRunII/Ntuplizer/flatTuple_ZeroBias_2018lowPU_new.root', help="The input root file containing the data events for unfolding")
    parser.add_argument('--MCtag', default='norminal', help="The tag for MC")
    parser.add_argument('--analysis', default='LowPU2018', help="The dataset name")
    parser.add_argument('--era', default='13TeV', help="The era name")
    parser.add_argument('--nonnormalize',action="store_true",default=False,help="unset the normalization")
    parser.add_argument('--cut', default='(PV_N_good==1&&PV_isgood&&(Instanton_N_Trk_highPurity_pt05>2))', help="The event selection")
    args = parser.parse_args()
    with open(args.varunfold, 'r') as fjson:
        info_var = json.load(fjson)

    fin = ROOT.TFile(args.inputfilesim,"READ")
    tree = fin.Get("ntuplizer/tree")
    fin_data = ROOT.TFile(args.inputfiledata,"READ")
    tree_data = fin_data.Get("ntuplizer/tree")
    if not os.path.exists(args.outputdir):
      os.makedirs(args.outputdir)
    WriteHist(tree,tree_data,args.outputdir+"/"+args.analysis+"_"+args.era+"_"+args.MCtag+"_"+args.var+".root",args.cut,info_var[args.var]["reco"],info_var[args.var]["gen"],"genWeight",info_var[args.var]["nbinsreco"],info_var[args.var]["minreco"],info_var[args.var]["maxreco"],info_var[args.var]["nbinsgen"],info_var[args.var]["mingen"],info_var[args.var]["maxgen"],args.nonnormalize)
    WriteDatacard(args.analysis, args.era,info_var[args.var]["nbinsreco"], info_var[args.var]["nbinsgen"], args.outputdir+"/"+args.analysis+"_"+args.era+"_"+args.MCtag+"_"+args.var+".root", args.outputdir+"/shape_"+args.analysis+"_"+args.era+"_"+args.MCtag+"_"+args.var+".root", args.outputdir+"/datacard_"+args.analysis+"_"+args.era+"_"+args.MCtag+"_"+args.var+".txt")
    WriteScript("datacard_"+args.analysis+"_"+args.era+"_"+args.MCtag+"_"+args.var+".txt",args.outputdir+"/do_fit_"+args.analysis+"_"+args.era+"_"+args.MCtag+"_"+args.var+".sh",info_var[args.var]["nbinsgen"],args.var,args.MCtag)

