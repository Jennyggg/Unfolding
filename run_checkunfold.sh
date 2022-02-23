#!/bin/bash
python checkunfold.py --var "nparticle_eta2p4pt05_pur"
python checkunfold.py --var "nparticle_eta2p4pt05_pur" --inputfilesim '/work/jinw/CMSSW_11_1_0/src/EXOVVNtuplizerRunII/Ntuplizer/flatTuple_MB_TuneUp_trk_noPU_new.root' --MCtag 'tuneup'
python checkunfold.py --var "nparticle_eta2p4pt05_pur" --inputfilesim '/work/jinw/CMSSW_11_1_0/src/EXOVVNtuplizerRunII/Ntuplizer/flatTuple_MB_TuneDown_trk_noPU_new.root' --MCtag 'tunedown'
python checkunfold.py --var "spherocity_eta2p4pt05_pur"
python checkunfold.py --var "spherocity_eta2p4pt05_pur" --inputfilesim '/work/jinw/CMSSW_11_1_0/src/EXOVVNtuplizerRunII/Ntuplizer/flatTuple_MB_TuneUp_trk_noPU_new.root' --MCtag 'tuneup'
python checkunfold.py --var "spherocity_eta2p4pt05_pur" --inputfilesim '/work/jinw/CMSSW_11_1_0/src/EXOVVNtuplizerRunII/Ntuplizer/flatTuple_MB_TuneDown_trk_noPU_new.root' --MCtag 'tunedown'
python checkunfold.py --var "thrust_eta2p4pt05_pur"
python checkunfold.py --var "thrust_eta2p4pt05_pur" --inputfilesim '/work/jinw/CMSSW_11_1_0/src/EXOVVNtuplizerRunII/Ntuplizer/flatTuple_MB_TuneUp_trk_noPU_new.root' --MCtag 'tuneup'
python checkunfold.py --var "thrust_eta2p4pt05_pur" --inputfilesim '/work/jinw/CMSSW_11_1_0/src/EXOVVNtuplizerRunII/Ntuplizer/flatTuple_MB_TuneDown_trk_noPU_new.root' --MCtag 'tunedown'
python checkunfold.py --var "broaden_eta2p4pt05_pur"
python checkunfold.py --var "broaden_eta2p4pt05_pur" --inputfilesim '/work/jinw/CMSSW_11_1_0/src/EXOVVNtuplizerRunII/Ntuplizer/flatTuple_MB_TuneUp_trk_noPU_new.root' --MCtag 'tuneup'
python checkunfold.py --var "broaden_eta2p4pt05_pur" --inputfilesim '/work/jinw/CMSSW_11_1_0/src/EXOVVNtuplizerRunII/Ntuplizer/flatTuple_MB_TuneDown_trk_noPU_new.root' --MCtag 'tunedown'
python checkunfold.py --var "transversespherocity_eta2p4pt05_pur"
python checkunfold.py --var "transversespherocity_eta2p4pt05_pur" --inputfilesim '/work/jinw/CMSSW_11_1_0/src/EXOVVNtuplizerRunII/Ntuplizer/flatTuple_MB_TuneUp_trk_noPU_new.root' --MCtag 'tuneup'
python checkunfold.py --var "transversespherocity_eta2p4pt05_pur" --inputfilesim '/work/jinw/CMSSW_11_1_0/src/EXOVVNtuplizerRunII/Ntuplizer/flatTuple_MB_TuneDown_trk_noPU_new.root' --MCtag 'tunedown'




