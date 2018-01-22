#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(
        "file:myMicroAODOutputFile.root"        
        )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test.root")
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )



process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
if os.environ["CMSSW_VERSION"].count("CMSSW_7_6"):
    process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'
elif os.environ["CMSSW_VERSION"].count("CMSSW_7_4"):
    process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v4' 
elif os.environ["CMSSW_VERSION"].count("CMSSW_8_0"):
    process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2'
else:
    raise Exception,"Could not find a sensible CMSSW_VERSION for default globaltag"

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test.root")
)

process.load("flashgg.Taggers.flashggUpdatedIdMVADiPhotons_cfi") 
from flashgg.Taggers.flashggPreselectedDiPhotons_cfi import flashggPreselectedDiPhotons
process.kinPreselDiPhotons = flashggPreselectedDiPhotons.clone(
cut=cms.string(
        "leadingPhoton.pt > 40 && subLeadingPhoton.pt > 30"
        " && abs(leadingPhoton.superCluster.eta)<2.5 && abs(subLeadingPhoton.superCluster.eta)<2.5 "
        " && ( abs(leadingPhoton.superCluster.eta)<1.4442 || abs(leadingPhoton.superCluster.eta)>1.566)"
        " && ( abs(subLeadingPhoton.superCluster.eta)<1.4442 || abs(subLeadingPhoton.superCluster.eta)>1.566)"
        )
                                                              )



process.load("flashgg.Taggers.bRegressionDumper_cfi")
import flashgg.Taggers.dumperConfigTools as cfgTools
from flashgg.Taggers.bRegressionDumper_cfi import bRegressionDumper
bRegressionDumper.dumpTrees = True


#process.bRegressionDumper.src = "kinPreselDiPhotons"

cfgTools.addCategories(bRegressionDumper,
                       [
#        ("Reject", "diPhoton.mass < 50 || diPhoton.mass > 130", -1),
        ("All", "1", 0)
        # ("EBHighR9", "abs(getProbe.superCluster.eta)<1.4442 && getProbe.full5x5_r9>0.94", 0),
        # ("EBLowR9", "abs(getProbe.superCluster.eta)<1.4442 && getProbe.full5x5_r9<=0.94", 0),
        # ("EEHighR9", "abs(getProbe.superCluster.eta)>1.566 && getProbe.full5x5_r9>0.94", 0),
        # ("EELowR9", "abs(getProbe.superCluster.eta)>1.566 && getProbe.full5x5_r9<=0.94", 0)
        ],
                       variables=[ "jetPt                   :=pt"],
                       histograms=[]
                       )    


process.load("flashgg.Taggers.flashggbRegressionProducer_cfi")
from flashgg.Taggers.flashggbRegressionProducer_cfi import flashggbRegressionProducer


from flashgg.Taggers.flashggTags_cff import flashggUnpackedJets

bRegSequence = cms.Sequence(flashggUnpackedJets+flashggbRegressionProducer+bRegressionDumper)

process.p = cms.Path(bRegSequence)