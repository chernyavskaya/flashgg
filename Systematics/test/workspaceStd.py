#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing
from flashgg.Systematics.SystematicDumperDefaultVariables import minimalVariables,minimalHistograms,minimalNonSignalVariables,systematicVariables
from flashgg.Systematics.SystematicDumperDefaultVariables import minimalVariablesHTXS,systematicVariablesHTXS
import os
import flashgg.Systematics.settings as settings

# SYSTEMATICS SECTION
dropVBFInNonGold = False  # for 2015 only!

process = cms.Process("FLASHggSyst")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("Configuration.StandardSequences.GeometryDB_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
if os.environ["CMSSW_VERSION"].count("CMSSW_7_6"):
    process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'
elif os.environ["CMSSW_VERSION"].count("CMSSW_7_4"):
    process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v4' 
elif os.environ["CMSSW_VERSION"].count("CMSSW_8_0"):
    process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2'
elif os.environ["CMSSW_VERSION"].count("CMSSW_9_4"):
    #process.GlobalTag.globaltag = '94X_mc2017_realistic_v12'
    process.GlobalTag.globaltag = '94X_mc2017_realistic_v14'
else:
    raise Exception,"Could not find a sensible CMSSW_VERSION for default globaltag"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )


systlabels = [""]
phosystlabels = []
metsystlabels = []
jetsystlabels = []
elesystlabels = []
musystlabels = []

from flashgg.MetaData.JobConfig import customize
customize.options.register('tthTagsOnly',
                           False,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'tthTagsOnly'
                           )
customize.options.register('doubleHTagsOnly',
                           True,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'doubleHTagsOnly'
                           )
customize.options.register('doubleHReweight',
                           -1,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.int,
                           'doubleHReweight'
                           )
customize.options.register('year',
                           '2016',
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.string,
                           'year'
                           )
customize.options.register('doDoubleHTag',
                           True,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'doDoubleHTag'
                           )
customize.options.register('doDoubleHttHKiller',
                           True,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'doDoubleHttHKiller'
                           )
customize.options.register('ttHKillerInputVariables',
                           False,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'ttHKillerInputVariables'
                           )
customize.options.register('doDoubleHGenAnalysis',
                           False,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'doDoubleHGenAnalysis'
                           )
customize.options.register('doBJetRegression',
                           True,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'doBJetRegression'
                           )
customize.options.register('doHTXS',
                           False,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'doHTXS'
                           )
customize.options.register('doMuFilter',
                           True,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'doMuFilter'
                           )
customize.options.register('doFiducial',
                           False,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'doFiducial'
                           )
customize.options.register('acceptance',
                           'NONE',
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.string,
                           'acceptance'
                           )
customize.options.register('doSystematics',
                           False,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'doSystematics'
                           )
customize.options.register('doPdfWeights',
                           False,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'doPdfWeights'
                           )
customize.options.register('dumpTrees',
                           False,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'dumpTrees'
                           )
customize.options.register('dumpWorkspace',
                           True,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'dumpWorkspace'
                           )
customize.options.register('verboseTagDump',
                           False,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'verboseTagDump'
                           )
customize.options.register('verboseSystDump',
                           False,
                           VarParsing.VarParsing.multiplicity.singleton,
                           VarParsing.VarParsing.varType.bool,
                           'verboseSystDump'
                           )



print "Printing defaults"
print 'doFiducial '+str(customize.doFiducial)
print 'acceptance '+str(customize.acceptance)
print 'tthTagsOnly '+str(customize.tthTagsOnly)
# import flashgg customization to check if we have signal or background
from flashgg.MetaData.JobConfig import customize
# set default options if needed
customize.setDefault("maxEvents",10000)
customize.setDefault("targetLumi",1.00e+3)
customize.parse()

process_type = 'Sim'
if customize.processId == 'Data' : process_type = customize.processId
settings.init(customize.year,process_type)
year = settings.year
process_type = settings.process_type
if year == "2016" and process_type != "Data":
    process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'
elif year == "2016" and process_type == "Data":
    process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7'
elif year == "2017" and process_type != "Data":  
    process.GlobalTag.globaltag = '94X_mc2017_realistic_v14'
elif year == "2017" and process_type == "Data":  
    process.GlobalTag.globaltag = '94X_dataRun2_ReReco_EOY17_v6'

MUON_ID = "Medium" #["Tight", "Medium" , "Loose", "Soft", "HighPt", "MediumPrompt", "TrkHighPt"]
MUON_ISO = "LooseRel" #{ LooseID : ["LooseRel"],MediumID:["LooseRel", "TightRel"] , TrkHighPtID:["LooseRelTk", "TightRelTk"], TightIDandIPCut:["LooseRel", "TightRel"], HighPtIDandIPCut:["LooseRelTk", "TightRelTk"] }

from flashgg.Systematics.SystematicsCustomize import *
jetSystematicsInputTags = createStandardSystematicsProducers(process , MUON_ID=MUON_ID , MUON_ISO=MUON_ISO)
if dropVBFInNonGold:
    process.flashggVBFTag.SetArbitraryNonGoldMC = True
    process.flashggVBFTag.DropNonGoldData = True
modifyTagSequenceForSystematics(process,jetSystematicsInputTags)

print "Printing options"
print 'doFiducial '+str(customize.doFiducial)
print 'acceptance '+str(customize.acceptance)
print 'tthTagsOnly '+str(customize.tthTagsOnly)
print 'doMuFilter '+str(customize.doMuFilter)

if customize.doFiducial:
    import flashgg.Systematics.fiducialCrossSectionsCustomize as fc
    fc.leadCut = 1./3.
    fc.subLeadCut = 1./4.
    fc.isoCut = 10.
    fc.etaCut = 2.5
    matchCut = "leadingPhoton.hasMatchedGenPhoton() && subLeadingPhoton.hasMatchedGenPhoton()"
    phoIDcut = '(leadingView().phoIdMvaWrtChosenVtx() >0.320 && subLeadingView().phoIdMvaWrtChosenVtx() >0.320)'
    accCut   = fc.getAccRecoCut()
    
    print process.flashggPreselectedDiPhotons.cut

    if customize.acceptance == 'IN':
        process.flashggPreselectedDiPhotons.cut = cms.string(str(process.flashggPreselectedDiPhotons.cut)[12:-2] +' && '+ str(matchCut)+ ' && '+ str(phoIDcut) +' && ' + str(accCut))

    if customize.acceptance == 'OUT':
        process.flashggPreselectedDiPhotons.cut = cms.string(str(process.flashggPreselectedDiPhotons.cut)[12:-2] +' && '+ str(matchCut)+ ' && '+ str(phoIDcut) +' && !' + str(accCut))
        
    if customize.acceptance == 'NONE':
        process.flashggPreselectedDiPhotons.cut = cms.string(str(process.flashggPreselectedDiPhotons.cut)[12:-2] +' && '+ str(phoIDcut))
    print "Here we print the preslection cut"
    print process.flashggPreselectedDiPhotons.cut


process.load("flashgg/Taggers/flashggTagSequence_cfi")
process.load("flashgg/MicroAOD/flashggDiPhotons_cfi")

# needed for 0th vertex from microAOD
if customize.tthTagsOnly:
    process.flashggDiPhotons.whichVertex = cms.uint32(0)
    process.flashggDiPhotons.useZerothVertexFromMicro = cms.bool(True)

process.load("flashgg/Taggers/flashggTagSequence_cfi")
print 'here we print the tag sequence before'
print process.flashggTagSequence
if customize.doFiducial:
    from PhysicsTools.PatAlgos.tools.helpers import cloneProcessingSnippet,massSearchReplaceAnyInputTag
    process.flashggTagSequence.remove(process.flashggVBFTag)
    process.flashggTagSequence.remove(process.flashggTTHDiLeptonTag)
    process.flashggTagSequence.remove(process.flashggTTHLeptonicTag)
    process.flashggTagSequence.remove(process.flashggTTHHadronicTag)
   # process.flashggTagSequence.remove(process.flashggTTHeptonTag)
     #haven't tested VH tags with fiducial cross-section measurement yet
    process.flashggTagSequence.remove(process.flashggVHEtTag)
    process.flashggTagSequence.remove(process.flashggVHLooseTag)
    process.flashggTagSequence.remove(process.flashggVHTightTag)
    process.flashggTagSequence.remove(process.flashggVHMetTag)
    process.flashggTagSequence.remove(process.flashggWHLeptonicTag)
    process.flashggTagSequence.remove(process.flashggZHLeptonicTag)
    process.flashggTagSequence.remove(process.flashggVHLeptonicLooseTag)
    process.flashggTagSequence.remove(process.flashggVHHadronicTag)
    process.flashggTagSequence.replace(process.flashggUntagged, process.flashggSigmaMoMpToMTag)

if customize.tthTagsOnly:
    process.flashggTagSequence.remove(process.flashggVBFTag)
    process.flashggTagSequence.remove(process.flashggVHMetTag)
    process.flashggTagSequence.remove(process.flashggWHLeptonicTag)
    process.flashggTagSequence.remove(process.flashggZHLeptonicTag)
    process.flashggTagSequence.remove(process.flashggVHLeptonicLooseTag)
    process.flashggTagSequence.remove(process.flashggVHHadronicTag)
    process.flashggTagSequence.remove(process.flashggUntagged)
    process.flashggTagSequence.remove(process.flashggVBFMVA)
    process.flashggTagSequence.remove(process.flashggVBFDiPhoDiJetMVA)

if customize.doDoubleHTag:
    import flashgg.Systematics.doubleHCustomize as hhc
    hhc.customizeTagSequence(  customize, process )
    minimalVariables += hhc.variablesToDump(customize)
#####################################Only if needed, tmp fix for Francesco#########################
    if customize.processId == "Data":
         minimalNonSignalVariables += hhc.variablesToDumpData(customize)
###################################################################################################
    print "saving variables:", minimalVariables

print 'here we print the tag sequence after'
print process.flashggTagSequence

if customize.doFiducial:
    print 'we do fiducial and we change tagsorter'
    process.flashggTagSorter.TagPriorityRanges = cms.VPSet(     cms.PSet(TagName = cms.InputTag('flashggSigmaMoMpToMTag')) )

if customize.doubleHTagsOnly:
    process.flashggTagSorter.TagPriorityRanges = cms.VPSet(     cms.PSet(TagName = cms.InputTag('flashggDoubleHTag')) )

if customize.tthTagsOnly:
    process.flashggTagSorter.TagPriorityRanges = cms.VPSet(   cms.PSet(TagName = cms.InputTag('flashggTTHDiLeptonTag')),
        cms.PSet(TagName = cms.InputTag('flashggTTHLeptonicTag')),
        cms.PSet(TagName = cms.InputTag('flashggTTHHadronicTag')) )

    print "customize.processId:",customize.processId

    print "Removing FracRV from syst and adding  PixelSeed"
    
    newvpset = cms.VPSet()
    for pset in process.flashggDiPhotonSystematics.SystMethods:
        if not pset.Label.value().count("FracRVWeight")and not pset.Label.value().count("FracRVNvtxWeight") :
            print  pset.Label.value()
            newvpset += [pset]
    from flashgg.Systematics.flashggDiPhotonSystematics_cfi import PixelSeedWeight
    newvpset += [ PixelSeedWeight ]
    
    process.flashggDiPhotonSystematics.SystMethods = newvpset

print "customize.processId:",customize.processId
# load appropriate scale and smearing bins here
# systematics customization scripts will take care of adjusting flashggDiPhotonSystematics
#process.load("flashgg.Systematics.escales.escale76X_16DecRereco_2015")

# Or use the official tool instead
useEGMTools(process)

# Only run systematics for signal events
# convention: ggh vbf wzh (wh zh) tth
signal_processes = ["ggh_","vbf_","wzh_","wh_","zh_","bbh_","thq_","thw_","tth_","HHTo2B2G","Acceptance"]
is_signal = reduce(lambda y,z: y or z, map(lambda x: customize.processId.count(x), signal_processes))
#if customize.processId.count("h_") or customize.processId.count("vbf_") or customize.processId.count("Acceptance") or customize.processId.count("hh_"): 
if is_signal:
    print "Signal MC, so adding systematics and dZ"
    if customize.doHTXS:
        variablesToUse = minimalVariablesHTXS
    else:
        variablesToUse = minimalVariables
    if customize.doFiducial:
        variablesToUse.extend(fc.getGenVariables(True))
        variablesToUse.extend(fc.getRecoVariables(True))
        variablesToUse.append("genLeadGenIso := ? diPhoton().leadingPhoton().hasMatchedGenPhoton() ? diPhoton().leadingPhoton().userFloat(\"genIso\") : -99")
        variablesToUse.append("decorrSigmarv := diPhotonMVA().decorrSigmarv")
        variablesToUse.append("leadmva := diPhotonMVA().leadmva")
        variablesToUse.append("subleadmva := diPhotonMVA().subleadmva")
    
    if customize.doSystematics:
        for direction in ["Up","Down"]:
            phosystlabels.append("MvaShift%s01sigma" % direction)
#            phosystlabels.append("MvaLinearSyst%s01sigma" % direction)
            phosystlabels.append("SigmaEOverEShift%s01sigma" % direction)
            phosystlabels.append("MaterialCentralBarrel%s01sigma" % direction)
            phosystlabels.append("MaterialOuterBarrel%s01sigma" % direction)
            phosystlabels.append("MaterialForward%s01sigma" % direction)
            phosystlabels.append("FNUFEB%s01sigma" % direction)
            phosystlabels.append("FNUFEE%s01sigma" % direction)
            phosystlabels.append("MCScaleGain6EB%s01sigma" % direction)
            phosystlabels.append("MCScaleGain1EB%s01sigma" % direction)
            jetsystlabels.append("JEC%s01sigma" % direction)
            jetsystlabels.append("JER%s01sigma" % direction)
            jetsystlabels.append("PUJIDShift%s01sigma" % direction)
            metsystlabels.append("metJecUncertainty%s01sigma" % direction)
            metsystlabels.append("metJerUncertainty%s01sigma" % direction)
            metsystlabels.append("metPhoUncertainty%s01sigma" % direction)
            metsystlabels.append("metUncUncertainty%s01sigma" % direction)
            variablesToUse.append("UnmatchedPUWeight%s01sigma[1,-999999.,999999.] := weight(\"UnmatchedPUWeight%s01sigma\")" % (direction,direction))
            variablesToUse.append("MvaLinearSyst%s01sigma[1,-999999.,999999.] := weight(\"MvaLinearSyst%s01sigma\")" % (direction,direction))
            variablesToUse.append("LooseMvaSF%s01sigma[1,-999999.,999999.] := weight(\"LooseMvaSF%s01sigma\")" % (direction,direction))
            variablesToUse.append("PreselSF%s01sigma[1,-999999.,999999.] := weight(\"PreselSF%s01sigma\")" % (direction,direction))
            variablesToUse.append("electronVetoSF%s01sigma[1,-999999.,999999.] := weight(\"electronVetoSF%s01sigma\")" % (direction,direction))
            variablesToUse.append("TriggerWeight%s01sigma[1,-999999.,999999.] := weight(\"TriggerWeight%s01sigma\")" % (direction,direction))
            variablesToUse.append("FracRVWeight%s01sigma[1,-999999.,999999.] := weight(\"FracRVWeight%s01sigma\")" % (direction,direction))
            variablesToUse.append("FracRVNvtxWeight%s01sigma[1,-999999.,999999.] := weight(\"FracRVNvtxWeight%s01sigma\")" % (direction,direction))
            variablesToUse.append("ElectronWeight%s01sigma[1,-999999.,999999.] := weight(\"ElectronWeight%s01sigma\")" % (direction,direction))
            if os.environ["CMSSW_VERSION"].count("CMSSW_8_0"):
                variablesToUse.append("MuonWeight%s01sigma[1,-999999.,999999.] := weight(\"MuonWeight%s01sigma\")" % (direction,direction))
                variablesToUse.append("MuonMiniIsoWeight%s01sigma[1,-999999.,999999.] := weight(\"MuonMiniIsoWeight%s01sigma\")" % (direction,direction))
            elif os.environ["CMSSW_VERSION"].count("CMSSW_9_4"):
                variablesToUse.append("MuonIDWeight%s01sigma[1,-999999.,999999.] := weight(\"Muon%sIDWeight%s01sigma\")" % (direction,MUON_ID,direction))
                variablesToUse.append("MuonIsoWeight%s01sigma[1,-999999.,999999.] := weight(\"Muon%sISOWeight%s01sigma\")" % (direction,MUON_ISO,direction))
	    variablesToUse.append("JetBTagCutWeight%s01sigma[1,-999999.,999999.] := weight(\"JetBTagCutWeight%s01sigma\")" % (direction,direction))
            variablesToUse.append("JetBTagReshapeWeight%s01sigma[1,-999999.,999999.] := weight(\"JetBTagReshapeWeight%s01sigma\")" % (direction,direction))
            for r9 in ["HighR9","LowR9"]:
                for region in ["EB","EE"]:
                    phosystlabels.append("ShowerShape%s%s%s01sigma"%(r9,region,direction))
#                    phosystlabels.append("MCSmear%s%s%s01sigma" % (r9,region,direction))
                    phosystlabels.append("MCScale%s%s%s01sigma" % (r9,region,direction))
                    for var in ["Rho","Phi"]:
                        phosystlabels.append("MCSmear%s%s%s%s01sigma" % (r9,region,var,direction))
        systlabels += phosystlabels
        systlabels += jetsystlabels
        systlabels += metsystlabels
    customizeSystematicsForSignal(process)
elif customize.processId == "Data":
    print "Data, so turn off all shifts and systematics, with some exceptions"
    variablesToUse = minimalNonSignalVariables
    if customize.doFiducial:
        variablesToUse.extend(fc.getRecoVariables(True))
    customizeSystematicsForData(process)
else:
    print "Background MC, so store mgg and central only"
    variablesToUse = minimalNonSignalVariables
    customizeSystematicsForBackground(process)

if customize.doubleHTagsOnly:
    variablesToUse = minimalVariables
    if customize.processId == "Data":
        variablesToUse = minimalNonSignalVariables

print "--- Systematics  with independent collections ---"
print systlabels
print "-------------------------------------------------"
print "--- Variables to be dumped, including systematic weights ---"
print variablesToUse
print "------------------------------------------------------------"




#from flashgg.Taggers.globalVariables_cff import globalVariables
#globalVariables.extraFloats.rho = cms.InputTag("rhoFixedGridAll")

#cloneTagSequenceForEachSystematic(process,systlabels,phosystlabels,jetsystlabels,jetSystematicsInputTags)
cloneTagSequenceForEachSystematic(process,systlabels,phosystlabels,metsystlabels,jetsystlabels,jetSystematicsInputTags)

# Dump an object called NoTag for untagged events in order to track QCD weights
# Will be broken if it's done for non-central values, so turn this on only for the non-syst tag sorter
process.flashggTagSorter.CreateNoTag = True # MUST be after tag sequence cloning

###### Dumper section

from FWCore.ParameterSet.VarParsing import VarParsing
from flashgg.MetaData.samples_utils import SamplesManager

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(
#"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall17-3_1_0/3_1_0/DoubleEG/RunIIFall17-3_1_0-3_1_0-v0-Run2017F-31Mar2018-v1/180611_135216/0001/myMicroAODOutputFile_1382.root",
#"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall17-3_1_0/3_1_0/DoubleEG/RunIIFall17-3_1_0-3_1_0-v0-Run2017F-31Mar2018-v1/180611_135216/0000/myMicroAODOutputFile_685.root"
#"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall17-3_0_0/3_0_0/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8/RunIIFall17-3_0_0-3_0_0-v0-RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/180325_164819/0000/myMicroAODOutputFile_1.root"
#"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_1-25ns_Moriond17/2_4_1/THQ_HToGG_13TeV-madgraph-pythia8_TuneCUETP8M1/RunIISummer16-2_4_1-25ns_Moriond17-2_4_1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170114_100016/0000/myMicroAODOutputFile_9.root"
#"root://eoscms.cern.ch//eos/cms//store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_1-25ns_Moriond17/2_4_1/bbHToGG_M-125_4FS_ybyt_13TeV_amcatnlo/RunIISummer16-2_4_1-25ns_Moriond17-2_4_1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170114_095140/0000/myMicroAODOutputFile_9.root"
#"root://eoscms.cern.ch//eos/cms//store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_1-25ns_Moriond17/2_4_1/bbHToGG_M-125_4FS_yb2_13TeV_amcatnlo/RunIISummer16-2_4_1-25ns_Moriond17-2_4_1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170114_095013/0000/myMicroAODOutputFile_1.root"
#"root://eoscms.cern.ch//eos/cms//store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_1-25ns_Moriond17/2_4_1/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8/RunIISummer16-2_4_1-25ns_Moriond17-2_4_1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/170113_234241/0000/myMicroAODOutputFile_1.root"
#"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/sethzenz/flashgg/ReMiniAOD-03Feb2017-2_5_0-test/2_5_0/DoubleEG/ReMiniAOD-03Feb2017-2_5_0-test-2_5_0-v0-Run2016G-03Feb2017-v1/170210_054444/0000/myMicroAODOutputFile_264.root"
#"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_1-25ns_Moriond17/2_4_1/VBFHToGG_M-125_13TeV_powheg_pythia8/RunIISummer16-2_4_1-25ns_Moriond17-2_4_1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170114_092754/0000/myMicroAODOutputFile_10.root"
#"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_3_0-25ns_Moriond17_MiniAODv2/2_3_0/DoubleEG/RunIISpring16DR80X-2_3_0-25ns_Moriond17_MiniAODv2-2_3_0-v0-Run2016B-23Sep2016-v2/161114_162452/0000/myMicroAODOutputFile_10.root"
#"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_1-25ns_Moriond17/2_4_1/VHToGG_M120_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16-2_4_1-25ns_Moriond17-2_4_1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170114_094407/0000/myMicroAODOutputFile_19.root"
#"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_3_0-25ns_Moriond17_MiniAODv2/2_3_0/DoubleEG/RunIISpring16DR80X-2_3_0-25ns_Moriond17_MiniAODv2-2_3_0-v0-Run2016B-23Sep2016-v2/161114_162452/0000/myMicroAODOutputFile_10.root"
#"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_1-25ns_Moriond17/2_4_1/VHToGG_M120_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16-2_4_1-25ns_Moriond17-2_4_1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170114_094407/0000/myMicroAODOutputFile_19.root"
#"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_1-25ns_Moriond17/2_4_1/ttHToGG_M125_13TeV_powheg_pythia8_v2/RunIISummer16-2_4_1-25ns_Moriond17-2_4_1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170114_093929/0000/myMicroAODOutputFile_1.root"
#"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISummer16-2_4_1-25ns_Moriond17/2_4_1/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16-2_4_1-25ns_Moriond17-2_4_1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170114_094103/0000/myMicroAODOutputFile_1.root"
#"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2/2_2_0/VBFHToGG_M125_13TeV_amcatnlo_pythia8/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2-2_2_0-v0-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext2-v1/160707_150558/0000/myMicroAODOutputFile_25.root"
#"file:/afs/cern.ch/work/s/sethzenz/fromscratch107/CMSSW_8_0_8_patch1/src/flashgg/myMicroAODOutputFile.root"
#"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_0_0-25ns/2_0_0/VBFHToGG_M-125_13TeV_powheg_pythia8/RunIISpring16DR80X-2_0_0-25ns-2_0_0-v0-RunIISpring16MiniAODv1-PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1/160524_093752/0000/myMicroAODOutputFile_1.root"
#	  "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall17-3_0_0/3_0_0/ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17-3_0_0-3_0_0-v0-RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/180325_144315/0000/myMicroAODOutputFile_1.root"
#         "/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall17-3_0_0/3_0_0/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/RunIIFall17-3_0_0-3_0_0-v0-RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/180325_204512/0000/myMicroAODOutputFile_1.root"
#         "/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall17-3_0_1/3_0_1/GluGluHToGG_M-70_13TeV_powheg_pythia8/RunIIFall17-3_0_1-3_0_1-v0-RunIISummer17MiniAOD-NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/180405_124501/0000/myMicroAODOutputFile_19.root"
#         "/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall17-3_0_1/3_0_1/DoubleEG/RunIIFall17-3_0_1-3_0_1-v0-Run2017F-17Nov2017-v1/180331_074338/0000/myMicroAODOutputFile_1.root"
#        "file:myMicroAODOutputFile.root"
        #        "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-ReMiniAOD-1_1_0-25ns/1_1_0/VBFHToGG_M-125_13TeV_powheg_pythia8/RunIISpring15-ReMiniAOD-1_1_0-25ns-1_1_0-v0-RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/160105_224017/0000/myMicroAODOutputFile_1.root"
#        "root://eoscms.cern.ch//eos/cms//store/group/phys_higgs/cmshgg/szenz/flashgg/RunIISpring15-ReReco74X-Rerun-1_1_0-25ns/1_2_0/DoubleEG/RunIISpring15-ReReco74X-Rerun-1_1_0-25ns-1_2_0-v0-Run2015D-04Dec2015-v2/160117_214114/0000/myMicroAODOutputFile_10.root"
#        "root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-ReMiniAOD-1_1_0-25ns/1_1_0/ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring15-ReMiniAOD-1_1_0-25ns-1_1_0-v0-RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/160105_224456/0000/myMicroAODOutputFile_2.root"
        #"root://eoscms.cern.ch//eos/cms//store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIISpring15-ReMiniAOD-1_1_0-25ns/1_1_0/VHToGG_M120_13TeV_amcatnloFXFX_madspin_pythia8/RunIISpring15-ReMiniAOD-1_1_0-25ns-1_1_0-v0-RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/160105_224138/0000/myMicroAODOutputFile_1.root"
#        "root://eoscms.cern.ch//eos/cms//store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall15DR76-1_3_0-25ns/1_3_0/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIIFall15DR76-1_3_0-25ns-1_3_0-v0-RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/160126_090235/0000/myMicroAODOutputFile_16.root"
#        "root://eoscms.cern.ch//eos/cms//store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIIFall15DR76-1_3_0-25ns_ext1/1_3_1/ttHJetToGG_M120_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall15DR76-1_3_0-25ns_ext1-1_3_1-v0-RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/160127_024939/0000/myMicroAODOutputFile_1.root"
#"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIIFall15DR76-1_3_0-25ns_ext1/1_3_1/DoubleEG/RunIIFall15DR76-1_3_0-25ns_ext1-1_3_1-v0-Run2015D-16Dec2015-v2/160127_022911/0000/myMicroAODOutputFile_100.root"
#"root://eoscms.cern.ch//eos/cms//store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIIFall15DR76-1_3_0-25ns_ext1/1_3_1/VBFHToGG_M-120_13TeV_powheg_pythia8/RunIIFall15DR76-1_3_0-25ns_ext1-1_3_1-v0-RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/160210_045711/0000/myMicroAODOutputFile_1.root"
#"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIIFall15DR76-1_3_0-25ns_ext1/1_3_1/QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIIFall15DR76-1_3_0-25ns_ext1-1_3_1-v0-RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/160127_023721/0000/myMicroAODOutputFile_1.root"
#"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIIFall15DR76-1_3_0-25ns_ext1/1_3_1/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIIFall15DR76-1_3_0-25ns_ext1-1_3_1-v0-RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/160210_050208/0000/myMicroAODOutputFile_1.root"
#"root://eoscms.cern.ch//eos/cms/store/group/phys_higgs/cmshgg/ferriff/flashgg/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2/2_2_0/ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_v2/RunIISpring16DR80X-2_2_0-25ns_ICHEP16_MiniAODv2-2_2_0-v0-RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/160707_152047/0000/myMicroAODOutputFile_8.root"
#"root://cms-xrd-global.cern.ch//store/user/micheli/HHbbgg/MicroAod/RunIIMoriond17_HHbbgg_breg_extra4/1/GluGluToHHTo2B2G_node_box_13TeV-madgraph/RunIIMoriond17_HHbbgg_breg_extra4-1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180528_185053/0000/myMicroAODOutputFile_1.root"
#"root://xrootd-cms.infn.it//store/user/micheli/HHbbgg/MicroAod/RunIIMoriond17_HHbbgg_breg_v2/1/GluGluToHHTo2B2G_node_SM_13TeV-madgraph/RunIIMoriond17_HHbbgg_breg-1-v1-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180412_131554/0000/myMicroAODOutputFile_1.root"
#"root://xrootd-cms.infn.it//store/group/phys_higgs/cmshgg/micheli/flashgg/LegacyReReco-20180629/1/DoubleEG/LegacyReReco-20180629-1-v0-Run2016E-07Aug17-v1/180629_143351/0000/myMicroAODOutputFile_44.root"
#"root://xrootd-cms.infn.it//store/group/phys_higgs/cmshgg/micheli/flashgg/LegacyReReco-20180629/1/DoubleEG/LegacyReReco-20180629-1-v0-Run2016E-07Aug17-v1/180629_143351/0000/myMicroAODOutputFile_247.root"
#"root://cms-xrd-global.cern.ch//store/user/micheli/HHbbgg/MicroAod/RunIIMoriond17_HHbbgg_breg_extra4/1/GluGluToHHTo2B2G_node_13_13TeV-madgraph/RunIIMoriond17_HHbbgg_breg_extra4-1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180528_183109/0000/myMicroAODOutputFile_1.root"
#"root://cms-xrd-global.cern.ch////store/group/phys_higgs/HiggsExo/HH_bbgg/RunIIFall17-3_1_0/3_1_0/GluGluToHHTo2B2G_node_SM_13TeV-madgraph/RunIIFall17-3_1_0-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/180830_190720/0000/myMicroAODOutputFile_1.root"
#"root://cms-xrd-global.cern.ch////store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall17-3_1_0/3_1_0/DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa/RunIIFall17-3_1_0-3_1_0-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/180706_100515/0000/myMicroAODOutputFile_253.root"
#"root://cms-xrd-global.cern.ch////store/group/phys_higgs/cmshgg/spigazzi/flashgg/RunIIFall17-3_2_0/RunIIFall17-3_2_0/GluGluToHHTo2B2G_node_7_13TeV-madgraph_correctedcfg/RunIIFall17-3_2_0-RunIIFall17-3_2_0-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/181011_162430/0000/myMicroAODOutputFile_3.root"
#"root://cms-xrd-global.cern.ch///store/group/phys_higgs/cmshgg/spigazzi/flashgg/RunIIFall17-3_2_0/RunIIFall17-3_2_0/GluGluToHHTo2B2G_node_12_13TeV-madgraph_correctedcfg/RunIIFall17-3_2_0-RunIIFall17-3_2_0-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/181011_161750/0000/myMicroAODOutputFile_1.root"
#"root://xrootd-cms.infn.it///store/user/micheli/HHbbgg/MicroAod/RunIIMoriond17_HHbbgg_breg_v2/1/GluGluToHHTo2B2G_node_SM_13TeV-madgraph/RunIIMoriond17_HHbbgg_breg-1-v1-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180412_131554/0000/myMicroAODOutputFile_3.root"
#"root://xrootd-cms.infn.it///store/user/micheli/HHbbgg/MicroAod/HHbbgg_Signal_SM_20181120/1/GluGluToHHTo2B2G_node_SM_13TeV-madgraph/HHbbgg_Signal_SM_20181120-1-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/181120_113738/0000/myMicroAODOutputFile_1.root"
#"root://cms-xrd-global.cern.ch///store/group/phys_higgs/cmshgg/sethzenz/flashgg/RunIIFall17-3_1_0/3_1_0/DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa/RunIIFall17-3_1_0-3_1_0-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/180706_100515/0000/myMicroAODOutputFile_478.root"
#"root://cms-xrd-global.cern.ch///store/user/micheli/HHbbgg/MicroAod/RunIIMoriond17_HHbbgg_breg_extra_5/1/VBFHToGG_M-125_13TeV_powheg_pythia8/RunIIMoriond17_HHbbgg_breg_extra_5-1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180530_092552/0000/myMicroAODOutputFile_6.root"
#"root://cms-xrd-global.cern.ch///store/group/phys_higgs/cmshgg/spigazzi/flashgg/RunIIFall17-3_2_0/RunIIFall17-3_2_0/DoubleEG/RunIIFall17-3_2_0-RunIIFall17-3_2_0-v0-Run2017E-31Mar2018-v1/181008_110407/0001/myMicroAODOutputFile_1661.root"
#"root://cms-xrd-global.cern.ch///store/user/micheli/MicroAOD/ReMiniAOD2016-DeepCSV-bRegression/ReMiniAOD2016-DeepCSV-bRegression/prod-uAOD-300-11-g22a116d/DoubleEG/ReMiniAOD2016-DeepCSV-bRegression-prod-uAOD-300-11-g22a116d-v0-Run2016C-03Feb2017-v1/180928_153054/0000/myMicroAODOutputFile_420.root"
#"root://cms-xrd-global.cern.ch//store/user/micheli/MicroAOD/ReMiniAOD2016-DeepCSV-bRegression/ReMiniAOD2016-DeepCSV-bRegression/prod-uAOD-300-11-g22a116d/DoubleEG/ReMiniAOD2016-DeepCSV-bRegression-prod-uAOD-300-11-g22a116d-v0-Run2016G-03Feb2017-v1/180928_153415/0000/myMicroAODOutputFile_57.root"
#"root://cms-xrd-global.cern.ch//store/user/micheli/HHbbgg/MicroAod/RunIIMoriond17_HHbbgg_breg_extra4/1/GluGluToHHTo2B2G_node_4_13TeV-madgraph/RunIIMoriond17_HHbbgg_breg_extra4-1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180528_183747/0000/myMicroAODOutputFile_1.root"
#"root://cms-xrd-global.cern.ch////store/user/micheli/HHbbgg/MicroAod/RunIIMoriond17_HHbbgg_breg_extra_5/1/ttHToGG_M125_13TeV_powheg_pythia8_v2/RunIIMoriond17_HHbbgg_breg_extra_5-1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/180530_093447/0000/myMicroAODOutputFile_2.root"
#"/store/group/phys_higgs/HiggsExo/HH_bbgg/RunIIFall17-3_1_0/3_1_0/GluGluToHHTo2B2G_node_SM_13TeV-madgraph/RunIIFall17-3_1_0-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/180830_190720/0000/myMicroAODOutputFile_2.root"
#"root://cms-xrd-global.cern.ch//store/user/micheli/HHbbgg/MicroAod/HHbbgg_Signal_SM_20181120_splitting/1/GluGluToHHTo2B2G_node_SM_13TeV-madgraph/HHbbgg_Signal_SM_20181120-1-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/190117_130615/0000/myMicroAODOutputFile_21.root"
#"root://cms-xrd-global.cern.ch//store/user/micheli/HHbbgg/MicroAod/HHbbgg_Signal_SM_20181120_splitting/1/GluGluToHHTo2B2G_node_SM_13TeV-madgraph/HHbbgg_Signal_SM_20181120-1-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/190117_130615/0000/myMicroAODOutputFile_15.root"
"/store/user/micheli/HHbbgg/MicroAod/RunIIMoriond17_HHbbgg_breg_extra_8/1/GluGluToHHTo2B2G_node_10_13TeV-madgraph/RunIIMoriond17_HHbbgg_breg_extra_8-1-v0-RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/190228_195437/0000/myMicroAODOutputFile_2.root"
))

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test.root"))

process.extraDumpers = cms.Sequence()
process.load("flashgg.Taggers.diphotonTagDumper_cfi") ##  import diphotonTagDumper 
import flashgg.Taggers.dumperConfigTools as cfgTools


process.tagsDumper.className = "DiPhotonTagDumper"
process.tagsDumper.src = "flashggSystTagMerger"
#process.tagsDumper.src = "flashggTagSystematics"
process.tagsDumper.processId = "test"
process.tagsDumper.dumpTrees = customize.dumpTrees
process.tagsDumper.dumpWorkspace = customize.dumpWorkspace
process.tagsDumper.dumpHistos = False
process.tagsDumper.quietRooFit = True
process.tagsDumper.nameTemplate = cms.untracked.string("$PROCESS_$SQRTS_$CLASSNAME_$SUBCAT_$LABEL")
process.tagsDumper.splitPdfByStage0Cat = cms.untracked.bool(customize.doHTXS)

if customize.options.WeightName :
    lheProduct = customize.dataset[1]["LHESourceName"].split("_")
    #print lheProduct
    process.tagsDumper.LHEEventProduct = cms.untracked.InputTag( str(lheProduct[1]) , str(lheProduct[2]) , str(lheProduct[3]) )
    #print process.tagsDumper.LHEEventProduct
    process.tagsDumper.LHEWeightName = cms.untracked.string(customize.options.WeightName)


if(customize.doFiducial):
#    if customize.processId == "Data":
#        fc.addRecoGlobalVariables(process, process.tagsDumper)
#    else:
    fc.addObservables(process, process.tagsDumper, customize.processId )

#tagList=[
#["UntaggedTag",4],
#["VBFTag",2],
#["VHTightTag",0],
#["VHLooseTag",0],
#["VHEtTag",0],
#["VHHadronicTag",0],
#["TTHHadronicTag",0],
##["TTHLeptonicTag",0]
#]


if customize.doFiducial:
    tagList=[["SigmaMpTTag",3]]
elif customize.tthTagsOnly:
    tagList=[
        ["NoTag",0],
        ["TTHHadronicTag",3],
        ["TTHLeptonicTag",2],
        ["TTHDiLeptonTag",0]
        ]
elif customize.doubleHTagsOnly:
    tagList = hhc.tagList(customize,process)
    print "taglist is:"
    print tagList
else:
    tagList=[
        ["NoTag",0],
        ["UntaggedTag",4],
        ["VBFTag",3],
        ["ZHLeptonicTag",0],
        ["WHLeptonicTag",0],
        ["VHLeptonicLooseTag",0],
        ["VHMetTag",0],
        ["VHHadronicTag",0],
        ["TTHHadronicTag",3],
        ["TTHLeptonicTag",2],
        ["TTHDiLeptonTag",0]
        ]

definedSysts=set()
process.tagsDumper.NNLOPSWeightFile=cms.FileInPath("flashgg/Taggers/data/NNLOPS_reweight.root")
process.tagsDumper.reweighGGHforNNLOPS = cms.untracked.bool(bool(customize.processId.count("ggh")))
process.tagsDumper.classifierCfg.remap=cms.untracked.VPSet()
for tag in tagList: 
  tagName=tag[0]
  tagCats=tag[1]
  # remap return value of class-based classifier
  process.tagsDumper.classifierCfg.remap.append( cms.untracked.PSet( src=cms.untracked.string("flashgg%s"%tagName), dst=cms.untracked.string(tagName) ) )
  for systlabel in systlabels:
      if not systlabel in definedSysts:
          # the cut corresponding to the systematics can be defined just once
          cutstring = "hasSyst(\"%s\") "%(systlabel)
          definedSysts.add(systlabel)
      else:
          cutstring = None
      if systlabel == "":
          currentVariables = variablesToUse
      else:
          if customize.doHTXS:
              currentVariables = systematicVariablesHTXS
          else:    
              currentVariables = systematicVariables
      if tagName == "NoTag":
          if customize.doHTXS:
              currentVariables = ["stage0cat[72,9.5,81.5] := tagTruth().HTXSstage0cat"]
          else:
              currentVariables = []
      isBinnedOnly = (systlabel !=  "")
     # if ( customize.doPdfWeights or customize.doSystematics ) and ( (customize.datasetName() and customize.datasetName().count("HToGG")) or customize.processId.count("h_") or customize.processId.count("vbf_") ) and (systlabel ==  "") and not (customize.processId == "th_125" or customize.processId == "bbh_125"):
      if ( customize.doPdfWeights) and (customize.doSystematics ) and ( (customize.datasetName() and customize.datasetName().count("HToGG")) or customize.processId.count("h_") or customize.processId.count("vbf_") ) and (systlabel ==  "") and not (customize.processId == "th_125" or customize.processId == "bbh_125"):
          print "Signal MC central value, so dumping PDF weights"
          dumpPdfWeights = True
          nPdfWeights = 60
          nAlphaSWeights = 2
          nScaleWeights = 9
      else:
          print "Data, background MC, or non-central value, or no systematics: no PDF weights"
          dumpPdfWeights = False
          nPdfWeights = -1
          nAlphaSWeights = -1
          nScaleWeights = -1
      cfgTools.addCategory(process.tagsDumper,
                           systlabel,
                           classname=tagName,
                           cutbased=cutstring,
                           subcats=tagCats, 
                           variables=currentVariables,
                           histograms=minimalHistograms,
                           binnedOnly=isBinnedOnly,
                           dumpPdfWeights=dumpPdfWeights,
                           nPdfWeights=nPdfWeights,
                           nAlphaSWeights=nAlphaSWeights,
                           nScaleWeights=nScaleWeights,
                           splitPdfByStage0Cat=customize.doHTXS
                           )

# Require standard diphoton trigger
trigger_data = 'HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*'
if year=='2016':
    trigger_data = 'HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*'
elif year=='2017':
    trigger_data = 'HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v*'

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring("%s"%trigger_data,
#                                                                "HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v1",
#                                                                "HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v1"
                                                                ))

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# ee bad supercluster filter on data
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.EERecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits") # Saved MicroAOD Collection (data only)
# Bad Muon filter LOADS WRONG IN 8_0_28, FIX LATER
#process.load('RecoMET.METFilters.badGlobalMuonTaggersMiniAOD_cff')
#process.badGlobalMuonTaggerMAOD.muons = cms.InputTag("flashggSelectedMuons")
#process.cloneGlobalMuonTaggerMAOD.muons = cms.InputTag("flashggSelectedMuons")
process.dataRequirements = cms.Sequence()
if customize.processId == "Data":
        process.dataRequirements += process.hltHighLevel
        process.dataRequirements += process.eeBadScFilter
#        if customize.doMuFilter:
#            process.dataRequirements += process.noBadGlobalMuonsMAOD

# Split WH and ZH
process.genFilter = cms.Sequence()
if ((customize.processId.count("wh") or customize.processId.count("zh")) and not (customize.processId.count("powheg"))) and not customize.processId.count("wzh") :
    print "enabling vh filter!!!!!"
    process.load("flashgg/Systematics/VHFilter_cfi")
    process.genFilter += process.VHFilter
    process.VHFilter.chooseW = bool(customize.processId.count("wh"))
    process.VHFilter.chooseZ = bool(customize.processId.count("zh"))

if (customize.processId == "th_125" or customize.processId == "bbh_125"):
    process.load("flashgg/Systematics/CentralHiggsFilter_cfi")
    process.genFilter += process.CentralHiggsFilter

#pythia8 has an unanticipated EM showering feature, check have two photons from hard scatter
process.penultimateFilter= cms.Sequence()
if customize.processId == "th_125": # for this sample the filter removes also H -> ZG
    process.load("flashgg/Systematics/HardProcessFinalStateFilter_cfi")
#    process.HardProcessFinalStateFilter.debug = True
    process.penultimateFilter += process.HardProcessFinalStateFilter

# Split out prompt-fake or fake-fake
process.finalFilter = cms.Sequence()
if (customize.processId.count("qcd") or customize.processId.count("gjet")) and customize.processId.count("fake"):
    process.load("flashgg/Systematics/PromptFakeFilter_cfi")
    process.finalFilter += process.PromptFakeFilter
    if (customize.processId.count("promptfake")):
        process.PromptFakeFilter.doPromptFake = cms.bool(True)
        process.PromptFakeFilter.doFakeFake =cms.bool(False)
    elif (customize.processId.count("fakefake")):
        process.PromptFakeFilter.doPromptFake =cms.bool(False)
        process.PromptFakeFilter.doFakeFake =cms.bool(True)
    else:
        raise Exception,"Mis-configuration of python for prompt-fake filter"


if customize.tthTagsOnly:
    process.p = cms.Path(process.dataRequirements*
                         process.genFilter*
                         process.flashggDiPhotons* # needed for 0th vertex from microAOD
                         process.flashggUpdatedIdMVADiPhotons*
                         process.flashggDiPhotonSystematics*
                         process.flashggMetSystematics*
                         process.flashggMuonSystematics*process.flashggElectronSystematics*
                         (process.flashggUnpackedJets*process.jetSystematicsSequence)*
                         (process.flashggTagSequence*process.systematicsTagSequences)*
                         process.flashggSystTagMerger*
                         process.penultimateFilter*
                         process.finalFilter*
                         process.tagsDumper)
else :
    process.p = cms.Path(process.dataRequirements*
                         process.genFilter*
                         process.flashggUpdatedIdMVADiPhotons*
                         process.flashggDiPhotonSystematics*
                         process.flashggMetSystematics*
                         process.flashggMuonSystematics*process.flashggElectronSystematics*
                         (process.flashggUnpackedJets*process.jetSystematicsSequence)*
                         (process.flashggTagSequence*process.systematicsTagSequences)*
                         process.flashggSystTagMerger*
                         process.penultimateFilter*
                         process.finalFilter*
                         process.tagsDumper)

#if customize.doBJetRegression:
#
#    bregProducers = []
#    bregJets = []
#    
#    from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag
#    from flashgg.Taggers.flashggbRegressionProducer_cfi import flashggbRegressionProducer
#    recoJetCollections = UnpackedJetCollectionVInputTag
#
#
#    for icoll,coll in enumerate(recoJetCollections):
#        print "doing icoll "+str(icoll)
#
#        producer = flashggbRegressionProducer.clone(JetTag = coll)
#
#        setattr(process,"bRegProducer%d" %icoll,producer)
#        bregProducers.append(producer)
#        bregJets.append("bRegProducer%d" %icoll)
#            
#    process.bregProducers = cms.Sequence(reduce(lambda x,y: x+y, bregProducers))
##    process.bbggtree.inputTagJets=cms.VInputTag(bregJets)
#    process.p.replace(process.jetSystematicsSequence,process.jetSystematicsSequence*process.flashggUnpackedJets+process.bregProducers)
##    process.p.replace(process.jetSystematicsSequence,process.flashggUnpackedJets*process.jetSystematicsSequence+process.bregProducers)  #tried but no difference
#

if customize.doBJetRegression:

    bregProducers = []
    doubleHTagProducers = []
    
    from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag
    from flashgg.Taggers.flashggbRegressionProducer_cfi import flashggbRegressionProducer
    recoJetCollections = UnpackedJetCollectionVInputTag

    from flashgg.Taggers.flashggDoubleHTag_cfi import flashggDoubleHTag


    for jetsyst in ["","JECDown01sigma","JECUp01sigma","JERDown01sigma","JERUp01sigma"]: 
    #for jetsyst in ['']:
    # for jetsyst in [systlabel[0]]+jetsystlabels:
       jetCollection = cms.VInputTag()
       for icoll,coll in enumerate(recoJetCollections):
            jetCollection.append(cms.InputTag(coll.moduleLabel,jetsyst))
       if jetsyst != "" : producer = flashggbRegressionProducer.clone(JetTags = jetCollection)
       else : producer = flashggbRegressionProducer.clone(JetTags = recoJetCollections)

       setattr(process,"bRegProducer%s" %(jetsyst),producer)
       bregProducers.append(producer)
           
    process.bregProducers = cms.Sequence(reduce(lambda x,y: x+y, bregProducers))
    process.p.replace(process.jetSystematicsSequence,process.jetSystematicsSequence*process.flashggUnpackedJets+process.bregProducers)
        
    if jetsystlabels!=[]:
       # for jetsyst in [systlabel[0]]+jetsystlabels:
        for jetsyst in ["","JECDown01sigma","JECUp01sigma","JERDown01sigma","JERUp01sigma"]: 
            jetTagsSystematics = cms.VInputTag("bRegProducer%s" %(jetsyst))
            getattr(process, "flashggDoubleHTag"+jetsyst).JetTags = jetTagsSystematics
            
        print "New tag sequence  : ", process.flashggTagSequence 
        print "New syst merger src  : ", process.flashggSystTagMerger.src 
        print "New systematics tag sequence  : ", process.systematicsTagSequences



if customize.doFiducial:
    if ( customize.doPdfWeights or customize.doSystematics ) and ( (customize.datasetName() and customize.datasetName().count("HToGG")) 
                                                                   or customize.processId.count("h_") or customize.processId.count("vbf_") ) and (systlabel ==  ""):
          print "Signal MC central value, so dumping PDF weights"
          dumpPdfWeights = True
          nPdfWeights = 60
          nAlphaSWeights = 2
          nScaleWeights = 9
    else:
          print "Data, background MC, or non-central value, or no systematics: no PDF weights"
          dumpPdfWeights = False
          nPdfWeights = -1
          nAlphaSWeights = -1
          nScaleWeights = -1
    if not customize.processId == "Data":
        fc.addGenOnlyAnalysis(process,customize.processId,customize.acceptance,tagList,systlabels,pdfWeights=(dumpPdfWeights,nPdfWeights,nAlphaSWeights,nScaleWeights))

if customize.doubleHReweight>0:
    hhc.addNodesReweighting(customize,process)

if customize.doDoubleHGenAnalysis:
    hhc.addGenAnalysis(customize,process,tagList)


if( not hasattr(process,"options") ): process.options = cms.untracked.PSet()
process.options.allowUnscheduled = cms.untracked.bool(True)

print "--- Dumping modules that take diphotons as input: ---"
mns = process.p.moduleNames()
for mn in mns:
    module = getattr(process,mn)
    if hasattr(module,"src") and type(module.src) == type(cms.InputTag("")) and module.src.value().count("DiPhoton"):
        print str(module),module.src
    elif hasattr(module,"DiPhotonTag"):
        print str(module),module.DiPhotonTag
print
printSystematicInfo(process)

# Detailed tag interpretation information printout (blinded)
process.flashggTagSorter.StoreOtherTagInfo = True
process.flashggTagSorter.BlindedSelectionPrintout = True

#### BELOW HERE IS MOSTLY DEBUGGING STUFF

#####################################################################
## Memory and timing, n.b. igprof is very complicated to interpret ##
##################################################################### 

#from Validation.Performance.TimeMemoryInfo import customise as TimeMemoryCustomize
#TimeMemoryCustomize(process)
#process.MessageLogger.cerr.threshold = 'WARNING'
#
#process.load("IgTools.IgProf.IgProfTrigger")
#process.igprof.reportEventInterval     = cms.untracked.int32(1000)
#process.igprof.reportToFileAtBeginJob  = cms.untracked.string("|gzip -c>igprof.begin-job.gz")
#process.igprof.reportToFileAtEvent     = cms.untracked.string("|gzip -c>igprof.%I.%E.%L.%R.event.gz")
#process.p += process.igprof
#
################################
## Dump merged tags to screen ##
################################

if customize.verboseTagDump:
    # crashes right now, dunno why - 02 May 2018
    pass
#    process.load("flashgg/Taggers/flashggTagTester_cfi")
#    process.flashggTagTester.TagSorter = cms.InputTag("flashggSystTagMerger")
#    process.flashggTagTester.ExpectMultiples = cms.untracked.bool(True)
#    process.p += process.flashggTagTester

############################################
## Additional details on tag sorter steps ##
############################################

if customize.verboseTagDump:
    process.flashggUpdatedIdMVADiPhotons.Debug = True
    process.flashggTagSorter.Debug = True
    customize.maxEvents = 10
                           
if customize.verboseSystDump:
    turnOnAllSystematicsDebug(process)
    customize.maxEVents = 10

##############
## Dump EDM ##
##############

#process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('CustomizeWillChangeThisAnyway.root'),
#                               outputCommands = cms.untracked.vstring('keep *') # dump everything! small tests only!
#                               )
#process.e = cms.EndPath(process.out)

############################
## Dump the output Python ##
############################
#print process.dumpPython()
#processDumpFile = open('processDump.py', 'w')
#print >> processDumpFile, process.dumpPython()

# call the customization
customize(process)
