import FWCore.ParameterSet.Config as cms

import flashgg.Taggers.dumperConfigTools as cfgTools

from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag
import flashgg.Taggers.PUJID_wps as pujid

from flashgg.Taggers.globalVariables_cff import globalVariables
from flashgg.MicroAOD.flashggJets_cfi import  maxJetCollections


jetID = ''

MaxJetEta = 2.5

flashggcHTag = cms.EDProducer("FlashggcHTagProducer",
                                   DiPhotonName = cms.string('flashggPreselectedDiPhotons'), # 
                                   DiPhotonSuffixes = cms.vstring(''), #nominal and systematic variations 
                                   JetsName = cms.string("flashggUnpackedJets"), # 
                                   JetsCollSize = cms.uint32(maxJetCollections), #
                                   JetsSuffixes = cms.vstring(''), #nominal and systematic variations 
                                   GenParticleTag = cms.InputTag( "flashggPrunedGenParticles" ), # to compute MC-truth info
                                   
                                   VetoConeSize   = cms.double(0.4),
                                   MinLeadPhoPt   = cms.double(1./3.),
                                   MinSubleadPhoPt   = cms.double(0.25),
                                   ScalingPtCuts = cms.bool(True),
                                   DoSigmaMDecorr =cms.untracked.uint32(1),#transformation of sigmaM/M
                                   SigmaMDecorrFile = cms.untracked.FileInPath("flashgg/Taggers/data/diphoMVA_sigmaMoMdecorr_split_Mgg40_180.root"),
                                   ApplyEGMPhotonID = cms.untracked.bool(False),
                                   PhotonIDCut = cms.double(0.2),#this is loose id for 2016
                                   PhotonElectronVeto =cms.untracked.vint32(1, 1), #0: Pho1, 1: Pho2

                                   MinJetPt   = cms.double(20.),
                                   MaxJetEta   = cms.double(MaxJetEta),
                                   CTagTypeCvsL = cms.vstring('mini_pfDeepFlavourJetTags:probc','mini_pfDeepFlavourJetTags:probuds','mini_pfDeepFlavourJetTags:probg'), #string for ctag algorithm
                                   CTagTypeCvsB  = cms.vstring('mini_pfDeepFlavourJetTags:probc','mini_pfDeepFlavourJetTags:probb','mini_pfDeepFlavourJetTags:probbb','mini_pfDeepFlavourJetTags:problepb'), #string for ctag algorithm
                                   UseJetID = cms.bool(True),
                                   JetIDLevel = cms.string(jetID),
                                   globalVariables=globalVariables,
                                  ) 

