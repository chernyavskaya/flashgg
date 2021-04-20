import FWCore.ParameterSet.Config as cms
from  flashgg.Systematics.flashggJetSystematics_cfi import jetSystematicsCustomize

class cHCustomize():
    """
    cH process customizaton class
    """
    
    def __init__(self, process, customize, metaConditions):
        self.process = process
        self.customize = customize
        self.metaConditions = metaConditions
        self.tagList = [ ["cHTag",1] ]
        self.customizeTagSequence()


    def variablesToDump(self):
        var_workspace = []
        variables = []
        if(self.customize.cHTagsOnly):
            var_workspace += [
                "eventNumber := eventNumber()",
                "leadingJet_pt := leadJet().pt",
            ]
            if self.customize.processId != "Data":
                var_workspace += [
                    'btagReshapeWeight := weight("JetBTagReshapeWeightCentral")',
                ]
                variables += [
                    "leadingJet_hflav := leadJet().hadronFlavour()",
                    "leadingJet_pflav := leadJet().partonFlavour()",
                    'btagReshapeWeight := weight("JetBTagReshapeWeightCentral")',
                ]
            variables += [
                "leadingJet_DeepFlavour_CvsL := leadJet().bDiscriminator('mini_pfDeepFlavourJetTags:probc')/(leadJet().bDiscriminator('mini_pfDeepFlavourJetTags:probc')+leadJet().bDiscriminator('mini_pfDeepFlavourJetTags:probuds')+leadJet().bDiscriminator('mini_pfDeepFlavourJetTags:probg'))",#FIXME make the ctag type configurable?
                "leadingJet_DeepFlavour_CvsB := leadJet().bDiscriminator('mini_pfDeepFlavourJetTags:probc')/(leadJet().bDiscriminator('mini_pfDeepFlavourJetTags:probc')+leadJet().bDiscriminator('mini_pfDeepFlavourJetTags:probb')+leadJet().bDiscriminator('mini_pfDeepFlavourJetTags:probbb')+leadJet().bDiscriminator('mini_pfDeepFlavourJetTags:problepb'))",#FIXME make the ctag type configurable?
                "leadingJet_puJetIdMVA := leadJet().puJetIdMVA()",
                "leadingJet_puJetIdMVA := leadJet().puJetIdMVA()",
                "diphotonCandidatePtOverdiHiggsM := diphotonPtOverM()",
                "customLeadingPhotonIDMVA := diPhoton.leadingView.phoIdMvaWrtChosenVtx",
                "customSubLeadingPhotonIDMVA := diPhoton.subLeadingView.phoIdMvaWrtChosenVtx",
                "EGMLeadingPhotonIDMVA := diPhoton.leadingPhoton.userFloat('EGMPhotonMVA')",
                "EGMSubLeadingPhotonIDMVA := diPhoton.subLeadingPhoton.userFloat('EGMPhotonMVA')",
                "leadingPhotonSigOverE := diPhoton.leadingPhoton.sigEOverE",
                "subleadingPhotonSigOverE := diPhoton.subLeadingPhoton.sigEOverE",
                "sigmaMOverM := sqrt(0.5*(diPhoton.leadingPhoton.sigEOverE*diPhoton.leadingPhoton.sigEOverE + diPhoton.subLeadingPhoton.sigEOverE*diPhoton.subLeadingPhoton.sigEOverE))",
                "sigmaMOverMDecorr := getSigmaMDecorr()",
                "PhoJetMinDr := getPhoJetMinDr()",#up to here input variables to MVA
                "leadingPhoton_pt := diPhoton.leadingPhoton.pt",
                "leadingPhoton_eta := diPhoton.leadingPhoton.eta",
                "leadingPhoton_phi := diPhoton.leadingPhoton.phi",
                "subleadingPhoton_pt := diPhoton.subLeadingPhoton.pt",
                "subleadingPhoton_eta := diPhoton.subLeadingPhoton.eta",
                "subleadingPhoton_phi := diPhoton.subLeadingPhoton.phi",
                
                "leadingJet_pt := leadJet().pt",
                "leadingJet_eta := leadJet().eta",
                "leadingJet_phi := leadJet().phi",
                "leadingJet_mass := leadJet().p4().M()",
                
                "nJets := nJets()",
            ]
        if  self.customize.dumpWorkspace :
            return var_workspace
        else :
            return variables


    def systematicVariables(self):
      systematicVariables=["CMS_hgg_mass[160,100,180]:=diPhoton().mass","eventNumber[40,0.,1000000.]:=eventNumber()",'btagReshapeWeight[100,-10.,10]:=weight("JetBTagReshapeWeightCentral")',"nJets[100,0.,10] := nJets()","leadingJet_pt[100,0,1000] := leadJet().pt"]
      
      return systematicVariables


    def variablesToDumpData(self):
        variables = [ #placeholder
        ]

        if not  self.customize.dumpWorkspace:
            return self.variablesToDump()

        return variables


    def customizeSystematics(self,systlabels,jetsystlabels,metsystlabels):
       for s in metsystlabels:
          systlabels.remove(s)
       metsystlabels = []
       return systlabels,jetsystlabels,metsystlabels

    def customizeTagSequence(self):
        self.process.load("flashgg.Taggers.flashggcHTag_cff")

        self.process.flashggcHTag.JetIDLevel=cms.string(str(self.metaConditions["doubleHTag"]["jetID"]))
        self.process.flashggcHTag.MaxJetEta = cms.double(self.metaConditions["bTagSystematics"]["eta"])

        ## add cH tag to the tag sequence
        ## remove single Higgs tags
        if self.customize.cHTagsOnly:
            self.process.flashggTagSequence.remove(self.process.flashggVBFTag)
            self.process.flashggTagSequence.remove(self.process.flashggTTHLeptonicTag)
            self.process.flashggTagSequence.remove(self.process.flashggTTHHadronicTag)
            self.process.flashggTagSequence.remove(self.process.flashggVHEtTag)
            self.process.flashggTagSequence.remove(self.process.flashggVHLooseTag)
            self.process.flashggTagSequence.remove(self.process.flashggVHTightTag)
            self.process.flashggTagSequence.remove(self.process.flashggVHMetTag)
            self.process.flashggTagSequence.remove(self.process.flashggWHLeptonicTag)
            self.process.flashggTagSequence.remove(self.process.flashggZHLeptonicTag)
            self.process.flashggTagSequence.remove(self.process.flashggVHLeptonicLooseTag)
            self.process.flashggTagSequence.remove(self.process.flashggVHHadronicTag)
            self.process.flashggTagSequence.remove(self.process.flashggVBFMVA)
            self.process.flashggTagSequence.remove(self.process.flashggVBFDiPhoDiJetMVA)
            self.process.flashggTagSequence.remove(self.process.flashggTTHDiLeptonTag)
            self.process.flashggTagSequence.remove(self.process.flashggUntagged)
            self.process.flashggTagSequence.remove(self.process.flashggUntagged)
            self.process.flashggTagSequence.remove(self.process.flashggTHQLeptonicTag)
 
    def cHTagMerger(self,systlabels=[]):
        '''
        Construct the actual tag sequence for the cH analysis. setting up the
        merging step taking into account that different syst variations are produced by the same producer in the case of the H tags
        '''

        self.process.p.remove(self.process.flashggTagSorter) 
        self.process.p.replace(self.process.flashggSystTagMerger,self.process.flashggcHTagSequence*self.process.flashggTagSorter*self.process.flashggSystTagMerger)

        for systlabel in systlabels:
            if systlabel!='':
                self.process.p.remove(getattr(self.process,'flashggTagSorter'+systlabel))
                self.process.p.replace(self.process.flashggSystTagMerger,getattr(self.process, 'flashggTagSorter'+systlabel)*self.process.flashggSystTagMerger)            
            setattr(getattr(self.process, 'flashggTagSorter'+systlabel), 'TagPriorityRanges', cms.VPSet( cms.PSet(TagName = cms.InputTag('flashggcHTag', systlabel))))
                 

    def cHTagRunSequence(self,systlabels,jetsystlabels,phosystlabels):
        if self.customize.cHTagsOnly: 
            self.cHTagMerger(systlabels)

        if len(systlabels)>1 :
            getattr(self.process, "flashggcHTag").JetsSuffixes = cms.vstring([systlabels[0]]+jetsystlabels)
            getattr(self.process, "flashggcHTag").DiPhotonSuffixes = cms.vstring([systlabels[0]]+phosystlabels)


