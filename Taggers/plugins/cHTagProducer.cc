#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"

#include "flashgg/DataFormats/interface/DiPhotonMVAResult.h"
#include "flashgg/DataFormats/interface/cHTag.h"
#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "flashgg/Taggers/interface/LeptonSelection.h"
#include "flashgg/MicroAOD/interface/MVAComputer.h"
#include "flashgg/DataFormats/interface/DoubleHttHTagger.h"

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

#include <vector>
#include <algorithm>
#include "TGraph.h"

using namespace std;
using namespace edm;

namespace flashgg {

    class cHTagProducer : public EDProducer
    {

    public:
        typedef math::XYZPoint Point;

        cHTagProducer( const ParameterSet & );
    private:
        void produce( Event &, const EventSetup & ) override;
        
        std::string inputJetsName_;
        std::vector<std::string> inputJetsSuffixes_;
        unsigned int inputJetsCollSize_;
        std::vector<edm::EDGetTokenT<edm::View<flashgg::Jet> > > jetTokens_;
        std::string inputDiPhotonName_;
        std::vector<std::string> inputDiPhotonSuffixes_;
       // EDGetTokenT<View<DiPhotonCandidate> > diPhotonToken_;
        std::vector<edm::EDGetTokenT<edm::View<DiPhotonCandidate> > > diPhotonTokens_;

        EDGetTokenT<View<reco::GenParticle> > genParticleToken_;

        std::vector< std::string > systematicsLabels;

        double minLeadPhoPt_, minSubleadPhoPt_;
        bool scalingPtCuts_, doPhotonId_;
        double photonIDCut_;
        double vetoConeSize_;         
        unsigned int doSigmaMDecorr_;
        edm::FileInPath sigmaMDecorrFile_;
        std::vector<int> photonElectronVeto_;


        DecorrTransform* transfEBEB_;
        DecorrTransform* transfNotEBEB_;

        double minJetPt_;
        double maxJetEta_;
        vector<std::string> cTagTypeCvsL_;
        vector<std::string> cTagTypeCvsB_;
        bool       useJetID_;
        string     JetIDLevel_;        

        ConsumesCollector cc_;
        GlobalVariablesComputer globalVariablesComputer_;
            

    };

    cHTagProducer::cHTagProducer( const ParameterSet &iConfig ) :
      //  diPhotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
        genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleTag" ) ) ),
        minLeadPhoPt_( iConfig.getParameter<double> ( "MinLeadPhoPt" ) ),
        minSubleadPhoPt_( iConfig.getParameter<double> ( "MinSubleadPhoPt" ) ),
        scalingPtCuts_( iConfig.getParameter<bool> ( "ScalingPtCuts" ) ),
        vetoConeSize_( iConfig.getParameter<double> ( "VetoConeSize" ) ),
        minJetPt_( iConfig.getParameter<double> ( "MinJetPt" ) ),
        maxJetEta_( iConfig.getParameter<double> ( "MaxJetEta" ) ),
        cTagTypeCvsL_( iConfig.getParameter<vector<std::string>>( "CTagTypeCvsL") ),
        cTagTypeCvsB_( iConfig.getParameter<vector<std::string>>( "CTagTypeCvsB") ),
        useJetID_( iConfig.getParameter<bool>   ( "UseJetID"     ) ),
        JetIDLevel_( iConfig.getParameter<string> ( "JetIDLevel"   ) ),
        cc_( consumesCollector() ),
        globalVariablesComputer_(iConfig.getParameter<edm::ParameterSet>("globalVariables"), cc_)
    {


      //  diPhotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
        inputDiPhotonName_= iConfig.getParameter<std::string > ( "DiPhotonName" );
        inputDiPhotonSuffixes_= iConfig.getParameter<std::vector<std::string> > ( "DiPhotonSuffixes" );
        std::vector<edm::InputTag>  diPhotonTags;
        for (auto & suffix : inputDiPhotonSuffixes_){ 
            systematicsLabels.push_back(suffix);
            std::string inputName = inputDiPhotonName_;
            inputName.append(suffix);
            if (!suffix.empty()) diPhotonTags.push_back(edm::InputTag(inputName));
            else  diPhotonTags.push_back(edm::InputTag(inputDiPhotonName_));
        }
        for( auto & tag : diPhotonTags ) { diPhotonTokens_.push_back( consumes<edm::View<flashgg::DiPhotonCandidate> >( tag ) ); }

        inputJetsName_= iConfig.getParameter<std::string> ( "JetsName" );
        inputJetsCollSize_= iConfig.getParameter<unsigned int> ( "JetsCollSize" );
        inputJetsSuffixes_= iConfig.getParameter<std::vector<std::string> > ( "JetsSuffixes" );
        std::vector<edm::InputTag>  jetTags;
        for (auto & suffix : inputJetsSuffixes_) {
            if (!suffix.empty()) systematicsLabels.push_back(suffix);  //nominal is already put in the diphoton loop
            for (unsigned int i = 0; i < inputJetsCollSize_ ; i++) {
                jetTags.push_back(edm::InputTag(inputJetsName_, std::to_string(i)));
            }         
        }
        for( auto & tag : jetTags ) { jetTokens_.push_back( consumes<edm::View<flashgg::Jet> >( tag ) ); }

        doPhotonId_ = iConfig.getUntrackedParameter<bool>("ApplyEGMPhotonID");        
        photonIDCut_ = iConfig.getParameter<double>("PhotonIDCut");

        photonElectronVeto_=iConfig.getUntrackedParameter<std::vector<int > >("PhotonElectronVeto");

        doSigmaMDecorr_ = iConfig.getUntrackedParameter<unsigned int>("DoSigmaMDecorr");
        if(doSigmaMDecorr_){
            sigmaMDecorrFile_ = iConfig.getUntrackedParameter<edm::FileInPath>("SigmaMDecorrFile");
            TFile* f_decorr = new TFile((sigmaMDecorrFile_.fullPath()).c_str(), "READ");
            TH2D* h_decorrEBEB_ = (TH2D*)f_decorr->Get("hist_sigmaM_M_EBEB"); 
            TH2D* h_decorrNotEBEB_ = (TH2D*)f_decorr->Get("hist_sigmaM_M_notEBEB");

            if(h_decorrEBEB_ && h_decorrNotEBEB_){
                transfEBEB_ = new DecorrTransform(h_decorrEBEB_ , 125., 1, 0);
                transfNotEBEB_ = new DecorrTransform(h_decorrNotEBEB_ , 125., 1, 0);

            } else {
                throw cms::Exception( "Configuration" ) << "The file "<<sigmaMDecorrFile_.fullPath()<<" provided for sigmaM/M decorrelation does not contain the expected histograms."<<std::endl;
            }
        }



       // produces<vector<cHTag>>();
        for (auto & systname : systematicsLabels) {
            produces<vector<cHTag>>(systname);
        }
        produces<vector<TagTruthBase>>();
    }


    void cHTagProducer::produce( Event &evt, const EventSetup & )
    {
        // update global variables
        globalVariablesComputer_.update(evt);

        // prepare output
        std::unique_ptr<vector<TagTruthBase> > truths( new vector<TagTruthBase> );
        edm::RefProd<vector<TagTruthBase> > rTagTruth = evt.getRefBeforePut<vector<TagTruthBase> >();

        // MC truth
        TagTruthBase truth_obj;
        if( ! evt.isRealData() ) {
            Handle<View<reco::GenParticle> > genParticles;
            std::vector<edm::Ptr<reco::GenParticle> > selHiggses;
            evt.getByToken( genParticleToken_, genParticles );
            Point higgsVtx(0.,0.,0.);
            for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
                int pdgid = genParticles->ptrAt( genLoop )->pdgId(); 
                if( pdgid == 25 || pdgid == 22 ) {
                    higgsVtx = genParticles->ptrAt( genLoop )->vertex();
                    break;
                }
            }
            truth_obj.setGenPV( higgsVtx );
            truths->push_back( truth_obj );
        }

      // read diphotons
      for (unsigned int diphoton_idx = 0; diphoton_idx < diPhotonTokens_.size(); diphoton_idx++) {//looping over all diphoton systematics
        Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
        evt.getByToken( diPhotonTokens_[diphoton_idx], diPhotons );
        
        unsigned int loopOverJets = 1;
        if (inputDiPhotonSuffixes_[diphoton_idx].empty()) loopOverJets = inputJetsSuffixes_.size();
        for (unsigned int jet_col_idx = 0; jet_col_idx < loopOverJets; jet_col_idx++) {//looping over all jet systematics, only for nominal diphotons
        std::unique_ptr<vector<cHTag> > tags( new vector<cHTag> );

        // loop over diphotons
        for( unsigned int candIndex = 0; candIndex < diPhotons->size() ; candIndex++ ) {
            edm::Ptr<flashgg::DiPhotonCandidate> dipho = diPhotons->ptrAt( candIndex );

            // kinematic cuts on diphotons
            auto leadPho = dipho->leadingPhoton();
            auto subleadPho = dipho->subLeadingPhoton();

            double leadPt = leadPho->pt();
            double subleadPt = subleadPho->pt();
            if( scalingPtCuts_ ) {
                leadPt /= dipho->mass();
                subleadPt /= dipho->mass();
            }
            if( leadPt <= minLeadPhoPt_ || subleadPt <= minSubleadPhoPt_ ) { continue; }
            //apply egm photon id with given working point
            if(doPhotonId_){
                if(leadPho->userFloat("EGMPhotonMVA")<photonIDCut_ || subleadPho->userFloat("EGMPhotonMVA")<photonIDCut_){
                    continue;
                }
            }
            //electron veto
            if(leadPho->passElectronVeto()<photonElectronVeto_[0] || subleadPho->passElectronVeto()<photonElectronVeto_[1]){
                continue;
            }


            // find vertex associated to diphoton object
            size_t vtx = (size_t)dipho->jetCollectionIndex();
            // and read corresponding jet collection
    

            edm::Handle<edm::View<flashgg::Jet> > jets;
            evt.getByToken( jetTokens_[jet_col_idx*inputJetsCollSize_+vtx], jets);  //take the corresponding vertex of current systematic

            // photon-jet cross-cleaning and pt/eta/btag/jetid cuts for jets
            std::vector<edm::Ptr<flashgg::Jet> > cleaned_jets;
            for( size_t ijet=0; ijet < jets->size(); ++ijet ) {//jets are ordered in pt
                auto jet = jets->ptrAt(ijet);
                if (jet->pt()<minJetPt_ || fabs(jet->eta())>maxJetEta_)continue;

                if( useJetID_ ){
                    if( JetIDLevel_ == "Loose" && !jet->passesJetID  ( flashgg::Loose ) ) continue;
                    if( JetIDLevel_ == "Tight" && !jet->passesJetID  ( flashgg::Tight ) ) continue;
                    if( JetIDLevel_ == "Tight2017" && !jet->passesJetID  ( flashgg::Tight2017 ) ) continue;
                    if( JetIDLevel_ == "Tight2018" && !jet->passesJetID  ( flashgg::Tight2018 ) ) continue;
                }
                if( reco::deltaR( *jet, *(dipho->leadingPhoton()) ) > vetoConeSize_ && reco::deltaR( *jet, *(dipho->subLeadingPhoton()) ) > vetoConeSize_ ) {
                    cleaned_jets.push_back( jet );
                }
            }
            int nCleanJets = cleaned_jets.size();
            if( nCleanJets < 1 ) { continue; }
            double ctagCvsL_ref = -999;
            double ctagCvsB_ref = -999;
            edm::Ptr<flashgg::Jet>  jet1, jet2;
            for( size_t ijet=0; ijet < cleaned_jets.size();++ijet){
                auto jet_1 = cleaned_jets[ijet];
                double ctagCvsL=0.;
                double ctagCvsB=0.;
                double ctag_probc = jet_1->bDiscriminator(cTagTypeCvsL_[0]); //prob c
                for (unsigned int ctag_num=0;ctag_num<cTagTypeCvsL_.size();ctag_num++)
                    ctagCvsL+=jet_1->bDiscriminator(cTagTypeCvsL_[ctag_num]); 
                for (unsigned int ctag_num=0;ctag_num<cTagTypeCvsB_.size();ctag_num++)
                    ctagCvsB+=jet_1->bDiscriminator(cTagTypeCvsB_[ctag_num]);
                if ((ctag_probc==-1000) && (ctagCvsL==-3000) && (ctagCvsB==-4000)) {
                    ctagCvsL = -1;
                    ctagCvsB = -1;
                } else { 
                    ctagCvsL = ctag_probc/ctagCvsL;
                    ctagCvsB = ctag_probc/ctagCvsB;
                }

                //if (ctagCvsL > ctagCvsL_ref) {
                  //  ctagCvsL_ref = ctagCvsL; //choosing as leading jet the one with the best separation against gluon and light
                if (ctagCvsB > ctagCvsB_ref) {
                    ctagCvsB_ref = ctagCvsB; //choosing as leading jet the one with the best separation against b
                    jet1 = jet_1;
                }
            }
            auto & leadJet = jet1; 
 
            // prepare tag object
            cHTag tag_obj( dipho, leadJet);
            tag_obj.setDiPhotonIndex( candIndex );
            if (loopOverJets == 1) 
                tag_obj.setSystLabel( inputDiPhotonSuffixes_[diphoton_idx] );
            else  
                tag_obj.setSystLabel( inputJetsSuffixes_[jet_col_idx]);

            
            if(doSigmaMDecorr_){
                tag_obj.setSigmaMDecorrTransf(transfEBEB_,transfNotEBEB_);
            }


            tag_obj.setEventNumber(evt.id().event() );
            tag_obj.nJets_ = nCleanJets;
           
            tag_obj.setCategoryNumber( 0 );
            tag_obj.includeWeights( *dipho );
          //  tag_obj.includeWeightsByLabel( *leadJet ,"JetBTagReshapeWeight");



             tags->push_back( tag_obj );
             // link mc-truth
            if( ! evt.isRealData() ) {
                  tags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );
            }                 
        }
        if (loopOverJets == 1) 
            evt.put( std::move( tags ),inputDiPhotonSuffixes_[diphoton_idx] );
        else  
            evt.put( std::move( tags ),inputJetsSuffixes_[jet_col_idx] );
        }
        }   
        evt.put( std::move( truths ) );
    }
    
    
    
}

typedef flashgg::cHTagProducer FlashggcHTagProducer;
DEFINE_FWK_MODULE( FlashggcHTagProducer );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
