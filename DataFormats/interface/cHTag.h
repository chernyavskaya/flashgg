#ifndef flashgg_cHTag
#define flashgg_cHTag

#include "TLorentzVector.h"

#include "flashgg/DataFormats/interface/DiPhotonTagBase.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "flashgg/Taggers/interface/FunctionHelpers.h"

namespace flashgg {

    class cHTag: public DiPhotonTagBase, public reco::LeafCandidate
    {
    public:
        cHTag();
        ~cHTag();

        cHTag( edm::Ptr<DiPhotonCandidate>, edm::Ptr<flashgg::Jet> );
        virtual cHTag *clone() const override;

        double diphotonPtOverM() const {return diPhoton()->pt()/mass(); }
        const flashgg::Jet & leadJet() const { return *leadJet_; } 
        //const flashgg::Jet & subleadJet() const { return *subleadJet_; } 
        
        float nJets() const { return nJets_; }

        float getPhoJetMinDr() const;
        float getSigmaMDecorr() const;
        void  setSigmaMDecorrTransf( DecorrTransform* transfEBEB, DecorrTransform* transfNotEBEB){ transfEBEB_= transfEBEB; transfNotEBEB_=transfNotEBEB;}
        void setEventNumber(double x) { eventNumber_ = x; }
        double eventNumber() const { return eventNumber_; }

        int nJets_; 

    private:
         long eventNumber_;
        edm::Ptr<flashgg::Jet> leadJet_;
     //   edm::Ptr<flashgg::Jet> subleadJet_;
        DecorrTransform* transfEBEB_;
        DecorrTransform* transfNotEBEB_;
        
    };
}

#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

