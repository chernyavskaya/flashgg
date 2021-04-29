#include "flashgg/DataFormats/interface/cHTag.h"

using namespace flashgg;

cHTag::cHTag() : DiPhotonTagBase::DiPhotonTagBase() 
{  
}

cHTag::~cHTag() {}

cHTag::cHTag( edm::Ptr<flashgg::DiPhotonCandidate> diPho, edm::Ptr<flashgg::Jet> leadJet )
    : leadJet_(leadJet)
{
    dipho_ = diPho;
}

cHTag *cHTag::clone() const
{
    cHTag *result = new cHTag( *this );
    return result;
}


float cHTag::getPhoJetMinDr() const
{
    float PhoJetMinDr = min(deltaR( diPhoton()->leadingPhoton()->p4(), leadJet().p4() ), deltaR( diPhoton()->subLeadingPhoton()->p4(), leadJet().p4()  ));
    
    return PhoJetMinDr;
}



float cHTag::getSigmaMDecorr() const
{
    double mass_sigma[2]={0.,0.};
    double dummy[1]={0.};
    mass_sigma[0]=diPhoton()->mass();
    mass_sigma[1] = 0.5*sqrt((diPhoton()->leadingPhoton()->sigEOverE()*diPhoton()->leadingPhoton()->sigEOverE() + diPhoton()->subLeadingPhoton()->sigEOverE()*diPhoton()->subLeadingPhoton()->sigEOverE()));
    float sigmaMOverMDecorr=-99;
    //Splitting EBEB and !EBEB, using cuts as in preselection
    if(abs(diPhoton()->leadingPhoton()->superCluster()->eta())<1.4442 && abs(diPhoton()->subLeadingPhoton()->superCluster()->eta())<1.4442){
        sigmaMOverMDecorr = (*transfEBEB_)(mass_sigma,dummy);
    }
    else{
        sigmaMOverMDecorr = (*transfNotEBEB_)(mass_sigma,dummy);
    }
    return sigmaMOverMDecorr;
}


// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

