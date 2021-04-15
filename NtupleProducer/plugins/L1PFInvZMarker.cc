#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <vector>

class L1PFInvZMarker : public edm::global::EDProducer<> {
public:
    explicit L1PFInvZMarker(const edm::ParameterSet&);
    ~L1PFInvZMarker() override;

private:
    void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
    
    edm::EDGetTokenT<std::vector<l1t::PFCandidate>> l1PFCandidates_;
    //edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_;
    edm::EDGetTokenT<edm::View<reco::Candidate> > genParticles_;
    edm::EDGetTokenT<edm::View<reco::Candidate> > genJets_;
    std::string genOut_="genParticlesInvZ";
    std::string genJetOut_="genJetsInvZ";
    std::string pfOut_="L1PFCandidatesInvZ";
    //std::string pfOut_="l1pfCandidates:PuppiInvZ";
    bool debug=false;
};

L1PFInvZMarker::L1PFInvZMarker(const edm::ParameterSet& cfg)
    : l1PFCandidates_(consumes<std::vector<l1t::PFCandidate>>(cfg.getParameter<edm::InputTag>("L1PFObjects"))),
      // genParticles_(consumes<std::vector<reco::GenParticle>>(cfg.getParameter<edm::InputTag>("genParticles"))) {
      genParticles_(consumes<edm::View<reco::Candidate> >(cfg.getParameter<edm::InputTag>("genParticles"))),
      genJets_(consumes<edm::View<reco::Candidate> >(cfg.getParameter<edm::InputTag>("genJets"))) {

    produces<std::vector<l1t::PFCandidate> >( pfOut_ );
    //produces<std::vector<reco::GenParticle> >( genOut_ );
    produces<reco::CandidatePtrVector>(genOut_);
    produces<reco::CandidatePtrVector>(genJetOut_);
}

L1PFInvZMarker::~L1PFInvZMarker() {}

void L1PFInvZMarker::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const {

    edm::Handle<std::vector<l1t::PFCandidate>> l1PFCandidates;
    //edm::Handle<std::vector<reco::GenParticle>> genParticles;
    edm::Handle<edm::View<reco::Candidate>> genParticles;
    edm::Handle<edm::View<reco::Candidate>> genJets;
    iEvent.getByToken(l1PFCandidates_, l1PFCandidates);
    iEvent.getByToken(genParticles_, genParticles);
    iEvent.getByToken(genJets_, genJets);

    std::unique_ptr<std::vector<l1t::PFCandidate> > outPF(new std::vector<l1t::PFCandidate>(0));
    //std::unique_ptr<std::vector<reco::GenParticle> > outGen(new std::vector<reco::GenParticle>(0));
    std::unique_ptr<reco::CandidatePtrVector> outGen(new reco::CandidatePtrVector);
    std::unique_ptr<reco::CandidatePtrVector> outGenJet(new reco::CandidatePtrVector);

    std::set<unsigned int> to_remove{};
    int pdgId, absPdgId;

    for(int ip=0; ip < int(genParticles->size()); ip++){
        //const auto& p = genParticles->at(ip);
        const auto& cand = genParticles->at(ip);
        const auto& p = *dynamic_cast<const reco::GenParticle*>( &cand );

        pdgId = p.pdgId();
        absPdgId = abs(pdgId);

        if(absPdgId==11 || absPdgId==13){
            // lep from tau decay is 1<<10, but we don't consider taus here
            if( p.statusFlags().isLastCopy() && (p.statusFlags().isHardProcessTauDecayProduct() || (p.statusFlags().isPrompt() && p.statusFlags().fromHardProcess()))){
                if(debug) std::cout << "found " << pdgId << ", pt = " << p.pt() << ", " << p.statusFlags().flags_ << endl;
                to_remove.insert(ip);
                // if( p.statusFlags().isHardProcessTauDecayProduct() ){
                //     if(debug) cout << "  #M " << p.numberOfMothers() << endl;
                //     for(int i=0;i<int(p.numberOfMothers());i++){
                //         if(debug) cout << "  m: " << p.mother(i)->pdgId() << " " << p.mother(i)->pt() << endl;
                //     }
                // }
            }
            // NB the visible genparticles are not killed for taus in this case... need to veto these evts or improve this
        } else if(absPdgId==15) {
            if( p.statusFlags().isLastCopy() && p.statusFlags().isPrompt() && p.statusFlags().fromHardProcess()){
                if(debug) std::cout << "found " << pdgId << ", pt = " << p.pt() << ", " << p.statusFlags().flags_ << endl;
                to_remove.insert(ip);
            }            
        }
        if (!to_remove.count(ip)){
            //outGen->push_back(p);
            outGen->push_back( genParticles->ptrAt(ip) );
        }
    }

    // check if pfCands are due to the Z decay and should be removed
    bool remove;
    for(int ip=0; ip < int(l1PFCandidates->size()); ip++){
        const auto& p = l1PFCandidates->at(ip);
        remove=false;

        for(unsigned int igen : to_remove){
            const auto& g = genParticles->at(igen);
            if( reco::deltaR(p, g) < ( abs(g.pdgId())==15 ? 0.3 : 0.1) ){
                remove=true;
                if(debug) cout << "removed " << p.id() << "  with pt " << p.pt() << endl;
            }
        }
        if(!remove) outPF->push_back(p);
    }

    // check if genJets are due to the Z decay and should be removed (instead of rebuilding them)
    for(int ijet=0; ijet < int(genJets->size()); ijet++){
        const auto& j = genJets->at(ijet);
        remove=false;

        for(unsigned int igen : to_remove){
            const auto& g = genParticles->at(igen);
            if( reco::deltaR(j, g) < 0.4 ){
                remove=true;
                if(debug) cout << "removed jet with pt " << j.pt() << endl;
            }
        }
        if(!remove) outGenJet->push_back( genJets->ptrAt(ijet) );
    }

    iEvent.put(std::move(outPF), pfOut_);
    iEvent.put(std::move(outGen), genOut_);
    iEvent.put(std::move(outGenJet), genJetOut_);


}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1PFInvZMarker);
