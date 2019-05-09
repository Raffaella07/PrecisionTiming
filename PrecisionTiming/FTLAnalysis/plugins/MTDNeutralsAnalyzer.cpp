
#ifndef _MTD_NEUTRALS_ANALIZER_
#define _MTD_NEUTRALS_ANALIZER_

#include "TMath.h"
#include "Math/PositionVector3D.h"


#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "FastSimulation/Particle/interface/RawParticle.h"
#include "FWCore/Utilities/interface/BranchType.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/Provenance.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHit.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/RefToBase.h"
//#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementError.h"
//#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Geometry/CommonTopologies/interface/Topology.h"
#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeomDetUnit.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"

#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruthFinder.h"
#include "RecoEgamma/EgammaMCTools/interface/PhotonMCTruth.h"

#include "RecoMTD/DetLayers/interface/MTDDetLayerGeometry.h"
#include "RecoMTD/DetLayers/interface/MTDTrayBarrelLayer.h"
#include "RecoMTD/DetLayers/interface/MTDDetTray.h"
#include "RecoMTD/DetLayers/interface/MTDRingForwardDoubleLayer.h"
#include "RecoMTD/DetLayers/interface/MTDDetRing.h"
#include "RecoMTD/Records/interface/MTDRecoGeometryRecord.h"

#include "PrecisionTiming/FTLAnalysis/interface/MTDNeutralsTree.h"
#include "PrecisionTiming/FTLAnalysis/interface/FTLTracksTree.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "Calibration/IsolatedParticles/interface/CaloPropagateTrack.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "FWCore/Framework/interface/data_default_record_trait.h"
#include "HepPDT/ParticleDataTable.hh"
#include "SimGeneral/HepPDTRecord/interface/PDTRecord.h"
using namespace std;                             

using namespace std;

class MTDHitMatchingInfo
{
public:
    MTDHitMatchingInfo()
        {
            hit = -1;
            estChi2 = std::numeric_limits<double>::max();
            timeChi2 = std::numeric_limits<double>::max();
        }

    MTDHitMatchingInfo(int this_hit, double this_estChi2, double this_timeChi2) :
        hit(this_hit),
        estChi2(this_estChi2),
        timeChi2(this_timeChi2)
        { }
    
    MTDHitMatchingInfo(const MTDHitMatchingInfo& a) :
        hit(a.hit),
        estChi2(a.estChi2),
        timeChi2(a.timeChi2)
        { }
    //Operator used to sort the hits while performing the matching step at the MTD
    inline bool operator<(const MTDHitMatchingInfo &m2) const {
	    //only for good matching in time use estChi2, otherwise use mostly time compatibility
	    if (timeChi2<5 && m2.timeChi2<5)
		    return chi2(5.) < m2.chi2(5.);
	    else
		    return chi2(10.) < m2.chi2(10.);
    }

    inline double chi2(float timeWeight=1.) const { return estChi2 + timeWeight*timeChi2; }

    int hit;
    double estChi2;
    double timeChi2;
};

class MTDNeutralsAnalyzer : public edm::EDAnalyzer
{
	public:                             

		typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> genXYZ;
		typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag> Point;

		explicit MTDNeutralsAnalyzer(const edm::ParameterSet& pSet);
		~MTDNeutralsAnalyzer() {};

		//---methods
		virtual void beginJob() override {};
		virtual void analyze(edm::Event const&, edm::EventSetup const&) override;
		virtual void endJob() override {};

		//--utils

		pair<float, float> getNeutralsMTDMatchingChi2s(reco::PFCandidate& cand, GlobalPoint& mtd_gp, const math::XYZTLorentzVectorD& genPV, float mtd_time); 
	private:
		//---inputs
		const MTDGeometry* mtdGeometry_;
		//--magnetic field 

		//---gen particles
		edm::Handle<reco::GenParticleCollection> genParticlesHandle_;
		edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;

		//---sim hit tracker, for EGamma mc truth tool
		edm::Handle<edm::SimTrackContainer> simTkHandle_;
		edm::EDGetTokenT<edm::SimTrackContainer> simTkToken_;
		edm::Handle<edm::SimVertexContainer> simVtxHandle_;
		edm::EDGetTokenT<edm::SimVertexContainer> simVtxToken_;


		//---sim MTD hits
		edm::Handle<std::vector<PSimHit> > simHitsBTLHandle_;
		edm::EDGetTokenT<std::vector<PSimHit> > simHitsBTLToken_;    
		edm::Handle<std::vector<PSimHit> > simHitsETLHandle_;
		edm::EDGetTokenT<std::vector<PSimHit> > simHitsETLToken_;    

		//---BTL clusters
		edm::Handle<FTLRecHitCollection> recHitsBTLHandle_;
		edm::EDGetTokenT<FTLRecHitCollection> recHitsBTLToken_;    
		edm::Handle<FTLClusterCollection> clustersBTLHandle_;
		edm::EDGetTokenT<FTLClusterCollection> clustersBTLToken_;    
		//---ETL clusters
		edm::Handle<FTLRecHitCollection> recHitsETLHandle_;
		edm::EDGetTokenT<FTLRecHitCollection> recHitsETLToken_;    
		edm::Handle<FTLClusterCollection> clustersETLHandle_;
		edm::EDGetTokenT<FTLClusterCollection> clustersETLToken_;    


		//electronstracks---Tracks
		edm::EDGetTokenT<edm::View<reco::Track> > tracksToken_;
		edm::Handle<edm::View<reco::Track> > tracksHandle_;
		edm::EDGetTokenT<edm::View<reco::Track> > extTracksToken_;
		edm::Handle<edm::View<reco::Track> > extTracksHandle_;

		edm::EDGetTokenT<edm::ValueMap<int> > generalToExtendedTrkMapToken_;
		edm::Handle<edm::ValueMap<int> > generalToExtendedTrkMapHandle_;
		edm::EDGetTokenT<edm::ValueMap<float> > extTracksMTDtimeToken_;
		edm::Handle<edm::ValueMap<float> > extTracksMTDtimeHandle_;


		//---PF
		edm::Handle<reco::PFCandidateCollection> pfCandidatesHandle_;
		edm::EDGetTokenT<reco::PFCandidateCollection> pfCandidatesToken_;


		//---PF clusters
		edm::Handle<reco::PFClusterCollection> pfClustersHandle_;
		edm::EDGetTokenT<reco::PFClusterCollection> pfClustersToken_;

		//---vertices
		edm::EDGetTokenT<genXYZ>                  genXYZToken_;
		edm::Handle<genXYZ>                       genXYZHandle_;
		edm::EDGetTokenT<float>                   genT0Token_;
		edm::Handle<float>                        genT0Handle_;



		edm::EDGetTokenT<vector<reco::Vertex> >   vtx3DToken_;
		edm::Handle<vector<reco::Vertex> >        vtx3DHandle_;
		edm::EDGetTokenT<vector<reco::Vertex> >   vtx4DToken_;
		edm::Handle<vector<reco::Vertex> >        vtx4DHandle_;    


		//---options
		BTLDetId::CrysLayout crysLayout_;

		//---workers
		PhotonMCTruthFinder photonMCTruthFinder_;


		//	edm::EDGetTokenT<vector<SimVertex> >                 genVtxToken_;
		//	edm::Handle<vector<SimVertex> >                      genVtxHandle_;    
		//---outputs
		MTDNeutralsTree genTree_;
		MTDNeutralsTree candTree_;
		MTDNeutralsTree mctTree_;
		MTDNeutralsTree outTreet_;
		MTDNeutralsTree dumbTree_;
		edm::Service<TFileService> fs_;  

};

MTDNeutralsAnalyzer::MTDNeutralsAnalyzer(const edm::ParameterSet& pSet):
	genParticlesToken_(consumes<reco::GenParticleCollection>(pSet.getUntrackedParameter<edm::InputTag>("genParticlesTag"))),
	simTkToken_(consumes<edm::SimTrackContainer>(pSet.getUntrackedParameter<edm::InputTag>("simTkTag"))),
	simVtxToken_(consumes<edm::SimVertexContainer>(pSet.getUntrackedParameter<edm::InputTag>("simVtxTag"))),
	simHitsBTLToken_(consumes<std::vector<PSimHit> >(pSet.getUntrackedParameter<edm::InputTag>("simHitsBTLTag"))), 
	simHitsETLToken_(consumes<std::vector<PSimHit> >(pSet.getUntrackedParameter<edm::InputTag>("simHitsETLTag"))),   
	clustersBTLToken_(consumes<FTLClusterCollection>(pSet.getUntrackedParameter<edm::InputTag>("clustersBTLTag"))),    
	clustersETLToken_(consumes<FTLClusterCollection>(pSet.getUntrackedParameter<edm::InputTag>("clustersETLTag"))),
	//electrontracks--comment if uninterested in misidentified early conversion electrons
	tracksToken_(consumes<edm::View<reco::Track> >(pSet.getUntrackedParameter<edm::InputTag>("generalTracksTag"))),
	extTracksToken_(consumes<edm::View<reco::Track> >(pSet.getUntrackedParameter<edm::InputTag>("extendedTracksTag"))),
//genVtxToken_(consumes<vector<SimVertex> >(pSet.getUntrackedParameter<edm::InputTag>("genVtxTag"))),    
	generalToExtendedTrkMapToken_(consumes<edm::ValueMap<int> >(pSet.getUntrackedParameter<edm::InputTag>("generalToExtendedTrkMapTag"))),    
   	extTracksMTDtimeToken_(consumes<edm::ValueMap<float> >(pSet.getUntrackedParameter<edm::InputTag>("extTracksMTDtimeTag"))),	
	pfCandidatesToken_(consumes<reco::PFCandidateCollection>(pSet.getUntrackedParameter<edm::InputTag>("pfCandidatesTag"))),    
	pfClustersToken_(consumes<reco::PFClusterCollection>(pSet.getUntrackedParameter<edm::InputTag>("pfClusterTag"))),    
	genXYZToken_(consumes<genXYZ>(pSet.getUntrackedParameter<edm::InputTag>("genXYZTag"))),
	genT0Token_(consumes<float>(pSet.getUntrackedParameter<edm::InputTag>("genT0Tag"))),
	vtx3DToken_(consumes<vector<reco::Vertex> >(pSet.getUntrackedParameter<edm::InputTag>("vtx3DTag"))),
	vtx4DToken_(consumes<vector<reco::Vertex> >(pSet.getUntrackedParameter<edm::InputTag>("vtx4DTag"))),
	crysLayout_((BTLDetId::CrysLayout)(pSet.getUntrackedParameter<int>("crysLayout"))),
	photonMCTruthFinder_()
{
	genTree_ = MTDNeutralsTree(pSet.getUntrackedParameter<string>("outTreeName1").c_str(), "4D TOFPID studies");
	mctTree_ = MTDNeutralsTree(pSet.getUntrackedParameter<string>("outTreeName4").c_str(), "4D TOFPID studies");
	candTree_ = MTDNeutralsTree(pSet.getUntrackedParameter<string>("outTreeName3").c_str(), "4D TOFPID studies");
	outTreet_ = MTDNeutralsTree(pSet.getUntrackedParameter<string>("outTreeName2").c_str(), "4D TOFPID studies");
	dumbTree_ = MTDNeutralsTree(pSet.getUntrackedParameter<string>("outTreeName5").c_str(), "4D TOFPID studies");
}

void MTDNeutralsAnalyzer::analyze(edm::Event const& event, edm::EventSetup const& setup)
{

	mctTree_.Reset();
	candTree_.Reset();
	genTree_.Reset();
	outTreet_.Reset();
	dumbTree_.Reset();

	//--load magnetic field
  // get Geometry, B-field, Topology
  	   edm::ESHandle<MagneticField> magfield;
 	   setup.get<IdealMagneticFieldRecord>().get(magfield);
 	   const MagneticField* mag = magfield.product();

	//---load gen particles
	event.getByToken(genParticlesToken_, genParticlesHandle_);
	auto genParticles = *genParticlesHandle_.product();

	//---gen mc-truth photons
	event.getByToken(simTkToken_, simTkHandle_);
	event.getByToken(simVtxToken_, simVtxHandle_);
	auto mcTruthPhotons = photonMCTruthFinder_.find(*simTkHandle_.product(), *simVtxHandle_.product());

	//---load sim hits
	event.getByToken(simHitsBTLToken_, simHitsBTLHandle_);
	auto simHitsBTL = *simHitsBTLHandle_.product();
	event.getByToken(simHitsETLToken_, simHitsETLHandle_);
	auto simHitsETL = *simHitsETLHandle_.product();


	//---get the MTD geometry
	edm::ESHandle<MTDGeometry> geoHandle;
	setup.get<MTDDigiGeometryRecord>().get(geoHandle);
	mtdGeometry_ = geoHandle.product();

	edm::ESHandle<MTDDetLayerGeometry> layerGeo;
	setup.get<MTDRecoGeometryRecord>().get(layerGeo);

	//---MTD clusters
	event.getByToken(clustersBTLToken_, clustersBTLHandle_);
	auto clustersBTL = *clustersBTLHandle_.product();
	event.getByToken(clustersETLToken_, clustersETLHandle_);
	auto clustersETL = *clustersETLHandle_.product();

	//---load PFCandidates
	event.getByToken(pfCandidatesToken_, pfCandidatesHandle_);
	auto pfCandidates = *pfCandidatesHandle_.product();    

	//---load PFClusters
	event.getByToken(pfClustersToken_, pfClustersHandle_);
	auto pfClusters = *pfClustersHandle_.product();    
	
	//---load general tracks
	event.getByToken(tracksToken_,tracksHandle_);
	auto tracks = *tracksHandle_.product();

	//---load extended tracks
	event.getByToken(extTracksToken_, extTracksHandle_);

	auto extTracks = *extTracksHandle_.product();
	
	//---load general to extended tracks map
	event.getByToken(generalToExtendedTrkMapToken_, generalToExtendedTrkMapHandle_);
	auto generalToExtendedTrkMap = *generalToExtendedTrkMapHandle_.product();

	//---load MTD time 
	event.getByToken(extTracksMTDtimeToken_, extTracksMTDtimeHandle_);
    	auto extTracksMTDtime = *extTracksMTDtimeHandle_.product();	

	//---load gen, sim and reco vertices
	// GEN
	event.getByToken(genXYZToken_, genXYZHandle_);
	event.getByToken(genT0Token_, genT0Handle_);
	auto xyz = genXYZHandle_.product();
	auto t = *genT0Handle_.product();
	auto v = math::XYZVectorD(xyz->x(), xyz->y(), xyz->z());
	auto genPV = SimVertex(v, t).position();
	// 3D
	event.getByToken(vtx3DToken_, vtx3DHandle_);
	auto vtxs3D = *vtx3DHandle_.product();
	// Full 4D
	event.getByToken(vtx4DToken_, vtx4DHandle_);
	auto vtxs4D = *vtx4DHandle_.product(); 

	/*	int counter =0;
		for(auto& gen_part : genParticles){

		std::cout << "Event  " << event.id().event() << "GenParticle  " << counter << "   : pdgId   " << gen_part.pdgId() << "    pt   " << gen_part.pt() << "   eta  " << gen_part.eta() << "   phi   " << gen_part.phi() << "  status   " << gen_part.status() << std::endl;
		counter++;	
		}*/

	//	outTree_.genpv_t->resize(30);
	//	int idx=0;

	//
	//	outTree_.genpv->pu

	//---fill global info
	outTreet_.run = event.id().run();
	outTreet_.event = event.id().event();
	outTreet_.lumi = event.id().luminosityBlock();
	genTree_.run = event.id().run();            
	genTree_.lumi = event.id().luminosityBlock();
	genTree_.event = event.id().event();            
	mctTree_.event = event.id().event();
	mctTree_.lumi = event.id().luminosityBlock();
	mctTree_.run = event.id().run();
	candTree_.event = event.id().event();
	candTree_.lumi = event.id().luminosityBlock();
	candTree_.run = event.id().run();
	dumbTree_.event = event.id().event();
	dumbTree_.lumi = event.id().luminosityBlock();
	dumbTree_.run = event.id().run();
	int l_counter=0;    
	int gen_counter=0; 

	for(auto& gen_part: genParticles)//loop on genparticles
	{
		genTree_.pt->push_back(gen_part.pt()); 
		genTree_.eta ->push_back( gen_part.eta());                
		genTree_.phi ->push_back( gen_part.phi());            
		genTree_.pdgId ->push_back( gen_part.pdgId());
		genTree_.label ->push_back(gen_counter);
		std::cout << "conter " << gen_counter << std::endl;
		gen_counter++;
		int sim_match=0;	
	//---sim hits
//		float min_time = 1e9;
		float minDR_simh_gen = 1e6;
		GlobalPoint matched_hit;
		
		float dr = 1e6;
		for(auto& simHit : simHitsBTL)
		{
	
			if(simHit.tof()<0 || simHit.tof()>25)
				continue;

			BTLDetId id = simHit.detUnitId();
			DetId geoId = id.geographicalId( crysLayout_ );
			const auto& det = mtdGeometry_ -> idToDet(geoId);
			const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
			const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());

			LocalPoint lp_entry(simHit.entryPoint().x()/10., simHit.entryPoint().y()/10., simHit.entryPoint().z()/10.);
			GlobalPoint gp_entry = det->toGlobal(topo.pixelToModuleLocalPoint(lp_entry,id.row(topo.nrows()),id.column(topo.nrows())));
			float dr_gen = deltaR(gp_entry.eta(), gp_entry.phi(), gen_part.eta(), gen_part.phi());
			if(dr_gen<dr){
				dr= dr_gen;
				matched_hit = gp_entry; 
				

			}
		}
		if(dr < 0.01){	
			candTree_.MTDx->push_back(matched_hit.x());
			candTree_.MTDy->push_back(matched_hit.y());
			candTree_.MTDz->push_back(matched_hit.z());
			candTree_.MTDphi->push_back(matched_hit.phi());
			candTree_.MTDeta->push_back(matched_hit.eta());
			candTree_.MTDdr->push_back(dr);
			sim_match++;
			}
/*		//---match with gen candidate
			if(dr_gen < minDR_simh_gen)
			{
				minDR_simh_gen = dr_gen;
				outTree_.simh_energy = simHit.energyLoss()*1000.;
				outTree_.simh_time = simHit.tof();
				outTree_.simh_tof = sqrt(pow(gp_entry.x()-genPV.x(), 2)+pow(gp_entry.y()-genPV.y(), 2)+pow(gp_entry.z()-genPV.z(), 2))/2.99792458e1;
				outTree_.simh_DR = dr_gen;
				if(min_chi2_pos > -1)// < outTree_.clus_neu_DR->size())
				{
					auto r_eta = outTree_.clus_eta->at(min_chi2_pos);
					auto r_phi = outTree_.clus_phi->at(min_chi2_pos);
					outTree_.simh_recoh_DR = deltaR(gp_entry.eta(), gp_entry.phi(), r_eta, r_phi);
					outTree_.simh_recoh_DPhi = deltaPhi(gp_entry.phi(), r_phi);                                
				}
			}*/



				reco::PFCluster matched_clus;
				 dr=1e6;
				
				for(auto& clus:pfClusters ){//loop on pf candidates
				float dr_clus = deltaR(gen_part.eta(), gen_part.phi(), clus.eta(), clus.phi());
				if(dr_clus < dr){
					dr = dr_clus;
					matched_clus = clus;


					}
				}
				if(dr< 0.3){
				candTree_.ecalx ->push_back (matched_clus.x());
				candTree_.ecaly ->push_back (matched_clus.y());
				candTree_.ecalz ->push_back (matched_clus.z());
				candTree_.ecalphi ->push_back (matched_clus.phi());
				candTree_.ecaleta ->push_back (matched_clus.eta());
				candTree_.dr ->push_back (dr);
			
				
				}

					
				
 					



		PhotonMCTruth* mc_truth_pho=NULL;
	

		float min_dr_mct = 1e6;
		for(auto& mcpho : mcTruthPhotons)//loop on photons: searches for the nearest photon to each gen particle
		{
			if(mcpho.fourMomentum().et() > 0.5  && (gen_part.pdgId() == 22|| gen_part.pdgId() == 111) && deltaR(mcpho.fourMomentum().eta(), mcpho.fourMomentum().phi(), gen_part.eta(), gen_part.phi()) < min_dr_mct)
			{
				min_dr_mct = deltaR(mcpho.fourMomentum().eta(), mcpho.fourMomentum().phi(), gen_part.eta(), gen_part.phi());
				if(min_dr_mct < 0.3){ 	
					mc_truth_pho = &mcpho;
					

			}


		}
		}
		if(mc_truth_pho){
					std::cout << "vertex z" << mc_truth_pho->vertex().z() << endl;
			
					genTree_.convRadius->push_back(mc_truth_pho->vertex().vect().rho());
				//	genTree_.sim_match->push_back(sim_match);
			if(mc_truth_pho->electrons().size()==0 || mc_truth_pho->vertex().vect().rho()>116){
				
					mctTree_.dr ->push_back(min_dr_mct);
					mctTree_.PId->push_back(4);
					mctTree_.label->push_back(l_counter);
					mctTree_.primary_x ->push_back( mc_truth_pho->primaryVertex().x() );
					mctTree_.primary_y ->push_back( mc_truth_pho->primaryVertex().y() );
					mctTree_.primary_z ->push_back( mc_truth_pho->primaryVertex().z() );
					mctTree_.primary_phi ->push_back( mc_truth_pho->primaryVertex().phi());
					mctTree_.primary_eta ->push_back( mc_truth_pho->primaryVertex().eta());
					mctTree_.convRadius->push_back(mc_truth_pho->vertex().vect().rho());
					mctTree_.convPhi->push_back(mc_truth_pho->vertex().vect().phi());
					mctTree_.convZ->push_back(mc_truth_pho->vertex().vect().z());            
					mctTree_.conv_x->push_back(mc_truth_pho->vertex().vect().x());
					mctTree_.conv_y->push_back(mc_truth_pho->vertex().vect().y());
					mctTree_.conv_z->push_back(mc_truth_pho->vertex().vect().z());
					mctTree_.conv_phi->push_back(mc_truth_pho->vertex().vect().phi());
					mctTree_.conv_eta->push_back(mc_truth_pho->vertex().vect().eta());
					mctTree_.energy->push_back(mc_truth_pho->fourMomentum().e());
					mctTree_.eta->push_back(mc_truth_pho->fourMomentum().eta());
					mctTree_.phi->push_back(mc_truth_pho->fourMomentum().phi());
					mctTree_.pt->push_back(mc_truth_pho->fourMomentum().et());
					auto conversion_legs = mc_truth_pho->electrons();
					mctTree_.nlegs->push_back(conversion_legs.size());
					mctTree_.plabel->push_back(gen_counter);
					CLHEP::HepLorentzVector mom;
					mom.setREtaPhi(130,mc_truth_pho->fourMomentum().eta(),mc_truth_pho->fourMomentum().phi());
					dumbTree_.PId->push_back(4);
					dumbTree_.ecalx->push_back(mom.x());
					dumbTree_.ecaly->push_back(mom.y());
					dumbTree_.ecalz->push_back(mom.z());
					dumbTree_.ecalphi->push_back(mom.phi());
					dumbTree_.ecaleta->push_back(mom.eta());
					mom.setREtaPhi(117,mc_truth_pho->fourMomentum().eta(),mc_truth_pho->fourMomentum().phi());
					dumbTree_.MTDx->push_back(mom.x());
					dumbTree_.MTDy->push_back(mom.y());
					dumbTree_.MTDz->push_back(mom.z());
					dumbTree_.MTDphi->push_back(mom.phi());
					dumbTree_.MTDeta->push_back(mom.eta());
					
					std::cout << "#######____" << mom.z() << std::endl;

	l_counter++;	
}else if(mc_truth_pho->electrons().size()==2 ){
					mctTree_.label->push_back(l_counter);
					mctTree_.plabel->push_back(gen_counter);
					mctTree_.dr ->push_back(min_dr_mct);
					mctTree_.PId->push_back(2);
					mctTree_.primary_x ->push_back( mc_truth_pho->primaryVertex().x() );
					mctTree_.primary_y ->push_back( mc_truth_pho->primaryVertex().y() );
					mctTree_.primary_z ->push_back( mc_truth_pho->primaryVertex().z() );
					mctTree_.primary_phi ->push_back( mc_truth_pho->primaryVertex().phi());
					mctTree_.primary_eta ->push_back( mc_truth_pho->primaryVertex().eta());
					mctTree_.convRadius->push_back(mc_truth_pho->vertex().vect().rho());
					mctTree_.convPhi->push_back(mc_truth_pho->vertex().vect().phi());
					mctTree_.convZ->push_back(mc_truth_pho->vertex().vect().z());            
					mctTree_.conv_x->push_back(mc_truth_pho->vertex().vect().x());
					mctTree_.conv_y->push_back(mc_truth_pho->vertex().vect().y());
					mctTree_.conv_z->push_back(mc_truth_pho->vertex().vect().z());
					mctTree_.conv_phi->push_back(mc_truth_pho->vertex().vect().phi());
					mctTree_.conv_eta->push_back(mc_truth_pho->vertex().vect().eta());
					mctTree_.energy->push_back(mc_truth_pho->fourMomentum().e());
					auto conversion_legs = mc_truth_pho->electrons();
					mctTree_.nlegs->push_back(conversion_legs.size());
					RawParticle part=RawParticle(conversion_legs[0].fourMomentum().x(),conversion_legs[0].fourMomentum().y(),conversion_legs[0].fourMomentum().z(),conversion_legs[0].fourMomentum().e());
			
					part.setVertex(conversion_legs[0].primaryVertex().x(),conversion_legs[0].primaryVertex().y(),conversion_legs[0].primaryVertex().z(),conversion_legs[0].primaryVertex().t());
					part.setCharge(conversion_legs[0].simTracks().charge());
				
					BaseParticlePropagator p=BaseParticlePropagator(part,0,0,0);
					GlobalPoint pos(conversion_legs[0].primaryVertex().x(),conversion_legs[0].primaryVertex().y(),conversion_legs[0].primaryVertex().z());
				//	bool prop;

					p.setMagneticField(mag->inTesla(pos).z());
					p.SetPx(conversion_legs[0].fourMomentum().x());
					p.SetPy(conversion_legs[0].fourMomentum().y());
					p.SetPz(conversion_legs[0].fourMomentum().z());
				//	prop = p.propagateToEcalEntrance();
					std::cout << mag->inTesla(pos).z() << std::endl;
								mctTree_.eta->push_back(conversion_legs[0].fourMomentum().eta());
								mctTree_.phi->push_back(conversion_legs[0].fourMomentum().phi());
								mctTree_.pt->push_back(conversion_legs[0].fourMomentum().et());
	//	mctTree_.eta->push_back(p.Eta());
	//	mctTree_.phi->push_back(p.Phi());
	//	mctTree_.pt->push_back(p.Pt());



	//		std::cout << "1_____" << prop << "____"<< p.Phi() <<   "_____"  <<conversion_legs[0].simTracks().charge() << std::endl;
					mctTree_.plabel->push_back(gen_counter);
					l_counter++;	
					mctTree_.label->push_back(l_counter);
					mctTree_.dr ->push_back(min_dr_mct);
					mctTree_.PId->push_back(2);
					mctTree_.primary_x ->push_back( mc_truth_pho->primaryVertex().x() );
					mctTree_.primary_y ->push_back( mc_truth_pho->primaryVertex().y() );
					mctTree_.primary_z ->push_back( mc_truth_pho->primaryVertex().z() );
					mctTree_.primary_phi ->push_back( mc_truth_pho->primaryVertex().phi());
					mctTree_.primary_eta ->push_back( mc_truth_pho->primaryVertex().eta());
					mctTree_.convRadius->push_back(mc_truth_pho->vertex().vect().rho());
					mctTree_.convPhi->push_back(mc_truth_pho->vertex().vect().phi());
					mctTree_.convZ->push_back(mc_truth_pho->vertex().vect().z());            
					mctTree_.conv_x->push_back(mc_truth_pho->vertex().vect().x());
					mctTree_.conv_y->push_back(mc_truth_pho->vertex().vect().y());
					mctTree_.conv_z->push_back(mc_truth_pho->vertex().vect().z());
					mctTree_.conv_phi->push_back(mc_truth_pho->vertex().vect().phi());
					mctTree_.conv_eta->push_back(mc_truth_pho->vertex().vect().eta());
					mctTree_.energy->push_back(mc_truth_pho->fourMomentum().e());
					mctTree_.nlegs->push_back(conversion_legs.size());
					part=RawParticle(conversion_legs[1].fourMomentum().x(),conversion_legs[1].fourMomentum().y(),conversion_legs[1].fourMomentum().z(),conversion_legs[1].fourMomentum().e());
			
					part.setVertex(conversion_legs[1].primaryVertex().x(),conversion_legs[1].primaryVertex().y(),conversion_legs[1].primaryVertex().z(),conversion_legs[1].primaryVertex().t());
					part.setCharge(conversion_legs[1].simTracks().charge());
					part.setMass(0.0005);
					p=BaseParticlePropagator(part,0,0,0);
					pos=GlobalPoint(conversion_legs[1].primaryVertex().x(),conversion_legs[1].primaryVertex().y(),conversion_legs[1].primaryVertex().z());
					p.setMagneticField(mag->inTesla(pos).z());
					p.SetPx(conversion_legs[1].fourMomentum().x());
					p.SetPy(conversion_legs[1].fourMomentum().y());
					p.SetPz(conversion_legs[1].fourMomentum().z());
				//	prop = p.propagateToEcalEntrance();
					std::cout << mag->inTesla(pos).z() << std::endl;
			mctTree_.eta->push_back(conversion_legs[1].fourMomentum().eta());
			mctTree_.phi->push_back(conversion_legs[1].fourMomentum().phi());
			mctTree_.pt->push_back(conversion_legs[1].fourMomentum().et());
	//				mctTree_.eta->push_back(p.Eta());
	//				mctTree_.phi->push_back(p.Phi());
	//				mctTree_.pt->push_back(p.Pt());

	l_counter++;	

	//	std::cout << "2_____" << prop << "____"<< p.Phi() <<   "_____"  <<conversion_legs[1].simTracks().charge() << std::endl;
					
			}
	}
			
	}
	std::cout << "n matches" << std::endl; 
	genTree_.GetTTreePtr()->Fill();
	mctTree_.GetTTreePtr()->Fill();
	dumbTree_.GetTTreePtr()->Fill();
			int i,j, index=-1,pindex=-1;
			std::cout << "ind def" << std::endl; 
			std::vector <int> indices;
			indices.resize(0);
			std::cout << "ind resized" << std::endl; 
			//std::cout << "ph vector size " << counter << endl;
			for(auto& cand:pfCandidates ){//loop on pf candidates

			float dr_min = 1e6;
			reco::PFCandidate* match = NULL;
			if(cand.pt()>0.5/* && (cand.particleId()==4 || cand.particleId()==2)*/){
				
			//	candTree_.ecalx ->push_back (cand.positionAtECALEntrance().X());
			//	candTree_.ecaly ->push_back (cand.positionAtECALEntrance().Y());
			//	candTree_.ecalz ->push_back (cand.positionAtECALEntrance().Z());
			//	candTree_.ecalphi ->push_back (cand.positionAtECALEntrance().phi());
			//	candTree_.ecaleta ->push_back (cand.positionAtECALEntrance().eta());
			//	std::cout << "before supercluster:" << std::endl;
			//	std::cout << *cand.superClusterRef() << std::endl;
			//	std::cout << *cand.superClusterRef()->seed() << std::endl;

				

				for(i=0;i<(int)mctTree_.plabel->size();i++){
				bool already=false;
				std::cout << "before if siz" << std::endl; 
				if(indices.size()>0){
				for (j=0;j<(int)indices.size();j++){
				std::cout << " in ind for" << std::endl; 
				if(mctTree_.label->at(i)==indices.at(j)){ 
				already=true;
				std::cout << " already =true" << std::endl; 
				break;
				}
				}
				}
				std::cout << "Pid size______________" << mctTree_.label->at(i) << "__________" << mctTree_.PId->at(i) << std::endl;
				float dr = deltaR(cand.eta(), cand.phi(), mctTree_.eta->at(i), mctTree_.phi->at(i));
				if(dr<dr_min && cand.particleId()>3 && already==false){
			
					dr_min = dr;
				std::cout << " doing matches---" << std::endl; 
					match = &cand;
					pindex = mctTree_.plabel->at(i);
					index = mctTree_.label->at(i);
					indices.push_back(index);
				}		
		
	
			}

}
			if(match){
				if(dr_min < 0.1/* && match->pt() >5) ||(dr_min < 0.3 && match->pt() <5)*/){
					candTree_.pt->push_back(cand.pt());
					candTree_.eta ->push_back(cand.eta());
					candTree_.phi ->push_back(cand.phi());
					candTree_.ecalEnergy ->push_back (cand.ecalEnergy());
					candTree_.hcalEnergy ->push_back(cand.hcalEnergy());
					candTree_.PId ->push_back(cand.particleId());
					candTree_.label->push_back(index);
					candTree_.plabel->push_back(pindex);
					std::cout << "______________" << index << "__________" << mctTree_.PId->at(index) << std::endl;
					
				}else{

					candTree_.pt->push_back(cand.pt());
					candTree_.eta ->push_back(cand.eta());
					candTree_.phi ->push_back(cand.phi());
					candTree_.ecalEnergy ->push_back (cand.ecalEnergy());
					candTree_.hcalEnergy ->push_back(cand.hcalEnergy());
					candTree_.PId ->push_back(cand.particleId());
					candTree_.label->push_back(99);
					candTree_.plabel->push_back(99);

				}
				}
}

	/*				float dr_ele1,dr_ele2;
					if(outTree_.mct_nlegs->at(index).back() >0) dr_ele1 = deltaR(cand.eta(), cand.phi(), mct_truth_pho->electrons()[0].eta, mct_truth_pho->electrons()[0].phi);
					else dr_ele1 = 1e6;
					if(outTree_.mct_nlegs->at(index).back() >1)  dr_ele2 = deltaR(cand.eta(), cand.phi(),mct_truth_pho->electrons()[1].eta, mct_truth_pho->electrons()[1].phi);
					else dr_ele2=1e6;

					if(dr_ele1 < dr_ele2 && dr_ele1 < 0.3){
					candTree_.label->push_back(counter);
					candTree_.plabel->push_back();
					mctTree_.PId->push_back(2);
					mctTree_.mct_convRadius->push_back(mc_truth_pho->vertex().vect().rho());
					mctTree_.mct_convPhi->push_back(mc_truth_pho->vertex().vect().phi());
					mctTree_.mct_convZ->push_back(mc_truth_pho->vertex().vect().z());            
					mctTree_.mct_energy->push_back(mc_truth_pho->fourMomentum().e());
					mctTree_.mct_eta->push_back(conversionlegs[0].fourMomentum().eta());
					mctTree_.mct_phi->push_back(conversionlegs[0]fourMomentum().phi());
					mctTree_.mct_pt->push_back(conversionlegs[0].fourMomentum().et());

					mctTree_.plabel->push_back(mctTree_.label->size());
					}else{ 



						if(dr_ele2 < 0.3){//checks for match with electron2
					mctTree_.label->push_back(l_counter);
					mctTree_.PId->push_back(2);
					mctTree_.mct_convRadius->push_back(mc_truth_pho->vertex().vect().rho());
					mctTree_.mct_convPhi->push_back(mc_truth_pho->vertex().vect().phi());
					mctTree_.mct_convZ->push_back(mc_truth_pho->vertex().vect().z());            
					mctTree_.mct_energy->push_back(mc_truth_pho->fourMomentum().e());
					mctTree_.mct_eta->push_back(conversionlegs[1].fourMomentum().eta()));
					mctTree_.mct_phi->push_back(conversionlegs[1]fourMomentum().phi());
					mctTree_.mct_pt->push_back(conversionlegs[1].fourMomentum().et());

					mctTree_.plabel->push_back(mctTree_.label->size());

						}else 
					}

				




			
			else{//fills the branches to let the size be always the same 
				outTree_.dr->at(index).push_back(-1);
				outTree_.label->at(index).push_back(-1);
				outTree_.pt->at(index).push_back(-1);
				outTree_.phi->at(index).push_back(-1);
				outTree_.eta->at(index).push_back(-1);
				outTree_.particleId->at(index).push_back(-1);
				outTree_.ecalEnergy->at(index).push_back(-1);
				outTree_.hcalEnergy->at(index).push_back(-1);
				outTree_.ele2_match->at(index).push_back(false);
				outTree_.ele1_match->at(index).push_back(false);


			}
			//	std::cout <<"_____" << outTree_.dr->at(index).size() << "_____" << outTree_.ele1_match->at(index).size() << "_____" << outTree_.ele2_match->at(index).size() << std::endl;
		




	
			if(mc_truth_pho->vertex().vect().rho()>)
			outTree_.mct_convRadius->at(counter).push_back(mc_truth_pho->vertex().vect().rho());
			outTree_.mct_convPhi->at(counter).push_back(mc_truth_pho->vertex().vect().phi());
			outTree_.mct_convZ->at(counter).push_back(mc_truth_pho->vertex().vect().z());            
			outTree_.mct_energy->at(counter).push_back(mc_truth_pho->fourMomentum().e());
			outTree_.mct_eta->at(counter).push_back(mc_truth_pho->fourMomentum().eta());
			outTree_.mct_phi->at(counter).push_back(mc_truth_pho->fourMomentum().phi());
			outTree_.mct_pt->at(counter).push_back(mc_truth_pho->fourMomentum().et());
			auto conversion_legs = mc_truth_pho->electrons();
			outTree_.mct_nlegs->at(counter).push_back(conversion_legs.size());
			if(conversion_legs.size()>0)
			{
				outTree_.mct_ele1_pt->at(counter).push_back(conversion_legs[0].fourMomentum().et());
				outTree_.mct_ele1_eta->at(counter).push_back(conversion_legs[0].fourMomentum().eta());
				outTree_.mct_ele1_phi->at(counter).push_back(conversion_legs[0].fourMomentum().phi());

				if(conversion_legs.size()>1)
				{

					outTree_.mct_ele2_pt->at(counter).push_back(conversion_legs[1].fourMomentum().et());
					outTree_.mct_ele2_eta->at(counter).push_back(conversion_legs[1].fourMomentum().eta());
					outTree_.mct_ele2_phi->at(counter).push_back(conversion_legs[1].fourMomentum().phi());
					outTree_.mct_eles_dr->at(counter).push_back(deltaR(outTree_.mct_ele1_eta->at(counter).back(), outTree_.mct_ele1_phi->at(counter).back(),
								outTree_.mct_ele2_eta->at(counter).back(), outTree_.mct_ele2_phi->at(counter).back()));





				}


			}

		
ccounter++;//moves to the other photon
	


	//std::cout <<"_____" << outTree_.mct_pt->size() << "_____" << outTree_.mct_ele1_pt->at(0).size() << "_____" << outTree_.mct_ele2_pt->at(1).size() << std::endl;
	//fill the first two entries of the candidates -just to help with visualization on root scan
	//fills tree*/

	std::set<MTDHitMatchingInfo> clus_match;                
	//---BTL clusters loop
	for(auto& clusIt : clustersBTL)
	{    
		DetId id = clusIt.detId();
		const auto& det = mtdGeometry_ -> idToDet(id);
		const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(det->topology());
		const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());    
		for(auto& cluster : clusIt)
		{
			float energy = cluster.energy();
			float time   = cluster.time();
			int size=cluster.size();
			int sizeX=cluster.sizeX();
			int sizeY=cluster.sizeY();
			//float seed_energy=cluster.seed().energy();
		//	float seed_time=cluster.seed().time();
			int seed_x=cluster.seed().x();
			int seed_y=cluster.seed().y();

			MeasurementPoint mp(cluster.x(),cluster.y());
			LocalPoint lp = topo.localPosition(mp);
			GlobalPoint gp = det->toGlobal(lp);

		//	auto dr_clus_neu = deltaR(cand.eta(), cand.phi(), gp.eta(), gp.phi());
	//		if(dr_clus_neu > 0.1)
	//			continue;

		//	auto mtd_neu_match_chi2s = getNeutralsMTDMatchingChi2s(cand, gp, genPV, time);
			//match.insert(MTDHitMatchingInfo(candTree_.clus_n, mtd_neu_match_chi2s.first, mtd_neu_match_chi2s.second));

			/*MTDDetId mtdId(id);
			int RR = 0;

			if ( mtdId.mtdSubDetector() == MTDDetId::BTL )
			{
				BTLDetId btlId(id);
				RR = btlId.mtdRR();
			} */                      

			candTree_.clus_n += 1;    
			candTree_.clus_size->push_back(size);
			candTree_.clus_size_x->push_back(sizeX);
			candTree_.clus_size_y->push_back(sizeY);
			candTree_.clus_energy->push_back(energy);
			candTree_.clus_time->push_back(time);
			candTree_.clus_z->push_back(gp.z());
			candTree_.clus_x->push_back(gp.x());
			candTree_.clus_y->push_back(gp.y());
			candTree_.clus_eta->push_back(gp.eta());
			candTree_.clus_phi->push_back(gp.phi());
			
			std::cout << "ECAL rad " << sqrt(lp.x()*lp.x()+lp.y()*lp.y()) << std::endl;
		}
	}


		
	

	
	candTree_.GetTTreePtr()->Fill();
	
	std::cout << "HEREEEEE" << std::endl;




	bool debug = false;	
	for(unsigned int itrack=0; itrack<tracks.size(); ++itrack){
	reco::Track cand = tracks[itrack];	
	reco::TrackBaseRef track_ref(tracksHandle_, itrack);
	reco::TrackBaseRef ext_track_ref;
        if(generalToExtendedTrkMap[track_ref] != -1)
            ext_track_ref = reco::TrackBaseRef(extTracksHandle_, generalToExtendedTrkMap[track_ref]);
		if(cand.pt()>0.5 && (cand.charge()==1 || cand.charge()==-1)){

			const reco::Track * trackAddress = &cand;
			const MagneticField* mag  = magfield.product();	

			//;
			//
			//		outTree_.Track_zref = cand.vz();
			outTreet_.Track_ECAL_x ->push_back( spr::propagateTrackToECAL(trackAddress,mag,debug).point.X());
			outTreet_.Track_ECAL_y ->push_back( spr::propagateTrackToECAL(trackAddress,mag,debug).point.Y());
			outTreet_.Track_ECAL_z ->push_back( spr::propagateTrackToECAL(trackAddress,mag,debug).point.Z());
			outTreet_.Track_ECAL_phi ->push_back( spr::propagateTrackToECAL(trackAddress,mag,debug).point.Phi());
			outTreet_.Track_ECAL_eta ->push_back( spr::propagateTrackToECAL(trackAddress,mag,debug).point.Eta());
			outTreet_.Track_pt  ->push_back(  cand.pt());
			outTreet_.Track_In_x  ->push_back(  cand.innerPosition().X());
			outTreet_.Track_In_y  ->push_back(  cand.innerPosition().Y());
			outTreet_.Track_In_z  ->push_back(  cand.innerPosition().Z());
			outTreet_.Track_In_phi  ->push_back(  cand.innerPosition().Phi());
			outTreet_.Track_In_eta  ->push_back(  cand.innerPosition().Eta());
			outTreet_.Track_charge  ->push_back( cand.charge());
			outTreet_.Track_nhits  ->push_back( cand.found());
			outTreet_.Track_pt_Iphi  ->push_back( cand.innerMomentum().phi());
			outTreet_.Track_pt_Ephi  ->push_back( cand.outerMomentum().phi());
			outTreet_.Track_pt_Ieta  ->push_back( cand.innerMomentum().eta());
			outTreet_.Track_pt_Eeta  ->push_back( cand.outerMomentum().eta());
			outTreet_.Track_mtdt -> push_back(!ext_track_ref.isNull() ? extTracksMTDtime[ext_track_ref] : -1);
			
			
			float DRMin = 1e9;
			for(unsigned int iPart=0; iPart<genParticles.size(); ++iPart)
			{
				reco::GenParticle genPart = genParticles[iPart];

				if( genPart.status() != 1 ) continue;
				//if( genPart.pdgId() == 22 ) continue;
	
			/*	float DR   = deltaR(cand.eta(), cand.phi(), genPart.eta(), genPart.phi());

				if( DR < DRMin )
				{
					DRMin = DR;

				*/
				outTreet_.Track_genVtx_t -> push_back(genPV.t());
				break;
				//}else outTreet_.Track_genVtx_t -> push_back(-1);

			}
		
			std::cout << " time 0  "<<  	outTreet_.Track_genVtx_t ->back() << "time MTD"<< outTreet_.Track_mtdt -> back()<< std::endl;
	}


		else{
			//	outTree_.Track_ECALx = 9999;
			//	outTree_.Track_ECALy = 9999;
			//	outTree_.Track_zref = 9999;
			//	outTree_.Track_etaref = cand.momentum().Eta();
			outTreet_.Track_ECAL_x ->push_back( 9999);
			outTreet_.Track_ECAL_y ->push_back( 9999);
			outTreet_.Track_ECAL_z ->push_back( 9999);
			outTreet_.Track_ECAL_phi ->push_back( 9999);
			outTreet_.Track_ECAL_eta ->push_back( 9999);
			outTreet_.Track_pt  ->push_back(  9999);
			outTreet_.Track_charge  ->push_back(  9999);
			outTreet_.Track_In_x  ->push_back(9999);
			outTreet_.Track_In_y  ->push_back(9999);
			outTreet_.Track_In_z  ->push_back(9999);
			outTreet_.Track_In_phi  ->push_back(9999);
			outTreet_.Track_In_eta  ->push_back(9999);
			outTreet_.Track_nhits  ->push_back(9999);
			outTreet_.Track_pt_Iphi  ->push_back( 9999);
			outTreet_.Track_pt_Ephi  ->push_back( 9999);
			outTreet_.Track_pt_Ieta  ->push_back( 9999);
			outTreet_.Track_pt_Eeta  ->push_back( 9999);

		} 
	}

	std::cout << "HERE_____2" << std::endl;
	outTreet_.GetTTreePtr()->Fill();              
}

DEFINE_FWK_MODULE(MTDNeutralsAnalyzer);

#endif
