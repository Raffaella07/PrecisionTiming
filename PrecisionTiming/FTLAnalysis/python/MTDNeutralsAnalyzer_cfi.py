import FWCore.ParameterSet.Config as cms  

MTDNeutralsAnalyzer = cms.EDAnalyzer('MTDNeutralsAnalyzer',
		genParticlesTag = cms.untracked.InputTag("genParticles"),
		simTkTag = cms.untracked.InputTag("g4SimHits"),
		simVtxTag = cms.untracked.InputTag("g4SimHits"),
		simHitsBTLTag = cms.untracked.InputTag("g4SimHits:FastTimerHitsBarrel"),
		simHitsETLTag = cms.untracked.InputTag("g4SimHits:FastTimerHitsEndcap"),
		clustersBTLTag = cms.untracked.InputTag("mtdClusters:FTLBarrel"),                                     
		clustersETLTag = cms.untracked.InputTag("mtdClusters:FTLEndcap"),
		pfCandidatesTag = cms.untracked.InputTag("particleFlow"),
		#pfCandidatesTag = cms.untracked.InputTag("particleFlowClusterECAL"),
		genXYZTag = cms.untracked.InputTag("genParticles", "xyz0"),
		genT0Tag = cms.untracked.InputTag("genParticles", "t0"),
		vtx3DTag = cms.untracked.InputTag("offlinePrimaryVertices", ""),
		vtx4DTag = cms.untracked.InputTag("offlinePrimaryVertices4D", ""),
		crysLayout = cms.untracked.int32(3),
#electrontracks                             
		generalTracksTag = cms.untracked.InputTag("generalTracks"),
		extendedTracksTag = cms.untracked.InputTag("trackExtenderWithMTD"),                             
               
		#genvtxtag = cms.untracked.inputtag("g4simhits"),
##########		
		outTreeName1 = cms.untracked.string("neu_tree"),
		outTreeName2 = cms.untracked.string("tracks_tree"),
		)
