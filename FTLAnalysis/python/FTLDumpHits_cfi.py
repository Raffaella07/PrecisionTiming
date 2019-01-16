import FWCore.ParameterSet.Config as cms

FTLDumpHits = cms.EDAnalyzer('FTLDumpHits',
                             genParticlesTag = cms.untracked.InputTag("genParticles"),
                             trackingParticlesTag = cms.untracked.InputTag("mix", "MergedTrackTruth"),
                             trackAndTrackingParticlesAssociatorMapTag = cms.untracked.InputTag("trackingParticleRecoTrackAsssociation"),
                             simHitsBTLTag = cms.untracked.InputTag("g4SimHits:FastTimerHitsBarrel"),
                             recHitsBTLTag = cms.untracked.InputTag("mtdRecHits:FTLBarrel"),
                             clustersBTLTag = cms.untracked.InputTag("mtdClusters:FTLBarrel"),
                             simHitsETLTag = cms.untracked.InputTag("g4SimHits:FastTimerHitsEndcap"),
                             recHitsETLTag = cms.untracked.InputTag("mtdRecHits:FTLEndcap"),
                             clustersETLTag = cms.untracked.InputTag("mtdClusters:FTLEndcap"),
                             generalTracksTag = cms.untracked.InputTag("generalTracks"),
                             extendedTracksTag = cms.untracked.InputTag("trackExtenderWithMTD"),                             
                             genVtxTag = cms.untracked.InputTag("g4SimHits"),
                             crysLayout = cms.untracked.int32(0),
                             track_hit_DRMax = cms.double(0.05),
                             track_hit_distMax = cms.double(99999.),
                             treeName = cms.untracked.string("DumpHits"),
                             verbosity = cms.bool(False),
                             dumpSimHits = cms.bool(False),
                             dumpRecHits = cms.bool(False)
                             )
