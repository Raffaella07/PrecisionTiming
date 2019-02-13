import FWCore.ParameterSet.Config as cms

MTD4DVertexingAnalyzer = cms.EDAnalyzer('MTD4DVertexingAnalyzer',
                                        genParticlesTag = cms.untracked.InputTag("genParticles"),
                                        trackingParticlesTag = cms.untracked.InputTag("mix", "MergedTrackTruth"),
                                        trackAndTrackingParticlesAssociatorMapTag = cms.untracked.InputTag("trackingParticleRecoTrackAsssociation"),
                                        generalTracksTag = cms.untracked.InputTag("generalTracks"),
                                        extendedTracksTag = cms.untracked.InputTag("trackExtenderWithMTD", "", "MTDRERECO"),
                                        btlMatchChi2Tag = cms.untracked.InputTag("trackExtenderWithMTD", "btlMatchChi2", "MTDRERECO"),
                                        btlMatchTimeChi2Tag = cms.untracked.InputTag("trackExtenderWithMTD", "btlMatchTimeChi2", "MTDRERECO"),
                                        etlMatchChi2Tag = cms.untracked.InputTag("trackExtenderWithMTD", "etlMatchChi2",  "MTDRERECO"),
                                        etlMatchTimeChi2Tag = cms.untracked.InputTag("trackExtenderWithMTD", "etlMatchTimeChi2", "MTDRERECO"),
                                        extTracksPathLengthTag = cms.untracked.InputTag("trackExtenderWithMTD", "pathLength", "MTDRERECO"),
                                        extTracksMTDtimeTag = cms.untracked.InputTag("trackExtenderWithMTD", "tmtd", "MTDRERECO"),
                                        t0TOFPIDTag = cms.untracked.InputTag("tofPID", "t0", "MTDRERECO"),
                                        sigmat0TOFPIDTag = cms.untracked.InputTag("tofPID", "sigmat0", "MTDRERECO"),
                                        probPiTOFPIDTag = cms.untracked.InputTag("tofPID", "probPi", "MTDRERECO"),
                                        probPTOFPIDTag = cms.untracked.InputTag("tofPID", "probP", "MTDRERECO"),
                                        probKTOFPIDTag = cms.untracked.InputTag("tofPID", "probK", "MTDRERECO"),
                                        genXYZTag = cms.untracked.InputTag("genParticles", "xyz0", "HLT"),
                                        genT0Tag = cms.untracked.InputTag("genParticles", "t0", "HLT"),
                                        simVtxTag = cms.untracked.InputTag("mix", "InitialVertices", "HLT"),
                                        vtx3DTag = cms.untracked.InputTag("offlinePrimaryVertices", "", "RECO"),
                                        vtx4DTag = cms.untracked.InputTag("offlinePrimaryVertices4D", "", "MTDRERECO"),                                        
                                        vtx4DNoPIDTag = cms.untracked.InputTag("offlinePrimaryVertices4DnoPID", ""),
                                        trackPUID_3DBDT_weights_file = cms.FileInPath("PrecisionTiming/FTLAnalysis/data/test_gist_clf3D.xml"),
                                        trackPUID_4DBDT_weights_file = cms.FileInPath("PrecisionTiming/FTLAnalysis/data/test_gist_clf4D.xml"),
                                        vtxsTreeName = cms.untracked.string("vtxs_tree"),
                                        trksTreeName = cms.untracked.string("trks_tree")
)
