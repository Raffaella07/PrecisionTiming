#ifndef _MTD_NEUTRALS__TREE_
#define _MTD_NEUTRALS__TREE_

#include "ExternalTools/DynamicTTree/interface/DynamicTTreeBase.h"

using namespace std;
//---Define the TTree branches
#define DYNAMIC_TREE_NAME MTDNeutralsTree
#define DATA_TABLE                              \
    DATA(int, run)                              \
    DATA(int, lumi)                             \
    DATA(int, event)                            \
    DATA(float, Track_ECALx)                     \
    DATA(float, Track_eta)                   \
    DATA(float, Track_pt)                       \
    DATA(float, Track_phi)                       \
    DATA(int, Track_charge)                      \
   /* DATA(float, Track_ECALy) */                    \

#define DATA_CLASS_TABLE                     \
    DATA(vector <float>, gen_DR)                         \
    DATA(vector <float>, genpv_z)                        \
    DATA(vector <float>, pt)                             \
    DATA(vector <float>, eta)                            \
    DATA(vector <float>, phi)                            \
    DATA(vector <float>, ecalEnergy)                     \
    DATA(vector <float>, hcalEnergy)                 \
    DATA(vector <float>, particleId)                     \
    DATA(vector <int>, clus_n)                           \
    DATA(vector < float>, gen_pt)                         \
    DATA(vector < float>, gen_eta)                        \
    DATA(vector < float>, gen_phi)                        \
    DATA(vector < int>, gen_pdgId)                        \
    DATA(vector < float>, mct_pt)                         \
    DATA(vector < float>, mct_eta)                        \
    DATA(vector < float>, mct_phi)                        \
    DATA(vector < float>, mct_energy)                     \
    DATA(vector < float>, mct_convRadius)                 \
    DATA(vector < float>, mct_convZ)                      \
    DATA(vector < float>, mct_convPhi)                    \
    DATA(vector < float>, mct_ele1_pt)                    \
    DATA(vector < float>, mct_ele1_eta)                   \
    DATA(vector < float>, mct_nlegs)                   \
    DATA(vector < float>, mct_ele1_phi)                   \
    DATA(vector < float>, ele1_match)                   \
    DATA(vector < float>, mct_ele2_pt)                    \
    DATA(vector < float>, mct_ele2_eta)                   \
    DATA(vector < float>, mct_ele2_phi)                   \
    DATA(vector < float>, ele2_match)                   \
    DATA(vector < float>, mct_eles_dr)                    \
    DATA(vector < float>, mct_ele1_dr)                    \
    DATA(vector < float>, mct_ele2_dr)                    \
 /*   DATA(vector <float>, minDR)                          \
    DATA(vector <float>, minDEta)                        \
    DATA(vector <float>,minDPhi)                        \
    DATA(vector <float>,tof)                            \
    DATA(vector <float>,mtdTime)                        \
    DATA(vector <float>,mtdEnergy)                      \
    DATA(vector <float>,simh_DR)                        \
    DATA(vector <float>,simh_recoh_DR)                  \
    DATA(vector <float>,simh_recoh_DPhi)                \
    DATA(vector <float>,simh_energy)                    \
    DATA(vector <float>,simh_time)                      \
    DATA(vector <float>,simh_tof)                       \
    DATA(vector <float>,simh_clus_DR)                   \
    DATA(vector <float>,simh_recoh_clus_DR)             \
    DATA(vector <float>,simh_recoh_clus_DPhi)           \
    DATA(vector <float>,simh_clus_energy)               \
    DATA(vector <float>,simh_clus_time)                 \
    DATA(vector <float>,simh_clus_tof)                  \
    DATA(vector <int>, mct_nlegs)                        \
    DATA(vector<float>, clus_size)              \
    DATA(vector<float>, clus_size_x)            \
    DATA(vector<float>, clus_size_y)            \
    DATA(vector<float>, clus_energy)            \
    DATA(vector<float>, clus_time)              \
    DATA(vector<float>, clus_rr)                \
    DATA(vector<float>, clus_eta)               \
    DATA(vector<float>, clus_phi)               \
    DATA(vector<float>, clus_seed_energy)       \
        DATA(vector<float>, clus_seed_time)     \
    DATA(vector<float>, clus_seed_x)            \
    DATA(vector<float>, clus_seed_y)            \
    DATA(vector<float>, clus_x)                 \
    DATA(vector<float>, clus_y)                 \
    DATA(vector<float>, clus_z)                 \
    DATA(vector<float>, clus_global_R)          \
    DATA(vector<float>, clus_global_dist)       \
    DATA(vector<float>, clus_neu_DPhi)          \
    DATA(vector<float>, clus_neu_DEta)          \
    DATA(vector<float>, clus_neu_DR)            \
    DATA(vector<float>, clus_cele1_DR)          \
    DATA(vector<float>, clus_cele2_DR)       */   \
    DATA(vector<float > , genpv_t,)
#include "ExternalTools/DynamicTTree/interface/DynamicTTreeInterface.h"

#endif
