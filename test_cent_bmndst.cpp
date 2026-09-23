#include "ROOT/RVec.hxx"
#include "RtypesCore.h"
#define BMNHYPNUCLPAIR_H

#include <ROOT/RDataFrame.hxx>
#include <Math/LorentzVector.h>
#include <Math/PtEtaPhiE4D.h>
#include <vector>
#include <string>
#include <typeinfo>
#include <cxxabi.h>

using namespace ROOT;
using namespace ROOT::Math;
using namespace ROOT::RDF;

constexpr auto BT_mask 	 = []( UInt_t map ){ return map & (1<<2); };
constexpr auto MBT_mask  = []( UInt_t map ){ return map & (1<<5); };
constexpr auto CCT1_mask = []( UInt_t map ){ return map & (1<<4); };
constexpr auto CCT2_mask = []( UInt_t map ){ return map & (1<<7); };

void RemoveBadRuns(std::vector<int> &_Runs, std::vector<int> _BadRUns){
	for(auto _run : _BadRUns){
		_Runs.erase(remove(_Runs.begin(), _Runs.end(), _run), _Runs.end());
	}
}

TChain* makeChain(string& filename, const char* treename) {
    cout << "Adding files to chain:" << endl;
    TChain *chain = new TChain(treename);
    if (filename.rfind(".root") < filename.size())
      chain->Add(filename.data());
    else {
      TFileCollection fc("fc", "", filename.c_str());
      chain->AddFileInfoList((TCollection*)fc.GetList());
    }
    chain->ls();
    return chain;
}

//pile-up by Oleg Golosov
int getValue(RVec<int> values, int index)
try {
  if(values.size()==0)return -1;
  return values.at(index);
} catch (const std::exception& e) {
  std::cout << __func__ << ": " << e.what() << std::endl;
  throw;
}

int CentralHitIndexBC1S(BmnBC1hitInfo hit)
try {
  return hit.GetCentralHitIndexBC1S();
} catch (const std::exception& e) {
  std::cout << __func__ << ": " << e.what() << std::endl;
  throw;
}

int BC1hitClasses(BmnBC1hitInfo hit, int centralHitIndex)
try {
  if(centralHitIndex<0)return -1;
  return hit.GetBC1hitClasses().at(centralHitIndex)==BmnEventClass::k1;
} catch (const std::exception& e) {
  std::cout << __func__ << ": " << e.what() << std::endl;
  throw;
}

auto ClosestBC1hitsDt(BmnEventClass::id classId) {
  return [classId](BmnBC1hitInfo& hit) -> RVecF {
    try {
      return hit.GetClosestBC1hitsDt(classId);
    } catch (const std::exception& e) {
      std::cout << __func__ << ": " << e.what() << std::endl;
      throw;
    }
  };
};

// global track
enum class TrackParamConvert { Pq, P, Pt, Phi, Eta, Pz, Dedx, Q};

auto GlobalTrackParamFloat(TrackParamConvert param) {
  return [param](const RVec<BmnGlobalTrack> tracks) -> RVecF {
    RVecF result;
    result.reserve(tracks.size());
    
    for (auto track:tracks) {
      auto* par = track.GetParamFirst();
      if (!par) {
          result.push_back(-999.0f);
          continue;
      }
      switch (param) {
        case TrackParamConvert::Pq:
          result.push_back(1.0 / par->GetQp());
          break;
        case TrackParamConvert::P:
          result.push_back(std::abs(1.0 / par->GetQp()));
          break;
        case TrackParamConvert::Pt: {
          TVector3 mom;
          par->Momentum(mom);
          result.push_back(mom.Pt());
          break;
        }
        case TrackParamConvert::Phi: {
          TVector3 mom;
          par->Momentum(mom);
          result.push_back(mom.Phi());
          break;
        }
        case TrackParamConvert::Eta: {
          TVector3 mom;
          par->Momentum(mom);
          result.push_back(mom.Eta());
          break;
        }
        case TrackParamConvert::Pz: {
          TVector3 mom;
          par->Momentum(mom);
          result.push_back(mom.Pz());
          break;
        }
        case TrackParamConvert::Dedx: {
          result.push_back(track.GetdQdNLower());
          break;
        }
        case TrackParamConvert::Q: {
          result.push_back(track.GetZ());
          break;
        }
      }
    }
    return result;
  };
}

//DCA gl tr
enum class DcaTypeConvert { X, Y, Z, R };

const auto recDcaXYZR = [](DcaTypeConvert type) {
  return [type](const RVec<BmnGlobalTrack> tracks, const CbmVertex vtx) {
    RVecF dca;
    dca.reserve(tracks.size());
    
    float vx = vtx.GetX();
    float vy = vtx.GetY();
    float vz = vtx.GetZ();
    
    for (auto track : tracks) {
      auto* par = track.GetParamFirst();
      if (!par) {
        dca.push_back(-999.0f);
        continue;
      }
      switch (type) {
        case DcaTypeConvert::X: dca.push_back(par->GetX() - vx); break;
        case DcaTypeConvert::Y: dca.push_back(par->GetY() - vy); break;
        case DcaTypeConvert::Z: dca.push_back(par->GetZ() - vz); break;
        case DcaTypeConvert::R: {
          float dx = par->GetX() - vx;
          float dy = par->GetY() - vy;
          dca.push_back(std::hypot(dx, dy));
          break;
        }
      }
    }
    return dca;
  };
};







void test_cent_bmndst(
    std::string in_file_name="/eos/nica/bmn/exp/dst/run8/25.09.0/mpd_run_Top_8130_ev0_p16.root",
    std::string outPath="./",
    std::string in_file_vtx="./run8_25.09_corr_VtxXYZ.root", 
    std::string in_file_RefMultCorr="./run8_25.09_corrections_DifRefRunId.root")
{

    // In File
    const char* TreeName = "bmndata";
    TString name = gSystem->BaseName(in_file_name.c_str()); // mpd_run_Top_8130_ev0_p16.root
    name.ReplaceAll(".root", "");
    std::string only_name_iFileDst = name.Data();
    TChain *chain=makeChain(in_file_name, TreeName);
    
    int nEvents = chain->GetEntries();
    std::cout<<"nEvents:\t"<<nEvents<<std::endl;

    DstRunHeader* run_header = dynamic_cast<DstRunHeader*>( chain->GetCurrentFile()->Get("DstRunHeader") );
    if (run_header) {
      cout << "input dst file:\t"<<only_name_iFileDst<<".root"<<endl;
      cout << "\n|||||||||||||||| RUN SUMMARY |||||||||||||||" << endl;
      cout << "||\t\t\t\t\t  ||" << endl;
      cout << "||   Period:        " << run_header->GetPeriodId() << "\t\t\t  ||" << endl;
      cout << "||   Number:        " << run_header->GetRunId() << "\t\t  ||" << endl;
      cout << "||   Start Time:    " << run_header->GetStartTime().AsString("s") << "\t  ||" << endl;
      cout << "||   End Time:      " << run_header->GetFinishTime().AsString("s") << "\t  ||" << endl;
      cout << "||   Beam:          A = " << run_header->GetBeamA() << ", Z = " << run_header->GetBeamA() << "\t  ||" << endl;
      cout << "||   Beam energy:   " << run_header->GetBeamEnergy() << " GeV\t\t  ||" << endl;
      cout << "||   Target:        A = " << run_header->GetTargetA() << ", Z = " << run_header->GetTargetZ() << "\t  ||" << endl;
      //cout << "||   Field voltage: " << setprecision(4) << run_header->GetMagneticField() << " mV\t\t  ||" << endl;
      cout << "||\t\t\t\t\t  ||" << endl;
      cout << "||||||||||||||||||||||||||||||||||||||||||||\n" << endl;
    }
    
    auto run_id = run_header->GetRunId();
  
	auto file_fit_VtxXYZ = TFile::Open( in_file_vtx.c_str(), "READ" );
	file_fit_VtxXYZ->cd();
	
	TGraphErrors* g1_FitVtxX = nullptr;
	TGraphErrors* g1_FitVtxY = nullptr;
	TGraphErrors* g1_FitVtxZ = nullptr;
	TGraphErrors* f1_FitRunIdFactor_1 = nullptr;

	if(file_fit_VtxXYZ!=nullptr){
		file_fit_VtxXYZ->cd();
		g1_FitVtxX = file_fit_VtxXYZ->Get<TGraphErrors>("grGausNew_cct_nVtxMpd2_h2_RunId_vtx_x_mpd");
		g1_FitVtxY = file_fit_VtxXYZ->Get<TGraphErrors>("grNew_cct_nVtxMpd2_h2_RunId_vtx_y_mpd");
		g1_FitVtxZ = file_fit_VtxXYZ->Get<TGraphErrors>("grNew_cct_nVtxMpd2_h2_RunId_vtx_z_mpd");
	}else{
		std::cout<<"Warning! file_fit_VtxXYZ not found!"<<std::endl;
	}

	auto file_RefMultCorr = TFile::Open( in_file_RefMultCorr.c_str(), "READ" );
	file_RefMultCorr->cd();
	if(file_RefMultCorr!=nullptr){
		file_RefMultCorr->cd();
		f1_FitRunIdFactor_1 = file_RefMultCorr->Get<TGraphErrors>("RunId_corr_factor_h2_RunId_RefMult_gt_mpd_8185_8195");
		//f1_FitRunIdFactor_1 = file_RefMultCorr->Get<TGraphErrors>("RunId_corr_factor_h2_RunId_RefMult_gt_mpd_7643_7648");
		//f1_FitRunIdFactor_1 = file_RefMultCorr->Get<TGraphErrors>("RunId_corr_factor_h2_RunId_RefMult_gt_mpd_7386_7400");
	}

	///////// Warning part /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if(g1_FitVtxX==nullptr) std::cout<<"Warning! No corrections found for VtxX. These corrections will not be applied."<<std::endl;
	if(g1_FitVtxY==nullptr) std::cout<<"Warning! No corrections found for VtxY. These corrections will not be applied."<<std::endl;
	if(g1_FitVtxZ==nullptr) std::cout<<"Warning! No corrections found for VtxZ. These corrections will not be applied."<<std::endl;
	if(f1_FitRunIdFactor_1==nullptr) std::cout<<"Warning! No corrections found for RefMult. These corrections will not be applied."<<std::endl;

	//RefMult cuts
	float _ptMin = 0.05;
	float _ptMax = 2.0;
	float _etaMin = 0.7;
	float _etaMax = 2.5;
	float _dcaR = 1.;

	auto RefMult_gt = [](float ptMin, float ptMax, float etaMin, float etaMax, float dcaRCut){ 
		return [ptMin, ptMax, etaMin, etaMax, dcaRCut](RVecF tr_pt, RVecF tr_eta, RVecF tr_dca_r){
		  unsigned long Mult = 0;
		  for( int i=0; i<tr_pt.size(); ++i ){
			if(tr_pt.at(i) < ptMin) continue;
			if(tr_pt.at(i) > ptMax) continue;
			if(tr_eta.at(i) < etaMin) continue;
			if(tr_eta.at(i) > etaMax) continue;
			if(tr_dca_r.at(i) > dcaRCut) continue;
			Mult += 1;
		  }
		  return Mult;
		};
	  };

	  auto centrality_function = [](int NumVecMult) {
		return [NumVecMult](Double_t _refMult) -> int {
			// Границы для binMult = _refMult + 1 
			
			bool OutCentBin=true; 
			// если true вывод в формате 0 [0-10%], 1 [10-20%], 3 [20-30%], 4 [30-40%], 5 [40-50%], 6 [50-60%], -1 - событие выше анкер поинта, или не прошло по критериям отбора
			// если false вывод в формате центра бина 5 [0-10%], 15 [10-20%], 25 [20-30%],...

			std::vector<float> multiplicity_edges;			
			std::vector<float> centrality_percentage{ 0, 10, 20, 30, 40, 50, 60, 70, 100};
			
			if (NumVecMult == 1) multiplicity_edges = { 115, 57, 41, 29, 20, 13, 8, 4, 0}; // 8185-8195
			//if (NumVecMult == 2) multiplicity_edges = { 109, 51, 37, 26, 18, 12, 7, 4, 0}; // 7643-7648
			//if (NumVecMult == 3) multiplicity_edges = {  95, 46, 33, 24, 17, 11, 7, 4, 0}; // 7386-7400

			if (multiplicity_edges.empty()) return -1;
			
			if (_refMult < 0 || _refMult >= multiplicity_edges[0]) return -1;
			
			int idx = 0;
			float bin_edge = multiplicity_edges[idx];
			while (_refMult < bin_edge && idx < (int)multiplicity_edges.size() - 1) {
				idx++;
				bin_edge = multiplicity_edges[idx];
			}

			if(idx-1>=6) return -1; // cent > 60% нету
			
			if(OutCentBin){
				return (idx - 1);
			}else{
				return (centrality_percentage[idx-1] + centrality_percentage[idx])/2;
			}

			return -1;
		};
	};

    std::vector<int> physical_runs{6667, 6668, 6669, 6670, 6671, 6672, 6673, 6674, 6675, 6676, 6677, 6678, 6679, 6680, 6681, 6683, 6684, 6685, 6686, 6687, 6689, 6690, 6691, 6692, 6694, 6695, 6696, 6698, 6699, 6732, 6733, 6734, 6737, 6738, 6739, 6740, 6745, 6752, 6753, 6760, 6761, 6765, 6766, 6767, 6768, 6769, 6771, 6772, 6773, 6774, 6779, 6780, 6782, 6783, 6785, 6786, 6788, 6794, 6795, 6797, 6799, 6800, 6803, 6815, 6816, 6817, 6818, 6819, 6820, 6821, 6822, 6879, 6882, 6883, 6884, 6886, 6887, 6889, 6891, 6900, 6901, 6902, 6903, 6904, 6905, 6906, 6907, 6908, 6909, 6910, 6911, 6915, 6916, 6918, 6919, 6920, 6921, 6923, 6924, 6926, 6927, 6928, 6929, 6930, 6931, 6932, 6933, 6934, 6935, 6936, 6937, 6939, 6940, 6968, 6970, 6972, 6973, 6975, 6976, 6977, 6978, 6979, 6980, 6981, 6982, 6983, 6984, 6990, 6991, 6992, 6993, 6994, 6995, 6997, 6998, 6999, 7000, 7002, 7003, 7004, 7005, 7006, 7008, 7009, 7010, 7011, 7012, 7030, 7031, 7032, 7033, 7034, 7035, 7037, 7038, 7040, 7041, 7042, 7043, 7044, 7046, 7047, 7048, 7049, 7050, 7051, 7052, 7053, 7054, 7055, 7056, 7075, 7076, 7077, 7078, 7081, 7082, 7083, 7084, 7086, 7087, 7091, 7092, 7093, 7094, 7096, 7097, 7098, 7100, 7101, 7102, 7103, 7104, 7125, 7126, 7127, 7128, 7129, 7130, 7131, 7132, 7133, 7135, 7136, 7137, 7138, 7146, 7149, 7150, 7151, 7154, 7155, 7156, 7157, 7159, 7160, 7161, 7162, 7163, 7164, 7165, 7166, 7167, 7168, 7173, 7174, 7175, 7176, 7177, 7178, 7179, 7180, 7181, 7182, 7184, 7186, 7187, 7188, 7191, 7192, 7193, 7194, 7195, 7200, 7202, 7203, 7205, 7206, 7207, 7208, 7209, 7211, 7212, 7213, 7214, 7215, 7216, 7217, 7218, 7219, 7220, 7223, 7225, 7255, 7258, 7261, 7263, 7265, 7267, 7268, 7269, 7271, 7272, 7274, 7276, 7278, 7279, 7281, 7284, 7286, 7288, 7290, 7291, 7312, 7313, 7320, 7321, 7322, 7323, 7325, 7326, 7327, 7328, 7337, 7342, 7343, 7344, 7345, 7346, 7348, 7349, 7351, 7352, 7353, 7354, 7355, 7356, 7357, 7358, 7359, 7361, 7363, 7364, 7365, 7367, 7369, 7374, 7376, 7377, 7378, 7379, 7380, 7381, 7382, 7386, 7387, 7388, 7389, 7390, 7391, 7392, 7393, 7395, 7396, 7397, 7398, 7399, 7400, 7401, 7402, 7403, 7405, 7406, 7408, 7409, 7410, 7411, 7412, 7413, 7414, 7415, 7417, 7418, 7419, 7421, 7422, 7423, 7425, 7427, 7428, 7429, 7431, 7432, 7433, 7434, 7435, 7437, 7439, 7440, 7441, 7442, 7444, 7445, 7446, 7447, 7449, 7451, 7452, 7453, 7454, 7455, 7456, 7457, 7458, 7460, 7461, 7469, 7471, 7472, 7473, 7474, 7477, 7478, 7480, 7481, 7482, 7483, 7484, 7487, 7488, 7489, 7490, 7491, 7492, 7493, 7495, 7497, 7498, 7500, 7501, 7502, 7513, 7514, 7515, 7517, 7519, 7520, 7521, 7528, 7529, 7530, 7531, 7532, 7533, 7534, 7537, 7538, 7539, 7542, 7543, 7545, 7546, 7547, 7549, 7550, 7551, 7552, 7553, 7554, 7564, 7565, 7566, 7567, 7569, 7570, 7572, 7573, 7574, 7575, 7577, 7579, 7581, 7584, 7585, 7586, 7587, 7590, 7591, 7592, 7596, 7597, 7599, 7600, 7604, 7605, 7606, 7607, 7608, 7609, 7611, 7612, 7613, 7622, 7623, 7625, 7626, 7627, 7628, 7630, 7631, 7633, 7634, 7635, 7636, 7638, 7639, 7640, 7641, 7643, 7644, 7645, 7646, 7647, 7649, 7655, 7656, 7657, 7659, 7660, 7662, 7663, 7664, 7665, 7666, 7668, 7669, 7670, 7671, 7673, 7674, 7675, 7676, 7677, 7678, 7679, 7681, 7682, 7684, 7685, 7687, 7688, 7689, 7690, 7692, 7693, 7694, 7696, 7698, 7700, 7701, 7702, 7703, 7704, 7705, 7710, 7712, 7713, 7714, 7715, 7716, 7717, 7718, 7721, 7723, 7724, 7725, 7726, 7727, 7728, 7729, 7730, 7732, 7733, 7734, 7735, 7736, 7737, 7751, 7752, 7753, 7755, 7756, 7761, 7762, 7763, 7764, 7766, 7767, 7768, 7769, 7771, 7772, 7775, 7776, 7778, 7779, 7780, 7781, 7783, 7784, 7785, 7786, 7788, 7789, 7790, 7791, 7794, 7795, 7796, 7797, 7798, 7801, 7802, 7803, 7814, 7816, 7819, 7821, 7824, 7825, 7828, 7829, 7830, 7831, 7832, 7834, 7835, 7836, 7842, 7843, 7845, 7846, 7847, 7848, 7850, 7851, 7852, 7853, 7855, 7856, 7857, 7858, 7859, 7865, 7868, 7869, 7870, 7871, 7873, 7874, 7876, 7877, 7878, 7880, 7882, 7883, 7884, 7885, 7886, 7887, 7890, 7891, 7892, 7893, 7894, 7896, 7897, 7898, 7899, 7900, 7901, 7903, 7904, 7905, 7906, 7907, 7908, 7910, 7911, 7912, 7913, 7914, 7931, 7932, 7933, 7935, 7937, 7938, 7939, 7941, 7942, 7944, 7948, 7949, 7950, 7952, 7954, 7955, 7957, 7958, 7960, 7961, 7962, 7963, 7965, 7966, 7967, 7975, 7977, 7978, 7979, 7981, 7982, 7986, 7988, 7989, 7990, 7991, 7992, 7995, 7996, 7997, 7998, 7999, 8000, 8001, 8002, 8004, 8005, 8006, 8007, 8008, 8009, 8013, 8014, 8015, 8016, 8018, 8020, 8021, 8022, 8023, 8026, 8027, 8028, 8029, 8030, 8031, 8032, 8033, 8038, 8039, 8040, 8041, 8042, 8044, 8045, 8046, 8047, 8048, 8050, 8051, 8052, 8053, 8055, 8056, 8057, 8058, 8059, 8061, 8063, 8064, 8065, 8066, 8068, 8069, 8070, 8071, 8072, 8074, 8075, 8076, 8077, 8079, 8080, 8081, 8082, 8084, 8086, 8087, 8088, 8089, 8090, 8097, 8100, 8101, 8102, 8104, 8106, 8108, 8109, 8110, 8111, 8112, 8113, 8115, 8116, 8117, 8118, 8119, 8121, 8122, 8123, 8124, 8129, 8130, 8131, 8133, 8137, 8138, 8139, 8140, 8141, 8142, 8144, 8156, 8157, 8158, 8159, 8160, 8161, 8162, 8165, 8166, 8167, 8168, 8169, 8170, 8173, 8174, 8175, 8176, 8177, 8180, 8183, 8184, 8186, 8188, 8190, 8191, 8192, 8193, 8195, 8196, 8198, 8199, 8202, 8203, 8204, 8205, 8206, 8207, 8208, 8209, 8210, 8211, 8212, 8213, 8217, 8219, 8220, 8221, 8228, 8229, 8230, 8231, 8235, 8236, 8238, 8239, 8240, 8242, 8244, 8245, 8246, 8247, 8248, 8250, 8251, 8253, 8254, 8255, 8256, 8257, 8258, 8265, 8266, 8267, 8268, 8270, 8271, 8273, 8274, 8275, 8276, 8277, 8278, 8279, 8281, 8284, 8286, 8287, 8288, 8289, 8290, 8292, 8293, 8294, 8295, 8297, 8298, 8299, 8300, 8305, 8306};
    std::vector<int> bad_runs	  {6968, 6970, 6972, 6973, 6975, 6976, 6977, 6978, 6979, 6980, 6981, 6982, 6983, 6984, 7313, 7326, 7415, 7417, 7435,7469, 7471, 7517, 7519, 7520, 7537, 7538, 7542, 7543, 7545, 7546, 7547, 7573, 7575, 7657, 7659, 7679, 7681, 7843, 7847, 7848, 7850, 7851, 7852, 7853, 7855, 7856, 7857, 7858, 7859, 7865, 7868, 7869, 7907, 7910, 7932, 7933, 7935, 7937, 7954, 7955, 8018, 8028, 8029,8030, 8031, 8032, 8033, 8115, 8121, 8167, 8201, 8202, 8203, 8204, 8205, 8206, 8207, 8208, 8209, 8210, 8211, 8212, 8213, 8215, 8284, 8281, 8289};

	auto PhysRunCut = [&physical_runs]( UInt_t run_id ){return std::find( physical_runs.begin(), physical_runs.end(), run_id) != physical_runs.end();};

    // Vtx Z cirrections functions
    auto vtx_correction_generator = []( TGraphErrors* g1_calib ){
	   return [g1_calib](double _vtx, UInt_t _runId){
	   	return _vtx - (g1_calib!=nullptr ? g1_calib->Eval( static_cast<double>(_runId)) : 0.);
	   };
	};
	//RefMult
	auto ref_mult_corr = []( TGraphErrors* g1_calib ){
	   return [g1_calib](unsigned long _mult, UInt_t _runId){
	   	return (_mult * (g1_calib!=nullptr ? g1_calib->Eval( static_cast<double>(_runId)) : 1. ));
	   };
	};

	//Для физ анализа спользовать раны с 6925 (runId > 6925). До этого проблемы с gem.
	std::cout<<physical_runs.size()<<"\t"<<bad_runs.size()<<"\n";
  	RemoveBadRuns(physical_runs,bad_runs);
  	std::cout<<physical_runs.size()<<"\n";

    //Отбор событий: 1-прошли отбор, 0-нет.
    auto passSel = [](UInt_t runId, bool PhysRun, UInt_t CCT2, int noPileup, float vtxChi2Ndf, int vtxNtracks, double vtxRcorr, double vtxZcorr) -> bool
    {
        if (runId >= 8312) return false;
        if (runId <= 7050) return false;
        if (!PhysRun) return false;
        if (CCT2==0) return false;
        if (!(vtxChi2Ndf > std::numeric_limits<float>::min())) return false;
        if (!(vtxNtracks >= 2)) return false;
        if (!(vtxRcorr < 1.)) return false;
        if (!(std::abs(vtxZcorr) < 1.)) return false;
        if (noPileup != 1) return false;
        return true;
    };

    ROOT::RDataFrame d(*chain); 
	auto data = d
        .Define("runId",[run_id](){ return run_id; }, {} )
        .Define("evtId","DstEventHeader.fEventId")
        //trigger mask
        .Define("triggerMapBR","BmnTrigInfo.fInputsBR")
        .Define("triggerMapAR","BmnTrigInfo.fInputsAR")
        // pileup by Oleg Golosov
        .Define("centralHitIndex", CentralHitIndexBC1S,{"BmnBC1hitInfo."})
        .Define("centralHitClass", getValue, {"BmnBC1hitInfo.fBC1hitClasses","centralHitIndex"})
        .Define("k1central",  BC1hitClasses,{"BmnBC1hitInfo.","centralHitIndex"})
        .Define("k0closest",  ClosestBC1hitsDt(BmnEventClass::k0), {"BmnBC1hitInfo."})
        .Define("kV0closest", ClosestBC1hitsDt(BmnEventClass::kV0),{"BmnBC1hitInfo."})
        .Define("k1closest",  ClosestBC1hitsDt(BmnEventClass::k1), {"BmnBC1hitInfo."})
        .Define("kV1closest", ClosestBC1hitsDt(BmnEventClass::kV1),{"BmnBC1hitInfo."})
        .Define("k0cut", "k0closest.at(0)<-200 && k0closest.at(1)>200") 
        .Define("kV0cut", "kV0closest.at(0)<-200 && kV0closest.at(1)>200") 
        .Define("k1cut", "k1closest.at(0)<-1600 && k1closest.at(1)>600") 
        .Define("kV1cut", "kV1closest.at(0)<-1800 && kV1closest.at(1)>600")
        .Define("noPileup", "if(k1central==1 && k0cut && kV0cut && k1cut && kV1cut) return 1; else if(k1central==-1) return -1; else return 0;")
        //Vertex Corr (Run-by-Run)
        .Define("vtxNtracks","MpdVertex.fNTracks")
        .Define("vtxChi2","(Float_t)MpdVertex.fChi2")
        .Define("vtxNdf","MpdVertex.fNDF")
        .Define("vtxChi2Ndf","return (Float_t)MpdVertex.fChi2/MpdVertex.fNDF")
        .Define("vtxXcorr", vtx_correction_generator(g1_FitVtxX), {"MpdVertex.fX","runId"})
        .Define("vtxYcorr", vtx_correction_generator(g1_FitVtxY), {"MpdVertex.fY","runId"})
        .Define("vtxZcorr", vtx_correction_generator(g1_FitVtxZ), {"MpdVertex.fZ","runId"})
        .Define("vtxRcorr","return sqrt(vtxXcorr*vtxXcorr + vtxYcorr*vtxYcorr);" )
        //Track parametrs
        .Define("trDcaX",recDcaXYZR(DcaTypeConvert::X),{"BmnGlobalTrack","MpdVertex."})
        .Define("trDcaY",recDcaXYZR(DcaTypeConvert::Y),{"BmnGlobalTrack","MpdVertex."})
        .Define("trDcaR",recDcaXYZR(DcaTypeConvert::R),{"BmnGlobalTrack","MpdVertex."})
        .Define("trPt",GlobalTrackParamFloat(TrackParamConvert::Pt),{"BmnGlobalTrack"})
        .Define("trEta",GlobalTrackParamFloat(TrackParamConvert::Eta),{"BmnGlobalTrack"})
        //RefMult and Centrality
        .Define("RefMult", RefMult_gt(_ptMin,_ptMax,_etaMin,_etaMax,_dcaR),{"trPt","trEta","trDcaR"}) //0.05<pt<2 && 0.7<eta<2.7 && dca_R<1
		// Отбор событий        
        .Define("PhysRunCut",PhysRunCut,{"runId"} ) //PhysRuns
		.Define("CCT2_mask",CCT2_mask, {"triggerMapAR"}) //CCT
        .Define("passSel",passSel,{"runId","PhysRunCut","CCT2_mask","noPileup","vtxChi2Ndf","vtxNtracks","vtxRcorr","vtxZcorr"})
		.Define("RefMultCorr", [ref_mult_corr, f1_FitRunIdFactor_1](unsigned long refMult, UInt_t runId, bool pass) -> float {
            return pass ? ref_mult_corr(f1_FitRunIdFactor_1)(refMult, runId) : -1.0f;},
            {"RefMult", "runId", "passSel"})
        .Define("centrality_Bin", [centrality_function](float refMultCorr, bool pass) -> int {
            return pass ? centrality_function(1)(refMultCorr) : -1;},
            {"RefMultCorr", "passSel"})
	;

    std::vector<std::string> IncludeToWrite = {
        "runId",
        "evtId",
        "vtxXcorr",
        "vtxYcorr",
        "vtxZcorr",
        "noPileup",
        "PhysRunCut",
        "CCT2_mask",
        "passSel",
        "RefMult",
        "RefMultCorr",
        "centrality_Bin"
    };
    
    TString tstr(in_file_name);
    TString base = gSystem->BaseName(tstr);
    std::string out_file_name = outPath + "/cent_" + std::string(base.Data());
    
    data.Snapshot(TreeName, out_file_name, IncludeToWrite);

	file_fit_VtxXYZ->Close();
	file_RefMultCorr->Close();

	std::cout<<"cent_completed!"<<std::endl;
	std::cout<<"file:\t"<<out_file_name<<std::endl;

}


void UsageExample(
    std::string in_file_name_bmndst="/eos/nica/bmn/exp/dst/run8/25.09.0/mpd_run_Top_8130_ev0_p16.root",
    std::string in_file_name_cent="./cent_mpd_run_Top_8130_ev0_p16.root")
{

    TChain *chain=makeChain(in_file_name_bmndst, "bmndata");
    TChain *chainCent=makeChain(in_file_name_cent, "bmndata");
    chain->AddFriend(chainCent);

    //Далее работать с chain
    // Ветка centrality_Bin:
    //   0 :  0-10%
    //   1 : 10-20%
    //   2 : 20-30%
    //   3 : 30-40%
    //   4 : 40-50%
    //   5 : 50-60%
    //  -1 : все остальное (>60% cent, или события, для которых центральность не определена, или события не прошедшие критерии отбора 
    //       (CCT2,|VtxZ|<1,VtxR<1,MpdVertex.fNTracks>1, pileup cut, ))

    ROOT::RDataFrame dd(*chain);
    
    auto histoPtr = dd.Filter("passSel==1").Histo1D({"h1_Cent_10", ";cent bin", 11, -1.5, 9.5}, "centrality_Bin");

	std::vector<std::string> labels = {"err","0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%","null"};
    for (size_t i = 0; i < labels.size(); ++i) {
		histoPtr->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
	}

    TCanvas *c = new TCanvas("c", "c", 800, 600);
    histoPtr->Draw();
    c->SaveAs("./h1_Cent_10.png");

}