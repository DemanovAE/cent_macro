#include <string>
#include <iostream>
#include <fstream> 

using namespace ROOT;
using namespace ROOT::Math;
using namespace ROOT::RDF;

void RemoveBadRuns(std::vector<int> &_Runs, std::vector<int> _BadRUns){
	for(auto _run : _BadRUns){
		_Runs.erase(remove(_Runs.begin(), _Runs.end(), _run), _Runs.end());
	}
}

bool IsGoodFile(const char* filename) {
    std::unique_ptr<TFile> f(TFile::Open(filename, "READ"));
    if (!f || f->IsZombie() || f->TestBit(TFile::kRecovered)) {
        return false;
    }

    TTree* tree = (TTree*)f->Get("t"); // "bmndata"
    if (!tree) {
        return false;  // Дерева нет
    }
    
    Long64_t nEntries = tree->GetEntries();
    if (nEntries <= 0) {
        return false;  // Дерево пустое
    }
    
    return true;
}


void test_cent_new(std::string str_in_list, std::string out_file_name="qa.root", std::string in_file_vtx="./run8_25.09_corr_VtxXYZ.root", std::string in_file_RefMultCorr="./run8_25.09_corrections_DifRefRunId.root"){	

	//CorrRunId.root

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

	std::vector<int> physical_runs{6667, 6668, 6669, 6670, 6671, 6672, 6673, 6674, 6675, 6676, 6677, 6678, 6679, 6680, 6681, 6683, 6684, 6685, 6686, 6687, 6689, 6690, 6691, 6692, 6694, 6695, 6696, 6698, 6699, 6732, 6733, 6734, 6737, 6738, 6739, 6740, 6745, 6752, 6753, 6760, 6761, 6765, 6766, 6767, 6768, 6769, 6771, 6772, 6773, 6774, 6779, 6780, 6782, 6783, 6785, 6786, 6788, 6794, 6795, 6797, 6799, 6800, 6803, 6815, 6816, 6817, 6818, 6819, 6820, 6821, 6822, 6879, 6882, 6883, 6884, 6886, 6887, 6889, 6891, 6900, 6901, 6902, 6903, 6904, 6905, 6906, 6907, 6908, 6909, 6910, 6911, 6915, 6916, 6918, 6919, 6920, 6921, 6923, 6924, 6926, 6927, 6928, 6929, 6930, 6931, 6932, 6933, 6934, 6935, 6936, 6937, 6939, 6940, 6968, 6970, 6972, 6973, 6975, 6976, 6977, 6978, 6979, 6980, 6981, 6982, 6983, 6984, 6990, 6991, 6992, 6993, 6994, 6995, 6997, 6998, 6999, 7000, 7002, 7003, 7004, 7005, 7006, 7008, 7009, 7010, 7011, 7012, 7030, 7031, 7032, 7033, 7034, 7035, 7037, 7038, 7040, 7041, 7042, 7043, 7044, 7046, 7047, 7048, 7049, 7050, 7051, 7052, 7053, 7054, 7055, 7056, 7075, 7076, 7077, 7078, 7081, 7082, 7083, 7084, 7086, 7087, 7091, 7092, 7093, 7094, 7096, 7097, 7098, 7100, 7101, 7102, 7103, 7104, 7125, 7126, 7127, 7128, 7129, 7130, 7131, 7132, 7133, 7135, 7136, 7137, 7138, 7146, 7149, 7150, 7151, 7154, 7155, 7156, 7157, 7159, 7160, 7161, 7162, 7163, 7164, 7165, 7166, 7167, 7168, 7173, 7174, 7175, 7176, 7177, 7178, 7179, 7180, 7181, 7182, 7184, 7186, 7187, 7188, 7191, 7192, 7193, 7194, 7195, 7200, 7202, 7203, 7205, 7206, 7207, 7208, 7209, 7211, 7212, 7213, 7214, 7215, 7216, 7217, 7218, 7219, 7220, 7223, 7225, 7255, 7258, 7261, 7263, 7265, 7267, 7268, 7269, 7271, 7272, 7274, 7276, 7278, 7279, 7281, 7284, 7286, 7288, 7290, 7291, 7312, 7313, 7320, 7321, 7322, 7323, 7325, 7326, 7327, 7328, 7337, 7342, 7343, 7344, 7345, 7346, 7348, 7349, 7351, 7352, 7353, 7354, 7355, 7356, 7357, 7358, 7359, 7361, 7363, 7364, 7365, 7367, 7369, 7374, 7376, 7377, 7378, 7379, 7380, 7381, 7382, 7386, 7387, 7388, 7389, 7390, 7391, 7392, 7393, 7395, 7396, 7397, 7398, 7399, 7400, 7401, 7402, 7403, 7405, 7406, 7408, 7409, 7410, 7411, 7412, 7413, 7414, 7415, 7417, 7418, 7419, 7421, 7422, 7423, 7425, 7427, 7428, 7429, 7431, 7432, 7433, 7434, 7435, 7437, 7439, 7440, 7441, 7442, 7444, 7445, 7446, 7447, 7449, 7451, 7452, 7453, 7454, 7455, 7456, 7457, 7458, 7460, 7461, 7469, 7471, 7472, 7473, 7474, 7477, 7478, 7480, 7481, 7482, 7483, 7484, 7487, 7488, 7489, 7490, 7491, 7492, 7493, 7495, 7497, 7498, 7500, 7501, 7502, 7513, 7514, 7515, 7517, 7519, 7520, 7521, 7528, 7529, 7530, 7531, 7532, 7533, 7534, 7537, 7538, 7539, 7542, 7543, 7545, 7546, 7547, 7549, 7550, 7551, 7552, 7553, 7554, 7564, 7565, 7566, 7567, 7569, 7570, 7572, 7573, 7574, 7575, 7577, 7579, 7581, 7584, 7585, 7586, 7587, 7590, 7591, 7592, 7596, 7597, 7599, 7600, 7604, 7605, 7606, 7607, 7608, 7609, 7611, 7612, 7613, 7622, 7623, 7625, 7626, 7627, 7628, 7630, 7631, 7633, 7634, 7635, 7636, 7638, 7639, 7640, 7641, 7643, 7644, 7645, 7646, 7647, 7649, 7655, 7656, 7657, 7659, 7660, 7662, 7663, 7664, 7665, 7666, 7668, 7669, 7670, 7671, 7673, 7674, 7675, 7676, 7677, 7678, 7679, 7681, 7682, 7684, 7685, 7687, 7688, 7689, 7690, 7692, 7693, 7694, 7696, 7698, 7700, 7701, 7702, 7703, 7704, 7705, 7710, 7712, 7713, 7714, 7715, 7716, 7717, 7718, 7721, 7723, 7724, 7725, 7726, 7727, 7728, 7729, 7730, 7732, 7733, 7734, 7735, 7736, 7737, 7751, 7752, 7753, 7755, 7756, 7761, 7762, 7763, 7764, 7766, 7767, 7768, 7769, 7771, 7772, 7775, 7776, 7778, 7779, 7780, 7781, 7783, 7784, 7785, 7786, 7788, 7789, 7790, 7791, 7794, 7795, 7796, 7797, 7798, 7801, 7802, 7803, 7814, 7816, 7819, 7821, 7824, 7825, 7828, 7829, 7830, 7831, 7832, 7834, 7835, 7836, 7842, 7843, 7845, 7846, 7847, 7848, 7850, 7851, 7852, 7853, 7855, 7856, 7857, 7858, 7859, 7865, 7868, 7869, 7870, 7871, 7873, 7874, 7876, 7877, 7878, 7880, 7882, 7883, 7884, 7885, 7886, 7887, 7890, 7891, 7892, 7893, 7894, 7896, 7897, 7898, 7899, 7900, 7901, 7903, 7904, 7905, 7906, 7907, 7908, 7910, 7911, 7912, 7913, 7914, 7931, 7932, 7933, 7935, 7937, 7938, 7939, 7941, 7942, 7944, 7948, 7949, 7950, 7952, 7954, 7955, 7957, 7958, 7960, 7961, 7962, 7963, 7965, 7966, 7967, 7975, 7977, 7978, 7979, 7981, 7982, 7986, 7988, 7989, 7990, 7991, 7992, 7995, 7996, 7997, 7998, 7999, 8000, 8001, 8002, 8004, 8005, 8006, 8007, 8008, 8009, 8013, 8014, 8015, 8016, 8018, 8020, 8021, 8022, 8023, 8026, 8027, 8028, 8029, 8030, 8031, 8032, 8033, 8038, 8039, 8040, 8041, 8042, 8044, 8045, 8046, 8047, 8048, 8050, 8051, 8052, 8053, 8055, 8056, 8057, 8058, 8059, 8061, 8063, 8064, 8065, 8066, 8068, 8069, 8070, 8071, 8072, 8074, 8075, 8076, 8077, 8079, 8080, 8081, 8082, 8084, 8086, 8087, 8088, 8089, 8090, 8097, 8100, 8101, 8102, 8104, 8106, 8108, 8109, 8110, 8111, 8112, 8113, 8115, 8116, 8117, 8118, 8119, 8121, 8122, 8123, 8124, 8129, 8130, 8131, 8133, 8137, 8138, 8139, 8140, 8141, 8142, 8144, 8156, 8157, 8158, 8159, 8160, 8161, 8162, 8165, 8166, 8167, 8168, 8169, 8170, 8173, 8174, 8175, 8176, 8177, 8180, 8183, 8184, 8186, 8188, 8190, 8191, 8192, 8193, 8195, 8196, 8198, 8199, 8202, 8203, 8204, 8205, 8206, 8207, 8208, 8209, 8210, 8211, 8212, 8213, 8217, 8219, 8220, 8221, 8228, 8229, 8230, 8231, 8235, 8236, 8238, 8239, 8240, 8242, 8244, 8245, 8246, 8247, 8248, 8250, 8251, 8253, 8254, 8255, 8256, 8257, 8258, 8265, 8266, 8267, 8268, 8270, 8271, 8273, 8274, 8275, 8276, 8277, 8278, 8279, 8281, 8284, 8286, 8287, 8288, 8289, 8290, 8292, 8293, 8294, 8295, 8297, 8298, 8299, 8300, 8305, 8306};
	std::vector<int> bad_runs	  {6968, 6970, 6972, 6973, 6975, 6976, 6977, 6978, 6979, 6980, 6981, 6982, 6983, 6984, 7313, 7326, 7415, 7417, 7435,7469, 7471, 7517, 7519, 7520, 7537, 7538, 7542, 7543, 7545, 7546, 7547, 7573, 7575, 7657, 7659, 7679, 7681, 7843, 7847, 7848, 7850, 7851, 7852, 7853, 7855, 7856, 7857, 7858, 7859, 7865, 7868, 7869, 7907, 7910, 7932, 7933, 7935, 7937, 7954, 7955, 8018, 8028, 8029,8030, 8031, 8032, 8033, 8115, 8121, 8167, 8201, 8202, 8203, 8204, 8205, 8206, 8207, 8208, 8209, 8210, 8211, 8212, 8213, 8215, 8284, 8281, 8289};
	
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
			
			bool OutCentBin=false; 
			// если true вывод в формате 0 [0-10%], 1 [10-20%], 3 [20-30%]
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

	auto PhisRunCut = [&physical_runs]( UInt_t run_id ){return std::find( physical_runs.begin(), physical_runs.end(), run_id) != physical_runs.end();};

	auto CCT2_mask = []( ROOT::VecOps::RVec<unsigned int> map ){ return map[0] & (1<<7); }; // CCT2_mask 

    // Vtx Z cirrections functions
    auto vtx_correction_generator = []( TGraphErrors* g1_calib ){
	   return [g1_calib](float _vtx, UInt_t _runId){
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

	// Out File
	auto file_out = TFile::Open( out_file_name.c_str(), "RECREATE" );
	file_out->cd();

    TFileCollection collection("collection", "", str_in_list.c_str());
	auto* chain = new TChain("t");
	// Перебираем файлы вручную и добавляем ТОЛЬКО хорошие
	TIter next(collection.GetList());
	TFileInfo* info = nullptr;
	while ((info = (TFileInfo*)next())) {
		const char* fname = info->GetCurrentUrl()->GetUrl();
		if (IsGoodFile(fname)){
			chain->Add(fname);
		}else{
			std::cout<<"bad file: "<<fname<<std::endl;
		}
	}
	
	ROOT::RDataFrame d(*chain); 

	auto dd = d
		.Filter("runId < 8312")
		.Filter("runId > 7016")
		.Filter(PhisRunCut,{"runId"} ) //PhysRUns
		.Filter(CCT2_mask, {"triggerMapAR"}) //CCCT2
		//New vertex
		.Define("RefMult", RefMult_gt(_ptMin,_ptMax,_etaMin,_etaMax,_dcaR),{"trPt","trEta","trDcaR"}) //0.05<pt<2 && 0.7<eta<2.5 && dca_R<1
		.Define("vtxXcorr", vtx_correction_generator(g1_FitVtxX), {"vtxX","runId"})
		.Define("vtxYcorr", vtx_correction_generator(g1_FitVtxY), {"vtxY","runId"})
		.Define("vtxZcorr", vtx_correction_generator(g1_FitVtxZ), {"vtxZ","runId"})
		.Define("vtxRcorr", "return sqrt(vtxXcorr*vtxXcorr + vtxYcorr*vtxYcorr);" )
		//
		.Define( "vtxChi2Ndf", "return vtxChi2/vtxNdf;" )
		// Cuts 
		.Filter("vtxChi2Ndf > std::numeric_limits<float>::min()") // Chi2/ndf > 0.000000...01
		.Filter("vtxNtracks >= 2")
		.Filter("vtxRcorr < 1.") // sqrt(vtxXcorr*vtxXcorr+vtxYcorr*vtxYcorr)
		.Filter("abs(vtxZcorr) < 1.")
		.Filter("noPileup==1") // Кат на пайлап от Олега Голосова 
		.Define("RefMultCorr", ref_mult_corr( f1_FitRunIdFactor_1 ), {"RefMult", "runId"} )
		.Define("centrality", centrality_function(1), {"RefMultCorr"} )
	;

	vector <RResultPtr<::TH1D >> histo1;;
	histo1.push_back( dd.Histo1D({"h1_Cent", ";cent bin", 10, 0, 100}, "centrality" ));

	file_out->cd();
	std::for_each( histo1.begin(), histo1.end(), []( auto h1 ){ h1->Write(); } );

	file_fit_VtxXYZ->Close();
	file_RefMultCorr->Close();
	file_out->Close();

	std::cout<<"cent_completed!"<<std::endl;
	std::cout<<"file:\t"<<out_file_name<<std::endl;

}