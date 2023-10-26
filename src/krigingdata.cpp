#include "clusterclasses.h"
#include <boost/random/uniform_int_distribution.hpp>

void ClusterApp::krigingData::generate_grid() {
	int n_cl_range = upper_int[0] - lower_int[0];
	int n_ind_range = upper_int[1] - lower_int[1];
	int n_cl_step = n_cl_range / 20;
	int n_ind_step = n_ind_range / 20;
	for (int i = 0; i < 20; i++) {
		n_cl_grid[i] = lower_int[0] + i * n_cl_step;
		n_ind_grid[i] = lower_int[1] + i * n_ind_step;
		if (i % 2 == 0) {
			delete[] n_cl_grid_label[i];
			std::string label = std::format("{:.0f}", n_cl_grid[i]);
			n_cl_grid_label[i] = new char[label.length() + 1];
			strcpy(n_cl_grid_label[i], label.c_str());
			delete[] n_ind_grid_label[i];
			std::string label2 = std::format("{:.0f}", n_ind_grid[i]);
			n_ind_grid_label[i] = new char[label2.length() + 1];
			strcpy(n_ind_grid_label[i], label2.c_str());
		}
		else {
			n_cl_grid_label[i] = "";
			n_ind_grid_label[i] = "";
		}		
	}
}

void ClusterApp::krigingData::new_sample(int n) {
	n_ind.clear();
	n_cl.clear();
	power.clear();
	std::time_t now = std::time(0);
	boost::random::mt19937 rng{ static_cast<std::uint32_t>(now) };
	boost::random::uniform_int_distribution<> r_cl(lower_int[0], upper_int[0]);
	boost::random::uniform_int_distribution<> r_ind(lower_int[1], upper_int[1]);
	for (int i = 0; i < n; i++) {
		n_ind.push_back(r_ind(rng));
		n_cl.push_back(r_cl(rng));
	}
}

void ClusterApp::krigingData::generate_data() {
	int m = n_ind.size();
	int n = m - power.size();
	int start = power.size();
	if (n > 0) {
		std::vector<float> n_weights;
		std::vector<int> cl_vector;
		std::vector<int> ind_vector;
		ClusterApp::modelSummary summary(glmm.designs);
		// get weights
		n_weights.resize(glmm.designs.sequences);
		int J = glmm.designs.total_clusters();
		for (int i = 0; i < glmm.designs.sequences; i++) {
			n_weights[i] = (float)(*glmm.designs.n_clusters(i)) / (float)J;
		}

		cl_vector.resize(glmm.designs.sequences);
		for (int i = 0; i < glmm.designs.sequences; i++) {
			cl_vector[i] = *glmm.designs.n_clusters(i);
			for (int t = 0; t < glmm.designs.time; t++) {
				ind_vector.push_back(*glmm.designs.n(i, t));
			}
		}
		// now create new totals 
		std::vector<int> new_n_cl(n);
		for (int i = 0; i < n; i++) {
			new_n_cl = glmm.round_weights(n_weights, n_cl[i + start]);
			for (int k = 0; k < glmm.designs.sequences; k++) {
				*glmm.designs.n_clusters(k) = new_n_cl[k];
				for (int t = 0; t < glmm.designs.time; t++) {
					*glmm.designs.n(k, t) = n_ind[i + start];
				}
			} 
			glmm.update_model_data(updater.generate_data());
			glmm.power(summary);
			power.push_back(summary.power);
		}
		int counter = 0;
		for (int i = 0; i < glmm.designs.sequences; i++) {
			*glmm.designs.n_clusters(i) = cl_vector[i];
			for (int t = 0; t < glmm.designs.time; t++) {
				*glmm.designs.n(i, t) = ind_vector[counter];
				counter++;
			}
		} 
		initialised = true;
	}
	

}

void ClusterApp::krigingData::update(bool resample) {
	updating = true;
	Eigen::MatrixXf C(n_ind.size(),n_ind.size());
	Eigen::VectorXf Cx(n_ind.size());
	int n_cl_range = upper_int[0] - lower_int[0];
	int n_ind_range = upper_int[1] - lower_int[1];
	float dif0, dif1, dif2;
	for (int i = 0; i < n_ind.size(); i++) {
		C(i, i) = 1.0f;
	}
	for (int i = 0; i < n_ind.size()-1; i++) {
		for (int j = i; j < n_ind.size(); j++) {
			dif0 = (n_ind[i] - n_ind[j]) / (float)n_ind_range;
			dif1 = (n_cl[i] - n_cl[j]) / (float)n_cl_range;
			dif2 = exp(-1.0*sqrt(dif0 * dif0 + dif1 * dif1)/bandwidth);
			C(i, j) = dif2;
			C(j, i) = dif2;
;		}
	}
	
	Eigen::VectorXf ones = Eigen::VectorXf::Ones(n_ind.size());
	Eigen::MatrixXf Cinv = C.llt().solve(Eigen::MatrixXf::Identity(n_ind.size(), n_ind.size()));
	float v = ones.transpose() * Cinv * ones;
	if (power.size() != n_ind.size()) generate_data();
	Eigen::VectorXf y = Eigen::Map<Eigen::VectorXf>(power.data(), power.size());
	y.array() *= 0.01;
	mu = (1 / v) * ones.transpose() * Cinv * y;
	y.array() -= mu;
	if (resample) {
		std::pair<int, int> new_n;
		std::pair<int, int> keep_n;
		std::time_t now = std::time(0);
		boost::random::mt19937 rng{ static_cast<std::uint32_t>(now) };
		boost::random::uniform_int_distribution<> r_cl(lower_int[0], upper_int[0]);
		boost::random::uniform_int_distribution<> r_ind(lower_int[1], upper_int[1]);
		float best_tmse = 10;
		float tmse, v, new_m, sd_mean;
		for (int i = 0; i < resample_total; i++) {
			new_n.first = r_cl(rng);
			new_n.second = r_ind(rng);
			for (int j = 0; j < n_ind.size(); j++) {
				dif0 = (new_n.second - n_ind[j]) / (float)n_ind_range;
				dif1 = (new_n.first - n_cl[j]) / (float)n_cl_range;
				dif2 = exp(-1.0 * sqrt(dif0 * dif0 + dif1 * dif1) / bandwidth);
				Cx(j) = dif2;
			}
			new_m = mu + Cx.transpose() * Cinv * y;
			v = sqrt(1 - Cx.transpose() * Cinv * Cx);
			sd_mean = (new_m - threshold_power) / v;
			tmse = v * (1 / (sqrt(2 * M_PI) * v)) * exp(-0.5 * sd_mean * sd_mean);
			if (tmse < best_tmse) {
				keep_n = new_n;
			}
		}
		n_cl.push_back(keep_n.first);
		n_ind.push_back(keep_n.second);
		generate_data();
	}
	else {		
		for (int i = 0; i < 20; i++) {
			for (int j = 0; j < 20; j++) {
				for (int k = 0; k < n_ind.size(); k++) {
					dif0 = (n_ind_grid[j] - n_ind[k]) / (float)n_ind_range;
					dif1 = (n_cl_grid[i] - n_cl[k]) / (float)n_cl_range;
					dif2 = exp(-1.0 * sqrt(dif0 * dif0 + dif1 * dif1) / bandwidth);
					Cx(k) = dif2;
				}
				surface[j + (19-i) * 20] = (mu + Cx.transpose() * Cinv * y) * 100; // need to include in row major order from the bottom for the plot
			}			
		}
		surface_initialised = true;
	}
	updating = false;
}
