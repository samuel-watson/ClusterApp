#include "clusterclasses.h"

void ClusterApp::plotData::update_data() {

	ClusterApp::modelSummary summary(glmm.designs);
	float lower = x_axis_limits.first;
	float upper = x_axis_limits.second;
	float diff = (upper - lower) / 19;
	std::vector<int> int_vector;
	std::vector<int> int_vector_series;
	float fl_value;
	float fl_value_series = 0;
	std::vector<float> n_weights;

	if (xaxis == ClusterApp::XAxis::clusters || (multiple_series && series == ClusterApp::XAxis::clusters)) {
		// get weights
		n_weights.resize(glmm.designs.sequences);
		int J = glmm.designs.total_clusters();
		for (int i = 0; i < glmm.designs.sequences; i++) {
			n_weights[i] = (float)(*glmm.designs.n_clusters(i)) / (float)J;
		}		
	}

	for (int i = 0; i < 20; i++) {
		x_data[i] = lower + i * diff;
	}

	 //save existing values	for series variables
	if (multiple_series) {
		switch (series) {
		case ClusterApp::XAxis::clusters:
		{
			int_vector_series.resize(glmm.designs.sequences);
			for (int i = 0; i < glmm.designs.sequences; i++) int_vector_series[i] = *glmm.designs.n_clusters(i);
			break;
		}			
		case ClusterApp::XAxis::individual_n:
		{
			int_vector_series.clear();
			for (int i = 0; i < glmm.designs.sequences; i++) {
				for (int t = 0; t < glmm.designs.time; t++) {
					int_vector_series.push_back(*glmm.designs.n(i, t));
				}
			}
			break;
		}
		case ClusterApp::XAxis::icc:
		{
			fl_value_series = glmm.statmodel.ixx_pars[0];
			break;
		}
		case ClusterApp::XAxis::treatment_effect:
		{
			fl_value_series = glmm.statmodel.te_pars[0];
			break;
		}
		case ClusterApp::XAxis::baseline:
		{
			fl_value_series = glmm.statmodel.beta_pars[0];
			break;
		}
		case ClusterApp::XAxis::cac:
		{
			fl_value_series = glmm.statmodel.ixx_pars[1];
			break;
		}
			
		}
	}

	int n_times = multiple_series ? n_series : 1;

	for (int j = 0; j < n_times; j++) {

		// do series variables
		if (multiple_series) {
			switch (series) {
			case ClusterApp::XAxis::clusters:
			{
				std::vector<int> new_n_total = glmm.round_weights(n_weights, (int)x_series[j]);
				for (int k = 0; k < glmm.designs.sequences; k++) *glmm.designs.n_clusters(k) = new_n_total[k];
				break;
			}
			case ClusterApp::XAxis::individual_n:
			{
				for (int k = 0; k < glmm.designs.sequences; k++) {
					for (int t = 0; t < glmm.designs.time; t++) {
						*glmm.designs.n(k, t) = (int)x_series[j];
					}
				}
				break;
			}
			case ClusterApp::XAxis::icc:
			{
				glmm.statmodel.ixx_pars[0] = x_series[j];
				break;
			}
			case ClusterApp::XAxis::treatment_effect:
			{
				glmm.statmodel.te_pars[0] = x_series[j];
				break;
			}
			case ClusterApp::XAxis::baseline:
			{
				glmm.statmodel.beta_pars[0] = x_series[j];
				break;
			}
			case ClusterApp::XAxis::cac:
			{
				glmm.statmodel.ixx_pars[1] = x_series[j];
				break;
			}

			}
		}


		switch (xaxis) {
		case ClusterApp::XAxis::clusters:
		{
			int_vector.resize(glmm.designs.sequences);
			for (int i = 0; i < glmm.designs.sequences; i++) int_vector[i] = *glmm.designs.n_clusters(i);
			// now create new totals 
			std::vector<int> new_n_total(n_data_points);
			for (int i = 0; i < n_data_points; i++) {
				new_n_total = glmm.round_weights(n_weights, (int)x_data[i]);
				for (int k = 0; k < glmm.designs.sequences; k++) *glmm.designs.n_clusters(k) = new_n_total[k];
				glmm.update_model_data(updater.generate_data());
				extract_y(summary, i, j);
			}
			for (int i = 0; i < glmm.designs.sequences; i++) *glmm.designs.n_clusters(i) = int_vector[i];
			break;
		}
		case ClusterApp::XAxis::individual_n:
		{
			int_vector.clear();
			for (int i = 0; i < glmm.designs.sequences; i++) {
				for (int t = 0; t < glmm.designs.time; t++) {
					int_vector.push_back(*glmm.designs.n(i, t));
				}
			}
			for (int i = 0; i < n_data_points; i++) {
				for (int k = 0; k < glmm.designs.sequences; k++) {
					for (int t = 0; t < glmm.designs.time; t++) {
						*glmm.designs.n(k, t) = (int)x_data[i];
					}
				}
				glmm.update_parameters();
				extract_y(summary, i, j);
			}			
			int counter = 0;
			for (int i = 0; i < glmm.designs.sequences; i++) {
				for (int t = 0; t < glmm.designs.time; t++) {
					*glmm.designs.n(i, t) = int_vector[counter];
					counter++;
				}
			}
			break;
		}
		case ClusterApp::XAxis::icc:
		{
			fl_value = glmm.statmodel.ixx_pars[0];
			for (int i = 0; i < n_data_points; i++) {
				glmm.statmodel.ixx_pars[0] = x_data[i];
				glmm.update_parameters();
				extract_y(summary,i,j);
			}
			glmm.statmodel.ixx_pars[0] = fl_value;
			break;
		}
		case ClusterApp::XAxis::treatment_effect:
		{
			fl_value = glmm.statmodel.te_pars[0];
			for (int i = 0; i < n_data_points; i++) {
				glmm.statmodel.te_pars[0] = x_data[i];
				glmm.update_parameters();
				extract_y(summary, i, j);
			}
			glmm.statmodel.te_pars[0] = fl_value;
			break;
		}
		case ClusterApp::XAxis::baseline:
		{
			fl_value = glmm.statmodel.beta_pars[0];
			for (int i = 0; i < n_data_points; i++) {
				glmm.statmodel.beta_pars[0] = x_data[i];
				glmm.update_parameters();
				extract_y(summary, i, j);
			}
			glmm.statmodel.beta_pars[0] = fl_value;
			break;
		}
		case ClusterApp::XAxis::cac:
		{
			fl_value = glmm.statmodel.ixx_pars[1];
			for (int i = 0; i < n_data_points; i++) {
				glmm.statmodel.ixx_pars[1] = x_data[i];
				glmm.update_parameters();
				extract_y(summary, i, j);
			}
			glmm.statmodel.ixx_pars[1] = fl_value;
			break;
		}

		}
	}

	if (multiple_series) {
		switch (series) {
		case ClusterApp::XAxis::clusters:
		{
			for (int i = 0; i < glmm.designs.sequences; i++) *glmm.designs.n_clusters(i) = int_vector_series[i];
			break;
		}
		case ClusterApp::XAxis::individual_n:
		{
			int counter = 0;
			for (int i = 0; i < glmm.designs.sequences; i++) {
				for (int t = 0; t < glmm.designs.time; t++) {
					*glmm.designs.n(i, t) = int_vector_series[counter];
					counter++;
				}
			}
			break;
		}
		case ClusterApp::XAxis::icc:
		{
			glmm.statmodel.ixx_pars[0] = fl_value_series;
			break;
		}
		case ClusterApp::XAxis::treatment_effect:
		{
			glmm.statmodel.te_pars[0] = fl_value_series;
			break;
		}
		case ClusterApp::XAxis::baseline:
		{
			glmm.statmodel.beta_pars[0] = fl_value_series;
			break;
		}
		case ClusterApp::XAxis::cac:
		{
			glmm.statmodel.ixx_pars[1] = fl_value_series;
			break;
		}

		}
	}

	updater.update_data();
	y_axis_limits.first = min_y();
	y_axis_limits.second = max_y();
	initialised = true;
}

void ClusterApp::plotData::extract_y(ClusterApp::modelSummary& summary, int i, int s) {
	switch (yaxis) {
	case ClusterApp::YAxis::power:
	{
		glmm.power(summary);
		switch (s) {
		case 0:
			y_data_1[i] = summary.power;
			break;
		case 1:
			y_data_2[i] = summary.power;
			break;
		case 2:
			y_data_3[i] = summary.power;
			break;
		}
		break;
	}		
	case ClusterApp::YAxis::ci_width:
	{
		glmm.power(summary);
		switch (s) {
		case 0:
			y_data_1[i] = summary.ci_width;
			break;
		case 1:
			y_data_2[i] = summary.ci_width;
			break;
		case 2:
			y_data_3[i] = summary.ci_width;
			break;
		}
		break;
	}
	case ClusterApp::YAxis::power_bw:
	{
		glmm.power(summary);
		switch (s) {
		case 0:
			y_data_1[i] = summary.power_bw;
			break;
		case 1:
			y_data_2[i] = summary.power_bw;
			break;
		case 2:
			y_data_3[i] = summary.power_bw;
			break;
		}
		break;
	}
	case ClusterApp::YAxis::ci_width_bw:
	{
		glmm.power(summary);
		switch (s) {
		case 0:
			y_data_1[i] = summary.ci_width_bw;
			break;
		case 1:
			y_data_2[i] = summary.ci_width_bw;
			break;
		case 2:
			y_data_3[i] = summary.ci_width_bw;
			break;
		}
		break;
	}
	case ClusterApp::YAxis::power_kr:
	{
		glmm.power(summary);
		switch (s) {
		case 0:
			y_data_1[i] = summary.power_kr;
			break;
		case 1:
			y_data_2[i] = summary.power_kr;
			break;
		case 2:
			y_data_3[i] = summary.power_kr;
			break;
		}
		break;
	}
	case ClusterApp::YAxis::ci_width_kr:
	{
		glmm.power(summary);
		switch (s) {
		case 0:
			y_data_1[i] = summary.ci_width_kr;
			break;
		case 1:
			y_data_2[i] = summary.ci_width_kr;
			break;
		case 2:
			y_data_3[i] = summary.ci_width_kr;
			break;
		}
		break;
	}
	case ClusterApp::YAxis::power_de:
	{
		glmm.power(summary);
		switch (s) {
		case 0:
			y_data_1[i] = summary.power_de;
			break;
		case 1:
			y_data_2[i] = summary.power_de;
			break;
		case 2:
			y_data_3[i] = summary.power_de;
			break;
		}
		break;
	}
	case ClusterApp::YAxis::ci_width_de:
	{
		glmm.power(summary);
		switch (s) {
		case 0:
			y_data_1[i] = summary.ci_width_de;
			break;
		case 1:
			y_data_2[i] = summary.ci_width_de;
			break;
		case 2:
			y_data_3[i] = summary.ci_width_de;
			break;
		}
		break;
	}
	}

}

bool ClusterApp::plotData::check() {
	if (initialised) {
		CRC crc;
		crc(static_cast<std::underlying_type_t<ClusterApp::XAxis>>(xaxis));
		crc(static_cast<std::underlying_type_t<ClusterApp::YAxis>>(yaxis));
		crc(static_cast<std::underlying_type_t<ClusterApp::XAxis>>(series));
		crc(x_axis_limits.first);
		crc(x_axis_limits.second);
		crc(n_series);
		crc(multiple_series);
		bool has_changed = crc_val != crc.get();
		if (has_changed)crc_val = crc.get();
		return has_changed;
	}
	else {
		return true;
	}	
}

float ClusterApp::plotData::max_y() {
	float max = 0;
	for (int i = 0; i < n_data_points; i++) {
		if (y_data_1[i] > max)max = y_data_1[i];
		if (multiple_series) {
			if (n_series > 1 && y_data_2[i] > max)max = y_data_2[i];
			if (n_series == 3 && y_data_3[i] > max)max = y_data_3[i];
		}
	}
	return max;
}

float ClusterApp::plotData::min_y() {
	float min = 1000;
	for (int i = 0; i < n_data_points; i++) {
		if (y_data_1[i] < min)min = y_data_1[i];
		if (multiple_series) {
			if (n_series > 1 && y_data_2[i] < min)min = y_data_2[i];
			if (n_series == 3 && y_data_3[i] < min)min = y_data_3[i];
		}
	}
	return min;
}