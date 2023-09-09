#include "clusterclasses.h"

void ClusterApp::plotData::update_data() {
	double y = 1.0;
	(void)y;

	//ClusterApp::modelSummary summary(glmm.designs);
	//float lower = x_axis_limits.first;
	//float upper = x_axis_limits.second;
	//float diff = (upper - lower) / 19;
	//std::vector<int> int_vector;
	//std::vector<int> int_vector_series;
	//float fl_value;
	//float fl_value_series = 0;

	//for (int i = 0; i < 20; i++) {
	//	x_data[i] = lower + i * diff;
	//}

	//// save existing values	for series variables
	//if (n_series > 1) {
	//	switch (xaxis) {

	//	}
	//}

	//for (int j = 0; j < n_series; j++) {
	//	//update series variables if needed

	//	switch (xaxis) {
	//	case ClusterApp::XAxis::icc:
	//	{
	//		fl_value = glmm.statmodel.ixx_pars[0];
	//		for (int i = 0; i < 20; i++) {
	//			glmm.statmodel.ixx_pars[0] = x_data[i];
	//			glmm.update_parameters();
	//			extract_y(summary,i,j);
	//		}
	//		glmm.statmodel.ixx_pars[0] = fl_value;
	//	}


	//	}
	//}

	//replace series variables
}

//void ClusterApp::plotData::extract_y(ClusterApp::modelSummary& summary, int i, int series) {
//	switch (yaxis) {
//	case ClusterApp::YAxis::power:
//	{
//		glmm.power(summary);
//		switch (series) {
//		case 1:
//			y_data_1[i] = summary.power;
//			break;
//		case 2:
//			y_data_2[i] = summary.power;
//			break;
//		case 3:
//			y_data_3[i] = summary.power;
//			break;
//		}
//		break;
//	}		
//	case ClusterApp::YAxis::ci_width:
//	{
//		glmm.power(summary);
//		switch (series) {
//		case 1:
//			y_data_1[i] = summary.ci_width;
//			break;
//		case 2:
//			y_data_2[i] = summary.ci_width;
//			break;
//		case 3:
//			y_data_3[i] = summary.ci_width;
//			break;
//		}
//		break;
//	}
//		
//	}
//
//	// complete all the y axis options
//}