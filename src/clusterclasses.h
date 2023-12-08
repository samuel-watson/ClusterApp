#pragma once

#include <vector>
#include <map>
#include <memory>
#include <type_traits>
#include <chrono>
#include <numeric>
#include <random>
#include <ctime>
#include <cstdint>
#include <format>
#include <boost/random.hpp>
#include "modeltypes.h"
#include "glmmr.h"

#define DIV_ROUND_CLOSEST(n, d) ((((n) < 0) == ((d) < 0)) ? (((n) + (d)/2)/(d)) : (((n) - (d)/2)/(d)))

namespace ClusterApp {

    enum class PowerType {
        GLS = 0,
        BW = 1,
        Sat = 2,
        KR = 3,
        DesignEff = 4,
        KR2 = 5
    };

    class sequencePeriod {
    public:
        bool active = true;
        int n = 10;
        bool intervention = false;
        bool intervention_2 = false;
        float dose = 1.0;
        sequencePeriod() {};
        sequencePeriod(bool active_) : active(active_) {};
        sequencePeriod(bool active_, int n_, bool intervention_) : active(active_), n(n_), intervention(intervention_) {};
        sequencePeriod(bool active_, int n_, bool intervention_, bool intervention2_) : active(active_), n(n_), intervention(intervention_), intervention_2(intervention2_) {};
        sequencePeriod(const ClusterApp::sequencePeriod& sp);
        void set_active(bool active_);
        void set_n(int n_);
        void set_intervention(bool intervention_);
        void set_intervention_2(bool intervention_);
        bool get_active();
        bool get_intervention();
        bool get_intervention_2();
        int get_n();
        void copy_properties(const ClusterApp::sequencePeriod& sp);
        bool operator< (const ClusterApp::sequencePeriod& sp) const;
        bool operator== (const ClusterApp::sequencePeriod& sp) const;
        virtual ~sequencePeriod() = default;
    };

    class design {
        std::vector<std::vector<ClusterApp::sequencePeriod > > periods;
        std::vector<int> n_per_sequence;
        int default_clusters = 10;
        void initialise_data();
    public:
        int sequences = 2;
        int time = 1;
        int crc_val = 0;
        design();
        void add_sequence();
        void add_sequence(int i);
        void remove_sequence(int i);
        void add_period();
        void add_period(int t);
        void remove_period(int t);
        bool* active(int i, int t);
        bool* intervention(int i, int t);
        bool* intervention_2(int i, int t);
        int* n(int i, int t);
        int* n_clusters(int i);
        float* dose(int i, int t);
        int seq_by_cluster(int i);
        int total_clusters();
        int total_cluster_periods();
        int total_n();
        double mean_n();
        void set_parallel(const int t, const int n, const int J);
        void set_parallel_with_baseline(const int t, const int n, const int J);
        void set_stepped_wedge(const int t, const int n, const int J);
        void set_crossover(const int n, const int J);
        void set_staircase(const int t, const int n, const int J);
        void set_factorial(const int t, const int n, const int J);
        void split_sequences();
        void combine_sequences();
        ~design();
        bool check(bool update_if_changed = true);
        void apply_design(std::vector<std::vector<int> >& n);
        int active_time_periods();
        void swap_cells(int source, int target);
        void copy_cells(int source, int target);
        void move_cells(int source, int target);
        float randomisation_ratio(int intervention_arm = 1, bool cluster=true);
    };

    class statisticalModel {
    public:
        Family family = Family::gaussian;
        Link link = Link::identity;
        Covariance covariance = Covariance::exchangeable;
        IndividualCovariance ind_covariance = IndividualCovariance::exchangeable;
        LinearPredictor linearpredictor = LinearPredictor::time_fixed_effects;
        Sampling sampling = Sampling::cross_sectional;
        int include_intercept = 1;
        float sigma = 1;
        std::vector<float> te_pars = { 0.5f,0.5f,0.5f };
        std::vector<float> ixx_pars = { 0.05f,0.8f,0.5f };
        std::vector<float> cov_pars = std::vector<float>(5, 0.5f);
        std::vector<float> beta_pars = std::vector<float>(1, 0.0f);
        std::vector<float> c_vals = { 1.0f,1.0f,1.0f };
        statisticalModel() {};
        int crc_val = 0;
        int crc_val_pars = 0;
        std::pair<bool, bool> check();
        void update_beta(ClusterApp::design& design);
        void set_beta_random(const double m, const double s);
        float alpha = 0.05;
    };

    struct modelSummary {
    public:
        double power = 0;
        double power_kr = 0;
        double power_kr2 = 0;
        double power_bw = 0;
        double power_box = 0;
        double power_sat = 0;
        double ci_width = 0;
        double ci_width_kr = 0;
        double ci_width_kr2 = 0;
        double ci_width_bw = 0;
        double ci_width_sat = 0;
        double se = 1;
        double se_kr = 1;
        double se_kr2 = 1;
        double dof = 1;
        double dof_kr = 1;
        double dof_bw = 1;
        double dof_box = 1;
        int total_n = 20;

        double power_2 = 0;
        double power_kr_2 = 0;
        double power_kr2_2 = 0;
        double power_bw_2 = 0;
        double power_box_2 = 0;
        double power_sat_2 = 0;
        double ci_width_2 = 0;
        double ci_width_kr_2 = 0;
        double ci_width_kr2_2 = 0;
        double ci_width_bw_2 = 0;
        double ci_width_sat_2 = 0;
        double se_2 = 1;
        double se_kr_2 = 1;
        double se_kr2_2 = 1;
        double dof_2 = 1;
        double dof_kr_2 = 1;
        double dof_bw_2 = 1;
        double dof_box_2 = 1;

        double power_12 = 0;
        double power_kr_12 = 0;
        double power_kr2_12 = 0;
        double power_bw_12 = 0;
        double power_box_12 = 0;
        double power_sat_12 = 0;
        double ci_width_12 = 0;
        double ci_width_kr_12 = 0;
        double ci_width_kr2_12 = 0;
        double ci_width_bw_12 = 0;
        double ci_width_sat_12 = 0;
        double se_12 = 1;
        double se_kr_12 = 1;
        double se_kr2_12 = 1;
        double dof_12 = 1;
        double dof_kr_12 = 1;
        double dof_bw_12 = 1;
        double dof_box_12 = 1;

        double design_effect = 0;
        double individual_se = 0;
        double power_de = 0;
        double ci_width_de = 0;
        double se_de = 1;
        double individual_n = 1;
        double individual_var = 1;

        modelSummary(ClusterApp::design& design) : total_n(design.total_n()) {};
    };

    class glmmModel {
    public:
        std::unique_ptr<glmm> model;
        ClusterApp::statisticalModel& statmodel;
        ClusterApp::options& option;
        ClusterApp::design& designs;
        std::string formula = "int+(1|gr(cl))";
        std::string family = "gaussian";
        std::string link = "identity";
        const std::vector<std::string> colnames = { "cl","t","n","int","int2","int12" };
        boost::math::normal norm = boost::math::normal(0.0, 1.0);
        double zcutoff = boost::math::quantile(norm, 0.975);
        int dof = 1;
        std::vector<double> optimal_weights = { 0.5, 0.5 };
        glmmModel(ClusterApp::statisticalModel& statmodel_, ClusterApp::options& option_, ClusterApp::design& designs_) : statmodel(statmodel_), option(option_), designs(designs_) {};
        ~glmmModel() = default;
        void update_formula();
        void update_parameters();
        void update_model_data(const Eigen::ArrayXXd& data);
        std::vector<double> sim_data();
        void power(ClusterApp::modelSummary& summary);
        void power_kr(ClusterApp::modelSummary& summary);
        void power_box(ClusterApp::modelSummary& summary);
        void power_bw(ClusterApp::modelSummary& summary);
        void optimum(int N);
        float individual_n();
        double design_effect();
        void power_de(ClusterApp::modelSummary& summary, int type);
        double mean_individual_variance(bool weighted = true);
        std::pair<double,double> mean_outcome();
        std::vector<int> round_weights(std::vector<float> w, int n);
    };

    class modelUpdater {
    public:
        ClusterApp::design& designs;
        ClusterApp::statisticalModel& model;
        ClusterApp::modelSummary& summary;
        ClusterApp::glmmModel& glmm;
        bool update = false;
        bool manual_n_optim = false;
        int de_mode = 0;
        Eigen::ArrayXXd data = Eigen::ArrayXXd::Constant(1, 6, 1); // order of columns cl, t, n, int1, int2, int1*int2
        std::vector<std::vector<double> > optimum_data = { {0.5},{0.5} };
        std::vector<std::vector<int> > optimum_n = { {10},{10} };
        modelUpdater(ClusterApp::design& designs_,
            ClusterApp::statisticalModel& model_,
            ClusterApp::modelSummary& summary_,
            ClusterApp::glmmModel& glmm_);
        ~modelUpdater() = default;
        Eigen::ArrayXXd generate_data();
        void update_data();
        void update_formula();
        void update_parameters();
        void update_summary_statistics();
        void update_optimum();
    };

    class plotData {
    public:
        ClusterApp::glmmModel& glmm;
        ClusterApp::modelUpdater& updater;
        ClusterApp::XAxis xaxis = ClusterApp::XAxis::icc;
        ClusterApp::YAxis yaxis = ClusterApp::YAxis::power;
        ClusterApp::XAxis series = ClusterApp::XAxis::cac;
        int n_series = 1;
        bool multiple_series = false;
        int n_data_points = 20;
        plotData(ClusterApp::glmmModel& glmm_, ClusterApp::modelUpdater& updater_) : glmm(glmm_), updater(updater_) {};
        std::pair<float, float> x_axis_limits;
        std::pair<float, float> y_axis_limits;
        float x_data[20];
        float y_data_1[20];
        float y_data_2[20];
        float y_data_3[20];
        float x_series[3];
        void update_data();
        bool initialised = false;
        bool check();
        int crc_val = 0;
        void extract_y(ClusterApp::modelSummary& summary, int i, int series);
        float max_y();
        float min_y();
        bool updating = false;
        int lower_int[2] = { 10,10 };
        int upper_int[2] = { 40, 100 };
        float upper_float[4] = { 0.25, 1.0, 1.0, 0.9 };
        float lower_float[4] = { 0.01, 0.0, 0.0, 0.1 };
    };

    class krigingData {
    public: 
        ClusterApp::glmmModel& glmm;
        ClusterApp::modelUpdater& updater;
        std::vector<int> n_ind;
        std::vector<int> n_cl;
        std::vector<float> power;
        void new_sample(int n = 25);
        void generate_data();
        void update(bool resample = true);
        int lower_int[2] = { 10,10 };
        int upper_int[2] = { 40, 100 };
        float bandwidth = 1;
        float threshold_power = 0.8;
        krigingData(ClusterApp::glmmModel& glmm_, ClusterApp::modelUpdater& updater_) : glmm(glmm_), updater(updater_) { 
            generate_grid();
            new_sample();
        };
        bool initialised = false;
        bool surface_initialised = false;
        bool updating = false;
        bool start = false;
        float surface[400];
        float n_ind_grid[20];
        float n_cl_grid[20];
        char* n_ind_grid_label[20];
        char* n_cl_grid_label[20];
        float mu;
        void generate_grid();
        int resample_total = 20;
        void set_power_type(PowerType type_);
    private:
        PowerType type = PowerType::GLS;
    };

    class modelChecker {
    public:
        ClusterApp::design& designs;
        ClusterApp::statisticalModel& model;
        ClusterApp::modelUpdater& updater;
        ClusterApp::plotData& plot;
        ClusterApp::options& option;
        ClusterApp::krigingData& krig;
        double update_interval = 1000;
        std::chrono::steady_clock clock;
        std::chrono::time_point< std::chrono::steady_clock> t0;
        modelChecker(ClusterApp::design& designs_,
            ClusterApp::statisticalModel& model_,
            ClusterApp::modelUpdater& updater_,
            ClusterApp::plotData& plot_,
            ClusterApp::krigingData& krig_,
            ClusterApp::options& option_,
            double interval = 1000);
        void check();
        void check_time();
    };


}

