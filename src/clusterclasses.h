#pragma once

#include <vector>
#include <map>
#include <memory>
#include <type_traits>
#include <chrono>
#include <thread>
#include <functional>
#include <numeric>
#include <boost/math/distributions/normal.hpp>
#include "modeltypes.h"
#include "glmmr/model.hpp"


namespace ClusterApp {

    class sequencePeriod {
    public:
        bool active = true;
        int n = 10;
        bool intervention = false;
        bool intervention_2 = false;
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
        int default_clusters = 1;
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
        int seq_by_cluster(int i);
        int total_clusters();
        int total_cluster_periods();
        int total_n();
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
    };

    struct modelSummary {
    public:
        double power = 0;
        double power_kr = 0;
        double power_bw = 0;
        double ci_width = 0;
        double ci_width_kr = 0;
        double ci_width_bw = 0;
        double se = 1;
        double se_kr = 1;
        double dof = 1;
        double dof_kr = 1;
        double dof_bw = 1;
        int total_n = 20;

        double power_2 = 0;
        double power_kr_2 = 0;
        double power_bw_2 = 0;
        double ci_width_2 = 0;
        double ci_width_kr_2 = 0;
        double ci_width_bw_2 = 0;
        double se_2 = 1;
        double se_kr_2 = 1;
        double dof_2 = 1;
        double dof_kr_2 = 1;
        double dof_bw_2 = 1;

        double power_12 = 0;
        double power_kr_12 = 0;
        double power_bw_12 = 0;
        double ci_width_12 = 0;
        double ci_width_kr_12 = 0;
        double ci_width_bw_12 = 0;
        double se_12 = 1;
        double se_kr_12 = 1;
        double dof_12 = 1;
        double dof_kr_12 = 1;
        double dof_bw_12 = 1;

        modelSummary(ClusterApp::design& design) : total_n(design.total_n()) {};
    };

    class glmmModel {
    public:
        std::unique_ptr<glmmr::ModelBits> modelbits;
        std::unique_ptr<glmmr::Model> model;
        ClusterApp::statisticalModel& statmodel;
        ClusterApp::options& option;
        std::string formula = "int+(1|gr(cl))";
        std::string family = "gaussian";
        std::string link = "identity";
        const std::vector<std::string> colnames = { "cl","t","n","int","int2","int12" };
        boost::math::normal norm = boost::math::normal(0.0, 1.0);
        double zcutoff = boost::math::quantile(norm, 0.975);
        int dof = 1;
        std::vector<double> optimal_weights = { 0.5, 0.5 };
        glmmModel(ClusterApp::statisticalModel& statmodel_, ClusterApp::options& option_) : statmodel(statmodel_), option(option_) {};
        ~glmmModel() = default;
        void update_formula();
        void update_parameters();
        void update_model_data(const Eigen::ArrayXXd& data);
        void power(ClusterApp::modelSummary& summary);
        void power_kr(ClusterApp::modelSummary& summary);
        void power_bw(ClusterApp::modelSummary& summary);
        void optimum(int N);
    };

    class modelUpdater {
    public:
        ClusterApp::design& designs;
        ClusterApp::statisticalModel& model;
        ClusterApp::modelSummary& summary;
        ClusterApp::glmmModel& glmm;
        bool update = false;
        Eigen::ArrayXXd data = Eigen::ArrayXXd::Constant(1, 6, 1); // order of columns cl, t, n, int1, int2, int1*int2
        std::vector<std::vector<double> > optimum_data = { {0.5},{0.5} };
        std::vector<std::vector<int> > optimum_n = { {10},{10} };
        modelUpdater(ClusterApp::design& designs_,
            ClusterApp::statisticalModel& model_,
            ClusterApp::modelSummary& summary_,
            ClusterApp::glmmModel& glmm_);
        ~modelUpdater() = default;
        void update_data();
        void update_formula();
        void update_parameters();
        void update_summary_statistics();
        void update_optimum();
    };

    class modelChecker {
    public:
        ClusterApp::design& designs;
        ClusterApp::statisticalModel& model;
        ClusterApp::modelUpdater& updater;
        unsigned int update_interval = 1000;
        modelChecker(ClusterApp::design& designs_,
            ClusterApp::statisticalModel& model_,
            ClusterApp::modelUpdater& updater_,
            unsigned int interval = 1000);
        void check();
        void timer_start(std::function<void(void)> func, unsigned int interval);
    };


}
