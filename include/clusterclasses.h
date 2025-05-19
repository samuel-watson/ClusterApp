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
#include "imgui_internal.h"
#include "modeltypes.h"
#include "glmmr.h"

#define DIV_ROUND_CLOSEST(n, d) ((((n) < 0) == ((d) < 0)) ? (((n) + (d)/2)/(d)) : (((n) - (d)/2)/(d)))

namespace ClusterApp {

    enum class PowerType {
        GLS = 0,
        BW = 1,
        Sat = 2,
        KR = 3,
        Sand = 4,
        SandBW = 5,
        DesignEffGLM = 6,
        DesignEffNP = 7,
        Box = 8
    };

    class sequencePeriod {
    public:
        bool    active = true;
        int     n = 10;
        bool    intervention = false;
        bool    intervention_2 = false;
        float   dose = 1.0;
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
        std::vector<std::vector<ClusterApp::sequencePeriod > >  periods;
        std::vector<int>                                        n_per_sequence;
        int                                                     default_clusters = 10;
        void initialise_data();
    public:
        int     sequences = 2;
        int     time = 1;
        int     crc_val = 0;
        design();
        void    add_sequence();
        void    add_sequence(int i);
        void    remove_sequence(int i);
        void    add_period();
        void    add_period(int t);
        void    remove_period(int t);
        bool*   active(int i, int t);
        bool*   intervention(int i, int t);
        bool*   intervention_2(int i, int t);
        int*    n(int i, int t);
        int*    n_clusters(int i);
        float*  dose(int i, int t);
        int     seq_by_cluster(int i);
        int     total_clusters();
        int     total_cluster_periods();
        int     total_n();
        double  mean_n();
        void    set_parallel(const int t, const int n, const int J);
        void    set_parallel_with_baseline(const int t, const int n, const int J);
        void    set_stepped_wedge(const int t, const int n, const int J, bool implement_period = false);
        void    set_stepped_wedge(const int t, const int n, const int J);
        void    set_crossover(const int n, const int J);
        void    set_staircase(const int t, const int n, const int J);
        void    set_factorial(const int t, const int n, const int J);
        void    split_sequences();
        void    combine_sequences();
        //std::vector<int>    periods_with_comparisons();
        ~design();
        bool    check(bool update_if_changed = true);
        void    apply_design(std::vector<std::vector<int> >& n);
        int     active_time_periods();
        void    swap_cells(int source, int target);
        void    copy_cells(int source, int target);
        void    move_cells(int source, int target);
        float   randomisation_ratio(int intervention_arm = 1, bool cluster=true);
    };

    class statisticalModel {
    public:
        Family                  family = Family::gaussian;
        Link                    link = Link::identity;
        Covariance              covariance = Covariance::exchangeable;
        IndividualCovariance    ind_covariance = IndividualCovariance::exchangeable;
        LinearPredictor         linearpredictor = LinearPredictor::time_fixed_effects;
        Sampling                sampling = Sampling::cross_sectional;
        ClusterApp::options&    option;
        int                     include_intercept = 1;
        float                   sigma = 1;
        std::array<float,3>     te_pars = { 0.5f,0.5f,0.5f };
        std::array<float,4>     ixx_pars = { 0.05f,0.8f,0.5f, 0.05f };
        std::array<float,7>     cov_pars = { 0.5f,0.5f,0.5f, 0.5f,0.5f,0.5f, 0.5f };
        std::vector<float>      beta_pars = std::vector<float>(1, 0.0f);
        std::vector<float>      c_vals = { 1.0f,1.0f,1.0f };
        int                     crc_val = 0;
        int                     crc_val_pars = 0;
        std::pair<bool, bool>   check();
        statisticalModel(ClusterApp::options& option_) : option(option_) {};
        void                    update_beta(ClusterApp::design& design);
        void                    set_beta_random(const double m, const double s);
        float                   alpha = 0.05;
        ClusterApp::PowerType   powertype = ClusterApp::PowerType::GLS;
        float                   target_power = 80.0;
        float                   cv_size = 0;
        float                   cv_size_within = 0;
        int                     mean_size = 10;
        //float                   quantile = 0.5;
    };

    struct modelSummary {
    public:
        double power = 0;
        double ci_width = 0;
        double se = 1;
        double dof = 1;
        double min_eff = 0;
        int total_n = 20;

        double power_2 = 0;
        double ci_width_2 = 0;
        double se_2 = 1;
        double dof_2 = 1;
        double min_eff_2 = 0;

        double power_12 = 0;
        double ci_width_12 = 0;
        double se_12 = 1;
        double dof_12 = 1;
        double min_eff_12 = 0;

        double design_effect = 0;
        double individual_se = 0;
        double individual_n = 1;
        double individual_var = 1;

        modelSummary(ClusterApp::design& design) : total_n(design.total_n()) {};
    };

    struct AppLog
    {
    public:
        ImGuiTextBuffer     Buf;
        ImGuiTextFilter     Filter;
        ImVector<int>       LineOffsets; // Index to lines offset. We maintain this with AddLog() calls.
        bool                AutoScroll = true;  // Keep scrolling if already at the bottom.
        bool                ShowMatrix = true;
        const char* cat[3] = { "info", "warn", "error" };
        AppLog(){};
        void Clear();
        void AddLog(const char* fmt, ...) IM_FMTARGS(2);
        void Draw(const char* title, bool* p_open = NULL);
    };

    class glmmModel {
    public:
        std::unique_ptr<glmm>           model;
        ClusterApp::statisticalModel&   statmodel;
        ClusterApp::options&            option;
        ClusterApp::design&             designs;
        ClusterApp::AppLog&             logger;
        std::string                     formula = "int+(1|gr(cl))";
        std::string                     family = "gaussian";
        std::string                     link = "identity";
        const std::vector<std::string>  colnames = { "cl","t","n","int","int2","int12","control" };
        boost::math::normal             norm = boost::math::normal(0.0, 1.0);
        double                          zcutoff = boost::math::quantile(norm, 0.975);
        int                             dof = 1;
        std::vector<double>             optimal_weights = { 0.5, 0.5 };
        glmmModel(ClusterApp::statisticalModel& statmodel_, ClusterApp::options& option_, ClusterApp::design& designs_, ClusterApp::AppLog& log_) : statmodel(statmodel_), option(option_), designs(designs_), logger(log_) {};
        ~glmmModel() = default;
        void                    update_formula();
        void                    update_parameters();
        void                    update_model_data(const Eigen::ArrayXXd& data);
        std::vector<double>     sim_data();
        void                    power(ClusterApp::modelSummary& summary, const ClusterApp::PowerType& powertype);
        void                    optimum(int N);
        float                   individual_n();
        double                  design_effect();
        double                  mean_individual_variance(bool weighted = true);
        std::pair<double,double> mean_outcome();
        std::vector<int>        round_weights(std::vector<float> w, int n);
    private:
        void                    power_(ClusterApp::modelSummary& summary, const ClusterApp::PowerType& powertype, const float var_inflate);
    };

    class modelUpdater {
    public:
        ClusterApp::design&             designs;
        ClusterApp::statisticalModel&   model;
        ClusterApp::modelSummary&       summary;
        ClusterApp::glmmModel&          glmm;
        ClusterApp::AppLog&             log;
        bool                            update = false;
        bool                            manual_n_optim = false;
        int                             de_mode = 0;
        Eigen::ArrayXXd                 data = Eigen::ArrayXXd::Constant(1, 7, 1); // order of columns cl, t, n, int1, int2, int1*int2
        std::vector<std::vector<double> >   optimum_data = { {0.5},{0.5} };
        std::vector<std::vector<int> >      optimum_n = { {10},{10} };
        bool                            requires_update = false;
        bool                            is_updating = false;
        bool                            plot_requires_update = false;
        bool                            initialized = false;

        modelUpdater(ClusterApp::design& designs_,
            ClusterApp::statisticalModel& model_,
            ClusterApp::modelSummary& summary_,
            ClusterApp::glmmModel& glmm_,
            ClusterApp::AppLog& log_);
        ~modelUpdater() = default;
        Eigen::ArrayXXd     generate_data();
        void                update_data();
        void                update_formula();
        void                update_parameters();
        void                update_summary_statistics();
        void                update_optimum();
        int                 sample_size_search(const bool& clusters, const ClusterApp::PowerType& powertype);
    };

    class plotData {
    public:
        ClusterApp::glmmModel&      glmm;
        ClusterApp::modelUpdater&   updater;
        ClusterApp::XAxis           xaxis = ClusterApp::XAxis::icc;
        ClusterApp::YAxis           yaxis = ClusterApp::YAxis::power;
        ClusterApp::XAxis           series = ClusterApp::XAxis::cac;
        int                         n_series = 1;
        bool                        multiple_series = false;
        int                         n_data_points = 20;
        std::pair<float, float>     x_axis_limits;
        std::pair<float, float>     y_axis_limits;
        float                       x_data[20];
        float                       y_data_1[20];
        float                       y_data_2[20];
        float                       y_data_3[20];
        float                       x_series[3];
        bool                        initialised = false;
        bool                        updating = false;
        int                         lower_int[2] = { 10,10 };
        int                         upper_int[2] = { 40, 100 };
        float                       upper_float[4] = { 0.25, 1.0, 1.0, 0.9 };
        float                       lower_float[4] = { 0.01, 0.0, 0.0, 0.1 };
        int                         crc_val = 0;
        plotData(ClusterApp::glmmModel& glmm_, ClusterApp::modelUpdater& updater_) : glmm(glmm_), updater(updater_) {};
        void    update_data();
        bool    check();
        void    extract_y(ClusterApp::modelSummary& summary, int i, int series);
        float   max_y();
        float   min_y();
    };

    class krigingData {
    public: 
        ClusterApp::glmmModel&      glmm;
        ClusterApp::modelUpdater&   updater;
        std::vector<int>            n_ind;
        std::vector<int>            n_cl;
        std::vector<float>          power;
        int                         lower_int[2] = { 10,10 };
        int                         upper_int[2] = { 40, 100 };
        float                       bandwidth = 1;
        float                       threshold_power = 0.8;
        bool                        initialised = false;
        bool                        surface_initialised = false;
        bool                        updating = false;
        bool                        start = false;
        float                       surface[400];
        float                       n_ind_grid[20];
        float                       n_cl_grid[20];
        char*                       n_ind_grid_label[20];
        char*                       n_cl_grid_label[20];
        float                       mu;
        int                         resample_total = 20;
        void    new_sample(int n = 25);
        void    generate_data();
        void    update(bool resample = true);
        krigingData(ClusterApp::glmmModel& glmm_, ClusterApp::modelUpdater& updater_) : glmm(glmm_), updater(updater_) { 
            generate_grid();
            new_sample();
        };
        void    generate_grid();
        void    set_power_type(PowerType type_);
    private:
        PowerType type = PowerType::GLS;
    };

    class modelChecker {
    public:
        ClusterApp::design&             designs;
        ClusterApp::statisticalModel&   model;
        ClusterApp::modelUpdater&       updater;
        ClusterApp::plotData&           plot;
        ClusterApp::options&            option;
        ClusterApp::krigingData&        krig;
        double                          update_interval = 1000;
        std::chrono::steady_clock       clock;
        std::chrono::time_point< std::chrono::steady_clock> t0;
        bool                            dcheck = false;
        std::pair<bool, bool>           mcheck = {false, false};
        bool                            pcheck = false;
        bool                            plotopen = false;
        modelChecker(ClusterApp::design& designs_,
            ClusterApp::statisticalModel& model_,
            ClusterApp::modelUpdater& updater_,
            ClusterApp::plotData& plot_,
            ClusterApp::krigingData& krig_,
            ClusterApp::options& option_,
            double interval = 1000);
        void check();
        void update();
        void check_time();
    };

    class appModel {
    public:
        ClusterApp::options&            option;
        ClusterApp::AppLog&             logger;
        ClusterApp::design              designs;
        ClusterApp::statisticalModel    model;
        ClusterApp::modelSummary        summary;
        ClusterApp::glmmModel           glmm;
        ClusterApp::modelUpdater        updater;
        ClusterApp::plotData            plot;
        ClusterApp::krigingData         krig;
        ClusterApp::modelChecker        checker;
        const int                       id;
        // specific values for menu bars
        int     family_item_current = 0;
        int     estimator_item_current = 0;
        int     link_item_current = 0;
        int     outcome_item_current = 0;
        bool    designopen = false;
        bool    modelopen = false;
        bool    optimalopen = false;
        bool    kriggeropen = false;
        bool    openkrig = false;
        bool    samplesizeopen = false;
        int     target_ind_size = 0;
        int     target_cl_size  = 0;
        bool    find_period_size_trigger = false;
        bool    find_clusters_trigger = false;

        appModel(ClusterApp::options& option_, ClusterApp::AppLog& logger_, const int id_);
    };
    
    /*
    class Model {
        public:
            ClusterApp::options            option;
            ClusterApp::AppLog             logger;            

            int size();
            ClusterApp::appModel& operator[](const int index);

            void add();
            void remove(const int index);
            bool active(const int index);

            Model();

        private:
            std::array<bool, 3>  isin{true,false,false};
            ClusterApp::appModel model1 = ClusterApp::appModel(option, logger, 0);
            ClusterApp::appModel model2 = ClusterApp::appModel(option, logger, 1);
            ClusterApp::appModel model3 = ClusterApp::appModel(option, logger, 2);
            //std::array<ClusterApp::appModel, 3> models_ = {(option, logger), (option, logger), (option, logger)};
    };*/

    // Prints a matrix to the log
    inline void AddMatrixToLog(const Eigen::MatrixXd& mat, ClusterApp::AppLog& log) {
        int rows = mat.rows();
        int cols = mat.cols();
        for (int i = 0; i < rows; i++) {
            std::string out = "[";
            for (int j = 0; j < cols; j++) {
                out += std::format("{:.3f}", mat(i, j));
                if (j < (cols - 1))out += " ";
            }
            out += "]";
            log.AddLog("[%05d] [%s] %s \n", ImGui::GetFrameCount(), log.cat[0], out.c_str());
        }
    }

    template<typename T>
    inline void AddVectorToLog(const std::vector<T>& vec, ClusterApp::AppLog& log, int max_print = 0) {
        int size = vec.size();
        int vals_to_print;
        if (max_print == 0) {
            vals_to_print = size;
        }
        else {
            vals_to_print = max_print > size ? size : max_print;
        }
        std::string out = "[";
        for (int j = 0; j < vals_to_print; j++) {
            out += std::format("{:.3f}", vec[j]);
            if (j < (size - 1))out += " ";
        }
        if (vals_to_print < size)out += " ...";
        out += "]";
        log.AddLog("[%05d] [%s] %s \n", ImGui::GetFrameCount(), log.cat[0], out.c_str());
    }



}

