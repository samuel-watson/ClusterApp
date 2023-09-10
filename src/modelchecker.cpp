#include "clusterclasses.h"

ClusterApp::modelChecker::modelChecker(ClusterApp::design& designs_,
    ClusterApp::statisticalModel& model_,
    ClusterApp::modelUpdater& updater_,
    ClusterApp::plotData& plot_,
    double interval) : designs(designs_), model(model_), updater(updater_), plot(plot_),
    update_interval(interval), t0(clock.now()) {};


void ClusterApp::modelChecker::check() {
    bool dcheck = designs.check();
    if (dcheck) model.update_beta(designs);
    std::pair<bool, bool> mcheck = model.check();    
    if (dcheck || mcheck.first || mcheck.second) {
        updater.update = true;
        if (dcheck || mcheck.first) {
            updater.update_formula();
            updater.update_data();
        }
        if (mcheck.second && !mcheck.first)updater.update_parameters();
        updater.update = false;
    }
    bool pcheck = plot.check();
    if (pcheck) plot.update_data();
}


void ClusterApp::modelChecker::check_time() {
    auto t1 = clock.now();
    std::chrono::duration<double, std::ratio<1, 1000> > diff = t1 - t0;
    if (diff.count() > update_interval) {
        check();
        t0 = clock.now();
    }
}
