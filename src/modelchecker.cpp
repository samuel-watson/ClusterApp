#include "clusterclasses.h"

ClusterApp::modelChecker::modelChecker(ClusterApp::design& designs_,
    ClusterApp::statisticalModel& model_,
    ClusterApp::modelUpdater& updater_,
    unsigned int interval) : designs(designs_), model(model_), updater(updater_), update_interval(interval) {
    timer_start([this]() {this->check(); }, interval);
};

void ClusterApp::modelChecker::check() {
    bool dcheck = designs.check();
    std::pair<bool, bool> mcheck = model.check();
    if (dcheck || mcheck.first || mcheck.second) {
        updater.update = true;
        if (dcheck || mcheck.first) {
            updater.update_formula();
            updater.update_data();
        }
        if (mcheck.second && !mcheck.first)updater.update_parameters();
        std::this_thread::sleep_for(std::chrono::milliseconds(250));
        updater.update = false;
    }
}

void ClusterApp::modelChecker::timer_start(std::function<void(void)> func, unsigned int interval) {
    std::thread([func, interval]() {
        while (true)
        {
            func();
            std::this_thread::sleep_for(std::chrono::milliseconds(interval));
        }
        }).detach();
};
