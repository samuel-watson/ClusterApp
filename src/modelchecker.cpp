#include "clusterclasses.h"

ClusterApp::modelChecker::modelChecker(ClusterApp::design& designs_,
    ClusterApp::statisticalModel& model_,
    ClusterApp::modelUpdater& updater_,
    double interval) : designs(designs_), model(model_), updater(updater_), update_interval(interval), t0(clock.now()) {};

//timer_start([this]() {this->check(); }, interval);

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
        updater.update = false;
    }
}

//void ClusterApp::modelChecker::timer_start(std::function<void(void)> func, unsigned int interval) {
//    std::thread([func, interval]() {
//        while (true)
//        {
//            func();
//            std::this_thread::sleep_for(std::chrono::milliseconds(interval));
//        }
//        }).detach();
//};

void ClusterApp::modelChecker::check_time() {
    auto t1 = clock.now();
    std::chrono::duration<double, std::ratio<1, 1000> > diff = t1 - t0;
    if (diff.count() > update_interval) {
        check();
        t0 = clock.now();
    }
}
