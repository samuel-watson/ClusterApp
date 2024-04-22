#include "clusterclasses.h"

ClusterApp::modelChecker::modelChecker(ClusterApp::design& designs_,
    ClusterApp::statisticalModel& model_,
    ClusterApp::modelUpdater& updater_,
    ClusterApp::plotData& plot_,
    ClusterApp::krigingData& krig_,
    ClusterApp::options& option_,
    double interval) : designs(designs_), model(model_), updater(updater_), plot(plot_), krig(krig_), option(option_),
        update_interval(interval), t0(clock.now()) {};


void ClusterApp::modelChecker::check() {
    // if (!plot.updating && !krig.updating) {        
    // }  
    if(!dcheck) dcheck = designs.check();
    static std::pair<bool, bool> mcheck_new;
    if(!(mcheck.first && mcheck.second)){
        mcheck_new = model.check();
    }
    if(!mcheck.first) mcheck.first = mcheck_new.first;
    if(!mcheck.second) mcheck.second = mcheck_new.second;
    if(!pcheck) pcheck = plot.check();
    if(!updater.requires_update && (dcheck || mcheck.first || mcheck.second))
    {
        updater.requires_update = true;
        if(dcheck) updater.log.AddLog("[%05d] [%s] %s \n", ImGui::GetFrameCount(), updater.log.cat[0], "Design change detected.");
        if(mcheck.first || mcheck.second) updater.log.AddLog("[%05d] [%s] %s \n", ImGui::GetFrameCount(), updater.log.cat[0], "Statistical model change detected.");
    }  
    if(!updater.plot_requires_update && pcheck)
    {
        updater.plot_requires_update = true;
        if(pcheck) updater.log.AddLog("[%05d] [%s] %s \n", ImGui::GetFrameCount(), updater.log.cat[0], "Plot change detected.");
    }   
}

void ClusterApp::modelChecker::update() {
    updater.log.AddLog("[%05d] [%s] %s [%s, %s, %s] \n", ImGui::GetFrameCount(), updater.log.cat[0], "Running updater: ",dcheck ? "DESIGN" : ".", mcheck.first ? "MODEL 1": ".", mcheck.second ? "MODEL 2": ".");
    // if (!plot.updating && !krig.updating) {                    
    // } 
    updater.requires_update = false;  
    if (dcheck) model.update_beta(designs);
    if (dcheck || mcheck.first || mcheck.second) {
        updater.log.AddLog("[%05d] [%s] %s \n", ImGui::GetFrameCount(), updater.log.cat[0], "Updating model.");
        updater.update = true;
        if (dcheck || mcheck.first) {
            updater.update_formula();
            updater.update_data();
        }
        if (mcheck.second && !mcheck.first) updater.update_parameters();
        updater.update = false;
        updater.log.AddLog("[%05d] [%s] GLMM checksum: %d \n", ImGui::GetFrameCount(), updater.log.cat[0], designs.crc_val);
        updater.log.AddLog("[%05d] [%s] Model checksum: %d \n", ImGui::GetFrameCount(), updater.log.cat[0], model.crc_val);
        updater.log.AddLog("[%05d] [%s] Parameter checksum: %d \n", ImGui::GetFrameCount(), updater.log.cat[0], model.crc_val_pars);
    }
    if ((plot.initialised && (pcheck || dcheck || mcheck.first || mcheck.second) && option.plotter) || option.plotter && !plot.initialised){
        plot.update_data();
        bool update_crc_plot = plot.check();
        updater.plot_requires_update = false; 
        pcheck = false;
    } 
    dcheck = false;
    mcheck.first = false;
    mcheck.second = false;
    bool update_crc = designs.check();
    std::pair<bool, bool> mcheck_new = model.check();
}


void ClusterApp::modelChecker::check_time() {
    auto t1 = clock.now();
    std::chrono::duration<double, std::ratio<1, 1000> > diff = t1 - t0;
    if (diff.count() > update_interval) {
        check();
        if(option.auto_update || !updater.initialized){
            update();
            updater.initialized = true;
        } 
        t0 = clock.now();
    }
}
