#include "clusterclasses.h"

ClusterApp::sequencePeriod::sequencePeriod(const ClusterApp::sequencePeriod& sp) {
    active = sp.active;
    n = sp.n;
    intervention = sp.intervention;
    intervention_2 = sp.intervention_2;
}

void ClusterApp::sequencePeriod::set_active(bool active_) {
    active = active_;
};

void ClusterApp::sequencePeriod::set_n(int n_) {
    n = n_;
};

void ClusterApp::sequencePeriod::set_intervention(bool intervention_) {
    intervention = intervention_;
}

void ClusterApp::sequencePeriod::set_intervention_2(bool intervention_) {
    intervention_2 = intervention_;
}

bool ClusterApp::sequencePeriod::get_active() { return active; }

bool ClusterApp::sequencePeriod::get_intervention() { return intervention; }

bool ClusterApp::sequencePeriod::get_intervention_2() { return intervention_2; }

int ClusterApp::sequencePeriod::get_n() { return n; }

void ClusterApp::sequencePeriod::copy_properties(const ClusterApp::sequencePeriod& sp) {
    // use in case of drag & drop
    active = sp.active;
    n = sp.n;
    intervention = sp.intervention;
    intervention_2 = sp.intervention_2;
};

bool ClusterApp::sequencePeriod::operator< (const ClusterApp::sequencePeriod& sp) const {
    bool isless = false;
    if (sp.active) {
        if (!active) {
            isless = true;
        }
        else {
            if (sp.intervention) {
                if (!intervention) {
                    isless = true;
                }
                else {
                    if (sp.intervention_2) {
                        if (!intervention_2)isless = true;
                    }
                }
            }
        }
    }
    return isless;
}

bool ClusterApp::sequencePeriod::operator== (const ClusterApp::sequencePeriod& sp) const {
    return sp.active == active && sp.intervention == intervention && sp.intervention_2 == intervention_2;
}
