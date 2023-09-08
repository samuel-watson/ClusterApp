#include "clusterclasses.h"

void ClusterApp::design::initialise_data() {

    for (int i = 0; i < sequences; i++) {
        periods.emplace_back();
        for (int t = 0; t < time; t++) {
            periods[i].push_back(ClusterApp::sequencePeriod());
        }
    }
    periods[1][0].set_intervention(true);
    n_per_sequence = { default_clusters,default_clusters };
}

ClusterApp::design::design() {
    initialise_data();
};

void ClusterApp::design::add_sequence() {
    periods.emplace_back();
    int len = periods.size() - 1;
    for (int t = 0; t < time; t++) {
        periods[len].push_back(ClusterApp::sequencePeriod(true));
    }
    n_per_sequence.push_back(default_clusters);
    sequences++;
};

void ClusterApp::design::add_sequence(int i) {
    std::vector<ClusterApp::sequencePeriod > sqvec;
    for (int t = 0; t < time; t++) {
        sqvec.push_back(ClusterApp::sequencePeriod(true));
    }
    periods.insert(periods.begin() + i, sqvec);
    n_per_sequence.insert(n_per_sequence.begin() + i, default_clusters);
    sequences++;
};

void ClusterApp::design::remove_sequence(int i) {
    sequences--;
    periods.erase(periods.begin() + i);
    n_per_sequence.erase(n_per_sequence.begin() + i);
}

void ClusterApp::design::add_period() {
    for (int i = 0; i < sequences; i++) {
        periods[i].push_back(ClusterApp::sequencePeriod(true));
    }
    time++;
}

void ClusterApp::design::add_period(int t) {
    for (int i = 0; i < sequences; i++) {
        periods[i].insert(periods[i].begin() + t, ClusterApp::sequencePeriod(true));
    }
    time++;
}

void ClusterApp::design::remove_period(int t) {
    time--;
    for (int i = 0; i < sequences; i++) {
        periods[i].erase(periods[i].begin() + t);
    }
}

bool* ClusterApp::design::active(int i, int t) {
    return &periods[i][t].active;
}

bool* ClusterApp::design::intervention(int i, int t) {
    return &periods[i][t].intervention;
}

bool* ClusterApp::design::intervention_2(int i, int t) {
    return &periods[i][t].intervention_2;
}

int* ClusterApp::design::n(int i, int t) {
    return &periods[i][t].n;
}

int* ClusterApp::design::n_clusters(int i) {
    return &n_per_sequence[i];
}

int ClusterApp::design::seq_by_cluster(int i) {
    int count = 0;
    int seq = 0;
    bool found = false;
    while (!found && seq < sequences) {
        found = i >= count && i < count + n_per_sequence[seq];
        if (!found) {
            count += n_per_sequence[seq];
            seq++;            
        }
    }
    return seq;
};

int ClusterApp::design::total_clusters() {
    int total = 0;
    for (int i = 0; i < sequences; i++) {
        bool any_active = false;
        for (int t = 0; t < time; t++) {
            if (*active(i, t)) {
                any_active = true;
                break;
            }
        }
        if(any_active)total += n_per_sequence[i];
    }
    return total;
}

int ClusterApp::design::total_cluster_periods() {
    int total = 0;
    for (int i = 0; i < sequences; i++) {
        for (int t = 0; t < time; t++) {
            if (*active(i, t)) {
                total += n_per_sequence[i];
            }
        }
    }
    return total;
}

int ClusterApp::design::total_n() {
    int total = 0;
    for (int i = 0; i < sequences; i++) {
        for (int t = 0; t < time; t++) {
            if (*active(i, t)) {
                total += n_per_sequence[i] * (*n(i,t));
            }
        }
    }
    return total;
}

double ClusterApp::design::mean_n() {
    int totaln = total_n();
    int totalcl_periods = total_cluster_periods();
    return totaln / (double)totalcl_periods;
}

void ClusterApp::design::set_parallel(const int t, const int n, const int J) {
    sequences = 1;
    time = 1;
    periods.clear();
    n_per_sequence.clear();
    periods.emplace_back();
    for (int i = 0; i < t; i++) {
        periods[0].push_back(ClusterApp::sequencePeriod(true, n, false));
    }
    n_per_sequence.push_back(J);
    periods.emplace_back();
    for (int i = 0; i < t; i++) {
        periods[1].push_back(ClusterApp::sequencePeriod(true, n, true));
    }
    n_per_sequence.push_back(J);
    sequences = 2;
    time = t;
}

void ClusterApp::design::set_parallel_with_baseline(const int t, const int n, const int J) {
    sequences = 1;
    time = 1;
    periods.clear();
    n_per_sequence.clear();
    periods.emplace_back();
    for (int i = 0; i < t + 1; i++) {
        periods[0].push_back(ClusterApp::sequencePeriod(true, n, false));
    }
    n_per_sequence.push_back(J);
    periods.emplace_back();
    for (int i = 0; i < t + 1; i++) {
        periods[1].push_back(ClusterApp::sequencePeriod(true, n, true));
    }
    periods[1][0].set_intervention(false);
    n_per_sequence.push_back(J);
    sequences = 2;
    time = t + 1;
}

void ClusterApp::design::set_stepped_wedge(const int t, const int n, const int J) {
    sequences = 1;
    time = 1;
    periods.clear();
    n_per_sequence.clear();
    for (int s = 0; s < (t - 1); s++) {
        periods.emplace_back();
        for (int u = 0; u < t; u++) {
            bool has_intervention = u > s;
            periods[s].push_back(ClusterApp::sequencePeriod(true, n, has_intervention));
        }
        n_per_sequence.push_back(J);
    }
    sequences = t - 1;
    time = t;
}

void ClusterApp::design::set_crossover(const int n, const int J) {
    sequences = 1;
    time = 1;
    periods.clear();
    n_per_sequence.clear();
    periods.emplace_back();
    periods[0].push_back(ClusterApp::sequencePeriod(true, n, false));
    periods[0].push_back(ClusterApp::sequencePeriod(true, n, true));
    n_per_sequence.push_back(J);
    periods.emplace_back();
    periods[1].push_back(ClusterApp::sequencePeriod(true, n, true));
    periods[1].push_back(ClusterApp::sequencePeriod(true, n, false));
    n_per_sequence.push_back(J);
    sequences = 2;
    time = 2;
}

void ClusterApp::design::set_staircase(const int t, const int n, const int J) {
    sequences = 1;
    time = 1;
    periods.clear();
    n_per_sequence.clear();
    for (int s = 0; s < (t + 1); s++) {
        periods.emplace_back();
        for (int u = 0; u < t; u++) {
            bool is_active = (u == s || u == s - 1);
            bool has_intervention = u == s;
            periods[s].push_back(ClusterApp::sequencePeriod(is_active, n, has_intervention));
        }
        n_per_sequence.push_back(J);
    }
    sequences = t + 1;
    time = t;
}

void ClusterApp::design::set_factorial(const int t, const int n, const int J) {
    sequences = 1;
    time = 1;
    periods.resize(4);
    n_per_sequence.resize(4);
    for (int i = 0; i < 4; i++)n_per_sequence[i] = J;
    periods[0].resize(t);
    periods[1].resize(t);
    periods[2].resize(t);
    periods[3].resize(t);

    for (int i = 0; i < t; i++) {
        periods[0][i] = ClusterApp::sequencePeriod(true, n, false, false);
        periods[1][i] = ClusterApp::sequencePeriod(true, n, true, false);
        periods[2][i] = ClusterApp::sequencePeriod(true, n, false, true);
        periods[3][i] = ClusterApp::sequencePeriod(true, n, true, true);
    }
    sequences = 4;
    time = t;
}

void ClusterApp::design::split_sequences() {
    int idx = 0;
    int curr_sequences = sequences;
    for (int i = 0; i < curr_sequences; i++) {
        if (n_per_sequence[idx] > 1) {
            int new_n = n_per_sequence[idx];  //2, idx = 0
            n_per_sequence[idx] = 1; //n[0] = 1;
            for (int j = 0; j < (new_n - 1); j++) {
                add_sequence(idx);
                for (int t = 0; t < time; t++) {
                    *active(idx, t) = *active(idx + 1, t);
                    *intervention(idx, t) = *intervention(idx + 1, t);
                    *intervention_2(idx, t) = *intervention_2(idx + 1, t);
                }
                n_per_sequence[idx] = 1;
            }
            idx += new_n;
        }
        else {
            idx++;
        }
    }
}

void ClusterApp::design::combine_sequences() {
    std::map<std::vector<ClusterApp::sequencePeriod >, int> vecs;
    int idx = 0;
    int curr_sequences = sequences;
    for (int i = 0; i < curr_sequences; i++) {
        if (vecs.find(periods[idx]) != vecs.end()) {
            int newint = 1;
            vecs[periods[idx]]++;
            remove_sequence(idx);
        }
        else {
            vecs[periods[idx]] = 1;
            idx++;
        }
    }
    for (int i = 0; i < sequences; i++) {
        n_per_sequence[i] = vecs[periods[i]];
    }
}

ClusterApp::design::~design() {
    for (int i = 0; i < sequences; i++) {
        periods[i].clear();
    }
};

bool ClusterApp::design::check(bool update_if_changed) {
    CRC crc;
    for (int i = 0; i < sequences; i++) {
        for (int t = 0; t < time; t++) {
            crc(*active(i, t));
            crc(*intervention(i, t));
            crc(*intervention_2(i, t));
            crc(*n(i, t));
        }
        crc(n_per_sequence[i]);
    }
    bool has_changed = crc_val != crc.get();
    if (has_changed && update_if_changed)crc_val = crc.get();
    return has_changed;
}

void ClusterApp::design::apply_design(std::vector<std::vector<int> >& newN) {
    for (int i = 0; i < newN.size(); i++) {
        int seq = seq_by_cluster(i);
        for (int t = 0; t < time; t++) {
            if (*active(seq, t)) {
                if (newN[i][t] > 0) {
                    *n(seq, t) = (int)newN[i][t];
                }
                else {
                    *active(seq, t) = false;
                    *n(seq, t) = 0;
                }
            }
        }
    }
}

int ClusterApp::design::active_time_periods() {
    std::vector<int> active_period(time);
    std::fill(active_period.begin(), active_period.end(), 0);
    for (int i = 0; i < sequences; i++) {
        for (int t = 0; t < time; t++) {
            if (*active(i, t))active_period[t] = 1;
        }
    }
    int totalT = std::reduce(active_period.begin(),active_period.end());
    return totalT;
}

void ClusterApp::design::swap_cells(int source, int target) {
    int i_source = source / time;
    int t_source = source - i_source * time;
    int i_target = target / time;
    int t_target = target - i_target * time;

    ClusterApp::sequencePeriod sp_target(periods[i_target][t_target]);
    periods[i_target][t_target].copy_properties(periods[i_source][t_source]);
    periods[i_source][t_source].copy_properties(sp_target);
}

void ClusterApp::design::copy_cells(int source, int target) {
    int i_source = source / time;
    int t_source = source - i_source * time;
    int i_target = target / time;
    int t_target = target - i_target * time;
    periods[i_target][t_target].copy_properties(periods[i_source][t_source]);
}

void ClusterApp::design::move_cells(int source, int target) {
    int i_source = source / time;
    int t_source = source - i_source * time;
    int i_target = target / time;
    int t_target = target - i_target * time;
    periods[i_target][t_target].copy_properties(periods[i_source][t_source]);
    *active(i_source, t_source) = false;
}