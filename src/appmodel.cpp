#include "clusterclasses.h"

// main holding classes

ClusterApp::appModel::appModel(ClusterApp::options& option_, ClusterApp::AppLog& logger_, const int id_) : option(option_), logger(logger_),
     designs(), model(option), summary(designs), glmm(model, option, designs, logger), updater(designs, model, summary, glmm, logger),
     plot(glmm, updater), krig(glmm, updater), checker(designs, model, updater, plot, krig, option), id(id_) {};

     /*
ClusterApp::Model::Model() : option(), logger() {};

int ClusterApp::Model::size(){
    int total = 0;
    for(int i = 0; i < 3; i++) if(isin[i])total++;
    return total;
}

ClusterApp::appModel& ClusterApp::Model::operator[](const int index){
    if(index == 0){
        return model1;
    } else if (index == 1){
        return model2; 
    } else if (index == 2){
        return model3;
    } else {
        if(index >= 3)throw std::runtime_error("Out of range");
    }
    // return (models_[index]);
}

void ClusterApp::Model::add(){
    for(int i = 0; i < 3; i++){
        if(!isin[i]){
            // models_[i] = std::make_shared<ClusterApp::appModel>(option, logger);
            isin[i] = true;
            break;
        }
    }
}

void ClusterApp::Model::remove(const int index){
    // models_[index].reset();
    if(size() > 1){
        isin[index] = false;
    }   
}

bool ClusterApp::Model::active(const int index){
    return isin[index];  
}
     */