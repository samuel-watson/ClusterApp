#pragma once

#include "general.h"
#include "modelbits.hpp"
#include "matrixw.hpp"
#include "randomeffects.hpp"
#include "openmpheader.h"
#include "maths.h"
#include "matrixfield.h"

namespace glmmr {

    using namespace Eigen;

    template<typename modeltype>
    class ModelMatrix {
    public:
        modeltype& model;
        glmmr::MatrixW<modeltype> W;
        glmmr::RandomEffects<modeltype>& re;
        ModelMatrix(modeltype& model_, glmmr::RandomEffects<modeltype>& re_) : model(model_), W(model_), re(re_) { gen_sigma_blocks(); };
        MatrixXd information_matrix();
        MatrixXd Sigma(bool inverse = false);
        MatrixXd observed_information_matrix();
        MatrixXd sandwich_matrix();
        std::vector<MatrixXd> sigma_derivatives();
        MatrixXd information_matrix_theta();
        kenward_data kenward_roger();
        MatrixXd linpred();
        vector_matrix b_score();
        vector_matrix re_score();
        matrix_matrix hess_and_grad();
        VectorXd log_gradient(const VectorXd& v, bool beta = false);
        std::vector<glmmr::SigmaBlock> get_sigma_blocks();

    private:
        std::vector<glmmr::SigmaBlock> sigma_blocks;
        void gen_sigma_blocks();
        MatrixXd sigma_block(int b, bool inverse = false);
        MatrixXd sigma_builder(int b, bool inverse = false);
        MatrixXd information_matrix_by_block(int b);
    };

}

template<typename modeltype>
inline std::vector<glmmr::SigmaBlock> glmmr::ModelMatrix<modeltype>::get_sigma_blocks() {
    return sigma_blocks;
}

template<typename modeltype>
inline MatrixXd glmmr::ModelMatrix<modeltype>::information_matrix_by_block(int b) {
    ArrayXi rows = Map<ArrayXi, Unaligned>(sigma_blocks[b].RowIndexes.data(), sigma_blocks[b].RowIndexes.size());
    MatrixXd X = glmmr::Eigen_ext::submat(model.linear_predictor.X(), rows, ArrayXi::LinSpaced(model.linear_predictor.P(), 0, model.linear_predictor.P() - 1));
    MatrixXd S = sigma_block(b, true);
    MatrixXd M = X.transpose() * S * X;
    return M;
}

template<typename modeltype>
inline MatrixXd glmmr::ModelMatrix<modeltype>::information_matrix() {
    W.update();
    MatrixXd M = MatrixXd::Zero(model.linear_predictor.P(), model.linear_predictor.P());
    for (unsigned int i = 0; i < sigma_blocks.size(); i++) {
        M += information_matrix_by_block(i);
    }
    return M;
}

template<typename modeltype>
inline void glmmr::ModelMatrix<modeltype>::gen_sigma_blocks() {
    int block_counter = 0;
    intvec2d block_ids(model.n());
    int block_size;
    sparse Z = model.covariance.Z_sparse();
    int i, j, k;
    auto it_begin = Z.Ai.begin();
    for (int b = 0; b < model.covariance.B(); b++) {
        block_size = model.covariance.block_dim(b);
        for (i = 0; i < block_size; i++) {
#pragma omp parallel for shared(it_begin, i)
            for (j = 0; j < model.n(); j++) {
                auto it = std::find(it_begin + Z.Ap[j], it_begin + Z.Ap[j + 1], (i + block_counter));
                if (it != (it_begin + Z.Ap[j + 1])) {
                    block_ids[j].push_back(b);
                }
            }
        }
        block_counter += block_size;
    }
    block_counter = 0;
    intvec idx_matches;
    int n_matches;
    for (i = 0; i < model.n(); i++) {
        if (block_counter == 0) {
            glmmr::SigmaBlock newblock(block_ids[i]);
            newblock.add_row(0);
            sigma_blocks.push_back(newblock);
        }
        else {
            for (j = 0; j < block_counter; j++) {
                if (sigma_blocks[j] == block_ids[i]) {
                    idx_matches.push_back(j);
                }
            }
            n_matches = idx_matches.size();
            if (n_matches == 0) {
                glmmr::SigmaBlock newblock(block_ids[i]);
                newblock.add_row(i);
                sigma_blocks.push_back(newblock);
            }
            else if (n_matches == 1) {
                sigma_blocks[idx_matches[0]].add(block_ids[i]);
                sigma_blocks[idx_matches[0]].add_row(i);
            }
            else if (n_matches > 1) {
                std::reverse(idx_matches.begin(), idx_matches.end());
                for (k = 0; k < (n_matches - 1); k++) {
                    sigma_blocks[idx_matches[n_matches - 1]].merge(sigma_blocks[idx_matches[k]]);
                    sigma_blocks.erase(sigma_blocks.begin() + idx_matches[k]);
                }
            }
        }
        idx_matches.clear();
        block_counter = sigma_blocks.size();
    }
}

template<typename modeltype>
inline MatrixXd glmmr::ModelMatrix<modeltype>::Sigma(bool inverse) {
    W.update();
    MatrixXd S = sigma_builder(0, inverse);
    return S;
}

template<typename modeltype>
inline MatrixXd glmmr::ModelMatrix<modeltype>::sigma_block(int b,
    bool inverse) {
    sparse ZLs = submat_sparse(model.covariance.ZL_sparse(), sigma_blocks[b].RowIndexes);
    MatrixXd ZL = sparse_to_dense(ZLs, false);
    MatrixXd S = ZL * ZL.transpose();
    for (int i = 0; i < S.rows(); i++) {
        S(i, i) += 1 / W.W()(sigma_blocks[b].RowIndexes[i]);
    }
    if (inverse) {
        S = S.llt().solve(MatrixXd::Identity(S.rows(), S.cols()));
    }
    return S;
}

template<typename modeltype>
inline MatrixXd glmmr::ModelMatrix<modeltype>::sigma_builder(int b,
    bool inverse) {
    int B_ = sigma_blocks.size();
    if (b == B_ - 1) {
        return sigma_block(b, inverse);
    }
    else {
        MatrixXd mat1 = sigma_block(b, inverse);
        MatrixXd mat2;
        if (b == B_ - 2) {
            mat2 = sigma_block(b + 1, inverse);
        }
        else {
            mat2 = sigma_builder(b + 1, inverse);
        }
        int n1 = mat1.rows();
        int n2 = mat2.rows();
        MatrixXd dmat = MatrixXd::Zero(n1 + n2, n1 + n2);
        dmat.block(0, 0, n1, n1) = mat1;
        dmat.block(n1, n1, n2, n2) = mat2;
        return dmat;
    }
}

template<typename modeltype>
inline MatrixXd glmmr::ModelMatrix<modeltype>::observed_information_matrix() {
    // this works but its too slow doing all the cross partial derivatives
    //MatrixXd XZ(n_,P_+Q_);
    //int iter = zu_.cols();
    //XZ.leftCols(P_) = linp_.X();
    //XZ.rightCols(Q_) = sparse_to_dense(ZL_,false);
    //MatrixXd result = MatrixXd::Zero(P_+Q_,P_+Q_);
    //MatrixXd I = MatrixXd::Identity(P_+Q_,P_+Q_);
    //dblvec params(P_+Q_);
    //std::copy_n(linp_.parameters_.begin(),P_,params.begin());
    //for(int i = 0; i < iter; i++){
    //  for(int j = 0; j < Q_; j++){
    //    params[P_+j] = u_(j,i);
    //  }
    //  matrix_matrix hess = vcalc_.jacobian_and_hessian(params,XZ,Map<MatrixXd>(offset_.data(),offset_.size(),1));
    //  result += hess.mat1;
    //}
    //result *= (1.0/iter);
    //return result;
    W.update();
    MatrixXd XtXW = (model.linear_predictor.X()).transpose() * W.W_.asDiagonal() * model.linear_predictor.X();
    MatrixXd ZL = model.covariance.ZL();
    MatrixXd XtWZL = (model.linear_predictor.X()).transpose() * W.W_.asDiagonal() * ZL;
    MatrixXd ZLWLZ = ZL.transpose() * W.W_.asDiagonal() * ZL;
    ZLWLZ += MatrixXd::Identity(model.covariance.Q(), model.covariance.Q());
    MatrixXd infomat(model.linear_predictor.P() + model.covariance.Q(), model.linear_predictor.P() + model.covariance.Q());
    infomat.topLeftCorner(model.linear_predictor.P(), model.linear_predictor.P()) = XtXW;
    infomat.topRightCorner(model.linear_predictor.P(), model.covariance.Q()) = XtWZL;
    infomat.bottomLeftCorner(model.covariance.Q(), model.linear_predictor.P()) = XtWZL.transpose();
    infomat.bottomRightCorner(model.covariance.Q(), model.covariance.Q()) = ZLWLZ;
    return infomat;
}

template<typename modeltype>
inline MatrixXd glmmr::ModelMatrix<modeltype>::sandwich_matrix() {
    MatrixXd infomat = observed_information_matrix();
    infomat = infomat.llt().solve(MatrixXd::Identity(model.linear_predictor.P() + model.covariance.Q(), model.linear_predictor.P() + model.covariance.Q()));
    infomat.conservativeResize(model.linear_predictor.P(), model.linear_predictor.P());
    MatrixXd zuOffset = re.Zu();
    zuOffset.colwise() += model.data.offset;
    MatrixXd J = model.calc.jacobian(model.linear_predictor.parameters, model.linear_predictor.Xdata, zuOffset);
    MatrixXd sandwich = infomat * (J * J.transpose()) * infomat;
    return sandwich;
}

template<typename modeltype>
inline std::vector<MatrixXd> glmmr::ModelMatrix<modeltype>::sigma_derivatives() {
    std::vector<MatrixXd> derivs;
    model.covariance.derivatives(derivs, 2);
    return derivs;
}

template<typename modeltype>
inline MatrixXd glmmr::ModelMatrix<modeltype>::information_matrix_theta() {
    int n = model.n();
    std::vector<MatrixXd> derivs;
    model.covariance.derivatives(derivs, 1);
    int R = model.covariance.npar();
    int Rmod = model.family.family == "gaussian" ? R + 1 : R;
    MatrixXd M = information_matrix();
    M = M.llt().solve(MatrixXd::Identity(M.rows(), M.cols()));
    MatrixXd SigmaInv = Sigma(true);
    MatrixXd Z = model.covariance.Z();
    MatrixXd X = model.linear_predictor.X();
    MatrixXd SigX = SigmaInv * X;
    MatrixXd M_theta = MatrixXd::Zero(Rmod, Rmod);
    MatrixXd partial1(model.n(), model.n());
    MatrixXd partial2(model.n(), model.n());
    glmmr::MatrixField<MatrixXd> S;
    glmmr::MatrixField<MatrixXd> P;
    glmmr::MatrixField<MatrixXd> Q;
    int counter = 0;
    for (int i = 0; i < Rmod; i++) {
        if (i < R) {
            partial1 = Z * derivs[1 + i] * Z.transpose();
        }
        else {
            partial1 = MatrixXd::Identity(n, n);
            if ((model.data.weights != 1).any())partial1 = model.data.weights.inverse().matrix().asDiagonal();
        }
        P.add(-1 * SigX.transpose() * partial1 * SigX);
        for (int j = i; j < Rmod; j++) {
            if (j < R) {
                partial2 = Z * derivs[1 + j] * Z.transpose();
            }
            else {
                partial2 = MatrixXd::Identity(n, n);
                if ((model.data.weights != 1).any())partial2 = model.data.weights.inverse().matrix().asDiagonal();
            }
            Q.add(X.transpose() * SigmaInv * partial1 * SigmaInv * partial2 * SigX);
            S.add(SigmaInv * partial1 * SigmaInv * partial2);
        }
    }
    counter = 0;
    for (int i = 0; i < Rmod; i++) {
        for (int j = i; j < Rmod; j++) {
            M_theta(i, j) = 0.5 * (S(counter).trace()) - (M * Q(counter)).trace() + 0.5 * ((M * P(i) * M * P(j)).trace());
            if (i != j)M_theta(j, i) = M_theta(i, j);
            counter++;
        }
    }
    M_theta = M_theta.llt().solve(MatrixXd::Identity(Rmod, Rmod));
    return M_theta;
}

template<typename modeltype>
inline kenward_data glmmr::ModelMatrix<modeltype>::kenward_roger() {
    int n = model.n();
    std::vector<MatrixXd> derivs;
    model.covariance.derivatives(derivs, 2);
    int R = model.covariance.npar();
    int Rmod = model.family.family == "gaussian" ? R + 1 : R;

    MatrixXd M = information_matrix();
    M = M.llt().solve(MatrixXd::Identity(M.rows(), M.cols()));
    MatrixXd M_new(M);
    MatrixXd SigmaInv = Sigma(true);
    MatrixXd Z = model.covariance.Z();
    MatrixXd X = model.linear_predictor.X();
    MatrixXd SigX = SigmaInv * X;
    MatrixXd M_theta = MatrixXd::Zero(Rmod, Rmod);
    MatrixXd partial1(model.n(), model.n());
    MatrixXd partial2(model.n(), model.n());
    MatrixXd meat = MatrixXd::Zero(SigX.cols(), SigX.cols());
    glmmr::MatrixField<MatrixXd> P;
    glmmr::MatrixField<MatrixXd> Q;
    glmmr::MatrixField<MatrixXd> RR;
    glmmr::MatrixField<MatrixXd> S;
    int counter = 0;
    for (int i = 0; i < Rmod; i++) {
        if (i < R) {
            partial1 = Z * derivs[1 + i] * Z.transpose();
        }
        else {
            partial1 = MatrixXd::Identity(n, n);
            if ((model.data.weights != 1).any())partial1 = model.data.weights.inverse().matrix().asDiagonal();
        }
        P.add(-1 * SigX.transpose() * partial1 * SigX);
        for (int j = i; j < Rmod; j++) {
            if (j < R) {
                partial2 = Z * derivs[1 + j] * Z.transpose();
            }
            else {
                partial2 = MatrixXd::Identity(n, n);
                if ((model.data.weights != 1).any())partial2 = model.data.weights.inverse().matrix().asDiagonal();
            }
            S.add(SigmaInv * partial1 * SigmaInv * partial2);
            Q.add(X.transpose() * SigmaInv * partial1 * SigmaInv * partial2 * SigX);
            if (i < R && j < R) {
                int scnd_idx = i + j * (R - 1) - j * (j - 1) / 2;
                RR.add(SigX.transpose() * Z * derivs[R + 1 + scnd_idx] * Z.transpose() * SigX);
            }
        }
    }
    counter = 0;
    for (int i = 0; i < Rmod; i++) {
        for (int j = i; j < Rmod; j++) {
            M_theta(i, j) = 0.5 * (S(counter).trace()) - (M * Q(counter)).trace() + 0.5 * ((M * P(i) * M * P(j)).trace());
            if (i != j)M_theta(j, i) = M_theta(i, j);
            counter++;
        }
    }
    M_theta = M_theta.llt().solve(MatrixXd::Identity(Rmod, Rmod));
    for (int i = 0; i < (Rmod - 1); i++) {
        if (i < R) {
            partial1 = Z * derivs[1 + i] * Z.transpose();
        }
        else {
            partial1 = MatrixXd::Identity(n, n);
            if ((model.data.weights != 1).any())partial1 = model.data.weights.inverse().matrix().asDiagonal();
        }
        for (int j = (i + 1); j < Rmod; j++) {
            if (j < R) {
                partial2 = Z * derivs[1 + j] * Z.transpose();
            }
            else {
                partial2 = MatrixXd::Identity(n, n);
            }
            int scnd_idx = i + j * (Rmod - 1) - j * (j - 1) / 2;
            meat += M_theta(i, j) * (Q(scnd_idx) + Q(scnd_idx).transpose() - P(i) * M * P(j) - P(j) * M * P(i));//(SigX.transpose()*partial1*PG*partial2*SigX);//
            if (i < R && j < R) {
                scnd_idx = i + j * (R - 1) - j * (j - 1) / 2;
                meat -= 0.5 * M_theta(i, j) * (RR(scnd_idx));
            }
        }
    }
    for (int i = 0; i < Rmod; i++) {
        if (i < R) {
            partial1 = Z * derivs[1 + i] * Z.transpose();
        }
        else {
            partial1 = MatrixXd::Identity(n, n);
            if ((model.data.weights != 1).any())partial1 = model.data.weights.inverse().matrix().asDiagonal();
        }
        int scnd_idx = i + i * (Rmod - 1) - i * (i - 1) / 2;
        meat += M_theta(i, i) * (Q(scnd_idx) - P(i) * M * P(i));
        if (i < R) {
            scnd_idx = i + i * (R - 1) - i * (i - 1) / 2;
            meat -= 0.25 * M_theta(i, i) * RR(scnd_idx);
        }
    }
    M_new = M + 2 * M * meat * M;

    // degrees of freedom correction

    kenward_data out(model.linear_predictor.P(), model.linear_predictor.P(), Rmod, Rmod);
    out.vcov_beta = M_new;
    out.vcov_theta = M_theta;

    double a1, a2, B, g, c1, c2, c3, v0, v1, v2, rhotop, rho;
    int mult = 1;
    VectorXd L = VectorXd::Zero(model.linear_predictor.P());
    MatrixXd Theta(model.linear_predictor.P(), model.linear_predictor.P());
    for (int p = 0; p < L.size(); p++) {
        L.setZero();
        L(p) = 1;
        double vlb = L.transpose() * M * L;
        Theta = (1 / vlb) * (L * L.transpose());
        Theta = Theta * M;
        a1 = 0;
        a2 = 0;
        for (int i = 0; i < Rmod; i++) {
            for (int j = i; j < Rmod; j++) {
                mult = i == j ? 1 : 2;
                a1 += mult * M_theta(i, j) * (Theta * P(i) * M).trace() * (Theta * P(j) * M).trace();
                a2 += mult * M_theta(i, j) * (Theta * P(i) * M * Theta * P(j) * M).trace();
            }
        }
        B = (a1 + 6 * a2) * 0.5;
        g = (2 * a1 - 5 * a2) / (3 * a2);
        c1 = g / (3 + 2 * (1 - g));
        c2 = (1 - g) / (3 + 2 * (1 - g));
        c3 = (3 - g) / (3 + 2 * (1 - g));
        v0 = abs(1 + c1 * B) < 1e-10 ? 0 : 1 + c1 * B;
        v1 = 1 - c2 * B;
        v2 = 1 / (1 - c3 * B);
        rhotop = abs(1 - a2) < 1e-10 && abs(v1) < 1e-10 ? 1.0 : (1 - a2) / v1;
        rho = rhotop * rhotop * v0 * v2;
        out.dof(p) = 4 + 3 / (rho - 1);
        out.lambda(p) = (1 - a2) * out.dof(p) / (out.dof(p) - 2);
    }

    return out;
}

template<typename modeltype>
inline MatrixXd glmmr::ModelMatrix<modeltype>::linpred() {
    return (re.zu_.colwise() + (model.linear_predictor.xb() + model.data.offset));
}

template<typename modeltype>
inline vector_matrix glmmr::ModelMatrix<modeltype>::b_score() {
    MatrixXd zuOffset = re.Zu();
    zuOffset.colwise() += model.data.offset;
    matrix_matrix hess = model.calc.jacobian_and_hessian(model.linear_predictor.parameters, model.linear_predictor.Xdata, zuOffset);
    vector_matrix out(hess.mat1.rows());
    out.mat = hess.mat1;
    out.mat *= -1.0;
    out.vec = hess.mat2.rowwise().sum();
    return out;
}

template<typename modeltype>
inline matrix_matrix glmmr::ModelMatrix<modeltype>::hess_and_grad() {
    MatrixXd zuOffset = re.Zu();
    zuOffset.colwise() += model.data.offset;
    matrix_matrix hess = model.calc.jacobian_and_hessian(model.linear_predictor.parameters, model.linear_predictor.Xdata, zuOffset);
    return hess;
}

template<typename modeltype>
inline vector_matrix glmmr::ModelMatrix<modeltype>::re_score() {
    VectorXd xbOffset = model.linear_predictor.xb() + model.data.offset;
    matrix_matrix hess = model.vcalc.jacobian_and_hessian(dblvec(re.u(false).col(0).data(), re.u(false).col(0).data() + re.u(false).rows()), sparse_to_dense(re.ZL, false), Map<MatrixXd>(xbOffset.data(), xbOffset.size(), 1));
    vector_matrix out(model.covariance.Q());
    hess.mat1 *= -1.0;
    out.mat = hess.mat1 + MatrixXd::Identity(model.covariance.Q(), model.covariance.Q());
    out.vec = hess.mat2.rowwise().sum();
    out.vec -= re.u(false).col(0);
    return out;
}

template<typename modeltype>
inline VectorXd glmmr::ModelMatrix<modeltype>::log_gradient(const VectorXd& v,
    bool beta) {
    ArrayXd size_n_array = model.xb();
    ArrayXd size_q_array = ArrayXd::Zero(model.covariance.Q());
    ArrayXd size_p_array = ArrayXd::Zero(model.linear_predictor.P());
    sparse ZLt = re.ZL;
    ZLt.transpose();
    size_n_array += (re.ZL * v).array();

    if (beta) {
        VectorXd zuOffset = re.ZL * v;
        zuOffset += model.data.offset;
        MatrixXd J = model.calc.jacobian(model.linear_predictor.parameters, model.linear_predictor.Xdata, zuOffset);
        size_p_array = J.transpose().rowwise().sum().array();
    }
    else {
        switch (model.family.flink) {
        case 1:
        {
            size_n_array = size_n_array.exp();
            if (!beta) {
                size_n_array = model.data.y.array() - size_n_array;
                size_q_array = ZLt * size_n_array - v.array();
            }
            else {
                size_p_array += (model.linear_predictor.X().transpose() * (model.data.y - size_n_array.matrix())).array();
            }
            break;
        }
        case 2:
        {
            size_n_array = size_n_array.inverse();
            size_n_array = model.data.y.array() * size_n_array;
            size_n_array -= ArrayXd::Ones(model.n());
            if (beta) {
                size_p_array += (model.linear_predictor.X().transpose() * size_n_array.matrix()).array();
            }
            else {
                size_q_array = ZLt * size_n_array - v.array();
            }
            break;
        }
        case 3: case 13:
        {
            ArrayXd logitxb = model.xb().array().exp();
            logitxb += 1;
            logitxb = logitxb.inverse();
            logitxb *= model.xb().array().exp();
            size_n_array = model.data.y.array() * (ArrayXd::Constant(model.n(), 1) - logitxb) - (model.data.variance - model.data.y.array()) * logitxb;
            if (beta) {
                size_p_array += (model.linear_predictor.X().transpose() * size_n_array.matrix()).array();
            }
            else {
                size_q_array = ZLt * size_n_array - v.array();
            }
            break;
        }
        case 4: case 14:
        {
            ArrayXd logitxb = model.xb().array().exp();
            logitxb += 1;
            logitxb = logitxb.inverse();
            logitxb *= model.xb().array().exp();
            size_n_array = (model.data.variance - model.data.y.array()) * logitxb;
            size_n_array += model.data.y.array();
            if (beta) {
                size_p_array += (model.linear_predictor.X().transpose() * size_n_array.matrix()).array();
            }
            else {
                size_q_array = ZLt * size_n_array - v.array();
            }
            break;
        }
        case 5: case 15:
        {
            size_n_array = size_n_array.inverse();
            size_n_array *= model.data.y.array();
            ArrayXd n_array2 = ArrayXd::Constant(model.n(), 1.0) - model.xb().array();
            n_array2 = n_array2.inverse();
            n_array2 *= (model.data.variance - model.data.y.array());
            size_n_array -= n_array2;
            if (beta) {
                size_p_array += (model.linear_predictor.X().transpose() * size_n_array.matrix()).array();
            }
            else {
                size_q_array = ZLt * size_n_array - v.array();
            }
            break;
        }
        case 6: case 16:
        {
            ArrayXd n_array2(model.n());
            boost::math::normal norm(0, 1);
#pragma omp parallel for    
            for (int i = 0; i < model.n(); i++) {
                size_n_array(i) = (double)pdf(norm, size_n_array(i)) / ((double)cdf(norm, size_n_array(i)));
                n_array2(i) = -1.0 * (double)pdf(norm, size_n_array(i)) / (1 - (double)cdf(norm, size_n_array(i)));
            }
            size_n_array = model.data.y.array() * size_n_array + (model.data.variance - model.data.y.array()) * n_array2;
            if (beta) {
                size_p_array += (model.linear_predictor.X().transpose() * size_n_array.matrix()).array();
            }
            else {
                size_q_array = ZLt * size_n_array - v.array();
            }
            break;
        }
        case 7:
        {
            if (beta) {
                size_n_array -= model.data.y.array();
                size_n_array *= -1;
                size_n_array *= model.data.weights;
                size_p_array += ((1.0 / (model.data.var_par)) * (model.linear_predictor.X().transpose() * size_n_array.matrix())).array();
            }
            else {
                size_n_array = model.data.y.array() - size_n_array;
                size_n_array *= model.data.weights;
                size_q_array = (ZLt * size_n_array) - v.array();
                size_q_array *= 1.0 / (model.data.var_par);
            }
            break;
        }
        case 8:
        {
            if (beta) {
                size_n_array -= model.data.y.array();
                size_n_array *= -1;
                size_n_array *= model.data.weights;
                size_p_array += ((1.0 / (model.data.var_par)) * (model.linear_predictor.X().transpose() * (model.data.y - size_n_array.matrix()))).array();
            }
            else {
                size_n_array = model.data.y.array() - size_n_array;
                size_n_array *= model.data.weights;
                size_q_array = ZLt * size_n_array - v.array();
                size_q_array *= 1.0 / (model.data.var_par);
            }
            break;
        }
        case 9:
        {
            size_n_array *= -1.0;
            size_n_array = size_n_array.exp();
            if (beta) {
                size_p_array += (model.linear_predictor.X().transpose() * (model.data.y.array() * size_n_array - 1).matrix() * model.data.var_par).array();
            }
            else {
                size_n_array *= model.data.y.array();
                size_q_array = ZLt * size_n_array - v.array();
                size_q_array *= model.data.var_par;
            }
            break;
        }
        case 10:
        {
            size_n_array = size_n_array.inverse();
            if (beta) {
                size_p_array += (model.linear_predictor.X().transpose() * (size_n_array.matrix() - model.data.y) * model.data.var_par).array();
            }
            else {
                size_n_array -= model.data.y.array();
                size_q_array = ZLt * size_n_array - v.array();
                size_q_array *= model.data.var_par;
            }
            break;
        }
        case 11:
        {
            size_n_array = size_n_array.inverse();
            if (beta) {
                size_p_array += (model.linear_predictor.X().transpose() * ((model.data.y.array() * size_n_array * size_n_array).matrix() - size_n_array.matrix()) * model.data.var_par).array();
            }
            else {
                size_n_array *= (model.data.y.array() * size_n_array - ArrayXd::Ones(model.n()));
                size_q_array = ZLt * size_n_array - v.array();
                size_q_array *= model.data.var_par;
            }
            break;
        }
        case 12:
        {
#pragma omp parallel for 
            for (int i = 0; i < model.n(); i++) {
                size_n_array(i) = exp(size_n_array(i)) / (exp(size_n_array(i)) + 1);
                size_n_array(i) = (size_n_array(i) / (1 + exp(size_n_array(i)))) * model.data.var_par * (log(model.data.y(i)) - log(1 - model.data.y(i)) - boost::math::digamma(size_n_array(i) * model.data.var_par) + boost::math::digamma((1 - size_n_array(i)) * model.data.var_par));
            }
            if (beta) {
                size_p_array += (model.linear_predictor.X().transpose() * size_n_array.matrix()).array();
            }
            else {
                size_q_array = ZLt * size_n_array - v.array();
            }
            break;
        }

        }
    }
    // we can use autodiff here, but the above method is faster
    // else {
    //   VectorXd xbOffset_ = linpred_.xb() + offset_;
    //   MatrixXd J = vcalc_.jacobian(dblvec(v.data(),v.data()+v.size()),
    //                                               sparse_to_dense(ZL_,false),
    //                                               xbOffset_);
    //   size_q_array = (J.transpose().rowwise().sum() - v).array();
    // }
    return beta ? size_p_array.matrix() : size_q_array.matrix();
}
