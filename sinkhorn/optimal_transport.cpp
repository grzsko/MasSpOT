/*
 * =====================================================================================
 *
 *       Filename:  optimal_transport.cpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  11/04/19 10:55:59
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (),
 *   Organization:
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <iostream>
#include <Eigen/Dense>
#include "optimal_transport.h"

using namespace Eigen;
typedef Map<VectorXd> MapT;

// TODO how to do it passable via cffi?
// const double lambda = 1;
// const double epsilon = 0.1;
// const double tol = 1e-05;
// const double threshold = 1e+02;
// const double max_iter = 1000;

VectorXd proxdiv_TV(VectorXd s, VectorXd u, VectorXd p, double lambda,
                    double epsilon) {
    ArrayXd e_lam = ((lambda - u.array()) / epsilon).exp();
    ArrayXd e_minus_lam = (-(lambda + u.array()) / epsilon).exp();
    return e_lam.cwiseMin(e_minus_lam.cwiseMax(p.array() / s.array())).matrix();
}

VectorXd proxdiv_i(VectorXd s, VectorXd u, VectorXd p, double lambda,
                   double epsilon) {
    return (p.array() / s.array()).matrix();
}

MatrixXd calc_distance_matrix(VectorXd mzs1, VectorXd mzs2) {
    return (mzs1 * VectorXd::Ones(mzs2.rows()).transpose() -
                VectorXd::Ones(mzs1.rows()) * mzs2.transpose()).cwiseAbs();
}

MatrixXd calc_transport_plan(VectorXd ints1, VectorXd ints2, MatrixXd dists,
                             double lambda, double epsilon, double tol,
                             double threshold, double max_iter) {
    int n1 = ints1.rows();
    int n2 = ints2.rows();
    MatrixXd K_0 = (-dists / epsilon).array().exp().matrix();
    MatrixXd K = K_0;

    VectorXd b = VectorXd::Ones(n2);
    VectorXd a = VectorXd::Zero(n1);
    VectorXd u = VectorXd::Zero(n1);
    VectorXd v = VectorXd::Zero(n2);
	VectorXd a_new, b_new;
    double coverging_val;

    int tick = 0;
    bool coverged = false;

    do {
        a_new = proxdiv_TV(K * b, u, ints1, lambda, epsilon); // * is matrix-matrix
        b_new = proxdiv_TV(K.transpose() * a_new, v, ints2, lambda, epsilon);

        if ((((a_new.array().log().abs() - threshold) > 0).any()) or
            (((b_new.array().log().abs() - threshold) > 0).any())) {
            // Stabilizing
            u += epsilon * a_new.array().log().matrix();
            v += epsilon * b_new.array().log().matrix();
            // below is exp((u_i + v_j - dists_{ij}) / eps) forall i,j
            K = (u / epsilon).array().exp().matrix().asDiagonal() * K_0 *
                (v / epsilon).array().exp().matrix().asDiagonal();
        }

        coverging_val = std::max((a_new - a).maxCoeff(),
                                 (b_new - b).maxCoeff());
        tick++;
        coverged = coverging_val < tol or tick >= max_iter;
        a = a_new;
        b = b_new;
    } while (not coverged);
    // below is (a_i * K_{ij} * b_j)_{ij}
    return a.asDiagonal() * K * b.asDiagonal();
}

double calc_distance_from_plan(MatrixXd transport_plan, MatrixXd distances,
                               VectorXd ints1, VectorXd ints2, double lambda) {
    double transport = (transport_plan.array() * distances.array()).sum();
    // ^^^ coeff-wise multiplication ^^^

    double trash = 0;
    trash += lambda * ((transport_plan.rowwise().sum() - ints1).cwiseAbs().sum());
    trash += lambda * ((transport_plan.colwise().sum().transpose() - ints2).cwiseAbs().sum());
    // ^^ colwise sum is row, not column
    return transport + trash;
}

double calc_distance_cpp(double* mzs1, double* ints1, int len1, double* mzs2,
                         double* ints2, int len2, double lambda,
                         double epsilon, double tol, double threshold,
                         double max_iter) {
    /*  TODO description
     *
     *  ints1, and ints2 should be normalized */
    /*
     * TODO sharing memory causes segfault, so maybe copying to matrix is
     * better?
     */

    MapT mzs1_map(mzs1, len1);
    MapT mzs2_map(mzs2, len2);
    MapT ints1_map(ints1, len1);
    MapT ints2_map(ints2, len2);
    MatrixXd dists = calc_distance_matrix(mzs1_map, mzs2_map);
    MatrixXd transport_plan = calc_transport_plan(ints1_map, ints2_map, dists,
            lambda, epsilon, tol, threshold, max_iter);
    return calc_distance_from_plan(transport_plan, dists, ints1_map, ints2_map,
                                   lambda);
}


double calc_distance_c(double* mzs1, double* ints1, int len1, double* mzs2,
                       double* ints2, int len2, double lambda, double epsilon,
                       double tol, double threshold, double max_iter) {
    return calc_distance_cpp(mzs1, ints1, len1, mzs2, ints2, len2, lambda,
                             epsilon, tol, threshold, max_iter);
}

