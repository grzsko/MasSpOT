/*
 * =====================================================================================
 *
 *       Filename:  optimal_transport.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/04/19 10:54:46
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#ifndef optimal_transport_h__
#define optimal_transport_h__
extern "C" {
    double calc_distance_c(double*, double*, int, double*, double*, int,
            double, double, double, double, double);
}


#endif
