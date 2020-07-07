//
// Created by Alan Pryor on 2/1/17.
//

#ifndef INTERP_INTERP_H_H
#define INTERP_INTERP_H_H
#include <iostream>
#include <cmath>
#include "Array2D.h"

// Declarations
using namespace std;
using namespace GENFIRE;

// Interpolate position x,y from 2D array stored in data, which is assumed to lie on an integer coordinate system.
// Values of x and y that lie outside of the dimensions of data are set to 0.
template <typename T>
void interp2(const Array2D<T>& data, const double* x, const double* y, size_t N, typename T::value_type* result);

template <typename T>
void interp2(const T* data,const size_t& nrows, const size_t& ncols,  const double* x, const double* y, size_t N, T* result);


// Implementations
template <typename T>
void interp2(const Array2D<T>& data, const double* x, const double* y, size_t N, typename T::value_type* result){
    for (auto i = 0; i < N; ++i) {
        const int x_1 = (int) std::floor(x[i]);
        const int x_2 = x_1 + 1;
        const int y_1 = (int) std::floor(y[i]);
        const int y_2 = y_1 + 1;
        const typename T::value_type &f_11 = data.at(x_1, y_1);
        const typename T::value_type &f_12 = data.at(x_1, y_2);
        const typename T::value_type &f_21 = data.at(x_2, y_1);
        const typename T::value_type &f_22 = data.at(x_2, y_2);

        typename T::value_type w_x1 = x_2 - x[i];
        typename T::value_type w_x2 = x[i] - x_1;
        typename T::value_type w_y1 = y_2 - y[i];
        typename T::value_type w_y2 = y[i] - y_1;

        double a;
        double b;
        a = f_11 * w_x1 + f_21 * w_x2;
        b = f_12 * w_x1 + f_22 * w_x2;
        result[i] = a * w_y1 + b * w_y2;
    }
}


template <typename T>
void interp2(const T* data, const size_t& nrows, const size_t& ncols, const double* x, const double* y,size_t N, T* result){
for (auto i = 0; i < N; ++i) {
const int x_1 = (int) std::floor(x[i]);
const int x_2 = x_1 + 1;
const int y_1 = (int) std::floor(y[i]);
const int y_2 = y_1 + 1;

const T& f_11 = data[x_1 + y_1*nrows];
const T& f_12 = data[x_1 + y_2*nrows];
const T& f_21 = data[x_2 + y_1*nrows];
const T& f_22 = data[x_2 + y_2*nrows];

 T w_x1 = x_2 - x[i];
 T w_x2 = x[i] - x_1;
 T w_y1 = y_2 - y[i];
 T w_y2 = y[i] - y_1;

double a;
double b;
a = f_11 * w_x1 + f_21 * w_x2;
b = f_12 * w_x1 + f_22 * w_x2;
result[i] = a * w_y1 + b * w_y2;
}
}


#endif //INTERP_INTERP_H_H
