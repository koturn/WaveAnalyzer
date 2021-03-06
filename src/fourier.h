/*!
 * @brief fourier.cのヘッダファイル
 *
 * @author  koturn 0;
 * @date    2013 05/01
 * @file    fourier.h
 * @version 0.1
 */
#ifndef FOURIER_H
#define FOURIER_H

#include "compatibility.h"


void auto_dft(double *restrict X_real, double *restrict X_imag, const double *restrict x_real, const double *restrict x_imag, unsigned int N);
void dft(double *restrict X_real, double *restrict X_imag, const double *restrict x_real, const double *restrict x_imag, unsigned int N);
void fft(double *restrict X_real, double *restrict X_imag, const double *restrict x_real, const double *restrict x_imag, unsigned int N);
void auto_idft(double *restrict x_real, double *restrict x_imag, const double *restrict X_real, const double *restrict X_imag, unsigned int N);
void idft(double *restrict x_real, double *restrict x_imag, const double *restrict X_real, const double *restrict X_imag, unsigned int N);
void ifft(double *restrict x_real, double *restrict x_imag, const double *restrict X_real, const double *restrict X_imag, unsigned int N);
void calc_power_spectrum(double *restrict power, const double *restrict X_real, const double *restrict X_imag, unsigned int N);


#endif
