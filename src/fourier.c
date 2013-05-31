/*!
 * @brief フーリエ変換のための関数を提供する
 *
 * @author  koturn 0;
 * @date    2013 05/30
 * @file    fourier.c
 * @version 0.2
 */
#ifdef _MSC_VER
#  define _USE_MATH_DEFINES
#endif
#include <math.h>
#include <memory.h>
#include <stdlib.h>
#include "compatibility.h"
#include "fourier.h"


#ifdef __GNUC__
#  define __DO__          (          //!< GCC用のdo ~ while (0) のdoに当たる部分
#  define __WHILE_ZERO__  )          //!< GCC用のdo ~ while (0) のwhile (0) に当たる部分
#else
#  define __DO__          do         //!< do ~ while (0) の do (gcc以外で使用)
#  define __WHILE_ZERO__  while (0)  //!< do ~ while (0) の while (0) (gcc以外で使用)
#endif

#define SWAP(type, a, b)                  \
__DO__ {                                  \
  register type __tmp_swap_var__ = *(a);  \
  *(a) = *(b);                            \
  *(b) = __tmp_swap_var__;                \
} __WHILE_ZERO__

#define CALLOC(type, n)  ((type *)calloc((n), sizeof(type)))
#define FREE(ptr)        (free(ptr), (ptr) = NULL)


static void dft_ri(double *restrict X_real, double *restrict X_imag, const double *restrict x_real, const double *restrict x_imag, unsigned int N);
static void dft_r(double *restrict X_real, double *restrict X_imag, const double *restrict x_real, unsigned int N);
static void dft_i(double *restrict X_real, double *restrict X_imag, const double *restrict x_imag, unsigned int N);
static void idft_ri(double *restrict x_real, double *restrict x_imag, const double *restrict X_real, const double *restrict X_imag, unsigned int N);
static void idft_r(double *restrict x_real, double *restrict x_imag, const double *restrict X_real, unsigned int N);
static void idft_i(double *restrict x_real, double *restrict x_imag, const double *restrict X_imag, unsigned int N);
static void dbl_array_copy(double *restrict dst, const double *restrict src, unsigned int n);
__attribute__((const)) inline static unsigned int my_log2(unsigned int n);
__attribute__((const)) inline static unsigned int my_pow2(unsigned int n);


/*!
 * @brief 離散フーリエ変換を行う
 *
 * 標本化周波数により、高速フーリエ変換を行う
 * @param [out] X_real  周波数信号の実部(結果出力先)
 * @param [out] X_imag  周波数信号の虚部(結果出力先)
 * @param [in]  x_real  時間信号の実部(入力)
 * @param [in]  x_imag  時間信号の虚部(入力)
 * @param [in]  N       標本数
 * @see fft
 * @see dft
 */
void auto_ft(double *restrict X_real, double *restrict X_imag, const double *restrict x_real, const double *restrict x_imag, unsigned int N) {
  if (!(N & (N - 1))) {  // 2の累乗なら
    fft(X_real, X_imag, x_real, x_imag, N);
  } else {
    dft(X_real, X_imag, x_real, x_imag, N);
  }
}


/*!
 * @brief 通常の離散フーリエ変換を行う
 *
 * 入力信号にNULLを許容し、NULLの入力信号があった場合、<br>
 * NULLの信号を全て0の信号として扱い、計算処理を軽くした離散フーリエ変換を行う。
 * @param [out] X_real  周波数信号の実部(結果出力先)
 * @param [out] X_imag  周波数信号の虚部(結果出力先)
 * @param [in]  x_real  時間信号の実部(入力)
 * @param [in]  x_imag  時間信号の虚部(入力)
 * @param [in]  N       標本数
 * @see auto_ft
 * @see dft_ri
 * @see dft_r
 * @see dft_i
 */
void dft(double *restrict X_real, double *restrict X_imag, const double *restrict x_real, const double *restrict x_imag, unsigned int N) {
  if (x_real != NULL && x_imag != NULL) {
    dft_ri(X_real, X_imag, x_real, x_imag, N);
  } else if (x_imag == NULL) {
    dft_r(X_real, X_imag, x_real, N);
  } else if (x_real == NULL) {
    dft_i(X_real, X_imag, x_imag, N);
  }
}


/*!
 * @brief 実部信号と虚部信号がある場合の離散フーリエ変換を行う
 * @param [out] X_real  周波数信号の実部(結果出力先)
 * @param [out] X_imag  周波数信号の虚部(結果出力先)
 * @param [in]  x_real  時間信号の実部(入力)
 * @param [in]  x_imag  時間信号の虚部(入力)
 * @param [in]  N       標本数
 * @see dft
 */
static void dft_ri(double *restrict X_real, double *restrict X_imag, const double *restrict x_real, const double *restrict x_imag, unsigned int N) {
  unsigned int k;
  #pragma omp parallel for
  for (k = 0; k < N; k++) {
    unsigned int n;
    for (n = 0; n < N; n++) {
      double w_real = cos(2.0 * M_PI * k * n / N);
      double w_imag = sin(2.0 * M_PI * k * n / N);
      X_real[k] += w_real * x_real[n] + w_imag * x_imag[n];  // X(k)の実数部
      X_imag[k] += w_real * x_imag[n] - w_imag * x_real[n];  // X(k)の虚数部
    }
  }
}

/*!
 * @brief 実部信号のみがある場合の離散フーリエ変換を行う
 * @param [out] X_real  周波数信号の実部(結果出力先)
 * @param [out] X_imag  周波数信号の虚部(結果出力先)
 * @param [in]  x_real  時間信号の実部(入力)
 * @param [in]  N       標本数
 * @see dft
 */
static void dft_r(double *restrict X_real, double *restrict X_imag, const double *restrict x_real, unsigned int N) {
  unsigned int k;
  #pragma omp parallel for
  for (k = 0; k < N; k++) {
    unsigned int n;
    for (n = 0; n < N; n++) {
      double w_real = cos(2.0 * M_PI * k * n / N);
      double w_imag = sin(2.0 * M_PI * k * n / N);
      X_real[k] +=  w_real * x_real[n];  // X(k)の実数部
      X_imag[k] += -w_imag * x_real[n];  // X(k)の虚数部
    }
  }
}


/*!
 * @brief 虚部信号のみがある場合の離散フーリエ変換を行う
 * @param [out] X_real  周波数信号の実部(結果出力先)
 * @param [out] X_imag  周波数信号の虚部(結果出力先)
 * @param [in]  x_imag  時間信号の虚部(入力)
 * @param [in]  N       標本数
 * @see dft
 */
static void dft_i(double *restrict X_real, double *restrict X_imag, const double *restrict x_imag, unsigned int N) {
  unsigned int k;
  #pragma omp parallel for
  for (k = 0; k < N; k++) {
    unsigned int n;
    for (n = 0; n < N; n++) {
      double w_real = cos(2.0 * M_PI * k * n / N);
      double w_imag = sin(2.0 * M_PI * k * n / N);
      X_real[k] += w_imag * x_imag[n];  // X(k)の実数部
      X_imag[k] += w_real * x_imag[n];  // X(k)の虚数部
    }
  }
}


/*!
 * @brief 高速離散フーリエ変換を行う
 * @param [out] X_real  周波数信号の実部(結果出力先)
 * @param [out] X_imag  周波数信号の虚部(結果出力先)
 * @param [in]  x_real  時間信号の実部(入力)
 * @param [in]  x_imag  時間信号の虚部(入力)
 * @param [in]  N       標本数
 * @see auto_dft
 */
void fft(double *restrict X_real, double *restrict X_imag, const double *restrict x_real, const double *restrict x_imag, unsigned int N) {
  unsigned int  k, stage;
  unsigned int *restrict index;
  unsigned int  number_of_stage = my_log2(N);  // FFTの段数

  // Copy
  #pragma omp parallel sections
  {
    #pragma omp parallel section
    dbl_array_copy(X_real, x_real, N);
    #pragma omp parallel section
    dbl_array_copy(X_imag, x_imag, N);
  }

  // バタフライ計算(順番が関わるので、並列化できない)
  for (stage = 1; stage <= number_of_stage; stage++) {
    unsigned int i;
    unsigned int p1 = my_pow2(stage - 1);
    for (i = 0; i < p1; i++) {
      unsigned int j;
      unsigned int p2 = my_pow2(number_of_stage - stage);
      for (j = 0; j < p2; j++) {
        int n = my_pow2(number_of_stage - stage + 1) * i + j;
        int m = my_pow2(number_of_stage - stage) + n;
        int r = my_pow2(stage - 1) * j;
        double a_real = X_real[n];
        double a_imag = X_imag[n];
        double b_real = X_real[m];
        double b_imag = X_imag[m];
        if (stage < number_of_stage) {
          double w_real = cos((2.0 * M_PI * r) / N);
          double w_imag = sin((2.0 * M_PI * r) / N);
          X_real[n] = a_real + b_real;
          X_imag[n] = a_imag + b_imag;
          X_real[m] = (a_real - b_real) * w_real + (a_imag - b_imag) * w_imag;
          X_imag[m] = (a_imag - b_imag) * w_real - (a_real - b_real) * w_imag;
        } else {
          X_real[n] = a_real + b_real;
          X_imag[n] = a_imag + b_imag;
          X_real[m] = a_real - b_real;
          X_imag[m] = a_imag - b_imag;
        }
      }
    }
  }

  // インデックスの並び替えのためのテーブルの作成
  index = CALLOC(unsigned int, N);
  for (stage = 1; stage <= number_of_stage; stage++) {
    unsigned int i;
    for (i = 0; i < my_pow2(stage - 1); i++) {
      index[my_pow2(stage - 1) + i] = index[i] + my_pow2(number_of_stage - stage);
    }
  }

  // インデックスの並び替え
  #pragma omp parallel for
  for (k = 0; k < N; k++) {
    if (index[k] > k) {
      SWAP(double, &X_real[index[k]], &X_real[k]);
      SWAP(double, &X_imag[index[k]], &X_imag[k]);
    }
  }
  FREE(index);
}


/*!
 * @brief 逆離散フーリエ変換を行う
 *
 * 標本化周波数により、逆高速フーリエ変換を行う
 * @param [out] x_real  時間信号の実部(結果出力先)
 * @param [out] x_imag  時間信号の虚部(結果出力先)
 * @param [in]  X_real  周波数信号の実部(入力)
 * @param [in]  X_imag  周波数信号の虚部(入力)
 * @param [in]  N       標本数
 * @see ifft
 * @see idft
 */
void auto_ift(double *restrict x_real, double *restrict x_imag, const double *restrict X_real, const double *restrict X_imag, unsigned int N) {
  if (X_real != NULL && X_imag != NULL) {
    idft_ri(x_real, x_imag, X_real, X_imag, N);
  } else if (X_imag == NULL) {
    idft_r(x_real, x_imag, X_real, N);
  } else if (x_real == NULL) {
    idft_i(x_real, x_imag, X_imag, N);
  }
}


/*!
 * @brief 通常の逆離散フーリエ変換を行う
 *
 * 入力信号にNULLを許容し、NULLの入力信号があった場合、<br>
 * NULLの信号を全て0の信号として扱い、計算処理を軽くした逆離散フーリエ変換を行う。
 * 標本化周波数により、逆高速フーリエ変換を行う
 * @param [out] x_real  時間信号の実部(結果出力先)
 * @param [out] x_imag  時間信号の虚部(結果出力先)
 * @param [in]  X_real  周波数信号の実部(入力)
 * @param [in]  X_imag  周波数信号の虚部(入力)
 * @param [in]  N       標本数
 * @see auto_ift
 * @see idft_ri
 * @see idft_r
 * @see idft_i
 */
void idft(double *restrict x_real, double *restrict x_imag, const double *restrict X_real, const double *restrict X_imag, unsigned int N) {
  if (X_real != NULL && X_imag != NULL) {
    idft_ri(x_real, x_imag, X_real, X_imag, N);
  } else if (X_imag == NULL) {
    idft_r(x_real, x_imag, X_real, N);
  } else if (X_real == NULL) {
    idft_i(x_real, x_imag, X_imag, N);
  }
}


/*!
 * @brief 実部信号と虚部信号がある場合の逆離散フーリエ変換を行う
 * @param [out] x_real  時間信号の実部(結果出力先)
 * @param [out] x_imag  時間信号の虚部(結果出力先)
 * @param [in]  X_real  周波数信号の実部(入力)
 * @param [in]  X_imag  周波数信号の虚部(入力)
 * @param [in]  N       標本数
 * @see idft
 */
static void idft_ri(double *restrict x_real, double *restrict x_imag, const double *restrict X_real, const double *restrict X_imag, unsigned int N) {
  unsigned int n;
  #pragma omp parallel for
  for (n = 0; n < N; n++) {
    unsigned int k;
    for (k = 0; k < N; k++) {
      double w_real = cos(2.0 * M_PI * k * n / N);
      double w_imag = sin(2.0 * M_PI * k * n / N);
      x_real[n] += (w_real * X_real[k] - w_imag * X_imag[k]) / N;  // x(n)の実数部
      x_imag[n] += (w_real * X_imag[k] + w_imag * X_real[k]) / N;  // x(n)の虚数部
    }
  }
}


/*!
 * @brief 実部信号のみがある場合の逆離散フーリエ変換を行う
 * @param [out] x_real  時間信号の実部(結果出力先)
 * @param [out] x_imag  時間信号の虚部(結果出力先)
 * @param [in]  X_real  周波数信号の実部(入力)
 * @param [in]  N       標本数
 * @see idft
 */
static void idft_r(double *restrict x_real, double *restrict x_imag, const double *restrict X_real, unsigned int N) {
  unsigned int n;
  #pragma omp parallel for
  for (n = 0; n < N; n++) {
    unsigned int k;
    for (k = 0; k < N; k++) {
      double w_real = cos(2.0 * M_PI * k * n / N);
      double w_imag = sin(2.0 * M_PI * k * n / N);
      x_real[n] += w_real * X_real[k] / N;  // x(n)の実数部
      x_imag[n] += w_imag * X_real[k] / N;  // x(n)の虚数部
    }
  }
}


/*!
 * @brief 虚部信号のみがある場合の逆離散フーリエ変換を行う
 * @param [out] x_real  時間信号の実部(結果出力先)
 * @param [out] x_imag  時間信号の虚部(結果出力先)
 * @param [in]  X_imag  周波数信号の虚部(入力)
 * @param [in]  N       標本数
 * @see idft
 */
static void idft_i(double *restrict x_real, double *restrict x_imag, const double *restrict X_imag, unsigned int N) {
  unsigned int n;
  #pragma omp parallel for
  for (n = 0; n < N; n++) {
    unsigned int k;
    for (k = 0; k < N; k++) {
      double w_real = cos(2.0 * M_PI * k * n / N);
      double w_imag = sin(2.0 * M_PI * k * n / N);
      x_real[n] += w_imag * X_imag[k] / N;  // x(n)の実数部
      x_imag[n] += w_real * X_imag[k] / N;  // x(n)の虚数部
    }
  }
}


/*!
 * @brief 逆高速離散フーリエ変換を行う
 * @param [out] x_real  時間信号の実部(結果出力先)
 * @param [out] x_imag  時間信号の虚部(結果出力先)
 * @param [in]  X_real  周波数信号の実部(入力)
 * @param [in]  X_imag  周波数信号の虚部(入力)
 * @param [in]  N       標本数
 * @see auto_idft
 */
void ifft(double *restrict x_real, double *restrict x_imag, const double *restrict X_real, const double *restrict X_imag, unsigned int N) {
  unsigned int  k, stage;
  unsigned int *restrict index;
  unsigned int  number_of_stage = my_log2(N);  // FFTの段数

  // Copy
  #pragma omp parallel sections
  {
    #pragma omp parallel section
    dbl_array_copy(x_real, X_real, N);
    #pragma omp parallel section
    dbl_array_copy(x_imag, X_imag, N);
  }

  // バタフライ計算
  for (stage = 1; stage <= number_of_stage; stage++) {
    unsigned int i;
    unsigned int p1 = my_pow2(stage - 1);
    for (i = 0; i < p1; i++) {
      unsigned int j;
      unsigned int p2 = my_pow2(number_of_stage - stage);
      for (j = 0; j < p2; j++) {
        int n = my_pow2(number_of_stage - stage + 1) * i + j;
        int m = my_pow2(number_of_stage - stage) + n;
        int r = my_pow2(stage - 1) * j;
        double a_real = x_real[n];
        double a_imag = x_imag[n];
        double b_real = x_real[m];
        double b_imag = x_imag[m];
        if (stage < number_of_stage) {
          double w_real = cos((2.0 * M_PI * r) / N);
          double w_imag = sin((2.0 * M_PI * r) / N);
          x_real[n] = a_real + b_real;
          x_imag[n] = a_imag + b_imag;
          x_real[m] = (a_real - b_real) * w_real - (a_imag - b_imag) * w_imag;
          x_imag[m] = (a_imag - b_imag) * w_real + (a_real - b_real) * w_imag;
        } else {
          x_real[n] = a_real + b_real;
          x_imag[n] = a_imag + b_imag;
          x_real[m] = a_real - b_real;
          x_imag[m] = a_imag - b_imag;
        }
      }
    }
  }

  // インデックスの並び替えのためのテーブルの作成
  index = CALLOC(unsigned int, N);
  for (stage = 1; stage <= number_of_stage; stage++) {
    unsigned int i;
    for (i = 0; i < my_pow2(stage - 1); i++) {
      index[my_pow2(stage - 1) + i] = index[i] + my_pow2(number_of_stage - stage);
    }
  }

  /* インデックスの並び替え */
  #pragma omp parallel for
  for (k = 0; k < N; k++) {
    if (index[k] > k) {
      SWAP(double, &x_real[index[k]], &x_real[k]);
      SWAP(double, &x_imag[index[k]], &x_imag[k]);
    }
  }

  // 計算結果をNで割る
  #pragma omp parallel for
  for (k = 0; k < N; k++) {
    x_real[k] /= N;
    x_imag[k] /= N;
  }

  FREE(index);
}


/*!
 * @brief 2の対数を計算する
 * @param [in] n  log2の引数(ただし、正数)
 * @return log2(n) (ただし、正数)
 */
inline static unsigned int my_log2(unsigned int n) {
  unsigned int log2_val = 0;
  while (n > 1) {
    n >>= 1;
    log2_val++;
  }
  return log2_val;
}


/*!
 * @brief 2の累乗を計算する
 * @param [in] n  指数部(ただし、正数)
 * @return 2 ** n
 */
inline static unsigned int my_pow2(unsigned int n) {
  return n == 0 ? 1 : (2 << (n - 1));
}


/*!
 * @brief double型の配列のコピーを行う
 *
 * コピー元配列がNULLの場合、出力先配列にゼロをセットする
 * @param [out] dst  コピー先double型配列
 * @param [in]  src  コピー元double型配列
 * @param [in]  n    コピーする要素数
 */
static void dbl_array_copy(double *restrict dst, const double *restrict src, unsigned int n) {
  if (src != NULL) {
    memcpy(dst, src, sizeof(double) * n);
  } else {
    unsigned int i;
    for (i = 0; i < n; i++) {
      dst[i] = 0;
    }
  }
}
