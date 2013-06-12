/*!
 * @brief フーリエ変換を行い、結果をgnuplotに表示する
 *
 * @author  koturn 0;
 * @date    2013 06/01
 * @file    main.c
 * @version 0.5
 */
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "compatibility.h"
#include "fourier.h"
#include "gplot.h"
#include "wav.h"


#define CALLOC(type, n)  ((type *)calloc((n), sizeof(type)))
#define FREE(ptr_p)      (free(*(ptr_p)), *(ptr_p) = NULL)
#define SQ(x)            ((x) * (x))

static unsigned int str2int(const char *restrict nstr);
static void print_result(const POINT *restrict points, unsigned int N);


/*!
 * @brief プログラムのエントリポイント
 *
 * 使い方<br>
 *   $ ./main [wavefile-name] N
 * @param [in] argc コマンドライン引数の数
 * @param [in] argv コマンドライン引数
 * @return 終了コード
 */
int main(int argc, char *argv[]) {
  unsigned int           N;
  unsigned int           k;
  double                 frq_rate;
  double       *restrict f_real;
  double       *restrict f_imag;
  WAVE_DATA    *restrict wd;
  POINT        *restrict points;
  RANGE                  range = {{0.0, 0.0}, {0.0, 0.0}};

  if (argc != 3) {
    fprintf(stderr, "bad argument!\n");
    return EXIT_FAILURE;
  }
  wd = read_wave_file(argv[1]);
  N  = str2int(argv[2]);
  print_wave_header(wd);

  f_real = CALLOC(double, N);
  f_imag = CALLOC(double, N);
  points = CALLOC(POINT, N);
  if (f_real == NULL || f_imag == NULL || points == NULL) {
    fputs("メモリ確保エラーです\n", stderr);
    return EXIT_FAILURE;
  }
  auto_dft(f_real, f_imag, wd->s_data1, NULL, N);

  frq_rate = (double)wd->samples_per_sec / N;
  for (k = 0; k < N; k++) {
    points[k].x = frq_rate * k;
    points[k].y = sqrt(SQ(f_real[k]) + SQ(f_imag[k]));
    if (range.max.y < points[k].y) range.max.y = points[k].y;
  }
  range.max.x = (double)wd->samples_per_sec / 2;
  range.max.y = (range.max.y == 0.0 ? 1.0 : range.max.y * 1.2);

  // print_result(points, N);
  if (show_gplot(points, N, range) == -1) {
    fputs("gnuplotを起動できません\n", stderr);
    fputs("gnuplotにPATHが通っているか、確認してください\n", stderr);
  }

  FREE(&f_real);
  FREE(&f_imag);
  FREE(&points);
  free_wave_data(wd);
  return EXIT_SUCCESS;
}


/*!
 * @brief 引数の文字列を数値に変換する。
 * @param [in] nstr  数値に変換する文字列
 * @return 変換した数値
 */
static unsigned int str2int(const char *restrict nstr) {
  char *check;
  int   num = strtol(nstr, &check, 10);  // char * -> long
  if (*check != '\0') {
    fputs("文字列に数値以外がありました\n", stderr);
    exit(EXIT_FAILURE);
  }
  if (num <= 0) {
    fputs("0以下の値を指定しないでください\n", stderr);
    exit(EXIT_FAILURE);
  } else if (num == INT_MAX) {
    fputs("値が大きすぎます\n", stderr);
    exit(EXIT_FAILURE);
  }
  return (unsigned int)num;
}


/*!
 * @brief 周波数毎のパワースペクトルを表示する
 * @param [in] points  周波数軸上の座標の配列
 * @param [in] N       標本数
 */
static void print_result(const POINT *restrict points, unsigned int N) {
  unsigned int k;
  printf("%12s = %-12s\n", "Frequency", "Power");
  for (k = 0; k < N; k++) {
    printf("%10.2fHz = %10.2f\n", points[k].x, points[k].y);
  }
}
