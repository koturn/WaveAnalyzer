/*!
 * @brief gnuplotをC言語から使用するラッパー関数を提供する
 * @author  koturn 0;
 * @date    2013 06/01
 * @file    gplot.c
 * @version 0.6
 */
#include <stdio.h>
#include <stdlib.h>
#include "compatibility.h"
#include "gplot.h"


#if defined(WIN16) || defined(_WIN16) || defined(__WIN16) || defined(__WIN16__)   \
  || defined(WIN32) || defined(_WIN32) || defined(__WIN32) || defined(__WIN32__)  \
  || defined(WIN64) || defined(_WIN64) || defined(__WIN64) || defined(__WIN64__)
#  define IS_WINDOWS
#endif

#ifdef _MSC_VER
#  define popen(cmd, mode)  _popen(cmd, mode)
#  define pclose(fp)        _pclose(fp)
#endif


/*!
 * @brief 入力座標をgnuplotで表示する
 * @param [in] points  周波数軸上の座標の配列
 * @param [in] N       標本数
 * @param [in] range   gnuplotの描画範囲
 * @return 成功時に0を、失敗時に-1を返す
 */
int show_gplot(const POINT *restrict points, unsigned int N, RANGE range) {
  unsigned int i;
  FILE *gp;
  if ((gp = popen("gnuplot -persist", "w")) == NULL) {
    return -1;
  }

  fprintf(gp, "set xrange [%f:%f]\n", range.min.x, range.max.x);
  fprintf(gp, "set yrange [%f:%f]\n", range.min.y, range.max.y);
  fprintf(gp, "plot '-' with lines linetype 1 title \"Power Spectrum\"\n");
  for (i = 0; i <= N; i++) {
    fprintf(gp, "%f\t%f\n", points[i].x, points[i].y);  // データの書き込み
  }
  fprintf(gp, "e\n");
  fflush(gp);  // バッファに格納されているデータを吐き出す（必須）

#ifdef IS_WINDOWS
  system("pause");
#else
  system("read");
#endif

  pclose(gp);
  return 0;
}
