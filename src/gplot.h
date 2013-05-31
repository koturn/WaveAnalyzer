/*!
 * @brief gplot.cのヘッダファイル
 *
 * @author  koturn 0;
 * @date    2013 05/30
 * @file    gplot.h
 * @version 0.1
 */
#ifndef GPLOT_H
#define GPLOT_H

#include "compatibility.h"


typedef struct {
  double x;
  double y;
} POINT;

typedef struct {
  POINT min;
  POINT max;
} RANGE;

int show_gplot(const POINT *restrict points, unsigned int N, RANGE range);


#endif
