/*
 * affine.h
 *
 *  Created on: 2017-9-23
 *      Author: Eddie Zhou
 */

#ifndef AFFINE_H_
#define AFFINE_H_

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <util/util.h>
#include <ti/mathlib/mathlib.h>
#include <ti/vlib/vlib.h>
#include <ti/dsplib/dsplib.h>
#include "para.h"
#include "utility.h"

// 卡尔曼滤波
extern VLIB_F32 Z[4];
extern VLIB_F32 Residual[4];
/**
 * 估计全局运动放射变换
 */
int estimateTransformation(uint8_t count, uint16_t *X, uint16_t *Y, uint16_t *newX, uint16_t *newY, float *MM, uint8_t fullAffine);

tracePara traceSmooth(float *MM, tracePara *p, VLIB_kalmanFilter_4x8_F32 *KF);
void affinePara(AffineWarpParamq_t *para, tracePara *p);

#endif /* AFFINE_H_ */
