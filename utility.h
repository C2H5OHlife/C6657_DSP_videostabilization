/*
 * utility.h
 *
 *  Created on: 2017-8-21
 *      Author: Eddie Zhou
 */

#ifndef UTILITY_H_
#define UTILITY_H_

#include "para.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ti/vlib/vlib.h>
#include <ti/mathlib/mathlib.h>
#include <ti/dsplib/dsplib.h>
#include <math.h>
#include <float.h>

#define SW_BREAKPOINT   asm(" SWBP 0 ");

#define 	MAX(a,b)	((a) > (b) ? (a) : (b))

/**
 * 金字塔指针初始化
 */
void pyramidInit(uint8_t **pyr, uint8_t *pyrbuf, uint8_t *im, int width, int height);

/**
 * 金字塔生成（六层）
 */
void imagePyramid8_6Level(uint8_t *im, uint8_t *pyrbuf, int width, int height);

/**
 * 卡尔曼滤波器初始化
 */
void kalmanInit_4x8_F32(VLIB_kalmanFilter_4x8_F32 *KF);

/**
 * 数据读写
 */
int readdata(char *fname, uint8_t *im, int height, int width, int dim);
int writedata(char*fname, uint8_t *im, int height, int width, int dim);

/**
 * 特征点筛选
 */
int goodTrackingPoints(int16_t *score, int width, int height, uint8_t *feat, uint16_t *X, uint16_t *Y, int thresh);

/**
 * 匹配点筛选
 */
uint8_t goodMatchSelect(uint16_t *status, int nFeatures, uint16_t *X, uint16_t *Y, uint16_t *newX, uint16_t *newY, int thr);

/**
 * 求解方程组
 */
int solveAffMatrix(float *A, float *B, float *MM, int orderA, int method);

#endif /* UTILITY_H_ */
