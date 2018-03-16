/*
 * para.h
 *
 *  Created on: 2017-8-21
 *      Author: Eddie Zhou
 */

#ifndef PARA_H_
#define PARA_H_

// 视频参数
#define		HEIGHT				480
#define		WIDTH 				640
#define 	DIM 				1
#define		START_FRAME			1
#define		END_FRAME			10

// 角点匹配
#define		MAX_POINTS			50
#define		HARRIS_K			0.08
#define		MAX_ITER			10
#define		QUALITY_LEVEL		0.85

// 运动估计
#define		SOLVE_QR			1
#define		SOLVE_LU			2
#define		SOLVE_SVD			3

#define		METHOD				SOLVE_LU

// RANSAC
#define 	RANSAC_MAX_ITERS 	500
#define 	RANSAC_SIZE0 		3
#define		RANSAC_GOOD_RATIO	0.7
#define		RANSAC_THRESHOLD	0.02

// KALMAN
#define		INIT_ERROR_COV		1		// errorCov
#define		NOICE_COV_M_T		1		// measurementNoiseCov for translation
#define		NOICE_COV_M_R		1		// measurementNoiseCov for rotation
#define 	NOICE_COV_S_T		0.001	// processNoiseCov for translation
#define 	NOICE_COV_S_R		0.001	// processNoiseCov for translation
//#define		NOICE_COV_S_TV		1
//#define		NOICE_COV_S_RV		1

typedef struct trace
{
	float x;
	float y;
	float a;
}tracePara;



#endif /* PARA_H_ */
