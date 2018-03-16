/*
 * affine.c
 *
 *  Created on: 2017-9-23
 *      Author: Eddie Zhou
 */
#include "affine.h"

#pragma DATA_SECTION(A, ".ddr_mem");
#pragma DATA_ALIGN(A,8);
float A[36];
#pragma DATA_SECTION(B, ".ddr_mem");
#pragma DATA_ALIGN(B,8);
float B[12];

int getAffMatrix(uint16_t *X, uint16_t *Y, uint16_t *newX, uint16_t *newY, float *MM, int count, uint8_t fullAffine)
{
	int i, status;
	//初始化各个矩阵
	memset(A , 0, 36*sizeof(float));
	memset(B , 0, 12*sizeof(float));
	memset(MM, 0, 12*sizeof(float));

	if(fullAffine)
	{
		for(i = 0; i < count; i++)
		{
			A[0] += (X[i] >> 4) * (X[i] >> 4);
			A[1] += (Y[i] >> 4) * (X[i] >> 4);
			A[2] += X[i] >> 4;

			A[6] += (X[i] >> 4) * (Y[i] >> 4);
			A[7] += (Y[i] >> 4) * (Y[i] >> 4);
			A[8] += Y[i] >> 4;

			A[12] += X[i] >> 4;
			A[13] += Y[i] >> 4;
			A[14] += 1.0;

			B[0] += (X[i] >> 4) * (newX[i] / 16.0);
			B[1] += (Y[i] >> 4) * (newX[i] / 16.0);
			B[2] += newX[i] / 16.0;
			B[3] += (X[i] >> 4) * (newY[i] / 16.0);
			B[4] += (Y[i] >> 4) * (newY[i] / 16.0);
			B[5] += newY[i] / 16.0;
		}

		A[21] = A[0];
		A[22] = A[1];
		A[23] = A[2];
		A[27] = A[6];
		A[28] = A[7];
		A[29] = A[8];
		A[33] = A[12];
		A[34] = A[13];
		A[35] = A[14];

		// A B 计算完毕

		// LU分解法求矩阵方程的解 A * MM = B
		status = solveAffMatrix(A, B, MM, 6, METHOD);
	}
	else
	{
		for(i = 0; i < count; i++)
		{
			A[0] += (X[i] >> 4) * (X[i] >> 4) + (Y[i] >> 4) * (Y[i] >> 4);
			A[2] += X[i] >> 4;
			A[3] += Y[i] >> 4;
			A[10]+= 1.0;

			B[0] += (X[i] >> 4) * (newX[i] / 16.0) + (Y[i] >> 4) * (newY[i] / 16.0);
			B[1] += (X[i] >> 4) * (newY[i] / 16.0) - (Y[i] >> 4) * (newX[i] / 16.0);
			B[2] += newX[i] / 16.0;
			B[3] += newY[i] / 16.0;
		}

		A[5]  = A[0];
		A[6]  =-A[3];
		A[7]  = A[2];
		A[8]  = A[2];
		A[9]  =-A[3];
		A[12] = A[3];
		A[13] = A[2];
		A[15] = A[10];

		// A B 计算完毕

		// LU分解法求矩阵方程的解 A * MM = B
//		DSPF_sp_lud(4, A, L, U, P);             // 先做LU分解
//		DSPF_sp_lud_inverse(4, P, L, U, inv_A); // 利用L U P矩阵求矩阵的逆 求得A^-1
//		DSPF_sp_mat_mul(inv_A, 4, 4, B, 2, MM); // MM = A^-1*B
		status = solveAffMatrix(A, B, MM, 4, METHOD);
//		MM[8] = MM[0];
//		MM[10]= MM[6];
//		MM[6] = MM[2];
//		MM[2] =-MM[2];
		MM[4] = MM[0];
		MM[5] = MM[3];
		MM[3] = MM[1];
		MM[1] =-MM[1];
	}
	return status;
}

int RANSAC(uint16_t *X, uint16_t *Y, uint16_t *newX, uint16_t *newY, int count, int fullAffine)
{
	if(count < 3) return 0;

	int i = 0, j = 0, k = 0, k1 = 0, flag = 1, goodCount = 0;
	uint16_t idx[RANSAC_SIZE0] ,ax[RANSAC_SIZE0], ay[RANSAC_SIZE0], bx[RANSAC_SIZE0], by[RANSAC_SIZE0], goodIdx[count];
	float MM[12], bias;
	// 判断是否为内点的阈值
	float th = MAX((DSP_maxval(newX, count) - DSP_minval(newX, count)), (DSP_maxval(newY, count) - DSP_minval(newY, count))) / 16.0;
	for( k = 0; k < RANSAC_MAX_ITERS; k++)
	{

		memset(ax, 0, sizeof(ax));
		memset(bx, 0, sizeof(bx));
		memset(ay, 0, sizeof(ay));
		memset(by, 0, sizeof(by));

		for(i = 0; i < RANSAC_SIZE0; i++) // 取点开始
		{
			for(k1 = 0; k1 < RANSAC_MAX_ITERS; k1++) // 取每个点的迭代次数
			{
				flag = 1;
				idx[i] = rand() % count;

				// check that the points are not very close one each other
				for(j = 0; flag && j < i; j++)
				{
					flag = idx[j] == idx[i] ? 0 : 1;
					flag = (fabs((X[idx[i]] - X[idx[j]]) / 16.0) + fabs((Y[idx[i]] - Y[idx[j]]) / 16.0)) < FLT_EPSILON ? 0 : 1;
					flag = (fabs((newX[idx[i]] - newX[idx[j]]) / 16.0) + fabs((newY[idx[i]] - newY[idx[j]]) / 16.0)) < FLT_EPSILON ? 0 : 1;
				}

				if(!flag) continue; // 找点不符合条件，这个点重新找

				if(i + 1 == RANSAC_SIZE0) // 找够点了
				{
					//printf("\n idx: %d %d %d\n",idx[0], idx[1], idx[2]);
					// 要求三点不共线
					ax[0] = X[idx[0]];
					ax[1] = X[idx[1]];
					ax[2] = X[idx[2]];
					ay[0] = Y[idx[0]];
					ay[1] = Y[idx[1]];
					ay[2] = Y[idx[2]];

					bx[0] = newX[idx[0]];
					bx[1] = newX[idx[1]];
					bx[2] = newX[idx[2]];
					by[0] = newY[idx[0]];
					by[1] = newY[idx[1]];
					by[2] = newY[idx[2]];

					float dax1 = (ax[1] - ax[0]) >> 4, day1 = (ay[1] - ay[0]) >> 4;
					float dax2 = (ax[2] - ax[0]) >> 4, day2 = (ay[2] - ay[0]) >> 4;
					float dbx1 = (bx[1] - bx[0]) / 16.0, dby1 = (by[1] - by[0]) / 16.0;
					float dbx2 = (bx[2] - bx[0]) / 16.0, dby2 = (by[2] - by[0]) / 16.0;
					const float eps = 0.01;
					if (fabs(dax1 * day2 - day1 * dax2) < eps * sqrtsp(dax1 * dax1 + day1 * day1) * sqrtsp(dax2 * dax2 + day2 * day2)||
						fabs(dbx1 * dby2 - dby1 * dbx2) < eps * sqrtsp(dbx1 * dbx1 + dby1 * dby1) * sqrtsp(dbx2 * dbx2 + dby2 * dby2))
						continue; // 第三个点不满足共线，重新找第三个点 k
				}
				break; // 当前找的点满足<不相同且不共线>条件（跳出k1循环，i++）
			}

			if(k1 >= RANSAC_MAX_ITERS) break; // k1 到达最大迭代次数还是未能找到符合要求的第i个点，从第一个点重新开始（跳出i循环）
		}
		if(i < RANSAC_SIZE0) continue; // k1 没有找齐三个点，直接重新开始找第一个点（k+1 i=0）

		getAffMatrix(ax, ay, bx, by, MM, 3, fullAffine); // 找齐三个点了

		for(i = 0, goodCount = 0; i < count; i++)
		{
			bias = fabs(MM[0] * (X[i] >> 4) + MM[1] * (Y[i] >> 4) + MM[2] - (newX[i] / 16.0))
					+ fabs(MM[3] * (X[i] >> 4) + MM[4] * (Y[i] >> 4) + MM[5] - (newY[i] / 16.0));
			if(bias < th * RANSAC_THRESHOLD)
			{
				goodIdx[goodCount++] = i;
			}
		}

		if(goodCount >= count * RANSAC_GOOD_RATIO) break; // 这三个点满足RANSAC要求
	}

	 if(k >= RANSAC_MAX_ITERS)
		 return 0; // RANSAC失败

	// RANSAC成功,把RANSAC筛出的点取出
	for(i = 0; i < goodCount; i++)
	{
		j = goodIdx[i];
		X[i] = X[j];
		Y[i] = Y[j];
		newX[i] = newX[j];
		newY[i] = newY[j];
	}

	return goodCount;
}


int estimateTransformation(uint8_t count, uint16_t *X, uint16_t *Y, uint16_t *newX, uint16_t *newY, float *MM, uint8_t fullAffine)
{
	// 不够三对匹配点，直接估计失败

	// 能进RANSAC一定有三对以上的匹配点
	int goodCount = RANSAC(X, Y, newX, newY, count, fullAffine);
	printf("%d features\n",goodCount);
	// RANSAC 失败时放射参数MM没有被改变
	if(!goodCount) return 1;

	// 用筛选后的匹配点做仿射估计
	int status = getAffMatrix(X, Y, newX, newY, MM, goodCount, fullAffine);
	if(status == -1) return 1;
	return 0;
}

#pragma DATA_SECTION(Z, ".ddr_mem");			// 平滑
#pragma DATA_ALIGN(Z, 8);
VLIB_F32 Z[4];
#pragma DATA_SECTION(Residual, ".ddr_mem");
#pragma DATA_ALIGN(Residual, 8);
VLIB_F32 Residual[4];

tracePara traceSmooth(float *MM, tracePara *p, VLIB_kalmanFilter_4x8_F32 *KF)
{
	float dx, dy, da;
	VLIB_F32 Z[4],Residual[4];
	dx = MM[2];
	dy = MM[5];
	da = atan2sp(MM[3],MM[0]);
	printf("\tT:\t%f\t%f\t%f\n", dx, dy, da);
	VLIB_kalmanFilter_4x8_Predict_F32(KF);
	Z[0] = p -> x += dx;
	Z[1] = p -> y += dy;
	Z[2] = p -> a += da;
	Z[3] = 0; // 其实只有三个measurement，但是接口是4个
	// 卡尔曼滤波
	VLIB_kalmanFilter_4x8_Correct_F32(KF, Z, Residual);
	// Residual(k) = Z(k) - H * X(k|k-1)(非state的输出结果而是A * X(k-1))

//	float diff_x = - dx - KF -> state[0] + Z[0];
//	float diff_y = - dy - KF -> state[1] + Z[1];
//	float diff_a = - da - KF -> state[2] + Z[2];
	float diff_x = - KF -> state[0] + Z[0];
	float diff_y = - KF -> state[1] + Z[1];
	float diff_a = - KF -> state[2] + Z[2];

	printf("\tZ:\t%f\t%f\t%f\n",Z[0],Z[1],Z[2]);
	printf("\tstate:\t%f\t%f\t%f\n",KF -> state[0],KF -> state[1],KF -> state[2]);
//	p -> x += dx;
//	p -> y += dy;
//	p -> a += da;

	tracePara T = {diff_x, diff_y, diff_a};
	printf("\t\tdx = %f\tdy=%f\tda=%f\n",diff_x,diff_y,diff_a);

	return T;
}

void affinePara(AffineWarpParamq_t *para, tracePara *p)
{
	float co,si,dx,dy;
	co = cossp(p -> a);
	si = sinsp(p -> a);
	dx = p -> x;
	dy = p -> y;
	float shift = powsp(2.0,16.0);
	para -> xstep_c = co * shift;
	para -> xstep_r =-si * shift;
	para -> xshift  = dx * shift;
	para -> ystep_c = si * shift;
	para -> ystep_r = co * shift;
	para -> yshift  = dy * shift;
}


