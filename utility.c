/*
 * utility.c
 *
 *  Created on: 2017-8-21
 *      Author: Eddie Zhou
 */
#include "utility.h"

void pyramidInit(uint8_t **pyr, uint8_t *pyrbuf, uint8_t *im, int width, int height)
{
    pyr[0] = im;
    pyr[1] = pyrbuf;
    pyr[2] = pyrbuf + width / 2 * height / 2;
    pyr[3] = pyrbuf + width / 2 * height / 2 + width / 4 * height / 4;
    pyr[4] = pyrbuf + width / 2 * height / 2 + width / 4 * height / 4 + width / 8 * height / 8;
    pyr[5] = pyrbuf + width / 2 * height / 2 + width / 4 * height / 4 + width / 8 * height / 8 + width / 16 * height / 16;
    pyr[6] = pyrbuf + width / 2 * height / 2 + width / 4 * height / 4 + width / 8 * height / 8 + width / 16 * height / 16 + width / 32 * height / 32;
}

void imagePyramid8_6Level(uint8_t *im, uint8_t *pyrbuf, int width, int height)
{
	VLIB_imagePyramid8(im, width, height, pyrbuf); // 三层
	VLIB_imagePyramid8(pyrbuf + width / 2 * height / 2 + width / 4 * height / 4, width / 8, height / 8, pyrbuf + width / 2 * height / 2 + width / 4 * height / 4 + width / 8 * height / 8);

}

void kalmanInit_4x8_F32(VLIB_kalmanFilter_4x8_F32 *KF)
{
	// X(k)=A X(k-1) + B U(k) + W(k) 状态方程
	// Z(k)=H X(k) + V(k) 测量方程
	memset(KF, 0, sizeof(VLIB_kalmanFilter_4x8_F32));
	int i, index;
	for( i = 0; i < sD_4x8; i++ )
	{
		index = i * sD_4x8 + i;
		KF->errorCov[index] = INIT_ERROR_COV; // 最优角度估计的协方差P(k|k)
		KF->transition[index] = 1; // transition 对角元素填1
		if( i < mD_4x8 )
		{
			KF->measurement[index] = 1; // measurement 是 矩阵H[4x8] 对角矩阵
		}
	}

//	KF->transition[4]  = 1;
//	KF->transition[13] = 1;
//	KF->transition[22] = 1; // transition 是 矩阵A[8x8]

	KF->measurementNoiseCov[0]  = NOICE_COV_M_T; // 测量噪声 R
	KF->measurementNoiseCov[5]  = NOICE_COV_M_T;
	KF->measurementNoiseCov[10] = NOICE_COV_M_R;
	KF->measurementNoiseCov[15] = 1;

	KF->processNoiseCov[0]  = NOICE_COV_S_T; // 过程噪声 Q
	//KF->processNoiseCov[4]  = 0.25;
	KF->processNoiseCov[9]  = NOICE_COV_S_T;
	//KF->processNoiseCov[13] = 0.25;
	KF->processNoiseCov[18] = NOICE_COV_S_R;
	KF->processNoiseCov[27] = 1;
//	KF->processNoiseCov[36] = NOICE_COV_S_TV;
//	KF->processNoiseCov[45] = NOICE_COV_S_TV;
//	KF->processNoiseCov[54] = NOICE_COV_S_RV;
	KF->processNoiseCov[36] = 1;
	KF->processNoiseCov[45] = 1;
	KF->processNoiseCov[54] = 1;
	KF->processNoiseCov[63] = 1;


	/* set initial state to all 0s */
	KF->state[0] = 0; // X(k) in and out
	KF->state[1] = 0;
	KF->state[2] = 0;
	KF->state[3] = 0;
	KF->state[4] = 0;
	KF->state[5] = 0;
	KF->state[6] = 0;
	KF->state[7] = 0;

	KF->scaleFactor = 1;
}


int readdata(char *fname, uint8_t *im, int height, int width, int dim)
{

	FILE *fp;

	fp = fopen(fname, "rb");
	if (fp == NULL)
	{
		printf("can't open the image.\n ");
		exit(0);
	}
	else
	{
		fread(im, sizeof(uint8_t), height * width * dim, fp);
	}
	// printf("Read successfully.\n");
	fclose(fp);
	return 0;
}

int writedata(char *fname, uint8_t *im, int height, int width, int dim)
{
	FILE *fp;

	fp = fopen(fname, "wb");
	if (fp == NULL)
	{
		printf("can't create image file.\n");
		exit(0);
	}
	else
	{
		fwrite(im, sizeof(uint8_t), height * width * dim, fp);
	}
	// printf("Write successfully.\n");
	fclose(fp);
	return 0;
}

int goodTrackingPoints(int16_t *score, int width, int height, uint8_t *feat, uint16_t *X, uint16_t *Y, int thresh)
{
	int i,j;
	int nFeatures = 0, flag = 1;
	memset(feat, 0, height * width * sizeof(uint8_t));
	VLIB_nonMaxSuppress_7x7_S16(score, width, height, thresh, feat + 3 * width + 3); // 和梯度一样，输出也需要指向非边缘（边缘为两个像素）的第一个像素的位置
	// border：梯度1 + Harris角点3 + 非极大值抑制3 = 7


//	printf("\n");
	for (i = 7; flag && i < (height - 7); i++)
	{
		for(j = 7; flag && j < (width - 7); j++)
		{
			if(feat[i * width + j] == 255)
			{

				X[nFeatures] = j << 4;
				Y[nFeatures] = i << 4;
//				// 左移四位又右移三位，为了输入到光流金字塔的最低级
//				newX[i] = j << 1;
//				newY[i] = i << 1;
				nFeatures++;
			}
			flag = nFeatures >= MAX_POINTS ? 0 : 1;
		}
	}
	//printf("\n");
	if(nFeatures%2 != 0)
	{
		nFeatures--;   // 特征点数量需要是2的倍数
	}
	return nFeatures;
}

uint8_t goodMatchSelect(uint16_t *status, int nFeatures, uint16_t *X, uint16_t *Y, uint16_t *newX, uint16_t *newY, int thr)
{
	int i, goodPoints = 0;
	for(i = 0; i < nFeatures; i++)
	{
		if(status[i] < thr)
		{
			X[goodPoints] = X[i];
			Y[goodPoints] = Y[i];
			newX[goodPoints] = newX[i];
			newY[goodPoints] = newY[i];
			goodPoints++;
		}
	}
	return goodPoints;
}


#pragma DATA_SECTION(inv_A, ".ddr_mem");		// 逆
#pragma DATA_ALIGN(inv_A,8);
float inv_A[36];

#pragma DATA_SECTION(L, ".ddr_mem");      		// 分解LU
#pragma DATA_ALIGN(L,8);
float L[36];
#pragma DATA_SECTION(U, ".ddr_mem");
#pragma DATA_ALIGN(U,8);
float U[36];
#pragma DATA_SECTION(P, ".ddr_mem");
#pragma DATA_ALIGN(P,8);
unsigned short P[36];

#pragma DATA_SECTION(Q,".ddr_mem")				// 分解QR
#pragma DATA_ALIGN(Q,8);
float Q[36];
#pragma DATA_SECTION(R,".ddr_mem")
#pragma DATA_ALIGN(R,8);
float R[36];
#pragma DATA_SECTION(u,".ddr_mem")
#pragma DATA_ALIGN(u,8);
float u[6];
#pragma DATA_SECTION(y,".ddr_mem")
#pragma DATA_ALIGN(y,8);
float y[12];

#pragma DATA_SECTION(V,".ddr_mem")				// 分解SVD
#pragma DATA_ALIGN(V,8)
float V[36];
#pragma DATA_SECTION(U1,".ddr_mem")
#pragma DATA_ALIGN(U1,8)
float U1[36];
#pragma DATA_SECTION(transU,".ddr_mem")
#pragma DATA_ALIGN(transU,8)
float transU[36];
#pragma DATA_SECTION(diag,".ddr_mem")
#pragma DATA_ALIGN(diag,8)
float diag[6];
#pragma DATA_SECTION(rediag,".ddr_mem")
#pragma DATA_ALIGN(rediag,8)
float rediag[6];
#pragma DATA_SECTION(D,".ddr_mem")
#pragma DATA_ALIGN(D,8)
float D[36];
#pragma DATA_SECTION(VD,".ddr_mem")
#pragma DATA_ALIGN(VD,8)
float VD[36];
#pragma DATA_SECTION(superdiag,".ddr_mem")
#pragma DATA_ALIGN(superdiag,8)
float superdiag[6];

int solveAffMatrix(float *A, float *B, float *MM, int orderA, int method)
{
	if(method == SOLVE_QR)
	{
		memset(Q, 0, 36 * sizeof(float));
		memset(R, 0, 36 * sizeof(float));
		memset(u, 0, 6 * sizeof(float));
		memset(y, 0, 12 * sizeof(float));

		int status = DSPF_sp_qrd(orderA, orderA, A, Q, R, u);
		if(status!=-1)
		{
			status = DSPF_sp_qrd_solver(orderA, orderA, Q, R, B, y, MM);
		}

		return status;
	}
	else if(method == SOLVE_LU)
	{
		memset(inv_A, 0, 36 * sizeof(float));
		memset(L, 0, 36 * sizeof(float));
		memset(U, 0, 36 * sizeof(float));
		memset(P, 0, 36 * sizeof(float));

		int i = orderA - 1;
		for(; i>=0; i--)
		{
			B[2 * i] = B[i];
			B[2 * i +1] = 0;
		}

		DSPF_sp_lud(orderA, A, L, U, P);             // 先做LU分解
		DSPF_sp_lud_inverse(orderA, P, L, U, inv_A); // 利用L U P矩阵求矩阵的逆 求得A^-1
		DSPF_sp_mat_mul(inv_A, orderA, orderA, B, 2, MM); // MM = A^-1*B

		for(i = 0; i < orderA; i++)
		{
			MM[i] = MM[2 * i];
		}

		return 0;
	}
	else if(method == SOLVE_SVD)
	{
		int i;
		memset(inv_A, 0, 36 * sizeof(float));
		memset(V, 0, 36 * sizeof(float));
		memset(U, 0, 36 * sizeof(float));
		memset(transU, 0, 36 * sizeof(float));
		memset(D, 0, 36 * sizeof(float));
		memset(VD, 0, 36 * sizeof(float));
		memset(U1, 0, 36 * sizeof(float));
		memset(diag, 0, 6 * sizeof(float));
		memset(superdiag, 0, 6 * sizeof(float));

		int status = DSPF_sp_svd(orderA, orderA, A, U, V, U1, diag, superdiag); // A = U * D * V'
		if(status == -1)
		{
			printf("decomposition fail!\n");
		}

		// singular values less than a tolerance are treated as zero
		float tolerance = FLT_EPSILON * orderA * DSPF_sp_maxval(diag, orderA);
		DSPF_sp_vecrecip(diag, rediag, 4); // D^-1 = 1 / diag
		memset(D, 0, 36 * sizeof(float));

		for(i = 0; i < orderA; i++)
		{
//			D[i * 4 + i] = diag[i] > tolerance ? rediag[i] : 0;
			D[i * 4 + i] = rediag[i];
		}

		// inv_A = V * D^-1 * U'
		DSPF_sp_mat_mul(V, 4, 4, D, 4, VD); // V * D^-1
		DSPF_sp_mat_trans(U, 4, 4, transU); // U'
		DSPF_sp_mat_mul(VD, 4, 4, transU, 4, inv_A);

		for(i = 3; i>=0; i--)
		{
			B[2 * i] = B[i];
			B[2 * i +1] = 0;
		}

		DSPF_sp_mat_mul(inv_A, 4, 4, B, 2, MM); // MM = A^-1*B

		for(i = 0; i < orderA; i++)
		{
			MM[i] = MM[2 * i];
		}

		return 0;
	}
	else
		return -1;
}
