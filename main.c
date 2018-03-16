/*
 * main.c
 */
#include "utility.h"
#include "affine.h"
#include <ti/csl/csl_cacheAux.h>
#include <c6x.h>
#include <time.h>

// �ڴ����
#pragma DATA_SECTION(imPrev, ".ddr_mem");       // ͼ��
uint8_t imPrev[HEIGHT * WIDTH * DIM];
#pragma DATA_SECTION(imCur, ".ddr_mem");
uint8_t imCur[HEIGHT * WIDTH * DIM];
#pragma DATA_SECTION(imTrans, ".ddr_mem");
uint8_t imTrans[HEIGHT * WIDTH * DIM];

#pragma DATA_SECTION(gradx, ".ddr_mem");  		// �ݶ�
int16_t gradx[HEIGHT * WIDTH];
#pragma DATA_SECTION(grady, ".ddr_mem");
int16_t grady[HEIGHT * WIDTH];

#pragma DATA_SECTION(score, ".ddr_mem");  		// �÷�
int16_t score[HEIGHT * WIDTH];
#pragma DATA_SECTION(feat, ".ddr_mem");   		// �ǵ�
uint8_t feat[HEIGHT * WIDTH];
#pragma DATA_SECTION(status, ".ddr_mem"); 		// ���
uint16_t status[MAX_POINTS];

#pragma DATA_SECTION(scratch, ".sharebuf");		// ����
#pragma DATA_ALIGN(scratch, 2);
uint8_t scratch[893];
#pragma DATA_SECTION(buffer, ".ddr_mem");
#pragma DATA_ALIGN(buffer, 4);
uint8_t buffer[96 * WIDTH];

#pragma DATA_SECTION(X, ".sharebuf");      		// ����
#pragma DATA_ALIGN(X, 4);
uint16_t X[MAX_POINTS];
#pragma DATA_SECTION(Y, ".sharebuf");
#pragma DATA_ALIGN(Y, 4);
uint16_t Y[MAX_POINTS];
#pragma DATA_SECTION(newX, ".sharebuf");
#pragma DATA_ALIGN(newX, 4);
uint16_t newX[MAX_POINTS];
#pragma DATA_SECTION(newY, ".sharebuf");
#pragma DATA_ALIGN(newY, 4);
uint16_t newY[MAX_POINTS];
#pragma DATA_SECTION(pointIndex, ".sharebuf");
#pragma DATA_ALIGN(pointIndex, 4);
uint8_t pointIndex[MAX_POINTS];

#pragma DATA_SECTION(M, ".ddr_mem");			// ����
float M[12];


#pragma DATA_SECTION(oldpyrbuf,".ddr_mem");		// ���������
uint8_t oldpyrbuf[HEIGHT * WIDTH * 21 / 64];
#pragma DATA_SECTION(newpyrbuf,".ddr_mem");
uint8_t newpyrbuf[HEIGHT * WIDTH * 21 / 64];
#pragma DATA_SECTION(pyrgradx, ".ddr_mem");		// �������ݶ�
int16_t pyrgradx[HEIGHT * WIDTH];
#pragma DATA_SECTION(pyrgrady, ".ddr_mem");
int16_t pyrgrady[HEIGHT * WIDTH];
#pragma DATA_SECTION(pyramidX, ".ddr_mem");		// ����������
uint16_t pyramidX[200];
#pragma DATA_SECTION(pyramidY, ".ddr_mem");
uint16_t pyramidY[200];

#pragma DATA_SECTION(_KF, ".ddr_mem");			// ƽ��
#pragma DATA_ALIGN(_KF, 8);
VLIB_kalmanFilter_4x8_F32 _KF;

int main(void) {
	clock_t start, finish;
	double time = 0.0;

	// ͼƬԴ��ַ
	const char img_name[512]="E:\\Dataset\\Stabilization\\BlurFace_bin\\%03d.bin";
	const char dst_name[512]="E:\\Dataset\\Stabilization\\BlurFace_result\\%d.bin";
	char img_path[512];
	char dst_path[512];

	// ������õ��Ĳ���
	int i,j;
	int nFeatures = 0;
	uint8_t goodPoints = 0;  // �õ�ƥ��������ĸ���
	tracePara currentTrace = {0.0, 0.0, 0.0};

	// ����任�ṹ��
	ImgBlock8_t outBlock = {imTrans, HEIGHT, WIDTH, WIDTH*sizeof(uint8_t)};  // ����任�ò��� ע�����任�����ӿڵ��߼�
	ImgBlock8_t inBlock  = {imCur, HEIGHT, WIDTH, WIDTH*sizeof(uint8_t)};
	AffineWarpParamq_t affPara;


	// �ɡ��� ������ָ�� 3����ԭͼ�� 4��ָ��ֱ�ָ��ԭͼlevel0��level1/2/3
    uint8_t *oldpyr[7];
    uint8_t *newpyr[7];
    pyramidInit(oldpyr, oldpyrbuf, imPrev, WIDTH, HEIGHT);
    pyramidInit(newpyr, newpyrbuf,  imCur, WIDTH, HEIGHT);

    // �������˲���ʼ��
    VLIB_kalmanFilter_4x8_F32 *KF = &_KF;
    kalmanInit_4x8_F32(KF);

    // --��ʼ
	int frame = START_FRAME;
	sprintf(img_path, img_name, frame);
	readdata(img_path, imPrev, HEIGHT, WIDTH, DIM);

	// ����ǰһ֡�Ľ�����(3��)
    VLIB_imagePyramid8(imPrev, WIDTH, HEIGHT, oldpyrbuf);

	// �����һ֡imTrans = imPrev
	sprintf(dst_path, dst_name, frame);
	writedata(dst_path, imPrev, HEIGHT, WIDTH, DIM);

	int16_t k = HARRIS_K * powsp(2.0, 15.0);  // 0.04 -> UQ 1.15 (1310) С����kҪ��һ��������������������ʾ

	TSCH = 0;
	TSCL = 0;
	start = _itoll(TSCH,TSCL);
	printf("\n");
	for (frame = START_FRAME + 1; frame <= END_FRAME; frame++)
	{
		printf("%d.", frame);
		// Step.1 ������ƥ��
		// --ǰһ֡�ݶȼ���
		memset(gradx, 0, sizeof(int16_t) * HEIGHT * WIDTH);
		memset(grady, 0, sizeof(int16_t) * HEIGHT * WIDTH);
		VLIB_xyGradients(imPrev, gradx + WIDTH + 1, grady + WIDTH + 1, WIDTH, HEIGHT - 1); // gradָ��ƫ�Ʒ�ֹ��������������ƫһ������

		// --��ǰһ֡�����ǵ���
		memset(score, 0, HEIGHT*WIDTH*sizeof(int16_t));
		VLIB_harrisScore_7x7(gradx, grady, WIDTH, HEIGHT, score, k, buffer);  // ������ͼscore


		// --ɸѡ���ٵ� ��X,Y������4λ��
		nFeatures = goodTrackingPoints(score, WIDTH, HEIGHT, feat, X, Y, DSP_maxval(score, WIDTH * HEIGHT) * QUALITY_LEVEL);
		printf("\t%d corners -> ",nFeatures);

		// --��ȡ��֡
		sprintf(img_path, img_name, frame);
		readdata(img_path, imCur, HEIGHT, WIDTH, DIM);

		VLIB_imagePyramid8(imCur, WIDTH, HEIGHT, newpyrbuf);

		// Since we are starting from the level 3 in pyramid, the initial estimates of newX, newY are same as the
		// features of level 3 image of imPrev

		for(i = 0; i<nFeatures; i++) {
			newX[i] = X[i] >> 3;
			newY[i] = Y[i] >> 3;
		}

		// --��������������������ƥ��
		for (i = 3; i > 0; i--) // ����ͷֱ��ʿ�ʼ����,����������������
		{
			// pyramidX, pyramidY ����level i of imPrev�ϵ�������
			for (j = 0; j < nFeatures; j++)
			{
				pyramidX[j] = X[j] >> i;
				pyramidY[j] = Y[j] >> i;
			}

			// Estimates newX, newY are updated at level i. Input features are pyramidX, pyramidY
			memset(pyrgradx, 0, sizeof(int16_t) * HEIGHT * WIDTH);
			memset(pyrgrady, 0, sizeof(int16_t) * HEIGHT * WIDTH);
			VLIB_xyGradients(oldpyr[i], pyrgradx + (WIDTH >> i) + 1, pyrgrady + (WIDTH >> i) + 1, WIDTH >> i, (HEIGHT >> i) - 1);
			VLIB_trackFeaturesLucasKanade_7x7(oldpyr[i], newpyr[i], pyrgradx, pyrgrady, WIDTH >> i, HEIGHT >> i, nFeatures,
								  pyramidX, pyramidY, newX, newY, NULL, MAX_ITER, 0, scratch);

			// newX, newY refined at level i are scaled to become estimates for next iteration
			// newX��newY������������ʾ��һ���ֱ����ϵĳ�ʼ��������
			for (j = 0; j < nFeatures; j++)
			{
				newX[j] = newX[j] << 1;
				newY[j] = newY[j] << 1;
			}
		}

		// Fine tune newX,newY fourth time with original resolution images
		VLIB_trackFeaturesLucasKanade_7x7(imPrev, imCur, gradx, grady, WIDTH, HEIGHT, nFeatures,
										  X, Y, newX, newY, status, MAX_ITER, 1e-4, scratch);

		// --ƥ��������ɸѡ X, Y, newX, newY are modified. nFeatures -> goodPoints.
		goodPoints = goodMatchSelect(status, nFeatures, X, Y, newX, newY, 10000);
		printf("%d features -> ",goodPoints);

		// Step.2 ȫ���˶�����
		// --���Ʒ���仯
		uint8_t err = estimateTransformation(goodPoints, X, Y, newX, newY, M, 0); // ���һ�����������Ƿ�Ϊȫ���䣨���߶ȣ�
		if(err)
		{
			printf("error when estimate the transformation parameters\n");
		}
		tracePara T = traceSmooth(M, &currentTrace, KF);
		affinePara(&affPara, &T);
		util_affineWarp8q(&inBlock, &affPara, &outBlock);
		// ���imTrans
//		sprintf(dst_path, dst_name, frame);
//		writedata(dst_path, imTrans, HEIGHT, WIDTH, DIM);

		// ��ǰ֡��ǰһ֡
		memcpy(imPrev, imCur, WIDTH * HEIGHT);
		memcpy(oldpyrbuf, newpyrbuf, WIDTH * HEIGHT * 21 / 64);
	}

	finish = _itoll(TSCH,TSCL);
	time = 1000000.0 * (END_FRAME - START_FRAME) / (finish - start);
	time = time * 1000.0;
	printf("%f fps\n",time);
	printf("ָ����������%d\n",(finish - start));
//	printf("ָ����������%d\n",time0);

	SW_BREAKPOINT
	return 0;
}
