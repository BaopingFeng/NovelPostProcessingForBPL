
/***************************************************************************
* see paper "A Novel Post-Processing Method for BPL decoding of Polar Codes"
* submit to IEEE Communication Letter, 2021.1.5
* fengbaoping@buaa.edu.cn
*
* Name: NBPL.cpp
* 
* Function: A Novel post-processing method for BPL decoding
*
***************************************************************************/
#include "mex.h"
#include "matrix.h"
#include "math.h"
#include<iostream>
#include<fstream>
#include<string>
#define INF 64.0
#define maxvalue 64.0
using namespace std;
double boxOperate(double a, double b)
{
    double c = log((1+exp(a+b))/(exp(a)+exp(b)));
    if(c>maxvalue)
        c=maxvalue;
    if(c<-maxvalue)
        c=-maxvalue;
    return c;
}
bool crcDector(int* input, int K, int * poly, int crclen,int * crcResult)
{
	int* crc = new int[crclen];
	for (int i = 0; i < crclen; i++)
	{
		crc[i] = input[K - crclen + i];
		input[K - crclen + i] = 0;
	}

	int * re = new int[crclen + 1];
	int str = 0;
	for (int i = 0; i < K; i++)
	{
		if (input[i] == 1)
		{
			str = i;
			break;
		}
	}
	for (int i = 0; i < crclen + 1; i++)
	{
		re[i] = poly[i] ^ input[i + str];
	}
	for (int i = str + crclen + 1; i < K;)
	{
		int ix = 0;
		for (int j = 0; j < crclen + 1; j++)
		{
			if (re[j] != 0)
			{
				ix = j;
				break;
			}
		}
		if (i + ix - K > 0)
		{
			for (int j = 0; j < crclen - (K - i) + 1; j++)
			{
				re[j] = re[j + K - i];
			}
			for (int j = 0; j < K - i; j++)
			{
				re[crclen + 1 - (K - i) + j] = 0;
			}
			break;
		}
		else {
			if (i + ix - K == 0)
			{

				for (int j = 0; j < crclen + 1 - ix; j++)
				{
					re[j] = re[j + ix];
				}
				for (int j = 0; j < ix; j++)
				{
					re[crclen - j] = 0;
				}
				for (int id = 0; id < crclen + 1; id++)
				{
					if (re[id] != poly[id])
					{
						for (int j = 0; j < crclen + 1; j++)
						{
							re[j] = re[j] ^ poly[j];
						}
						break;
					}
					if (id == crclen)
					{
						for (int j = 0; j < crclen + 1; j++)
						{
							re[j] = 0;
						}
					}
				}
				break;
			}
		}
		if (ix == 0)
		{
			while (input[i] == 0) i++;
			if (K - i <= crclen)
			{
				break;
			}
			for (int j = 0; j < crclen + 1; j++)
			{
				re[j] = input[i++];
			}
		}
		else
		{
			for (int j = 0; j < crclen + 1 - ix; j++)
			{
				re[j] = re[ix + j];
			}
			for (int j = 0; j < ix; j++)
			{
				re[crclen + 1 - ix + j] = input[i++];
			}
		}
		for (int j = 0; j < crclen + 1; j++)
		{
			re[j] = re[j] ^ poly[j];
		}
	}
	bool flag = true;
	for (int i = 0; i < crclen; i++)
	{
		if (re[i + 1] != crc[i])
		{
			flag = false;
			break;
		}
	}
	for (int i = 0; i < crclen; i++)
	{
		crcResult[i] = re[i + 1];
	}
	delete re;
	delete crc;
	return flag;
}

void BP_decoding(double* LLR,int * Out, int M, int * frozen,int iternum,int * poly,int crclen,int* frozen_ori,int * inx,int *flg)
{
	int m = (int)log2(M);
	double ** L = new double*[m+1];
	for (int i = 0; i < m+1; i++)
	{
		L[i] = new double[M];
	}
	double ** R = new double*[m+1];
	for (int i = 0; i < m+1; i++)
	{
		R[i] = new double[M];
	}
	for (int i = 0; i < m+1; i++)
	{
		for (int j = 0; j < M; j++)
		{
			R[i][j] = 0;
		}
	}
	int infoLen = M;
	for (int i = 0; i < M; i++)
	{
		L[m][i] = LLR[i];
		if (frozen[i] == 0)
		{
			R[0][i] = INF;
			infoLen = infoLen - 1;
		}
	}
	
	int * input = new int[infoLen];
	int * input2 = new int[M];
	
	double * minValue = new double[4];        // |LLR|
	int * minIndex = new int[4];             //  positions
	double * mf = new double[4];             //  sign of LLR
	int * LminIndex = new int[4];            // last convergence positions
	for(int i=0;i<4;i++)
	{
		minValue[i] = 1000;
		minIndex[i] = -1;
		mf[i] = 0;
		LminIndex[i] = 0;
	}
	int *LdeCRC = new int[crclen];
	int *CdeCRC = new int[crclen];
	for (int i = 0; i < crclen; i++)
	{
		LdeCRC[i] = 0;               // last iteration CRC
		CdeCRC[i] = 0;               // current iteration CRC
	}

	bool conFlag = false;
	double reverse = 1;
	int cfactor = 0;        // convergence factor
	double mval = 20;       // M_0
		
	for (int ite = 0; ite < iternum; ite++)
	{	
		if (conFlag)     //First-level Detection
		{
			for (int i = 0; i < M; i++)
			{
				L[m][i] = LLR[i];
			}
			conFlag = false;
			if (minIndex[0] == LminIndex[0] && minIndex[1] == LminIndex[1] && minIndex[2] == LminIndex[2] && minIndex[3] == LminIndex[3])        // Second-level Detection
			{
				reverse = -1;
			}
			else
			{
				reverse = 1;
			}
			cfactor = cfactor + 1;
            //===================Two-Level Post-Processing Operations==============================//
			if (cfactor % 4 == 0)
			{
				L[m][minIndex[0]] = reverse*mf[0]*mval;
				L[m][minIndex[1]] = mf[1]*mval;
				L[m][minIndex[2]] = mf[2]*mval;
				L[m][minIndex[3]] = mf[3]*mval;
				
			}
			if (cfactor % 4 == 1)
			{
				L[m][minIndex[0]] = reverse*mf[0]*mval;	
				L[m][minIndex[1]] = reverse*mf[1]*mval;							
				L[m][minIndex[2]] = mf[2]*mval;
				L[m][minIndex[3]] = mf[3]*mval;
			}
			if (cfactor % 4 == 2)
			{
				L[m][minIndex[0]] = reverse*mf[0]*mval;
				L[m][minIndex[1]] = reverse*mf[1]*mval;
				L[m][minIndex[2]] = reverse*mf[2]*mval;						
				L[m][minIndex[3]] = mf[3]*mval;
			}
			if (cfactor % 4 == 3)
			{
				L[m][minIndex[0]] = reverse*mf[0]*mval;
				L[m][minIndex[1]] = reverse*mf[1]*mval;
				L[m][minIndex[2]] = reverse*mf[2]*mval;
				L[m][minIndex[3]] = reverse*mf[3]*mval;			
			}
			LminIndex[0] = minIndex[0];
			LminIndex[1] = minIndex[1];
			LminIndex[2] = minIndex[2];
			LminIndex[3] = minIndex[3];			
		}
		for (int i = m-1; i >=0; i--)				
		{
			int pj = 1 << (i);
			for (int j = 0; j < (1 << (m-i-1)); j++)
			{
				
				for (int k = 0; k < (1 << (i)); k++)
				{
					L[i][k + j*pj*2] = boxOperate(L[i+1][k + j*pj * 2], L[i+1][k+pj+j*pj * 2] + R[i][k + pj + j*pj * 2]);
					L[i][k+pj + j*pj*2] = L[i + 1][k + pj + j*pj * 2] + boxOperate(L[i + 1][k + j*pj * 2], R[i][k + j*pj * 2]);
				}
			}
		}
		for (int i = 0; i < m; i++)
		{
			int pj = 1 << (i);
			for (int j = 0; j < (1 << (m - i - 1)); j++)
			{
				
				for (int k = 0; k < (1 << (i)); k++)
				{

					R[i+1][k + j*pj * 2] = boxOperate(L[i+1][k + pj + j*pj * 2] + R[i][k +pj+ j*pj*2], R[i][k + j*pj * 2]);

					R[i+1][k + pj + j*pj * 2] = R[i][k + pj + j*pj * 2] + boxOperate(L[i+1][k + j*pj * 2], R[i][k + j*pj * 2]);
				}
			}
		}
		
		for (int i = 0; i < M; i++)
		{
			if (frozen[i] == 1)
			{
				input2[inx[i] - 1] = L[0][i] >= 0 ? 0 : 1;
			}
			else
			{
				input2[inx[i] - 1] = 0;
			}
		}
		
		int ix = 0;
		for (int i = 0; i < M; i++)
		{
			if (frozen_ori[i] == 1)
			{
				input[ix++] = input2[i];
			}
		}
		// CRC Check and Compute
		bool CFlag = crcDector(input, infoLen, poly, crclen,CdeCRC);
		if (CFlag)
		{
			flg[0] = 0;
			flg[1] = ite;
			break;
		}
		else
		{
			bool crcflag = false;
			for (int i = 0; i < crclen; i++)
			{				
				if (LdeCRC[i] != CdeCRC[i])
				{
					crcflag = true;
					break;
				}
			}			
			if (!crcflag)
			{
				conFlag = true; 				
				//======================= select the four minimum values===========================//
				//======divide N LLRs into four groups uniformly,===================================//
				//======and search the minimum in each group.    ==================================//
				for(int i=0;i<4;i++)
				{
					minValue[i] = 1000;
					minIndex[i] = -1;
					mf[i] = 0;
				}
				for (int i = 0; i < M / 4; i++)
				{
					double tmp = fabs(L[m][i] + R[m][i]);
					if (minValue[0] > tmp)
					{
						minValue[0] = tmp;
						minIndex[0] = i;
						mf[0] = (L[m][i] + R[m][i]) >= 0 ? -1 : 1;
					}
				}
				for (int i = M / 4; i < M/2; i++)
				{
					double tmp = fabs(L[m][i] + R[m][i]);
					if (minValue[1] > tmp)
					{
						minValue[1] = tmp;
						minIndex[1] = i;
						mf[1] = (L[m][i] + R[m][i]) >= 0 ? -1 : 1;
					}
				}
				for (int i = M/2; i < 3*M / 4; i++)
				{
					double tmp = fabs(L[m][i] + R[m][i]);
					if (minValue[2] > tmp)
					{
						minValue[2] = tmp;
						minIndex[2] = i;
						mf[2] = (L[m][i] + R[m][i]) >= 0 ? -1 : 1;
					}
				}
				for (int i = 3*M / 4; i < M; i++)
				{
					double tmp = fabs(L[m][i] + R[m][i]);
					if (minValue[3] > tmp)
					{
						minValue[3] = tmp;
						minIndex[3] = i;
						mf[3] = (L[m][i] + R[m][i]) >= 0 ? -1 : 1;
					}
				}
							
				//============sort m0,m1,m2,m3;===================================================//
				double tmpM[4];
				int tmpIndex[4];
				for(int i=0;i<4;i++)
				{
					tmpM[i] = minValue[i];
					tmpIndex[i] = i;
				}
				if(tmpM[0]>tmpM[1])
				{
					int Temp = tmpIndex[0];
					tmpIndex[0] = tmpIndex[1];
					tmpIndex[1] = Temp;
					
					double TempValue = tmpM[0];
					tmpM[0] = tmpM[1];
					tmpM[1] = TempValue;
				}
				if(tmpM[2]>tmpM[3])
				{
					int Temp = tmpIndex[2];
					tmpIndex[2] = tmpIndex[3];
					tmpIndex[3] = Temp;
					
					double TempValue = tmpM[2];
					tmpM[2] = tmpM[3];
					tmpM[3] = TempValue;
				}
				if(tmpM[0]>tmpM[2])
				{
					int Temp = tmpIndex[0];
					tmpIndex[0] = tmpIndex[2];
					tmpIndex[2] = Temp;
					
					double TempValue = tmpM[0];
					tmpM[0] = tmpM[2];
					tmpM[2] = TempValue;
				}
				if(tmpM[1]>tmpM[3])
				{
					int Temp = tmpIndex[1];
					tmpIndex[1] = tmpIndex[3];
					tmpIndex[3] = Temp;
					
					double TempValue = tmpM[1];
					tmpM[1] = tmpM[3];
					tmpM[3] = TempValue;
				}
				if(tmpM[1]>tmpM[2])
				{
					int Temp = tmpIndex[1];
					tmpIndex[1] = tmpIndex[2];
					tmpIndex[2] = Temp;
					
					double TempValue = tmpM[1];
					tmpM[1] = tmpM[2];
					tmpM[2] = TempValue;
				}
				
				double Temp2Value[4];
				int Temp2Index[4];
				double Temp2Flag[4];
				for(int i=0;i<4;i++)
				{
					Temp2Value[i] = minValue[tmpIndex[i]];
					Temp2Index[i] = minIndex[tmpIndex[i]];
					Temp2Flag[i] = mf[tmpIndex[i]];
				}
				
				//=============select m0=======================================================// 
				minValue[0] = Temp2Value[0];
				minIndex[0] = Temp2Index[0];
				mf[0] = Temp2Flag[0];
				
				//=============select m1=======================================================//
				for(int i=tmpIndex[0]*M/4;i<(tmpIndex[0]+1)*M/4;i++)
				{
					double tmp = fabs(L[m][i] + R[m][i]);
					if(Temp2Value[1]>tmp && i!=minIndex[0])
					{
						Temp2Value[1] = tmp;
						Temp2Index[1] = i;
						Temp2Flag[1] = (L[m][i] + R[m][i]) >= 0 ? -1 : 1;
					}
				}
				minValue[1] = Temp2Value[1];
				minIndex[1] = Temp2Index[1];
				mf[1] = Temp2Flag[1];
				
				//============select m2=========================================================//
				for(int i=tmpIndex[1]*M/4;i<(tmpIndex[1]+1)*M/4;i++)
				{
					double tmp = fabs(L[m][i] + R[m][i]);
					if(Temp2Value[2]>tmp && i!=minIndex[1])
					{
						Temp2Value[2] = tmp;
						Temp2Index[2] = i;
						Temp2Flag[2] = (L[m][i] + R[m][i]) >= 0 ? -1 : 1;
					}
				}
				for(int i=tmpIndex[0]*M/4;i<(tmpIndex[0]+1)*M/4;i++)
				{
					double tmp = fabs(L[m][i] + R[m][i]);
					if(Temp2Value[2]>tmp && i!=minIndex[0] && i!=minIndex[1])
					{
						Temp2Value[2] = tmp;
						Temp2Index[2] = i;
						Temp2Flag[2] = (L[m][i] + R[m][i]) >= 0 ? -1 : 1;
					}
				}
				minValue[2] = Temp2Value[2];
				minIndex[2] = Temp2Index[2];
				mf[2] = Temp2Flag[2];
				
				//============ select m3 =======================================================//
				for(int i=tmpIndex[2]*M/4;i<(tmpIndex[2]+1)*M/4;i++)
				{
					double tmp = fabs(L[m][i] + R[m][i]);
					if(Temp2Value[3]>tmp && i!=minIndex[2])
					{
						Temp2Value[3] = tmp;
						Temp2Index[3] = i;
						Temp2Flag[3] = (L[m][i] + R[m][i]) >= 0 ? -1 : 1;
					}
				}
				for(int i=tmpIndex[1]*M/4;i<(tmpIndex[1]+1)*M/4;i++)
				{
					double tmp = fabs(L[m][i] + R[m][i]);
					if(Temp2Value[3]>tmp && i!=minIndex[2] && i!=minIndex[1])
					{
						Temp2Value[3] = tmp;
						Temp2Index[3] = i;
						Temp2Flag[3] = (L[m][i] + R[m][i]) >= 0 ? -1 : 1;
					}
				}
				for(int i=tmpIndex[0]*M/4;i<(tmpIndex[0]+1)*M/4;i++)
				{
					double tmp = fabs(L[m][i] + R[m][i]);
					if(Temp2Value[3]>tmp && i!=minIndex[2] && i!=minIndex[1] && i!=minIndex[0])
					{
						Temp2Value[3] = tmp;
						Temp2Index[3] = i;
						Temp2Flag[3] = (L[m][i] + R[m][i]) >= 0 ? -1 : 1;
					}
				}				
				minValue[3] = Temp2Value[3];
				minIndex[3] = Temp2Index[3];
				mf[3] = Temp2Flag[3];	
                //============================================================================//				
			}
			else
			{
				for (int i = 0; i < crclen; i++)
				{
					LdeCRC[i] = CdeCRC[i];
				}
			}
		}
	}
	for (int i = 0; i < M; i++)
	{
		if (L[0][i] + R[0][i] > 0)
			Out[inx[i] - 1] = 0;
		else
			Out[inx[i] - 1] = 1;
	}
	
	delete minValue;
	delete minIndex;
	delete mf;
	delete LminIndex;	
	delete input;
	delete input2;
	delete LdeCRC;
	delete CdeCRC;
	for (int i = 0; i < m+1; i++)
	{
		delete R[i];
	}
	delete R;

	for (int i = 0; i < m+1; i++)
	{
		delete L[i];
	}
	delete L;
}

void mexFunction(int output_size, mxArray *output[], int input_size, const mxArray *input[])
{
	double *LLR_msg = mxGetPr(input[0]);                 // channel LLR
	double *frozen_set = mxGetPr(input[1]);              // Frozen set
	double *dN = mxGetPr(input[2]);                      // code length
	double *diter = mxGetPr(input[3]);                   // max iteration number
	double *dpoly = mxGetPr(input[4]);                   // crc polynomial   
	double *dcrclen = mxGetPr(input[5]);                 // crc length
	double *dindex = mxGetPr(input[6]);                  // index permutation
	double * fro_ori = mxGetPr(input[7]);                // Original frozen set  
	//////////////////////////////////////
	int N = (int)(*dN);
	int iternum = (int)(*diter);
	int crclen = (int)(*dcrclen);
	/////////////////////////////////////  
	output[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
	output[1] = mxCreateDoubleMatrix(2, 1, mxREAL);
	double *outData = mxGetPr(output[0]);
	double *outflg = mxGetPr(output[1]);
	
	int *hd_dec = new int[N];
	int * frozen = new int[N];
	int *flg = new int[2];
	flg[0] = N;
	flg[1] = iternum;
	int* index = new int[N];
	for (int i = 0; i < N; i++)
	{
		index[i] = (int)dindex[i];
	}
	int * frozen_ori = new int[N];
	for (int i = 0; i < N; i++)
	{
		frozen_ori[i] = (int)fro_ori[i];
	}
	for (int i = 0; i < N; i++)
	{
		frozen[i] = (int)frozen_set[i];
	}
	double * LLR = new double[N];
	for (int i = 0; i < N; i++)
	{
		LLR[i] = LLR_msg[i];
	}
	int * poly = new int[crclen + 1];
	for (int i = 0; i < crclen + 1; i++)
	{
		poly[i] = (int)dpoly[i];
	}

    BP_decoding(LLR, hd_dec, N, frozen, iternum, poly, crclen,frozen_ori,index,flg);

	for (int i = 0; i<N; i++)
	{
		outData[i] = (double)hd_dec[i];
	}
	outflg[0] = flg[0];
	outflg[1] = flg[1];
	delete frozen_ori;
	delete index;
	delete poly;
	delete flg;
	delete LLR;
	delete frozen;
	delete hd_dec;
}