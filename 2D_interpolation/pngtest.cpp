#include "PngProc.h"
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <fstream>

using namespace std;

void ResizeBilinear(const unsigned char* pbIn, int lWidthIn,
	int lHeightIn, unsigned char* pbOut, int lWidthOut, int lHeightOut)
{
	for (int i = 0; i < lHeightOut; ++i)
	{
		double yy = (double)i * (double)lHeightIn / lHeightOut;
		int y = floor(yy);	// целая часть yy 
		double u = yy - (double)y; // дробная часть yy 
		for (int j = 0; j < lWidthOut; ++j)
		{
			double xx = (double)j * (double)lWidthIn / lWidthOut;
			int x = floor(xx);	// целая часть xx 
			double v = xx - (double)x; // дробная часть xx
			int lP00 = pbIn[y * lWidthIn + x],
				lP01 = pbIn[y * lWidthIn + x + 1],
				lP10 = pbIn[(y + 1) * lWidthIn + x],
				lP11 = pbIn[(y + 1) * lWidthIn + x + 1];
			pbOut[i * lWidthOut + j] = (unsigned char)((1. - u) * (1. - v) * lP00 + u * (1. - v) * lP10 + v * (1. - u) * lP01 +
				u * v * lP11);
		}
	}
	return;
}

void ResizeNearest(const unsigned char* pbIn, int lWidthIn,
	int lHeightIn, unsigned char* pbOut, int lWidthOut, int lHeightOut)
{
	for (int i = 0; i < lHeightOut; ++i)
	{
		double yy = (double)i * (double)lHeightIn / lHeightOut;
		int y = floor(yy);
		double u = yy - (double)y;
		for (int j = 0; j < lWidthOut; ++j)
		{
			double xx = (double)j * (double)lWidthIn / lWidthOut;
			int x = floor(xx);
			double v = xx - (double)x;
			int lP00 = pbIn[y * lWidthIn + x],
				lP01 = pbIn[y * lWidthIn + x + 1],
				lP10 = pbIn[(y + 1) * lWidthIn + x],
				lP11 = pbIn[(y + 1) * lWidthIn + x + 1];
			if (((int)(xx + 0.5) - x) > 0) {
				if (((int)(yy + 0.5) - y) > 0) {
					pbOut[i * lWidthOut + j] = lP11;
				}
				else {
					pbOut[i * lWidthOut + j] = lP10;
				}
			}
			else {
				if (((int)(yy + 0.5) - y) > 0) {
					pbOut[i * lWidthOut + j] = lP01;
				}
				else {
					pbOut[i * lWidthOut + j] = lP00;
				}
			}
		}
	}
	return;
}


int main(int argc, char* argv[])
{
	class CBitsPtrGuard
	{
	public:
		CBitsPtrGuard(unsigned char** pB) : m_ppBits(pB) { }
		~CBitsPtrGuard() { if (*m_ppBits) delete* m_ppBits, * m_ppBits = 0; }
	protected:
		unsigned char** m_ppBits;
	};

	// parse input parameters
	char	szInputFileName[256];
	char    szOutputFileName[256];
	if (argc < 2)
		printf("\nformat: pngtest <input_file> [<output_file>]");
	else
	{
		strcpy(szInputFileName, argv[1]);
		if (argc > 2)
			strcpy(szOutputFileName, argv[2]);
		else
		{
			strcpy(szOutputFileName, szInputFileName);
			strcat(szOutputFileName, "_out.png");
		}
	}


	size_t nReqSize = NPngProc::readPngFile(szInputFileName, 0, 0, 0, 0);
	if (nReqSize == NPngProc::PNG_ERROR)
	{
		printf("\nError ocured while pngfile was read");
		return -1;
	}


	unsigned char* pInputBits = new unsigned char[nReqSize];
	if (!pInputBits)
	{
		printf("\nCan't allocate memory for image, required size is %u", nReqSize);
		return -1;
	}
	CBitsPtrGuard InputBitsPtrGuard(&pInputBits);


	size_t nWidth, nHeight;
	unsigned int nBPP;

	size_t nRetSize = NPngProc::readPngFileGray(szInputFileName, pInputBits, &nWidth, &nHeight/*, &nBPP*/);
	nBPP = 8;

	//new w & h
	int outWidth = nWidth * 4;
	int outHeight = nHeight * 4;
	const int OutSize = outWidth * outHeight;

	unsigned char* pOutputBits = new unsigned char[OutSize];
	if (!pOutputBits)
	{
		printf("\nCan't allocate memory for image, required size is %u", nReqSize);
		return -1;
	}


	CBitsPtrGuard OutputBitsPtrGuard(&pOutputBits);

//раскомментить это для интерполяии ближайшим соседом
	//ResizeNearest(pInputBits, nWidth, nHeight, pOutputBits, outWidth, outHeight);

//раскомментить это для билинейной интерполяии
	//ResizeBilinear(pInputBits, nWidth, nHeight, pOutputBits, outWidth, outHeight);

	if (NPngProc::writePngFile(szOutputFileName, pOutputBits, outWidth, outHeight, nBPP) == NPngProc::PNG_ERROR)
	{
		printf("\nError ocuured during png file was written");
		return -1;
	}

	return 0;
}
