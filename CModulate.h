#ifndef CMODULATE_H
#define CMODULATE_H
#include"mkl.h"
#include<cstdint>
#include<iostream>
#include"CTool.h"
#include<cstdlib>
#include<cmath>
#include"./Constants/Constants_SSE.h"
#include"CChannel.h"
#include"CLDPC.h"

using namespace std;
/* CModulate */
struct ModStatistic
{
	unsigned long ErrorSymbol;
	unsigned long ErrorBits;
	unsigned long ErrorFrame;
	char padding[48];// padding data no use
};

class CModulate
{
public:

	float * DemodSeq;//DemodSeq
	int8_t* InterLeaveSeq;
	float * DeInterLeaveSeq;
	int8_t* ILSeq;
	float * DILSeq;
	MKL_Complex8 *ModSeq;//Higher Order modulation
	float* BPSKModSeq;//BPSK modualtion symbols
	unsigned long SourceLen;// The length of the modualtion bits
	unsigned long SymbolLen;// The length of the modulation symbols
	int ModulationType;// The moduation type
	int InterleaveModType;
	char padding[20];// padding chars


	CModulate();
	~CModulate();
	void Initial(unsigned long codelen);//Initialization

	void BPSKModulation(int8_t*SourceSeq);//BPSK modulation
	void BPSKDemodulation(float* ReceivedSeq);//BPSK demodulation
	void Modulation(int8_t* ReceivedSeq);// High order modulation
	void Demodulation(MKL_Complex8* ReceivedSeq);// High order demodulation

	void BeforeModulationInterleaver(int8_t *codeseq);
	void AfterDeModulationDeInterleaver();

	ModStatistic ModCalErr(float *demodseq, int8_t *encodeseq);
};


#endif // !CSIMULATE_H
