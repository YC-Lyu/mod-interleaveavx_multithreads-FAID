#ifndef CSIMULATE_H
#define CSIMULATE_H

#include"CLDPC.h"
#include"CChannel.h"
#include"CModulate.h"
#include"CTool.h"
#include<cstring>
#include<pthread.h>
#include<cstdlib>
using namespace std;
/*
	class CSimulate:

*/
class CSimulate
{
public:
	CLDPC *ldpc;
	CChannel *channel;
	CModulate *modulate;
	unsigned long ErrorBits;
	unsigned long TestFrame;
	unsigned long ErrorFrame;
	unsigned long ModErrorBits;
	unsigned long ModErrorSymbol;
	unsigned long ModErrorFrame;
	unsigned long LT3ErrBitFrame;
	float snr;
	float sigma;
	float scale;
	int decode_method;
	int nb_iteration;
	int cpu_id;
	pthread_t pid;
	char padding[48];

	CSimulate();
	~CSimulate();
	void Initial(Parameter_Simulation &p, int cpu_index);//Initial parameters
	void Configure(float Eb_N0, const int _decode_method);//set Eb_N0
	void Run();// Simualate 
	void set_cpu(int i);//set cpu affinity
	static void *start_thread(void * arg);//start simulation
	int Start();
	void End();// join end

	/*void ReadProfile();*/
};

#endif // !CSIMULATE_H
