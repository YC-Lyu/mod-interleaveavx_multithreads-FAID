#include "CSimulate.h"

#define BF_ITER_COUNT 1
#define FAKE_ENCODE 0

/*
        Seed is used to initial the Gaussian Noise seed that used in traditional
        C function
*/

static int seed[] = { 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337,
    347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467,
    479, 487, 491, 499, 503, 507, 521, 523, 541, 547, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631,
    641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787,
    797, 809, 811, 821, 823, 827, 829, 839, 953, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947,
    953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019 };
CSimulate::CSimulate()
{
    ldpc = nullptr;
    channel = nullptr;
    modulate = nullptr;
    TestFrame = 0;
    ErrorFrame = 0;
    ErrorBits = 0;
    ModErrorBits = 0;
    ModErrorSymbol = 0;
    ModErrorFrame = 0;
    LT3ErrBitFrame = 0;
    scale = 0;
}
CSimulate::~CSimulate()
{
    delete ldpc;
    delete channel;
    delete modulate;
}
/*
        Initial Module
*/
void CSimulate::Initial(Parameter_Simulation& p, int index)
{
    scale = p.scale;
    ldpc = new CLDPC();
    channel = new CChannel();
    modulate = new CModulate();
    cpu_id = index;
    ldpc->Initial(p.nb_frames, p.Max_Iteration);
    modulate->ModulationType = p.mod_type;
    modulate->InterleaveModType = p.interleavemod_type;
    /*
    DATE:20190306
    Modualte:the modulated bits number are the realtransmited bits
    速率匹配：打孔和shorten实际没有参与
    */
    modulate->Initial(ldpc->m_frame * (ldpc->m_N - ldpc->m_PunLen - ldpc->m_ShortenLen));
    channel->RandomSeed = seed[index];
    channel->Initial(modulate->SymbolLen, index);
}

/*
        FunctionName: Configure
        Function:set  Eb_N0

*/

void CSimulate::Configure(float Eb_N0, const int _decode_method)
{
    snr = Eb_N0;
    if (modulate->ModulationType == 1) {
        sigma = (float)(1.0 / sqrt(2.0 * ldpc->m_Rate * modulate->ModulationType * pow(10.0, 0.1 * snr)));
    } else {
        sigma = (float)(1.0 / sqrt(ldpc->m_Rate * modulate->ModulationType * pow(10.0, 0.1 * snr)));
    }
    decode_method = _decode_method;
    /*
            Coding and modulation
            All zero coding, so the Modulation symbols are same
            Encoding and modulation only need initialize once
    */
    // ldpc->FakeEncoder();
    // modulate->BPSKModulation(ldpc->outputBits);
    // Reset the statistics
    TestFrame = 0;
    ErrorFrame = 0;
    ErrorBits = 0;
    ModErrorBits = 0;
    ModErrorSymbol = 0;
    ModErrorFrame = 0;
    LT3ErrBitFrame = 0;
}
void CSimulate::Run()
{
    int i;
    Statistic Test;
    ModStatistic ModTest;
#if BF_ITER_COUNT == 1
    ofstream iterOut;
    int BFiters_[51] = { 0 }; // 统计 BF 迭代次数
#endif
    int BFiter;

#if FAKE_ENCODE == 1
    ldpc->FakeEncoder();
#else
    ldpc->GenMsgSeq();
    ldpc->Encode();
#endif
    /*编码是固定的，因此只需要调制及交织1次*/
    //编码后打孔的序列存入outputbits进入交织器，存入modulate->interleaveseq，outputbit存放格式32信息+32校验
    if (modulate->ModulationType == 1) {
        modulate->BPSKModulation(ldpc->outputBits);
    } else {
        modulate->BeforeModulationInterleaver(ldpc->outputBits);
        modulate->Modulation(modulate->InterLeaveSeq); //调制输入为交织后的interleaveseq
    }
    for (i = 0; i < 50; ++i) {

        TestFrame += 32;

        if (modulate->ModulationType == 1) {
            // 设计思想是：初始化阶段完成编码调制以后，由于是BPSK调制，因此我们只是取幅度，所以不需要解调了。
            channel->BPSKAWGNChannel(modulate->BPSKModSeq, sigma); //一个一个出
            ldpc->float2LimitChar_4bit(ldpc->fixInput, channel->BPSKSymbol, scale, BitsOverChannel * 32);
        } else {
            channel->AWGNChannel(modulate->ModSeq, sigma / sqrt(2)); // 复信号虚实部各加噪声方差为原来的1/2
            modulate->Demodulation(channel->SymbolSeq);
            modulate->AfterDeModulationDeInterleaver();
            // ModTest = modulate->ModCalErr(modulate->DeInterLeaveSeq, ldpc->outputBits); // 32*bitsoverchannel
            // ModErrorSymbol += ModTest.ErrorSymbol;
            // ModErrorBits += ModTest.ErrorBits;
            ldpc->float2LimitChar_4bit(ldpc->fixInput, modulate->DeInterLeaveSeq, scale, BitsOverChannel * 32);
        }

        // BitsOverChannel=Novar-puncture-shorten
        switch (decode_method) {
        case 0:
            ldpc->Decode(); // NMS
            break;
        case 1:
            ldpc->Decode_OMS(); // OMS
            break;
        case 2:
            ldpc->Decode_FAID(); // FAID+DTBF
            break;
        case 3:
            BFiter = ldpc->Decode_OMSBF(); // OMS+BF
#if BF_ITER_COUNT == 1
            BFiters_[BFiter]++;
#endif
            break;
        case 4:
            BFiter = ldpc->Decode_OMS_DTBF(); // OMS+DTBF
#if BF_ITER_COUNT == 1
            BFiters_[BFiter]++;
#endif
            break;
        case 5:
            ldpc->Decode_FAID_2B1C(); // FAID+2B1C
            break;
        default:
            ldpc->Decode(); // NMS
            break;
        }

        Test = ldpc->CalculateErrors(modulate->DeInterLeaveSeq, ldpc->fixInput, collectflag); // 32个码字的串行比较
        ErrorFrame += Test.ErrorFrame;
        ErrorBits += Test.ErrorBits;
        LT3ErrBitFrame += Test.LT3ErrBitFrame;
    }
#if BF_ITER_COUNT == 1
    iterOut.open("iterCount.txt", std::ios::app);
    for (int i = 1; i <= 50; i++) {
        if (BFiters_[i] != 0) {
            iterOut << i << ": " << BFiters_[i] << std::endl;
        }
    }
    iterOut.close();
#endif
}
///*
// date:20190415
//如果出现某一帧错误比特小于3（自行设定的固定值）
//则将该线程的迭代次数增加到60(decode1()中自行设定）
//*/
//       if (LT3ErrBitFrame > 0)
//		{
//			break;
//		}
//	}
//
//	   if (LT3ErrBitFrame > 0)
//	   {
//
//		   TestFrame = 0;
//		   ErrorFrame = 0;
//		   ErrorBits = 0;
//		   LT3ErrBitFrame = 0;
//		   for (i = 0; i < 10000; ++i)
//			{
//				TestFrame += 32;
//				channel->BPSKAWGNChannel(modulate->BPSKModSeq, sigma);
//				ldpc->float2LimitChar_6bit(ldpc->fixInput, channel->BPSKSymbol, scale;, BitsOverChannel
//* 32); 				ldpc->Decode1();//如果出现某一帧错误小于3比特，则将迭代次数增加为60，重新译码
// Test = ldpc->CalculateErrors(channel->BPSKSymbol, ldpc->fixInput);//32个码字的串行比较
// ErrorFrame += Test.ErrorFrame; 				ErrorBits += Test.ErrorBits;
// LT3ErrBitFrame
// += Test.LT3ErrBitFrame;
//			}
//	    }
//}

/*
        FunctionName: set_cpu
        Function: Set CPU affinity

*/
void CSimulate::set_cpu(int i)
{
    // cout << "set cpu" << endl;
    // getchar();

    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(i, &mask);
    if (pthread_setaffinity_np(pthread_self(), sizeof(mask), &mask) < 0) {
        fprintf(stderr, "set thread affinity failed\n");
        getchar();
        exit(EXIT_FAILURE);
    }
}

/*
        FunctionName: start_thread
        Function: Run
*/

void* CSimulate::start_thread(void* arg)
{
    // cout << "start_thread" << endl;
    // getchar();

    CSimulate* ptr = (CSimulate*)arg;
    // cout << "CSimulate* ptr = (CSimulate*)arg; end" << endl;
    // getchar();

    // cout << "start_thread ptr->Run();" << endl;
    // getchar();
    ptr->Run();
    return NULL;
}
void CSimulate::End() { pthread_join(pid, NULL); }
// TestFrame and ErrorFrame are set to zero for each start and the
// value are store in Main function
int CSimulate::Start()
{
    TestFrame = 0;
    ErrorFrame = 0;
    ErrorBits = 0;
    LT3ErrBitFrame = 0;
    ModErrorBits = 0;
    ModErrorSymbol = 0;
    ModErrorFrame = 0;

    // cout << "start set_cpu" << endl;
    // getchar();
    set_cpu(cpu_id);

    // cout << "start  start_thread" << endl;
    // getchar();
    if (pthread_create(&pid, NULL, start_thread, (void*)this) != 0) //创建线程
    {
        std::cout << "pthread_create error\n";
        getchar();
        return -1;
    }
    return 0;
}
