#include "CLDPC.h"
#include "Codeword.h"
#include <fstream>
#include <iostream>
#include <omp.h>
#include <stdio.h>

int vSAT_NEG_MSG = SAT_NEG_MSG;
int vSAT_POS_MSG = SAT_POS_MSG;

int vSAT_NEG_VAR = SAT_NEG_VAR;
int vSAT_POS_VAR = SAT_POS_VAR;

int vSAT_NEG_LLR = SAT_NEG_MSG;
int vSAT_POS_LLR = SAT_POS_MSG;

CLDPC::CLDPC()
{

    // scale=0;//for fixedpoint decoding
    inputBits = nullptr;
    outputBits = nullptr;
    var_nodes = nullptr;
    p_vn_adr = nullptr;
    var_msgs = nullptr;
    encoder_varnodes = nullptr; // for parallel encoding

    decodedBits = nullptr;
    fixInput = nullptr; // fixed LLR for the length is N*m_frames;
    // errorindex = nullptr;
    m_M = 0; // check
    m_N = 0; // code length
    m_K = 0; // Information bits
    m_PunLen = 0; // Puncture
    m_ShortenLen = 0; // Shorten
    m_Rate = 0; // code ratio
    m_frame = 0; // no of frames for parallel
}
CLDPC::~CLDPC()
{
    free(inputBits);
    free(outputBits);
    free(var_nodes);
    free(p_vn_adr);
    free(var_msgs);
    free(encoder_varnodes);
    free(decodedBits);
    free(fixInput);
    free(errorbitblock);
    free(errorbitindex);
    free(errorcheckblock);
    free(errorcheckindex);
    free(errorfloat);
    free(errorchar);
}
/*
    GenMsgSeq:
        Random information Bits
*/
void CLDPC::GenMsgSeq()
{
    int i;
    for (i = 0; i < (NmoinsK - _ShortenBits) * m_frame; ++i) {
        inputBits[i] = rand() % 2;
    }
}

void CLDPC::Encode()
{
    TYPE* encodedBits = encoder_varnodes;
    // set all zeros so the shorten can be assigned
    memset(encodedBits, 0, sizeof(TYPE) * NOEUD);
    int i, j, z;

    // Get The information Bits
    if ((NmoinsK - _ShortenBits) % 32 == 0) {
        uchar_transpose_avx((TYPE*)inputBits, (TYPE*)encodedBits, (NmoinsK - _ShortenBits));
    } else {
        char* ptSource = (char*)encodedBits;
        for (i = 0; i < (NmoinsK - _ShortenBits); ++i) {
            for (z = 0; z < 32; ++z) {
                ptSource[32 * i + z] = inputBits[z * (NmoinsK - _ShortenBits) + i];
            }
        }
    }

    // Gaussian elimination encoding
    const unsigned short* p = GenMatrix;
    for (i = 0; i < _NoCheck; ++i) {
        short int weight = *(p++);
        for (j = 0; j < weight; ++j) {
            encodedBits[i + NmoinsK] = _mm256_xor_si256(encodedBits[i + NmoinsK], encodedBits[*(p++)]);
        }
    }
    // the puncturebits and the shorten bits are not transmitted over the channel
    // The information bits
    // The frame structure is :The 32 information Bits The 32 Check Bits;
    //	Which maybe different from the the Version_2
    //  modified date:20170923
    if ((NmoinsK - _PunctureBits - _ShortenBits) % 32 == 0) {
        uchar_itranspose_avx(
            (TYPE*)(encodedBits + _PunctureBits), (TYPE*)outputBits, (NmoinsK - _PunctureBits - _ShortenBits));
    } else {
        char* ptr = (char*)(encodedBits + _PunctureBits); //把指针的起始位置放在puncture后
        for (i = 0; i < (NmoinsK - _PunctureBits - _ShortenBits); i += 1) {
            for (j = 0; j < 32; j += 1) {
                outputBits[j * (NmoinsK - _PunctureBits - _ShortenBits) + i] = (ptr[32 * i
                    + j]); // outputBits为编码后去除打孔和shorten，在outputbits中填充去除打孔和shorten后的编码信息
            }
        }
    }
    // The check Bits
    int8_t* pp = (outputBits + (NmoinsK - _PunctureBits - _ShortenBits) * 32);
    if (_NoCheck % 32 == 0) {
        uchar_itranspose_avx((TYPE*)(encodedBits + NmoinsK), (TYPE*)(pp), _NoCheck);
    } else {
        char* ptr = (char*)(encodedBits + NmoinsK);
        for (i = 0; i < _NoCheck; ++i)

        {
            for (j = 0; j < 32; ++j) {

                pp[j * _NoCheck + i] = ptr[32 * i + j];
            }
        }
    }
    ofstream enout;

    //	int* encodeout;
    //	encodeout = (int*)malloc(sizeof(int)*BitsOverChannel * 32);
    //	int l = 0;
    //	for(int i = 0 ; i < 32 ; i++ )
    //{
    //	enout.open("encodeout.txt", std::ios::app);
    //	for (j = 0; j < BitsOverChannel; j++)
    //	{
    //		encodeout[l] = (int)outputBits[i*BitsOverChannel + j];
    //		l++;
    //	}
    //	enout << "encodeout: ";//通过awgn信道后的值，information+check
    //	for (j = 0; j < l; j++)
    //	{
    //		enout << encodeout[j] << "\t";
    //
    //	}
    //	enout << endl;
    //	enout.close();
    //	for (j = 0; j <= l; j++)
    //	{
    //		encodeout[j] = 0;
    //	}
    //	l = 0;
    // }
    //	free(encodeout);
}

/**
 * @brief 读取固定码字，outputbits 中是先交织存储信息位，再交织存储校验位
 * @date 2022-02-22
 * @author ZhangYifan
 *
 */
void CLDPC::FakeEncoder()
{
    int i, j, z;
    TYPE* encodedBits = encoder_varnodes;
    memset(encodedBits, 0, sizeof(TYPE) * NOEUD);
    for (i = 0; i < NOEUD; ++i) {
        encodedBits[i] = VECTOR_SET1(CodeWord_sym[i]);
    }

    if ((NmoinsK - _ShortenBits) % 32 == 0) {
        uchar_itranspose_avx((TYPE*)encodedBits, (TYPE*)inputBits, (NmoinsK - _ShortenBits));
    } else {
        char* ptSource = (char*)encodedBits;
        for (i = 0; i < (NmoinsK - _ShortenBits); ++i) {
            for (z = 0; z < 32; ++z) {
                inputBits[z * (NmoinsK - _ShortenBits) + i] = ptSource[32 * i + z];
            }
        }
    }
    if ((NmoinsK - _PunctureBits - _ShortenBits) % 32 == 0) {
        uchar_itranspose_avx(
            (TYPE*)(encodedBits + _PunctureBits), (TYPE*)outputBits, (NmoinsK - _PunctureBits - _ShortenBits));
    } else {
        char* ptr = (char*)(encodedBits + _PunctureBits); //把指针的起始位置放在puncture后
        for (i = 0; i < (NmoinsK - _PunctureBits - _ShortenBits); i += 1) {
            for (j = 0; j < 32; j += 1) {
                outputBits[j * (NmoinsK - _PunctureBits - _ShortenBits) + i] = (ptr[32 * i + j]);
                // outputBits为编码后去除打孔和shorten，在outputbits中填充去除打孔和shorten后的编码信息
            }
        }
    }

    // The check Bits
    int8_t* pp = (outputBits + (NmoinsK - _PunctureBits - _ShortenBits) * 32);
    if (_NoCheck % 32 == 0) {
        uchar_itranspose_avx((TYPE*)(encodedBits + NmoinsK), (TYPE*)(pp), _NoCheck);
    } else {
        char* ptr = (char*)(encodedBits + NmoinsK);
        for (i = 0; i < _NoCheck; ++i) {
            for (j = 0; j < 32; ++j) {
                pp[j * _NoCheck + i] = ptr[32 * i + j];
            }
        }
    }
}

/*
        Decode
        parallel decoing 32 codewords
        迭代次数读取Profile
*/
void CLDPC::Decode()
{
    Parameter_Simulation p_simulation;
    ReadProfile(&p_simulation);

    int8_t* ptr;
    int8_t* pp;
    int factor_1 = p_simulation.Factor_1;
    int factor_2 = p_simulation.Factor_2;
    const TYPE zero = VECTOR_ZERO;
    //#define SAT_NEG_MSG(-(0x0001 << (NB_BITS_MESSAGES - 1)) + 1) -31(6bit量化）
    const TYPE minMesg = VECTOR_SET1(SAT_NEG_MSG);
    size_t i, j, z;

    /* Lmn Initialization */
    for (i = 0; i < MESSAGE; ++i) {
        var_msgs[i] = zero;
    }
    // The information part 交织 解调初始似然比排列格式
    // 第一个码字第一个对数似然比，第二个码字第一个，……第32个码字第一个，第一个码字第二个……
    if ((NmoinsK - _PunctureBits - _ShortenBits) % 32 == 0) {
        uchar_transpose_avx(
            (TYPE*)fixInput, (TYPE*)(var_nodes + _PunctureBits), (NmoinsK - _PunctureBits - _ShortenBits));
    } else {
        ptr = (int8_t*)(var_nodes + _PunctureBits); //指针起始位置在puncture后，即打孔在最前
        for (i = 0; i < (NmoinsK - _PunctureBits - _ShortenBits); ++i) {
            for (j = 0; j < 32; ++j) {
                //除去puncture和shorten的原始似然比
                ptr[i * 32 + j] = fixInput[j * (NmoinsK - _PunctureBits - _ShortenBits) + i];
            }
        }
    }
    // The check part
    pp = fixInput + (NmoinsK - _PunctureBits - _ShortenBits) * 32;
    // pp = fixInput + (NmoinsK - _PunctureBits - _ShortenBits);
    if (_NoCheck % 32 == 0) {
        uchar_transpose_avx((TYPE*)(pp), (TYPE*)(var_nodes + NmoinsK), _NoCheck);
    } else {
        ptr = (int8_t*)(var_nodes + NmoinsK);
        for (i = 0; i < _NoCheck; ++i) {
            for (j = 0; j < 32; ++j) {
                ptr[i * 32 + j] = pp[j * _NoCheck + i];
            }
        }
    }

    // // puncture 最前 llr=0
    // for (i = 0; i < _PunctureBits; ++i) {
    //     var_nodes[i] = zero;
    // }

    // // shorten 最后 llr=minMesg
    // for (i = NmoinsK - _ShortenBits; i < NmoinsK; ++i) {
    //     var_nodes[i] = minMesg;
    // }

    for (i = 0; i < 384; ++i) {
        var_nodes[_NoVar-1-i] = zero;
    }

    int nombre_iterations = nb_iteration;
    // Fixed iterations
    while (nombre_iterations--) { //一次迭代开始，没有动态终止
        TYPE* p_msg1r = var_msgs; // read Lmn
        TYPE* p_msg1w = var_msgs; // write Lmn
#if PETIT == 1
        TYPE** p_indice_nod1 = p_vn_adr; // p_vn_adr: Stores the En,read
        TYPE** p_indice_nod2 = p_vn_adr; // write

#else
        const unsigned short* p_indice_nod1 = PosNoeudsVariable;
        const unsigned short* p_indice_nod2 = PosNoeudsVariable;
#endif

        const TYPE min_var = VECTOR_SET1(vSAT_NEG_VAR); // semimimmum value 8bit -127
        const TYPE max_msg = VECTOR_SET1(vSAT_POS_MSG); // maximum value 6bit +31
        const TYPE max_var = VECTOR_SET1(vSAT_POS_VAR); // maximum value 8 bit 127
        // DEG_1_COMPUTATIONS: The number of the rows of degree DEG_1
        /**************************************DEG_1****************************************************************/
        for (i = 0; i < DEG_1_COMPUTATIONS; i++) {
            TYPE tab_vContr[DEG_1]; // DEG1
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR); // 8bit maximum value  +127
            TYPE min2 = min1; // min2=min1=127

#if (DEG_1 & 0x01) == 1
            const unsigned char sign8 = 0x80; // sign8 =128=1000000B
            // isign =12*16=11000000B,之所以使用的是X100,0000,是为了如果他的符号为是0，则在进行
            // _mm_sign_epi8运算时，符合要求。
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80;
            const unsigned char isign8b = 0x40; // isign8b=01000000B=64
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

#if PETIT == 1
#if MANUAL_PREFETCH == 1
            _mm_prefetch((const char*)(p_indice_nod1[DEG_1]), _MM_HINT_T0);
            _mm_prefetch((const char*)(&p_msg1r[DEG_1]), _MM_HINT_T0);
#endif
#endif

#pragma unroll(DEG_1)
            // Kernal 2 in algorithm
            for (j = 0; j < DEG_1; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1); // vNoeud:in algorithm:En
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r); // vMessg in algorithm: Lmn
                // The substraction operation does not achieve the maximum value just could achieve the minimum value
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var); // update Lnm
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8); // get the sign bit of the message
                sign = VECTOR_XOR(sign, cSign); // sign initial vector zero
                TYPE vAbs = VECTOR_ABS(vContr); // the magnitude Lnm the maximum is 8bit
                tab_vContr[j] = vContr;
                TYPE vTemp = min1; // min1 the mimimum value Lnm
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2); // the second minimum value of the Lnm
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

            TYPE norm_1 = VECTOR_SET2(factor_1);
            // 个人理解，先乘后除，为了提高精度。
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1); // a8 00 a9 00...a15 00
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1); // a1 00 a2 00...a7 00
                                                   //
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2); // shift 5 bits to right
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2); // saturation

            TYPE norm_2 = VECTOR_SET2(factor_2); // factor_1=factor_2 note that the NMS factor can be different
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1); // 16 bits

            cste_1 = VECTOR_MIN_1(cste_1, max_msg); // Lmn is 6bit value
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);

#if PETIT == 1
#if MANUAL_PREFETCH == 1
            for (j = 0; j < DEG_1; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_1]), _MM_HINT_T0);
#endif
#endif

#if (DEG_1 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8); // xor operation
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif
            // immediately update En the operations in Algorithm from line 18 to 21
#pragma unroll(DEG_1)
            for (j = 0; j < DEG_1; j++) {
                TYPE vContr = tab_vContr[j]; // Lnm in the ith iteration
                TYPE vAbs = VECTOR_ABS(vContr); // Lmn_new test if saturates
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8)); // sign is in Algorithm Line 15
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig); // compute the Lmn in the next iteration
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var); // update the En in Algorithm 21
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St); // Store Lmn
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr); // Store En
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
#if PETIT == 1
#if MANUAL_PREFETCH == 1
            _mm_prefetch((const char*)(*p_indice_nod2), _MM_HINT_T0);
#endif
#endif
        }

        ////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////

#if NB_DEGRES >= 2
        /********************************************DEG_2******************************************************************/
        for (int i = 0; i < DEG_2_COMPUTATIONS; i++) {

#if (DEG_2 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_2];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_2)
            for (int j = 0; j < DEG_2; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_2; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_2]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_2 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_2)
            for (int j = 0; j < DEG_2; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                                             // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 3
        /****************************************************DEG_3****************************************************************/
        for (int i = 0; i < DEG_3_COMPUTATIONS; i++) {

#if (DEG_3 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_3];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_3)
            for (int j = 0; j < DEG_3; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_3; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_3]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_3 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_3)
            for (int j = 0; j < DEG_3; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 4
        /*************************************DEG_4*************************************************/
        for (int i = 0; i < DEG_4_COMPUTATIONS; i++) {

#if (DEG_4 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_4];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_4)
            for (int j = 0; j < DEG_4; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_4; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_4]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_4 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_4)
            for (int j = 0; j < DEG_4; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 5
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_5_COMPUTATIONS; i++) {

#if (DEG_5 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_5];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_5)
            for (int j = 0; j < DEG_5; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_5; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_5]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_5 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_5)
            for (int j = 0; j < DEG_5; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 6
        /*************************************DEG_6*************************************************/
        for (int i = 0; i < DEG_6_COMPUTATIONS; i++) {

#if (DEG_6 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_6];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_6)
            for (int j = 0; j < DEG_6; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_6; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_6]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_6 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_6)
            for (int j = 0; j < DEG_6; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 7
        /*************************************DEG_7*************************************************/
        for (int i = 0; i < DEG_7_COMPUTATIONS; i++) {

#if (DEG_7 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_7];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_7)
            for (int j = 0; j < DEG_7; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_7; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_7]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_7 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_7)
            for (int j = 0; j < DEG_7; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 8
        /*************************************DEG_8*************************************************/
        for (int i = 0; i < DEG_8_COMPUTATIONS; i++) {

#if (DEG_8 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_8];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_8)
            for (int j = 0; j < DEG_8; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_8; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_8]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_8 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_8)
            for (int j = 0; j < DEG_8; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 9
        /*************************************DEG_9*************************************************/
        for (int i = 0; i < DEG_9_COMPUTATIONS; i++) {

#if (DEG_9 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_9];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_9)
            for (int j = 0; j < DEG_9; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_9; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_9]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_9 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_9)
            for (int j = 0; j < DEG_9; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 10
        /*************************************DEG_10*************************************************/
        for (int i = 0; i < DEG_10_COMPUTATIONS; i++) {

#if (DEG_10 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_10];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_10)
            for (int j = 0; j < DEG_10; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_10; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_10]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_10 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_10)
            for (int j = 0; j < DEG_10; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 11
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_11_COMPUTATIONS; i++) {

#if (DEG_11 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_11];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_11)
            for (int j = 0; j < DEG_11; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_11; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_11]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_11 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_11)
            for (int j = 0; j < DEG_11; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 12
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_12_COMPUTATIONS; i++) {

#if (DEG_12 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_12];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_12)
            for (int j = 0; j < DEG_12; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_12; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_12]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_12 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_12)
            for (int j = 0; j < DEG_12; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif

#if NB_DEGRES >= 13
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_13_COMPUTATIONS; i++) {

#if (DEG_13 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_13];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_13)
            for (int j = 0; j < DEG_13; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_13; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_13]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_13 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_13)
            for (int j = 0; j < DEG_13; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 14
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_14_COMPUTATIONS; i++) {

#if (DEG_14 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_14];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_14)
            for (int j = 0; j < DEG_14; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_14; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_14]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_14 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_14)
            for (int j = 0; j < DEG_14; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 15
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_15_COMPUTATIONS; i++) {

#if (DEG_15 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_15];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_15)
            for (int j = 0; j < DEG_15; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_15; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_15]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_15 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_15)
            for (int j = 0; j < DEG_15; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 16
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_16_COMPUTATIONS; i++) {

#if (DEG_16 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_16];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_16)
            for (int j = 0; j < DEG_16; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_16; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_16]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_16 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_16)
            for (int j = 0; j < DEG_16; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 17
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_17_COMPUTATIONS; i++) {

#if (DEG_17 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_17];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_17)
            for (int j = 0; j < DEG_17; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_17; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_17]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_17 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_17)
            for (int j = 0; j < DEG_17; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 18
        /*************************************DEG_18*************************************************/
        for (int i = 0; i < DEG_18_COMPUTATIONS; i++) {

#if (DEG_18 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_18];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_18)
            for (int j = 0; j < DEG_18; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_18; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_18]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_18 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_18)
            for (int j = 0; j < DEG_18; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 19
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_19_COMPUTATIONS; i++) {

#if (DEG_19 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_19];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_19)
            for (int j = 0; j < DEG_19; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_19; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_19]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_19 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_19)
            for (int j = 0; j < DEG_19; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 20
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_20_COMPUTATIONS; i++) {

#if (DEG_20 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_20];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_20)
            for (int j = 0; j < DEG_20; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_20; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_20]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_20 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_20)
            for (int j = 0; j < DEG_20; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
        ////shorten
        // for (int i = NmoinsK - _ShortenBits;i < NmoinsK;++i)
        //{
        //	ldpc->var_nodes[i] = maxMesg;
        //}
    }
    // DATE:20181028
    /*
                The Program is modified for LDPC codes without puncture and shorten and the decoder output
                is the whole codeword,and the statistic result of BER and FER is the whole codeword
                it is different from the Program for 5G platform
        */
    if ((NOEUD) % 32 == 0) //解交织
    {
        uchar_itranspose_avx((TYPE*)var_nodes, (TYPE*)decodedBits, (NOEUD));
    } else {
        char* ptr = (char*)var_nodes;
        for (int i = 0; i < (NOEUD); i += 1) {
            for (int j = 0; j < 32; j += 1) {
                decodedBits[j * (NOEUD) + i] = (ptr[32 * i + j] > 0);
            }
        }
    }
    ////DATE:20190306
    ///*
    // The program is modified for 5G LDPC codes, the decodedBits are the information bits
    // and the FER and BER are compared with the information bits
    //*/
    //
    // if ((NmoinsK - _ShortenBits) % 32 == 0)
    //{
    //	uchar_itranspose_avx((TYPE*)var_nodes, (TYPE*)decodedBits, (NmoinsK - _ShortenBits));
    //}
    // else
    //{
    //	char* ptr = (char*)var_nodes;
    //	for (int i = 0; i < (NmoinsK - _ShortenBits); i += 1)
    //	{
    //		for (int j = 0; j < 32; j += 1)
    //		{
    //			decodedBits[j *(NmoinsK - _ShortenBits) + i] = (ptr[32 * i + j] >
    // 0);//varnode0存第一帧第一个信息 varnode1存第二帧第一个信息 varnode32存第一帧第二个信息
    //		}
    //	}
    //}
    //
}
void CLDPC::Decode1() //固定迭代次数为60的译码过程
{
    Parameter_Simulation p_simulation;
    ReadProfile(&p_simulation);

    int8_t* ptr;
    int8_t* pp;
    int factor_1 = p_simulation.Factor_1;
    int factor_2 = p_simulation.Factor_2;
    const TYPE zero = VECTOR_ZERO;
    const TYPE minMesg = VECTOR_SET1(SAT_NEG_MSG);
    size_t i, j, z;
    /* Lmn Initialization */
    for (i = 0; i < MESSAGE; ++i) {
        var_msgs[i] = zero;
    }
    // The information part
    if ((NmoinsK - _PunctureBits - _ShortenBits) % 32 == 0) {
        uchar_transpose_avx(
            (TYPE*)fixInput, (TYPE*)(var_nodes + _PunctureBits), (NmoinsK - _PunctureBits - _ShortenBits));
    } else {
        ptr = (int8_t*)(var_nodes + _PunctureBits); //指针起始位置在puncture后，即打孔在最前
        for (i = 0; i < (NmoinsK - _PunctureBits - _ShortenBits); ++i) {
            for (int j = 0; j < 32; ++j) {
                ptr[i * 32 + j]
                    = fixInput[j * (NmoinsK - _PunctureBits - _ShortenBits) + i]; //除去puncture和shorten的原始似然比
            }
        }
    }
    // The check part
    pp = fixInput + (NmoinsK - _PunctureBits - _ShortenBits) * 32;
    if (_NoCheck % 32 == 0) {

        uchar_transpose_avx((TYPE*)(pp), (TYPE*)(var_nodes + NmoinsK), _NoCheck);
    } else {
        ptr = (int8_t*)(var_nodes + NmoinsK);
        for (i = 0; i < _NoCheck; ++i) {
            for (int j = 0; j < 32; ++j) {
                ptr[i * 32 + j] = pp[j * _NoCheck + i];
            }
        }
    }

    // puncture 最前 llr=0
    for (i = 0; i < _PunctureBits; ++i) {
        var_nodes[i] = zero;
    }
    // shorten 最后 llr=minMesg
    for (i = NmoinsK - _ShortenBits; i < NmoinsK; ++i) {
        var_nodes[i] = minMesg;
    }

    int nombre_iterations = nb_iteration;
    // Fixed iterations
    while (nombre_iterations--) { //一次迭代开始，没有动态终止
        TYPE* p_msg1r = var_msgs; // read Lmn
        TYPE* p_msg1w = var_msgs; // write Lmn
#if PETIT == 1
        TYPE** p_indice_nod1 = p_vn_adr; // p_vn_adr: Stores the En,read
        TYPE** p_indice_nod2 = p_vn_adr; // write

#else
        const unsigned short* p_indice_nod1 = PosNoeudsVariable;
        const unsigned short* p_indice_nod2 = PosNoeudsVariable;
#endif

        const TYPE min_var = VECTOR_SET1(vSAT_NEG_VAR); // semimimmum value 8bit -127
        const TYPE max_msg = VECTOR_SET1(vSAT_POS_MSG); // maximum value 6bit +31
        const TYPE max_var = VECTOR_SET1(vSAT_POS_VAR); // maximum value 8 bit 127
                                                        // DEG_1_COMPUTATIONS: The number of the rows of degree DEG_1
        /**************************************DEG_1****************************************************************/
        for (i = 0; i < DEG_1_COMPUTATIONS; i++) {
            TYPE tab_vContr[DEG_1]; // DEG1
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR); // 8bit maximum value  +127
            TYPE min2 = min1; // min2=min1=127

#if (DEG_1 & 0x01) == 1
            const unsigned char sign8 = 0x80; // sign8 =128=1000000B
            const unsigned char isign8
                = 0xC0; // isign =12*16=11000000B,之所以使用的是X100,0000,是为了如果他的符号为是0，则在进行
                        //_mm_sign_epi8运算时，符合要求。
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80;
            const unsigned char isign8b = 0x40; // isign8b=01000000B=64
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

#if PETIT == 1
#if MANUAL_PREFETCH == 1
            _mm_prefetch((const char*)(p_indice_nod1[DEG_1]), _MM_HINT_T0);
            _mm_prefetch((const char*)(&p_msg1r[DEG_1]), _MM_HINT_T0);
#endif
#endif

#pragma unroll(DEG_1)
            // Kernal 2 in algorithm
            for (int j = 0; j < DEG_1; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1); // vNoeud:in algorithm:En
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r); // vMessg in algorithm: Lmn
                                                    // The substraction operation does not achieve the maximum value
                                                    // just could achieve the minimum value
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var); // update Lnm
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8); // get the sign bit of the message
                sign = VECTOR_XOR(sign, cSign); // sign initial vector zero
                TYPE vAbs = VECTOR_ABS(vContr); // the magnitude Lnm the maximum is 8bit
                tab_vContr[j] = vContr;
                TYPE vTemp = min1; // min1 the mimimum value Lnm
                min1 = VECTOR_MIN_1(vAbs, min1);
                //   #define VECTOR_MIN_2(val,old_min1,min2) \
															                 //(VECTOR_MIN(min2,VECTOR_MAX(val,old_min1)))
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2); // the second minimum value of the Lnm
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

            TYPE norm_1 = VECTOR_SET2(factor_1);
            // 个人理解，先乘后除，为了提高精度。
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1); // a8 00 a9 00...a15 00
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1); // a1 00 a2 00...a7 00
                                                   //
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2); // shift 5 bits to right
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2); // saturation

            TYPE norm_2 = VECTOR_SET2(factor_2); // factor_1=factor_2 note that the NMS factor can be different
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1); // 16 bits

            cste_1 = VECTOR_MIN_1(cste_1, max_msg); // Lmn is 6bit value
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);

#if PETIT == 1
#if MANUAL_PREFETCH == 1
            for (int j = 0; j < DEG_1; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_1]), _MM_HINT_T0);
#endif
#endif

#if (DEG_1 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8); // xor operation
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif
            // immediately update En the operations in Algorithm from line 18 to 21
#pragma unroll(DEG_1)
            for (int j = 0; j < DEG_1; j++) {
                TYPE vContr = tab_vContr[j]; // Lnm in the ith iteration
                TYPE vAbs = VECTOR_ABS(vContr); // Lmn_new test if saturates
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8)); // sign is in Algorithm Line 15
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig); // compute the Lmn in the next iteration
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var); // update the En in Algorithm 21
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St); // Store Lmn
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr); // Store En
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
#if PETIT == 1
#if MANUAL_PREFETCH == 1
            _mm_prefetch((const char*)(*p_indice_nod2), _MM_HINT_T0);
#endif
#endif
        }

        ////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////

#if NB_DEGRES >= 2
        /********************************************DEG_2******************************************************************/
        for (int i = 0; i < DEG_2_COMPUTATIONS; i++) {

#if (DEG_2 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_2];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_2)
            for (int j = 0; j < DEG_2; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_2; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_2]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_2 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_2)
            for (int j = 0; j < DEG_2; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                                             // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 3
        /****************************************************DEG_3****************************************************************/
        for (int i = 0; i < DEG_3_COMPUTATIONS; i++) {

#if (DEG_3 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_3];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_3)
            for (int j = 0; j < DEG_3; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_3; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_3]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_3 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_3)
            for (int j = 0; j < DEG_3; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 4
        /*************************************DEG_4*************************************************/
        for (int i = 0; i < DEG_4_COMPUTATIONS; i++) {

#if (DEG_4 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_4];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_4)
            for (int j = 0; j < DEG_4; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_4; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_4]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_4 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_4)
            for (int j = 0; j < DEG_4; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 5
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_5_COMPUTATIONS; i++) {

#if (DEG_5 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_5];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_5)
            for (int j = 0; j < DEG_5; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_5; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_5]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_5 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_5)
            for (int j = 0; j < DEG_5; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 6
        /*************************************DEG_6*************************************************/
        for (int i = 0; i < DEG_6_COMPUTATIONS; i++) {

#if (DEG_6 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_6];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_6)
            for (int j = 0; j < DEG_6; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_6; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_6]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_6 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_6)
            for (int j = 0; j < DEG_6; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 7
        /*************************************DEG_7*************************************************/
        for (int i = 0; i < DEG_7_COMPUTATIONS; i++) {

#if (DEG_7 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_7];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_7)
            for (int j = 0; j < DEG_7; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_7; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_7]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_7 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_7)
            for (int j = 0; j < DEG_7; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 8
        /*************************************DEG_8*************************************************/
        for (int i = 0; i < DEG_8_COMPUTATIONS; i++) {

#if (DEG_8 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_8];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_8)
            for (int j = 0; j < DEG_8; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_8; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_8]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_8 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_8)
            for (int j = 0; j < DEG_8; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 9
        /*************************************DEG_9*************************************************/
        for (int i = 0; i < DEG_9_COMPUTATIONS; i++) {

#if (DEG_9 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_9];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_9)
            for (int j = 0; j < DEG_9; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_9; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_9]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_9 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_9)
            for (int j = 0; j < DEG_9; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 10
        /*************************************DEG_10*************************************************/
        for (int i = 0; i < DEG_10_COMPUTATIONS; i++) {

#if (DEG_10 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_10];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_10)
            for (int j = 0; j < DEG_10; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_10; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_10]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_10 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_10)
            for (int j = 0; j < DEG_10; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 11
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_11_COMPUTATIONS; i++) {

#if (DEG_11 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_11];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_11)
            for (int j = 0; j < DEG_11; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_11; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_11]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_11 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_11)
            for (int j = 0; j < DEG_11; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 12
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_12_COMPUTATIONS; i++) {

#if (DEG_12 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_12];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_12)
            for (int j = 0; j < DEG_12; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_12; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_12]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_12 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_12)
            for (int j = 0; j < DEG_12; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif

#if NB_DEGRES >= 13
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_13_COMPUTATIONS; i++) {

#if (DEG_13 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_13];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_13)
            for (int j = 0; j < DEG_13; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_13; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_13]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_13 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_13)
            for (int j = 0; j < DEG_13; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 14
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_14_COMPUTATIONS; i++) {

#if (DEG_14 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_14];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_14)
            for (int j = 0; j < DEG_14; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_14; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_14]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_14 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_14)
            for (int j = 0; j < DEG_14; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 15
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_15_COMPUTATIONS; i++) {

#if (DEG_15 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_15];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_15)
            for (int j = 0; j < DEG_15; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_15; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_15]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_15 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_15)
            for (int j = 0; j < DEG_15; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 16
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_16_COMPUTATIONS; i++) {

#if (DEG_16 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_16];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_16)
            for (int j = 0; j < DEG_16; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_16; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_16]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_16 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_16)
            for (int j = 0; j < DEG_16; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 17
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_17_COMPUTATIONS; i++) {

#if (DEG_17 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_17];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_17)
            for (int j = 0; j < DEG_17; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_17; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_17]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_17 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_17)
            for (int j = 0; j < DEG_17; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 18
        /*************************************DEG_18*************************************************/
        for (int i = 0; i < DEG_18_COMPUTATIONS; i++) {

#if (DEG_18 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_18];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_18)
            for (int j = 0; j < DEG_18; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_18; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_18]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_18 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_18)
            for (int j = 0; j < DEG_18; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 19
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_19_COMPUTATIONS; i++) {

#if (DEG_19 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_19];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_19)
            for (int j = 0; j < DEG_19; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_19; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_19]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_19 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_19)
            for (int j = 0; j < DEG_19; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
#if NB_DEGRES >= 20
        /*************************************DEG_5*************************************************/
        for (int i = 0; i < DEG_20_COMPUTATIONS; i++) {

#if (DEG_20 & 0x01) == 1
            const unsigned char sign8 = 0x80;
            const unsigned char isign8 = 0xC0;
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8 = VECTOR_SET1(isign8);
#else
            const unsigned char sign8 = 0x80; // 1000 0000
            const unsigned char isign8b = 0x40; // 0100 0000
            const TYPE msign8 = VECTOR_SET1(sign8);
            const TYPE misign8b = VECTOR_SET1(isign8b);
#endif

            TYPE tab_vContr[DEG_20];
            TYPE sign = VECTOR_ZERO;
            TYPE min1 = VECTOR_SET1(vSAT_POS_VAR);
            TYPE min2 = min1;

#pragma unroll(DEG_20)
            for (int j = 0; j < DEG_20; j++) {
#if PETIT == 1
                TYPE vNoeud = VECTOR_LOAD(*p_indice_nod1);
#else
                TYPE vNoeud = VECTOR_LOAD(&var_nodes[(*p_indice_nod1)]);
#endif
                TYPE vMessg = VECTOR_LOAD(p_msg1r);
                TYPE vContr = VECTOR_SUB_AND_SATURATE_VAR_8bits(vNoeud, vMessg, min_var);
                TYPE cSign = VECTOR_GET_SIGN_BIT(vContr, msign8);
                sign = VECTOR_XOR(sign, cSign);
                TYPE vAbs = VECTOR_ABS(vContr); // modify
                tab_vContr[j] = vContr;
                TYPE vTemp = min1;
                min1 = VECTOR_MIN_1(vAbs, min1);
                min2 = VECTOR_MIN_2(vAbs, vTemp, min2);
                p_indice_nod1 += 1;
                p_msg1r += 1;
            }

#if PETIT == 1
            for (int j = 0; j < DEG_20; j++) {
                _mm_prefetch((const char*)(p_indice_nod1[j]), _MM_HINT_T0);
            }
            _mm_prefetch((const char*)(p_indice_nod1[DEG_20]), _MM_HINT_T0);
#endif

            TYPE norm_1 = VECTOR_SET2(factor_1);
            TYPE h_min1 = VECTOR_UNPACK_HIGH(min1);
            TYPE l_min1 = VECTOR_UNPACK_LOW(min1);
            TYPE h_cste_2 = VECTOR_MUL(h_min1, norm_1);
            h_cste_2 = VECTOR_DIV32(h_cste_2);
            TYPE l_cste_2 = VECTOR_MUL(l_min1, norm_1);
            l_cste_2 = VECTOR_DIV32(l_cste_2);
            TYPE cste_2 = VECTOR_PACK(h_cste_2, l_cste_2);

            TYPE norm_2 = VECTOR_SET2(factor_2);
            TYPE h_min2 = VECTOR_UNPACK_HIGH(min2);
            TYPE l_min2 = VECTOR_UNPACK_LOW(min2);
            TYPE h_cste_1 = VECTOR_MUL(h_min2, norm_2);
            h_cste_1 = VECTOR_DIV32(h_cste_1);
            TYPE l_cste_1 = VECTOR_MUL(l_min2, norm_2);
            l_cste_1 = VECTOR_DIV32(l_cste_1);
            TYPE cste_1 = VECTOR_PACK(h_cste_1, l_cste_1);
            cste_1 = VECTOR_MIN_1(cste_1, max_msg);
            cste_2 = VECTOR_MIN_1(cste_2, max_msg);
#if (DEG_20 & 0x01) == 1
            sign = VECTOR_XOR(sign, misign8);
#else
            sign = VECTOR_XOR(sign, misign8b);
#endif

#pragma unroll(DEG_20)
            for (int j = 0; j < DEG_20; j++) {
                TYPE vContr = tab_vContr[j];
                TYPE vAbs = VECTOR_ABS(vContr);
                // TYPE vRes = VECTOR_CMOV(vAbs, min1, cste_1, cste_2);
                TYPE z = VECTOR_EQUAL(vAbs, min1); // a==min1?== 0xff else 0x00
                TYPE g = VECTOR_AND(cste_1, z); //
                TYPE h = VECTOR_ANDNOT(z, cste_2); // cste_1 the second minimum value
                TYPE vRes = VECTOR_OR(g, h); //
                TYPE vSig = VECTOR_XOR(sign, VECTOR_GET_SIGN_BIT(vContr, msign8));
                TYPE v2St = VECTOR_invSIGN2(vRes, vSig);
                TYPE v2Sr = VECTOR_ADD_AND_SATURATE_VAR_8bits(vContr, v2St, min_var);
                v2Sr = VECTOR_MIN(v2Sr, max_var);
                VECTOR_STORE(p_msg1w, v2St);
#if PETIT == 1
                VECTOR_STORE(*p_indice_nod2, v2Sr);
#else
                VECTOR_STORE(&var_nodes[(*p_indice_nod2)], v2Sr);
#endif
                p_msg1w += 1;
                p_indice_nod2 += 1;
            }
        }
#endif
        ////shorten
        // for (int i = NmoinsK - _ShortenBits;i < NmoinsK;++i)
        //{
        //	ldpc->var_nodes[i] = maxMesg;
        //}
    }
    // DATE:20181028
    /*
                The Program is modified for LDPC codes without puncture and shorten and the decoder output
                is the whole codeword,and the statistic result of BER and FER is the whole codeword
                it is different from the Program for 5G platform
        */
    if ((NOEUD) % 32 == 0) {
        uchar_itranspose_avx((TYPE*)var_nodes, (TYPE*)decodedBits, (NOEUD));
    } else {
        char* ptr = (char*)var_nodes;
        for (int i = 0; i < (NOEUD); i += 1) {
            for (int j = 0; j < 32; j += 1) {
                decodedBits[j * (NOEUD) + i] = (ptr[32 * i + j] > 0);
            }
        }
    }

    ////DATE:20190306
    ///*
    // The program is modified for 5G LDPC codes, the decodedBits are the information bits
    // and the FER and BER are compared with the information bits
    //*/
    // if ((NmoinsK - _ShortenBits) % 32 == 0)
    //{
    //	uchar_itranspose_avx((TYPE*)var_nodes, (TYPE*)decodedBits, (NmoinsK - _ShortenBits));
    //}
    // else
    //{
    //	char* ptr = (char*)var_nodes;
    //	for (int i = 0; i < (NmoinsK - _ShortenBits); i += 1)
    //	{
    //		for (int j = 0; j < 32; j += 1)
    //		{
    //			decodedBits[j *(NmoinsK - _ShortenBits) + i] = (ptr[32 * i + j] > 0);
    //		}
    //	}
    //}
}

void CLDPC::float2LimitChar_6bit(int8_t* output, const float* input, float scale,
    int length) // float2LimitChar_8bit(去除shorten和puncture原始似然比,
                // channel->BPSKSymbol, scale, BitsOverChannel * 32);
// change dataform from float to int8_t

{
    int qbits = 2;
    int limit = 31;

    int i;
    const __m128i limit0_128i = _mm_set1_epi8(31);
    const __m128i limit1_128i = _mm_set1_epi8(-31);
    // float	scale =5.8;
    // (float)(1 << qbits);
    __m256 scale_256;
    __m256 llra_in_256, llrb_in_256;
    __m256i llra_in_256i, llrb_in_256i;
    __m128i* output_128i = (__m128i*)output; //强制类型转换
    __m128i llra_in_128i, llrb_in_128i, llrab_in_128i;

    int8_t output_tmp;
    // float absmax_float;

    // CChannel *channel;
    // int Blength = 32 * BitsOverChannel;
    // for (int j = 0; j <Blength; j++)
    //{
    //	cout << "in floattochar" << (float)channel->BPSKSymbol[j] << endl;
    //	cout << "pressenter" << endl;
    //	getchar();
    //}

    for (i = 0; i + 15 < length; i += 16) // i=0 16 32……
    {

        llra_in_256 = _mm256_loadu_ps(
            input + i); // float是32位比特，这个是两个并行的计算 a存8float=8*32位=8*（32/8）字节=32字节
        // llra -0.79 0.61 -1.7 0.39 0.97 0.46 0.21 0.95
        llrb_in_256 = _mm256_loadu_ps(input + i + 8); // b存8float
        // llrb 0.93 0.26 0.66 -0.52 1.48 -0.8 0.56 -0.35
        // AVX版本
        // load 16 float data

        // absmax_float = absmax(llra_in_256, llrb_in_256);
        scale_256 = _mm256_set1_ps(scale); // scale 是一个标量，用来放大信号的 存了8个float 值为32

        llra_in_256 = _mm256_mul_ps(llra_in_256, scale_256); // llra 8个float 与scale_256 8个float执行SIMD乘法。
        // llra -25.43 19.54 -57.50
        llrb_in_256 = _mm256_mul_ps(llrb_in_256, scale_256); //一个8

        //执行乘法之后0.5 -0.5都会被量化为0 四舍五入
        llra_in_256i = _mm256_cvtps_epi32(llra_in_256); // float to int32 把8个float转为8个int32

        // m256i_i32 -25 19 -57
        llrb_in_256i = _mm256_cvtps_epi32(llrb_in_256);
        //_mm256_cvttps_epi32

        /*__m128i _mm_packs_epi32 (__m128i a, __m128i b);将a和b中的8个int32打包成8个int16 a0 a1 a2 a3 b0 b1 b2 b3*/
        llra_in_128i = _mm_packs_epi32(
            _mm256_extractf128_si256(llra_in_256i, 0), _mm256_extractf128_si256(llra_in_256i, 1)); // in32 to int16
        // llra m256i_i16 -25 19 -57
        llrb_in_128i
            = _mm_packs_epi32(_mm256_extractf128_si256(llrb_in_256i, 0), _mm256_extractf128_si256(llrb_in_256i, 1));
        // llrb m256i_i16 29 8 21
        llrab_in_128i = _mm_packs_epi16(llra_in_128i, llrb_in_128i); // int16 to int8（超过8位=255，否则本身）
        // llrab m128i_i8 -25 19 -57……29 8 21……一共16个数

        llrab_in_128i = _mm_min_epi8(llrab_in_128i, limit0_128i);
        llrab_in_128i = _mm_max_epi8(llrab_in_128i, limit1_128i); // set limitation(6bit)

        _mm_storeu_si128(output_128i++, llrab_in_128i);
    }
    for (; i < length; i++) {
        output_tmp = (int8_t)(*(input + i) * scale);
        output_tmp = (output_tmp > limit) ? limit : output_tmp;
        *(output + i) = (output_tmp < -limit) ? -limit : output_tmp; // fractional data
    }
}

void CLDPC::float2LimitChar_5bit(int8_t* output, const float* input, float scale,
    int length) // float2LimitChar_8bit(去除shorten和puncture原始似然比, channel->BPSKSymbol, scale, BitsOverChannel *
                // 32); change dataform from float to int8_t

{
    int qbits = 2;
    int limit = 31;
    int i;
    const __m128i limit0_128i = _mm_set1_epi8(15);
    const __m128i limit1_128i = _mm_set1_epi8(-16);
    // float	scale =5.8;
    // (float)(1 << qbits);
    __m256 scale_256;
    __m256 llra_in_256, llrb_in_256;
    __m256i llra_in_256i, llrb_in_256i;
    __m128i* output_128i = (__m128i*)output; //强制类型转换
    __m128i llra_in_128i, llrb_in_128i, llrab_in_128i;

    int8_t output_tmp;
    // float absmax_float;
    for (i = 0; i + 15 < length; i += 16) // i=0 16 32……
    {
        llra_in_256 = _mm256_loadu_ps(
            input + i); // float是32位比特，这个是两个并行的计算 a存8float=8*32位=8*（32/8）字节=32字节

        llrb_in_256 = _mm256_loadu_ps(input + i + 8); // b存8float
                                                      // AVX版本
                                                      // load 16 float data

        scale_256 = _mm256_set1_ps(scale); // scale 是一个标量，用来放大信号的 存了8个float 值为32

        llra_in_256 = _mm256_mul_ps(llra_in_256, scale_256); // llra 8个float 与scale_256 8个float执行SIMD乘法。

        llrb_in_256 = _mm256_mul_ps(llrb_in_256, scale_256); //一个8

        llra_in_256i = _mm256_cvttps_epi32(llra_in_256); // float to int32 把8个float转为8个int32

        llrb_in_256i = _mm256_cvttps_epi32(llrb_in_256);

        /*__m128i _mm_packs_epi32 (__m128i a, __m128i b);将a和b中的8个int32打包成8个int16 a0 a1 a2 a3 b0 b1 b2 b3*/
        llra_in_128i = _mm_packs_epi32(
            _mm256_extractf128_si256(llra_in_256i, 0), _mm256_extractf128_si256(llra_in_256i, 1)); // in32 to int16

        llrb_in_128i
            = _mm_packs_epi32(_mm256_extractf128_si256(llrb_in_256i, 0), _mm256_extractf128_si256(llrb_in_256i, 1));

        llrab_in_128i = _mm_packs_epi16(llra_in_128i, llrb_in_128i); // int16 to int8（超过8位=255，否则本身）

        llrab_in_128i = _mm_min_epi8(llrab_in_128i, limit0_128i);
        llrab_in_128i = _mm_max_epi8(llrab_in_128i, limit1_128i); // set limitation(5bit)

        _mm_storeu_si128(output_128i++, llrab_in_128i);
    }
    for (; i < length; i++) {
        output_tmp = (int8_t)(*(input + i) * scale);
        output_tmp = (output_tmp > limit) ? limit : output_tmp;
        *(output + i) = (output_tmp < -limit) ? -limit : output_tmp; // fractional data
    }
}

void CLDPC::float2LimitChar_4bit(int8_t* output, const float* input, float scale,
    int length) // float2LimitChar_8bit(去除shorten和puncture原始似然比, channel->BPSKSymbol, scale, BitsOverChannel *
                // 32); change dataform from float to int8_t

{
    int qbits = 2;
    int limit = 7;
    int i;
    const __m128i limit0_128i = _mm_set1_epi8(7);
    const __m128i limit1_128i = _mm_set1_epi8(-7); // -7 会影响性能！
    // float	scale =5.8;
    // (float)(1 << qbits);
    __m256 scale_256;
    __m256 llra_in_256, llrb_in_256;
    __m256i llra_in_256i, llrb_in_256i;
    __m128i* output_128i = (__m128i*)output; //强制类型转换
    __m128i llra_in_128i, llrb_in_128i, llrab_in_128i;

    int8_t output_tmp;
    // float absmax_float;
    for (i = 0; i + 15 < length; i += 16) // i=0 16 32……
    {
        llra_in_256 = _mm256_loadu_ps(
            input + i); // float是32位比特，这个是两个并行的计算 a存8float=8*32位=8*（32/8）字节=32字节

        llrb_in_256 = _mm256_loadu_ps(input + i + 8); // b存8float
                                                      // AVX版本
                                                      // load 16 float data

        scale_256 = _mm256_set1_ps(scale); // scale 是一个标量，用来放大信号的 存了8个float 值为32

        llra_in_256 = _mm256_mul_ps(llra_in_256, scale_256); // llra 8个float 与scale_256 8个float执行SIMD乘法。

        llrb_in_256 = _mm256_mul_ps(llrb_in_256, scale_256); //一个8

        llra_in_256i = _mm256_cvttps_epi32(llra_in_256); // float to int32 把8个float转为8个int32

        llrb_in_256i = _mm256_cvttps_epi32(llrb_in_256);

        /*__m128i _mm_packs_epi32 (__m128i a, __m128i b);将a和b中的8个int32打包成8个int16 a0 a1 a2 a3 b0 b1 b2 b3*/
        llra_in_128i = _mm_packs_epi32(
            _mm256_extractf128_si256(llra_in_256i, 0), _mm256_extractf128_si256(llra_in_256i, 1)); // in32 to int16

        llrb_in_128i
            = _mm_packs_epi32(_mm256_extractf128_si256(llrb_in_256i, 0), _mm256_extractf128_si256(llrb_in_256i, 1));

        llrab_in_128i = _mm_packs_epi16(llra_in_128i, llrb_in_128i); // int16 to int8（超过8位=255，否则本身）

        llrab_in_128i = _mm_min_epi8(llrab_in_128i, limit0_128i);
        llrab_in_128i = _mm_max_epi8(llrab_in_128i, limit1_128i); // set limitation(5bit)

        _mm_storeu_si128(output_128i++, llrab_in_128i);
    }
    for (; i < length; i++) {
        output_tmp = (int8_t)(*(input + i) * scale);
        output_tmp = (output_tmp > limit) ? limit : output_tmp;
        *(output + i) = (output_tmp < -limit) ? -limit : output_tmp; // fractional data
    }
}

void CLDPC::float2LimitChar_3bit(int8_t* output, const float* input, float scale,
    int length) // float2LimitChar_8bit(去除shorten和puncture原始似然比, channel->BPSKSymbol, scale, BitsOverChannel *
                // 32); change dataform from float to int8_t

{
    int qbits = 2;
    int limit = 31;
    int i;
    const __m128i limit0_128i = _mm_set1_epi8(3);
    const __m128i limit1_128i = _mm_set1_epi8(-4);
    // float	scale =5.8;
    // (float)(1 << qbits);
    __m256 scale_256;
    __m256 llra_in_256, llrb_in_256;
    __m256i llra_in_256i, llrb_in_256i;
    __m128i* output_128i = (__m128i*)output; //强制类型转换
    __m128i llra_in_128i, llrb_in_128i, llrab_in_128i;

    int8_t output_tmp;
    // float absmax_float;
    for (i = 0; i + 15 < length; i += 16) // i=0 16 32……
    {
        llra_in_256 = _mm256_loadu_ps(
            input + i); // float是32位比特，这个是两个并行的计算 a存8float=8*32位=8*（32/8）字节=32字节

        llrb_in_256 = _mm256_loadu_ps(input + i + 8); // b存8float
                                                      // AVX版本
                                                      // load 16 float data

        scale_256 = _mm256_set1_ps(scale); // scale 是一个标量，用来放大信号的 存了8个float 值为32

        llra_in_256 = _mm256_mul_ps(llra_in_256, scale_256); // llra 8个float 与scale_256 8个float执行SIMD乘法。

        llrb_in_256 = _mm256_mul_ps(llrb_in_256, scale_256); //一个8

        llra_in_256i = _mm256_cvttps_epi32(llra_in_256); // float to int32 把8个float转为8个int32

        llrb_in_256i = _mm256_cvttps_epi32(llrb_in_256);

        /*__m128i _mm_packs_epi32 (__m128i a, __m128i b);将a和b中的8个int32打包成8个int16 a0 a1 a2 a3 b0 b1 b2 b3*/
        llra_in_128i = _mm_packs_epi32(
            _mm256_extractf128_si256(llra_in_256i, 0), _mm256_extractf128_si256(llra_in_256i, 1)); // in32 to int16

        llrb_in_128i
            = _mm_packs_epi32(_mm256_extractf128_si256(llrb_in_256i, 0), _mm256_extractf128_si256(llrb_in_256i, 1));

        llrab_in_128i = _mm_packs_epi16(llra_in_128i, llrb_in_128i); // int16 to int8（超过8位=255，否则本身）

        llrab_in_128i = _mm_min_epi8(llrab_in_128i, limit0_128i);
        llrab_in_128i = _mm_max_epi8(llrab_in_128i, limit1_128i); // set limitation(5bit)

        _mm_storeu_si128(output_128i++, llrab_in_128i);
    }
    for (; i < length; i++) {
        output_tmp = (int8_t)(*(input + i) * scale);
        output_tmp = (output_tmp > limit) ? limit : output_tmp;
        *(output + i) = (output_tmp < -limit) ? -limit : output_tmp; // fractional data
    }
}

void CLDPC::float2LimitChar_2bit(int8_t* output, const float* input, float scale, int length)
// change dataform from float to int8_t

{
    int qbits = 2;
    int limit = 31;
    int i;
    const __m128i limit0_128i = _mm_set1_epi8(1);
    const __m128i limit1_128i = _mm_set1_epi8(-2);
    const __m128i zero = _mm_set1_epi8(0);
    const __m128i bit_1_xor = _mm_set1_epi8(63);
    // float	scale =5.8;
    // (float)(1 << qbits);
    __m256 scale_256;
    __m256 llra_in_256, llrb_in_256;
    __m256i llra_in_256i, llrb_in_256i;
    __m128i* output_128i = (__m128i*)output; //强制类型转换
    __m128i llra_in_128i, llrb_in_128i, llrab_in_128i;

    int8_t output_tmp;
    // float absmax_float;
    for (i = 0; i + 15 < length; i += 16) {
        llra_in_256 = _mm256_loadu_ps(input + i); // float是32位比特，这个是两个并行的计算
        llrb_in_256 = _mm256_loadu_ps(input + i + 8);
        // AVX版本
        // load 16 float data

        // absmax_float = absmax(llra_in_256, llrb_in_256);
        scale_256 = _mm256_set1_ps(scale); // scale 是一个标量，用来放大信号的。
        llra_in_256 = _mm256_mul_ps(llra_in_256, scale_256);
        llrb_in_256 = _mm256_mul_ps(llrb_in_256, scale_256); //一个8

        llra_in_256i = _mm256_cvttps_epi32(llra_in_256); // float to int32
        llrb_in_256i = _mm256_cvttps_epi32(llrb_in_256);

        llra_in_128i = _mm_packs_epi32(
            _mm256_extractf128_si256(llra_in_256i, 0), _mm256_extractf128_si256(llra_in_256i, 1)); // in32 to int16
        llrb_in_128i
            = _mm_packs_epi32(_mm256_extractf128_si256(llrb_in_256i, 0), _mm256_extractf128_si256(llrb_in_256i, 1));

        llrab_in_128i = _mm_packs_epi16(llra_in_128i, llrb_in_128i); // int16 to int8

        // llrab_in_128i = _mm_cmpgt_epi8(zero,
        // llrab_in_128i);//分别比较a的每个8bits整数是否大于b的对应位置的8bits整数，若大于，则返回-1(ff)，否则返回0。
        // llrab_in_128i = _mm_xor_si128(llrab_in_128i, bit_1_xor);//00000000xor00111111=00111111(63)
        // 11111111xor00111111=11000000(-64)

        llrab_in_128i = _mm_min_epi8(llrab_in_128i, limit0_128i); //正数符号位为0 01111111
        llrab_in_128i = _mm_max_epi8(llrab_in_128i, limit1_128i); // set limitation +-31

        _mm_storeu_si128(output_128i++, llrab_in_128i);
    }
    for (; i < length; i++) {
        output_tmp = (int8_t)(*(input + i) * scale);
        output_tmp = (output_tmp >= 0) ? limit : -limit;
        *(output + i) = output_tmp; // fractional data
    }
}

void CLDPC::float2LimitChar_1bit(int8_t* output, const float* input, float scale, int length)
// change dataform from float to int8_t

{
    int qbits = 2;
    int limit = 31;
    int i;
    const __m128i limit0_128i = _mm_set1_epi8(31);
    const __m128i limit1_128i = _mm_set1_epi8(-31);
    const __m128i zero = _mm_set1_epi8(0);
    const __m128i bit_1_xor = _mm_set1_epi8(-64);
    // const __m128i bit_1_xor = _mm_set1_epi8(63);
    // float	scale =5.8;
    // (float)(1 << qbits);
    __m256 scale_256; // 1bit需要非常大的scale值
    __m256 llra_in_256, llrb_in_256;
    __m256i llra_in_256i, llrb_in_256i;
    __m128i* output_128i = (__m128i*)output; //强制类型转换
    __m128i llra_in_128i, llrb_in_128i, llrab_in_128i;

    int8_t output_tmp;
    // float absmax_float;
    for (i = 0; i + 15 < length; i += 16) {
        llra_in_256 = _mm256_loadu_ps(input + i); // float是32位比特，这个是两个并行的计算
        llrb_in_256 = _mm256_loadu_ps(input + i + 8);
        // AVX版本
        // load 16 float data

        // absmax_float = absmax(llra_in_256, llrb_in_256);
        scale_256 = _mm256_set1_ps(scale); // scale 是一个标量，用来放大信号的。
        llra_in_256 = _mm256_mul_ps(llra_in_256, scale_256);
        llrb_in_256 = _mm256_mul_ps(llrb_in_256, scale_256); //一个8

        llra_in_256i = _mm256_cvttps_epi32(llra_in_256); // float to int32
        llrb_in_256i = _mm256_cvttps_epi32(llrb_in_256);

        llra_in_128i = _mm_packs_epi32(
            _mm256_extractf128_si256(llra_in_256i, 0), _mm256_extractf128_si256(llra_in_256i, 1)); // in32 to int16
        llrb_in_128i
            = _mm_packs_epi32(_mm256_extractf128_si256(llrb_in_256i, 0), _mm256_extractf128_si256(llrb_in_256i, 1));

        llrab_in_128i = _mm_packs_epi16(llra_in_128i, llrb_in_128i); // int16 to int8
        /*>=0=+31*/
        llrab_in_128i = _mm_cmpgt_epi8(llrab_in_128i,
            zero); //分别比较a的每个8bits整数是否大于b的对应位置的8bits整数，若大于，则返回-1(ff)，否则返回0。
        llrab_in_128i = _mm_xor_si128(
            llrab_in_128i, bit_1_xor); // 00000000xor11000000=1100000(-64) 11111111xor11000000=00111111(63)

        /*>0=+31*/
        // llrab_in_128i =
        // _mm_cmpgt_epi8(zero,llrab_in_128i);//分别比较a的每个8bits整数是否大于b的对应位置的8bits整数，若大于，则返回-1(ff)，否则返回0。
        // llrab_in_128i = _mm_xor_si128(llrab_in_128i, bit_1_xor);//00000000xor11000000=00111111(63)
        // 11111111xor00111111=11000000(63)

        llrab_in_128i = _mm_min_epi8(llrab_in_128i, limit0_128i); //正数符号位为0 01111111
        llrab_in_128i = _mm_max_epi8(llrab_in_128i, limit1_128i); // set limitation +-31
        //或者 llrab_in_128i = _mm_min_epi8(llrab_in_128i, 0);
        // llrab_in_128i = _mm_max_epi8(llrab_in_128i, -1);//即>0=00000000 <0=11111111
        // llrab_in_128i = _mm_xor_si128(llrab_in_128i, 31);//00000000xor00011111=00011111(31)
        // 11111111xor00011111=11100000(-32)//负数值取反加一为绝对值
        _mm_storeu_si128(output_128i++, llrab_in_128i);
    }
    for (; i < length; i++) {
        output_tmp = (int8_t)(*(input + i) * scale);
        output_tmp = (output_tmp >= 0) ? limit : -limit;
        *(output + i) = output_tmp; // fractional data
    }
}

void CLDPC::Initial(int nb_frame, int MaxItertion)
{
    m_M = _NoCheck;
    m_K = NmoinsK;
    m_N = _NoVar;
    m_PunLen = _PunctureBits;
    m_ShortenLen = _ShortenBits;
    // m_Rate = (float)(float)(_NoVar - _NoCheck - _ShortenBits) / (BitsOverChannel);
    m_Rate = 0.8444444;
    m_frame = nb_frame; // default 16
    nb_iteration = MaxItertion;

    inputBits = (int8_t*)vec_malloc(sizeof(int8_t) * (NmoinsK - _ShortenBits) * m_frame); // 32帧的全都要存
    outputBits = (int8_t*)vec_malloc(sizeof(int8_t) * BitsOverChannel * m_frame);
    ; // new int8_t[BitsOverChannel*m_frame];

    encoder_varnodes = (__m256i*)vec_malloc(sizeof(__m256i) * m_N); // new __m256i[m_N];
    var_nodes = (__m256i*)vec_malloc(sizeof(__m256i) * m_N);
    ; // new __m256i[m_N];只存一帧的

    var_msgs = (__m256i*)vec_malloc(sizeof(__m256i) * _NoOnes); // new __m256i[_NoOnes];
    p_vn_adr = (__m256i**)vec_malloc(sizeof(__m256i*) * _NoOnes); // new TYPE *[_NoOnes];
                                                                  // DATE:20190508
                                                                  // /
    decodedBits
        = (int8_t*)vec_malloc(sizeof(int8_t) * (NOEUD)*m_frame); //存所有的东西，包含punture information shorten check
    // decodedBits, the length is the length of the information bits
    // date:20190306
    // decodedBits = (int8_t*)vec_malloc(sizeof(int8_t)*(NmoinsK - _ShortenBits)*m_frame);//new int8_t[(NmoinsK -
    // _ShortenBits)*m_frame];
    fixInput = (int8_t*)vec_malloc(sizeof(int8_t) * BitsOverChannel * m_frame); // new int8_t[BitsOverChannel*m_frame];
    VN_weight_ = (int8_t*)malloc(sizeof(int8_t) * _NoVar);
    VN_weight_count();

    errorbitblock = (int*)malloc(sizeof(int) * BitsOverChannel * m_frame); //记录bit错误的块
    errorbitindex = (int*)malloc(sizeof(int) * BitsOverChannel * m_frame); //记录bit错误所在的块的位置
    errorcheckblock = (int*)malloc(sizeof(int) * BitsOverChannel * m_frame); //记录check错误的块
    errorcheckindex = (int*)malloc(sizeof(int) * BitsOverChannel * m_frame); //记录check错误所在的块的位置
    errorfloat = (float*)malloc(sizeof(float) * NOEUD * m_frame); //记录对应过信道后BPSKSymbol的float值
    errorchar = (int*)malloc(sizeof(int) * NOEUD * m_frame); //记录对应量化后FixInput的值

    int i = 0;
    for (i = 0; i < _NoOnes; ++i) {
        p_vn_adr[i] = &var_nodes[PosNoeudsVariable[i]];
    }
}

Statistic CLDPC::CalculateErrors(float* bpskinput, int8_t* charinput, int collectflag)
{
    Statistic Test;
    CLDPC* ldpc;
    CChannel* channel;
    Parameter_Simulation p_simulation;
    ofstream eout;
    ofstream nout;
    ofstream dout;
    ReadProfile(&p_simulation);
    Test.ErrorBits = 0;
    Test.ErrorFrame = 0;
    Test.LT3ErrBitFrame = 0;
    unsigned long errorBits = 0;
    unsigned long errorCheck = 0;
    int collect_err = collectflag;
    int Z = p_simulation.Z;
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int m = 0;

    for (i = 0; i < m_frame; ++i) // 32帧串行比较
    {
        errorBits = 0;
        m = 0;
        //统计BER FER 仅统计puncture+information 输出biterror(k)
        for (j = 0; j < ((NmoinsK - _ShortenBits)); ++j) {
            // decodedbits中存储了所有块的信息，但只比较信息比特，毕竟 inputBits 里只有信息比特
            if (decodedBits[i * (NOEUD) + j] != inputBits[i * (NmoinsK - _ShortenBits) + j]) {
                errorBits++;
                int bitblocknum = j / Z;
                errorbitblock[k] = bitblocknum + 1;
                errorbitindex[k] = j - Z * bitblocknum;
                k++;
            }
        }
        //统计check错误 输出checkerror(m)
        for (j = NmoinsK; j < _NoVar; ++j) //统计check
        {
            if (decodedBits[i * (NOEUD) + j]
                != outputBits[32 * (NmoinsK - _PunctureBits - _ShortenBits) + i * _NoCheck + (j - NmoinsK)]) // check
            {
                int checkblocknum = j / Z;
                errorcheckblock[m] = checkblocknum + 1;
                errorcheckindex[m] = j - Z * checkblocknum;
                m++;
            }
        }
        if (errorBits > 0) {
            Test.ErrorBits += errorBits;
            Test.ErrorFrame++;
            // date20190508
            //统计错误数小于3比特的帧，如果存在则增加迭代次数
            if (errorBits < 3) {
                Test.LT3ErrBitFrame++;
            }
            if (collect_err == 1) //是否输出错误信息开关
            {

                eout.open("errorindex.txt", std::ios::app);
                nout.open("errorfloat.txt", std::ios::app);
                dout.open("errordecode.txt", std::ios::app);

                // errorindex.txt
                eout << "ErrorFrame: " << i << endl; // 32帧的第几帧发生错误

                // information部分
                eout << "ErrorBit Num: " << k << endl; //该帧错误的比特数

                eout << "Errorbit Block: "; //那个块发生错误
                for (j = 0; j < k; j++)

                {
                    eout << errorbitblock[j] << "\t";
                }
                eout << endl;

                eout << "Errobit Index: "; //错误发生在块的哪一个位置
                for (j = 0; j < k; j++) {
                    eout << errorbitindex[j] << "\t";
                }
                eout << endl;

                // check部分
                eout << "Errorcheck Num: " << m << endl; //该帧错误的比特数
                eout << "Errorcheck Block: "; //那个块发生错误
                for (j = 0; j < m; j++) {
                    eout << errorcheckblock[j] << "\t";
                }
                eout << endl;

                eout << "Errorcheck Index: "; //错误发生在块的哪一个位置
                for (j = 0; j < m; j++) {
                    eout << errorcheckindex[j] << "\t";
                }
                eout << endl;
                eout.close();

                // ErrorNum.txt
                //输出过信道之后的float值和量化后的limitchar
                // date20190629 只需要bitsoverchannel
                for (j = 0; j < (NmoinsK - _PunctureBits); j++) //信息比特部分
                {
                    errorfloat[l] = (float)bpskinput[i * (NmoinsK - _PunctureBits) + j];
                    errorchar[l] = (int)charinput[i * (NmoinsK - _PunctureBits) + j];
                    l++;
                }

                for (j = 0; j < _NoCheck; j++) // check部分
                {
                    errorfloat[l] = (float)bpskinput[(NmoinsK - _PunctureBits) * 32 + i * _NoCheck + j];
                    errorchar[l] = (int)charinput[(NmoinsK - _PunctureBits) * 32 + i * _NoCheck + j];
                    l++;
                }

                // nout << "ErrorFrame Num: " << i << endl;//32帧的第几帧发生错误
                nout << "ErrorFloat=[ "; //通过awgn信道后的值，information+check
                for (j = 0; j < l; j++) {
                    nout << errorfloat[j] << "\t";
                }
                nout << "];" << endl;

                nout << "ErrorChar=[";
                for (j = 0; j < l; j++) {
                    nout << errorchar[j] << "\t";
                }
                nout << "];" << endl;
                nout << endl;
                nout.close();

                // errordecode.txt
                dout << "Decodedbits=[";
                for (j = 0; j < NOEUD; j++) {
                    dout << (int)decodedBits[i * NOEUD + j] << "\t";
                }
                dout << "];";
                dout << endl;
                dout << "inputbits=[";
                for (j = 0; j < (NmoinsK - _ShortenBits); j++) {
                    dout << (int)inputBits[i * (NmoinsK - _ShortenBits) + j] << "\t";
                }
                dout << "];";
                dout << endl;

                dout << "outputbits=[";
                for (j = 0; j < (NmoinsK - _PunctureBits); j++) {
                    dout << (int)outputBits[i * (NmoinsK - _PunctureBits) + j] << "\t";
                }
                for (j = 0; j < _NoCheck; j++) {
                    dout << (int)outputBits[32 * (NmoinsK - _PunctureBits) + i * _NoCheck + j] << "\t";
                }
                dout << "];";
                dout << endl;
                dout << endl;
                dout.close();
                for (j = 0; j <= k; j++) {
                    errorbitindex[j] = 0;
                    errorbitblock[j] = 0;
                }
                k = 0;
                for (j = 0; j <= m; j++) {
                    errorcheckindex[j] = 0;
                    errorcheckblock[j] = 0;
                }
                m = 0;
                for (j = 0; j <= l; j++) {
                    errorfloat[j] = 0;
                    errorchar[j] = 0;
                }
                l = 0;
            }
        }
    }
    return Test;
}

// count the weight of each VN
void CLDPC::VN_weight_count()
{
    for (int i = 0; i < _NoOnes; i++) {
        VN_weight_[PosNoeudsVariable[i]]++;
    }
}
