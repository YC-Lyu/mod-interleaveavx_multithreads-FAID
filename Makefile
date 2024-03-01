objects=main.o CTool.o CSimulate.o CModulate.o CChannel.o CLDPC.o CDecoder_OMS.o CDecoder_FAID.o CDecoder_OMSBF.o CDecoder_OMS_DTBF.o CDecoder_FAID_2B1C.o
main: $(objects)
	icc -std=c++11 -o main  $(objects) -lpthread -mkl -fopenmp
CTool.o: CTool.h
	icc -std=c++11 -c CTool.cpp 
main.o: CChannel.h CTool.h CLDPC.h CSimulate.h CModulate.h
	icc -std=c++11  -c main.cpp
CSimulate.o: CLDPC.h CChannel.h CModulate.h CTool.h CSimulate.h
	icc -std=c++11  -c CSimulate.cpp
CModulate.o: CTool.h CModulate.h
	icc -std=c++11  -c CModulate.cpp
CLDPC.o: CLDPC.h CTool.h ./Constants/Constants_SSE.h
	icc -std=c++11  -c CLDPC.cpp
CDecoder_OMS.o: CLDPC.h CTool.h ./Constants/Constants_SSE.h
	icc -std=c++11  -c CDecoder_OMS.cpp
CDecoder_FAID.o: CLDPC.h CTool.h ./Constants/Constants_SSE.h
	icc -std=c++11  -c CDecoder_FAID.cpp
CDecoder_OMSBF.o: CLDPC.h CTool.h ./Constants/Constants_SSE.h
	icc -std=c++11  -c CDecoder_OMSBF.cpp
CDecoder_OMS_DTBF.o: CLDPC.h CTool.h ./Constants/Constants_SSE.h
	icc -std=c++11  -c CDecoder_OMS_DTBF.cpp
CDecoder_FAID_2B1C.o: CLDPC.h CTool.h ./Constants/Constants_SSE.h
	icc -std=c++11  -c CDecoder_FAID_2B1C.cpp
CChannel.o: CTool.h CChannel.h
	icc -std=c++11 -c CChannel.cpp
clean:
	rm main $(objects)


