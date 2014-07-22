GPUDIR   = /usr/local/cuda-4.2
SDKDIR   = $(GPUDIR)/NVIDIA_GPU_Computing_SDK/C
MEXDIR   = /usr/local/matlab2012a

NVCC     = $(GPUDIR)/bin/nvcc
MWCC	 = $(MEXDIR)/bin/mex -largeArrayDims
INCLUDES = -I$(MEXDIR)/extern/include -I$(SDKDIR)/common/inc
LIBS     = -L$(GPUDIR)/lib64 -L$(SDKDIR)/lib -lcudart
NVFLAGS  = -Xcompiler -fPIC --ptxas-options -v -arch sm_20 #-maxrregcount 16

EXE_NAME = xcorr1d
HDRFILES = 
CPPFILES =
CU_FILES = main.cu
OBJFILES = main.o

all: $(CU_FILES) $(CPPFILES) $(OBJFILES)
	$(MWCC) $(INCLUDES) $(LIBS) $(OBJFILES) -cxx -output $(EXE_NAME)
	rm *.o

%.o: %.cu $(HDRFILES)
	$(NVCC) $(INCLUDES) $(NVFLAGS) $(LIBS) $< -c -o $@
