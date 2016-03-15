#Makefile

# Supported platforms
#       Unix / Linux               	LNX
#       Mac                        	MAC
# Compilation options
#       link to LAPACK              WITH_LAPACK
#       32-bit binary        		FORCE_32BIT

# Set this variable to either LNX or MAC
SYS = LNX
# Leave blank after "=" to disable; put "= 1" to enable
# Disable WITH_LAPACK option can slow computation speed significantly and is not recommended
# Disable WITH_ARPACK option only disable -apprx option in the software
WITH_LAPACK = 1
FORCE_32BIT = 
DIST_NAME = Estep_mcmc

# --------------------------------------------------------------------
# Edit below this line with caution
# --------------------------------------------------------------------


BIN_DIR  = ./bin

SRC_DIR  = ./src

CPP = g++

CPPFLAGS = -ggdb -Wall -O3 -I./libStatGen/include/ -I./zlib/ -D__ZLIB_AVAILABLE__ -D_FILE_OFFSET_BITS=64 -D__STDC_LIMIT_MACROS #-pg

LIBS = -lgsl -lgslcblas -pthread -lz -lm ./libStatGen/libStatGen.a

OUTPUT = $(BIN_DIR)/Estep_mcmc

SOURCES = $(SRC_DIR)/main.cpp

HDR = 

# Detailed library paths, D for dynamic and S for static

LIBS_LNX_D_LAPACK = -llapack
LIBS_MAC_D_LAPACK = -framework Veclib
LIBS_LNX_S_LAPACK = /usr/lib/lapack/liblapack.a -lgfortran  /usr/lib/atlas-base/libatlas.a /usr/lib/libblas/libblas.a -Wl,--allow-multiple-definition 

# Options


  SOURCES += $(SRC_DIR)/param.cpp $(SRC_DIR)/sfba.cpp $(SRC_DIR)/io.cpp $(SRC_DIR)/lm.cpp $(SRC_DIR)/lmm.cpp  $(SRC_DIR)/mvlmm.cpp $(SRC_DIR)/bvsrm.cpp $(SRC_DIR)/prdt.cpp $(SRC_DIR)/mathfunc.cpp $(SRC_DIR)/gzstream.cpp $(SRC_DIR)/ReadVCF.cpp $(SRC_DIR)/compress.cpp
  HDR += $(SRC_DIR)/param.h $(SRC_DIR)/sfba.h $(SRC_DIR)/io.h $(SRC_DIR)/lm.h $(SRC_DIR)/lmm.h  $(SRC_DIR)/mvlmm.h $(SRC_DIR)/bvsrm.h $(SRC_DIR)/prdt.h $(SRC_DIR)/mathfunc.h $(SRC_DIR)/gzstream.h $(SRC_DIR)/ReadVCF.h $(SRC_DIR)/compress.h


ifdef WITH_LAPACK
  OBJS += $(SRC_DIR)/lapack.o
  CPPFLAGS += -DWITH_LAPACK

ifeq ($(SYS), MAC)
  LIBS += $(LIBS_MAC_D_LAPACK)
else
  LIBS += $(LIBS_LNX_S_LAPACK)
endif

  SOURCES += $(SRC_DIR)/lapack.cpp
  HDR += $(SRC_DIR)/lapack.h
endif

ifdef FORCE_32BIT
  CPPFLAGS += -m32
else
  CPPFLAGS += -m64
endif


# all
OBJS = $(SOURCES:.cpp=.o)

all: $(OUTPUT)

$(OUTPUT): $(OBJS)
	$(CPP) $(CPPFLAGS) $(OBJS) $(LIBS) -o $(OUTPUT)

$(OBJS) : $(HDR)

.cpp.o: 
	$(CPP) $(CPPFLAGS) $(HEADERS) -c $*.cpp -o $*.o
.SUFFIXES : .cpp .c .o $(SUFFIXES)


clean:
	rm -rf ${SRC_DIR}/*.o ${SRC_DIR}/*~ *~ ${SRC_DIR}/*_float.* $(OUTPUT)

DIST_COMMON = COPYING.txt README.txt Makefile
DIST_SUBDIRS = src doc example bin

tar:
	mkdir -p ./$(DIST_NAME)
	cp $(DIST_COMMON) ./$(DIST_NAME)/
	cp -r $(DIST_SUBDIRS) ./$(DIST_NAME)/
	tar cvzf $(DIST_NAME).tar.gz ./$(DIST_NAME)/
	rm -r ./$(DIST_NAME)
