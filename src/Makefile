include ../make.inc

KBLAS_ROOT= ..
KBLAS_SRC =

KBLAS_SRC = strmm.cpp dtrmm.cpp ctrmm.cpp ztrmm.cpp Xtr_common.cpp
KBLAS_SRC+= strsm.cpp dtrsm.cpp ctrsm.cpp ztrsm.cpp

OBJ_DIR=./obj
ALL_OBJ = $(patsubst %.cpp,$(OBJ_DIR)/%.coo,$(KBLAS_SRC))

INC=-I${KBLAS_ROOT}/include
DEPS=../make.inc Xtr_common.ch ${KBLAS_ROOT}/include/kblas.h

$(OBJ_DIR)/%tr_common.coo: %tr_common.cpp $(DEPS)
	$(CXX) $(COPTS) $(INC) -c $< -o $@

$(OBJ_DIR)/%trmm.coo: %trmm.cpp $(DEPS)
	$(CXX) $(COPTS) $(INC) -c $< -o $@

$(OBJ_DIR)/%trsm.coo: %trsm.cpp $(DEPS)
	$(CXX) $(COPTS) $(INC) -c $< -o $@


KBLAS_DIR=$(KBLAS_ROOT)/lib
LIBKBLAS=$(KBLAS_DIR)/lib$(LIB_KBLAS_NAME).a

all: $(LIBKBLAS)

$(LIBKBLAS): $(ALL_OBJ) 
	ar -cru $(LIBKBLAS) $(ALL_OBJ)
	ranlib $@
	
#cp $(LIBKBLAS) ${KBLAS_ROOT}/lib
#rm -f ${KBLAS_ROOT}/testing/*.o
	
clean:
	rm -f $(OBJ_DIR)/*.coo

