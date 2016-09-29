include ../make.inc

# ROOTS
_KBLAS_ROOT_=..

# include and lib paths
INCLUDES=
INCLUDES+= -I${_KBLAS_ROOT_}/include

LIB_PATH=
LIB_PATH+= -L${_KBLAS_ROOT_}/lib

# libraries to link against
LIB= -lm -l$(LIB_KBLAS_NAME)
LIB+= $(BLAS_LIB)

KBLAS_LIB=../lib/lib${LIB_KBLAS_NAME}.a

C_SRC=

C_SRC+=	              \
	test_strmm.c  \
	test_dtrmm.c  \
	test_ctrmm.c  \
	test_ztrmm.c  \
	test_strsm.c  \
	test_dtrsm.c  \
	test_ctrsm.c  \
	test_ztrsm.c
	
#C_SRC+=	              \
#	test_dgemm.c



INC_SRC=testing_utils.h

OBJ_DIR=./obj
C_OBJ = $(patsubst %.c,$(OBJ_DIR)/%.o,$(C_SRC))

BIN_DIR=./bin
EXE=$(patsubst %.c,$(BIN_DIR)/%,$(C_SRC))

$(OBJ_DIR)/test_%trmm.o: test_%trmm.c test_Xtrmm.ch testing_utils.h $(INC_SRC) $(KBLAS_LIB)
	$(CXX) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/test_%trsm.o: test_%trsm.c test_Xtrsm.ch testing_utils.h $(INC_SRC) $(KBLAS_LIB)
	$(CXX) $(INCLUDES) -c $< -o $@

#$(OBJ_DIR)/test_%gemm.o: test_%gemm.c test_Xgemm.ch testing_utils.h $(INC_SRC) $(KBLAS_LIB)
#	$(CXX) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/%.o: %.c $(INC_SRC) $(KBLAS_LIB)
	$(CC) $(COPTS) $(INCLUDES) -c $< -o $@

$(BIN_DIR)/test_%: $(OBJ_DIR)/test_%.o 
	$(CC) $< -o $@ $(LIB_PATH) $(LIB)

all: $(EXE)

$(EXE): $(C_OBJ)

clean:
	rm -f $(OBJ_DIR)/*.o $(EXE)

