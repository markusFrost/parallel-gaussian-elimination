COMP    = icc
CFLAGS  = -O2 -pthread
TARGET  = gauss_pthread
SRC     = src/Gaussian.c src/Gaussian_parallel.c src/Gaussian_Block_Parallel.c src/Gaussian_Reuse_Parallel.c 
DEP     = src/Gaussian.h src/Gaussian_parallel.h src/Gaussian_Block_Parallel.h src/Gaussian_Reuse_Parallel.h

all: $(TARGET)

$(TARGET): $(DEP)
	$(COMP) $(CFLAGS) -o $@ $(SRC)

clean:
	rm -f $(TARGET)

