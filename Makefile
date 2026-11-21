NVCC ?= nvcc
ARCH ?= sm_86
NVCCFLAGS ?= -O3 -std=c++17 -arch=$(ARCH) --extended-lambda --expt-relaxed-constexpr -lineinfo -Xcompiler=-Wall,-Wextra
# deposit.cu is compiled once per dimension (-DSIC_DIM), so that make -j builds the three in parallel
PARTS = deposit1 deposit2 deposit3 fields mesh
HEADERS = src/sic.hpp Makefile
DOUBLE = $(PARTS:%=build/double/%.o)
SINGLE = $(PARTS:%=build/single/%.o)

.PHONY: all clean
all: build/sic build/sic-single
build/double/deposit%.o: src/deposit.cu $(HEADERS)
	@mkdir -p $(@D)
	$(NVCC) $(NVCCFLAGS) -DSIC_DIM=$* -c $< -o $@
build/single/deposit%.o: src/deposit.cu $(HEADERS)
	@mkdir -p $(@D)
	$(NVCC) $(NVCCFLAGS) -DSIC_SINGLE -DSIC_DIM=$* -c $< -o $@
build/double/%.o: src/%.cu $(HEADERS)
	@mkdir -p $(@D)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@
build/single/%.o: src/%.cu $(HEADERS)
	@mkdir -p $(@D)
	$(NVCC) $(NVCCFLAGS) -DSIC_SINGLE -c $< -o $@
build/sic: src/main.cu $(DOUBLE) $(HEADERS)
	$(NVCC) $(NVCCFLAGS) $< $(DOUBLE) -o $@
build/sic-single: src/main.cu $(SINGLE) $(HEADERS)
	$(NVCC) $(NVCCFLAGS) -DSIC_SINGLE $< $(SINGLE) -o $@
build/test: tests/test.cu $(DOUBLE) $(HEADERS) | build/sic  # (the tests also run the program)
	$(NVCC) $(NVCCFLAGS) -Isrc $< $(DOUBLE) -o $@
build/test-single: tests/test.cu $(SINGLE) $(HEADERS) | build/sic-single
	$(NVCC) $(NVCCFLAGS) -DSIC_SINGLE -Isrc $< $(SINGLE) -o $@
clean:
	rm -rf build
