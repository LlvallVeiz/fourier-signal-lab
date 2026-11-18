# Makefile -- fourier-signal-lab
# C11, no external dependencies (only libc/libm). `make test` builds and
# runs the assert-based unit test suite plus the numpy cross-check;
# `make examples` regenerates the example outputs under examples/.

CC       ?= gcc
CFLAGS   ?= -std=c11 -O2 -Wall -Wextra -Wpedantic -Isrc
LDFLAGS  ?=
LDLIBS   := -lm

SRC      := $(wildcard src/*.c)
LIB_SRC  := $(filter-out src/main.c,$(SRC))
TEST_SRC := $(wildcard tests/*.c)

TARGET      := fourier-lab
TEST_TARGET := run_tests
BENCH_TARGET := /tmp/benchmark_fft

.PHONY: all clean test examples bench

all: $(TARGET)

$(TARGET): $(SRC) $(wildcard src/*.h)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS) $(LDLIBS)

$(TEST_TARGET): $(LIB_SRC) $(TEST_SRC) $(wildcard src/*.h tests/*.h)
	$(CC) $(CFLAGS) $(LIB_SRC) $(TEST_SRC) -o $(TEST_TARGET) $(LDFLAGS) $(LDLIBS)

test: $(TEST_TARGET) $(TARGET)
	./$(TEST_TARGET)
	python3 tests/validate_numpy.py
	bash tests/run_tests.sh

examples: $(TARGET)
	bash examples/run_examples.sh

bench: $(TARGET)
	$(CC) $(CFLAGS) examples/benchmark_fft.c src/fft.c src/signal.c \
		-o $(BENCH_TARGET) $(LDFLAGS) $(LDLIBS)
	$(BENCH_TARGET)

clean:
	rm -f $(TARGET) $(TEST_TARGET)
	rm -rf build examples/outputs
