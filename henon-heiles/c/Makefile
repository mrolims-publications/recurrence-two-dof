CC      ?= gcc
CSTD    ?= c99

SRC_DIR   := src
PROG_DIR  := programs
INC_DIR   := include
OBJ_DIR   := build/obj
BIN_DIR   := bin

WARN      := -Wall -Wextra -Wpedantic
OPT       := -O3

CPPFLAGS  := -I$(INC_DIR)
CFLAGS    := $(WARN) $(OPT) -std=$(CSTD)
LDFLAGS   :=
LDLIBS    := -lm

USE_OMP   ?= 0

SRC       := $(wildcard $(SRC_DIR)/*.c)
OBJS      := $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(SRC))
DEPS      := $(OBJS:.o=.d)

PROGS_SRC := $(wildcard $(PROG_DIR)/*.c)
PROGS     := $(patsubst $(PROG_DIR)/%.c,$(BIN_DIR)/%,$(PROGS_SRC))

# Detect compiler family
CC_VERSION := $(shell $(CC) --version 2>/dev/null | head -n 1)

IS_ICX         := $(findstring icx,$(CC))
IS_APPLE_CLANG := $(findstring Apple clang,$(CC_VERSION))
IS_CLANG       := $(findstring clang,$(CC_VERSION))
IS_GCC         := $(findstring GCC,$(CC_VERSION))

# LTO
ifeq ($(IS_ICX),icx)
    CFLAGS  += -ipo
    LDFLAGS += -ipo
else
    CFLAGS  += -flto
    LDFLAGS += -flto
endif

# OpenMP
ifeq ($(USE_OMP),1)
    ifeq ($(IS_ICX),icx)
        CFLAGS  += -qopenmp
        LDFLAGS += -qopenmp

    else ifneq ($(IS_APPLE_CLANG),)
        # Apple clang: plain -fopenmp does not work
        # Assumes libomp installed, e.g. via Homebrew
        OMP_PREFIX ?= /opt/homebrew/opt/libomp

        CPPFLAGS  += -I$(OMP_PREFIX)/include
        CFLAGS    += -Xpreprocessor -fopenmp
        LDFLAGS   += -L$(OMP_PREFIX)/lib
        LDLIBS    += -lomp

    else ifneq ($(IS_CLANG),)
        # Upstream LLVM clang
        CFLAGS  += -fopenmp
        LDFLAGS += -fopenmp

    else
        # GCC
        CFLAGS  += -fopenmp
        LDFLAGS += -fopenmp
    endif
endif

.PHONY: all programs clean
all: $(OBJS)

programs: $(PROGS)

%: $(BIN_DIR)/%
	@:

$(OBJ_DIR):
	mkdir -p $@

$(BIN_DIR):
	mkdir -p $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -MMD -MP -c $< -o $@

$(BIN_DIR)/%: $(PROG_DIR)/%.c $(OBJS) | $(BIN_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(OBJS) $(LDFLAGS) $(LDLIBS) -o $@.x

clean:
	rm -rf build $(BIN_DIR)

-include $(DEPS)
