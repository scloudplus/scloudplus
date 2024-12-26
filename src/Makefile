CC = gcc
ifeq "$(CC)" "gcc"
	COMPILER = gcc
else ifeq "$(CC)" "clang"
	COMPILER = clang
endif

CFLAGS = -O3 -march=native -g

LDFLAGS = -lm -lcrypto

SOURCES = aes_ni.c encode.c random.c fips202.c sample.c pke.c matrix.c util.c kem.c
        
HEADERS = aes.h encode.h random.h fips202.h sample.h pke.h matrix.h util.h kem.h ds_benchmark.h config.h

ifeq ($(OS),Windows_NT)
	UNAME_P := $(PROCESSOR_ARCHITECTURE)
	RM = cmd /C del
	RMDIR = cmd /C rmdir /S /Q
	TARGET := $(TARGET).exe

else
	UNAME_P := $(shell uname -p)
	RM = rm -f
	RMDIR = rm -rf
endif

ifeq ($(UNAME_P),AMD64)
	CFLAGS += -DAMD
else ifeq ($(UNAME_P),arm64)
	CFLAGS += -DARM
else ifeq ($(UNAME_P),x86_64)
	CFLAGS += -DAMD
else ifeq ($(UNAME_P),arm)
	CFLAGS += -DARM
endif



.PHONY: all speed clean

all: scloudplus

scloudplus:\
	scloudplus128_aes\
	scloudplus192_aes\
	scloudplus256_aes\



scloudplus128_aes: $(HEADERS) $(SOURCES) test.c
	$(CC) $(CFLAGS) -Dscloudplus_l=128  test.c -o $@ $(SOURCES) $(LDFLAGS)

scloudplus192_aes: $(HEADERS) $(SOURCES) test.c
	$(CC) $(CFLAGS) -Dscloudplus_l=192  test.c -o $@ $(SOURCES) $(LDFLAGS)

scloudplus256_aes: $(HEADERS) $(SOURCES) test.c
	$(CC) $(CFLAGS) -Dscloudplus_l=256  test.c -o $@ $(SOURCES) $(LDFLAGS)


clean:
	-$(RM) scloudplus128*
	-$(RM) scloudplus192*
	-$(RM) scloudplus256*

