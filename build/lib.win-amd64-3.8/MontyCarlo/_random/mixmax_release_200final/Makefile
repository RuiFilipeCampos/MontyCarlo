ifndef CC
CC=gcc
endif

ifndef CXX
CXX=g++
endif

ifdef _N
NFLAG=-D_N=${_N}
endif

ifndef PREFIX
PREFIX=/usr/local
endif

ifeq ($(CC), icc)
CFLAGS=-std=c99 -O3 -funroll-loops  -Wall
else
CFLAGS=-std=c99 -O3 -funroll-loops -Wall -I ${PREFIX}/include/
endif

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
DYNLIBFLAGS=-dynamiclib -Wl,-dead_strip_dylibs -Wl,-headerpad_max_install_names,-undefined,dynamic_lookup,-compatibility_version,1.0,-current_version,1.0,-install_name,/usr/local/lib/libmixmax.dylib
dynlib=libmixmax.dylib
else
CFLAGS+= -fPIC -lc
DYNLIBFLAGS=-shared -fPIC -Wl,-dead_strip_dylibs
dynlib=libmixmax.so
endif


objects = mixmax.o
mainobj = driver_main.o
skipobj = driver_skip_charpoly.o
verobj = driver_verification.o

.c.o:
	${CC} ${CFLAGS}  ${NFLAG} -c $<

mixmax: $(objects) $(mainobj) mixmax.c driver_main.c mixmax.h
	${CC} ${CFLAGS} -o mixmax $(objects) $(mainobj)

# multithreaded demo
multi: mixmax.c mixmax.h driver_threads.c 
	${CC} ${CFLAGS}  ${NFLAG}  -o multi mixmax.c driver_threads.c -lpthread

all: mixmax verify testU01 lib check

verify: $(objects) $(verobj) mixmax.c mixmax.h driver_verification.c 
	${CC} ${CFLAGS}  ${NFLAG} -o verify $(objects) $(verobj)
	
check:  clean mixmax verify
	@ echo ""
	@ echo [Test 1] - Testing correct floating point output
	@ if [ -e N240_1200.diff ]; then rm N240_1200.diff; fi
	@ echo "1200" | ./mixmax 2> /dev/null > N240_1200_local.out
	@ echo "Outputting the first 1200 numbers into N240_1200_local.out"
	@ diff -w N240_1200.out N240_1200_local.out > N240_1200.diff
	@ echo "The diff file N240_1200.diff should be empty."
	@ /bin/ls -1s N240_1200.diff
	@ if [ -s N240_1200.diff ]; then  echo "<<Files differ!>> Test failed."; else echo "Test 1 passed"; fi
	@ echo ""
	@ echo "[Test 2] - Writing and reading the state from file"
	${CC} ${CFLAGS}  ${NFLAG}  -o testIO driver_iotest.c mixmax.c
	@ ./testIO
	@ echo ""
	@ echo "[Test 3] - Checking the explicit MIXMAX matrix"
	@ if [ -e mathematica.txt.gz ]; then gunzip mathematica.txt.gz; fi
	@ echo 100 | ./verify 2> /dev/null > mathematica.local.txt
	@ diff -w mathematica.N240.txt mathematica.local.txt > mathematica.txt.diff
	@ echo "The diff file mathematica.txt.diff should be empty."
	@ /bin/ls -1s mathematica.txt.diff
	@ if [ -s mathematica.txt.diff ]; then  echo "<<Files differ!>> Test failed."; else echo "Test 3 passed"; fi
	@ echo ""
	@ echo "[Test 4] - Optional"
	@ echo "You may now run the contents of mathematica.local.txt in Mathematica (yourself) to get an independent check of correctness."

gsl: mixmax.c mixmax.h driver_gsl.c
	${CC} ${CFLAGS} driver_gsl.c mixmax.c -lgsl -L${PREFIX}/lib -DHOOKUP_GSL=1  -o gsl
	@ echo "Run it with:\n GSL_RNG_SEED=123 GSL_RNG_TYPE=MIXMAX ./gsl"

testU01: mixmax.c mixmax.h driver_testU01.c 
	${CC} ${CFLAGS}  ${NFLAG} -o testU01 driver_testU01.c -ltestu01  mixmax.c -mtune=native

lib: mixmax.o mixmax.h
	@ echo "Dynamic library is ${dynlib}"
	ar rcs libmixmax.a mixmax.o
	gcc ${DYNLIBFLAGS} -o ${dynlib} mixmax.o

install: lib
	cp mixmax.h ${PREFIX}/include/
	cp libmixmax.a ${PREFIX}/lib/
	cp ${dynlib} ${PREFIX}/lib/

asm: mixmax.c mixmax.h
	${CC} ${CFLAGS}  ${NFLAG} -S mixmax.c

src: mixmax.c mixmax.h
	gcc-mp-5 -std=c99 -O3 -funroll-loops -Wall -I /usr/local/include/   mixmax.c  -mtune=native -S
	zip -9 src.zip mixmax.c mixmax.s Makefile mixmax.c mixmax.h driver_testU01.c \
	driver_main.c mixmax_skip_N240.c mixmax_skip_N17.c \
	driver_verification.c driver_threads.c driver_gsl.c driver_iotest.c driver_time.c README.pdf N240_1200.out mathematica.N240.txt mixmax.hpp example.cpp
	
clean:
	- rm $(objects) $(mainobj) $(skipobj) $(verobj) verify N256_1200_local.out N256_1200.diff mixmax gsl \
	testU01 libmixmax.* a.out

plusplus:
	${CXX} -std=c++0x -O3 -funroll-loops example.cpp -o mixmax.plusplus
