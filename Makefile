CC=			gcc
CXX=		g++
CFLAGS=		-std=c99 -g -Wall -O3
CXXFLAGS=	$(CFLAGS)
CPPFLAGS=
INCLUDES=
OBJS=		kagraph.o dict.o cleanup.o mrna.o read.o format.o
PROG=		gffio
LIBS=		-lz -lm

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address -ldl -lpthread
endif

.SUFFIXES:.c .cpp .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

gffio:$(OBJS) main.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr *.o a.out $(PROG) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

cleanup.o: gfpriv.h gffio.h kagraph.h
dict.o: gfpriv.h gffio.h khashl.h
format.o: gfpriv.h gffio.h
kagraph.o: kagraph.h
main.o: gffio.h ketopt.h
mrna.o: gfpriv.h gffio.h ksort.h
read.o: gfpriv.h gffio.h kseq.h
