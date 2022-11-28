CC=			gcc
CXX=		g++
CFLAGS=		-std=c99 -g -Wall -O3
CXXFLAGS=	$(CFLAGS)
CPPFLAGS=
INCLUDES=
OBJS=		kagraph.o dict.o regular.o query.o read.o format.o
PROG=		minigff
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

minigff:$(OBJS) main.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr *.o a.out $(PROG) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

dict.o: mgf-priv.h minigff.h khashl.h
format.o: mgf-priv.h minigff.h
kagraph.o: kagraph.h
main.o: minigff.h ketopt.h
query.o: mgf-priv.h minigff.h ksort.h
read.o: mgf-priv.h minigff.h kseq.h
regular.o: mgf-priv.h minigff.h kagraph.h
