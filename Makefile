HOMEL=$(PWD)    
BIN=$(HOMEL)//bin
IDIR=-I.   -Inicksrc

ND = ./nicksrc
NLIB = $(ND)/libnick.a

override LDLIBS +=  -lm -lgsl -lopenblas -lfftw3 $(NLIB)  
override CFLAGS += -c -g -p -Wimplicit  $(IDIR)  

override CFLAGS += -DHAVE_CONFIG_H -Wimplicit-int -D_IOLIB=2

all:   dates_expfit dates simpjack2 dowtjack grabpars

$(NLIB): 
	$(MAKE) -C $(ND) 

dates_expfit:    $(NLIB) dates_expfit.o   regsubs.o  fitexp.o   gslfit.o
dates: $(NLIB) dates.o qpsubs.o mcio.o ldsubs.o admutils.o egsubs.o regsubs.o fftsubs.o
dowtjack: $(NLIB) dowtjack.o
simpjack2: $(NLIB) simpjack2.o 
grabpars: $(NLIB) grabpars.o 

clean:  
	rm -f  *.o  nicksrc/*.o nicklib

clobber:  clean
	rm -f  dates_expfit dates simpjack2 grabpars dowtjack nicksrc/libnick.a bin/*

install: all
	mkdir -p ./bin
	cp dates_expfit dates grabpars simpjack2 dowtjack perlsrc/*  ./bin/.
	chmod 777 bin/*
