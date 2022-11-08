# Unix makefile for libmq.a, qcalc, libstats.a qmtst

CPU = $(shell uname -m)
OSNAME=$(shell uname -s)
ASM=qasm.o shift.o qsquare.o

ifeq ($(CPU),x86_64)
ifeq ($(OSNAME),Darwin)
#ASM=$(CPU)/qasm-$(CPU)_$(OSNAME).o
echo "Do not use this directory for OS X-x86"
echo "Change directory to:"
echo "cd ./x86_64/mac"
exit
endif
endif

#ifeq ($(CPU),arm64)
F128=f128.o
#endif

OPTIM=-O3
#LDFLAGS=-fprofile-instr-generate
#PROFILE=-fcoverage-mapping -fprofile-instr-generate

ARITH = qsqrt.o qflti.o qfltbi.o bsr.o $(ASM)
#
# Set compiler
#
CC= gcc
#CC=~/lcc/book/lcc

ifeq ($(CC),gcc)
CFLAGS= -Wall -I. -g $(OPTIM)  $(PROFILE)
endif

ifeq ($(CC),~/lcc/book/lcc)
CFLAGS=-O -g -c -S -z
LCCLIBS=/usr/local/lcc/lib/lcclib.a
=======
CFLAGS= -Wall -I. -g $(OPTIM) $(PROFILE)
endif

ifeq ($(CC),~/lcc/book/lcc)
CFLAGS=-O3 -g -c
endif

HEADERS= mconf.h qcalc.h qhead.h qtens.h
LD=gcc
STATSLIBSRC= statslib/bdtrl.c statslib/betainv.c statslib/btdtrl.c statslib/chdtrl.c statslib/fdtrl.c statslib/igami.c\
statslib/igamil.c statslib/incbetl.c statslib/incbil.c statslib/kolmogor.c statslib/nbdtrl.c statslib/ndtril.c\
statslib/pdtrl.c statslib/psi.c statslib/simplestats.c statslib/stdtrl.c statslib/tbeta.c statslib/tbetainv.c\
statslib/test.c statslib/weibull.c
CFLAGS += -D$(CPU) -D$(OSNAME) 

ifeq ($(OSNAME),Darwin)
ifeq ($(CC),gcc)
CFLAGS += -fno-stack-protector
endif
endif

#for profiling
#CFLAGS +=-gline-tables-only

LIBOBJS= qacosh.o qairy.o qasin.o qasinh.o qatanh.o qatn.o qbeta.o \
qcbrt.o qcgamma.o qcmplx.o qchyp1f1.o qcerf.o qconst.o qcos.o qcosh.o \
qcpolylog.o qdawsn.o qei.o qellie.o qellik.o qellpe.o qellpj.o qellpk.o \
qerf.o qerfc.o qeuclid.o qexp.o qexp10.o qexp2.o qexpn.o qfac.o qhypot.o \
qfresf.o qgamma.o qhy2f1.o qhyperg.o qigam.o qigami.o qin.o qinv_fact.o \
qincb.o qine.o qjn.o qjypn.o qjyqn.o qagm.o qcatalan.o qnthroot.o \
qk0.o qkn.o qkne.o qlog.o qlog1.o qlog10.o qndtr.o qndtri.o qplanck.o \
qpolylog.o qpow.o qpsi.o qrand.o qremain.o qshici.o qsici.o qprob.o \
qsimq.o qsin.o qsindg.o qsinh.o qspenc.o qstudt.o qtan.o qtanh.o \
qincbi.o qyn.o qzetac.o qfloor.o mtherr.o $(F128) $(ARITH)


all: libmq.a qcalc qmtst qparanoi tsqrt statslib/libstats.a qtime

check: libmq.a qmtst qparanoi
	-qparanoi > temp.tmp
	-diff temp.tmp qparanoi.exp
	qmtst > temp.tmp
	-diff temp.tmp qmtst.exp

libmq.a: $(LIBOBJS) $(HEADERS) $(ASM) 
	rm -f libmq.a
	ar -r libmq.a $(LIBOBJS) $(ASM) 
	-ranlib libmq.a


qcalc: qcalc.o $(HEADERS) incbet.o incbi.o igami.o igam.o \
powi.o polevl.o ndtri.o \
const.o libmq.a statslib/libstats.a
	$(LD) $(LDFLAGS) -o qcalc qcalc.o incbet.o incbi.o ndtri.o \
igami.o igam.o powi.o polevl.o  \
const.o statslib/libstats.a  $(ASM)  \
libmq.a -lreadline -lcurses -lm $(LCCLIBS)
#libmq.a -lreadline -lm
#libmq.a -lm

qcalc.o: qcalc.c
	$(CC) $(CFLAGS) -DUSE_READLINE=1 -c qcalc.c
#	$(CC) $(CFLAGS) -c qcalc.c

mtherr.o: mtherr.c $(HEADERS)

const.o: const.c  $(HEADERS)

incbet.o: incbet.c  $(HEADERS)

incbi.o: incbi.c  $(HEADERS)

gamma.o: gamma.c  $(HEADERS)

igami.o: igami.c  $(HEADERS)

igam.o: igam.c  $(HEADERS)

exp.o: exp.c  $(HEADERS)

sin.o: sin.c $(HEADERS)

pow.o: pow.c $(HEADERS)

powi.o: powi.c $(HEADERS)

sqrt.o: sqrt.c  $(HEADERS)
f128.o: f128.c  $(HEADERS)



qccalc: $(HEADERS) qccalc.o cmplx.o const.o $(HEADERS) libmq.a
	$(LD) $(LDFLAGS) -o qccalc qccalc.o cmplx.o const.o \
libmq.a -lreadline -lm $(LCCLIBS)
#libmq.a -lreadline -lm
#libmq.a -lm

qccalc.o: qccalc.c  $(HEADERS)
	$(CC) $(CFLAGS) -DUSE_READLINE=1 -c qccalc.c
#	$(CC) $(CFLAGS) -c qccalc.c

cmplx.o: cmplx.c $(HEADERS)  $(HEADERS)



qmtst: qmtst.o ndtri.o polevl.o const.o \
drand.o $(HEADERS) libmq.a
	$(LD) $(LDFLAGS) -o qmtst qmtst.o ndtri.o  \
polevl.o const.o drand.o libmq.a -lm $(LCCLIBS)

qmtst.o: qmtst.c
ndtri.o: ndtri.c $(HEADERS)

polevl.o: polevl.c $(HEADERS)

drand.o: drand.c $(HEADERS)

log.o: log.c

qacosh.o: qacosh.c $(HEADERS)
qairy.o: qairy.c $(HEADERS)
qasin.o: qasin.c $(HEADERS)
qasinh.o: qasinh.c $(HEADERS)
qatanh.o: qatanh.c $(HEADERS)
qatn.o: qatn.c $(HEADERS)
qbeta.o: qbeta.c $(HEADERS)
qcbrt.o: qcbrt.c $(HEADERS)
qcgamma.o: qcgamma.c $(HEADERS)
qcmplx.o: qcmplx.c $(HEADERS)
qconst.o: qconst.c $(HEADERS)
qcos.o: qcos.c $(HEADERS)
qcosh.o: qcosh.c $(HEADERS)
qdawsn.o: qdawsn.c $(HEADERS)
qellie.o: qellie.c $(HEADERS)
qellik.o: qellik.c $(HEADERS)
qellpe.o: qellpe.c $(HEADERS)
qellpj.o: qellpj.c $(HEADERS)
qellpk.o: qellpk.c $(HEADERS)
qerf.o: qerf.c $(HEADERS)
qerfc.o: qerfc.c $(HEADERS)
qeuclid.o: qeuclid.c $(HEADERS)
qexp.o: qexp.c $(HEADERS) qinv_fact.c
qexp10.o: qexp10.c $(HEADERS)
qexp2.o: qexp2.c $(HEADERS)
qexpn.o: qexpn.c $(HEADERS)
qfac.o: qfac.c $(HEADERS)
qfresf.o: qfresf.c $(HEADERS)
qgamma.o: qgamma.c $(HEADERS)
qhy2f1.o: qhy2f1.c $(HEADERS)
qhyp.o: qhyp.c $(HEADERS)
qigam.o: qigam.c $(HEADERS)
qigami.o:qigami.c $(HEADERS)
qin.o: qin.c $(HEADERS)
qincb.o: qincb.c $(HEADERS)
qincbi.o:qincbi.c $(HEADERS)
qincg.o: qincg.c $(HEADERS)
qincgs.o: qincgs.c $(HEADERS)
qine.o: qine.c $(HEADERS)
qjn.o: qjn.c $(HEADERS)
qjypn.o: qjypn.c $(HEADERS)
qjyqn.o: qjyqn.c $(HEADERS)
qk0.o: qk0.c $(HEADERS)
qkn.o: qkn.c $(HEADERS)
qkne.o: qkne.c $(HEADERS)
qlog.o: qlog.c $(HEADERS)
qlog1.o: qlog1.c $(HEADERS)
qlog10.o: qlog10.c $(HEADERS)
qndtr.o: qndtr.c $(HEADERS)
qndtri.o: qndtri.c $(HEADERS)
qpow.o: qpow.c $(HEADERS)
qprob.o: qprob.c $(HEADERS)
qpsi.o: qpsi.c $(HEADERS)
qrand.o: qrand.c $(HEADERS)
qremain.o: qremain.c $(HEADERS)
qshici.o: qshici.c $(HEADERS)
qsici.o: qsici.c $(HEADERS)
qsimq.o: qsimq.c $(HEADERS)
qsin.o: qsin.c $(HEADERS) qinv_fact.c
qsindg.o: qsindg.c $(HEADERS)
qsinh.o: qsinh.c $(HEADERS)
qspenc.o: qspenc.c $(HEADERS)
qstudt.o: qstudt.c $(HEADERS)
qtan.o: qtan.c $(HEADERS)
qtanh.o: qtanh.c $(HEADERS)
qyn.o: qyn.c $(HEADERS)
qzetac.o: qzetac.c $(HEADERS)
qfloor.o: qfloor.c $(HEADERS)
mtherr.o: mtherr.c $(HEADERS)
qsqrt.o: qsqrt.c $(HEADERS)
qagm.o: qagm.c $(HEADERS)
qnthroot.o: qnthroot.c $(HEADERS)
qinv_fact.o: qinv_fact.c $(HEADERS)
qhypot.o: qhypot.c $(HEADERS)

#qparanoi: qparanoi.o qflt.o qflta.o qsqrta.o libmq.a
#	$(CC) $(CFLAGS) -o qparanoi qparanoi.o qflt.o qflta.o \
#qsqrta.o libmq.a -lm

# This will test the arithmetic that is actually in the library.
qparanoi: qparanoi.o libmq.a $(HEADERS)
	$(LD) $(LDFLAGS) -o qparanoi qparanoi.o libmq.a -lm $(LCCLIBS)

#qflt.o: qflt.c $(HEADERS)
#	$(CC) $(CFLAGS) -DSTICKY=1 -c qflt.c

qflta.o: qflta.c $(HEADERS)
	$(CC) $(CFLAGS) -c qflta.c

qflti.o: qflti.c $(HEADERS)
	$(CC) $(CFLAGS) $(DEFS) -c qflti.c

qparanoi.o: qparanoi.c $(HEADERS)
	$(CC) $(CFLAGS) -Wno-implicit -c qparanoi.c
#	$(CC) -O -c qparanoi.c

qsqrta.o: qsqrta.c $(HEADERS)

qtime:	qtime.o libmq.a 
	$(LD) $(LDFLAGS) -o qtime qtime.o libmq.a -lm $(LCCLIBS) 

qtime.o:	qtime.c
	$(CC) $(CFLAGS) -c qtime.c

# i386, coff version (DJGPP)
qfltb386.o: qfltbi.386
	as -o qfltb386.o qfltbi.386

# i386, ELF version (linux)
qfltbelf.o: qfltbelf.386
	as -o qfltbelf.o qfltbelf.386

qfltbi.o: qfltbi.c
	$(CC) $(CFLAGS) -c qfltbi.c

qf68k.o: qf68k.a
	as -o qf68k.o qf68k.a
bsr.o:	$(CPU)/bsr-$(CPU)_$(OSNAME).s
	as -c -o bsr.o $(CPU)/bsr-$(CPU)_$(OSNAME).s
tsqrt:	tsqrt.o libmq.a
	$(LD) $(LDFLAGS) -o tsqrt tsqrt.o libmq.a -lm
tsqrt.o:	tsqrt.c 
	$(CC) $(CFLAGS) -c tsqrt.c
statslib/libstats.a: $(STATSLIBSRC)
	cd statslib;make;cd ..

x86_64/qasm-x86_64_Darwin.s:	x86_64/qasm-x86_64.s
	sed -f sed.cmds x86_64/qasm-x86_64.s >x86_64/qasm-x86_64_Darwin.s
	as -c -g -o qasm.o x86_64/qasm-x86_64_Darwin.s

x86_64/bsr-x86_64_Darwin.o:	x86_64/bsr-x86_64.s
	sed -e "s/bsr64/_bsr64/" x86_64/bsr-x86_64.s >x86_64/bsr-x86_64_Darwin.s
	as -c -g -o x86_64/bsr-x86_64_Darwin.o x86_64/bsr-x86_64_Darwin.s

qasm.o:	$(CPU)/qasm-$(CPU)_$(OSNAME).s
	as -c -g -o qasm.o $(CPU)/qasm-$(CPU).s

shift.o: $(CPU)/shift-$(CPU).s
	as -g -c -o shift.o $(CPU)/shift-$(CPU).s
qsquare.o: $(CPU)/qsquare-$(CPU).s
	as -c -g -o qsquare.o $(CPU)/qsquare-$(CPU).s

clean:
	rm -f *.o *.lil *.s libmq.a qcalc qparanoi qmtst qtime tsqrt statslib/*.o $(CPU)/*.o
