LINKER=gcc
LIB=qlib.a
# 32-bit words, C language.
# qfltbi requires compiler support of 64-bit long long ints.
#ARITH = qsqrt.o qflti.o qasm.o 

# 32-bit integer array elements, i386 assembly language
#ARITH = qsqrt.o qflti.o qfltbi.o pack.o a128mul.o a128shl.o a128shr.o au128div.o bsr.o
ARITH= qsqrt.o qflti.o qasm.o

CC=gcc
CFLAGS= -I. -DSTICKY=1 -DINFINITIES=1 -DNANS=1 -DDENORMAL=1 -DMINUSZERO=1  -g
ASFLAGS=
INCS= qcalc.h qhead.h qtens.h

LIBOBJS= qacosh.o qasin.o qatn.o qatanh.o qasinh.o qremquo.o \
qconst.o qcos.o qcosh.o qtanh.o qcbrt.o qsinh.o qgamma.o qigam.o qpsi.o \
qexp.o qexp10.o qexp2.o qfac.o qerf.o qerfc.o qrand.o qjn.o qyn.o \
qlog.o qlog1.o qlog10.o  qpow.o qremain.o qcmplx.o \
qsin.o qtan.o  qfloor.o mtherr.o qexpm1.o qhyp.o qairy.o \
qstudt.o qndtri.o qndtr.o qincbi.o qpolylog.o qzetac.o \
qincb.o qellpk.o qhy2f1.o qexpn.o $(ARITH)

EXTENDED=igami.o igam.o qfresf.o qellpj.o qshici.o qeuclid.o \
qigami.o qprob.o qei.o qsindg.o btdtrl.o \
qcpolylog.o qplanck.o qcgamma.o qsici.o \
qellpe.o polevll.o incbetl.o \
qin.o qine.o qjypn.o qjyqn.o qk0.o ndtril.o \
qkn.o qchyp1f1.o qcerf.o qellie.o qellik.o betainv.o bdtrl.o incbil.o \

all: qfloat.a qcalc qmtst qparanoi tsqrt

qfloat.a:	$(LIBOBJS) $(INCS) makefile $(EXTENDED)
	rm qfloat.a;ar q qfloat.a $(LIBOBJS) $(EXTENDED)

qcalc: qcalc.o $(INCS)  qfloat.a
	gcc -o qcalc qcalc.o  qfloat.a 

qcalc.o: qcalc.c
	$(CC) $(CFLAGS) -c qcalc.c

mtherr.o: mtherr.c $(INCS)

const.o: const.c
	$(CC) $(CFLAGS) -c const.c

incbet.o: incbet.c

betainv.o:	betainv.c

incbi.o: incbi.c

gamma.o: gamma.c

igami.o: igami.c

igam.o: igam.c

exp.o: exp.c

sin.o: sin.c

pow.o: pow.c

powi.o: powi.c


qccalc: qccalc.o qdawsn.o cmplx.o const.o qkne.o $(INCS) qfloat.a
	gcc -o qccalc qccalc.o cmplx.o polevl.o const.o $(EXTENDED) ./qfloat.a 

qccalc.o: qccalc.c
	$(CC) $(CFLAGS) -DUSE_READLINE=1 -c qccalc.c

cmplx.o: cmplx.c $(INCS)

qmtst: qmtst.o polevl.o \
const.o drand.o $(INCS) $(EXTENDED) ./qfloat.a
	gcc -o qmtst qmtst.o drand.o ./qfloat.a

qmtst.o: qmtst.c
ndtri.o: ndtri.c $(INCS)

polevl.o: polevl.c $(INCS)

drand.o: drand.c $(INCS)


qacosh.o: qacosh.c
qairy.o: qairy.c
qasin.o: qasin.c
qasinh.o: qasinh.c
qatanh.o: qatanh.c
qatn.o: qatn.c
qbeta.o: qbeta.c
qcbrt.o: qcbrt.c
qcgamma.o: qcgamma.c
qcmplx.o: qcmplx.c
qconst.o: qconst.c
qcos.o: qcos.c
qcosh.o: qcosh.c
qdawsn.o: qdawsn.c
qellie.o: qellie.c
qellik.o: qellik.c
qellpe.o: qellpe.c
qellpj.o: qellpj.c
qellpk.o: qellpk.c
qerf.o: qerf.c
qerfc.o: qerfc.c
qeuclid.o: qeuclid.c
qexp.o: qexp.c
qexp10.o: qexp10.c
qexp2.o: qexp2.c
qexpn.o: qexpn.c
qfac.o: qfac.c
qfresf.o: qfresf.c
qgamma.o: qgamma.c
qhy2f1.o: qhy2f1.c
qhyp.o: qhyp.c
qigam.o: qigam.c
qigami.o:qigami.c
qin.o: qin.c
qincb.o: qincb.c
qincbi.o:qincbi.c
qine.o: qine.c
qjn.o: qjn.c
qjypn.o: qjypn.c
qjyqn.o: qjyqn.c
qk0.o: qk0.c
qkn.o: qkn.c
qkne.o: qkne.c
qlog1.o: qlog1.c
qlog10.o: qlog10.c
qndtr.o: qndtr.c
qndtri.o: qndtri.c
qpow.o: qpow.c
qprob.o: qprob.c
qpsi.o: qpsi.c
qrand.o: qrand.c
qremain.o: qremain.c
qshici.o: qshici.c
qsici.o: qsici.c
qsimq.o: qsimq.c
qsin.o: qsin.c
qsindg.o: qsindg.c
qsinh.o: qsinh.c
qspenc.o: qspenc.c
qstudt.o: qstudt.c
qtan.o: qtan.c
qtanh.o: qtanh.c
qyn.o: qyn.c
qzetac.o: qzetac.c
qfloor.o: qfloor.c
mtherr.o: mtherr.c
qsqrt.o: qsqrt.c


# This will test the arithmetic that is actually in the library.
qparanoi: qparanoi.o qfloat.a 
	gcc -o qparanoi -g qparanoi.o qfloat.a

qflt.o: qflt.c $(INCS)
	$(CC) $(CFLAGS) -DSTICKY=1 -c qflt.c

qflta.o: qflta.c $(INCS)
	$(CC) $(CFLAGS) -c qflta.c

qflti.o: qflti.c $(INCS)
	$(CC) $(CFLAGS) $(DEFS) -c qflti.c

qparanoi.o: qparanoi.c $(INCS)
	$(CC) $(CFLAGS) -c qparanoi.c

qsqrta.o: qsqrta.c $(INCS)


#qfltbi.o: qfltbi.c
#	$(CC) $(CFLAGS) qfltbi.c
qasm.o:	qasm.s
	as -g $(ASFLAGS) qasm.s -o qasm.o
a128shl.o:	a128shl.asm
	lcc a128shl.asm

a128shr.o:	a128shr.asm
	lcc a128shr.asm
au128div.o:	au128div.asm
	lcc -c au128div.asm
bsr.o:	bsr.asm
	lcc bsr.asm
tsqrt:	tsqrt.o qfloat.a
	gcc -o tsqrt tsqrt.o qfloat.a
tsqrt.o:	tsqrt.c
	$(CC) $(CFLAGS) -c tsqrt.c
clean:
	rm -f *.o qfloat.lib qccalc.exe qparanoi qmtst qcalc tsqrt qfloat.a
