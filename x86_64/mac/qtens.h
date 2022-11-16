#if WORDSIZE == 32
#define NTEN 13
#define MINNTEN -8192
#define MAXNTEN 8192

#if NQ == 9
extern QELT qtens[NTEN+1][NQ]; 
extern QELT qmtens[NTEN+1][NQ];
#endif


#if NQ < 9
#define NTT 8
extern QELT qtens[NTEN+1][NTT];
extern QELT qmtens[NTEN+1][NTT];
#else

#define NTT NQ
extern QELT qtens[NTEN+1][NTT];
extern QELT qmtens[NTEN+1][NTT];
#endif
#endif
