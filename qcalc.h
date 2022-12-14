/*		calc.h
 * include file for calc.c
 */
 
/* Set LARGEMEM 1 for 32-bit wide memory addresses
 * or to 0 for 16 bit addresses (small memory model on a PC)
 */
#ifndef LARGEMEM
#define LARGEMEM 1
#endif

/* data structure of symbol table */
struct symbol
	{
	char *spel;
	unsigned short attrib;
#if LARGEMEM
	long sym;
#else
	short sym;
#endif
	};

struct funent
	{
	char *spel;
	unsigned short attrib;
	int (*fun )();
	};

struct varent
        {
	char *spel;
	unsigned short attrib;
	Qfloatp value;
        };

struct strent
	{
	char *spel;
	unsigned short attrib;
	char *string;
	};


/*	general symbol attributes:	*/
#define OPR 0x8000
#define	VAR 0x4000
#define CONST 0x2000
#define FUNC 0x1000
#define ILLEG 0x800
#define BUSY 0x400
#define TEMP 0x200
#define STRING 0x100
#define COMMAN 0x80
#define MENU 0x40
#define IND 0x1

/* attributes of operators (ordered by precedence): */
#define BOL 1
#define EOL 2
/* end of expression (comma): */
#define EOE 3
#define EQU 4
#define PLUS 5
#define MINUS 6
#define MULT 7
#define DIV 8
#define UMINUS 9
#define LPAREN 10
#define RPAREN 11
#define COMP 12
#define MOD 13
#define LAND 14
#define LOR 15
#define LXOR 16


extern struct funent funtbl[];
/*extern struct symbol symtbl[];*/
extern struct varent indtbl[];

