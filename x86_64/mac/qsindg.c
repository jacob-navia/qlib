/*	qsindg.c
* sin, cos, tan in degrees
*/

#include "qhead.h"

/* 1.745329251994329576923690768488612713442871888541725456E-2 */
static Qfloat pi180[1] = {
	0,EXPONE-6,0x8efa351294e9c8aeULL,0x0ec5f66e9485c4d9ULL,0x00b7aef501b5e6b8ULL,
	0xe502a9b4c94c8512ULL,0xb6f6116781911487ULL,0x10c50c969d5140c9ULL,0x60d4a6b49598f1edULL};


void qsindg(Qfloatp x,Qfloatp y)
{
	Qfloat w[1];

	qmul( x, pi180, w );
	qfsin( w, y );
}


void qcosdg(Qfloatp x,Qfloatp y)
{
	Qfloat w[1];

	qmul( x, pi180, w );
	qfcos( w, y );
}


void qtandg(Qfloatp x,Qfloatp y)
{
	Qfloat w[1];

	qmul( x, pi180, w );
	qftan( w, y );
}

