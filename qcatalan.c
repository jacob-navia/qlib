#include <stdio.h>
#include <string.h>
/*
1,1,2,5,14,42,132,429,1430,4862,16796,58786,208012,742900,2674440,
9694845,35357670,129644790,477638700,1767263190,6564120420,
24466267020,91482563640,343059613650,1289904147324,4861946401452,
18367353072152,69533550916004,263747951750360,1002242216651368,3814986502092304

With 128 bit integers we arrive at catalan(62).
*/
 
void print_uint128(unsigned __int128 n,char *buffer) 
{
  char str[45] = {0}; // log10(1 << 128) + '\0'
  char *s = str + sizeof(str) - 2; // start at the end
  char c;
  if (n == 0) {
     buffer[0]='0';
     buffer[1]=0;
	return ;
  } 

  while (n != 0) {

    c  = '0' + (n % 10); // save last digit
	*s-- = c;
    n /= 10;                     // drop it
  }
  strcpy(buffer,++s);
}
// Returns value of Binomial Coefficient C(n, k)
unsigned __int128  binomialCoeff(unsigned __int128 n,unsigned  __int128 k)
{
    unsigned __int128  res = 1;
 
    // Since C(n, k) = C(n, n-k)
    if (k > n - k)
        k = n - k;
 
    // Calculate value of [n*(n-1)*---*(n-k+1)] / [k*(k-1)*---*1]
    for (int i = 0; i < k; ++i)
    {
        res *= (n - i);
        res /= (i + 1);
    }
 
    return res;
}
 
// A Binomial coefficient based function to find nth catalan
// number in O(n) time
unsigned __int128  catalan(unsigned __int128 n)
{
    // Calculate value of 2nCn
   unsigned  __int128  c = binomialCoeff(2*n, n);
 
    // return 2nCn/(n+1)
    return c/(n+1);
}

#include "mconf.h"
void qcatalan(Qfloatp n,Qfloatp result)
{
	Qfloat nm1[1],ip1[1],r[1];
	unsigned __int128 n128,n128x2,res;
	long long i;

	n128 = qtoi128(n);
	if ((__int128)n128 < 0) {
		qclear(result);
		return;
	}
	if (n128 < 2) {
		qmov(qone,result);
		return;
	}

	if (n128 < 63) {
		n128 = catalan(n128);
		i128toq(n128,result);
		return;
	}
	res = 1;
	n128x2 = n128 * 2;
	for (i=0;i<10;i++) {
		res *= (n128x2-i);
		res /= (i+1);
	}
	i128toq(res,r);
	//qmov(qone,r);
	for (; i<n128;i++) {
		i128toq(n128x2-i,nm1);
		i128toq(i+1,ip1);
		qmul(nm1,r,r);
		qdiv(ip1,r,r);
		//qfloor(r,r);
	}
	i128toq(n128+1,ip1);
	qdiv(ip1,r,r);
	qfloor(r,result);
}
 
#ifdef STANDALONE
// Driver program to test above functions
int main(int argc,char *argv[])
{
	int start,top;
	char buf[60],buf1[60];
	if (argc == 2) {
		start = atoi(argv[1]);
		top = start+1;
	}
	else {
		top = 63;
		start = 1;
	}
    for (int i = start; i < top; i++) {
		print_uint128(catalan(i),buf);
		print_uint128(binomialCoeff(2*i,i),buf1);
        printf("[%3d], %s     %s\n",i,buf,buf1);
	}
    return 0;
}
#endif
