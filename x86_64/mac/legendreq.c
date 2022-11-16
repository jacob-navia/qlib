#include <qfloat.h>
qfloat EXPORT legendreq( int n, qfloat x)
{
  int i=0;
  qfloat cx0,cx1,cxi;

  if ( n < 0 )
  {
    return 0.0;
  }
  if (n == 0.0)
    return 1.0Q;
  if (n == 1)
    return x;

  cx0 = 1.0;

  cx1 = x;

  for ( i = 2; i <= n; i++ )
  {
    cxi = ( ( qfloat ) ( 2 * i - 1 ) * x * cx1 
            + ( qfloat ) (   - i + 1 )     * cx0 ) 
            / ( qfloat ) (     i     );
    cx0=cx1;
    cx1=cxi;

  }

  return cxi;
}

