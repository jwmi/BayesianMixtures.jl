/* Combined congruential and Tauseworthe generators from SuperDuper
 * package. Should work on machines with unsigned long of at least 32
 * bits. JC and JT must be initialized to values with 0 < JC < 2^32 and
 * 0 < JT < 2^32. JC must be odd.
 * References: Marsaglia, Ananthanarayanan & Paul, 1973,
 * Learmonth & Lewis, 1973, and Dudewicz, 1976)
 */

/* Compilation flags:
   -DSUBR sdrand is used as a subroutine sdrand(u) 
          instead of a function u = sdrand()
   -DSUNF generate Sun Fortran compatible entry names
          (with trailing _ symbol)
   -DDOUB sdrand or its argument are double precision
   -DBSD  solves problem with Sun Fortran, pre-Solaris versions
          (i.e. BSD), with single precision function option.
	  Do not use this option for Solaris.
   -DRETS returns effective seed in argument to sdrni;
          if used, it is essential to use a variable, not
          a constant, in calling sdrni.
   -DLOG  prints seed value in file "rnilog".

   Examples:
	cc -c -o sd.o -DSUNF -DRETS sd.c
   (single precision Sun Fortran function, returning seed value from sdrni)

	cc -c -o sdc.o sd.c
   (single precision C function)
*/

#ifdef LOG
#include<stdio.h>
#endif

#ifdef SUNF
#ifndef DOUB
#include<math.h>
#endif
#endif

static unsigned long JC, JT;
static double Norm=4.656612873E-10;

#ifdef SUNF
#ifdef SUBR
  void sdrand_(u)
#else
#ifdef DOUB
  double sdrand_()
#else
#ifdef BSD
  FLOATFUNCTIONTYPE sdrand_()
#else
  float sdrand_()
#endif
#endif
#endif
#else
#ifdef SUBR
  void sdrand(u)
#else
#ifdef DOUB
  double sdrand()
#else
  float sdrand()
#endif
#endif
#endif

#ifdef SUBR
#ifdef DOUB
  double *u;
#else
  float *u;
#endif
#endif

{
  JC = (JC * 69069) & 037777777777;	/* congruential part */
  JT ^= JT >> 15;	/* tausworthe part */
  JT ^= (JT << 17) & 037777777777;

#ifdef SUBR
#ifdef DOUB
  *u = ((JT ^ JC) >> 1) * Norm;
#else
  *u = (float)(((JT ^ JC) >> 1) * Norm);
#endif
#else
#ifdef DOUB
  return(((JT ^ JC) >> 1) * Norm);
#else
#ifdef SUNF
#ifdef BSD
  RETURNFLOAT((float)((JT ^ JC) >> 1) * Norm);
#else
  return((float)((JT ^ JC) >> 1) * Norm);
#endif
#else
  return((float)((JT ^ JC) >> 1) * Norm);
#endif
#endif
#endif

}

#ifdef SUNF
  void sdrni_(i)
#else
  void sdrni(i)
#endif

unsigned long *i;
{
#ifdef LOG
	FILE *stream;
#endif
	unsigned long k=*i;
	if(k==0) k=time(0);
	JT = k/65536; JC = k-65536*JT;
	JT = 65536*JT+1; JC = 32768*JC+1;
#ifdef LOG
	stream = fopen("rnilog","a+");
	fprintf(stream,"%12d\n",k);
	fclose(stream);
#endif
#ifdef RETS
	*i = k;
#endif
}
