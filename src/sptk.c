#include <stdio.h>


void next_line(FILE *filep){
	 /*
      *  Skips to next input line.
      */

  int dummy;
  while( (dummy=getc(filep)) != '\n');
}


int ncord(int l, int i, int j, int ix, int iy)
{
   int nc,i1,i2;

   i1=i+ix;
   i2=j+iy;

/* Periodic boundary conditions */
   if (i1>=l) i1-=l;
   if (i1<0) i1+=l;
   if (i2>=l) i2-=l;
   if (i2<0) i2+=l;

/* No-flux boundary conditions */
//    if (i1>=l) i1=l-1;
//    if (i1<0) i1=0;
//    if (i2>=l) i2=l-1;
//    if (i2<0) i2=0;

   nc=i1+i2*l;

   return nc;
}

