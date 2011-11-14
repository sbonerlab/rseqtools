#include <bios/format.h>
#include <bios/log.h>
#include <bios/intervalFind.h>

/** 
 *   \file intervalSize.c Module to compute the size of the intervals in nucleotides.
 *         Usage: intervalSize <file.interval> \n
 *         Takes Interval format from fileName.
 */


int main ( int argc, char* argv[] ) {
  int i,j;
  Array intervals =  intervalFind_parseFile (argv[1],0);
  for( i=0; i<arrayMax(intervals); i++) {
    Interval* currInterval = arrp(intervals, i, Interval);
    int size=0;
    for( j=0; j<arrayMax(currInterval->subIntervals); j++) {
      SubInterval* si = arrp( currInterval->subIntervals, j, SubInterval );
      size += (si->end - si->start );
    }	
    printf( "%s\t%d\n", currInterval->name, size);
  }
  arrayDestroy( intervals );
  return 0;
} 
