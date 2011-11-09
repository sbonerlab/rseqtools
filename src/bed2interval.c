#include <bios/format.h>
#include <bios/log.h>
#include <bios/intervalFind.h>
#include <bios/bedParser.h>
#include <bios/common.h>



/** 
 *   \file bed2interval.c Module to convert BED format into Interval format.
 *         Takes BED from STDIN \n
 *         Note: The Interval and BED formats are zero-based and half-open. \n
 */



int main (int argc, char *argv[]) 
{
  Array beds;
  Bed *currBed;
  int i,j;

  bedParser_initFromFile ("-");
  j=1;
  while( currBed = bedParser_nextEntry() ) {
    Interval* currInterval;
    AllocVar( currInterval );
    currInterval->chromosome = hlr_strdup( currBed->chromosome );
    currInterval->start = currBed->start;
    currInterval->end = currBed->end;
    if( currBed->extended ) {
      currInterval->name = hlr_strdup( currBed->name );
      currInterval->strand = currBed->strand;
      currInterval->subIntervalCount = currBed->blockCount;
      currInterval->subIntervals = arrayCreate( currBed->blockCount, SubInterval );
      for( i=0; i<currBed->blockCount; i++) { 
	SubInterval *currSubInterval = arrayp( currInterval->subIntervals, arrayMax( currInterval->subIntervals), SubInterval);
	SubBlock* currSubBlock = arrp( currBed->subBlocks, i, SubBlock );
	currSubInterval->start = currSubBlock->start + currBed->start;
	currSubInterval->end = currSubInterval->start + currSubBlock->size;
      }
    } else {
      Stringa tmp = stringCreate( 10 );
      stringPrintf( tmp, "BED_%d", j);
      currInterval->name = hlr_strdup( string( tmp )  );
      stringDestroy( tmp );
      currInterval->strand = '.';
      currInterval->subIntervalCount = 1;
      currInterval->subIntervals = arrayCreate( 1, SubInterval );
      SubInterval *currSubInterval = arrayp( currInterval->subIntervals, arrayMax( currInterval->subIntervals), SubInterval);
      currSubInterval->start = currBed->start;
      currSubInterval->end = currBed->end;
    }
    printf( "%s\n", intervalFind_writeInterval( currInterval ) );
    freeMem( currInterval );
    j++;
  }
  return 0;
}
