#include <bios/format.h>
#include <bios/log.h>
#include <bios/bgrParser.h>
#include <bios/linestream.h>
#include <bios/intervalFind.h>
#include <bios/common.h>


static void getValuesForInterval( Array bgrs, Interval* currInterval, int* length, double* value) 
{
  Array entries;
  int i,j;
  *value = 0.0;
  *length = 0;
  for( i=0; i < arrayMax( currInterval->subIntervals ); i++ ) {
    SubInterval* currSubInterval = arrp( currInterval->subIntervals, i, SubInterval ); 
    *length += currSubInterval->end - currSubInterval->start;
    entries = bgrParser_getValuesForRegion( bgrs, currInterval->chromosome, currSubInterval->start, currSubInterval->end);
    for( j=0; j<arrayMax( entries ); j++) 
      *value += arru( entries, j, double );
  }
  arrayDestroy( entries );
}

/**
 * \file bgrQuantifier <annotation.interval>.
 * \pre: it requires a BedGraph file from STDIN normalized by the number of mapped nucleotides
 */
int main( int argc, char* argv[] ) {
  Array bgrs;
  Array intervals;
  int i, length;
  double value;
  if( argc < 2 ) {
    usage("%s <annotation.interval>\n%s requires a BedGraph from STDIN", argv[0], argv[0]);
  }
  bgrs = arrayCreate( 1000, BedGraph );
  bgrParser_initFromFile ( "-" );
  bgrs = bgrParser_getAllEntries ();
  bgrParser_deInit();
  arraySort( bgrs, (ARRAYORDERF) bgrParser_sort );
  
  intervalFind_addIntervalsToSearchSpace ( argv[1], 0 );
  intervals = intervalFind_getAllIntervals ();
  
  for( i=0; i<arrayMax(intervals); i++ ) {
    Interval *currInterval = arrp( intervals, i, Interval );
    getValuesForInterval( bgrs, currInterval, &length, &value);
    
    printf("%s\t%s:%d-%d\t%f\n", currInterval->name, 
	   currInterval->chromosome, 
	   currInterval->start+1, 
	   currInterval->end, 
	   value /= length / 1000.0 );       
  }
  arrayDestroy( intervals );
  return 0;
}
