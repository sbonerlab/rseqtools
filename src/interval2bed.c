#include <bios/format.h>
#include <bios/log.h>
#include <bios/intervalFind.h>



/** 
 *   \file interval2bed.c Module to convert Interval format to BED format.
 *         Usage: interval2bed <trackName> \n
 *         Takes Interval format from stdin.
 */


int main (int argc, char *argv[])
{
  int i,j;
  Array intervals;
  Interval *currInterval;
  SubInterval *currSubInterval;
  Stringa sizes =NULL;
  Stringa starts=NULL;

  if (argc < 2) {
    usage ("%s <trackName> [simple]",argv[0]);
  }
  if( (argc==3) && !strEqual( argv[2],"simple") ) {
    usage("%s <trackName> [simple]",argv[0]);
  }
  intervalFind_addIntervalsToSearchSpace ("-",0);
  intervals = intervalFind_getAllIntervals ();
  puts ("browser hide all");
  printf ("track name=\"%s\" visibility=2\n",argv[1]);
  for (i = 0; i < arrayMax (intervals); i++) {
    currInterval = arrp (intervals,i,Interval);
    if( argc == 3 ) {
      for (j = 0; j < arrayMax (currInterval->subIntervals); j++) {
	currSubInterval = arrp (currInterval->subIntervals,j,SubInterval);
	printf ("%s\t%d\t%d\t%s\t900\t%c\n",
		currInterval->chromosome,currSubInterval->start,currSubInterval->end,currInterval->name,currInterval->strand);
      }
    } else {
      stringCreateClear( starts, 10);
      stringCreateClear( sizes,   10);
      for( j = 0; j < arrayMax (currInterval->subIntervals); j++) {
	currSubInterval = arrp (currInterval->subIntervals,j,SubInterval);
	stringAppendf( sizes, "%d", currSubInterval->end - currSubInterval->start );
        stringAppendf( starts, "%d", currSubInterval->start - currInterval->start );
	if( j<arrayMax( currInterval->subIntervals) ) {
	  stringAppendf( sizes, "," );
	  stringAppendf(starts, "," );
	}
      }
      printf ("%s\t%d\t%d\t%s\t900\t%c\t%d\t%d\t.\t%d\t%s\t%s\n",
	      currInterval->chromosome,currInterval->start,currInterval->end, currInterval->name,currInterval->strand, currInterval->start, currInterval->end, currInterval->subIntervalCount, string(sizes), string(starts) );
    }
  }
  return 0;
}
