#include <stdlib.h>
#include <time.h>

#include <bios/log.h>
#include <bios/format.h>
#include <mrf/mrf.h>



/** 
 *   \file mrfSampler.c Module to sample reads from MRF.
 *         Usage: mrfSampler <proportionOfReadsToSample> \n
 *         Takes MRF from STDIN. \n
 */



int main (int argc, char *argv[])
{
  MrfEntry *currEntry;
  double proportion;

  if (argc != 2) {
    usage ("%s <proportionOfReadsToSample>",argv[0]);
  }
  proportion = atof (argv[1]);
  srand (time (0));
  mrf_init ("-"); 
  puts (mrf_writeHeader ());
  while (currEntry = mrf_nextEntry ()) {
    if ((1.0 * rand () / RAND_MAX) > proportion) {
      continue;
    }  
    puts (mrf_writeEntry (currEntry));
  }
  mrf_deInit (); 
  return 0;
}
