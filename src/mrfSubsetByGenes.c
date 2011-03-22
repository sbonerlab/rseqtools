#include <stdio.h>

#include <bios/log.h>
#include <bios/format.h>
#include <bios/intervalFind.h>

#include <mrf/mrf.h>



/** 
 *   \file mrfSubsetByGenes.c Module to subset MRF file by a set of genes.
 *         Usage: mrfSubsetByGenes <genes.interval> <samplePrefix>\n
 */



typedef struct {
  FILE *fp;
  Interval *gene;
  char *fileName;
} Entry;



static int sortEntries (Entry *a, Entry *b) 
{
  return a->gene - b->gene;
}



static int sortIntervals (Interval **a, Interval **b) 
{
  return *a - *b;
}



static void processAlignmentBlocks (Array intervals, Array alignmentBlocks)
{
  MrfBlock *currBlock;
  int i,j;
  Array genes;
  

  for (i = 0; i < arrayMax (alignmentBlocks); i++) {
    currBlock = arrp (alignmentBlocks,i,MrfBlock);
    genes = intervalFind_getOverlappingIntervals (currBlock->targetName,currBlock->targetStart - 1,currBlock->targetEnd);
    for (j = 0; j < arrayMax (genes); j++) {
      array (intervals,arrayMax (intervals),Interval*) = arru (genes,j,Interval*);
    }
  }
}



static void processEntry (Array entries, MrfEntry *currMrfEntry, char *prefix)
{
  int i;
  static Array intervals = NULL;
  static Stringa buffer = NULL;
  int index;
  Entry testEntry,*currEntry;

  if (intervals == NULL) {
    intervals = arrayCreate (1000,Interval*);
  }
  else {
    arrayClear (intervals);
  }
  processAlignmentBlocks (intervals,currMrfEntry->read1.blocks);
  if (currMrfEntry->isPairedEnd == 1) {
    processAlignmentBlocks (intervals,currMrfEntry->read2.blocks);
  }
  arraySort (intervals,(ARRAYORDERF)sortIntervals);
  arrayUniq (intervals,NULL,(ARRAYORDERF)sortIntervals);
  for (i = 0; i < arrayMax (intervals); i++) {
    testEntry.gene = arru (intervals,i,Interval*);
    testEntry.fp = NULL;
    if (arrayFindInsert (entries,&testEntry,&index,(ARRAYORDERF)sortEntries) == 1) {
      // new gene inserted
      currEntry = arrp (entries,index,Entry);
      stringCreateClear (buffer,100);
      stringPrintf (buffer,"%s_%s_%s:%d-%d.mrf",prefix,testEntry.gene->name,
                    testEntry.gene->chromosome,testEntry.gene->start,testEntry.gene->end);
      if (!(currEntry->fp = fopen (string (buffer),"w"))) {
        die ("Unable to open file: %s",string (buffer));
      }
      currEntry->fileName = hlr_strdup (string (buffer));
      fprintf (currEntry->fp,"%s\n",mrf_writeHeader ());
    }
    else {
      // gene is already present      
      ;
    }
    currEntry = arrp (entries,index,Entry);
    fprintf (currEntry->fp,"%s\n",mrf_writeEntry (currMrfEntry));
  }
}



int main (int argc, char *argv[])
{
  MrfEntry *currMRF;
  Entry *currEntry;
  Array entries;
  int i;

  if (argc != 3) {
    usage ("%s <genes.interval> <samplePrefix>",argv[0]);
  }
  intervalFind_addIntervalsToSearchSpace (argv[1],0);
  entries = arrayCreate (1000,Entry);
  mrf_init ("-");
  while (currMRF = mrf_nextEntry ()) {
    processEntry (entries,currMRF,argv[2]);
  }
  mrf_deInit ();
  for (i = 0; i < arrayMax (entries); i++) {
    currEntry = arrp (entries,i,Entry);
    fclose (currEntry->fp);
    warn ("Done with: %s",currEntry->fileName);
  }
  return 0;
}
