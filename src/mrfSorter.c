#include <stdlib.h>
#include <time.h>

#include <bios/log.h>
#include <bios/format.h>
#include <bios/common.h>
#include <mrf/mrf.h>



/** 
 *   \file mrfSorter.c Module to sample reads from MRF.
 *         Usage: mrfSampler <proportionOfReadsToSample> \n
 *         Takes MRF from STDIN. \n
 */

void destroyBlocks( MrfRead *read) {
  int j;
  MrfBlock *currBlock;
  for( j=0;j<arrayMax(read->blocks); j++ ) {
    currBlock = arrp( read->blocks, j, MrfBlock );
    hlr_free(currBlock->targetName);
  }
}

int copyBlock(MrfRead *src, MrfRead *dest ) {
  int i;
  MrfBlock *currDestBlock, *currSrcBlock;
  if( dest->blocks!=NULL )  {
    destroyBlocks( dest );
    arrayDestroy ( dest->blocks );
  }
  dest->blocks = arrayCreate( arrayMax(src->blocks), MrfBlock);
    
  for( i=0; i<arrayMax(src->blocks); i++ ) {
    currSrcBlock = arrp( src->blocks, i, MrfBlock);
    currDestBlock = arrayp( dest->blocks, i, MrfBlock );
    currDestBlock->targetName  = hlr_strdup( currSrcBlock->targetName );
    currDestBlock->strand      = currSrcBlock->strand;
    currDestBlock->targetStart = currSrcBlock->targetStart;
    currDestBlock->targetEnd   = currSrcBlock->targetEnd;
    currDestBlock->queryStart  = currSrcBlock->queryStart;
    currDestBlock->queryEnd    = currSrcBlock->queryEnd;
  }
  return 0;
}

int copyRead( MrfRead *src, MrfRead *dest) {
  if( src->sequence ) dest->sequence = hlr_strdup(src->sequence);
  if( src->qualityScores )  dest->qualityScores = hlr_strdup(src->qualityScores);
  //if( src->queryId )  dest->queryId = hlr_strdup(src->queryId); // no need to swap ID
  copyBlock( src, dest );  
}  

/**
   \function swapReads: change the order of the reads
   [in] read1, read2 of type MrfReads
 */
int swapReads( MrfRead *read1, MrfRead *read2 ) {
   MrfRead* tmpRead;
   AllocVar( tmpRead );
   tmpRead->blocks=NULL;
   copyRead( read1, tmpRead );
   copyRead( read2, read1 );
   copyRead( tmpRead, read2 );
   free(tmpRead);
   return 0;
}

int main (int argc, char *argv[])
{
  MrfEntry *currEntry;
  MrfBlock *currReadCoord1,  *currReadCoord2;

  mrf_init ("-"); 
  puts (mrf_writeHeader ());
  while (currEntry = mrf_nextEntry ()) {
    if( currEntry->isPairedEnd ) {
      currReadCoord1 = arrp(currEntry->read1.blocks, 0, MrfBlock);
      currReadCoord2 = arrp(currEntry->read2.blocks, 0, MrfBlock);
      if( strcmp(currReadCoord1->targetName, currReadCoord2->targetName)>0 )
	swapReads( &currEntry->read1, &currEntry->read2 );
      else {
	if( strEqual( currReadCoord1->targetName, currReadCoord2->targetName) ) {
	  if( currReadCoord1->targetStart > currReadCoord2->targetStart )
	    swapReads( &currEntry->read1, &currEntry->read2 );
	}
      } 
    }
    puts (mrf_writeEntry ( currEntry )) ;
  }
  mrf_deInit (); 
  return 0;
}
