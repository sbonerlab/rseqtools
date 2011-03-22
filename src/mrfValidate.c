#include <bios/log.h>
#include <bios/format.h>
#include <bios/fasta.h>
#include <bios/seq.h>

#include <mrf/mrf.h>


/** 
 * \file mrfValidate.c Module to validate MRF files.
 * \author Andrea Sboner (andrea.sboner@yale.edu)
 * \brief Validation of MRF files.
 * \details It validates the correctness of MRF files. 
 *
 * Usage: mrfValidate <errors|noerros> <gaps|nogaps>\n
 * \par
 * \arg \c errors: prints the wrong MRF entries
 * \arg \c noerrors: prints the correct MRF entries
 * \arg \c gaps: the alignment tool allows gaps in the query sequence
 * \arg \c nogaps: the alignment tool does not allow gaps in the query sequence
 *
 * \par
 * Error codes:
 * \arg   0 : success
 * \arg   1 : query end      <= query start
 * \arg   2 : target end     <= target start
 * \arg   4 : query size     != target size
 * \arg   8 : target name 1  != target name 2
 * \arg  16 : strand 1       != strand 2
 * \arg  32 : query start 2  <= query start 1  OR query end 2 <= query end 1
 * \arg  64 : target start 2 <= target start 1  OR target end 2 <= target end 2
 * \arg 128 : ( query start 2 - query end 1 ) != 1 
 * \note
 * error codes greater than 4 occur only if the read has multiple blocks. Error code 128 can occur if \b nogaps is specified.
 */

int processBlock ( MrfBlock* currBlock, MrfBlock* prevBlock, int gaps ) {
  int errorCode = 0;
  if( currBlock->queryEnd < currBlock->queryStart )
    errorCode += 1;
  if( currBlock->targetEnd < currBlock->targetStart )
    errorCode += 2;
  int querySize  = currBlock->queryEnd  - currBlock->queryStart  + 1;
  int targetSize = currBlock->targetEnd - currBlock->targetStart + 1;
  if( targetSize != querySize )
    errorCode += 4;

  if( prevBlock != NULL ) {
    if( !strEqual( prevBlock->targetName, currBlock->targetName) )
      errorCode += 8;
    if( prevBlock->strand != currBlock->strand )
      errorCode += 16;
    if( currBlock->queryStart <= prevBlock->queryStart || currBlock->queryEnd <= prevBlock->queryEnd )
      errorCode += 32;
    if( currBlock->targetStart <= prevBlock->targetStart || currBlock->targetEnd <= prevBlock->targetEnd )
      errorCode += 64;
    if( gaps == 0 ) {
      if( (currBlock->queryStart - prevBlock->queryEnd ) != 1 )
	errorCode += 128;    
    }
  }
  return errorCode;
}

int main (int argc, char *argv[]) 
{
  if( argc != 3 ) {
    usage( "%s <errors|noerrors> <gaps|nogaps>", argv[0]);
    return -1;
  }
  MrfEntry *currEntry;
  MrfBlock *currBlock, *prevBlock;
  int i, error, errorCode;
  int errorCodesSummary[255];
  for( i=0; i<255; i++) 
    errorCodesSummary[i]=0;
  int gaps = 0;
  if( strEqual( "gaps", argv[2] ) )
    gaps = 1;
  seq_init();
  mrf_init ("-");
  printf( "%s\n", mrf_writeHeader());
  while (currEntry = mrf_nextEntry ()) {
    error = 0;
    for( i=0; i<arrayMax(currEntry->read1.blocks); i++) {   
      currBlock=arrp( currEntry->read1.blocks, i, MrfBlock);
      if( i==0 ) {
	errorCode = processBlock ( currBlock, NULL, gaps );
      } else {
	errorCode = processBlock ( currBlock, arrp( currEntry->read1.blocks, i-1, MrfBlock ), gaps );
      }
      if(  errorCode != 0 ) { 
	if( i==0 )
	warn( "Error code (read1): %d\t%s:%c:%d:%d:%d:%d", errorCode, currBlock->targetName, currBlock->strand, currBlock->targetStart, 
	      currBlock->targetEnd, currBlock->queryStart, currBlock->queryEnd);	
	else {
	  prevBlock=arrp( currEntry->read1.blocks, i-1, MrfBlock );
	  warn( "Error code (read1): %d\tblocks: [%d-%d]\t%s:%c:%d:%d:%d:%d,%s:%c:%d:%d:%d:%d",  errorCode, i-1, i, prevBlock->targetName, prevBlock->strand, prevBlock->targetStart, 
		prevBlock->targetEnd, prevBlock->queryStart, prevBlock->queryEnd, currBlock->targetName, currBlock->strand, currBlock->targetStart, 
		currBlock->targetEnd, currBlock->queryStart, currBlock->queryEnd);	
	}
	error++;
      }
      errorCodesSummary[errorCode]++; 
    }

    if (currEntry->isPairedEnd & error == 0) {
      for( i=0; i<arrayMax(currEntry->read2.blocks); i++) {
	currBlock=arrp( currEntry->read2.blocks, i, MrfBlock);
	if( i==0 ) {
	  errorCode = processBlock ( currBlock, NULL , gaps ); 
	} else {
	  errorCode = processBlock (currBlock, arrp( currEntry->read2.blocks, i-1, MrfBlock) , gaps);
	}
	if( errorCode !=0 ) {
	  if( i==0 ) {
	    warn( "Error code (read2): %d\t%s:%c:%d:%d:%d:%d", errorCode, currBlock->targetName, currBlock->strand, currBlock->targetStart, 
		  currBlock->targetEnd, currBlock->queryStart, currBlock->queryEnd);
	  }
	  else {
	    prevBlock=arrp( currEntry->read2.blocks, i-1, MrfBlock );
	    warn( "Error code (read1): %d\tblocks [%d-%d]\t%s:%c:%d:%d:%d:%d,%s:%c:%d:%d:%d:%d", errorCode, i-1, i, prevBlock->targetName, prevBlock->strand, prevBlock->targetStart, 
		  prevBlock->targetEnd, prevBlock->queryStart, prevBlock->queryEnd, currBlock->targetName, currBlock->strand, currBlock->targetStart, 
		  currBlock->targetEnd, currBlock->queryStart, currBlock->queryEnd);	
	  }
	  error++;
	}
	errorCodesSummary[errorCode]++; 
      }
    }    
   
    if( error > 0 && strEqual(argv[1], "errors")) {
      printf( "%s\n", mrf_writeEntry( currEntry ));
      error=0;
    } 
    if( error == 0 && strEqual(argv[1], "noerrors")) 
      printf( "%s\n", mrf_writeEntry( currEntry ));
  }
  mrf_deInit ();
  for( i = 0; i<255; i++ ) {
    if( errorCodesSummary[i]>0 ) 
      warn("%s: Error code (%d):\t%d", argv[0], i, errorCodesSummary[i] );
  }
  warn("%s: done", argv[0]);
  return 0;
}
