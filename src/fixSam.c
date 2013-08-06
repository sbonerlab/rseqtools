#include <stdlib.h>
#include <string.h>

#include <bios/format.h>
#include <bios/log.h>
#include <bios/linestream.h>
#include <bios/common.h>
#include <bios/seq.h>

#include <mrf/mrf.h>
#include <mrf/sam.h>

/** 
 *   \file fixSam.c Program to fix SAM files
 *   @pre The sam file should have already excluded unmapped and non primary alignments and also be sorted by read name.
 */



#define R_FIRST		0
#define R_SECOND	1

void fixSingleEnd( SamEntry *currSamE ) {
  if( ( currSamE->flags & S_QUERY_UNMAPPED )  && !(strEqual( currSamE->rname, "*") ) )
    currSamE->flags -= S_QUERY_UNMAPPED;
  if( ( currSamE->flags & S_MATE_UNMAPPED )  && !(strEqual( currSamE->mrnm, "*") ) )
    currSamE->flags -= S_MATE_UNMAPPED;
					       
  if( strEqual( currSamE->mrnm, "*" ) && !(currSamE->flags & S_MATE_UNMAPPED) ) {
    currSamE->flags += S_MATE_UNMAPPED;
    if( currSamE->flags & S_PAIR_MAPPED ) currSamE->flags -= S_PAIR_MAPPED ; 
    if( currSamE->flags & S_MATE_STRAND ) currSamE->flags -= S_MATE_STRAND;
  }
  char *pos = strchr( currSamE->qname, '/'); // removing /1 or /2 in the read name
  if( pos != NULL ) {
      if( ( currSamE->flags & S_FIRST) && ( currSamE->flags & S_SECOND ) ) {    
	if( strEqual( pos, "/1") ) 
	  currSamE->flags -= S_SECOND ;
	else 
	  currSamE->flags -= S_FIRST ;
      }
      if( !( currSamE->flags & S_FIRST) && !( currSamE->flags & S_SECOND ) ) {
 	if( strEqual( pos, "/1") ) 
	  currSamE->flags -= S_SECOND ;
	else 
	  currSamE->flags -= S_FIRST ;
      }
      *pos='\0';
  }
  if( ( currSamE->flags & S_QUERY_UNMAPPED ) && ( currSamE->flags & S_NOT_PRIMARY ) ) currSamE->flags -= S_NOT_PRIMARY; 
}

void fixPairedEnds( SamEntry *r1, SamEntry *r2 ) {
  if( ( r1->flags & S_MATE_UNMAPPED) && !( r2->flags & S_QUERY_UNMAPPED) ) {
    r1->flags -= S_MATE_UNMAPPED;
    r1->mrnm = hlr_strdup( r2->rname );
    r2->mrnm = hlr_strdup( r1->rname );
    r1->mpos = r2->pos;
    r2->mpos = r1->pos;
  } 
 
  if( !( r1->flags & S_QUERY_UNMAPPED ) && ( r2->flags & S_MATE_UNMAPPED) ) {
    r2->flags -= S_MATE_UNMAPPED;
    r1->mrnm = hlr_strdup( r2->rname );
    r2->mrnm = hlr_strdup( r1->rname );
    r1->mpos = r2->pos;
    r2->mpos = r1->pos;
  }

  if( ( ((r1->flags & S_MATE_STRAND )>0) != ((r2->flags & S_QUERY_STRAND )>0 ) ) ) 
    (  r2->flags & S_QUERY_STRAND ) ? ( r1->flags += S_MATE_STRAND) : ( r1->flags -= S_MATE_STRAND );
  if( ( ((r2->flags & S_MATE_STRAND )>0) != ((r1->flags & S_QUERY_STRAND )>0 ) ) ) 
    (r1->flags & S_QUERY_STRAND ) ? (r2->flags += S_MATE_STRAND) : ( r2->flags -= S_MATE_STRAND );

  if( strEqual( r1->rname, r2->rname ) ) {
    r1->isize = r1->pos - r2->pos;
    r2->isize = -1*r1->isize;
  }
  if( ! ( r1->flags & S_READ_PAIRED ) )  r1->flags += S_READ_PAIRED;
  if( ! ( r2->flags & S_READ_PAIRED ) )  r2->flags += S_READ_PAIRED;
  if(   ( r1->flags & S_FIRST )  && !( r2->flags & S_SECOND ) ) r2->flags += S_SECOND;
  if(   ( r1->flags & S_SECOND ) && !( r2->flags & S_FIRST ) )  r2->flags += S_FIRST;
  if(   ( r2->flags & S_FIRST )  && !( r1->flags & S_SECOND ) ) r1->flags += S_SECOND;
  if(   ( r2->flags & S_SECOND ) && !( r1->flags & S_FIRST ) )  r1->flags += S_FIRST; 
}

int main (int argc, char **argv)
{
  
  char* pos;
  SamEntry *currSamE = NULL;
  SamEntry *mateSamE = NULL;
  SamEntry *tmpSamE = NULL;
  SamEntry *r1 = NULL;
  SamEntry *r2 = NULL;
  SamEntry *r3 = NULL;
  short unsigned int flag=1;

  samParser_initFromFile("-");
  while( tmpSamE = samParser_nextEntry() ) {
    samParser_freeEntry( currSamE );currSamE = NULL;
    samParser_copyEntry( &currSamE, tmpSamE );
    fixSingleEnd( currSamE );
    tmpSamE = samParser_nextEntry();
    if( !tmpSamE ) break;
    samParser_freeEntry( mateSamE );mateSamE = NULL; 
    samParser_copyEntry( &mateSamE, tmpSamE );
    fixSingleEnd( mateSamE ); 
    if(  strEqual( currSamE->qname, mateSamE->qname ) ) {
      fixPairedEnds( currSamE, mateSamE );
      printf("%s\n", samParser_writeEntry( currSamE ) );
      printf("%s\n", samParser_writeEntry( mateSamE ) );
      if( r1 ) samParser_freeEntry( r1 ); r1 = NULL;
      if( r2 ) samParser_freeEntry( r2 ); r2 = NULL;
      if( r3 ) samParser_freeEntry( r3 ); r3 = NULL;
      continue;
    } else {
      if( r1==NULL ) {
	samParser_copyEntry( &r1, mateSamE );
	if( r3 ) {
	  if( strEqual( r3->qname, currSamE->qname ) ) {
	    fixPairedEnds( r3, currSamE );
	    printf("%s\n", samParser_writeEntry( r3 ) );
	    printf("%s\n", samParser_writeEntry( currSamE ) );
	    samParser_freeEntry( r3 );r3=NULL;
	    continue;
	  } else {
	    samParser_freeEntry( r3 ); r3=NULL;
	  }
	}
	continue;
      } else if( r2==NULL ) {
	samParser_copyEntry( &r2, currSamE );
	samParser_copyEntry( &r3, mateSamE );
      }    
    }
    if( r1->qname && r2->qname ) {
      if( strEqual( r1->qname, r2->qname ) ) {
	fixPairedEnds( r1, r2 );
	printf("%s\n", samParser_writeEntry( r1 ) );
	printf("%s\n", samParser_writeEntry( r2 ) );
	samParser_freeEntry( r1 );r1=NULL;
	samParser_freeEntry( r2 );r2=NULL;
      }
    }	   
  }
  if( r1 ) { samParser_freeEntry( r1 ); r1=NULL; }
  if( r2 ) { samParser_freeEntry( r2 ); r2=NULL; }
  samParser_deInit();
}
