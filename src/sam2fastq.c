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
 *   \file sam2fastq.c Program to convert SAM paired end files to FASTQ
 */



#define R_FIRST		0
#define R_SECOND	1




void writeSamToFastq( SamEntry *r1, int end ) {
  
}

int main (int argc, char **argv)
{

  char* pos;
  SamEntry *currSamE = NULL;

  samParser_initFromFile("-");
  while( currSamE = samParser_nextEntry() ) {
    char *pos = strchr( currSamE->qname, '/'); // checking /1 or /2 in the read name
    if( pos != NULL ) {
      if( !( currSamE->flags & S_FIRST) && !( currSamE->flags & S_SECOND ) ) {
	if( strEqual( pos, "/1") ) {
	  writeSamToFastq( currEntry, 1);
	} else 	  {
	  writeSamToFastq( currEntry, 2);
	}
      }
    } else die("Error");
  }
  samParser_deInit();
}
