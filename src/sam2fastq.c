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

void writeSamToFastq(SamEntry *r1, int end) {
  
}

int main (int argc, char **argv) {
  SamParser* parser = samParser_initFromFile("-");
  for (SamEntry* entry = NULL; entry = samparser_next_entry(parser); ) {
    char *pos = strchr(entry->qname, '/'); // checking /1 or /2 in the read name
    if (pos != NULL) {
      if (!(entry->flags & S_FIRST) && !(entry->flags & S_SECOND)) {
        if (strEqual(pos, "/1") ) {
          writeSamToFastq(currEntry, 1);
        } else {
          writeSamToFastq(currEntry, 2);
        }
      }
    } else {
      die("Error");
    }
  }
  samparser_free(parser);
}
