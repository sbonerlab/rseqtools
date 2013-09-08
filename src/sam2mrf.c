#define _GNU_SOURCE

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
 *   \file sam2mrf.c Module to convert SAM to MRF.
 */

#define R_FIRST     0
#define R_SECOND    1

static void printMrfAlignBlocks(SamEntry *entry, int strand) {
  char mrf_strand = '.';
  Texta tokens = textFieldtokP(entry->cigar, "MN");
  if (strand == R_FIRST) {
    if (entry->flags & S_QUERY_STRAND) {
      mrf_strand = '-';
    } else {
      mrf_strand = '+';
    }
  } else {
    if (entry->flags & S_MATE_STRAND) {
      mrf_strand = '-';
    } else {
      mrf_strand = '+';
    }
  }

  // Process first item in cigar
  int len = atoi(textItem(tokens, 0));
  int pos = entry->pos;
  int q = 1;
  printf("%s:%c:%d:%d:%d:%d",
         entry->rname, mrf_strand, entry->pos, pos + len - 1, 1, len);
  pos += len - 1;
  q += len;

  // Process rest of cigar
  if (arrayMax(tokens) > 2) {
    for (int i = 2; i < arrayMax(tokens) - 1; i += 2) {
      int len = atoi(textItem(tokens, i));
      int intronic = atoi(textItem(tokens, i - 1));
      pos += intronic + 1;
      printf(",%s:%c:%d:%d:%d:%d",
             entry->rname, mrf_strand, pos, pos + len - 1, q, q + len - 1);
      pos += len - 1;
      q += len;
    }
  }

  textDestroy(tokens);
}

int generateSamEntry(Texta tokens, SamEntry *currSamE, int* hasSeqs, 
                     int* hasQual) {
  currSamE->qname = strdup(textItem(tokens, 0));
  currSamE->flags = atoi(textItem(tokens, 1));
  currSamE->rname = strdup(textItem(tokens, 2));
  currSamE->pos   = atoi(textItem(tokens, 3));
  currSamE->mapq  = atoi(textItem(tokens, 4));
  currSamE->cigar = strdup(textItem(tokens, 5));
  currSamE->mrnm  = strdup(textItem(tokens, 6));
  currSamE->mpos  = atoi(textItem(tokens, 7));
  currSamE->isize = atoi(textItem(tokens, 8));
  currSamE->seq   = NULL;
  currSamE->qual  = NULL;
  currSamE->tags  = NULL;

#if 0
  // Debug only
  printf("READ_PAIRED:\t%s\nPAIR_MAPPED:\t%s\nQUERY_UNMAPPED:\t%s\nMATE_UNMAPPED:\t%s\nQUERY_STRAND:\t%s\nMATE_STRAND:\t%s\nFIRST:\t\t%s\nSECOND:\t\t%s\nNOT_PRIMARY:\t%s\nFAILS_CHECK:\t%s\nDUPLICATE:\t%s\n", 
         (currSamE->flags & S_READ_PAIRED)?"paired":"single",
         (currSamE->flags & S_PAIR_MAPPED)?"yes":"no",
         (currSamE->flags & S_QUERY_UNMAPPED)?"unmapped":"mapped",
         (currSamE->flags & S_MATE_UNMAPPED)?"unmapped":"mapped",
         (currSamE->flags & S_QUERY_STRAND)?"forward":"reverse",
         (currSamE->flags & S_MATE_STRAND)?"forward":"reverse",
         (currSamE->flags & S_FIRST)?"yes":"no",
         (currSamE->flags & S_SECOND)?"yes":"no",
         (currSamE->flags & S_NOT_PRIMARY)?"yes":"no",
         (currSamE->flags & S_FAILS_CHECKS)?"yes":"no",
         (currSamE->flags & S_DUPLICATE)?"yes":"no" );
#endif

  // Skip if unmapped or fails platform/vendor checks
  if (currSamE->flags & S_QUERY_UNMAPPED ||
      currSamE->flags & S_MATE_UNMAPPED ||
      currSamE->flags & S_FAILS_CHECKS || 
      currSamE->flags & S_NOT_PRIMARY) {
    return 0;
  }

  // Get tokens
  if (arrayMax(tokens) > 11) {
    Stringa tags = stringCreate(10);
    for (int j = 11; j < arrayMax(tokens); j++) {
      if (j > 11) {
        stringAppendf(tags, "\t");
      }
      stringAppendf(tags, "%s", textItem(tokens, j));
    }
    currSamE->tags = strdup(string(tags));
    stringDestroy(tags);
  }

  if (strcmp(textItem(tokens, 9),  "*") != 0) {
    *hasSeqs = 1;
    currSamE->seq = strdup(textItem(tokens, 9));
  }
  if (strcmp(textItem(tokens, 10), "*") != 0) {
    *hasQual = 1;
    currSamE->qual = strdup(textItem(tokens, 10));
  }
  return 1;
}

void destroySamEntry(SamEntry* currSamE) {
  free(currSamE->qname);
  free(currSamE->rname);
  free(currSamE->cigar);
  free(currSamE->mrnm);
  if (currSamE->seq) {
    free(currSamE->seq);
  }
  if (currSamE->qual) {
    free(currSamE->qual);
  }
  if (currSamE->tags) {
    free(currSamE->tags);
  }
}

int isMateUnmapped(SamEntry* samE) {
  return samE->flags & S_MATE_UNMAPPED;
}

int isPaired(SamEntry* samE) {
  return samE->flags & S_READ_PAIRED;
}

int isValidSamLine(Texta tokens) {
  if (arrayMax(tokens) < 11) {
    textDestroy(tokens);
    return 0;
  } else {
    return 1;
  }
}

int haveSameName(SamEntry *query, SamEntry *mate, char delim) {
  char* pos = strchr(query->qname, delim);
  if (*pos != '\0') {
    *pos='\0';
  }
  pos = strchr(mate->qname, delim);
  if (*pos != '\0') {
    *pos='\0';
  }
  return strEqual(query->qname, mate->qname);
}

int main (int argc, char **argv) {
  char delim = '\0';
  if (argc == 2) {
    delim = argv[1][0];
  }
 
  LineStream ls = ls_createFromFile("-");
  ls_bufferSet(ls, 1);
  char *line = NULL;
  int start = 1;
  while (line = ls_nextLine(ls)) {
    // Put all the lines of the SAM header in comments
    if (line[0] == '@') {
      printf("# %s\n", line);
      continue;
    }
    // Parse each SAM entry and store into array   
    Texta tokens = textFieldtokP(line, "\t");
    if (isValidSamLine(tokens) != 1) {
      ls_destroy(ls);
      die("Invalid SAM entry: %s", line);
    }
 
    SamEntry *currSamE = NULL;
    SamEntry *mateSamE = NULL;
    AllocVar(currSamE);

    int hasQual = 0;
    int hasSeqs = 0;
    int ret = generateSamEntry(tokens, currSamE, &hasSeqs, &hasQual);
    textDestroy(tokens);
    if (ret == 0) {
      if (isPaired(currSamE)) { 
        // only if paired, discard next entry too (the mate) if same name
        line = ls_nextLine(ls); 
        if (ls_isEof(ls)) {
          destroySamEntry(currSamE);
          freeMem(currSamE);
          continue;         
        }
        tokens = textFieldtokP(line, "\t");
        if (isValidSamLine(tokens) != 1) {
          ls_destroy(ls);
          die("Invalid SAM entry: %s", line);
        }
        AllocVar(mateSamE);
        ret = generateSamEntry(tokens, mateSamE, &hasSeqs, &hasQual); 
        // it's part of another pair
        if (!haveSameName(currSamE, mateSamE, delim)) {
          ls_back(ls, 1); // reloading the line
        }
        destroySamEntry(mateSamE);
        freeMem(mateSamE);
      }
      destroySamEntry(currSamE);
      freeMem(currSamE);
      continue;
    }   
    if (isPaired(currSamE)) {
      int hasQual2, hasSeq2;
      AllocVar(mateSamE);
      Texta secondEnd = NULL;
      secondEnd = textFieldtok(ls_nextLine(ls), "\t");
      ret = generateSamEntry(secondEnd, mateSamE, &hasSeq2, &hasQual2);
      textDestroy(secondEnd);
      // This should have been taken care already with the analysis of the 
      // first end.
      if (ret == 0) { 
        destroySamEntry(currSamE);
        destroySamEntry(mateSamE);
        freeMem(currSamE);
        freeMem(mateSamE);
        continue;
      }
      
      if (!haveSameName(currSamE, mateSamE, delim)) {
        die("Please note that for paired-end data, sam2mrf requires the mate "
            "pairs to be on subsequent lines. You may want to sort the SAM "
            "file first.\nEx: sort -r file.sam | sam2mrf > file.mrf\n");
      }
    } 

    // Print MRF headers
    if (start) {
      printf("%s", MRF_COLUMN_NAME_BLOCKS);
      if (hasSeqs) {
        printf("\t%s", MRF_COLUMN_NAME_SEQUENCE);
      }
      if (hasQual) {
        printf("\t%s", MRF_COLUMN_NAME_QUALITY_SCORES);
      }
      printf("\t%s\n", MRF_COLUMN_NAME_QUERY_ID);
      start = 0;
    }
    
    // Print AlignmentBlocks   
    printMrfAlignBlocks(currSamE, R_FIRST);
    if (isPaired(currSamE)) {  
      printf("|");
      printMrfAlignBlocks(mateSamE, R_SECOND);
    }

    seq_init();
    // Print Sequence
    if (hasSeqs) {
      if (!currSamE->seq) {
        die ("Entry missing sequence column\n");
      }
      if (currSamE->flags & S_QUERY_STRAND) {
        seq_reverseComplement(currSamE->seq, strlen(currSamE->seq));
      }
      printf("\t%s", currSamE->seq);
      if (mateSamE) {
        if (!mateSamE->seq) {
          die("Entry missing sequence column\n");
        }
        if (mateSamE->flags & S_MATE_STRAND) {
          seq_reverseComplement(mateSamE->seq, strlen(mateSamE->seq));
        }
        printf("|%s", mateSamE->seq);
      }
    }
    // Print quality scores
    if (hasQual) {
      if (!currSamE->qual) {
        die("Entry missing quality scores column\n");
      }
      printf("\t%s", currSamE->qual);
      if (mateSamE) {
        if (!mateSamE->qual)
          die("Entry missing quality scores column\n");
        printf("|%s", mateSamE->qual);
      }
    }

    // Print queryID
    if (mateSamE) {
      // No need to print out both IDs, but need the pipe symbol for 
      // consistency.
      printf("\t%s|%s", currSamE->qname,"2"); 
    } else {
      printf("\t%s", currSamE->qname);
    }
    printf("\n");
    
    destroySamEntry(currSamE);
    freeMem(currSamE); 
    if (isPaired(currSamE)) {
      destroySamEntry(mateSamE);
      freeMem(mateSamE);
    }
  }
  // clean up
  ls_destroy(ls);
  return EXIT_SUCCESS;
}
