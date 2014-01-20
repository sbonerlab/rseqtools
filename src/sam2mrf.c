/// @file sam2mrf.c Module to convert SAM to MRF.

#define _GNU_SOURCE

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <bios/format.h>
#include <bios/log.h>
#include <bios/linestream.h>
#include <bios/common.h>
#include <bios/seq.h>

#include <mrf/mrf.h>
#include <mrf/sam.h>

static void print_mrf_alignment_blocks(SamEntry *sam_entry, int sam_strand) {
  // Parse CIGAR string.
  Array cigar_operations = sam_entry->cigar_ops;
  char mrf_strand = '.'; 

  // Print alignment blocks.
  if( sam_entry->flags & S_QUERY_STRAND ) mrf_strand = '-'; else mrf_strand = '+';
  int target_start = sam_entry->pos;
  int query_start = 1;
  int query_end = query_start;
  bool first = true;
  for (int i = 0; i < arrayMax(cigar_operations); ++i) {
    CigarOperation* op = arrp(cigar_operations, i, CigarOperation);
    switch (op->type) {
      case kCigarAlignmentMatch: {
        int target_end = target_start + op->length - 1;
        query_end = query_start + op->length - 1;
	printf("%s%s:%c:%d:%d:%d:%d",
	       (first == false) ? "," : "",
	       sam_entry->rname, mrf_strand, target_start, target_end,
	       query_start, query_end);
	first = false;
        target_start += op->length;
        query_start += op->length ;
        break;
      }
      case kCigarSkippedRegion: {
        target_start += op->length;
        break;
      }
      case kCigarSoftClipping: {
	if ( i == 0 )  query_start += op->length; else query_end = query_start ;
	break;
      }
    }
    
  }
  
}

void print_headers(bool has_seqs, bool has_qual) {
  printf("%s", MRF_COLUMN_NAME_BLOCKS);
  if (has_seqs) {
    printf("\t%s", MRF_COLUMN_NAME_SEQUENCE);
  }
  if (has_qual) {
    printf("\t%s", MRF_COLUMN_NAME_QUALITY_SCORES);
  }
  printf("\t%s\n", MRF_COLUMN_NAME_QUERY_ID);
}

int sam_to_mrf(char delim) {
  SamParser* parser = samparser_from_file("-");
  Array comments = samparser_get_comments(parser);
  if( comments != NULL ) {
    for (int i = 0; i < arrayMax(comments); ++i) {
      char* comment = arru(comments, i, char*);
      printf("# %s\n", comment);
    }
  }

  bool first = true;
  bool has_seqs = false;
  bool has_qual = false;
  for (SamEntry* entry = NULL; entry = samparser_next_entry(parser); ) {
    SamEntry* mate_entry = NULL;
    // If entry is not valid, skip it. If it is the first of a pair, then make
    // sure that the mate entry is skipped as well.
    if (!samentry_is_valid(entry)) { // valid = mapped and primary alignment
      if (!samentry_is_paired(entry)) {
        samentry_free(entry);
        continue;
      }
      mate_entry = samparser_next_entry(parser);
      if (mate_entry == NULL) {
        samentry_free(entry);
        break;
      }
      if (!samentry_same_name(entry, mate_entry, delim)) {
        samparser_back(parser, 1);
      }
      samentry_free(entry);
      samentry_free(mate_entry);
      continue;
    }
    
    // If entry is the first of a paired-end, then attempt to read then next 
    // entry and validate that it is the mate of the entry.
    if (samentry_is_paired(entry)) {
      mate_entry = samparser_next_entry(parser);
      if (mate_entry == NULL || !samentry_is_valid(mate_entry)) {
	samentry_free(entry);
	if (mate_entry != NULL) {
	    samentry_free(mate_entry);
	}
        continue;
      }
      if (!samentry_same_name(entry, mate_entry, delim)) {
        die("Please note that for paired-end data, sam2mrf requires the mate "
            "pairs to be on subsequent lines. You may want to sort the SAM "
            "file first.\nEx: sort -r file.sam | sam2mrf > file.mrf\n");
      }
    }
    
    // If this is the first entry, print the MRF headers. Note that we assume that all SAM files have sequence and quality scores
    if (first == true) {
      has_seqs = true; //samentry_has_seqs(entry);
      has_qual = true; //samentry_has_qual(entry);
      print_headers(has_seqs, has_qual);
      first = false;
    }

    // Print MRF AlignmentBlocks.
    print_mrf_alignment_blocks(entry, R_FIRST);
    if (samentry_is_paired(entry)) {
      printf("|");
      print_mrf_alignment_blocks(mate_entry, R_SECOND);
    }

    // Print sequence(s) if available.
    if (has_seqs == true) {
      seq_init();
      if (entry->flags & S_QUERY_STRAND) {
        seq_reverseComplement(entry->seq, strlen(entry->seq));
      }
      printf("\t%s", entry->seq);
      if (mate_entry != NULL) {
        if (mate_entry->flags & S_MATE_STRAND) {
          seq_reverseComplement(mate_entry->seq, strlen(mate_entry->seq));
        }
        printf("|%s", mate_entry->seq);
      }
    }

    // Print quality score(s) if available.
    if (has_qual == true) {
      printf("\t%s", entry->qual);
      if (mate_entry != NULL) {
        printf("|%s", mate_entry->qual);
      }
    }
    
    // Print queryID.
    if (mate_entry == NULL) {
      printf("\t%s", entry->qname);
    } else {
      // We don't need to print out both IDs, but we need the pipe symbol for
      // consistency.
      printf("\t%s|%s", entry->qname, "2");
    }

    printf("\n");
    samentry_free(entry);
    if (mate_entry != NULL) {
      samentry_free(mate_entry);
    }
  }
  samparser_free(parser);
  return EXIT_SUCCESS;
  }

int main(int argc, char **argv) {
  char delim = '\0';
  if (argc == 2) {
    delim = argv[1][0];
  } 
  return sam_to_mrf(delim);
}

