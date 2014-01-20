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
 * @file fixSam.c Program to fix SAM files
 * @pre The sam file should have already excluded unmapped and non primary 
 *      alignments and also be sorted by read name.
 */

void fix_single_end(SamEntry *entry) {
  if ((entry->flags & S_QUERY_UNMAPPED) &&
      !(strEqual(entry->rname, "*"))) {
    entry->flags -= S_QUERY_UNMAPPED;
  }
  if ((entry->flags & S_MATE_UNMAPPED) && 
      !(strEqual(entry->mrnm, "*"))) {
    entry->flags -= S_MATE_UNMAPPED;
  }
                                               
  if (strEqual(entry->mrnm, "*") && 
      !(entry->flags & S_MATE_UNMAPPED)) {
    entry->flags += S_MATE_UNMAPPED;
    if (entry->flags & S_PAIR_MAPPED) {
      entry->flags -= S_PAIR_MAPPED;
    }
    if (entry->flags & S_MATE_STRAND) {
      entry->flags -= S_MATE_STRAND;
    }
  }
  // removing /1 or /2 in the read name
  char *pos = strchr(entry->qname, '/');
  if (pos != NULL) {
    if ((entry->flags & S_FIRST) && 
        (entry->flags & S_SECOND)) {
      if (strEqual(pos, "/1")) {
        entry->flags -= S_SECOND ;
      } else {
        entry->flags -= S_FIRST ;
      }
    }
    if (!(entry->flags & S_FIRST) && 
        !(entry->flags & S_SECOND)) {
      if (strEqual( pos, "/1")) { 
        entry->flags -= S_SECOND;
      } else {
        entry->flags -= S_FIRST;
      }
      *pos='\0';
    }
  }
  if ((entry->flags & S_QUERY_UNMAPPED) && 
      (entry->flags & S_NOT_PRIMARY)) {
    entry->flags -= S_NOT_PRIMARY; 
  }
}

void fix_paired_ends(SamEntry *r1, SamEntry *r2) {
  if ((r1->flags & S_MATE_UNMAPPED) && !(r2->flags & S_QUERY_UNMAPPED)) {
    r1->flags -= S_MATE_UNMAPPED;
    r1->mrnm = hlr_strdup(r2->rname);
    r2->mrnm = hlr_strdup(r1->rname);
    r1->mpos = r2->pos;
    r2->mpos = r1->pos;
  } 
 
  if (!(r1->flags & S_QUERY_UNMAPPED) && (r2->flags & S_MATE_UNMAPPED)) {
    r2->flags -= S_MATE_UNMAPPED;
    r1->mrnm = hlr_strdup(r2->rname);
    r2->mrnm = hlr_strdup(r1->rname);
    r1->mpos = r2->pos;
    r2->mpos = r1->pos;
  }

  if ((((r1->flags & S_MATE_STRAND) > 0) != ((r2->flags & S_QUERY_STRAND) > 0))) {
    if (r2->flags & S_QUERY_STRAND) {
      r1->flags += S_MATE_STRAND;
    } else {
      r1->flags -= S_MATE_STRAND;
    }
  }
  if ((((r2->flags & S_MATE_STRAND) > 0) != ((r1->flags & S_QUERY_STRAND) > 0))) { 
    if (r1->flags & S_QUERY_STRAND) {
      r2->flags += S_MATE_STRAND;
    } else {
      r2->flags -= S_MATE_STRAND;
    }
  }

  if (strEqual(r1->rname, r2->rname)) {
    r1->isize = r1->pos - r2->pos;
    r2->isize = -1 * r1->isize;
  }
  if (!(r1->flags & S_READ_PAIRED)) {
    r1->flags += S_READ_PAIRED;
  }
  if (!(r2->flags & S_READ_PAIRED)) {
    r2->flags += S_READ_PAIRED;
  }
  if ((r1->flags & S_FIRST) && !(r2->flags & S_SECOND)) {
    r2->flags += S_SECOND;
  }
  if ((r1->flags & S_SECOND) && !(r2->flags & S_FIRST)) {
    r2->flags += S_FIRST;
  }
  if ((r2->flags & S_FIRST) && !(r1->flags & S_SECOND)) {
    r1->flags += S_SECOND;
  }
  if ((r2->flags & S_SECOND) && !(r1->flags & S_FIRST)) {
    r1->flags += S_FIRST; 
  }
}

int main (int argc, char **argv) {
  SamEntry *mate_entry = NULL;
  SamEntry *tmp_entry = NULL;
  SamEntry *r1 = NULL;
  SamEntry *r2 = NULL;
  SamEntry *r3 = NULL;

  SamParser* parser = samparser_from_file("-");
  for (SamEntry* tmp_entry = NULL; tmp_entry = samparser_next_entry(parser); ) {
    //samentry_free(tmp_entry);
    SamEntry* entry = NULL;
   
    samentry_copy(&entry, tmp_entry);
    fix_single_end(entry);
    tmp_entry = samparser_next_entry(parser);
    if (tmp_entry == NULL) {
      break;
    }
    samentry_free(mate_entry);
    mate_entry = NULL; 
    
    samentry_copy(&mate_entry, tmp_entry);
    fix_single_end(mate_entry);
    if (strEqual(entry->qname, mate_entry->qname)) {
      fix_paired_ends(entry, mate_entry);
      printf("%s\n", samentry_to_string(entry));
      printf("%s\n", samentry_to_string(mate_entry));
      if (r1) {
        samentry_free(r1); 
        r1 = NULL;
      }
      if (r2) {
        samentry_free(r2); 
        r2 = NULL;
      }
      if (r3) {
        samentry_free(r3); 
        r3 = NULL;
      }
      continue;
    } 
    if (r1 == NULL) {
      samentry_copy(&r1, mate_entry);
      if (r3) {
        if (strEqual(r3->qname, entry->qname)) {
          fix_paired_ends(r3, entry);
          printf("%s\n", samentry_to_string(r3));
          printf("%s\n", samentry_to_string(entry));
          samentry_free(r3);
          r3 = NULL;
          continue;
        } else {
          samentry_free(r3); 
          r3 = NULL;
        }
      }
      continue;
    } else if (r2 == NULL) {
      samentry_copy(&r2, entry);
      samentry_copy(&r3, mate_entry);
    }    
    if (r1->qname != NULL && r2->qname != NULL) {
      if (strEqual(r1->qname, r2->qname)) {
        fix_paired_ends(r1, r2);
        printf("%s\n", samentry_to_string(r1));
        printf("%s\n", samentry_to_string(r2));
        samentry_free(r1);
        r1 = NULL;
        samentry_free(r2);
        r2 = NULL;
      }
    }   
  }
  if (r1 != NULL) { 
    samentry_free(r1); 
    r1 = NULL;
  }
  if (r2) {
    samentry_free(r2); 
    r2 = NULL;
  }
  samparser_free(parser);
  return EXIT_SUCCESS;
}
