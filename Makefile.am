AM_MAKEFLAGS = --no-print-directory
AM_CPPFLAGS = \
	-I${top_srcdir}/src

bin_PROGRAMS = \
	src/psl2mrf \
	src/bowtie2mrf \
	src/singleExport2mrf \
	src/mrfSubsetByTargetName \
	src/mrfQuantifier \
	src/mrfAnnotationCoverage \
	src/mrf2wig \
	src/mrf2gff \
	src/mrfSampler \
	src/mrf2bgr \
	src/wigSegmenter \
	src/mrfMappingBias \
	src/mrfSelectRegion \
	src/mrfSelectSpliced \
	src/mrfSelectAnnotated \
	src/createSpliceJunctionLibrary \
	src/gff2interval \
	src/export2fastq \
	src/mergeTranscripts \
	src/interval2gff \
	src/interval2sequences \
	src/intervalSize \
	src/bed2interval \
	src/interval2bed \
	src/mrf2sam \
	src/sam2mrf \
	src/mrfValidate \
	src/bgrQuantifier \
	src/bgrSegmenter \
	src/mrfCountRegion \
	src/mrfSubsetByGenes \
	src/mrfSorter \
	src/fixSam	

src_psl2mrf_SOURCES = src/psl2mrf.c
src_psl2mrf_LDADD = -lmrf -lbios

src_bowtie2mrf_SOURCES = src/bowtie2mrf.c
src_bowtie2mrf_LDADD = -lmrf -lbios

src_singleExport2mrf_SOURCES = src/singleExport2mrf.c
src_singleExport2mrf_LDADD = -lmrf -lbios

src_mrfSubsetByTargetName_SOURCES = src/mrfSubsetByTargetName.c
src_mrfSubsetByTargetName_LDADD = -lmrf -lbios

src_mrfQuantifier_SOURCES = src/mrfQuantifier.c
src_mrfQuantifier_LDADD = -lmrf -lm -lbios

src_mrfAnnotationCoverage_SOURCES = src/mrfAnnotationCoverage.c
src_mrfAnnotationCoverage_LDADD = -lmrf -lm -lbios

src_mrf2wig_SOURCES = src/mrf2wig.c
src_mrf2wig_LDADD = -lmrf -lm -lbios

src_mrf2gff_SOURCES = src/mrf2gff.c
src_mrf2gff_LDADD = -lmrf -lbios

src_mrfSampler_SOURCES = src/mrfSampler.c
src_mrfSampler_LDADD = -lmrf -lbios

src_mrf2bgr_SOURCES = src/mrf2bgr.c
src_mrf2bgr_LDADD = -lmrf -lm -lbios

src_wigSegmenter_SOURCES = src/wigSegmenter.c
src_wigSegmenter_LDADD = -lmrf -lbios

src_mrfMappingBias_SOURCES = src/mrfMappingBias.c
src_mrfMappingBias_LDADD = -lmrf -lm -lbios

src_mrfSelectRegion_SOURCES = src/mrfSelectRegion.c
src_mrfSelectRegion_LDADD = -lmrf -lm -lbios

src_mrfSelectSpliced_SOURCES = src/mrfSelectSpliced.c
src_mrfSelectSpliced_LDADD = -lmrf -lbios

src_mrfSelectAnnotated_SOURCES = src/mrfSelectAnnotated.c
src_mrfSelectAnnotated_LDADD = -lmrf -lm -lbios

src_createSpliceJunctionLibrary_SOURCES = src/createSpliceJunctionLibrary.c
src_createSpliceJunctionLibrary_LDADD = -lmrf -lbios

src_gff2interval_SOURCES = src/gff2interval.c
src_gff2interval_LDADD = -lmrf -lbios

src_export2fastq_SOURCES = src/export2fastq.c
src_export2fastq_LDADD = -lmrf -lbios

src_mergeTranscripts_SOURCES = src/mergeTranscripts.c
src_mergeTranscripts_LDADD = -lmrf -lm -lbios

src_interval2gff_SOURCES = src/interval2gff.c
src_interval2gff_LDADD = -lmrf -lm -lbios

src_interval2sequences_SOURCES = src/interval2sequences.c
src_interval2sequences_LDADD = -lmrf -lm -lbios

src_intervalSize_SOURCES = src/intervalSize.c
src_intervalSize_LDADD = -lmrf -lm -lbios

src_bed2interval_SOURCES = src/bed2interval.c
src_bed2interval_LDADD = -lmrf -lm -lbios

src_interval2bed_SOURCES = src/interval2bed.c
src_interval2bed_LDADD = -lmrf -lm -lbios

src_mrf2sam_SOURCES = src/mrf2sam.c
src_mrf2sam_LDADD = -lmrf -lm -lbios

src_sam2mrf_SOURCES = src/sam2mrf.c
src_sam2mrf_LDADD = -lmrf -lbios

src_mrfValidate_SOURCES = src/mrfValidate.c
src_mrfValidate_LDADD = -lmrf -lbios

src_bgrQuantifier_SOURCES = src/bgrQuantifier.c
src_bgrQuantifier_LDADD = -lmrf -lbios

src_bgrSegmenter_SOURCES = src/bgrSegmenter.c
src_bgrSegmenter_LDADD = -lmrf -lbios

src_mrfCountRegion_SOURCES = src/mrfCountRegion.c
src_mrfCountRegion_LDADD = -lmrf -lm -lbios

src_mrfSubsetByGenes_SOURCES = src/mrfSubsetByGenes.c
src_mrfSubsetByGenes_LDADD = -lmrf -lm -lbios

src_mrfSorter_SOURCES = src/mrfSorter.c
src_mrfSorter_LDADD = -lmrf -lm -lbios

src_fixSam_SOURCES = src/fixSam.c
src_fixSam_LDADD = -lmrf -lm -lbios
