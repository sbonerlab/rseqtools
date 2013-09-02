AM_MAKEFLAGS = --no-print-directory
AM_CPPFLAGS = \
	-I${top_srcdir}/src

bin_PROGRAMS = \
	psl2mrf \
	bowtie2mrf \
	singleExport2mrf \
	mrfSubsetByTargetName \
	mrfQuantifier \
	mrfAnnotationCoverage \
	mrf2wig \
	mrf2gff \
	mrfSampler \
	mrf2bgr \
	wigSegmenter \
	mrfMappingBias \
	mrfSelectRegion \
	mrfSelectSpliced \
	mrfSelectAnnotated \
	createSpliceJunctionLibrary \
	gff2interval \
	export2fastq \
	mergeTranscripts \
	interval2gff \
	interval2sequences \
	intervalSize \
	bed2interval \
	interval2bed \
	mrf2sam \
	sam2mrf \
	mrfValidate \
	bgrQuantifier \
	bgrSegmenter \
	mrfCountRegion \
	mrfSubsetByGenes \
	mrfSorter \
	fixSam	

psl2mrf_SOURCES = src/psl2mrf.c
psl2mrf_LDADD = -lmrf -lbios

bowtie2mrf_SOURCES = src/bowtie2mrf.c
bowtie2mrf_LDADD = -lmrf -lbios

singleExport2mrf_SOURCES = src/singleExport2mrf.c
singleExport2mrf_LDADD = -lmrf -lbios

mrfSubsetByTargetName_SOURCES = src/mrfSubsetByTargetName.c
mrfSubsetByTargetName_LDADD = -lmrf -lbios

mrfQuantifier_SOURCES = src/mrfQuantifier.c
mrfQuantifier_LDADD = -lmrf -lm -lbios

mrfAnnotationCoverage_SOURCES = src/mrfAnnotationCoverage.c
mrfAnnotationCoverage_LDADD = -lmrf -lm -lbios

mrf2wig_SOURCES = src/mrf2wig.c
mrf2wig_LDADD = -lmrf -lm -lbios

mrf2gff_SOURCES = src/mrf2gff.c
mrf2gff_LDADD = -lmrf -lbios

mrfSampler_SOURCES = src/mrfSampler.c
mrfSampler_LDADD = -lmrf -lbios

mrf2bgr_SOURCES = src/mrf2bgr.c
mrf2bgr_LDADD = -lmrf -lm -lbios

wigSegmenter_SOURCES = src/wigSegmenter.c
wigSegmenter_LDADD = -lmrf -lbios

mrfMappingBias_SOURCES = src/mrfMappingBias.c
mrfMappingBias_LDADD = -lmrf -lm -lbios

mrfSelectRegion_SOURCES = src/mrfSelectRegion.c
mrfSelectRegion_LDADD = -lmrf -lm -lbios

mrfSelectSpliced_SOURCES = src/mrfSelectSpliced.c
mrfSelectSpliced_LDADD = -lmrf -lbios

mrfSelectAnnotated_SOURCES = src/mrfSelectAnnotated.c
mrfSelectAnnotated_LDADD = -lmrf -lm -lbios

createSpliceJunctionLibrary_SOURCES = src/createSpliceJunctionLibrary.c
createSpliceJunctionLibrary_LDADD = -lmrf -lbios

gff2interval_SOURCES = src/gff2interval.c
gff2interval_LDADD = -lmrf -lbios

export2fastq_SOURCES = src/export2fastq.c
export2fastq_LDADD = -lmrf -lbios

mergeTranscripts_SOURCES = src/mergeTranscripts.c
mergeTranscripts_LDADD = -lmrf -lm -lbios

interval2gff_SOURCES = src/interval2gff.c
interval2gff_LDADD = -lmrf -lm -lbios

interval2sequences_SOURCES = src/interval2sequences.c
interval2sequences_LDADD = -lmrf -lm -lbios

intervalSize_SOURCES = src/intervalSize.c
intervalSize_LDADD = -lmrf -lm -lbios

bed2interval_SOURCES = src/bed2interval.c
bed2interval_LDADD = -lmrf -lm -lbios

interval2bed_SOURCES = src/interval2bed.c
interval2bed_LDADD = -lmrf -lm -lbios

mrf2sam_SOURCES = src/mrf2sam.c
mrf2sam_LDADD = -lmrf -lm -lbios

sam2mrf_SOURCES = src/sam2mrf.c
sam2mrf_LDADD = -lmrf -lbios

mrfValidate_SOURCES = src/mrfValidate.c
mrfValidate_LDADD = -lmrf -lbios

bgrQuantifier_SOURCES = src/bgrQuantifier.c
bgrQuantifier_LDADD = -lmrf -lbios

bgrSegmenter_SOURCES = src/bgrSegmenter.c
bgrSegmenter_LDADD = -lmrf -lbios

mrfCountRegion_SOURCES = src/mrfCountRegion.c
mrfCountRegion_LDADD = -lmrf -lm -lbios

mrfSubsetByGenes_SOURCES = src/mrfSubsetByGenes.c
mrfSubsetByGenes_LDADD = -lmrf -lm -lbios

mrfSorter_SOURCES = src/mrfSorter.c
mrfSorter_LDADD = -lmrf -lm -lbios

fixSam_SOURCES = src/fixSam.c
fixSam_LDADD = -lmrf -lm -lbios
