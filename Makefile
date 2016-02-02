# http://www.ncbi.nlm.nih.gov/sra/?term=elc1
# 
#  Novel RNA polymerase II mutation suppresses transcriptional fidelity and oxidative stress sensitivity in rpb9D yeast
#     Hiroshi Koyama, Tomofumi Ueda, Takahiro Ito and Kazuhisa Sekimizu*
#   says spt4 delta mutants have no error rate increase


#  REQUIRES:
#  SRA toolkit to download from SRA (fastq-dump)
#  samtools
#  bamUtil (trimBam) http://genome.sph.umich.edu/wiki/BamUtil:_trimBam
#       https://github.com/statgen/bamUtil
#  bwa (for alignment)
#

M=matlabrun.pl -local -po -jvm  -m mpileup_count_errors 
MPFILES := $(shell ls *mpileup.gz | perl -pne 's/.gz//')
Q=qsub -b y -V -cwd 
YEASTGFF=/homes/users/lcarey/single_cell_behavior//Data/Yeast/genome.gff
YEASTGENOME=/homes/users/lcarey/single_cell_behavior//Data/Yeast/genome.fasta

clean:
	rm -f dl.[eo]* bwa.[epo]* mp.[oep]*


# 67-69 are wild-type ; 70-72 are elc1 mutants
#  download using fastq-dump from the NCBI SRA-toolkit
download:
	qsub -b y -cwd -V -N dl -l h_vmem=3G ~/single_cell_behavior/src/sratoolkit.2.3.5-2-centos_linux64/bin/fastq-dump --split-files --gzip SRR502967
	qsub -b y -cwd -V -N dl -l h_vmem=3G ~/single_cell_behavior/src/sratoolkit.2.3.5-2-centos_linux64/bin/fastq-dump --split-files --gzip SRR502968
	qsub -b y -cwd -V -N dl -l h_vmem=3G ~/single_cell_behavior/src/sratoolkit.2.3.5-2-centos_linux64/bin/fastq-dump --split-files --gzip SRR502969
	qsub -b y -cwd -V -N dl -l h_vmem=3G ~/single_cell_behavior/src/sratoolkit.2.3.5-2-centos_linux64/bin/fastq-dump --split-files --gzip SRR502970
	qsub -b y -cwd -V -N dl -l h_vmem=3G ~/single_cell_behavior/src/sratoolkit.2.3.5-2-centos_linux64/bin/fastq-dump --split-files --gzip SRR502971
	qsub -b y -cwd -V -N dl -l h_vmem=3G ~/single_cell_behavior/src/sratoolkit.2.3.5-2-centos_linux64/bin/fastq-dump --split-files --gzip SRR502972

_bwa:
	bwa aln -n 2 -t $(NSLOTS) -R ''@RG\\tID:spt\\tSM:bar''  $(YEASTGENOME) $(F)_1.fastq.gz -f $(F).sai
	bwa samse  -r '@RG\tID:'$(F)'\tSM:Elongation\tPL:ILLUMINA' $(YEASTGENOME) $(F).sai $(F)_1.fastq.gz \
	| samtools view -Shu -q 20 - \
	| samtools sort -O bam -m 3G -@ $(NSLOTS) -T .samsort_X_$F -o $(F).bam -
	bam trimBam $(F).bam - 10 | samtools sort -O bam -T .$(F).tmpbam -m 3G -@ $(NSLOTS) - -o $(F).clipped.bam
	samtools index $(F).clipped.bam
	du -sh $(F)_1.fastq.gz $(F).sai  $(F).bam  $(F).clipped.bam

bwa:
	qsub -hold_jid dl -b y -cwd -V -N bwa -pe threaded 3-5  -l h_vmem=4G 'make -k _bwa F=SRR502967'
	qsub -hold_jid dl -b y -cwd -V -N bwa -pe threaded 3-5  -l h_vmem=4G 'make -k _bwa F=SRR502968'
	qsub -hold_jid dl -b y -cwd -V -N bwa -pe threaded 3-5  -l h_vmem=4G 'make -k _bwa F=SRR502969'
	qsub -hold_jid dl -b y -cwd -V -N bwa -pe threaded 3-5  -l h_vmem=4G 'make -k _bwa F=SRR502970'
	qsub -hold_jid dl -b y -cwd -V -N bwa -pe threaded 3-5  -l h_vmem=4G 'make -k _bwa F=SRR502971'
	qsub -hold_jid dl -b y -cwd -V -N bwa -pe threaded 3-5  -l h_vmem=4G 'make -k _bwa F=SRR502972'

mp:
	qsub -hold_jid bwa -b y -cwd -V -N mp -l h_vmem=3G 'make _mp F=SRR502967'
	qsub -hold_jid bwa -b y -cwd -V -N mp -l h_vmem=3G 'make _mp F=SRR502968'
	qsub -hold_jid bwa -b y -cwd -V -N mp -l h_vmem=3G 'make _mp F=SRR502969'
	qsub -hold_jid bwa -b y -cwd -V -N mp -l h_vmem=3G 'make _mp F=SRR502970'
	qsub -hold_jid bwa -b y -cwd -V -N mp -l h_vmem=3G 'make _mp F=SRR502971'
	qsub -hold_jid bwa -b y -cwd -V -N mp -l h_vmem=3G 'make _mp F=SRR502972'

_mp:
	samtools  mpileup -f $(YEASTGENOME) -q 35 -d1000000  -C50 -Q 40 $(F).clipped.bam | gzip > $(F).q35.Q40.clipped.mpileup.gz
	#samtools  mpileup -f $(YEASTGENOME) -q 35 -d1000000  -C50 -Q 40 $(F).bam | gzip > $(F).q35.Q40.mpileup.gz



mpE:
	for S in $(MPFILES); do \
		$Q -hold_jid mp -l h_vmem=3G -N mpE make _mpE S=$$S ;\
	done
	sleep 15
	chgrp -R single_cell_behavior ../
	
O1 = --max_mismatches_per_position=1 --max_percent_mismatches_per_position=10
N1 = _MMPP=1_MPMPP=10_
O4 = --max_mismatches_per_position=4 --max_percent_mismatches_per_position=10
N4 = _MMPP=4_MPMPP=10_
O9 = --max_mismatches_per_position=1000 --max_percent_mismatches_per_position=1000
N9 = _MMPP=1000_MPMPP=1000_
_mpE:
	mpileupCountMismatches.pl $(O1) --mpileup_gzipped_file=$${S} > $${S}$(N1).tab
	mpileupCountMismatches.pl $(O4) --mpileup_gzipped_file=$${S} > $${S}$(N4).tab
	mpileupCountMismatches.pl $(O9) --mpileup_gzipped_file=$${S} > $${S}$(N9).tab


