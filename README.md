# ITDFinder
FLT3-ITD mutation detection software

perl ITDFinder.pl --bam bamfile --output outdir --ref FLT3-ITD_hg19 --tmp 1
perl convert2annovar.pl -includeinfo -format vcf4old FLT3-ITD.vcf > FLT3.vcf.avinput; 
perl table_annovar.pl FLT3.vcf.avinput Annovar/humandb -outfile FLT3.vcf -buildver hg19 -protocol refGene,clinvar_20220709,cosmic96 -operation g,f,f -arg ,, --otherinfo -nastring NULL; \
python stat_result.py --infile FLT3.vcf.hg19_multianno.txt --outfile predictQuality.txt
