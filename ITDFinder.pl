#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper qw(Dumper);
use List::Util qw(max);
use List::Util qw(min);
use List::Util qw(sum);
use POSIX qw(ceil);
use POSIX qw(floor);
use Math::Combinatorics qw(combine);

use Getopt::Long qw(GetOptions);


my $version = '1.0';
my $update_data = '2021-11-02';


my $logo='
############################################################################################
#ITDFinder usage:
############################################################################################';

my $usage = <<__EOUSAGE__;

$logo
#
#  Required:
#
#  To include running acornHMAS_ITD:
#
#      --bam <string>                    BWA-MEM BAM file
#
#      --output <string>                 Output path
#
#      --ref <string>                    reference genome in FASTA format (with fai index)
#

#  Optional:
#
#    --tmp <int>                         file for temporary files,chose 1 (default: 0)
#
#    --version <string>               	 report version ($version)
#
#
#######################
#
#   Quick guide to running:
#
#       perl ITDFinder.pl --bam bam file --ref ref file --tmp 1 --output ouddir
#
############################################################################################

__EOUSAGE__

    ;

my $help_flag; my $inputbam = ''; my $outpath =''; my $ref ='';my $tmp = '';my $REPORT_VERSION = 0;
&GetOptions ( 'help|h' => \$help_flag,
              'bam=s' => \$inputbam,
              'output=s' => \$outpath,
              'ref=s' => \$ref,
              'version|v' => \$REPORT_VERSION,
              'tmp|t=s' => \$tmp,
    );

if ($help_flag) {
    die $usage;
}
if ($REPORT_VERSION) {
    print "\n\nITDFinder version: $version\n\n";
    exit(0);
}

unless ($inputbam) {
    die "Error, need --bam BWA-MEM BAM file\n";
}

unless ($outpath) {
    die "Error, need --output Output path\n";
}

unless ($ref) {
    die "Error, need --ref reference genome in FASTA format (with fai index)\n";
}

## 软件路径设置
my $clustercmd = "usearch";
my $samtoolscmd = "samtools";
my $bamToFastqcmd = "bamToFastq";
my $bwacmd = "bwa";
my $python ='python';
my $similarityfile = 'similarity.py';

my @line_reads; my @sub_reads; my @a;
my $key; my $readname; my $cigar; my $seq; my $qa; my $i; my $j; my $k;
my $outbase; my $outdir;my $fastq1 = '';my $fastq2 = '';my $mappingprefix  = "remap";my $targetbed = "";
my $expand = 9; my $refstartpos = 28607438; my $refendpos = 28608937;
my $refseq = "TGAATTTTTTTAGAGACGGGATTTCACTGTGTTGCCTGGGCTGGTGTCTAATTCTTGGACTCAAGTGATCCTCCCATCTTGGCCTCCCAAAGTGCTGGGACTACAGGCGTGAGCCATCGTGCCTGGACCAGATGTAATTCTTGTAAATTTACCCTAGTAAATGGTTTTAAAGAATGGATTGAGAATAATCCGCAATTTTCTAGGGAGGAAATGTATCCTTTACAGGGACATTGCCTGATTGTCTGTGGGGAGTAGAGTATATGTAGAGTGGTTGTTAGGACTGAAAATGATTATTACTGAAACAGGATGTGAGAGATTATAATGAGTTGTCCACTATTTATAATGTCACACAGGAATTCTGTTTCATCGCTGAGTGACACTCTTTTGTTGCAGGCCCCTTCCCTTTCATCCAAGACAACATCTCATTCTATGCAACAATTGGTGTTTGTCTCCTCTTCATTGTCGTTTTAACCCTGCTAATTTGTCACAAGTACAAAAAGGTAAAAGCAAAGGTAAAAATTCATTATTCTTTCCTCTATCTGCAGAACTGCCTATTCCTAACTGACTCATCATTTCATCTCTGAAGCAATTTAGGTATGAAAGCCAGCTACAGATGGTACAGGTGACCGGCTCCTCAGATAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGATCTCAAATGGGAGTTTCCAAGAGAAAATTTAGAGTTTGGTAAGAATGGAATGTGCCAAATGTTTCTGCAGCATTTCTTTTCCATTGGAAAATCTTTAAAATGCACGTACTCACCATTTGTCTTTGCAGGGAAGGTACTAGGATCAGGTGCTTTTGGAAAAGTGATGAACGCAACAGCTTATGGAATTAGCAAAACAGGAGTCTCAATCCAGGTTGCCGTCAAAATGCTGAAAGGTACAGTATAGTGGAAGGACAGCAACAAAGATGCACAAAAATGGGAGGCACAGTTTCCCACCCATGCCTTCTTCTCTTTTCCATCCTTTTAATGGTTACTGTTTGCCATGTTTCAAGGCTAAAAATGGAGTGGATTGGGGTGTCAACCCAGTCATGAATACAAAATTCAAGTCATAAATACAAAACTCCACCTCTTCTTCCACTTCTCTCTTTGTTTATTTTTATTTTCTTTTATTTTTTTGAGACAGGGTCTCACTCTGTCACCTAGGCTGGAGTACAGTAGTGTGATCATGGCTCACTGCAGCCTCCACCCCCTGGGCTCAGGCTCCTCAGCAGCTGGGGCTACAGGCGTGCACCACCACGCTCAGCTAAGTTTTTTATTTTTGGTAGAGACGGGGTTTCACCATGTTGCCTAGGCTGATCTGGGACTCCTGAGCTCAAACTATCTGCCCGCCTCGGCCTCCCAAAGTGTTGGGATTACAGGCGTGAGCCACTGTGCCAGGCCTATTTAGTTTTTATAGGCTTTGATATTTGGGGATACAGTTGCCACAGGCAGAAGTTCTTTCTTACATTTGTGTTGTTTTAT";


########################第一步：寻找候选ITD##################################################
@line_reads = split(/\//, $outpath);
$outdir = $line_reads[scalar(@line_reads)-1]; # 180563978B3
$outbase = $outpath . "/" . $outdir;
$fastq1 = $outbase . ".R1.fastq";
$fastq2 = $outbase . ".R2.fastq";

if( system( sprintf( "%s view -H %s | grep -w 'SN:chr13\\|SN:13' > %s_header.txt", $samtoolscmd, $inputbam, $outbase ) ) ) {
  die "Error,Failed to extract chr13 from bam.\n";}
open(FI, $outbase . "_header.txt" ) or die $!;
while(<FI>) { if($_ =~ /SN:(\w+)/) { $targetbed = $1 . ":" . $refstartpos .  "-" . $refendpos; } }
close FI;
if( system( sprintf( "%s view -b -F 0x900 %s %s > %s.bam", $samtoolscmd, $inputbam, $targetbed, $outbase ) ) ||
	system( sprintf( "%s index %s.bam", $samtoolscmd, $outbase ) ) ||
	system( sprintf( "%s sort -n %s.bam -o %s.sorted.bam", $samtoolscmd, $outbase, $outbase ) ) ) {
  die "Error,Failed to extract FLT3-ITD locus from input bam.\n";
}
if( system( sprintf( "%s -i %s.sorted.bam -fq %s -fq2 %s",
	 $bamToFastqcmd, $outbase, $fastq1, $fastq2 ) ) ) {
die "Error,Failed to use bamToFastq.\n";
}
if( system( sprintf( "rm %s_header.txt", $outbase )) ) {
  die "Error,Failed to delete $outbase._header.txt.\n";}
if( system( sprintf( "rm %s.bam", $outbase )) ||
system( sprintf( "rm %s.bam.bai", $outbase )) ) {
  die "Error,Failed to delete $outbase.bam and $outbase.bam.bai .\n";
}
if( system( sprintf( "rm %s.sorted.bam", $outbase )) ) {
die "Error,Failed to delete $outbase.sorted.bam.\n";
}

# 寻找候选ITD的局部重新比对
if( system( sprintf( "%s mem -k 6 -M -O 6 -T 12 %s %s %s | grep FLT3-ITD_hg19 > %s_%s.sam",
	$bwacmd, $ref, $fastq1, $fastq2, $outbase, $mappingprefix ) ) ) {
  die "Error,Failed to Local alignments.\n";
}

my $read1orread2; my $strand;
my %acornFwdR2 = (); my %acornRevR2 = ();
my %acornprimarysam = (); my %acornsecondary = (); my %acornpriseq = (); my %acornpriqa = (); my %acornpristr = (); my %acornspos = ();
my %acornpricig = (); my %acornseccig = (); my %acornpriloc = (); my %acornsecloc = (); my %acornsecsstr = (); my %acornpripos = ();
my %acorndefineWT = (); my %acorndefineWTPaired = (); my %acornreadsPassed = (); my %pbase = ();
open( FI, $outbase . "_" . $mappingprefix . ".sam" ) or die $!;
while(<FI>){
    if( !/^@/ ) {
	@line_reads = split( /\t/, $_ );

	if( $line_reads[2] eq "*" ) { next; }

	if( ( $line_reads[1] & 64 ) == 64 ) {$read1orread2 = "_1";
	} elsif( ( $line_reads[1] & 128 ) == 128 ) {$read1orread2 = "_2";
	} else {next;}

	$readname = $line_reads[0] . $read1orread2;

	$cigar = $line_reads[5];
	if( $cigar =~ /^(\d+)H(.*)H$/ ) { next; }
	# For Archer, assume GSP2s are on R2
	# if( $read1orread2 eq "_2" ) { $readlengths{length( $line_reads[9] )} += 1; }
	my $acornrpos = $line_reads[3]-1;
	while( $cigar =~ /(\d+)([SHDIMN])(.*)/ ) {
	    $cigar = $3;
	    if( $2 eq "M" || $2 eq "N" || $2 eq "D" ) { $acornrpos += $1; }
	}

	if( ( $line_reads[1] & 16 ) == 16 ) { $strand = "-"; } else { $strand = "+"; }

	if( ( $line_reads[1] & 256 ) == 256 ) {
	    $acornsecondary{$readname} = $_;
	    $acornsecloc{$readname} = $line_reads[3]-1;
	    $acornseccig{$readname} = $line_reads[5];
	    $acornspos{$readname} = $acornrpos;
	    $acornsecsstr{$readname} = $strand;
	} elsif( ( $line_reads[1] & 2048 ) != 2048 ) {
	    $pbase{$line_reads[0]} = 1;
	    $acornprimarysam{$readname} = $_;
	    $acornpriloc{$readname} = $line_reads[3]-1;
	    $acornpricig{$readname} = $line_reads[5];
	    $acornpriseq{$readname} = $line_reads[9];
	    $acornpriqa{$readname} = $line_reads[10];
	    $acornpripos{$readname} = $acornrpos;
	    $acornpristr{$readname} = $strand;

	    if( $read1orread2 eq "_2" ) {
		if( ( $line_reads[1] & 16 ) == 0 ) {
		  $acornFwdR2{$line_reads[3]} += 1;
		} else {
		  $acornRevR2{$acornrpos} += 1;
		}
	    }
	}
    }
}
close FI;

if( !$tmp && system( "rm " . $outbase . "_" . $mappingprefix . ".sam" ) ) {
  die "Error removing local alignment file to FLT3 target locus. Exiting...";
}

my %acornlloc = (); my %acornrloc = ();
foreach my $read ( keys %acornprimarysam ) {
    # my $extSeq = ""; my $extLoc = "";
    if( exists( $acornsecondary{$read} ) ) {
		#有次要比对位置的直接加入$acornlloc，$acornrloc里候选ITDreads 集合。
	if( $acornpricig{$read} =~ /^(\d+)(S)(.*)/ && $acornpricig{$read} =~ /(\d+)(M)$/ && $acornseccig{$read} =~ /^(\d+)(M)(.*)/ ) {
	    $acornlloc{$read} = $acornsecloc{$read};
	    $acornrloc{$read} = $acornpripos{$read};
	} elsif( $acornpricig{$read} =~ /(\d+)(S)$/ && $acornpricig{$read} =~ /^(\d+)(M)(.*)/ && $acornseccig{$read} =~ /(\d+)(M)$/ ) {
	    $acornlloc{$read} = $acornpriloc{$read};
	    $acornrloc{$read} = $acornspos{$read};
	} else {
	    next;
	}
    } else {
		# 没有次要比对位置的，把净插入大于3bp也加入$acornlloc，$acornrloc里候选ITDreads 集合。novel ITD的出现。和elsif 就是野生型的判定
	$cigar = $acornpricig{$read};
	my $acornclipLength = 0;
	my $acorninsLength = 0;
	my $acorndelLength = 0;
	if( @a = $cigar =~ /(\d+)S/g ) { $acornclipLength = max(@a); }
	if( @a = $cigar =~ /(\d+)I/g ) { foreach my $x (@a) { $acorninsLength += $x; } }
	if( @a = $cigar =~ /(\d+)D/g ) { foreach my $x (@a) { $acorndelLength += $x; } }

	if( $acorninsLength-$acorndelLength >= 3 && $acornclipLength == 0 ) {  # net insertion of 3+ bp
	    $acornlloc{$read} = $acornpriloc{$read};
	    $acornrloc{$read} = $acornpripos{$read};
	} elsif( $acorninsLength <= 2 && $acorndelLength <= 2 && $acornclipLength <= 3 && ($acornprimarysam{$read} =~ /NM:i:(\d+)/ && $1 <= 5) ) {
	    $acorndefineWT{$read} = 1;  # reads aligning "stringently" to reference
	}
    }
}

foreach my $read (keys %acorndefineWT) {
    $read =~ s/\_1$//g; $read =~ s/\_2$//g;
    if( exists($acorndefineWT{$read."_1"}) && exists($acorndefineWT{$read."_2"}) && $acornpristr{$read."_1"} ne $acornpristr{$read."_2"} ) {
	$acorndefineWTPaired{$read} = 1;
	$acornreadsPassed{$read} = 1;
    }
}

my $acornfatsa = "";
open(FO1, ">" . $outbase . ".nowild.R1.fastq" ) or die $!;
open(FO2, ">" . $outbase . ".nowild.R2.fastq" ) or die $!;
foreach my $read (sort keys %pbase ) {
    if( exists($acornprimarysam{$read . "_1"}) && exists($acornprimarysam{$read . "_2"}) ) {
	@line_reads = split(/\t/, $acornprimarysam{$read . "_1"});
	$seq = $line_reads[9]; $qa = $line_reads[10];
       	if( ($line_reads[1] & 16) == 16 ) { $seq = reverse $seq; $seq =~ tr/ATGCatgc/TACGtacg/; $qa = reverse $qa; }
	$acornfatsa = "@" . $read . "/1\n" . $seq . "\n" . "+" . "\n" . $qa . "\n";
        if( !exists($acorndefineWTPaired{$read}) ) { print FO1 $acornfatsa; }
	@line_reads = split(/\t/, $acornprimarysam{$read . "_2"});
	$seq = $line_reads[9]; $qa = $line_reads[10];
       	if( ($line_reads[1] & 16) == 16 ) { $seq = reverse $seq; $seq =~ tr/ATGCatgc/TACGtacg/; $qa = reverse $qa; }
	$acornfatsa = "@" . $read . "/2\n" . $seq . "\n" . "+" . "\n" . $qa . "\n";
	if( !exists($acorndefineWTPaired{$read}) ) { print FO2 $acornfatsa; }
    }
}
close FO1;
close FO2;

######################################第三步：聚类和分组##########################################################
my $acornSeq = ''; my %acornByLength = (); my %acornByseq = ();
foreach my $read ( keys %acornlloc ) {
  my $acornitdsize = length($acornpriseq{$read}) - $acornrloc{$read} + $acornlloc{$read};
	###ITD长度的上限就是500bp。
  if( $acornitdsize > 0 && $acornitdsize <= 500 ) {
    $acornSeq = substr($refseq,0,$acornlloc{$read}) . $acornpriseq{$read} . substr($refseq,$acornrloc{$read});
    $acornByLength{$acornitdsize} += 1;
    $acornByseq{$acornitdsize}{$acornSeq} += 1;
  }
}
#'''$acornSeq 为理解的也是ref + ITD 长度的序列但不一定是ITD序列
# 打印的%acornByLength= {
#           '66' => 15,  # ITD长度为66的有多少条reads支持，整理多少个支持这个ITD的突变型reads
#           '45' => 381,
#           '21' => 639
#         };
# 打印的%acornByseq={   # 具体见软件下文件
#           '66' => {
#                     'TGAATTTTTTTAGAGACGGGATTTCACTGTGTTGCCTGGGCTGGTGTCTAATTCTTGGACTCAAGTGATCCTCCCATCTTGGCCTC
# CCAAAGTGCTGGGACTACAGGCGTGA'=> 15},
#           '21' => { 多个序列=>} => 639,
# 		  '45' => { 多个序列=>}=> 381,
# '''

my %usearchHash = (); my $usearchclust = ""; my %acornsizeHash = (); my %acornpreHash = ();
$i = 0; # 按照频率进行排序后整理在去聚类
foreach $j ( sort {$acornByLength{$b}<=>$acornByLength{$a}} keys %acornByLength ) {
    foreach ( sort{ $acornByseq{$j}{$b}<=>$acornByseq{$j}{$a} } keys %{$acornByseq{$j}} ) {
	if( $acornByseq{$j}{$_} >= 1 ){
	    for( $k=0; $k<$acornByseq{$j}{$_}; $k++ ) {
		$i += 1;
		$key = $i . "_length" . $j . "_size=" . $acornByseq{$j}{$_}; # size后面接的是这个ITD长度，这个序列的具体支持数
		$usearchHash{$key} = $_; # $_ 为具体序列：ref+ITD序列
		$acornsizeHash{$key} = length($_);
		$acornpreHash{$key} = $acornByseq{$j}{$_}; # num of initial reads generating the extended read
		$usearchclust .= ">" . $key . "\n" . $_ . "\n";
	    }
	}
    }
}

if( length($usearchclust) == 0 ) {
    print "no FLT3-ITD. \n";
    system( "touch " . $outbase . "_FLT3-ITD_negative.vcf" );
    if( !$tmp && ( system( "rm " . $outbase . ".nowild.R*.fastq" ) ||
        ( $inputbam ne "" && system( "rm " . $outbase . ".R*.fastq" ) ) ) ) {
      die "Error removing nowild.fastq.";
    }
    exit(0);
}

open (FO, ">" . $outbase . "_for_usearchclust.fasta" ) or die $!;
print FO $usearchclust;
close FO;

if( system( $clustercmd . " -cluster_fast " . $outbase . "_for_usearchclust.fasta -id 0.997 -sort size -centroids " . $outbase . "_clusters.fasta" .  " -uc " . $outbase . "_clusters.uc" ) ) {
  die "usearch command failed. Please find out the reason.\n";
}
###把聚类后的ITD加上序列信息，序列上是ref+ITD序列
open (FI, $outbase . "_clusters.fasta" ) or die $!;
open (FO, ">" . $outbase . "_candidate_clusters.fa" ) or die $!;
# my $lines_keys = ""; ——》这是声明局部变量。 our $lines_keys = "";——》这是声明全部变量。第一次出现的变量一定要指定是局部还是全局变量，不能不指定就使用，第二次使用的时候可以直接用，不用接my,our等
my $lines_keys = "";
while(<FI>) {
    chomp;
	if($_ =~ /^>/) {
		print FO $_ . "\n";
		$lines_keys = substr $_, 1;
		print FO $usearchHash{$lines_keys} . "\n";
	}
}
close FI;
close FO;

# Alignment of candidate ITD clusters against reference
if( system( sprintf( "%s mem -M -O 6 -T 20 %s %s_candidate_clusters.fa > %s_candidate_clusters_%s.sam",
	$bwacmd, $ref , $outbase, $outbase, $mappingprefix ) ) ) {
  die "Error，Fail to aligning FLT3 target region.";
}

my %acornHashinfos = ();
open (FI, $outbase . "_candidate_clusters_" . $mappingprefix . ".sam") or die $!;
while(<FI>) {
    if( !/^@/ ) {
	@line_reads = split( /\t/, $_ );

	my %acornerroacornrlocs = (); my $acorneacornrpos = $line_reads[3]-1;
	my $mdstring = substr($line_reads[12],5); #strip off MD:Z:
	while( $mdstring =~ /^(\d+)([ACGT]|\^[ACGT]+)(.*)/ ) {
	    $mdstring = $3;
	    if( substr( $2, 0, 1 ) eq "^" ) {
		$acorneacornrpos += $1 + length($2) - 1;
	    } else {
		$acorneacornrpos += $1 + 1; $acornerroacornrlocs{$acorneacornrpos} = $2;
	    }
	}

	my $acornspos = 1; my $acornrpos = $line_reads[3]; my $acornhpos = 0; my $i=0;
	$cigar = $line_reads[5];
	print $cigar . "\n";
	while( $cigar =~ /(\d+)([SHDIMN])(.*)/ ) {
	    my $mmkey = "";
	    $cigar = $3;
	    if( $2 eq "M" ) {
		if( $i==0 || $cigar eq "" || $1 >= 12 ) {
		    foreach $acorneacornrpos ( sort {$a<=>$b} keys %acornerroacornrlocs ) {
			if( $acorneacornrpos >= $acornrpos && $acorneacornrpos < $acornrpos+$1 ) {
			    if( length($mmkey) > 0 ) { $mmkey .= "+"; }
			    $mmkey .= $acorneacornrpos . ":" . $acornerroacornrlocs{$acorneacornrpos} . ">" .
				substr( $line_reads[9], $acorneacornrpos+$acornspos-$acornrpos-$acornhpos-1, 1);
			}
		    }
		    $acornHashinfos{$line_reads[0]}{$acornspos} = {r0=>$acornrpos+0, r1=>$acornrpos+$1-1, s0=>$acornspos+0, s1=>$acornspos+$1-1, mm=>$mmkey };
		}
		$acornspos += $1; $acornrpos += $1;
	    } elsif( $2 eq "I" || $2 eq "S" || $2 eq "H" ) {
		$acornspos += $1;
		if( $2 eq "H" ) { $acornhpos += $1; }
	    } elsif( $2 eq "D" ) {
		$acornrpos += $1;
	    }
	    $i++;
	}
    }
}
close FI;
# 打印的有错配的例子：%acornHashinfos  就是以前的%acornHashinfoss
# $VAR1 = {
#           '1_length33_size=590' => {
#                                      '1' => {
#                                               's1' => 673,
#                                               's0' => 1,
#                                               'r1' => 673,
#                                               'mm' => '',
#                                               'r0' => 1
#                                             },
#                                      '707' => {
#                                                 'mm' => '',
#                                                 'r1' => 1500,
#                                                 'r0' => 674,
#                                                 's1' => 1533,
#                                                 's0' => 707
#                                               }
#                                    },
#           '628_length33_size=1' => {
#                                      '713' => {
#                                                 's0' => 713,
#                                                 's1' => 1533,
#                                                 'r0' => 680,
#                                                 'mm' => '710:T>G+727:T>G+756:T>A+765:T>G+773:T>A',
#                                                 'r1' => 1500
#                                               },
#                                      '1' => {
#                                               's1' => 679,
#                                               's0' => 1,
#                                               'mm' => '',
#                                               'r1' => 679,
#                                               'r0' => 1
#                                             }
#                                    },
#           '631_length33_size=1' => {
#                                      '1' => {
#                                               'r1' => 692,
#                                               'mm' => '611:C>A+617:G>C+629:G>T+631:C>A+644:G>T+649:C>A',
#                                               'r0' => 1,
#                                               's1' => 692,
#                                               's0' => 1
#                                             },
#                                      '726' => {
#                                                 's0' => 726,
#                                                 's1' => 1533,
#                                                 'r0' => 693,
#                                                 'mm' => '704:G>T',
#                                                 'r1' => 1500
#                                               }
#                                    }
#         };
# '''
if( !$tmp && (
    system( "rm " . $outbase . "_for_usearchclust.fasta" ) ||
    system( "rm " . $outbase . "_clusters.fasta" ) ||
    system( "rm " . $outbase . "_clusters.uc" ) ||
	system( "rm " . $outbase . "_candidate_clusters.fa" ) ||
    system( "rm " . $outbase . "_candidate_clusters_" . $mappingprefix . ".sam" ) ) ) {
    die "Error,Fail to removing usearch files.";
}

open (FO, ">" . $outbase . "_candidate_clusters_middle.fa") or die $!; # 重新写一个新的文件。里面是ITD序列，是在负链上的实际reads序列和参考基因组反向互补。
foreach( keys %acornHashinfos ) {
    my @acornkeys = sort {$a<=>$b} keys %{$acornHashinfos{$_}};
    for( my $i=1; $i<scalar(@acornkeys); $i++ ) {
	if( ($acornHashinfos{$_}{$acornkeys[$i]}{s0}-$acornHashinfos{$_}{$acornkeys[$i-1]}{s1}-1) >= 10 ) {
		#ITD的终止位置-ITD起始位置-1 >= 10bp 如果不是不能证明有ITD插入，之前测试的时候一个是不符合>=10Bp,cigar值看着不像ITD
	    print FO ">" . $_ . "|start" . ($i-1) . "\n" .
		substr( $usearchHash{$_}, $acornHashinfos{$_}{$acornkeys[$i-1]}{s1},
			$acornHashinfos{$_}{$acornkeys[$i]}{s0}-$acornHashinfos{$_}{$acornkeys[$i-1]}{s1}-1 ) . "\n"; # #截取$extSeq 中 开始为：ITD起始位置，长度为ITD的长度。应该是最终报出ITD的反向互补序列，或者复杂情况就不是ITD序列
	}
    }
}
close FO;

#将ITD候选序列比对到参考基因组上
if( system( sprintf( "%s mem -k 10 -M -O 6 -T 10 %s %s_candidate_clusters_middle.fa > %s_candidate_clusters_middle_%s.sam",
	$bwacmd, $ref , $outbase, $outbase, $mappingprefix ) ) ) {
  die "Error,Fail to aligning FLT3 target .";
}

open (FI, $outbase . "_candidate_clusters_middle_" . $mappingprefix . ".sam") or die $!; # # ITD序列（是在负链上的实际reads序列和参考基因组反向互补）和ref比对的结果
while(<FI>) {
    if( !/^@/ ) {
	@line_reads = split( /\t/, $_ );
	if( ($line_reads[1] & 16) == 16 ) { next; }

	my %acornerroacornrlocs = (); my $acorneacornrpos = $line_reads[3]-1;
	my $mdstring = substr($line_reads[12],5); #strip off MD:Z:
	while( $mdstring =~ /^(\d+)([ACGT]|\^[ACGT]+)(.*)/ ) {
	    $mdstring = $3;
	    if( substr( $2, 0, 1 ) eq "^" ) {
		$acorneacornrpos += $1 + length($2) - 1;
	    } else {
		$acorneacornrpos += $1 + 1; $acornerroacornrlocs{$acorneacornrpos} = $2;
	    }
	}

	@sub_reads = split( /\|start/, $line_reads[0] );
	my @acornkeys = sort {$a<=>$b} keys %{$acornHashinfos{$sub_reads[0]}};

	my $acornspos = $acornHashinfos{$sub_reads[0]}{$acornkeys[$sub_reads[1]]}{s1}+1; my $acornrpos = $line_reads[3]; my $acornhpos = 0; my $i=0;
	$cigar = $line_reads[5];

	while( $cigar =~ /(\d+)([SHDIMN])(.*)/ ) {
	    my $mmkey = "";
	    $cigar = $3;
	    if( $2 eq "M" ) {
		foreach $acorneacornrpos ( sort {$a<=>$b} keys %acornerroacornrlocs ) {
		    if( $acorneacornrpos >= $acornrpos && $acorneacornrpos < $acornrpos+$1 ) {
			if( length($mmkey) > 0 ) { $mmkey .= "+"; }
			$mmkey .= $acorneacornrpos . ":" . $acornerroacornrlocs{$acorneacornrpos} . ">" .
			    substr( $line_reads[9], $acorneacornrpos+$acornspos-$acornHashinfos{$sub_reads[0]}{$acornkeys[$sub_reads[1]]}{s1}-$acornrpos-$acornhpos-1, 1);
		    }
		}
		$acornHashinfos{$sub_reads[0]}{$acornspos} = {r0=>$acornrpos+0, r1=>$acornrpos+$1-1, s0=>$acornspos+0, s1=>$acornspos+$1-1, mm=>$mmkey };
		$acornspos += $1; $acornrpos += $1;
	    } elsif( $2 eq "I" || $2 eq "S" || $2 eq "H" ) {
		$acornspos += $1;
		if( $2 eq "H" ) { $acornhpos += $1; }
	    } elsif( $2 eq "D" ) {
		$acornrpos += $1;
	    }
	    $i++;
	}
    }
}
### 到这次更新了ITD起始位置到哈希 %acornHashinfoss 就是现在的%acornHashinfos，%acornHashinfos三个键构建完毕，下面就是利用函数proc_alignment_dets_hash加信息
if( !$tmp && (
    system( "rm " . $outbase . "_candidate_clusters_middle.fa" ) ||
    system( "rm " . $outbase . "_candidate_clusters_middle_" . $mappingprefix . ".sam" ) ) ) {
  die "Error，Fail to removing candidate_clusters_middle (insert) files.";
}
###############################这一步是注释的核心，我们进行弱化为基本的信息整理#####################################
my %candidatedets = (); my %acornitd = ();

foreach( keys %acornHashinfos ) {
    my %acorncregion0 = add_dacornHashinfos_info( \%{$acornHashinfos{$_}}, $usearchHash{$_} );
    $candidatedets{$_} = \%acorncregion0;
    $acornitd{$_} = $candidatedets{$_}{posaltmis};
    # $rdotkeys{$_} = $candidatedets{$_}{keyrdot};
    #$cdotkeys{$_} = $candidatedets{$_}{keycdot};
    #$caltkeys{$_} = $candidatedets{$_}{keycdotalt};
}
if( scalar( keys %acornitd ) == 0 ) {
    print "NO ITD\n";
    system( "touch " . $outbase . "_FLT3-ITD_negative.vcf" );
    if( !$tmp && ( system( "rm " . $outbase . ".nowild.R*.fastq" ) ||
        ( $inputbam ne "" && system( "rm " . $outbase . ".R*.fastq" ) ) ) ) {
      die "Error,Fail to removing fastq.";
    }
    exit(0);
}

# my %nonacornitd = ();
# foreach( keys %acornitd ) {
#     if( $acornitd{$_} =~ /dup/ || $acornitd{$_} =~ /itd/ ) { next; }
#     $nonacornitd{$_} = $acornitd{$_};
# }

#### FILTER OUT RELATED ITDS (E.G. PREDICTED SEQUENCING ERRORS)
# Process each ITD size separately
# Then sort by underlying structure with greatest number of extended reads
# Then sort by representative candidate ITD with greatest number of extended reads and fewest mismatches

my %clustITD = ();  # sort into nearest in-frame sizes and remove predicted sequencing errors
my %clustITDTot = ();
my %keynamesMax = ();
foreach( keys %candidatedets ) {
    my $insize = $candidatedets{$_}{totlen} - ($candidatedets{$_}{totlen} % 3);
    if( ($candidatedets{$_}{totlen} % 3) == 2 ) { $insize += 3; }
    @line_reads = split( /\_/, $candidatedets{$_}{posaltmis} );
    my $acornsvkey = "";
    for( my $i=0; $i<scalar(@line_reads); $i++ ) {
	@sub_reads = split(/\(/, $line_reads[$i] );
	if( length($acornsvkey)>0 ) { $acornsvkey .= "_"; }
	$acornsvkey .= $sub_reads[0];
    }
    @line_reads = split( /size=/, $_ );
    $clustITD{$insize}{$acornsvkey}{$_} = ($line_reads[1]+0) + 1/($candidatedets{$_}{mismatches}+1);
    $clustITDTot{$insize}{$acornsvkey} += $line_reads[1];

    if( !exists( $keynamesMax{$insize} ) || $line_reads[1] > $keynamesMax{$insize}{max} ||
	( $line_reads[1] == $keynamesMax{$insize}{max} &&
	  $candidatedets{$_}{mismatches} < $keynamesMax{$insize}{mismatches} ) ) {
	$keynamesMax{$insize}{cname} = $_;
	$keynamesMax{$insize}{max} = $line_reads[1];
	$keynamesMax{$insize}{mismatches} = $candidatedets{$_}{mismatches};
    }
}


my %acornitdfilteredfenzu = (); my $acorn_oneitd_perlength = 0;my $acorn_singleread_itds =0;my $acorn_anylength_itds =0;
# $acorn_anylength_itds 的值改为1就没有ITD插入片段为3的倍数才输出；
my $similarity = 0;
foreach my $size ( keys %clustITDTot ) {
    if( $acorn_oneitd_perlength ) {
	$acornitdfilteredfenzu{ $keynamesMax{$size}{cname} } = 1;
    } elsif( scalar( keys %{$clustITD{$size}} ) > 1 ) {
	my $usearchclust = ""; my $i=0; my $j=0; my @totalitdinfoing = ();
	foreach my $sv ( sort{$clustITDTot{$size}{$b} <=> $clustITDTot{$size}{$a}} keys %{$clustITDTot{$size}} ) {
	    my @sorted = sort{$clustITD{$size}{$sv}{$b} <=> $clustITD{$size}{$sv}{$a}} keys %{$clustITD{$size}{$sv}};
	    $seq = "";
	    @line_reads = split(/\_/, $sv );
	    for( my $k=0; $k<scalar(@line_reads); $k++ ) {
		if( substr($line_reads[$k], 0, 3) eq "ins" ) {
		    $seq .= substr($line_reads[$k], 3);
		} elsif( substr($line_reads[$k], 0, 3) eq "ref" ) {
		    @sub_reads = split( /\-/, substr($line_reads[$k], 3) );
		    $seq .= substr($refseq, $sub_reads[0]-1, $sub_reads[1]-$sub_reads[0]+1);
		}
	    }
	    for( my $k=0; $k<sum( values %{$clustITD{$size}{$sv}} ); $k++ ) {
			#原来代码，修改为删掉长度是3的倍数的判断：if( ( $j==0 || $acornpreHash{$sorted[0]} > 1 || $proc_singleread_itds ) && ( ($acornsizeHash{$sorted[0]} % 3)==0 || $proc_anylength_itds ) )
		if( ( $j==0 || $acornpreHash{$sorted[0]} > 1 || $acorn_singleread_itds ) && ( ($acornsizeHash{$sorted[0]} % 3)==0 || $acorn_anylength_itds ) ) {
		    $i++;
			#print "打印数组@sorted\n";
			#print Dumper(\@sorted);
			#print "打印数组@sorted\n";
		    $usearchclust .= ">" . $sorted[0] . ":" . $i . "\n" . $seq . "\n";
		}
		$j++;
	    }
	}
	if( length($usearchclust)>0 ) {
	    open (FO, ">" . $outbase . "_for_usearchclust_len" . $size  . ".fasta" ) or die $!;
	    print FO $usearchclust;
	    close FO;
	    system( $clustercmd . " -cluster_fast " . $outbase . "_for_usearchclust_len" . $size  . ".fasta -id 0.995 -sort size -centroids " .
			$outbase . "_clusters_len.fasta" .  " -uc " . $outbase . "_clusters_len.uc" );

	    open (FI, $outbase . "_clusters_len.fasta" ) or die $!;
	    my @lines_keyss = ();
		while(<FI>) {
		chomp;
		if($_ =~ /^>/){
			$lines_keys = substr $_, 1;
			@lines_keyss = split( /\:/, $lines_keys );
			$acornitdfilteredfenzu{$lines_keyss[0]} = 1;
			my $itdinfo = acorn_vcf_from_posaltmis($acornitd{$lines_keyss[0]}, '.' );
			my @itdinfos = split( /\t/, $itdinfo );
			my @itdinfoings = [$itdinfos[1],$itdinfos[4],$lines_keyss[0]];
			push(@totalitdinfoing,@itdinfoings);  # 这是多维数组的构建
			# my @itdinfoingstest = [28608240,'AACTCCCATTTGAGATCATATTCATATTCTCCCCTGCCTGAGGG','726_length48_size=106'];
			# push(@totalitdinfoing,@itdinfoingstest);
		}
	    }
	    close FI;

	    if( !$tmp && (
		system( "rm " . $outbase . "_clusters_len.fasta" ) ||
		system( "rm " . $outbase . "_clusters_len.uc" ) ||
		system( "rm " . $outbase . "_for_usearchclust_len" . $size . ".fasta" ) ) ) {
		die "Error,Fail to removing usearch files.";
	    }
	}
		# 加的商检的分组模块
		#处理分组同一片段长度，有多个序列的先用聚类软件得到输出整理为：$acornitdfilteredfenzu{$lines_keyss[0]} = 1; 和@totalitdinfoing 进行二次分组处理
		# print "HHHHHHHHJJJJJJJJJXXXXXXXXXXXXXXXXXXXXXXJJJJJJJXXXXXXXXXXJJJJJJJXXXXXXXXXX\n";
		# print Dumper(\@totalitdinfoing);
		# print "\n";
		# my $kkk = scalar(@totalitdinfoing);
		# print "$kkk\n";
		# print "HHHHHHHHJJJJJJJJJXXXXXXXXXXXXXXXXXXXXXXJJJJJJJXXXXXXXXXXJJJJJJJXXXXXXXXXX\n";
		my $maxseq = ''; my $minseq = '';
		if( scalar(@totalitdinfoing) > 1) {
			foreach my $diedai (combine(2, @totalitdinfoing)){
				print "jjjjjjjjjj\n";
				print "$$diedai[0][0]\n";
				print "$$diedai[0][1]\n";
				print "$$diedai[1][0]\n";
				print "$$diedai[1][1]\n";
				print "jjjjjjjjjj\n";
				my $loc0 = $$diedai[0][0];
				my $item0 = $$diedai[0][1];
				my $loc1 = $$diedai[1][0];
				my $item1 = $$diedai[1][1];
				my $keyy0 = $$diedai[0][2];
				my $keyy1 = $$diedai[1][2];
				if(length($item0) > length($item1)){
					$maxseq = $item0;
					$minseq = $item1;
				} else{
					$maxseq = $item1;
					$minseq = $item0;
				}
				# 求二个序列相似值
				if( system( sprintf( "%s %s --A %s --B %s --out %s_similarity.txt", $python, $similarityfile, $minseq, $maxseq, $outbase ) ) ) {
					$similarity = 0;
      				die "Error,Failed to extract similarity value.\n";}
				open( FI, $outbase . "_similarity.txt") or die $!;
				while(<FI>) {
					chomp;
					print "相似值\n";
					$similarity = $_;
					print "$similarity\n";
					print "相似值\n";
					if($similarity >= 0.99 && $loc0 - $loc1 <= length($item0) + 10){
						print "删商检聚类的itd\n";
						my @line_cell0 = split( /\=/,$keyy0); # '1_length54_size=554
						my @line_cell1 = split( /\=/,$keyy1);
						my $keyyy0 = $line_cell0[1];
						my $keyyy1 = $line_cell1[1];
						print "$keyyy1\n";
						if($keyyy0 > $keyyy1){
							delete($acornitdfilteredfenzu{$keyy1});
						} else {delete($acornitdfilteredfenzu{$keyy0});
						}
					}
				}
			}
		}
    } else {
	@line_reads = keys %{$clustITD{$size}};
	my $sv = $line_reads[0];
	my @sorted = sort{$clustITD{$size}{$sv}{$b} <=> $clustITD{$size}{$sv}{$a}} keys %{$clustITD{$size}{$sv}};
	if( ($acornsizeHash{$sorted[0]} % 3)==0 || $acorn_anylength_itds ) { $acornitdfilteredfenzu{ $sorted[0] } = 1; }
    }
}

if( scalar( keys %acornitdfilteredfenzu ) == 0 ) {
    print "NO ITD.\n";
    system( "touch " . $outbase . "_FLT3-ITD_negative.vcf" );
    if( !$tmp && ( system( "rm " . $outbase . ".nowild.R*.fastq" ) ||
        ( $inputbam ne "" && system( "rm " . $outbase . ".R*.fastq" ) ) ) ) {
      die "Error,Fail to removing fastq files.";
    }
}

########################第6步：这里是评价ITD，二测要有10bp比对并且为定量准备哈希#######################################
my $number_itd = 0; my $headersam = "";my $xiugai = "";
my %defineMutPaired = (); my %acornprimarysamMut = (); my %acornrposMut = ();
# my %acornitdfilteredfenzuxiugai = ();
foreach $_ ( sort {$a cmp $b} keys %acornitdfilteredfenzu ) {
	@line_reads = split( /size=/ , $_ );
	$xiugai = $line_reads[1];
	# print "jjjjjjjHHHHHHHHHHH\n";
	$acornitdfilteredfenzu{$_} = $xiugai;
	# print Dumper(\%acornitdfilteredfenzu);
	# print "jjjjjjjHHHHHHHHHHH\n";
}
foreach my $itd ( sort {$acornpreHash{$b} <=> $acornpreHash{$a} or $b cmp $a} keys %acornitdfilteredfenzu) {
	print "\nheduiheduiehduiehduiehdiehdieudieuidiudieudieudieudeiudi\n";
	print Dumper(\%acornitdfilteredfenzu);
	print "\n";
	print "$itd\n";
	print "\nheduiheduiehduiehduiehdiehdieudieuidiudieudieudieudeiudi\n";
  $number_itd++;
  open(FO, ">" . $outbase . "_candidate_ITD_" . $number_itd . ".fa" ) or die $!;
  print FO ">" . $itd . "\n" . $usearchHash{$itd} . "\n";
  close FO;

  # Build bwa candidates index and perform bwa-mem against candidates
  if ( system( sprintf( "%s index -p %s_candidate_ITD_%s %s_candidate_ITD_%s.fa", $bwacmd, $outbase, $number_itd, $outbase, $number_itd ) ) ||
       system( sprintf( "%s mem -M -O 6 -T 20 %s_candidate_ITD_%s %s.nowild.R1.fastq %s.nowild.R2.fastq > %s_to_candidate_ITD_%s_%s.sam",
		 $bwacmd, $outbase, $number_itd, $outbase, $outbase, $outbase, $number_itd, $mappingprefix ) ) ) {
    die "Error,Fail to aligning nowild.fastq.";
  }

  my %defineMut = (); my %defineMutMJ = (); my %strandMut = ();
  open (FI, $outbase . "_to_candidate_ITD_" . $number_itd . "_" . $mappingprefix . ".sam" ) or die $!;
  while(<FI>) {
    if( /^@/ ) {
	if( /^\@SQ/ ) { $headersam .= $_ };
    } else {
	chomp;
	@line_reads = split( /\t/, $_ );
	if( exists( $defineMutPaired{$line_reads[0]} ) || $line_reads[2] eq "*" ||
	    ( $line_reads[1] & 256 ) == 256 || ( $line_reads[1] & 2048 ) == 2048 ) {
	  next;
	}
	if( ( $line_reads[1] & 64 ) == 64 ) {
	    $read1orread2 = "_1";
	} elsif( ( $line_reads[1] & 128 ) == 128 ) {
	    $read1orread2 = "_2";
	} else {
	    next;
	}

	$readname = $line_reads[0] . $read1orread2;
        $cigar = $line_reads[5];
	my $acornclipLength = 0;
	my $acorninsLength = 0;
	my $acorndelLength = 0;
	if( @a = $cigar =~ /(\d+)S/g ) { $acornclipLength = max(@a); }
	if( @a = $cigar =~ /(\d+)I/g ) { foreach my $x (@a) { $acorninsLength += $x; } }
	if( @a = $cigar =~ /(\d+)D/g ) { foreach my $x (@a) { $acorndelLength += $x; } }

	if( $acorninsLength <= 2 && $acorndelLength <= 2 && $acornclipLength <= 3 &&
		 ($_ =~ /NM:i:(\d+)/ && $1 <= 5) ) {
	  my $acornrpos = $line_reads[3]-1;
	  while( $cigar =~ /(\d+)([SHDIMN])(.*)/ ) {
	    $cigar = $3;
	    if( $2 eq "M" || $2 eq "N" || $2 eq "D" ) { $acornrpos += $1; }
	  }
	  $defineMut{$readname} = 1;  # reads aligning "stringently" to itd genome
	  $acornprimarysamMut{$readname} = $_;
	  $acornrposMut{$readname} = $acornrpos;
	  if( ( $line_reads[1] & 16 ) == 16 ) {
	    $strandMut{$readname} = "-";
	  } else {
	    $strandMut{$readname} = "+";
	  }
	  if( ( $candidatedets{$itd}{mutbk_s0} >= $line_reads[3]+$expand &&
		$candidatedets{$itd}{mutbk_s0} < $acornrpos-$expand ) ||
    	      ( $candidatedets{$itd}{mutbk_s1} > $line_reads[3]+$expand &&
		$candidatedets{$itd}{mutbk_s1} <= $acornrpos-$expand ) ) {
	    $defineMutMJ{$readname} = 1;
	  }
        }
    }
  }
  close FI;

  foreach my $read (keys %defineMut) {
    $read =~ s/\_1$//g; $read =~ s/\_2$//g;
    if( exists($defineMut{$read."_1"}) && exists($defineMut{$read."_2"}) &&
	$strandMut{$read."_1"} ne $strandMut{$read."_2"} &&
        ($defineMutMJ{$read."_1"} || $defineMutMJ{$read."_2"} ) ) {
	$defineMutPaired{$read} = 1;
	$acornreadsPassed{$read} = 1;
    }
  }

  if( !$tmp && (
    system( "rm " . $outbase . "_candidate_ITD_" . $number_itd. "*" ) ||
    system( "rm " . $outbase . "_to_candidate_ITD_" . $number_itd . "_" . $mappingprefix . ".sam" ) ) ) {
    die "Error,Fail to deleting ITD fasta, index, and sam.";
  }
}

if( !$tmp && ( system( "rm " . $outbase . ".nowild.R*.fastq" ) ||
  ( $inputbam ne "" && system( "rm " . $outbase . ".R*.fastq" ) ) ) ) {
  die "Error,Fail to deleting fastq.";
}

open (FO, ">" . $outbase . "_to_candidate_ITDs.filtered.sam" ) or die $!;
print FO $headersam;
foreach( sort keys %defineMutPaired ) {
    print FO $acornprimarysamMut{$_."_1"};
    print FO "\n";
    print FO $acornprimarysamMut{$_."_2"};
    print FO "\n";
}
close FO;

##################################最后一步定量:后期可能会优化，所以变量名变更的少###################################
my %readsMut = (); my %readsMutLeft = (); my %readsMutRight = ();
my %readsMutWt = (); my %readsMutWtLeft = (); my %readsMutWtRight = ();
# my %readsMutMut = (); my %readsMutMutLeft = (); my %readsMutMutRight = ();
my %readsWt = (); my %readsWtLeft = (); my %readsWtRight = ();
my %depthMut = (); my %depthMutLeft = (); my %depthMutRight = ();
my %depthMutWt = (); my %depthMutWtLeft = (); my %depthMutWtRight = ();
my %depthWt = ();  my %depthWtLeft = ();  my %depthWtRight = ();
my %readsMutByR2 = ();
# my $acornlposAdj; my $acornrposAdj;

open( FI, $outbase . "_to_candidate_ITDs.filtered.sam" ) or die $!;
while(<FI>) {
    if( /^@/ ) { next; }
    @line_reads = split( /\t/, $_ );

    $readname = $line_reads[0];
    my $cname = $line_reads[2];
    my $acornrpos=0;
    if( ( $line_reads[1] & 64 ) == 64 ) {
      $acornrpos = $acornrposMut{$readname."_1"};
    } elsif( ( $line_reads[1] & 128 ) == 128 ) {
      $acornrpos = $acornrposMut{$readname."_2"};
    } else {
      next;
    }

    # ADD LOGIC FOR NOBIND CASES (ALLOW RELAXED BUFFER?)

    if( $candidatedets{$cname}{mutbk_s0} >= $line_reads[3]+$expand &&
	$candidatedets{$cname}{mutbk_s0} < $acornrpos-$expand ) {
      if( ( $line_reads[1] & 128 ) == 128 ) { $readsMutByR2{$cname}{$readname} = 1; }
      if( !exists($readsMutLeft{$cname}{$readname}) ) {
	$readsMutLeft{$cname}{$readname} = 1;  #有用
	$readsMut{$cname}{$readname} = 1;
      }
    }
    if( $candidatedets{$cname}{mutbk_s1} > $line_reads[3]+$expand &&
	$candidatedets{$cname}{mutbk_s1} <= $acornrpos-$expand ) {
      if( ( $line_reads[1] & 128 ) == 128 ) { $readsMutByR2{$cname}{$readname} = 1; }
      if( !exists($readsMutRight{$cname}{$readname}) ) {
	$readsMutRight{$cname}{$readname} = 1; #有用
	$readsMut{$cname}{$readname} = 1;
      }
    }
    if( !exists($readsMutWtLeft{$cname}{$readname}) &&
	$candidatedets{$cname}{dup_r0} > $line_reads[3]+$expand &&
	$candidatedets{$cname}{dup_r0} <= $acornrpos-$expand ) {
	$readsMutWtLeft{$cname}{$readname} = 1; # 有用
	$readsMutWt{$cname}{$readname} = 1;
    }
    if( !exists($readsMutWtRight{$cname}{$readname}) &&
	$candidatedets{$cname}{dup_r1}+$candidatedets{$cname}{acornitdsize} >= $line_reads[3]+$expand &&
	$candidatedets{$cname}{dup_r1}+$candidatedets{$cname}{acornitdsize} < $acornrpos-$expand ) {
	$readsMutWtRight{$cname}{$readname} = 1; # 有用
	$readsMutWt{$cname}{$readname} = 1;
    }

    # if( $acornrpos <= $candidatedets{$cname}{dup_r1} ) {
    #   $acornrposAdj = $acornrpos;
    #   $acornlposAdj = $line_reads[3];
    # } elsif( $line_reads[3] >= $candidatedets{$cname}{dup_r1} ) {
    #   $acornrposAdj = $acornrpos - $candidatedets{$cname}{acornitdsize};
    #   $acornlposAdj = $line_reads[3] - $candidatedets{$cname}{acornitdsize};
    # } else {
    #   $acornrposAdj = $candidatedets{$cname}{dup_r1};
    #   $acornlposAdj = $candidatedets{$cname}{dup_r0};
    # }

    # Determine if mutant read also covers other mutant boundaries for depth puacornrposes
    # if( exists( $readsMut{$cname}{$readname} ) ) {
    #   foreach my $itd ( keys %acornitdfilteredfenzu ) {
	# if( $itd ne $cname ) {
    #       if( !exists($readsMutMutLeft{$cname}{$itd}{$readname}) &&
	#     $candidatedets{$itd}{dup_r0} > $acornlposAdj+$expand &&
	#     $candidatedets{$itd}{dup_r0} <= $acornrposAdj-$expand ) {
	#     $readsMutMutLeft{$cname}{$itd}{$readname} = 1;
	#     $readsMutMut{$cname}{$itd}{$readname} = 1;
	#   }
	#   if( !exists($readsMutMutRight{$cname}{$itd}{$readname}) &&
	#     $candidatedets{$cname}{dup_r1} >= $acornlposAdj+$expand &&
	#     $candidatedets{$cname}{dup_r1} < $acornrposAdj-$expand ) {
	#     $readsMutMutRight{$cname}{$itd}{$readname} = 1;
	#     $readsMutMut{$cname}{$itd}{$readname} = 1;
	#   }
	# }
    #   }
    # }
}
close FI;

# my @lens = sort {$readlengths{$b} <=> $readlengths{$a}} keys %readlengths;
# my $readlength = $lens[0];

foreach $readname ( sort keys %acorndefineWTPaired ) {
  foreach $read1orread2 ( ( "_1", "_2" ) ) {
    @line_reads = split( /\t/, $acornprimarysam{$readname . $read1orread2} );
    $cigar = $line_reads[5];
    my $acornrpos = $line_reads[3] - 1;
    while( $cigar =~ /(\d+)([SHDIMN])(.*)/ ) {
      $cigar = $3;
      if( $2 eq "M" || $2 eq "N" || $2 eq "D" ) { $acornrpos += $1; }
    }

    foreach my $cname ( keys %acornitdfilteredfenzu ) {
      if( !exists( $readsWtLeft{$cname}{$readname} ) &&
	  $candidatedets{$cname}{dup_r0} > $line_reads[3]+$expand &&
	  $candidatedets{$cname}{dup_r0} <= $acornrpos - $expand ) {
	$readsWtLeft{$cname}{$readname} = 1; # 有用
	$readsWt{$cname}{$readname} = 1;
      }
      if( !exists( $readsWtRight{$cname}{$readname} ) &&
	  $candidatedets{$cname}{dup_r1} >= $line_reads[3]+$expand &&
	  $candidatedets{$cname}{dup_r1} < $acornrpos - $expand ) {
	$readsWtRight{$cname}{$readname} = 1; # 有用
	$readsWt{$cname}{$readname} = 1;
      }
    }
  }
}

foreach( keys %acornitdfilteredfenzu ) {
  $depthMut{$_} = 0; $depthMutLeft{$_} = 0; $depthMutRight{$_} = 0;
  $depthMutWt{$_} = 0; $depthMutWtLeft{$_} = 0; $depthMutWtRight{$_} = 0;
  $depthWt{$_} = 0; $depthWtLeft{$_} = 0; $depthWtRight{$_} = 0;

  if(exists($readsMut{$_})) {
	  $depthMut{$_} = scalar( keys %{$readsMut{$_}} );
  }
  if(exists($readsMutLeft{$_})) {
	  $depthMutLeft{$_} = scalar( keys %{$readsMutLeft{$_}} );
  }
  if(exists($readsMutRight{$_})) {
	  $depthMutRight{$_} = scalar( keys %{$readsMutRight{$_}} );
  }
  if(exists($readsMutWt{$_})) {
	  $depthMutWt{$_} = scalar( keys %{$readsMutWt{$_}} );
  }
  if(exists($readsMutWtLeft{$_})) {
	  $depthMutWtLeft{$_} = scalar( keys %{$readsMutWtLeft{$_}} );
  }
  if(exists($readsMutWtRight{$_})) {
	  $depthMutWtRight{$_} = scalar( keys %{$readsMutWtRight{$_}} );
  }
  if(exists($readsWt{$_})) {
	  $depthWt{$_} = scalar( keys %{$readsWt{$_}} );
  }
  if(exists($readsWtLeft{$_})) {
	  $depthWtLeft{$_} = scalar( keys %{$readsWtLeft{$_}} );
  }
  if(exists($readsWtRight{$_})) {
	  $depthWtRight{$_} = scalar( keys %{$readsWtRight{$_}} );
  }
}

my %depthRadius = ();
foreach my $cname (keys %acornitdfilteredfenzu) {
    @line_reads = split( /\_/, $cname );
    $depthRadius{$cname} = [(0) x (substr($line_reads[1], 6)+length($refseq))];
}

open( FI, $outbase . "_to_candidate_ITDs.filtered.sam" ) or die $!;
open( FO, ">" . $outbase . "_to_candidate_ITDs.mutreads.sam" ) or die $!;
while(<FI>) {
    if( /^@/ ) {
	print FO $_;
    } else {
	@line_reads = split( /\t/, $_ );
	if( exists($readsMut{$line_reads[2]}{$line_reads[0]}) ) {
	    print FO $_;
	    for( my $i=$line_reads[3]; $i < $line_reads[3]+length($line_reads[9]); $i++ ) {
		$depthRadius{$line_reads[2]}[$i-1] = max($depthRadius{$line_reads[2]}[$i-1],
							  min($i-$line_reads[3]+1, $line_reads[3]+length($line_reads[9])-$i) );
	    }
	}
    }
}
close FO;
close FI;

if( !$tmp && system( "rm " . $outbase . "_to_candidate_ITDs.filtered.sam" ) ) {
  die "Error,Fail to deleting candidateITDs_filtered sam.";
}

# my $rAFsum = 0; my $uAFsum = 0; my $auAFsum = 0;
my %depthDets = (); my %depthTots = ();
foreach my $cname ( keys %acornitdfilteredfenzu ) {
    $depthDets{$cname} = { "0"=>0, "1-5"=>0, "6-10"=>0, "11-25"=>0, "26-50"=>0, ">50"=>0 };
    for( my $i = $candidatedets{$cname}{dup_r0};
	 $i <= 2 * $candidatedets{$cname}{dup_r1} - $candidatedets{$cname}{dup_r0} + $candidatedets{$cname}{inslen} + 1; $i++ ) {
	if( $depthRadius{$cname}[$i] == 0 ) {
	    $depthDets{$cname}{"0"} += 1;
	} elsif( $depthRadius{$cname}[$i] <= 5 ) {
	    $depthDets{$cname}{"1-5"} += 1;
	} elsif( $depthRadius{$cname}[$i] <= 10 ) {
	    $depthDets{$cname}{"6-10"} += 1;
	} elsif( $depthRadius{$cname}[$i] <= 25 ) {
	    $depthDets{$cname}{"11-25"} += 1;
	} elsif( $depthRadius{$cname}[$i] <= 50 ) {
	    $depthDets{$cname}{"26-50"} += 1;
	} else {
	    $depthDets{$cname}{">50"} += 1;
	}
    }
    my $mutCount = 0; my $mwCount = 0;
    # my $mutCountAdj = 0; my $mutCountUmiAdj = 0; my $mwCountAdj = 0; my $mwCountUmiAdj = 0;
	$mutCount = ($depthMutLeft{$cname}+$depthMutRight{$cname})/2;
	$mwCount = ($depthWtLeft{$cname}+$depthWtRight{$cname}+$depthMutWtLeft{$cname}+$depthMutWtRight{$cname})/2;

    $mwCount = max( $mutCount, $mwCount );
    $depthTots{$cname}{Mut} = $mutCount;
    $depthTots{$cname}{All}= $mwCount;
    if( $mwCount == 0 ) { $depthTots{$cname}{AF} = 0; } else { $depthTots{$cname}{AF} = $mutCount/$mwCount; }
}

# Remove itds with no mutant reads and flag ITDs where neither endpoints is within an exon + expand
my %itdother = ();
# 原来是：if( $depthMut{$_} == 0 || $depthTots{$_}{Mut} == 0 ) 之前定的阈值是0，如果支持突变的reads为0就把这个位点进行删除，不要这个ITD
foreach( keys %acornitdfilteredfenzu ) {
    if( $depthMut{$_} < 1 || $depthTots{$_}{Mut} <= 4 ) {
	delete $acornitdfilteredfenzu{$_};
    } elsif( !($candidatedets{$_}{dup_r0} > 585 && $candidatedets{$_}{dup_r0} < 721 ) &&
	     !($candidatedets{$_}{dup_r1} > 585 && $candidatedets{$_}{dup_r1} < 721 ) &&
	     !($candidatedets{$_}{dup_r0} > 808 && $candidatedets{$_}{dup_r0} < 916 ) &&
	     !($candidatedets{$_}{dup_r1} > 808 && $candidatedets{$_}{dup_r1} < 916 ) &&
	     !($candidatedets{$_}{dup_r0} > 393 && $candidatedets{$_}{dup_r0} < 501 ) &&
	     !($candidatedets{$_}{dup_r1} > 393 && $candidatedets{$_}{dup_r1} < 501 )
		) {
		print "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM\n";
	print "$candidatedets{$_}{dup_r0}\n";
	print "$candidatedets{$_}{dup_r1}\n";
		print "$depthTots{$_}{Mut}\n";
		print "$_\n";

	$itdother{$_} = 1;
    }
}

if( glob($outbase . "*.vcf") ) { system( "rm " . $outbase . "*.vcf" ); }
my $outvcf = ""; my $outother = "";
foreach( keys %acornitdfilteredfenzu ) {
    # my $nobind = "None"; my $dupbait = "None"; my $unambig = "None"; my @baits = ();
	my $ID = '.';
    my $vcfstring = acorn_vcf_from_posaltmis($acornitd{$_}, $ID) . "\t" .
	"SVLEN=" . $candidatedets{$_}{acornitdsize} . ";AF=" . sprintf( "%.4g", $depthTots{$_}{AF} ) .
	";SVTYPE=FLT3-ITD" . "\t" . "GT:AD" . "\t0\/1:";
	$vcfstring .= sprintf( "%.4g", $depthTots{$_}{All}-$depthTots{$_}{Mut}) . "," . sprintf( "%.4g", $depthTots{$_}{Mut} ) . "\n"; #alt深度
	# $vcfstring .= "\n";

    my $covLength = sum values %{$depthDets{$_}};
	if( !exists( $itdother{$_} ) && ( $candidatedets{$_}{acornitdsize} < 300 && $depthDets{$_}{"0"}/$covLength < 0.5 )
		&& ($depthTots{$_}{AF} >= 0.005 && $depthTots{$_}{Mut} >= 10) ) {
		$outvcf .= $vcfstring;
		print "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh\n";
		print "$depthTots{$_}{AF}\n";
		print "$depthTots{$_}{Mut}\n";
		print "$_\n";
	} elsif ( exists( $itdother{$_} ) ||
	( $candidatedets{$_}{acornitdsize} >= 300 && $covLength > 0 &&
	  $depthDets{$_}{"0"}/$covLength >= 0.5 ) ||  $depthTots{$_}{AF} <= 0.006 ) {
	$outother .= $vcfstring;
    } else {
		# 如果是其他类型的ITD就走了If语句，就没有走这个else,就没有追加其他类型的ITD到结果文件里
	$outvcf .= $vcfstring;
	# $outsummary .= $tabsummary;
    }
}
my $vcfheader = "##fileformat=VCFv4.2\n" .
	"##source=acornHMAS_ITD v1.0\n" .
	"##FILTER=<ID=PASS,Description=Accept as a confident somatic mutation>\n" .
	"##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n" .
	"##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n" .
	"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=The type of event, FLT3-ITD>\n" .
	"##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n" .
	"##FORMAT=<ID=AD,Number=.,Type=Integer,Description=Allelic depths for the ref and alt alleles in the order listed>\n" .
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";

if( length($outvcf) >= 1 ) {
  open( FO, ">" . $outbase . "_FLT3-ITD.vcf" ) or die $!;
  print FO $vcfheader . $outdir . "\n";
  print FO $outvcf;
  close FO;
} else {
  system( "touch " . $outbase . "_FLT3-ITD_negative.vcf" );
}

if( length($outother) >= 1 ) {
  open( FO, ">" . $outbase . "_filter_ITD.vcf" ) or die $!;
	print FO $vcfheader . $outdir . "\n";
  print FO $outother;
  close FO;
}

if( !$tmp && system( "rm " . $outbase . "_to_candidate_ITDs.mutreads.sam" ) ) {
    die "Error removing to_candidte_ITDs mutreads alignment files. Exiting..." ;
}

#################套用的函数##################

sub add_dacornHashinfos_info {
    # argument #1: %{} 是对$_[0]解引用。$acornHashinfoss里面的去掉一级键：640_length45_size=36，得到的：{$1+1就是ITD开始的位置+1: {{r0=>$rpos+0 $rpos是ITD开始的位置, r1=>$rpos+$1-1 $1是ITD序列长度, s0=>$spos+0 $spos就是ITD开始的位置, s1=>$spos+$1-1, mm=>$mmkey }}}
    my %acornHashinfos = %{$_[0]};

    # argument #2: # 得到的是$sumaHash{$key} = $extSeq 中的值$extSeq，整个ref序列，并且里面有ITD序列
    my $candidateseq = $_[1];
	my @acornkeys = sort {$a<=>$b} keys %acornHashinfos;
    my %acornHashinfosNew = ();
	my %acorndis  = ();
	my $keytmp;
    for( my $i=1; $i<scalar(@acornkeys); $i++ ) {
	if( ($acornHashinfos{$acornkeys[$i]}{s0}-$acornHashinfos{$acornkeys[$i-1]}{s1}) > 0 &&
	    ($acornHashinfos{$acornkeys[$i]}{s0}-$acornHashinfos{$acornkeys[$i-1]}{s1}) <= 4 &&
	    ($acornHashinfos{$acornkeys[$i]}{s0}-$acornHashinfos{$acornkeys[$i-1]}{s1}) == ($acornHashinfos{$acornkeys[$i]}{r0}-$acornHashinfos{$acornkeys[$i-1]}{r1}) ) {
	    $acorndis{$i}=1;
	}
    }
    for( my $i=0; $i<scalar(@acornkeys); $i++ ) {
	if( !exists($acorndis{$i}) ) {
	    # make a new copy to avoid overwriting original hash passed by reference
	    $keytmp = $i;
	    $acornHashinfosNew{$acornkeys[$keytmp]}{s0} = $acornHashinfos{$acornkeys[$keytmp]}{s0};
	    $acornHashinfosNew{$acornkeys[$keytmp]}{s1} = $acornHashinfos{$acornkeys[$keytmp]}{s1};
	    $acornHashinfosNew{$acornkeys[$keytmp]}{r0} = $acornHashinfos{$acornkeys[$keytmp]}{r0};
	    $acornHashinfosNew{$acornkeys[$keytmp]}{r1} = $acornHashinfos{$acornkeys[$keytmp]}{r1};
	    $acornHashinfosNew{$acornkeys[$keytmp]}{mm} = $acornHashinfos{$acornkeys[$keytmp]}{mm};
	} else {
	    $acornHashinfosNew{$acornkeys[$keytmp]}{s1} = $acornHashinfos{$acornkeys[$i]}{s1};
	    $acornHashinfosNew{$acornkeys[$keytmp]}{r1} = $acornHashinfos{$acornkeys[$i]}{r1};
	    my $offset = $acornHashinfos{$acornkeys[$i]}{s0} - $acornHashinfos{$acornkeys[$i-1]}{s1} - 1;
	    if( $offset > 0 ) {
                print "符合mm的条件打印打印打印打印打印打印打印打印打印\n";
		for( my $j=1; $j<=$offset; $j++ ) {
		    if( substr( $candidateseq, $acornHashinfos{$acornkeys[$i-1]}{s1}+$j-1, 1 ) ne
			substr( $refseq, $acornHashinfos{$acornkeys[$i-1]}{r1}+$j-1, 1 ) ) {
			if( length($acornHashinfosNew{$acornkeys[$keytmp]}{mm}) > 0 ) {
			    $acornHashinfosNew{$acornkeys[$keytmp]}{mm} .= "+";
			}
			$acornHashinfosNew{$acornkeys[$keytmp]}{mm} .= ($acornHashinfos{$acornkeys[$i-1]}{r1}+$j) . ":" .
			    substr( $refseq, $acornHashinfos{$acornkeys[$i-1]}{r1}+$j-1, 1 ) . ">" .
			    substr( $candidateseq, $acornHashinfos{$acornkeys[$i-1]}{s1}+$j-1, 1 );
		    }
		}
	    }
	    if( length($acornHashinfos{$acornkeys[$i]}{mm}) > 0 ) {
		if( length($acornHashinfosNew{$acornkeys[$keytmp]}{mm}) > 0 ) {
		    $acornHashinfosNew{$acornkeys[$keytmp]}{mm} .= "+";
		}
		$acornHashinfosNew{$acornkeys[$keytmp]}{mm} .= $acornHashinfos{$acornkeys[$i]}{mm};
	    }
	}
    }
	# 构建哈希：%cregion
    my %cregion = (posaltmis => "", duacornpriseq => "", insseq => "", duplen => 0, inslen => 0, totlen => length($candidateseq), acornitdsize => length($candidateseq)-1500,
	dup_r0 => -1, dup_r1 => -1, mutbk_s0 => -1, mutbk_s1 => -1, mismatches => 0);

    my $acorndiff0 = 0;
    my $acornmmtot = 0;
    my $acornrstartNew = 0;
    my @acornkeysNew = sort {$a<=>$b} keys %acornHashinfosNew;

    # process the full key
    my $altposmis = "ref" . $acornHashinfosNew{$acornkeysNew[0]}{r0} . "-" . $acornHashinfosNew{$acornkeysNew[0]}{r1};
    for( my $i=0; $i<scalar(@acornkeysNew); $i++ ) {
	if( $i>0 ) {
	    $acorndiff0 = $acornHashinfosNew{$acornkeysNew[$i]}{s0} - $acornHashinfosNew{$acornkeysNew[$i-1]}{s1} - 1;
	    $acornrstartNew = $acornHashinfosNew{$acornkeysNew[$i]}{r0} - min($acorndiff0, 0);

	    if( $acorndiff0 > 0 ) {
		$altposmis .= "_ins" . substr( $candidateseq, $acornHashinfosNew{$acornkeysNew[$i-1]}{s1}, $acorndiff0 );
	    }
	    $altposmis .= "_ref" . $acornrstartNew . "-" . $acornHashinfosNew{$acornkeysNew[$i]}{r1};
	}

	if( length($acornHashinfosNew{$acornkeysNew[$i]}{mm}) > 0 ) {
	    $altposmis .= "(" . $acornHashinfosNew{$acornkeysNew[$i]}{mm} . ")";
	    my @acorntmparr = split( /\+/, $acornHashinfosNew{$acornkeysNew[$i]}{mm} );
	    $acornmmtot += scalar( @acorntmparr );
	}
	# Need logic for when $acorndiff0 < 0 and mismatch is in this overlap region
    }
    $cregion{posaltmis} = $altposmis;
    $cregion{mismatches} = $acornmmtot;

    # create rdot, cdot, and pdot keys
    my $acornleftr0 = ""; my $acornleftr1 = ""; my $acornleftmm = "";
    my $acornright0 = ""; my $acornrightr1 = ""; my $acornrightmm = "";
    my @cells; my @subcells; my @subsubcells;

    @cells = split(/\_/, $altposmis);
    @subcells = split(/\(/, $cells[0]);
    @subsubcells = split(/\-/, $subcells[0]);
    if( substr($subsubcells[0],0,3) eq "ref" ) { $acornleftr0 = substr($subsubcells[0],3); $acornleftr1 = $subsubcells[1]; }
    if( scalar(@subcells) > 1 ) { $acornleftmm = substr($subcells[1],0,length($subcells[1])-1); }
    @subcells = split(/\(/, $cells[scalar(@cells)-1]);
    @subsubcells = split(/\-/, $subcells[0]);
    if( substr($subsubcells[0],0,3) eq "ref" ) { $acornright0 = substr($subsubcells[0],3); $acornrightr1 = $subsubcells[1]; }
    if( scalar(@subcells) > 1 ) { $acornrightmm = substr($subcells[1],0,length($subcells[1])-1); }

    #my @cdotL1 = proc_to_cdot_coords( $acornleftr1 );
    #my @cdotR0 = proc_to_cdot_coords( $acornright0 );
    $cregion{dup_r0} = $acornright0;
    $cregion{dup_r1} = $acornleftr1;
    $cregion{mutbk_s0} = $acornleftr1;

    if( $acornright0 > $acornleftr1 ) {
	$cregion{dup_r0} = $acornleftr1;
	$cregion{dup_r1} = $acornright0;
	$cregion{mutbk_s0} = $acornleftr1;
	$cregion{mutbk_s1} = $cregion{totlen} - (1500-$acornright0);

	@subcells = (); my $sub0 = -1; my $sub1 = -1; my $acornflag = 0;

	# should improve this to make more general
	if( scalar(@cells) > 3 && substr($cells[1],0,3) eq "ins" && substr($cells[2],0,3) eq "ref" ) {
	    @subcells = split(/\-/, substr($cells[2],3) );
	} elsif( scalar(@cells) > 2 && substr($cells[1],0,3) eq "ref" ) {
	    @subcells = split(/\-/, substr($cells[1],3) );
	}
	if( scalar(@subcells) > 0 ) {
	    $sub0 = $subcells[0];
	    @subsubcells = split(/\(/, $subcells[1]);
	    $sub1 = $subsubcells[0];

	    if( $sub0 < $acornleftr1 && $acornleftr1-$sub0+1 <= $cregion{totlen}-1500 ) {
		$cregion{dup_r0} = $sub0;
		$cregion{dup_r1} = $acornleftr1;
		$cregion{duplen} = $cregion{dup_r1}-$cregion{dup_r0}+1;
		$cregion{duacornpriseq} = substr( $refseq, $cregion{dup_r0}-1, $cregion{duplen} );
		if( substr($cells[1],0,3) eq "ins" ) {
		    $cregion{insseq} = substr($cells[1],3);
		    $cregion{inslen} = length($cregion{insseq});
		}
		$cregion{mutbk_s1} = $cregion{mutbk_s0} + $cregion{inslen} + 1;
		$acornflag = 1;
	    }
	}

	@subcells = (); $sub0 = -1; $sub1 = -1;
	if( scalar(@cells) > 3 && substr($cells[scalar(@cells)-2],0,3) eq "ins" &&
	    substr($cells[scalar(@cells)-3],0,3) eq "ref" ) {
	    @subcells = split(/\-/, substr($cells[scalar(@cells)-3],3) );
	} elsif( scalar(@cells) > 2 && substr($cells[scalar(@cells)-2],0,3) eq "ref" ) {
	    @subcells = split(/\-/, substr($cells[scalar(@cells)-2],3) );
	}
	if( scalar(@subcells) > 0 && !$acornflag ) {
	    $sub0 = $subcells[0];
	    @subsubcells = split(/\(/, $subcells[1]);
	    $sub1 = $subsubcells[0];
	    if( $sub1 > $acornright0 && $sub1-$acornright0 <= $cregion{totlen}-1500 && substr($cells[scalar(@cells)-2],0,3) eq "ins" ) {
		$cregion{dup_r0} = $acornright0;
		$cregion{dup_r1} = $sub1;
		$cregion{duplen} = $cregion{dup_r1}-$cregion{dup_r0}+1;
		$cregion{duacornpriseq} = substr( $refseq, $cregion{dup_r0}-1, $cregion{duplen} );
		if( substr($cells[scalar(@cells)-2],0,3) eq "ins" ) {
		    $cregion{insseq} = substr($cells[scalar(@cells)-2],3);
		    $cregion{inslen} = length($cregion{insseq});
		}
		$cregion{mutbk_s0} = $cregion{mutbk_s1} - $cregion{inslen} - 1;
	    }
	}
    } elsif( scalar(@cells) == 2 || ( scalar(@cells) == 3 && substr($cells[1],0,3) eq "ins" ) ) {
	my $acornrightmm0 = ""; my $acornrightmm1 = "";
	@subcells = split( /\+/, $acornrightmm );
	for( my $i=0;$i<scalar(@subcells); $i++ ) {
            print "YESYES\n";
	    @subsubcells = split( /:/, $subcells[$i] );
	    if( $subsubcells[0] <= $acornleftr1 ) {
		if( length($acornrightmm0)>0 ) { $acornrightmm0 .= "+"; }
		$acornrightmm0 .= $subcells[$i];
	    } else {
		if( length($acornrightmm1)>0 ) { $acornrightmm1 .= "+"; }
		$acornrightmm1 .= $subcells[$i];
	    }
	}
	my $acornleftmm0 = ""; my $acornleftmm1 = "";
	@subcells = split( /\+/, $acornleftmm );
	for( my $i=0;$i<scalar(@subcells); $i++ ) {
	    @subsubcells = split( /:/, $subcells[$i] );
	    if( $subsubcells[0] >= $acornright0 ) {
		if( length($acornleftmm1)>0 ) { $acornleftmm1 .= "+"; }
		$acornleftmm1 .= $subcells[$i];
	    } else {
		if( length($acornleftmm0)>0 ) { $acornleftmm0 .= "+"; }
		$acornleftmm0 .= $subcells[$i];
	    }
	}

	if( scalar(@cells) == 3 ) {
	    $cregion{insseq} = substr($cells[1],3);
	    $cregion{inslen} = length($cregion{insseq});
	    $cregion{mutbk_s1} = $acornleftr1 + $cregion{inslen} + 1;
	} else {
	    $cregion{mutbk_s1} = $acornleftr1 + 1;
	}
	$cregion{duplen} = $cregion{dup_r1}-$cregion{dup_r0}+1;
	$cregion{duacornpriseq} = substr( $refseq, $cregion{dup_r0}-1, $cregion{duplen} );
    } else {
		my $inacornsecsstr = "";
		my $inslen = 0;
		for (my $i = 1; $i < scalar(@cells) - 1; $i++) {
			if (length($inacornsecsstr) > 0) {$inacornsecsstr .= "/";}
			if (substr($cells[$i], 0, 3) eq "ins") {
				$inacornsecsstr .= substr($cells[$i], 3);
				$inslen += length(substr($cells[$i], 3));
			}
			elsif (substr($cells[$i], 0, 3) eq "ref") {
				@subcells = split(/\(/, substr($cells[$i], 3));
				@subsubcells = split(/\-/, $subcells[0]);
				$inslen += $subsubcells[1] - $subsubcells[0] + 1;
				#my @cdot0 = proc_to_cdot_coords($subsubcells[0]);
				#my @cdot1 = proc_to_cdot_coords($subsubcells[1]);

				$cregion{insseq} = "needs_manual_review";
				$cregion{inslen} = $inslen;
				$cregion{mutbk_s1} = $cregion{mutbk_s0} + $inslen + 1;
			}
		}
	}

    return( %cregion );
}

sub acorn_vcf_from_posaltmis {
    # argument #1: posaltmis
    my $keyf = $_[0];
	my $id = $_[1];

    my @regions = split(/\_/, $keyf); # my @regions = （ref1-676,ref656-1500）
    my $ref0 = $regions[0]; my $ref1 = $regions[scalar(@regions)-1];
    if( scalar(@regions) < 2 || $ref0 !~ /^ref/ || $ref1 !~ /^ref/ ) { return(""); }

    $ref0 =~ s/ref//; $ref1 =~ s/ref//;
    my $mm0 = ""; my $mm1 = "";
    my @cells = split( /\(/, $ref0 );
    my @coords0 = split( /\-/, $cells[0] );
    if( scalar(@cells) > 1 ) { $mm0 = $cells[1]; $mm0 =~ s/\)//; }

    @cells = split( /\(/, $ref1 );
    my @coords1 = split( /\-/, $cells[0] );
    if( scalar(@cells) > 1 ) { $mm1 = $cells[1]; $mm1 =~ s/\)//; }

    my $seq0 = $refseq; my $seq1 = $refseq; my @acorneacornrposs0 = (); my @acorneacornrposs1 = ();
    if( length($mm0) > 0 ) {
	@cells = split(/\+/, $mm0);
	foreach my $mm ( @cells ) {
	    my @subcells = split( /:/, $mm );
	    push @acorneacornrposs0, $subcells[0];
	    my @subsubcells = split (/>/, $subcells[1] );
	    substr( $seq0, $subcells[0]-1, 1 ) = $subsubcells[1];
	}
    }

    if( length($mm1) > 0 ) {
	@cells = split(/\+/, $mm1);
	foreach my $mm ( @cells ) {
	    my @subcells = split( /\:/, $mm );
	    push @acorneacornrposs1, $subcells[0];
	    my @subsubcells = split (/>/, $subcells[1] );
	    substr( $seq1, $subcells[0]-1, 1 ) = $subsubcells[1];
	}
    }

    my $inacornspos = $coords0[1];
    if( length($mm0) > 0 ) { $inacornspos = min(@acorneacornrposs0)-1; }
    my $dupstart = $coords1[0];
    if( length($mm1) > 0 ) { $dupstart = max(@acorneacornrposs1)+1; }

    my $ins = "";
    if( scalar(@regions) > 2 ) {
	for( my $i=1; $i < scalar(@regions)-1; $i++ ) {
	    my $tmacornpriseq = $regions[$i];
	    if( $tmacornpriseq =~ /^ins/ ) {
		$tmacornpriseq =~ s/ins//;
	    } elsif( $tmacornpriseq =~ /^ref/ ) {
		$tmacornpriseq =~ s/ref//;
		@cells = split( /\(/, $tmacornpriseq );
		my @coords = split( /\-/, $cells[0] ); my $seq2 = $refseq;
		if( scalar(@cells) > 1 ) {
		    my $mm2 = $cells[1]; $mm2 =~ s/\)//;
		    @cells = split(/\+/, $mm2);
		    foreach my $mm ( @cells ) {
			my @subcells = split( /\:/, $mm );
			my @subsubcells = split (/>/, $subcells[1] );
			substr( $seq2, $subcells[0]-1, 1 ) = $subsubcells[1];
		    }
		}
		$tmacornpriseq = substr( $seq2, $coords[0]-1, $coords[1]-$coords[0]+1 );
	    } else {
		return( "" );
	    }
	    $ins .= $tmacornpriseq;
	}
    }

    my $alt;
    my $ref;
    if( $inacornspos+1 >= $dupstart ) {
	$ref = substr( $refseq, $inacornspos, 1 );
	$alt = substr($seq0, $inacornspos, $coords0[1]-$inacornspos) . $ins .
	    substr($seq1, $coords1[0]-1, $inacornspos-$coords1[0]+1) . $ref;
    } else {
	$ref = substr( $refseq, $inacornspos, $dupstart-$inacornspos-1 );
	$alt = substr($seq0, $inacornspos, $coords0[1]-$inacornspos) . $ins .
	    substr($seq1, $coords1[0]-1, $dupstart-$coords1[0]);
    }
    $alt = reverse $alt; $alt =~ tr/ATGCatgc/TACGtacg/;
    $ref = reverse $ref; $ref =~ tr/ATGCatgc/TACGtacg/;

    return("chr13\t" . ($refendpos-$inacornspos) . "\t" . $id . "\t" . $ref . "\t" . $alt . "\t.\tPASS" );
}
