#!/usr/bin/perl -w
#
########################################
#
# File Name:
#   UMI_trim.pl
# 
# Description:
#   Trimming steps:
#   1. 5'end 6 bases, if qaulity score < 17, remove; otherwise, trim the 6 bases as UMI and record the UMI
#   2. after trimming, further trim up to 9 G's.
#   3. 3'end: remove A's and qaulity score <=2. if the resulted sequences < 25, remove the read.   
# 
# Usage:
#   
# 
# Author:
#   Xi Wang, Xi.Wang@dkfz.de
# 
# Date:
#   Tue Sep 20 20:15:41 CEST 2016
#
########################################

use strict;
use lib qw(/home/bq_xwang/tools/perl_package/share/perl/5.14 /home/bq_xwang/tools/perl_package/share/perl/5.14.2);
use Getopt::Long;
use IO::Zlib;
my $debug_lvl = 2;

##############
my $usage = "$0 [options] <in.fq.gz> <out.fq.gz.prefix>
Options:
  --phred64-quals  [default: phred32 quals]
  --UMI-min-quals  [default: 17]
  --UMI-length     [default: 6]
  --G-max-length   [default: 9]
  --min-read-len   [default: 25]
  --keep-3p-As     [default: remove]
  --keep-3p-lowQ   [default: remove]
  --lowQ-max       [default: 2]
  --help           [print this help info & exit]
\n";

#### options #####
my $phred64_quals = 0;
my $UMI_min_quals = 17;
my $UMI_length = 6;
my $G_max_length = 9; 
my $min_read_length = 25; 
my $keep_3p_As = 0; 
my $keep_3p_lowQ = 0; 
my $lowQ_max = 2; 
my $help = 0;

GetOptions ('phred64-quals' => \$phred64_quals, 'UMI-min-quals=i' => \$UMI_min_quals, 'UMI-length=i' => \$UMI_length, 'G-max-length=i' => \$G_max_length, 'min-read-len=i' => \$min_read_length, 'keep-3p-As' => \$keep_3p_As, 'keep-3p-lowQ' => \$keep_3p_lowQ, 'lowQ-max=i' => \$lowQ_max, 'help' => \$help) || die $usage;

my $in_file = shift || die $usage;
my $out_pre = shift || die $usage;
$help && die $usage;
my $qual_base = 32; 
$qual_base = 64 if($phred64_quals);

###################

my $IN = IO::Zlib->new($in_file, "rb") || die "Cannot open $in_file for reading.\n";
my $OUT1 = IO::Zlib->new("$out_pre.trimmed.fq.gz", "wb9") || die "Cannot open $out_pre.trimmed.fq.gz for writing\n";
my $OUT2 = IO::Zlib->new("$out_pre.discarded.fq.gz", "wb9") || die "Cannot open $out_pre.discarded.fq.gz for writing\n";

my $n=0;
my @read;
my %proc_result = (0 => 0, 1 => 0, 2 => 0);

while(<$IN>) {
  chomp;
  $read[$n] = $_; 
  $n++;
  if($n == 4) { # process the read
    $n = 0; 
    my $ret = &proc_read(\@read); 
    $proc_result{$ret} ++; 
  }
}

my $tot = $proc_result{0} + $proc_result{1} + $proc_result{2};
print "Processed: $tot\n";
print "Kept: $proc_result{0}\n";
print "UMI low quals: $proc_result{1}\n";
print "Too short: $proc_result{2}\n";

#################
sub proc_read { 
  my @read_ori = @{$_[0]};
  my @read = @read_ori; 
  ### check 5' quals 
  for(my $i=0; $i<$UMI_length; $i++) { 
    my $cur_q = substr($read_ori[3], $i, 1); 
    if((ord($cur_q) - $qual_base) < $UMI_min_quals) { 
      $read_ori[0] = $read_ori[0]." UMI_Low_Quals";
      print $OUT2 join "\n",@read_ori;
      print $OUT2 "\n";
      return(1);
    }
  }
  ### trim UMI
  my @a = split '\s', $read[0];
  $a[0] = $a[0].":UMI=".substr($read_ori[1], 0, $UMI_length);
  $read[0] = join " ", @a;
  $read[1] = substr($read_ori[1], $UMI_length);
  $read[3] = substr($read_ori[3], $UMI_length);

  ### trim G's from 5' 
  my $G_len = 0; 
  for(my $i=0; $i<$G_max_length; $i++) {
    last if("G" ne substr($read[1], $i, 1)); 
    $G_len ++;
  }
  if($G_len > 0) { 
    $read[1] = substr($read[1], $G_len);
    $read[3] = substr($read[3], $G_len);
  }
  ### check 3'end (remove A's, low qual bases)
  my $trim_3_len = 0; 
  my $read_len = length($read[1]);
  for(my $i = -1; $i > -$read_len; $i --) { 
    last if("A" ne substr($read[1], $i, 1) && (ord(substr($read[3], $i, 1))-$qual_base) > $lowQ_max);
    $trim_3_len ++; 
  }
  if($trim_3_len > 0) { 
    $read[1] = substr($read[1], 0, $read_len-$trim_3_len);
    $read[3] = substr($read[3], 0, $read_len-$trim_3_len);
  }
  
  ### final length check
  if(length($read[1]) < $min_read_length) { 
    $read_ori[0] = $read_ori[0]." Too_Short_After_Trimming";
    print $OUT2 join "\n",@read_ori;
    print $OUT2 "\n";
    return(2);
  }
  
  ### output
  print $OUT1 join "\n",@read;
  print $OUT1 "\n";
  return(0)
}

