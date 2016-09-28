#!/usr/bin/perl -w
#
########################################
#
# File Name:
#   UMI_collapse.pl
# 
# Description:
#   
# 
# Usage:
#   
# 
# Author:
#   Xi Wang, Xi.Wang@mdc-berlin.de
# 
# Date:
#   Mon Sep 26 17:38:46 CEST 2016
#
########################################

use strict;
use lib qw(/home/bq_xwang/tools/perl_package/share/perl/5.14 /home/bq_xwang/tools/perl_package/share/perl/5.14.2);
use Getopt::Long;

my $usage = "$0 [options] <in.bed.gz> <out.bed.gz>
Options:
  --min-map-qual  [default: 50]
  --help          [show this help information]
\n"; 

#### options #####
my $min_map_qual = 50;
my $help = 0;
GetOptions ('min-map-qual=i' => \$min_map_qual, 'help' => \$help) || die $usage;
$help && die $usage;

my $infile = shift || die $usage;
my $outfile = shift || die $usage;

system("zcat $infile | awk -vMIN=$min_map_qual '\$5 >= MIN {\$5=0; print}' | sort -k1,1V -k6,6 -k2,3n | uniq -c | awk -vOFS='\t' '{print \$2,\$3,\$4,\$5,\$1,\$7,\$8,\$9,0,\$11,\$12,\$13}' | gzip > $outfile");

