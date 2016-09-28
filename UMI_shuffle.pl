#!/usr/bin/perl -w
#
########################################
#
# File Name:
#   UMI_dedup.pl
# 
# Description:
#   deduplicate based on UMI
# 
# Usage:
#   
# 
# Author:
#   Xi Wang, Xi.Wang@mdc-berlin.de
# 
# Date:
#   Fri Sep 23 17:28:55 CEST 2016
#
########################################

use strict;
use lib qw(/home/bq_xwang/tools/perl_package/share/perl/5.14 /home/bq_xwang/tools/perl_package/share/perl/5.14.2 /home/bq_xwang/tools/perl_package/lib/perl/5.14.2/);
use Getopt::Long;
use IO::Zlib;
use List::Util qw(shuffle);
my $debug_lvl = 2;

#####################
my $usage = "$0 <in.bed.gz> <out.bed.gz>\n";
my $infile = shift || die $usage;
my $outfile = shift || die $usage;

my $IN = IO::Zlib->new("$infile", "rb") || die "Cannot open $infile for reading.\n";
my $OUT = IO::Zlib->new("$outfile", "wb9") || die "Cannot open $outfile for writing\n";

my (@UMIs, @UMIs_shuffle);
my @lines = <$IN>; 
for(my $i=0; $i<@lines; $i++) { 
  my @a= split /\t/, $lines[$i]; 
  push @UMIs, $a[3]; 
}
@UMIs_shuffle = shuffle(0..$#UMIs); 
for(my $i=0; $i<@lines; $i++) { 
  my @a= split /\t/, $lines[$i]; 
  $a[3] = $UMIs[$UMIs_shuffle[$i]]; 
  print $OUT join "\t", @a;
}
undef $OUT; 

