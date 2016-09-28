#!/usr/bin/perl -w
#
########################################
#
# File Name:
#   UMI_div_stat.pl
# 
# Description:
#   UMI diversity at the same genomic postion 
# 
#   modification: 
#     output two more stat: uni-UMI position percentage, dist within 1 postion percentage
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
use List::Util qw(min max sum);
use IO::Zlib;
my $debug_lvl = 2;

#####################
my $usage = "$0 [options] <in.bed.gz>
Options: 
  --min-UMI-cnt     [default: 1]
  --help            [show this help & exit]
\n";

#### options #####
my $min_UMI_cnt = 1;
my $help = 0;
GetOptions ('min-UMI-cnt=i' => \$min_UMI_cnt, 'help' => \$help) || die $usage;
$help && die $usage;
my $infile = shift || die $usage;

### 
my $IN = IO::Zlib->new("$infile", "rb") || die "Cannot open $infile for reading.\n";

## read $IN & cluster alignments
my $cur_chr = "chr1";
my $cur_start = 0;
my $cur_end = 0;
my $cur_strand = "-";
my $cur_s = 0;
my $cur_l = 0;
my %clusters;
my @cur_cluster; 
my $key;
while(<$IN>) {
  chomp;
  my @a = split; 
  next if($a[4] < $min_UMI_cnt);
  if($cur_chr eq $a[0] && $cur_strand eq $a[5] && $cur_start ==$a [1] && $cur_end == $a[2] && $cur_s eq $a[11] && $cur_l eq $a[10]) { # same with cur region
    push @cur_cluster, $_; 
  } else { # 
    ### put the cluster into hash
    $key = join "_", ($cur_chr, $cur_start, $cur_end, $cur_strand, $cur_l, $cur_s);
    $clusters{$key} = [ @cur_cluster ] unless $key eq "chr1_0_0_-_0_0"; ## 
    ### renew array
    @cur_cluster = (); 
    ### renew cur region
    $cur_chr = $a[0];
    $cur_strand = $a[5];
    $cur_start = $a[1];
    $cur_end = $a[2];
    $cur_l = $a[10];
    $cur_s = $a[11];
    push @cur_cluster, $_; 
  }
}
$key = join "_", ($cur_chr, $cur_start, $cur_end, $cur_strand, $cur_l, $cur_s);
$clusters{$key} = [ @cur_cluster ]; 

my %UMI_dist_cnt; 
my ($total_pos, $uni_UMI_pos, $dist_1_pos);
## in each cluster, select UMI and output 
foreach $key (sort keys %clusters) { 
  $total_pos ++; 
  #print $key."\n"; 
  #print join "\n", @{$clusters{$key}};
  #print "\n\n";
  my @UMI_dist = @{&same_pos_UMI_dist($clusters{$key})};
  my $max_dist = max(@UMI_dist); 
  $uni_UMI_pos ++ if(($max_dist) == 0);
  $dist_1_pos ++ if(($max_dist) == 1);
  for(my $i=0; $i<@UMI_dist; $i++) { 
    $UMI_dist_cnt{$UMI_dist[$i]} ++; 
  }
}

foreach my $d (sort keys %UMI_dist_cnt) { 
  print $d, "\t", $UMI_dist_cnt{$d}, "\n"; 
}

print "\n";
print "total\t$total_pos\n";
print "uni_UMI\t$uni_UMI_pos\n";
print "dist_1\t$dist_1_pos\n";

###########################
sub same_pos_UMI_dist { 
  my @lines = @{$_[0]};
  my %cluster_hash; 
  my $junc_pos;
  my @ret; 
  foreach my $l (@lines) { 
    my @a = split "\t", $l;
    if($a[9] == 1) {
      $junc_pos = 0;
    } else { 
      my @l = split /,/,$a[10]; 
      my @s = split /,/,$a[11]; 
      my @p; 
      for(my $i=1;$i<$a[9];$i++) { 
        push @p, $a[1] + $l[$i-1] + $s[$i-1];
        push @p, $a[1] + $s[$i]; 
      }
      $junc_pos = join ',',@p; 
    }
    my @UMI = ($a[1],$a[2],$a[3],$junc_pos); 
    my $k = join "\t", @UMI; 
    $cluster_hash{$k} = $a[4];
  }
 
  #foreach my $k (sort {$cluster_hash{$b} <=> $cluster_hash{$a}} keys %cluster_hash) { 
  #  print $k, "\t", $cluster_hash{$k}, "\n";
  #}
  my @UMI = sort {$cluster_hash{$b} <=> $cluster_hash{$a}} keys %cluster_hash; 
  if($#UMI == 0) { 
    push @ret, 0; 
  } else { 
    my @aa = split /\t/, $UMI[0];
    my @next = @{&next_kmers($aa[2])};
    my %next_hash = map { $_ => 1 } @next;
    for(my $i=1; $i<@UMI; $i++) { 
      my @bb = split /\t/, $UMI[$i];
      if(exists($next_hash{$bb[2]})) { 
        push @ret, 1; 
        next; 
      } 
      push @ret, &mismatch($aa[2], $bb[2]); 
    }
  }
  return(\@ret); 
}

sub next_kmers { 
  my $kmer = $_[0]; 
  #print $kmer, "\n";
  my %ret; 
  my @ab = ("A", "C", "G", "T"); 
  my @a = split //,$kmer; 
  for(my $i=0; $i<@a; $i++) { 
    for(my $j=0; $j<@ab; $j++) { 
      if($a[$i] ne $ab[$j]) { 
        my $temp = join '', (@a[0..($i-1)], $ab[$j], @a[($i+1)..$#a]);
        $ret{$temp} ++;
      }
    }
  }
  for(my $j=0; $j<@ab; $j++) { 
    my $temp = join '', (@a[1..$#a], $ab[$j]);
    $ret{$temp} ++;
    $temp = join '', ($ab[$j], @a[0..($#a-1)]);
    $ret{$temp} ++;
  }
  delete($ret{$kmer}) if(exists($ret{$kmer}));
  my @r = sort keys(%ret); 
  return(\@r); 
}

sub mismatch { 
  my $a = $_[0]; 
  my $b = $_[1]; 
  my $ret = 0; 
  my @aa = split //, $a; 
  my @ba = split //, $b; 
  for(my $i=0; $i<@aa; $i++) { 
    $ret ++ if($aa[$i] ne $ba[$i]); 
  }
  print "$a\t$b\n@aa\t@ba\n" if($ret > 6);
  return($ret); 
}

sub all_kmer { 
  my $k = $_[0]; 
  my @ret; 
  my @ab = qw(A C G T); 
  my $n_ab = scalar @ab;
  my $n_kmer = $n_ab ** $k; 
  for(my $i=0; $i<$n_kmer; $i++) { 
    my @tmp = (0) x $k; 
    my $j = $i; 
    my $p = $k - 1;
    while($j>0) { 
      $tmp[$p] = $j % $n_ab;
      $j = int($j/$n_ab); 
      $p --; 
    }
    my $mer = join "", @tmp; 
    $mer =~ tr/0123/ACTG/;
    push @ret, $mer;
  }
  return \@ret; 
}

