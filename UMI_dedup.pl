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
use List::Util qw(min max sum);
use Sort::Versions;
my $debug_lvl = 0;

#####################
my $usage = "$0 [options] <in.bed.gz> <out.bed.gz>
Options: 
  --min-UMI-cnt     [default: 1]
  --fuzz-pos-dist   [default: 5]
  --max-merge-dist  [default: 2]
  --add-num-to-UMI  [default: false]
  --help            [show this help & exit]
\n";

#### options
my $min_UMI_cnt = 1; 
my $fuzz_pos_dist = 5; 
my $max_merge_node_dist = 2;  #dist between nodes, if <= this value then nodes can be merged
my $add_no_to_UMI = 0; # add a uniq number to UMI
my $help = 0; 
## in-house
my $major_UMI_rel_frac = 0.5; # larger than this value, should keep (dist+=inf)
my $sub_major_UMI_rel_frac = 0.2; # larger than this value, should have more chance to keep (dist+=1)
my $remove_UMI_mean_frac = 0.01; # lower than this value, should be removed #TODO 
#################
GetOptions ('min-UMI-cnt=i' => \$min_UMI_cnt, 'fuzz-pos-dist=i' => \$fuzz_pos_dist, 'max-merge-dist=i' => \$max_merge_node_dist, 'add-num-to-UMI' => \$add_no_to_UMI, 'help' => \$help) || die $usage;
my $discard_UMI_cnt = $min_UMI_cnt - 1;

$help && die $usage;
my $infile = shift || die $usage;
my $outfile = shift || die $usage;

### 
my $IN = IO::Zlib->new("$infile", "rb") || die "Cannot open $infile for reading.\n";
my $OUT = IO::Zlib->new("$outfile", "wb9") || die "Cannot open $outfile for writing\n";


## read $IN & cluster alignments
my $cur_chr = "chr1";
my $cur_start = 0;
my $cur_end = 0;
my $cur_strand = "-";
my %clusters;
my @cur_cluster; 
my $key;
while(<$IN>) {
  chomp;
  my @a = split; 
  if($cur_chr eq $a[0] && $cur_strand eq $a[5] && $cur_end > $a[1]) { # overlapped with cur region
    $cur_end = $a[2] if($cur_end < $a[2]); 
    push @cur_cluster, $_; 
  } else { # not overlap
    ### put the cluster into hash
    $key = join "_", ($cur_chr, $cur_start, $cur_end, $cur_strand);
    $clusters{$key} = [ @cur_cluster ]; ## 
    ### renew array
    @cur_cluster = (); 
    ### renew cur region
    $cur_chr = $a[0];
    $cur_strand = $a[5];
    $cur_start = $a[1];
    $cur_end = $a[2];
    push @cur_cluster, $_; 
  }
}
$key = join "_", ($cur_chr, $cur_start, $cur_end, $cur_strand);
$clusters{$key} = [ @cur_cluster ]; 


## in each cluster, select UMI and output 
my $n = -1;
foreach $key (sort { versioncmp($a, $b) } keys %clusters) { 
  print STDERR $key."\n" if($debug_lvl > 1); 
  #print join "\n", @{$clusters{$key}};
  #print "\n\n";
  my @output = sort { versioncmp($a, $b) }  @{&parse_each_cluster($clusters{$key})};
  my @k_a = split /_/, $key; 
  for(my $i=0; $i<@output; $i++) { 
    print STDERR $output[$i], "\n" if($debug_lvl > 1);
    my @a = split /\t/, $output[$i]; 
    next unless($a[4] > $discard_UMI_cnt); 
    $n++;
    $a[2] = $a[2].":".$n if($add_no_to_UMI);
    if($a[3] eq "0") { 
      print $OUT join "\t", ($k_a[0], $a[0], $a[1], $a[2], $a[4], $k_a[3], $a[0], $a[1], "0", "1", $a[1]-$a[0], "0");
    } else {
      my @ba = split /,/, $a[3]; 
      my @b = ($a[0], @ba, $a[1]);
      my (@s, @l); 
      for(my $j=0; $j<@b; $j=$j+2) { 
        push @l, $b[$j+1] - $b[$j];
        push @s, $b[$j] - $b[0]; 
      }
      my $s_s = join ",", @s;
      my $l_s = join ",", @l;
      print $OUT join "\t", ($k_a[0], $a[0], $a[1], $a[2], $a[4], $k_a[3], $a[0], $a[1], "0", ($#b+1)/2, $l_s, $s_s);
    }
    print $OUT "\n";
  }
}
undef $OUT;

###########################
## subs
###########################
sub parse_each_cluster { 
  my @lines = @{$_[0]};
  my %UMI_cnt; 
  my $junc_pos;
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
    $UMI_cnt{$k} = $a[4];
  }
  
  my %nodes;     # key is the major node, value is an array storing merged nodes
  my %node_dist; # key is each node, value is the dist to major node
  my $full_dist = 999; 
  my $max_cnt = max(values %UMI_cnt);
  my $avg_cnt = mean(values %UMI_cnt);
  foreach my $k (sort {$UMI_cnt{$b} <=> $UMI_cnt{$a}} keys %UMI_cnt) { 
    my @major_nodes = sort {$UMI_cnt{$b} <=> $UMI_cnt{$a}} keys %nodes;
    if($#major_nodes == -1 || $UMI_cnt{$k} >= $max_cnt * $major_UMI_rel_frac) {
      push @{$nodes{$k}},$k; 
      $node_dist{$k} = 0; 
    } else { 
      last if($UMI_cnt{$k} < $avg_cnt * $remove_UMI_mean_frac);
      my $dist;
      my @aa = split /\t/, $k;
      for(my $i=0; $i<@major_nodes; $i++) { 
        my @all_nodes = @{$nodes{$major_nodes[$i]}}; 
        for(my $j=0; $j<@all_nodes; $j++) { 
          my @bb = split /\t/, $all_nodes[$j]; 
          $dist = $node_dist{$all_nodes[$j]}; 
          # count
          $dist += 1 if($UMI_cnt{$k} >= $max_cnt * $sub_major_UMI_rel_frac);
          # splicing pattern
          $dist += $full_dist if($aa[3] ne $bb[3]); 
          # mapped postion
          $dist += 1 if( ($aa[0] != $bb[0] && $aa[0] <= ($bb[0] + $fuzz_pos_dist) && $aa[0] >= ($bb[0] - $fuzz_pos_dist)) || ($aa[1] != $bb[1] && $aa[1] <= ($bb[1] + $fuzz_pos_dist) && $aa[1] >= ($bb[1] - $fuzz_pos_dist)) );
          $dist += $full_dist if( $aa[0] > ($bb[0] + $fuzz_pos_dist) || $aa[0] < ($bb[0] - $fuzz_pos_dist) || $aa[1] > ($bb[1] + $fuzz_pos_dist) || $aa[1] < ($bb[1] - $fuzz_pos_dist) );
          # UMI mismatch
          if($aa[2] ne $bb[2]) { 
            my @next = @{&next_kmers($aa[2])};
            my %next_hash = map { $_ => 1 } @next;
            if(exists($next_hash{$bb[2]})) {
              $dist += 1; 
            } else {
              $dist += $full_dist;
            }
          }
          if($dist <= $max_merge_node_dist) { 
            push @{$nodes{$major_nodes[$i]}}, $k; 
            $node_dist{$k} = $dist;
            last; 
          }
        }
        last if($dist <= $max_merge_node_dist);
      }
      if($dist > $max_merge_node_dist) { 
        push @{$nodes{$k}},$k;
        $node_dist{$k} = 0;
      }
    }
    print STDERR $k, "\t", $UMI_cnt{$k}, "\t", $node_dist{$k}, "\n" if($debug_lvl > 1);
  }
  ## from top count to bottom: 
  ## 1. push highest but indepdent UMI into new UMI node
  ## 2. if within postion distance and within allowed UMI barcode distance and with the same splice pattern, then the UMI is not independent. 
  
  my @ret;  
  foreach my $major (sort {$UMI_cnt{$b} <=> $UMI_cnt{$a}} keys %nodes) { 
    my $cnt = sum(@UMI_cnt{@{$nodes{$major}}});
    push @ret, $major."\t".$cnt; 
  }
  return(\@ret);
}

#my @test = @{&next_kmers("AAAAAA")};
#print join "\t", @test, "\n"; 

sub next_kmers { 
  my $kmer = $_[0]; 
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

sub mean {
   return @_ ? sum(@_) / @_ : 0;
}

## chr10_128035604_128036511_+
#
