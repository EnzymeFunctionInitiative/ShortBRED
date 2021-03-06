#!/usr/bin/env perl

#TODO: This is basically a copy of the remove_demuxed_nodes.pl script in the EST repo.  We need to centralize it.
#
# This removes all nodes from the input ID/cluster files that have been filtered out by the cd-hit
# process.

BEGIN {
    die "Please load efishared before runing this script" if not $ENV{EFISHARED};
    use lib $ENV{EFISHARED};
}

use strict;
use warnings;

use Getopt::Long;
use EFI::CdHit::Parser;


my ($accIn, $clusterIn, $accOut, $clusterOut, $cdhitFile);
my $result = GetOptions(
    "id-in=s"           => \$accIn,
    "cluster-in=s"      => \$clusterIn,
    "id-out=s"          => \$accOut,
    "cluster-out=s"     => \$clusterOut,
    "cdhit-file=s"      => \$cdhitFile,
);

die "-id-in parameter required" if not defined $accIn or not -f $accIn;
die "-cluster-in parameter required" if not defined $clusterIn or not -f $clusterIn;
die "-id-out parameter required" if not defined $accOut or not $accOut;
die "-cluster-out parameter required" if not defined $clusterOut or not $clusterOut;
die "-cdhit-file parameter $cdhitFile required" if not defined $cdhitFile or not -f $cdhitFile;



my $cp = new EFI::CdHit::Parser();

#parse cluster file to get parent/child sequence associations
open CDHIT, $cdhitFile or die "cannot open cdhit cluster file $cdhitFile\n";

print "Read in clusters\n";
my $line = "";
while (<CDHIT>) {
    $line = $_;
    chomp $line;
    $cp->parse_line($line);
}
$cp->finish;

close CDHIT;



my %remove;
foreach my $clusterId ($cp->get_clusters) {
    foreach my $child ($cp->get_children($clusterId)) {
        $remove{$child} = 1 if $child ne $clusterId;
    }
}

open INID, $accIn or die "Cannot open input ID list file $accIn: $!";
open OUTID, ">$accOut" or die "Cannot open output ID list file $accOut: $!";

while (<INID>) {
    chomp;
    print OUTID "$_\n" if not exists $remove{$_};
}

close OUTID;
close INID;


open INCLUSTER, $clusterIn or die "Cannot open input cluster list file $clusterIn: $!";
open OUTCLUSTER, ">$clusterOut" or die "Cannot open output cluster list file $clusterOut: $!";

while (<INCLUSTER>) {
    chomp;
    my ($clusterNum, $accId, @rest) = split(m/\t/);
    print OUTCLUSTER "$_\n" if not exists $remove{$accId};
}

close OUTCLUSTER;
close INCLUSTER;


