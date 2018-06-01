#!/usr/bin/env perl

use strict;
use warnings;

BEGIN {
    die "Please load efishared before runing this script" if not $ENV{EFISHARED};
    use lib $ENV{EFISHARED};
}

use Getopt::Long;
use FindBin;
use EFI::CdHitParser;

use lib $FindBin::Bin . "/lib";
use ShortBRED qw(getClusterMap);


my ($cdhitInput, $tableOutput, $clusterMapFile, $colorFile);


my $result = GetOptions(
    "cdhit-file=s"      => \$cdhitInput,
    "table-file=s"      => \$tableOutput,
    "cluster-map=s"     => \$clusterMapFile,
    "color-file=s"      => \$colorFile,
);

my $usage =
"$0 -cdhit-file path_to_input_cdhit_file -table-file path_to_output_table [-color-file path_to_input_colors_file]
";

die $usage if not defined $cdhitInput or not -f $cdhitInput or not defined $tableOutput or not $tableOutput or
              not defined $clusterMapFile or not -f $clusterMapFile;

my $colors = {};
if (defined $colorFile and -f $colorFile) {
    $colors = getColors($colorFile);
}

my $clusterMap = getClusterMap($clusterMapFile);


my $parser = new EFI::CdHitParser();
$parser->parse_file($cdhitInput);
my @clusters = $parser->get_clusters();

open OUTPUT, "> $tableOutput" or die "Unable to open table output $tableOutput: $!";

my @headers = ("CD-HIT Seed Sequence", "Protein", "ShortBRED Cluster Number");
push(@headers, "CD-HIT Seed Sequence Color") if keys %$colors;
print OUTPUT join("\t", @headers), "\n";

my $c = 1;
foreach my $cluster (@clusters) {
    my @children = $parser->get_children($cluster);
    my $color = exists $colors->{$c} ? $colors->{$c} : "";
    $c++;
    foreach my $child (@children) {
        my $clusterNum = exists $clusterMap->{$child} ? $clusterMap->{$child} : "N/A";
        my @vals = ($cluster, $child, $clusterNum);
        push(@vals, $color) if keys %$colors;
        print OUTPUT join("\t", @vals), "\n";
    }
}

close OUTPUT;




sub getColors {
    my $file = shift;

    my $colors = {};

    open FILE, $file or warn "Unable to open color file $file: $!";

    while (<FILE>) {
        chomp;
        my ($num, $color) = split(m/\t/);
        $colors->{$num} = $color;
    }

    close FILE;

    return $colors;
}


