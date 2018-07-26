#!/bin/env perl

BEGIN {
    die "Please load efishared before runing this script" if not $ENV{EFISHARED};
    use lib $ENV{EFISHARED};
}

use strict;
use warnings;

use XML::LibXML::Reader;
use Getopt::Long;
use EFI::Annotations;
use FindBin;
use Data::Dumper;

use lib $FindBin::Bin . "/lib";
use ShortBRED qw(expandMetanodeIds getClusterNumber);

my ($ssn, $accFile, $clusterFile, $useDefaultClusterNumbering );

my $result = GetOptions(
    "ssn=s"             => \$ssn,
    "accession-file=s"  => \$accFile,
    "cluster-file=s"    => \$clusterFile,
    "default-numbering" => \$useDefaultClusterNumbering,
);


my $usage =
"$0 -ssn=path_to_ssn -accession-file=output_accession_list -cluster-file=output_cluster_list";

die $usage if not defined $ssn or not -f $ssn or not defined $accFile or not $accFile or not defined $clusterFile or not $clusterFile;

$useDefaultClusterNumbering = 0 if not defined $useDefaultClusterNumbering;
$useDefaultClusterNumbering = 1 if defined $useDefaultClusterNumbering;


my $efiAnnoUtil = new EFI::Annotations;

my $reader = XML::LibXML::Reader->new(location => $ssn);

#my ($name, $nodes, $edges, $degrees) = getNodesAndEdges($reader);
my ($network, $nodeIds, $clusterNumbers) = getNodesAndEdges($reader);
my ($clusters, $constellations) = getClusters($network, $nodeIds);


# Sort by cluster size
my @clusterIds = sort { scalar(@{$clusters->{$b}}) <=> scalar(@{$clusters->{$a}}) } keys %$clusters;

open CLUSTER, "> $clusterFile" or die "Unable to write to cluster file $clusterFile: $!";
open ACCESSION, "> $accFile" or die "Unable to write to accession file $accFile: $!";

my $clusterCount = 0;
my $singleCount = 0;
my @singles;

#foreach my $clusterId (@clusterIds) {
#    my @ids = sort @{$clusters->{$clusterId}};
#    my $isSingle = scalar @ids == 1;
#    if ($isSingle) {
#        push @singles, $ids[0];
#        next;
#    }
#    $clusterCount++;
#    foreach my $id (@ids) {
#        my $clusterNumber = $clusterCount;
#        if (exists $clusterNumbers->{$id}) {
#            $clusterNumber = $clusterNumbers->{$id};
#        }
#        print CLUSTER join("\t", $clusterNumber, $id), "\n";
#        print ACCESSION "$id\n";
#    }
#}
#
#$singleCount = $clusterCount + 1;
#foreach my $singleId (@singles) {
#    $singleCount++;
#    my $clusterNumber = $singleCount;
#    if (exists $clusterNumbers->{$id}) {
#        $clusterNumber = $clusterNumbers->{$id};
#    }
#    print CLUSTER join("\t", "S$clusterNumber", $id), "\n";
#    print ACCESSION "$id\n";
#}

foreach my $clusterId (@clusterIds) {
    my @ids = sort @{$clusters->{$clusterId}};
    my $isSingle = scalar @ids == 1;
    if ($isSingle) {
        $singleCount++;
    } else {
        $clusterCount++;
    }
    foreach my $id (@ids) {
        my $clusterNumber = $clusterCount;
        if (exists $clusterNumbers->{$id}) {
            $clusterNumber = $clusterNumbers->{$id};
#        } elsif ($isSingle) {
#            $clusterNumber = $singleCount;
        }
        print CLUSTER join("\t", $clusterNumber, $id), "\n";
        print ACCESSION "$id\n";
    }
}

close ACCESSION;
close CLUSTER;














sub getNodesAndEdges{
    my $reader = shift @_;

    my @nodes;
    my @edges;
    my $parser = XML::LibXML->new();
    $reader->read();
    if ($reader->nodeType == 8) {
        $reader->read; #we do not want to start reading a comment
    }
    my $graphname = $reader->getAttribute('label');
    my $firstnode = $reader->nextElement();
    my $tmpstring = $reader->readOuterXml;
    my $tmpnode = $parser->parse_string($tmpstring);
    my $xmlNode = $tmpnode->firstChild;

    my @network;
    my @nodeIds;
    my $clusterNumbers = {};

    if ($reader->name() eq "node"){
        push @nodes, $xmlNode;
        my $nodeId = $xmlNode->getAttribute("label");
        push @nodeIds, $nodeId;
        my @expandedIds = expandMetanodeIds($nodeId, $xmlNode, $efiAnnoUtil);
        my $clusterNum = getClusterNumber($nodeId, $xmlNode);
        if ($clusterNum) {
            map { $clusterNumbers->{$_} = $clusterNum; } (@expandedIds, $nodeId);
        }
        push @nodeIds, @expandedIds;
    }

    my %degrees;
    while ($reader->nextSiblingElement()) {
        $tmpstring = $reader->readOuterXml;
        $tmpnode = $parser->parse_string($tmpstring);
        $xmlNode = $tmpnode->firstChild;
        if ($reader->name() eq "node") {
            push @nodes, $xmlNode;
            my $nodeId = $xmlNode->getAttribute("label");
            push @nodeIds, $nodeId;
            my @expandedIds = expandMetanodeIds($nodeId, $xmlNode, $efiAnnoUtil);
            my $clusterNum = getClusterNumber($nodeId, $xmlNode);
            if ($nodeId eq "A0A101VZA3") {
                print "CLUSTER: $clusterNum\n";
            }
            if ($clusterNum) {
                map { $clusterNumbers->{$_} = $clusterNum; } (@expandedIds, $nodeId);
            }
            push @nodeIds, @expandedIds;
        } elsif($reader->name() eq "edge") {
            push @edges, $xmlNode;
            
            my $label = $xmlNode->getAttribute("label");
            my ($source, $target);
            if (defined $label) {
                ($source, $target) = split(m/,/, $label);
            }
            if (not defined $source or not $source or not defined $target or not $target) {
                $source = $xmlNode->getAttribute("source");
                $target = $xmlNode->getAttribute("target");
            }

            die "$tmpstring: " . Dumper($xmlNode) if not $source;
            die "Target" if not $target;

            push @network, {source => $source, target => $target};
        }
    }
    
    return \@network, \@nodeIds, $clusterNumbers;
}


sub getClusters {
    my $edges = shift;
    my $nodeIds = shift;

    my %con;
    my %super;

    my $clusterId = 1;

    foreach my $edge (@$edges) {
        my $source = $edge->{source};
        my $target = $edge->{target};

        if (exists $con{$source}) {
            if (exists $con{$target}) {
                next if ($con{$target} eq $con{$source});
                push @{$super{$con{$source}}}, @{$super{$con{$target}}};
                delete $super{$con{$target}};

                my $oldTarget = $con{$target};
                foreach my $node (keys %con) {
                    if ($oldTarget == $con{$node}) {
                        $con{$node} = $con{$source};
                    }
                }
            } else {
                $con{$target} = $con{$source};
                push @{$super{$con{$source}}}, $target;
            }
        } elsif (exists $con{$target}) {
            $con{$source} = $con{$target};
            push @{$super{$con{$target}}}, $source;
        } else {
            $con{$source} = $con{$target} = $clusterId;
            push @{$super{$clusterId}}, $source, $target;
            $clusterId++;
        }
    }

    foreach my $id (@$nodeIds) {
        if (not exists $con{$id}) {
            push @{$super{$clusterId}}, $id;
            $con{$id} = $clusterId;
            $clusterId++;
        }
    }

    return (\%super, \%con);
}




