#!/bin/env perl


use strict;
use warnings;

use XML::LibXML::Reader;
use Getopt::Long;


my ($ssn, $accFile, $clusterFile);

my $result = GetOptions(
    "ssn=s"             => \$ssn,
    "accession-file=s"  => \$accFile,
    "cluster-file=s"    => \$clusterFile,
);




my $reader = XML::LibXML::Reader->new(location => $ssn);

#my ($name, $nodes, $edges, $degrees) = getNodesAndEdges($reader);
my ($network, $nodeIds) = getNodesAndEdges($reader);
my ($clusters, $constellations) = getClusters($network, $nodeIds);

my @nodeIds = sort { scalar(@{$clusters->{$b}}) <=> scalar(@{$clusters->{$a}}) } keys %$clusters;

open CLUSTER, "> $clusterFile" or die "Unable to write to cluster file $clusterFile: $!";
open ACCESSION, "> $accFile" or die "Unable to write to accession file $accFile: $!";

my $clusterCount = 1;
foreach my $clusterId (@nodeIds) {
    foreach my $id (sort @{$clusters->{$clusterId}}) {
        print CLUSTER join("\t", $clusterCount, $id), "\n";
        print ACCESSION "$id\n";
    }
    $clusterCount++;
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

    if ($reader->name() eq "node"){
        push @nodes, $xmlNode;
    }

    my @network;
    my @nodeIds;

    my %degrees;
    while ($reader->nextSiblingElement()) {
        $tmpstring = $reader->readOuterXml;
        $tmpnode = $parser->parse_string($tmpstring);
        $xmlNode = $tmpnode->firstChild;
        if ($reader->name() eq "node") {
            push @nodes, $xmlNode;
            push @nodeIds, $xmlNode->getAttribute("label");;
        } elsif($reader->name() eq "edge") {
            push @edges, $xmlNode;
            my $source = $xmlNode->getAttribute("source");
            my $target = $xmlNode->getAttribute("target");

            push @network, {source => $source, target => $target};
#
#            $network{$source} = {visited => 0, neighbors => []};
#            push(@{ $network{$source}->{neighbors} }, $target);
#            $network{$target} = {visited => 0, neighbors => []};
#            push(@{ $network{$target}->{neighbors} }, $source);
#
#            #$degrees{$source} = 0 if not exists $degrees{$source};
#            #$degrees{$target} = 0 if not exists $degrees{$target};
#            #$degrees{$source}++;
#            #$degrees{$target}++;
#        } else {
#            warn "not a node or an edge\n $tmpstring\n";
        }
    }
#    return ($graphname, \@nodes, \@edges, \%degrees);
    return \@network, \@nodeIds;
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

    foreach my $id (keys %$nodeIds) {
        if (not exists $con{$id}) {
            $super{$clusterId} = $id;
            $con{$id} = $clusterId;
            $clusterId++;
        }
    }

    return (\%super, \%con);
}





