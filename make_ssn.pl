#!/bin/env perl

BEGIN {
    die "Please load efishared before runing this script" if not $ENV{EFISHARED};
    use lib $ENV{EFISHARED};
}

use strict;
use warnings;

use XML::LibXML;
use XML::Writer;
use XML::LibXML::Reader;
use Getopt::Long;
use FindBin;

use lib $FindBin::Bin . "/lib";
use ShortBRED qw(getAbundanceData expandUniRefIds);
use EFI::Annotations;


my ($ssnIn, $ssnOut, $markerFile, $proteinFile, $clusterFile, $clusterMapFile, $isQuantify);
my $result = GetOptions(
    "ssn-in=s"              => \$ssnIn,
    "ssn-out=s"             => \$ssnOut,
    "marker-file=s"         => \$markerFile,
    "protein-file=s"        => \$proteinFile,
    "cluster-file=s"        => \$clusterFile,
    "cluster-map=s"         => \$clusterMapFile,
    "merge-ssn|quantify"    => \$isQuantify,
);

my $usage = <<USAGE;
$0 -ssn-in path_to_input_ssn -ssn-out path_to_output_ssn
    [-marker-file path_to_shortbred_marker_file -cluster-map path_to_cluster-protein_map_file]  <-- for identify step
    [-protein-file path_to_protein_abundance_file -cluster-file path_to_cluster_abundance_file -quantify] <-- for quantify step
USAGE

die $usage if not defined $ssnIn or not -f $ssnIn or not defined $ssnOut or not $ssnOut;
die $usage if not defined $isQuantify and (not defined $markerFile or not -f $markerFile);
die $usage if defined $isQuantify and (not defined $proteinFile or not -f $proteinFile or not defined $clusterFile or not -f $clusterFile);


$isQuantify = 0 if not defined $isQuantify;

my $efiAnnoUtil = new EFI::Annotations;

my $markerData = {};
my $clusterMap = {};
my $abData = {};

# Only get marker data and cluster map if we're generating the data from the initial step
if (not $isQuantify) {
    $markerData = getMarkerData($markerFile);
    $clusterMap = getClusterMap($clusterMapFile);
} else {
    $abData = getAbundanceData($proteinFile, $clusterFile, 0, 1); # cleanIds = no, use merged data (cluster num in separate column) = yes
}

my $reader = XML::LibXML::Reader->new(location => $ssnIn);
my ($title, $nodes, $edges) = getNodesAndEdges($reader);

my $output = new IO::File(">$ssnOut");
my $writer = new XML::Writer(DATA_MODE => 'true', DATA_INDENT => 2, OUTPUT => $output);

writeMarkerSsn($nodes, $edges, $writer, $title, $markerData, $clusterMap, $abData);






sub getMarkerData {
    my $file = shift;

    my $markerData = {};

    open FH, $file or die "Unable to open marker file $file: $!";

    while (<FH>) {
        chomp;
        if (m/^>/) {
            my $header = $_;
            (my $type = $header) =~ s/^.*_([TJQ]M)[0-9]*_.*$/$1/;
            $header =~ s/^>(tr|sp)\|/>/;
            if ($header =~ m/^>([A-Z0-9]{6,})/) {
                my $id = $1;
                $markerData->{$id} = {count => 0, type => $type} if not exists $markerData->{$id};
                $markerData->{$id}->{count}++;
            }
        }
    }

    close FH;

    return $markerData;
}


sub getClusterMap {
    my $file = shift;

    my $data = {};
    
    if (not defined $file or not -f $file) {
        return $data;
    }

    open FILE, $file;

    while (<FILE>) {
        chomp;
        my ($cluster, $protein) = split(m/\t/);
        $cluster = "N/A" if not defined $cluster or not length $cluster;
        $data->{$protein} = $cluster if defined $protein and $protein;
    }

    close FILE;

    return $data;
}



sub getNodesAndEdges{
    my $reader = shift;

    my @nodes;
    my @edges;
    my $parser = XML::LibXML->new();
    $reader->read();
    if ($reader->nodeType == 8) { #node type 8 is a comment
        print "XGMML made with ".$reader->value."\n";
        $reader->read; #we do not want to start reading a comment
    }
    my $graphname = $reader->getAttribute('label');
    my $firstnode = $reader->nextElement();
    my $tmpstring = $reader->readOuterXml;
    my $tmpnode = $parser->parse_string($tmpstring);
    my $node = $tmpnode->firstChild;
    
    if ($reader->name() eq "node") {
        push @nodes, $node;
    }

    while ($reader->nextSiblingElement()) {
        $tmpstring = $reader->readOuterXml;
        $tmpnode = $parser->parse_string($tmpstring);
        $node = $tmpnode->firstChild;
        if ($reader->name() eq "node") {
            push @nodes, $node;
        } elsif ($reader->name() eq "edge") {
            push @edges, $node;
        } else {
            warn "not a node or an edge\n $tmpstring\n";
        }
    }
    return ($graphname, \@nodes, \@edges);
}


sub writeMarkerSsn {
    my $nodes = shift;
    my $edges = shift;
    my $writer = shift;
    my $title = shift;
    my $markerData = shift;
    my $clusterMap = shift;
    my $abData = shift;

    $writer->startTag('graph', 'label' => $title . " ShortBRED markers", 'xmlns' => 'http://www.cs.rpi.edu/XGMML');
    writeMarkerSsnNodes($nodes, $writer, $markerData, $clusterMap, $abData);
    writeMarkerSsnEdges($edges, $writer, $markerData);
    $writer->endTag(); 
}


sub writeMarkerSsnNodes {
    my $nodes = shift;
    my $writer = shift;
    my $markerData = shift;
    my $clusterMap = shift;
    my $abd = shift; # Protein and cluster abundance dataa

    foreach my $node (@{$nodes}) {
        my $nodeId = $node->getAttribute('label');

        $writer->startTag('node', 'id' => $nodeId, 'label' => $nodeId);

        foreach my $attribute ($node->getChildnodes) {
            if ($attribute=~/^\s+$/) {
                #print "\t badattribute: $attribute:\n";
                #the parser is returning newline xml fields, this removes it
                #code will break if we do not remove it.
            } else {
                my $attrType = $attribute->getAttribute('type');
                my $attrName = $attribute->getAttribute('name');

                if ($attrType eq 'list') {
                    $writer->startTag('att', 'type' => $attrType, 'name' => $attrName);
                    foreach my $listelement ($attribute->getElementsByTagName('att')) {
                        $writer->emptyTag('att', 'type' => $listelement->getAttribute('type'),
                                          'name' => $listelement->getAttribute('name'),
                                          'value' => $listelement->getAttribute('value'));
                    }
                    $writer->endTag;
                } elsif ($attrName eq 'interaction') {
                    #do nothing
                    #this tag causes problems and it is not needed, so we do not include it
                } else {
                    if (defined $attribute->getAttribute('value')) {
                        $writer->emptyTag('att', 'type' => $attrType, 'name' => $attrName,
                                          'value' => $attribute->getAttribute('value'));
                    } else {
                        $writer->emptyTag('att', 'type' => $attrType, 'name' => $attrName);
                    }
                }
            }
        }

        writeResults($writer, $nodeId, $node, $markerData, $clusterMap, $abd);
        
        $writer->endTag(  );
    }
}


sub writeResults {
    my $writer = shift;
    my $nodeId = shift;
    my $node = shift;
    my $markerData = shift;
    my $clusterMap = shift;
    my $abd = shift; # Protein and cluster abundance dataa

    my @xids = expandUniRefIds($nodeId, $node, $efiAnnoUtil);
    push @xids, $nodeId; # xids = expanded IDs

    if ($isQuantify) {
        writeQuantifyResults($writer, $abd, @xids);
    } else {
        writeMarkerResults($writer, $markerData, $clusterMap, @xids);
    }
}


sub writeQuantifyResults {
    my $writer = shift;
    my $abd = shift;
    my @xids = @_;

    my $mgList = $abd->{metagenomes};

    my (%mgVals);
    foreach my $id (@xids) {
        next if not exists $abd->{proteins}->{$id};
        for (my $i = 0; $i <= $#$mgList; $i++) {
            my $mgId = $mgList->[$i];
            my $hasVal = exists $abd->{proteins}->{$id}->{$mgId};
            my $val = $hasVal ? $abd->{proteins}->{$id}->{$mgId} : "";
            push @{$mgVals{$mgId}}, $val;
        }
    }

    for (my $i = 0; $i <= $#$mgList; $i++) {
        my $mgId = $mgList->[$i];
        if (exists $mgVals{$mgId}) {
            writeGnnListField($writer, $mgId, 'real', $mgVals{$mgId});
        }
    }
}


sub writeMarkerResults {
    my $writer = shift;
    my $markerData = shift;
    my $clusterMap = shift;
    my @xids = @_;

    my (@markerTypeNames, @markerIsTrue, @markerCount, @markerIds, @markerClusters);
    foreach my $id (@xids) {
        next if not exists $markerData->{$id};
        my $markerType = $markerData->{$id}->{type};
        my $markerTypeName = $markerType eq "TM" ? "True" : $markerType eq "JM" ? "Junction" : $markerType eq "QM" ? "Quasi" : "";
        my $isTrue = $markerType eq "TM";
        my $mCount = $isTrue ? $markerData->{$id}->{count} : 0;
        push @markerTypeNames, $markerTypeName;
        push @markerIsTrue, $isTrue;
        push @markerCount, $mCount;
        push @markerIds, $id;
        my $cluster = exists $clusterMap->{$id} ? $clusterMap->{$id} : "N/A";
        push @markerClusters, $cluster;
    }

    if (scalar @markerIds) {
        writeGnnField($writer, 'Is Marker', 'string', "true");
        writeGnnListField($writer, 'Marker IDs', 'string', \@markerIds);
        writeGnnListField($writer, 'Marker Type', 'string', \@markerTypeNames);
        writeGnnListField($writer, 'True Marker Count', 'integer', \@markerCount);
        writeGnnListField($writer, 'ShortBRED SSN Cluster', 'integer', \@markerClusters);
    } else {
        writeGnnField($writer, 'Is Marker', 'string', "false");
    }
}


sub writeMarkerSsnEdges {
    my $edges = shift;
    my $writer = shift;

    foreach my $edge (@{$edges}) {
        $writer->startTag('edge', 'id' => $edge->getAttribute('id'), 'label' => $edge->getAttribute('label'), 'source' => $edge->getAttribute('source'), 'target' => $edge->getAttribute('target'));
        foreach my $attribute ($edge->getElementsByTagName('att')) {
            if ($attribute->getAttribute('name') eq 'interaction' or $attribute->getAttribute('name')=~/rep-net/) {
                #this tag causes problems and it is not needed, so we do not include it
            } else {
                $writer->emptyTag('att', 'name' => $attribute->getAttribute('name'), 'type' => $attribute->getAttribute('type'), 'value' =>$attribute->getAttribute('value'));
            }
        }
        $writer->endTag;
    }
}


sub writeGnnField {
    my $writer = shift;
    my $name = shift;
    my $type = shift;
    my $value = shift;

    unless ($type eq 'string' or $type eq 'integer' or $type eq 'real') {
        die "Invalid GNN type $type\n";
    }

    $writer->emptyTag('att', 'name' => $name, 'type' => $type, 'value' => $value);
}

sub writeGnnListField {
    my $writer = shift @_;
    my $name = shift @_;
    my $type = shift @_;
    my $valuesIn = shift @_;
    my $toSortOrNot = shift @_;

    unless($type eq 'string' or $type eq 'integer' or $type eq 'real'){
        die "Invalid GNN type $type\n";
    }
    $writer->startTag('att', 'type' => 'list', 'name' => $name);
    
    my @values;
    if (defined $toSortOrNot and $toSortOrNot) {
        @values = sort @$valuesIn;
    } else {
        @values = @$valuesIn;
    }

    foreach my $element (@values) {
        $writer->emptyTag('att', 'type' => $type, 'name' => $name, 'value' => $element);
    }
    $writer->endTag;
}

