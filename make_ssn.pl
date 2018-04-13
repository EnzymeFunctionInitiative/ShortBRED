#!/bin/env perl

use strict;
use warnings;

use XML::LibXML;
use XML::Writer;
use XML::LibXML::Reader;
use Getopt::Long;


my ($ssnIn, $ssnOut, $markerFile, $proteinFile, $clusterFile);
my $result = GetOptions(
    "ssn-in=s"          => \$ssnIn,
    "ssn-out=s"         => \$ssnOut,
    "marker-file=s"     => \$markerFile,
    "protein-file=s"    => \$proteinFile,
    "cluster-file=s"    => \$clusterFile,
);

my $usage = <<USAGE;
$0 -ssn-in path_to_input_ssn -ssn-out path_to_output_ssn -marker-file path_to_shortbred_marker_file
    [-protein-file path_to_protein_abundance_file -cluster-file path_to_cluster_abundance_file]
USAGE

die $usage if not defined $ssnIn or not -f $ssnIn or not defined $ssnOut or not $ssnOut or not defined $markerFile or not -f $markerFile;


my $markerData = getMarkerData($markerFile);
my $abData = getAbundanceData($proteinFile, $clusterFile);

my $reader = XML::LibXML::Reader->new(location => $ssnIn);
my ($title, $nodes, $edges) = getNodesAndEdges($reader);

my $output = new IO::File(">$ssnOut");
my $writer = new XML::Writer(DATA_MODE => 'true', DATA_INDENT => 2, OUTPUT => $output);

writeMarkerSsn($nodes, $edges, $writer, $title, $markerData, $abData);




sub getAbundanceData {
    my $protFile = shift;
    my $clustFile = shift;

    my $abd = {metagenomes => [], proteins => {}, clusters => {}};

    if (defined $protFile and -f $protFile) {
        open PROT, $protFile or die "Unable to open protein file $protFile: $!";

        my $header = <PROT>;
        chomp($header);
        my ($feat, @mg) = split(m/\t/, $header);
        push(@{$abd->{metagenomes}}, @mg);

        while (<PROT>) {
            chomp;
            my ($feature, @mgRes) = split(m/\t/);
            $feature =~ s/^.*\|//;
            for (my $i = $#mgRes; $i < $#mg; $i++) { # Ensure that there are the same amount of results as metagenome headers
                push(@mgRes, 0);
            }
            push(@{$abd->{proteins}->{$feature}}, @mgRes);
        }

        close PROT;
    }

    if (defined $clustFile and -f $clustFile) {
        open CLUST, $clustFile or die "Unable to open cluster file $clustFile: $!";

        my $header = <CLUST>;
        chomp($header);
        my ($feat, @mg) = split(m/\t/, $header);
        push(@{$abd->{metagenomes}}, @mg) if not scalar @{$abd->{metagenomes}};

        while (<CLUST>) {
            chomp;
            my ($feature, @mgRes) = split(m/\t/);
            for (my $i = $#mgRes; $i < $#mg; $i++) { # Ensure that there are the same amount of results as metagenome headers
                push(@mgRes, 0);
            }
            push(@{$abd->{clusters}->{$feature}}, @mgRes);
        }

        close CLUST;
    }

    return $abd;
}


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
    my $abData = shift;

    $writer->startTag('graph', 'label' => $title . " ShortBRED markers", 'xmlns' => 'http://www.cs.rpi.edu/XGMML');
    writeMarkerSsnNodes($nodes, $writer, $markerData, $abData);
    writeMarkerSsnEdges($edges, $writer, $markerData);
    $writer->endTag(); 
}


sub writeMarkerSsnNodes {
    my $nodes = shift;
    my $writer = shift;
    my $markerData = shift;
    my $abd = shift; # Protein and cluster abundance data

    foreach my $node (@{$nodes}) {
        my $nodeId = $node->getAttribute('label');

        $writer->startTag('node', 'id' => $nodeId, 'label' => $nodeId);

        my $isMarker = exists $markerData->{$nodeId} ? "true" : "false";
        my $markerType = exists $markerData->{$nodeId} ? $markerData->{$nodeId}->{type} : "";
        my $isTrueMarker = $markerType eq "TM";
        $markerType = $markerType eq "TM" ? "True" : $markerType eq "JM" ? "Junction" : $markerType eq "QM" ? "Quasi" : "";
        my $markerCount = not $isTrueMarker ? 0 : exists $markerData->{$nodeId} ? $markerData->{$nodeId}->{count} : "";

        my @protValues;
        if (exists $abd->{proteins}->{$nodeId}) {
            @protValues = @{$abd->{proteins}->{$nodeId}};
        }

        #TODO: include cluster data somehow

        my $savedAttrs = 0;

        foreach my $attribute ($node->getChildnodes) {
            if ($attribute=~/^\s+$/) {
                #print "\t badattribute: $attribute:\n";
                #the parser is returning newline xml fields, this removes it
                #code will break if we do not remove it.
            } else {
                my $attrType = $attribute->getAttribute('type');
                my $attrName = $attribute->getAttribute('name');
                #if ($attrName eq "Organism") { #TODO: need to make this a shared constant
                #    writeGnnField($writer, 'Is Marker', 'string', $isMarker);
                #    writeGnnField($writer, 'Marker Type', 'string', $markerType) if $markerType;
                #    writeGnnField($writer, 'True Marker Count', 'integer', $markerCount) if $isTrueMarker;
                #    $savedAttrs = 1;
                #}

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

        if (not $savedAttrs) {
            writeGnnField($writer, 'Is Marker', 'string', $isMarker);
            writeGnnField($writer, 'Marker Type', 'string', $markerType) if $markerType;
            writeGnnField($writer, 'True Marker Count', 'integer', $markerCount) if $isTrueMarker;

            for (my $i = 0; $i <= $#{ $abd->{metagenomes} }; $i++) {
                writeGnnField($writer, $abd->{metagenomes}->[$i], 'real', $protValues[$i]);
            }
        }

        $writer->endTag(  );
    }
}


sub writeMarkerSsnEdges {
    my $edges = shift;
    my $writer = shift;

    foreach my $edge (@{$edges}) {
        $writer->startTag('edge', 'id' => $edge->getAttribute('id'), 'label' => $edge->getAttribute('label'), 'source' => $edge->getAttribute('source'), 'target' => $edge->getAttribute('target'));
        foreach my $attribute ($edge->getElementsByTagName('att')) {
            if ($attribute->getAttribute('name') eq 'interaction' or $attribute->getAttribute('name')=~/rep-net/) {
                #print "do nothing\n";
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

