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
use File::Slurp;

use lib $FindBin::Bin . "/lib";
use ShortBRED qw(expandMetanodeIds getClusterNumber);

my ($ssn, $seqFullFile, $accUniqFile, $clusterFile, $cdhit100File, $cdhitSbFile, $minSeqLen, $maxSeqLen, $markerFile, $metadataFile);
my ($spClusterFile, $spSingleFile, $clusterSizeFile, $accFullFile);

my $result = GetOptions(
    "ssn=s"                 => \$ssn,
    "sequence-full=s"       => \$seqFullFile, # list of accession IDs after min/max length filters but before uniquing
    "accession-unique=s"    => \$accUniqFile, # list of accession IDs after min/max length filter and CD-HIT 100 uniquing
    "accession-full=s"      => \$accFullFile, # list of accession IDs after min/max length filter and CD-HIT 100 uniquing
    "cluster=s"             => \$clusterFile,
    "cdhit-100=s"           => \$cdhit100File, # output of CD-HIT 100 process (unique sequences)
    "cdhit-sb=s"            => \$cdhitSbFile, # output of ShortBRED CD-HIT process (# CD-HIT clusters)
    "min-seq-len=i"         => \$maxSeqLen, 
    "max-seq-len=i"         => \$maxSeqLen, 
    "markers=s"             => \$markerFile,
    "metadata=s"            => \$metadataFile,
    "cluster-size=s"        => \$clusterSizeFile,
    "swissprot-cluster=s"   => \$spClusterFile,
    "swissprot-singleton=s" => \$spSingleFile,
);

print "missing ssn\n" and exit(0) if not $ssn or not -f $ssn;
print "missing seqFullFile\n" and exit(0) if not $seqFullFile or not -f $seqFullFile;
print "missing accFullFile\n" and exit(0) if not $accFullFile or not -f $accFullFile;
print "missing accUniqFile\n" and exit(0) if not $accUniqFile or not -f $accUniqFile;
print "missing clusterFile\n" and exit(0) if not $clusterFile or not -f $clusterFile;
print "missing cdhit100File\n" and exit(0) if not $cdhit100File or not -f $cdhit100File;
print "missing cdhitSbFile\n" and exit(0) if not $cdhitSbFile or not -f $cdhitSbFile;
print "missing markerFile\n" and exit(0) if not $markerFile or not -f $markerFile;
print "missing metadataFile\n" and exit(0) if not $metadataFile;
print "missing clusterSizeFile\n" and exit(0) if not $clusterSizeFile;
print "missing spClusterFile\n" and exit(0) if not $spClusterFile;
print "missing spSingleFile\n" and exit(0) if not $spSingleFile;


$minSeqLen = "none" if not defined $minSeqLen or not $minSeqLen;
$maxSeqLen = "none" if not defined $maxSeqLen or not $maxSeqLen;

# Field for reading if node is SwissProt annotated.
my $spKey = "Swissprot Description";



my $metadata = {
    num_metanodes => 0,
    num_raw_accessions => 0,
    min_seq_len => $minSeqLen,
    max_seq_len => $maxSeqLen,
    num_filtered_seq => 0,
    num_unique_seq => 0,
    num_cdhit_clusters => 0,
    num_ssn_clusters => 0,
    num_ssn_singletons => 0,
    num_markers => 0,
};


my $efiAnnoUtil = new EFI::Annotations;
my $reader = XML::LibXML::Reader->new(location => $ssn);
my ($clusterSize, $spStatus) = countSsnAccessions($reader, $metadata);

my $numRawSeq = `wc -l < $accFullFile`;
chomp $numRawSeq;
$metadata->{num_raw_accessions} = $numRawSeq;

my $numFiltSeq = `grep '^>' $seqFullFile | wc -l`;
chomp $numFiltSeq;
$metadata->{num_filtered_seq} = $numFiltSeq;

my $numUniqSeq = `wc -l < $accUniqFile`;
chomp $numUniqSeq;
$metadata->{num_unique_seq} = $numUniqSeq;

my $numSbClusters = `grep '^>' $cdhitSbFile | wc -l`;
chomp $numSbClusters;
$metadata->{num_cdhit_clusters} = $numSbClusters;

my $numMarkers = `grep '^>' $markerFile | wc -l`;
chomp $numMarkers;
$metadata->{num_markers} = $numMarkers;



open METADATA, ">", $metadataFile or die "Unable to write to metadata file $metadataFile: $!";
foreach my $key (keys %$metadata) {
    print METADATA "$key\t", $metadata->{$key}, "\n";
}
close METADATA;



open CLUSTERSIZE, ">", $clusterSizeFile or die "Unable to write to cluster size file $clusterSizeFile: $!";
print CLUSTERSIZE "Cluster Number\tCluster Sequence Count\n";
#my @clusterIds = sort {$clusterSize->{$b} <=> $clusterSize->{$a}} keys %$clusterSize;
my @clusterIds = sort { $a <=> $b } keys %$clusterSize;
foreach my $id (@clusterIds) {
    print CLUSTERSIZE "$id\t" . $clusterSize->{$id} . "\n";
}
close CLUSTERSIZE;



open SPCLUSTER, ">", $spClusterFile or die "Unable to write to swissprot cluster file $spClusterFile: $!";
open SPSINGLE, ">", $spSingleFile or die "Unable to write to swissprot singleton file $spSingleFile: $!";

print SPCLUSTER join("\t", "Cluster Number", "Protein ID", "SwissProt Annotation"), "\n";
print SPSINGLE join("\t", "Cluster Number", "Protein ID", "SwissProt Annotation"), "\n";

@clusterIds = sort { (my $aa = $a) =~ s/\D//g; (my $bb = $b) =~ s/\D//g; $aa <=> $bb } keys %$spStatus;
foreach my $cid (@clusterIds) {
    my $fh = $cid =~ m/^\d/ ? \*SPCLUSTER : \*SPSINGLE;
    foreach my $id (sort keys %{ $spStatus->{$cid} }) {
        $fh->print(join("\t", $cid, $id, $spStatus->{$cid}->{$id}), "\n");
    }
}

close SPSINGLE;
close SPCLUSTER;









# Gets the node and edge objects, as well as writes any sequences in the XGMML to the sequence file.
sub countSsnAccessions {
    my $reader = shift;
    my $metadata = shift;

    my $parser = XML::LibXML->new();
    $reader->read();
    my $firstnode = $reader->nextElement();

#    my @nodeIds;
    my %spStatus;
    my %clusterSize;
    my $singleCount = 0;

    my %degrees;
    do {
        next if ($reader->nodeType == 8);

        my $tmpstring = $reader->readOuterXml;
        my $tmpnode = $parser->parse_string($tmpstring);
        my $xmlNode = $tmpnode->firstChild;
        if ($reader->name() eq "node") {
            my $nodeId = $xmlNode->getAttribute("label");
            my @expandedIds = expandMetanodeIds($nodeId, $xmlNode, $efiAnnoUtil);
            push @expandedIds, $nodeId if not grep /$nodeId/, @expandedIds;
            $metadata->{num_metanodes}++;
        
            $metadata->{is_uniref} = checkUniRef($xmlNode) if not exists $metadata->{is_uniref};
            
            my $clusterId = getClusterNumber($nodeId, $xmlNode);
            $singleCount++ if not $clusterId or $clusterId =~ m/^S/;
            $clusterSize{$clusterId} += scalar @expandedIds if $clusterId and $clusterId =~ m/^\d/;
            my $status = getSwissProtDesc($xmlNode);
            my $clSizeId = (not $clusterId) ? "Singletons" : $clusterId;
            $spStatus{$clSizeId}->{$nodeId} = $status if $status;
        }
    } while ($reader->nextSiblingElement());

# Somehow this doesn't work properly in all situations    $metadata->{num_raw_accessions} = scalar @nodeIds;

    my @clusterIds = keys %clusterSize;
    $metadata->{num_ssn_clusters} = scalar @clusterIds;
    $metadata->{num_ssn_singletons} = $singleCount;

    return \%clusterSize, \%spStatus;
}



# Returns the UniRef version (e.g. 50 or 90) that the SSN was generated with, or 0 if UniRef was not used.
sub checkUniRef {
    my $xmlNode = shift;

    my $urVersion = 0;

    my @annotations = $xmlNode->findnodes("./*");
    foreach my $annotation (@annotations) {
        my $attrName = $annotation->getAttribute("name");
        if ($attrName =~ m/UniRef(\d+)/) {
            $urVersion = $1;
            last;
        }
    }

    return $urVersion;
}


# Returns the SwissProt description, if any
sub getSwissProtDesc {
    my $xmlNode = shift;

    my $spStatus = "";

    my @annotations = $xmlNode->findnodes("./*");
    foreach my $annotation (@annotations) {
        my $attrName = $annotation->getAttribute("name");
        if ($attrName eq $spKey) {
            my $attrType = $annotation->getAttribute("type");

            if ($attrType and $attrType eq "list") {
                $spStatus = getSwissProtDesc($annotation);
            } else {
                my $val = $annotation->getAttribute("value");
                $spStatus = $val if $val and length $val > 3; # Skip NA and N/A
            }

            last if $spStatus;
        }
    }

    return $spStatus;
}

