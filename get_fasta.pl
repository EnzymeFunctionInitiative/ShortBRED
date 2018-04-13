#!/usr/bin/env perl



use strict;
use warnings;

use Capture::Tiny qw(:all);
use Getopt::Long;
use File::Which qw(which);

my ($result, $idFile, $outputFile, $blastDbPath);
$result = GetOptions(
    "id-file=s"     => \$idFile,
    "output=s"      => \$outputFile,
    "blast-db=s"    => \$blastDbPath,
);

my $usage = "$0 -id-file=path_to_file_containing_ids -output=output_fasta_file";

die $usage if not defined $idFile or not -f $idFile or
              not defined $blastDbPath or not $blastDbPath;


my $isBlast1 = which("fastacmd");
my $fastaCmd = $isBlast1 ? "fastacmd" : "blastdbcmd";
my $dbFlag = $isBlast1 ? "-d" : "-db";
my $entryFlag = $isBlast1 ? "-s" : "-entry";


my @ids = get_ids($idFile);

open OUTPUT, ">$outputFile" or die "Unable to open output file $outputFile: $!";


while (scalar @ids) {
    my $batchLine = join(",", splice(@ids, 0, 1000));
    my ($fastacmdOutput, $fastaErr) = capture {
        system($fastaCmd, $dbFlag, $blastDbPath, $entryFlag, $batchLine);
    };
    my @sequences = split /\n>/, $fastacmdOutput;
    $sequences[0] = substr($sequences[0], 1) if $#sequences >= 0 and substr($sequences[0], 0, 1) eq ">";
    foreach my $seq (@sequences) {
        if (($isBlast1 and $seq =~ s/^\w\w\|(\w{6,10})\|.*//) or (not $isBlast1 and $seq =~ s/^(\w{6,10})\s.*//)) {
            my $accession = $1;
            print OUTPUT ">$accession$seq\n\n";
        }
    }
}

close OUTPUT;












sub get_ids {
    my $file = shift;

    my @ids;

    open FILE, $file or die "Unable to read id-file $file: $!";

    while (<FILE>) {
        chomp;
        my ($id, @rest) = split(m/[\s,]+/);
        push @ids, $id;
    }

    close FILE;

    return @ids;
}



