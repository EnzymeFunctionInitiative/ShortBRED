#!/usr/bin/env perl

use strict;
use warnings;

BEGIN {
    die "Please load efishared before runing this script" if not $ENV{EFISHARED};
    use lib $ENV{EFISHARED};
}

use Getopt::Long;
use FindBin;
use EFI::SchedulerApi;
use EFI::Util qw(usesSlurm getLmod);
use EFI::CdHitParser;

use lib $FindBin::Bin . "/lib";
use ShortBRED;


die "Load efidb module in order to run this program." if not exists $ENV{EFIDBPATH};


my ($inputSsn, $outputDirName);
my ($np, $queue, $memQueue, $scheduler, $dryRun, $jobId, $outputSsnName, $cdhitFileName, $parentJobId);
my ($minSeqLen);


my $result = GetOptions(

    "ssn-in=s"          => \$inputSsn,
    "tmpdir=s"          => \$outputDirName,
    "ssn-out-name=s"    => \$outputSsnName,
    "cdhit-out-name=s"  => \$cdhitFileName,

    "parent-job-id=i"   => \$parentJobId,
    "min-seq-len=i"     => \$minSeqLen,

    "job-id=i"          => \$jobId,
    "np=i"              => \$np,
    "queue=s"           => \$queue,
    "mem-queue=s"       => \$memQueue,
    "scheduler=s"       => \$scheduler,
    "dryrun"            => \$dryRun,
);

my $usage =
"$0 -ssn-in=ssn_file -tmpdir=relative_output_dir -queue=cluster_queue -ssn-out-name output_name
    [-job-id=job_id -np=num_cpu -scheduler=scheduler_type -dryrun]

    -ssn-in         path to input SSN file (relative or absolute)
    -ssn-out-name   what to name the output xgmml file
    -cdhit-out-name what to name the output cdhit mapping table file
    -tmpdir         name of output directory (relative)
    -job-id         number of the job [optional]
    -np             number of CPUs to use [optional, defaults to 24]
    -queue          the name of the cluster queue to submit to
    -scheduler      the type of scheduler to use (torque or slurm) [optional, defaults to slurm]
    -dryrun         only output the commands to be executed, do not submit to cluster [optional]
";

die $usage if not defined $inputSsn or not -f $inputSsn or not defined $outputDirName or not $outputDirName or
              not defined $queue or not $queue or not defined $outputSsnName or not $outputSsnName;


$memQueue = $queue if not defined $memQueue or not $memQueue;
$np = 24 if not defined $np or not $np;
$scheduler = "slurm" if not defined $scheduler or not $scheduler;
$dryRun = 0 if not defined $dryRun;
$parentJobId = 0 if not defined $parentJobId;
$minSeqLen = 0 if not defined $minSeqLen;

my ($outputDir, $parentOutputDir, $scriptDir, $logDir) = initDirectoryStructure($parentJobId);
my ($S, $jobNamePrefix) = initScheduler($logDir);

# If we are linking to a parent job, most of the files we are using are in the parent output dir.
my $childOutputDir = $outputDir;
if (defined $parentJobId and $parentJobId) {
    $outputDir = $parentOutputDir;
}

my $pythonMod = getLmod("Python/2", "Python");
my $biopythonMod = getLmod("Biopython", "Biopython");
my $usearchMod = getLmod("USEARCH/9", "USEARCH");
my $muscleMod = getLmod("MUSCLE/3", "MUSCLE");
my $blastMod = getLmod("BLAST+", "BLAST+");
my $cdhitMod = getLmod("CD-HIT", "CD-HIT");
my $blastDbPath = "$ENV{EFIDBPATH}/combined.fasta";
my $dbSupport = "$ENV{EFIDBHOME}/support";
my $dbModule = $ENV{EFIDBMOD};
my $sbModule = $ENV{SHORTBRED_MOD};
my $sbParseSsnApp = $ENV{SHORTBRED_PARSESSN};
my $sbIdentifyApp = $ENV{SHORTBRED_IDENTIFY};

my $efiSbDir = $ENV{EFI_SHORTBRED_HOME};
my $ssnAccessionFile = "$outputDir/accession";
my $ssnClusterFile = "$outputDir/cluster";
my $ssnSequenceFile = "$outputDir/ssn-sequences.fa";
my $fastaFile = "$outputDir/sequences.fa";
my $sbOutputDir = "$outputDir/id-temp";
my $sbMarkerFile = "$outputDir/markers.faa";
my $cdhitFile = "$sbOutputDir/clust/clust.faa.clstr";
my $cdhitTableFile = (defined $cdhitFileName and $cdhitFileName) ? "$outputDir/$cdhitFileName" : "$outputDir/cdhit.txt";
my $colorFile = -f "$dbSupport/colors.tab" ? "$dbSupport/colors.tab" : "";
my $clusterSizeFile = "$outputDir/cluster.sizes";
my $ssnMarker = "$outputDir/$outputSsnName";
my $jobPrefix = (defined $jobId and $jobId) ? "${jobId}_" : "";
my $submitName = "";
my $submitFile = "";
my $depId = 0;

my $ssnErrorDir = $outputDir;
if ($parentJobId) {
    $ssnErrorDir = $childOutputDir;
    $ssnAccessionFile = "$childOutputDir/accession";
    $ssnClusterFile = "$childOutputDir/cluster";
    $cdhitTableFile = (defined $cdhitFileName and $cdhitFileName) ? "$childOutputDir/$cdhitFileName" : "$childOutputDir/cdhit.txt";
    $ssnMarker = "$childOutputDir/$outputSsnName";
}


my $B;


#######################################################################################################################
# Unzip the file if necessary
my $inputSsnZip = "";
if ($inputSsn =~ m/\.zip/i) {
    $inputSsnZip = $inputSsn;
    $inputSsn =~ s/\.zip/.xgmml/i;
}


#######################################################################################################################
# Get the clusters and accessions
my $minSeqLenArg = $minSeqLen ? "-min-seq-len $minSeqLen" : "";

$B = $S->getBuilder();
$submitName = "sb_get_clusters";
$B->setScriptAbortOnError(0); # grep causes the script to abort if we have set -e in the script.
$B->resource(1, 1, "70gb");
$B->addAction("module load $sbModule");
$B->addAction("$efiSbDir/unzip_file.pl -in $inputSsnZip -out $inputSsn") if $inputSsnZip =~ m/\.zip$/i;
$B->addAction("HASCLUSTERNUM=`head -2000 $inputSsn | grep -m1 -e \"Cluster Number\" -e \"Singleton Number\"`");
$B->addAction("if [[ \$HASCLUSTERNUM == \"\" ]]; then");
$B->addAction("    echo \"ERROR: Cluster Number is not present in SSN\"");
$B->addAction("    touch $ssnErrorDir/ssn_cl_num.failed");
$B->addAction("    exit 1");
$B->addAction("fi");
$B->addAction("$efiSbDir/get_clusters.pl -ssn $inputSsn -accession-file $ssnAccessionFile -cluster-file $ssnClusterFile -sequence-file $ssnSequenceFile $minSeqLenArg");
# Add this check because we disable set -e above for grep.
$B->addAction("if [ $? != 0 ]; then");
$B->addAction("    echo \"ERROR: in get_clusters.pl\"");
$B->addAction("    exit 1");
$B->addAction("fi");
$depId = doSubmit();


if (not $parentJobId) {
    # CD-HIT params
    my $lenDiffCutoff = "1";
    my $seqIdCutoff = "1";
    my $tempFasta = "$fastaFile.cdhit100";
    my $tempAcc = "$ssnAccessionFile.cdhit100";
    my $tempCluster = "$ssnClusterFile.cdhit100";
    my $sortedAcc = "$ssnAccessionFile.sorted";

    #######################################################################################################################
    # Get the FASTA files from the database
    
    $B = $S->getBuilder();
    $submitName = "sb_get_fasta";
    $B->resource(1, 1, "15gb");
    $B->addAction("module load $sbModule");
    $B->addAction("module load $dbModule");
    $B->addAction("sort $ssnAccessionFile > $sortedAcc");
    $B->addAction("$efiSbDir/get_fasta.pl -id-file $sortedAcc -output $fastaFile -blast-db $blastDbPath $minSeqLenArg");
    $B->addAction("SZ=`stat -c%s $ssnSequenceFile`");
    $B->addAction("if [[ \$SZ != 0 ]]; then");
    $B->addAction("    cat $ssnSequenceFile >> $fastaFile");
    $B->addAction("fi");
    $B->addAction("cd-hit -c $seqIdCutoff -s $lenDiffCutoff -i $fastaFile -o $tempFasta -M 14900");
    $B->addAction("$efiSbDir/remove_redundant_sequences.pl -id-in $sortedAcc -cluster-in $ssnClusterFile -id-out $tempAcc -cluster-out $tempCluster -cdhit-file $tempFasta.clstr");
    $B->addAction("mv $fastaFile $fastaFile.full");
    $B->addAction("mv $ssnAccessionFile $ssnAccessionFile.full");
    $B->addAction("mv $ssnClusterFile $ssnClusterFile.full");
    $B->addAction("mv $tempFasta $fastaFile");
    $B->addAction("mv $tempAcc $ssnAccessionFile");
    $B->addAction("mv $tempCluster $ssnClusterFile");
    $B->addAction("SZ=`stat -c%s $ssnAccessionFile`");
    $B->addAction("if [[ \$SZ == 0 ]]; then");
    $B->addAction("    echo \"Unable to find any FASTA sequences. Check input file.\"");
    $B->addAction("    touch $outputDir/get_fasta.failed");
    $B->addAction("    exit 1");
    $B->addAction("fi");
    $depId = doSubmit($depId);
    
    
    #######################################################################################################################
    # Run ShortBRED-Identify
    
    $B = $S->getBuilder();
    $submitName = "sb_identify";
    $B->resource(1, $np, "300gb");
    $B->addAction("module load $sbModule");
    $B->addAction("python $sbIdentifyApp --threads $np --goi $fastaFile --refdb $blastDbPath --markers $sbMarkerFile --tmp $sbOutputDir");
    $depId = doSubmit($depId);
}


#######################################################################################################################
# Build the XGMML file with the marker attributes added

my $colorFileArg = $colorFile ? "-color-file $colorFile" : "";
$B = $S->getBuilder();
$submitName = "sb_make_xgmml";
$B->queue($memQueue);
$B->resource(1, 1, "200gb");
$B->addAction("module load $sbModule");
#$B->addAction("$efiSbDir/find_cluster_sizes.pl -cluster-list-file $ssnClusterFile -output-file $clusterSizeFile");
$B->addAction("$efiSbDir/make_cdhit_table.pl -cdhit-file $cdhitFile -cluster-map $ssnClusterFile -table-file $cdhitTableFile $colorFileArg");
$B->addAction("$efiSbDir/make_ssn.pl -ssn-in $inputSsn -ssn-out $ssnMarker -marker-file $sbMarkerFile -cluster-map $ssnClusterFile -cdhit-file $cdhitTableFile");
$B->addAction("zip -j $ssnMarker.zip $ssnMarker");
$B->addAction("touch $childOutputDir/job.completed");
$depId = doSubmit($depId);





sub doSubmit {
    my $depJobId = shift;
    if (defined $depJobId and $depJobId) {
        $B->dependency(0, $depJobId);
    }

    $B->jobName("$jobPrefix$submitName");
    $submitFile = "$scriptDir/$submitName.sh";
    $B->renderToFile($submitFile);
    $depJobId = $S->submit($submitFile);
    
    chomp $depJobId;
    print "$submitName job ID is :\n $depJobId"; # Format is important so web scripts can read the job ID

    return $depJobId;
}



sub initScheduler {
    my ($logDir) = @_;

    # Set up the scheduler API so we can work with Torque or Slurm.
    my $schedType = "torque";
    $schedType = "slurm" if (defined($scheduler) and $scheduler eq "slurm") or (not defined($scheduler) and usesSlurm());
    
    my %schedArgs = (type => $schedType, queue => $queue, resource => [1, 1, "35gb"], dryrun => $dryRun);
    $schedArgs{output_base_dirpath} = $logDir if $logDir;
    my $S = new EFI::SchedulerApi(%schedArgs);
    my $jobNamePrefix = $jobId ? $jobId . "_" : ""; 
    
    return ($S, $jobNamePrefix);
}

sub initDirectoryStructure {
    my $parentId = shift;
    
    my $baseOutputDir = $ENV{PWD};
    my $outputDir = "$baseOutputDir/$outputDirName";
    mkdir $outputDir if not $dryRun;

    my $logDir = "$baseOutputDir/log";
    mkdir $logDir if not $dryRun;
    $logDir = "" if not -d $logDir;

    my $scriptDir = "$baseOutputDir/scripts";
    mkdir $scriptDir if not $dryRun;
    $scriptDir = $outputDir if not -d $scriptDir;

    my $parentOutputDir = "";
    if ($parentId) {
        ($parentOutputDir = $baseOutputDir) =~ s%/(\d+)/*$%/$parentId%;
        $parentOutputDir .= "/$outputDirName";
    }

    return ($outputDir, $parentOutputDir, $scriptDir, $logDir);
}



