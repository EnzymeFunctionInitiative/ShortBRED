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

use lib $FindBin::Bin . "/lib";
use ShortBRED;


die "Load efidb module in order to run this program." if not exists $ENV{EFIDBPATH};


my ($inputSsn, $outputDirName);
my ($np, $queue, $scheduler, $dryRun, $jobId, $outputSsnName);


my $result = GetOptions(

    "ssn-in=s"          => \$inputSsn,
    "tmpdir=s"          => \$outputDirName,
    "ssn-out-name=s"    => \$outputSsnName,

    "job-id=i"          => \$jobId,
    "np=i"              => \$np,
    "queue=s"           => \$queue,
    "scheduler=s"       => \$scheduler,
    "dryrun"            => \$dryRun,
);

my $usage =
"$0 -ssn-in=ssn_file -tmpdir=relative_output_dir -queue=cluster_queue -ssn-out-name output_name
    [-job-id=job_id -np=num_cpu -scheduler=scheduler_type -dryrun]

    -ssn-in         path to input SSN file (relative or absolute)
    -ssn-out-name   what to name the output xgmml file
    -tmpdir         name of output directory (relative)
    -job-id         number of the job [optional]
    -np             number of CPUs to use [optional, defaults to 24]
    -queue          the name of the cluster queue to submit to
    -scheduler      the type of scheduler to use (torque or slurm) [optional, defaults to slurm]
    -dryrun         only output the commands to be executed, do not submit to cluster [optional]
";

die $usage if not defined $inputSsn or not -f $inputSsn or not defined $outputDirName or not $outputDirName or
              not defined $queue or not $queue or not defined $outputSsnName or not $outputSsnName;


$np = 24 if not defined $np or not $np;
$scheduler = "slurm" if not defined $scheduler or not $scheduler;
$dryRun = 0 if not defined $dryRun;

my ($outputDir, $scriptDir, $logDir) = initDirectoryStructure();
my ($S, $jobNamePrefix) = initScheduler($logDir);

my $pythonMod = getLmod("Python/2", "Python");
my $biopythonMod = getLmod("Biopython", "Biopython");
my $usearchMod = getLmod("USEARCH/9", "USEARCH");
my $muscleMod = getLmod("MUSCLE/3", "MUSCLE");
my $blastMod = getLmod("BLAST+", "BLAST+");
my $cdhitMod = getLmod("CD-HIT", "CD-HIT");
my $blastDbPath = "$ENV{EFIDBPATH}/combined.fasta";
my $dbModule = $ENV{EFIDBMOD};
my $sbModule = $ENV{SHORTBRED_MOD};
my $sbParseSsnApp = $ENV{SHORTBRED_PARSESSN};
my $sbIdentifyApp = $ENV{SHORTBRED_IDENTIFY};

my $efiSbDir = $ENV{EFI_SHORTBRED_HOME};
my $ssnAccessionFile = "$outputDir/accession";
my $ssnClusterFile = "$outputDir/cluster";
my $fastaFile = "$outputDir/sequences.fa";
my $sbOutputDir = "$outputDir/id-temp";
my $sbMarkerFile = "$outputDir/markers.faa";
my $ssnMarker = "$outputDir/$outputSsnName";
my $jobPrefix = (defined $jobId and $jobId) ? "${jobId}_" : "";
my $submitName = "";
my $submitFile = "";
my $depId = 0;



#######################################################################################################################
# Get the clusters and accessions

my $B = $S->getBuilder();
$submitName = "sb_get_clusters";
$B->resource(1, 1, "5gb");
$B->addAction("module load $sbModule");
$B->addAction("$efiSbDir/get_clusters.pl -ssn $inputSsn -accession-file $ssnAccessionFile -cluster-file $ssnClusterFile");
$depId = doSubmit();


#######################################################################################################################
# Get the FASTA files from the database

$B = $S->getBuilder();
$submitName = "sb_get_fasta";
$B->resource(1, 1, "5gb");
$B->addAction("module load $sbModule");
$B->addAction("module load $dbModule");
$B->addAction("$efiSbDir/get_fasta.pl -id-file $ssnAccessionFile -output $fastaFile -blast-db $blastDbPath");
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
$B->resource(1, $np, "50gb");
$B->addAction("module load $sbModule");
#$B->addAction("module load $pythonMod");
#$B->addAction("module load $biopythonMod");
#$B->addAction("module load $usearchMod");
#$B->addAction("module load $muscleMod");
#$B->addAction("module load $blastMod");
#$B->addAction("module load $cdhitMod");
$B->addAction("python $sbIdentifyApp --threads $np --goi $fastaFile --refdb $blastDbPath --markers $sbMarkerFile --tmp $sbOutputDir");
$depId = doSubmit($depId);


#######################################################################################################################
# Build the XGMML file with the marker attributes added

$B = $S->getBuilder();
$submitName = "sb_make_xgmml";
$B->resource(1, 1, "10gb");
$B->addAction("module load $sbModule");
$B->addAction("$efiSbDir/make_ssn.pl -ssn-in $inputSsn -ssn-out $ssnMarker -marker-file $sbMarkerFile");
$B->addAction("touch $outputDir/job.completed");
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
    my $baseOutputDir = $ENV{PWD};
    my $outputDir = "$baseOutputDir/$outputDirName";
    mkdir $outputDir if not $dryRun;

    my $logDir = "$baseOutputDir/log";
    mkdir $logDir if not $dryRun;
    $logDir = "" if not -d $logDir;

    my $scriptDir = "$baseOutputDir/scripts";
    mkdir $scriptDir if not $dryRun;
    $scriptDir = $outputDir if not -d $scriptDir;

    return ($outputDir, $scriptDir, $logDir);
}



