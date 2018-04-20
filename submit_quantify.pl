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


my ($dbFiles, $metagenomeIdList, $idResDirName, $qResDirName);
my ($np, $queue, $scheduler, $dryRun, $jobId, $ssnName, $inputSsn, $proteinFileName, $clusterFileName);


my $result = GetOptions(
    "metagenome-db=s"   => \$dbFiles,
    "metagenome-ids=s"  => \$metagenomeIdList,
    
    "id-dir=s"          => \$idResDirName,
    "quantify-dir=s"    => \$qResDirName,
    "ssn-in=s"          => \$inputSsn,
    "ssn-out-name=s"    => \$ssnName,
    "protein-file=s"    => \$proteinFileName,
    "cluster-file=s"    => \$clusterFileName,

    "job-id=i"          => \$jobId,
    "np=i"              => \$np,
    "queue=s"           => \$queue,
    "scheduler=s"       => \$scheduler,
    "dryrun"            => \$dryRun,
);

my $usage =
"$0 -metagenome-db=path_to_dbs -metagenome-ids=metagenome_ids -tmpdir=relative_output_dir
    -queue=cluster_queue [-job-id=job_id -np=num_cpu -scheduler=scheduler_type -dryrun]

    -metagenome-db      path(s) to metagenome database directories
    -metagenome-ids     comma-separated list of metagenome IDs to use in analysis

    -id-dir             name of output directory containing identify results (relative)
    -quantify-dir       name of output directory containing quantify results (relative to id-dir)
    -ssn-in             path to input SSN (add the marker and abundance data to it)
    -ssn-out-name       name of output SSN (containing marker and protein and cluster info)
    -protein-file       name of output file containing protein abundances
    -cluster-file       name of output file containing cluster abundances
    
    -job-id             number of the job [optional]
    -np                 number of CPUs to use [optional, defaults to 24]
    -queue              the name of the cluster queue to submit to
    -scheduler          the type of scheduler to use (torque or slurm) [optional, defaults to slurm]
    -dryrun             only output the commands to be executed, do not submit to cluster [optional]
";

die $usage if not defined $dbFiles or not $dbFiles or not defined $metagenomeIdList or not $metagenomeIdList
              or not defined $idResDirName or not $idResDirName or not defined $qResDirName
              or not $qResDirName or not defined $queue or not $queue;

$np = 24 if not defined $np or not $np;
$scheduler = "slurm" if not defined $scheduler or not $scheduler;
$dryRun = 0 if not defined $dryRun;

my ($inputDir, $outputDir, $scriptDir, $logDir) = initDirectoryStructure();
my ($S, $jobNamePrefix) = initScheduler($logDir);

my $efiSbDir = $ENV{EFI_SHORTBRED_HOME};
my $sbModule = $ENV{SHORTBRED_MOD};
my $sbQuantifyApp = $ENV{SHORTBRED_QUANTIFY};
my $localMergeApp = "$efiSbDir/merge_shortbred.py";

my $ssnClusterFile = "$inputDir/cluster";
my $sbOutputDir = "$outputDir/quantify-temp";
my $sbMarkerFile = "$inputDir/markers.faa";
my $clusterResult = "$outputDir/$clusterFileName";
my $proteinResult = "$outputDir/$proteinFileName";
my $clusterMerged = "$inputDir/$clusterFileName";
my $proteinMerged = "$inputDir/$proteinFileName";
my $outputSsn = "$inputDir/$ssnName"; # Merge all results into one SSN, in the results/ directory
my $jobPrefix = (defined $jobId and $jobId) ? "${jobId}_" : "";
my $submitName = "";
my $submitFile = "";
my $depId = 0;


my @metagenomeIds = split(m/,/, $metagenomeIdList);
my $metagenomeInfo = getMetagenomeInfo($dbFiles, @metagenomeIds);




#######################################################################################################################
# Run ShortBRED-Quantify on the markers

my @resFiles;

my $B = $S->getBuilder();
$submitName = "sb_quantify";
$B->resource(1, $np, "50gb");
$B->addAction("module load $sbModule");
#$B->addAction("module load $pythonMod");
#$B->addAction("module load $biopythonMod");
#$B->addAction("module load $usearchMod");
#$B->addAction("module load $muscleMod");
#$B->addAction("module load $blastMod");
#$B->addAction("module load $cdhitMod");
foreach my $mgId (@metagenomeIds) {
    if (not exists $metagenomeInfo->{$mgId}) {
        warn "Metagenome file does not exists for $mgId.";
        next;
    }
    my $mgFile = $metagenomeInfo->{$mgId}->{file_path};
    my $resFile = "$outputDir/$mgId.txt";
    push(@resFiles, $resFile);
    $B->addAction("python $sbQuantifyApp --threads $np --markers $sbMarkerFile --wgs $mgFile --results $resFile --tmp $sbOutputDir-$mgId");
}
$depId = doSubmit($depId);


#######################################################################################################################
# Merge quantify outputs into one table.

my $resFileList = join(" ", @resFiles);

$B = $S->getBuilder();
$submitName = "sb_merge_quantify";
$B->resource(1, $np, "5gb");
$B->addAction("python $localMergeApp $resFileList -C $clusterResult -p $proteinResult -c $ssnClusterFile");
$depId = doSubmit($depId);


#######################################################################################################################
# Build the XGMML file with the marker attributes and the abundances added added

my $lockFile = "$inputDir/merge.lock";

$B = $S->getBuilder();
$B->setScriptAbortOnError(0); # Disable abort on error, so that we can disable the merged output lock.
$submitName = "sb_q_make_xgmml";
$B->resource(1, 1, "10gb");
$B->addAction("module load $sbModule");
$B->addAction("$efiSbDir/lock_merged_output.pl $lockFile lock");
$B->addAction("OUT=\$?");
$B->addAction("if [ \$OUT -ne 0 ]; then");
$B->addAction("    echo \"unable to lock output\"");
$B->addAction("    echo \$OUT > $outputDir/ssn.failed");
$B->addAction("    rm -f $lockFile");
$B->addAction("    exit 620");
$B->addAction("fi");
$B->addAction("$efiSbDir/merge_data.pl -input-dir $inputDir -quantify-dir-pat quantify- -merged-protein $proteinMerged -merged-cluster $clusterMerged -protein-name $proteinFileName -cluster-name $clusterFileName");
$B->addAction("$efiSbDir/make_ssn.pl -ssn-in $inputSsn -ssn-out $outputSsn -protein-file $proteinMerged -cluster-file $clusterMerged -quantify");
$B->addAction("OUT=\$?");
$B->addAction("if [ \$OUT -ne 0 ]; then");
$B->addAction("    echo \"make SSN failed.\"");
$B->addAction("    echo \$OUT > $outputDir/ssn.failed");
$B->addAction("    rm -f $lockFile");
$B->addAction("    exit 621");
$B->addAction("fi");
$B->addAction("$efiSbDir/lock_merged_output.pl $lockFile unlock");
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
    my $idOutputDir = "$baseOutputDir/$idResDirName";
    my $qOutputDir = "$idOutputDir/$qResDirName";

    mkdir $qOutputDir if not -d $qOutputDir;

    my $logDir = "$baseOutputDir/log";
    $logDir = "" if not -d $logDir;

    my $scriptDir = $qOutputDir;
    #my $scriptDir = "$baseOutputDir/scripts";
    #$scriptDir = $qOutputDir if not -d $scriptDir;

    return ($idOutputDir, $qOutputDir, $scriptDir, $logDir);
}


sub getMetagenomeInfo {
    my $dbList = shift;
    my @ids = @_;

    my $data = {};

    my @dbs = split(m/,/, $dbList);
    foreach my $dbFile (@dbs) {
        (my $dbDir = $dbFile) =~ s%^(.*)/[^/]+$%$1%;
        open DB, $dbFile or next;
        while (<DB>) {
            next if m/^#/;
            chomp;
            my ($id, $name, $desc, $file) = split(m/\t/);
            $file = "$dbDir/$id/$file" if $file;
            $data->{$id} = {name => $name, desc => $desc, file_path => $file};
        }
        close DB;
    }

    return $data;
}


