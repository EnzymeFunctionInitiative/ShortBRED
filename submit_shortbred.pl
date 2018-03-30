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


my ($ssn, $outputDirName);
my ($np, $queue, $scheduler, $dryRun, $jobId);


my $result = GetOptions(

    "ssn=s"             => \$ssn,
    "tmpdir=s"          => \$outputDirName,


    "job-id=i"          => \$jobId,
    "np=i"              => \$np,
    "queue=s"           => \$queue,
    "scheduler=s"       => \$scheduler,
    "dryrun"            => \$dryRun,
);

my $usage =
"$0 -ssn=ssn_file -tmpdir=relative_output_dir -queue=cluster_queue [-job-id=job_id -np=num_cpu
    -scheduler=scheduler_type -dryrun]

    -ssn        path to input SSN file (relative or absolute)
    -tmpdir     name of output directory (relative)
    -job-id     number of the job [optional]
    -np         number of CPUs to use [optional, defaults to 24]
    -queue      the name of the cluster queue to submit to
    -scheduler  the type of scheduler to use (torque or slurm) [optional, defaults to slurm]
    -dryrun     only output the commands to be executed, do not submit to cluster [optional]
";

die $usage if not defined $ssn or not -f $ssn or not defined $outputDirName or not $outputDirName or
              not defined $queue or not $queue;


$np = 24 if not defined $np or not $np;
$scheduler = "slurm" if not defined $scheduler or not $scheduler;
$dryRun = 0 if not defined $dryRun;

my $sbModule = $ENV{SHORTBRED_MOD};
my $sbParseSsnApp = $ENV{SHORTBRED_PARSESSN};
my $efiSbDir = $ENV{EFI_SHORTBRED_HOME};

my ($outputDir, $scriptDir, $logDir) = initDirectoryStructure();
my ($S, $jobNamePrefix) = initScheduler($logDir);


my $ssnAccessionFile = "$outputDir/accession";
my $ssnClusterFile = "$outputDir/cluster";
my $jobPrefix = (defined $jobId and $jobId) ? "${jobId}_" : "";
my $submitName = "";
my $submitFile = "";



my $B = $S->getBuilder();
$submitName = "get_clusters";
$submitFile = "$scriptDir/$submitName.sh";
$B->addAction("module load $sbModule");
$B->addAction("$efiSbDir/get_clusters.pl -ssn $ssn -accession-file $ssnAccessionFile -cluster-file $ssnClusterFile");
$B->jobName("$jobPrefix$submitName");
$B->renderToFile($submitFile);
$S->submit($submitFile);













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



