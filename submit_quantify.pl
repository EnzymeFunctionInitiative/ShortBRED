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
use ShortBRED qw(getMetagenomeInfo);


die "Load efidb module in order to run this program." if not exists $ENV{EFIDBPATH};


my ($dbFiles, $metagenomeIdList, $idResDirName, $qResDirName);
my ($np, $queue, $memQueue, $scheduler, $dryRun, $jobId, $ssnName, $inputSsn, $proteinFileName, $clusterFileName);
my ($proteinFileNameNormalized, $clusterFileNameNormalized);
my ($proteinFileNameGenomeNormalized, $clusterFileNameGenomeNormalized);
my ($parentIdentifyId, $parentQuantifyId);
my ($mergeAllRuns);


my $result = GetOptions(
    "metagenome-db=s"           => \$dbFiles,
    "metagenome-ids=s"          => \$metagenomeIdList,
    
    "id-dir=s"                  => \$idResDirName,
    "quantify-dir=s"            => \$qResDirName,
    "ssn-in=s"                  => \$inputSsn,
    "ssn-out-name=s"            => \$ssnName,
    "protein-file=s"            => \$proteinFileName,
    "cluster-file=s"            => \$clusterFileName,
    "protein-norm=s"            => \$proteinFileNameNormalized,
    "cluster-norm=s"            => \$clusterFileNameNormalized,
    "protein-genome-norm=s"     => \$proteinFileNameGenomeNormalized,
    "cluster-genome-norm=s"     => \$clusterFileNameGenomeNormalized,

    "parent-identify-id=i"      => \$parentIdentifyId,,
    "parent-quantify-id=i"      => \$parentQuantifyId,

    "global-merge"              => \$mergeAllRuns,
    "job-id=i"                  => \$jobId,
    "np=i"                      => \$np,
    "queue=s"                   => \$queue,
    "mem-queue=s"               => \$memQueue,
    "scheduler=s"               => \$scheduler,
    "dryrun"                    => \$dryRun,
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
    -protein-norm       name of output file containing normalized protein abundances
    -cluster-norm       name of output file containing normalized cluster abundances

    -parent-job-id      ID of parent IDENTIFY job; if specified the data from the parent job will be used
                        rather than re-running the time-consuming quantify process
    
    -job-id             number of the job [optional]
    -np                 number of CPUs to use [optional, defaults to 24]
    -queue              the name of the cluster queue to submit to
    -scheduler          the type of scheduler to use (torque or slurm) [optional, defaults to slurm]
    -dryrun             only output the commands to be executed, do not submit to cluster [optional]
";

die $usage if not defined $dbFiles or not $dbFiles or not defined $metagenomeIdList or not $metagenomeIdList
              or not defined $idResDirName or not $idResDirName or not defined $qResDirName
              or not $qResDirName or not defined $queue or not $queue;

$memQueue = $queue                              if not defined $memQueue or not $memQueue;
$np = 24                                        if not defined $np or not $np;
$scheduler = "slurm"                            if not defined $scheduler or not $scheduler;
$dryRun = 0                                     if not defined $dryRun;
$clusterFileName = "cluster_abundance.txt"      if not defined $clusterFileName or not $clusterFileName;
$proteinFileName = "protein_abundnace.txt"      if not defined $proteinFileName or not $proteinFileName;
$parentQuantifyId = 0                           if not defined $parentQuantifyId;
$parentIdentifyId = 0                           if not defined $parentIdentifyId;
$mergeAllRuns = 0                               if not defined $mergeAllRuns;

my ($inputDir, $parentInputDir, $outputDir, $parentOutputDir, $scriptDir, $logDir) = initDirectoryStructure($parentIdentifyId, $parentQuantifyId);
my ($S, $jobNamePrefix) = initScheduler($logDir);


my $efiSbDir = $ENV{EFI_SHORTBRED_HOME};
my $sbModule = $ENV{SHORTBRED_MOD};
my $sbQuantifyApp = $ENV{SHORTBRED_QUANTIFY};
my $localMergeApp = "$efiSbDir/merge_shortbred.py";
my $agsFilePath = exists $ENV{SHORTBRED_AGS_FILE} ? $ENV{SHORTBRED_AGS_FILE} : "";

($clusterFileNameNormalized = $clusterFileName) =~ s/\.txt$/_normalized.txt/ if not defined $clusterFileNameNormalized or not $clusterFileNameNormalized;
($proteinFileNameNormalized = $proteinFileName) =~ s/\.txt$/_normalized.txt/ if not defined $proteinFileNameNormalized or not $proteinFileNameNormalized;
($clusterFileNameGenomeNormalized = $clusterFileName) =~ s/\.txt$/_genome_normalized.txt/ if not defined $clusterFileNameGenomeNormalized or not $clusterFileNameGenomeNormalized;
($proteinFileNameGenomeNormalized = $proteinFileName) =~ s/\.txt$/_genome_normalized.txt/ if not defined $proteinFileNameGenomeNormalized or not $proteinFileNameGenomeNormalized;


my $ssnClusterFile = "$inputDir/cluster";
my $sbOutputDir = "$outputDir/quantify-temp";
my $sbMarkerFile = $parentIdentifyId ? "$parentInputDir/markers.faa" : "$inputDir/markers.faa";
my $cdhitFile = "$inputDir/cdhit.txt";

my $targetDir = $mergeAllRuns ? $inputDir : $outputDir; # If we merge all runs, the target dir is the input dir
my $clusterResult = "$outputDir/$clusterFileName";
my $proteinResult = "$outputDir/$proteinFileName";
my $clusterResultNormalized = "$outputDir/$clusterFileNameNormalized";
my $proteinResultNormalized = "$outputDir/$proteinFileNameNormalized";
my $clusterResultGenomeNormalized = "$outputDir/$clusterFileNameGenomeNormalized";
my $proteinResultGenomeNormalized = "$outputDir/$proteinFileNameGenomeNormalized";
my $clusterMerged = "$targetDir/$clusterFileName";
my $proteinMerged = "$targetDir/$proteinFileName";
my $clusterMergedNormalized = "$targetDir/$clusterFileNameNormalized";
my $proteinMergedNormalized = "$targetDir/$proteinFileNameNormalized";
my $clusterMergedGenomeNormalized = "$targetDir/$clusterFileNameGenomeNormalized";
my $proteinMergedGenomeNormalized = "$targetDir/$proteinFileNameGenomeNormalized";
my $outputSsn = "$targetDir/$ssnName"; # Merge all results into one SSN, in the results/ directory
my $jobPrefix = (defined $jobId and $jobId) ? "${jobId}_" : "";
my $submitName = "";
my $submitFile = "";
my $depId = 0;


my @metagenomeIds = split(m/,/, $metagenomeIdList);
my $metagenomeInfo = getMetagenomeInfo($dbFiles, @metagenomeIds);



# Get the list of result files.
my %resFiles;
my %mgFiles;
foreach my $mgId (@metagenomeIds) {
    if (not exists $metagenomeInfo->{$mgId}) {
        warn "Metagenome file does not exists for $mgId.";
        next;
    }
    my $mgFile = $metagenomeInfo->{$mgId}->{file_path};
    my $resFile = $parentQuantifyId ? "$parentOutputDir/$mgId.txt" : "$outputDir/$mgId.txt";
    $mgFiles{$mgId} = $mgFile;
    $resFiles{$mgId} = $resFile;
    #push(@resFiles, $resFile);
}


#######################################################################################################################
# Run ShortBRED-Quantify on the markers
#
# Don't run this if we are creating SSNs from another job's results.
my $B = $S->getBuilder();
my $useTasks = 1;
if (not $parentQuantifyId) {
    $submitName = "sb_quantify";
    $B->addAction("module load $sbModule");
    if (not $useTasks) {
        $B->resource(1, $np, "300gb");
        foreach my $mgId (@metagenomeIds) {
            my $mgFile = $mgFiles{$mgId};
            my $resFile = $resFiles{$mgId};
            $B->addAction("python $sbQuantifyApp --threads $np --markers $sbMarkerFile --wgs $mgFile --results $resFile --tmp $sbOutputDir-$mgId");
        }
    } else {
        my $numMg = scalar @metagenomeIds;
        my $numFiles = $numMg < 24 ? 1 : int($numMg / 24 + 1);
        my $maxTask = $numMg >= 24 ? 24 : $numMg;
        if ($maxTask == 24) {
            my $excessTask = int(($numFiles * 24 - $numMg) / $numFiles);
            $maxTask = $maxTask - $excessTask;
        }

        my $tmpMarker = "$outputDir/markers.faa.{JOB_ARRAYID}";

        $B->jobArray("1-$maxTask");
        $B->resource(1, 1, "13gb");
        my $c = 0;
        foreach my $mgId (@metagenomeIds) {
            my $mgFile = $mgFiles{$mgId};
            my $resFile = $resFiles{$mgId};
            if ($c % $numFiles == 0) {
                my $aid = int($c / $numFiles) + 1;
                if ($c > 0) {
                    $B->addAction("    rm $tmpMarker");
                    $B->addAction("fi");
                }
                $B->addAction("if [ {JOB_ARRAYID} == $aid ]; then");
                $B->addAction("    cp $sbMarkerFile $tmpMarker"); # Copy to possibly help performance out.
            }
            $B->addAction("    python $sbQuantifyApp --markers $tmpMarker --wgs $mgFile --results $resFile --tmp $sbOutputDir-$mgId");
            $c++;
        }
        if ($c > 1) {
            $B->addAction("    rm $tmpMarker");
            $B->addAction("fi");
        }
    }
    $depId = doSubmit($depId);
}


#######################################################################################################################
# Merge quantify outputs into one table.

my $resFileList = join(" ", sort values %resFiles);

$B = $S->getBuilder();
$submitName = "sb_merge_quantify";
$B->resource(1, $np, "20gb");
$B->addAction("python $localMergeApp $resFileList -C $clusterResult -p $proteinResult -c $ssnClusterFile");
$B->addAction("python $localMergeApp $resFileList -C $clusterResultNormalized -p $proteinResultNormalized -c $ssnClusterFile -n");
if ($agsFilePath) {
    $B->addAction("python $localMergeApp $resFileList -C $clusterResultGenomeNormalized -p $proteinResultGenomeNormalized -c $ssnClusterFile -g $agsFilePath");
}
$depId = doSubmit($depId);


#######################################################################################################################
# Build the XGMML file with the marker attributes and the abundances added added
my $lockFile = "$inputDir/merge.lock";

$B = $S->getBuilder();
$B->setScriptAbortOnError(0); # Disable abort on error, so that we can disable the merged output lock.
$submitName = "sb_q_make_xgmml";
$B->queue($memQueue);
$B->resource(1, 1, "200gb");
$B->addAction("MGIDS=\"$metagenomeIdList\"");
$B->addAction("MGDB=\"$dbFiles\"");
$B->addAction("module load $sbModule");

if ($mergeAllRuns) {
    $B->addAction("$efiSbDir/lock_merged_output.pl $lockFile lock");
    $B->addAction("OUT=\$?");
    $B->addAction("if [ \$OUT -ne 0 ]; then");
    $B->addAction("    echo \"unable to lock output\"");
    $B->addAction("    echo \$OUT > $outputDir/ssn.failed");
    $B->addAction("    rm -f $lockFile");
    $B->addAction("    exit 620");
    $B->addAction("fi");

    my $mergeArgsShared = "-cluster-list-file $ssnClusterFile -input-dir $inputDir -quantify-dir-pat quantify- -force-include $qResDirName";
    $B->addAction("$efiSbDir/merge_data.pl $mergeArgsShared -merged-protein $proteinMerged -merged-cluster $clusterMerged -protein-name $proteinFileName -cluster-name $clusterFileName");
    $B->addAction("$efiSbDir/merge_data.pl $mergeArgsShared -merged-protein $proteinMergedNormalized -merged-cluster $clusterMergedNormalized -protein-name $proteinFileNameNormalized -cluster-name $clusterFileNameNormalized");
    if ($agsFilePath) {
        $B->addAction("$efiSbDir/merge_data.pl $mergeArgsShared -merged-protein $proteinMergedGenomeNormalized -merged-cluster $clusterMergedGenomeNormalized -protein-name $proteinFileNameGenomeNormalized -cluster-name $clusterFileNameGenomeNormalized");
    }
}
$B->addAction("$efiSbDir/make_ssn.pl -ssn-in $inputSsn -ssn-out $outputSsn -protein-file $proteinMergedGenomeNormalized -cluster-file $clusterMergedGenomeNormalized -cdhit-file $cdhitFile -quantify -metagenome-db \$MGDB -metagenome-ids \$MGIDS");
$B->addAction("OUT=\$?");
$B->addAction("if [ \$OUT -ne 0 ]; then");
$B->addAction("    echo \"make SSN failed.\"");
$B->addAction("    echo \$OUT > $outputDir/ssn.failed");
$B->addAction("    rm -f $lockFile")                                        if $mergeAllRuns;
$B->addAction("    exit 621");
$B->addAction("fi");
$B->addAction("zip -j $outputSsn.zip $outputSsn");
$B->addAction("$efiSbDir/lock_merged_output.pl $lockFile unlock")           if $mergeAllRuns;
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
    my $parentIdentifyId = shift;
    my $parentQuantifyId = shift;

    my $baseOutputDir = $ENV{PWD};
    my $idOutputDir = "$baseOutputDir/$idResDirName";
    my $qOutputDir = "$idOutputDir/$qResDirName";

    mkdir $qOutputDir if not -d $qOutputDir;

    my $logDir = "$baseOutputDir/log";
    $logDir = "" if not -d $logDir;

    my $scriptDir = $qOutputDir;

    my $parentIdentifyOutputDir = "";
    if ($parentIdentifyId) {
        ($parentIdentifyOutputDir = $baseOutputDir) =~ s%/(\d+)/*$%/$parentIdentifyId%;
        $parentIdentifyOutputDir .= "/$idResDirName";
    }
    my $parentQuantifyOutputDir = "";
    if ($parentQuantifyId) {
        (my $qTemp = $qResDirName) =~ s/-\d+$/-$parentQuantifyId/;
        $parentQuantifyOutputDir = "$parentIdentifyOutputDir/$qTemp";
    }

    return ($idOutputDir, $parentIdentifyOutputDir, $qOutputDir, $parentQuantifyOutputDir, $scriptDir, $logDir);
    #return ($idOutputDir, $qOutputDir, $parentOutputDir, $scriptDir, $logDir);
}




