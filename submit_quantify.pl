#!/usr/bin/env perl

use strict;
use warnings;

BEGIN {
    die "Please load efishared before runing this script" if not $ENV{EFISHARED};
    use lib $ENV{EFISHARED};
}

use Getopt::Long qw(:config pass_through);
use FindBin;
use EFI::SchedulerApi;
use EFI::Util qw(usesSlurm getLmod);

use lib $FindBin::Bin . "/lib";
use ShortBRED qw(getMetagenomeInfo);


die "Load efidb module in order to run this program." if not exists $ENV{EFIDBPATH};


my ($dbFiles, $metagenomeIdList, $idResDirName, $qResDirName, $jobDir);
my ($np, $queue, $memQueue, $scheduler, $dryRun, $jobId, $ssnName, $inputSsn);
my ($proteinFileName, $clusterFileName);
my ($parentIdentifyId, $parentQuantifyId);
my ($mergeAllRuns);
my ($searchType);


my $result = GetOptions(
    "metagenome-db=s"           => \$dbFiles,
    "metagenome-ids=s"          => \$metagenomeIdList,
    
    "results-dir-name=s"        => \$idResDirName,
    "quantify-dir=s"            => \$qResDirName,
    "job-dir=s"                 => \$jobDir,
    "ssn-in=s"                  => \$inputSsn, # full path
    "ssn-out-name=s"            => \$ssnName,
    "protein-file=s"            => \$proteinFileName,
    "cluster-file=s"            => \$clusterFileName,

    "parent-identify-id=i"      => \$parentIdentifyId,,
    "parent-quantify-id=i"      => \$parentQuantifyId,

    "search-type=s"             => \$searchType,
    "global-merge"              => \$mergeAllRuns,
    "job-id=i"                  => \$jobId,
    "np=i"                      => \$np,
    "queue=s"                   => \$queue,
    "mem-queue=s"               => \$memQueue,
    "scheduler=s"               => \$scheduler,
    "dry-run"                   => \$dryRun,
);

my $usage =
"$0 -metagenome-db=path_to_dbs -metagenome-ids=metagenome_ids -tmpdir=relative_output_dir
    -queue=cluster_queue [-job-id=job_id -np=num_cpu -scheduler=scheduler_type -dryrun]

    -metagenome-db              path(s) to metagenome database directories
    -metagenome-ids             comma-separated list of metagenome IDs to use in analysis

    -id-dir                     name of output directory containing identify results (relative)
    -quantify-dir               name of output directory containing quantify results (relative to id-dir)
    -ssn-in                     path to input SSN (add the marker and abundance data to it)
    -ssn-out-name               name of output SSN (containing marker and protein and cluster info)

    -parent-identify-id         ID of parent IDENTIFY job; if specified the data from the parent job will be used
                                rather than re-running the time-consuming quantify process
    -parent-quantify-id         ID of parent IDENTIFY job; if specified the data from the parent job will be used
                                rather than re-running the time-consuming quantify process

    -job-id                     number of the job [optional]
    -np                         number of CPUs to use [optional, defaults to 24]
    -queue                      the name of the cluster queue to submit to
    -scheduler                  the type of scheduler to use (torque or slurm) [optional, defaults to slurm]
    -dryrun                     only output the commands to be executed, do not submit to cluster [optional]
    
    -protein-file               name of output file containing protein abundances (median)
    -cluster-file               name of output file containing cluster abundances (median)
    -protein-norm               name of output file containing normalized protein abundances (median)
    -cluster-norm               name of output file containing normalized cluster abundances (median)
    -protein-genome-norm        name of output file containing AGS normalized protein abundances (median)
    -cluster-genome-norm        name of output file containing AGS normalized cluster abundances (median)
  
 If these are not specified, then they are created with the same as the median options above, but with .mean appended to the filename.
    -protein-file-mean          name of output file containing protein abundances (mean)
    -cluster-file-mean          name of output file containing cluster abundances (mean)
    -protein-norm-mean          name of output file containing normalized protein abundances (mean)
    -cluster-norm-mean          name of output file containing normalized cluster abundances (mean)
    -protein-genome-norm-mean   name of output file containing AGS normalized protein abundances (mean)
    -cluster-genome-norm-mean   name of output file containing AGS normalized cluster abundances (mean)
";

die $usage if not $dbFiles or not $metagenomeIdList
              or not $idResDirName or not $qResDirName
              or not $queue or not $inputSsn;

$memQueue = $queue                                      if not defined $memQueue or not $memQueue;
$np = 24                                                if not defined $np or not $np;
$scheduler = "slurm"                                    if not defined $scheduler or not $scheduler;
$dryRun = 0                                             if not defined $dryRun;
$clusterFileName = "cluster_abundance"                  if not defined $clusterFileName or not $clusterFileName;
$proteinFileName = "protein_abundance"                  if not defined $proteinFileName or not $proteinFileName;
$parentQuantifyId = 0                                   if not defined $parentQuantifyId;
$parentIdentifyId = 0                                   if not defined $parentIdentifyId;
$mergeAllRuns = 0                                       if not defined $mergeAllRuns;

$searchType = (defined $searchType and $searchType eq "diamond") ? "diamond" : "usearch";

my ($inputDir, $parentInputDir, $outputDir, $parentOutputDir, $scriptDir, $logDir) = initDirectoryStructure($parentIdentifyId, $parentQuantifyId, $jobDir);
my ($S, $jobNamePrefix) = initScheduler($logDir);


my $efiSbDir = $ENV{EFI_SHORTBRED_HOME};
my $sbModule = $ENV{SHORTBRED_MOD};
my $sbQuantifyApp = $ENV{SHORTBRED_QUANTIFY};
my $localMergeApp = "$efiSbDir/merge_shortbred.py";
my $agsFilePath = exists $ENV{SHORTBRED_AGS_FILE} ? $ENV{SHORTBRED_AGS_FILE} : "";
my $extraPerl = $ENV{EFI_PERL_ENV} // "";

my $clusterFileNameMedian = "${clusterFileName}.txt";
my $proteinFileNameMedian = "${proteinFileName}.txt";
my $clusterFileNameNormalizedMedian = "${clusterFileName}_normalized.txt";
my $proteinFileNameNormalizedMedian = "${proteinFileName}_normalized.txt";
my $clusterFileNameGenomeNormalizedMedian = "${clusterFileName}_genome_normalized.txt";
my $proteinFileNameGenomeNormalizedMedian = "${proteinFileName}_genome_normalized.txt";
my $clusterFileNameMean = "${clusterFileName}_mean.txt";
my $proteinFileNameMean = "${proteinFileName}_mean.txt";
my $clusterFileNameNormalizedMean = "${clusterFileName}_normalized_mean.txt";
my $proteinFileNameNormalizedMean = "${proteinFileName}_normalized_mean.txt";
my $clusterFileNameGenomeNormalizedMean = "${clusterFileName}_genome_normalized_mean.txt";
my $proteinFileNameGenomeNormalizedMean = "${proteinFileName}_genome_normalized_mean.txt";


my $ssnClusterFile = "$inputDir/cluster";
my $sbOutputDir = "$outputDir/quantify-temp";
my $sbMarkerFile = $parentIdentifyId ? "$parentInputDir/markers.faa" : "$inputDir/markers.faa";
my $cdhitFile = "$inputDir/cdhit.txt";

my $targetDir = $mergeAllRuns ? $inputDir : $outputDir; # If we merge all runs, the target dir is the input dir

my $clusterResultMedian = "$outputDir/$clusterFileNameMedian";
my $proteinResultMedian = "$outputDir/$proteinFileNameMedian";
my $clusterResultNormalizedMedian = "$outputDir/$clusterFileNameNormalizedMedian";
my $proteinResultNormalizedMedian = "$outputDir/$proteinFileNameNormalizedMedian";
my $clusterResultGenomeNormalizedMedian = "$outputDir/$clusterFileNameGenomeNormalizedMedian";
my $proteinResultGenomeNormalizedMedian = "$outputDir/$proteinFileNameGenomeNormalizedMedian";

my $clusterResultMean = "$outputDir/$clusterFileNameMean";
my $proteinResultMean = "$outputDir/$proteinFileNameMean";
my $clusterResultNormalizedMean = "$outputDir/$clusterFileNameNormalizedMean";
my $proteinResultNormalizedMean = "$outputDir/$proteinFileNameNormalizedMean";
my $clusterResultGenomeNormalizedMean = "$outputDir/$clusterFileNameGenomeNormalizedMean";
my $proteinResultGenomeNormalizedMean = "$outputDir/$proteinFileNameGenomeNormalizedMean";

my $clusterMerged = "$targetDir/$clusterFileNameMedian";
my $proteinMerged = "$targetDir/$proteinFileNameMedian";
my $clusterMergedNormalized = "$targetDir/$clusterFileNameNormalizedMedian";
my $proteinMergedNormalized = "$targetDir/$proteinFileNameNormalizedMedian";
my $clusterMergedGenomeNormalized = "$targetDir/$clusterFileNameGenomeNormalizedMedian";
my $proteinMergedGenomeNormalized = "$targetDir/$proteinFileNameGenomeNormalizedMedian";

my $outputSsn = "$targetDir/$ssnName.xgmml"; # Merge all results into one SSN, in the results/ directory
my $metadataFile = "$targetDir/metadata.tab";
my $jobPrefix = (defined $jobId and $jobId) ? "${jobId}_" : "";
my $submitName = "";
my $submitFile = "";
my $depId = 0;


my @metagenomeIds = split(m/,/, $metagenomeIdList);
my ($metagenomeInfo, $mgMetadata) = getMetagenomeInfo($dbFiles);

# Get the specific avg genome size file, if present.
my ($mgDir) = split(m/,/, $dbFiles);
my $agsFileName = "AvgGenomeSize.txt";
if ($mgDir and -f "$mgDir/$agsFileName") {
    $agsFilePath = "$mgDir/$agsFileName";
}



# Get the list of result files.
my %resFilesMedian;
my %resFilesMean;
my %mgFiles;
foreach my $mgId (@metagenomeIds) {
    if (not exists $metagenomeInfo->{$mgId}) {
        warn "Metagenome file does not exists for $mgId.";
        next;
    }
    my $mgFile = $metagenomeInfo->{$mgId}->{file_path};
    my $resFileMedian = $parentQuantifyId ? "$parentOutputDir/$mgId.txt" : "$outputDir/$mgId.txt";
    my $resFileMean = "$resFileMedian.mean";
    $mgFiles{$mgId} = $mgFile;
    $resFilesMedian{$mgId} = $resFileMedian;
    $resFilesMean{$mgId} = $resFileMean;
}



#######################################################################################################################
# Run ShortBRED-Quantify on the markers
#
# Don't run this if we are creating SSNs from another job's results.

my $searchTypeArgs = $searchType ? "--search_program $searchType" : "";
    
my $B = $S->getBuilder();
my $useTasks = 1;
if (not $parentQuantifyId) {
    $submitName = "sb_quantify";
    $B->addAction("module load $sbModule");
    if (not $useTasks) {
        $B->resource(1, $np, "300gb");
        foreach my $mgId (@metagenomeIds) {
            my $mgFile = $mgFiles{$mgId};
            my $resFileMedian = $resFilesMedian{$mgId};
            my $resFileMean = $resFilesMean{$mgId};
            $B->addAction("python $sbQuantifyApp --threads $np $searchTypeArgs --markers $sbMarkerFile --wgs $mgFile --results $resFileMedian --results-mean $resFileMean --tmp $sbOutputDir-$mgId");
        }
    } else {
        my $numMg = scalar @metagenomeIds;
        my $numFiles = $numMg < $np ? 1 : int($numMg / $np + 1);
        my $maxTask = $numMg >= $np ? $np : $numMg;
        if ($maxTask == $np) {
            my $excessTask = int(($numFiles * $np - $numMg) / $numFiles);
            $maxTask = $maxTask - $excessTask;
        }

        #my $numIds = scalar @metagenomeIds;
        #my @sizeOrderedIds = @metagenomeIds; #sort { -s $mgFiles{$a} cmp -s $mgFiles{$b} } @metagenomeIds;
        #my @evenDistIds; # Due to the allocation algorithm we're using here, some of the entries in the array will be empty.
        #my $evenSize = $maxTask * $numFiles;
        #my %sizes;
        #for (my $i = 0; $i < $evenSize; $i++) {
        #    my $taskNum = ($i % $maxTask) * $numFiles;
        #    my $idx = int($i / $maxTask) + $taskNum;
        #    my $mgId = $idx < $numIds ? $sizeOrderedIds[$idx] : "";
        #    $evenDistIds[$idx] = $mgId;
        #    print "$i $taskNum $idx $mgId\n";
        #    
        #    $sizes{$taskNum} = 0 if not exists $sizes{$taskNum};
        #    if ($mgId) {
        #        $sizes{$taskNum} += -s $mgFiles{$mgId};
        #    }
        #}

        my $tmpMarker = "$outputDir/markers.faa.{JOB_ARRAYID}";

        $B->jobArray("1-$maxTask");
        $B->resource(1, 1, "13gb");
        my $c = 0;
        foreach my $mgId (@metagenomeIds) {
            my $mgFile = $mgFiles{$mgId};
            my $resFileMedian = $resFilesMedian{$mgId};
            my $resFileMean = $resFilesMean{$mgId};
            if ($c % $numFiles == 0) {
                my $aid = int($c / $numFiles) + 1;
                if ($c > 0) {
                    $B->addAction("    rm $tmpMarker");
                    $B->addAction("fi");
                }
                $B->addAction("if [ {JOB_ARRAYID} == $aid ]; then");
                $B->addAction("    cp $sbMarkerFile $tmpMarker"); # Copy to possibly help performance out.
            }
            $B->addAction("    python $sbQuantifyApp $searchTypeArgs --markers $tmpMarker --wgs $mgFile --results $resFileMedian --results-mean $resFileMean --tmp $sbOutputDir-$mgId");
            $c++;
        }
        if ($c > 1 or $numFiles == 1) {
            $B->addAction("    rm $tmpMarker");
            $B->addAction("fi");
        }
    }
    $depId = doSubmit($depId);
}


#######################################################################################################################
# Merge quantify outputs into one table.

# Sort the metagenome IDs according to a body site
my @sortedMgIds = sort mgSortFn @metagenomeIds;

my $resFileMedianList = join(" ", map { $resFilesMedian{$_} } @sortedMgIds);
my $resFileMeanList = join(" ", map { $resFilesMean{$_} } @sortedMgIds);

$B = $S->getBuilder();
$submitName = "sb_merge_quantify";
$B->resource(1, 1, "20gb");
$B->addAction("python $localMergeApp $resFileMedianList -C $clusterResultMedian -p $proteinResultMedian -c $ssnClusterFile");
$B->addAction("python $localMergeApp $resFileMedianList -C $clusterResultNormalizedMedian -p $proteinResultNormalizedMedian -c $ssnClusterFile -n");
if ($agsFilePath) {
    $B->addAction("python $localMergeApp $resFileMedianList -C $clusterResultGenomeNormalizedMedian -p $proteinResultGenomeNormalizedMedian -c $ssnClusterFile -g $agsFilePath");
}
$B->addAction("python $localMergeApp $resFileMeanList -C $clusterResultMean -p $proteinResultMean -c $ssnClusterFile");
$B->addAction("python $localMergeApp $resFileMeanList -C $clusterResultNormalizedMean -p $proteinResultNormalizedMean -c $ssnClusterFile -n");
if ($agsFilePath) {
    $B->addAction("python $localMergeApp $resFileMeanList -C $clusterResultGenomeNormalizedMean -p $proteinResultGenomeNormalizedMean -c $ssnClusterFile -g $agsFilePath");
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
$B->addAction("source $extraPerl");
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
    $B->addAction("$efiSbDir/merge_data.pl $mergeArgsShared -merged-protein $proteinMerged -merged-cluster $clusterMerged -protein-name $proteinFileNameMedian -cluster-name $clusterFileNameMedian");
    $B->addAction("$efiSbDir/merge_data.pl $mergeArgsShared -merged-protein $proteinMergedNormalized -merged-cluster $clusterMergedNormalized -protein-name $proteinFileNameNormalizedMedian -cluster-name $clusterFileNameNormalizedMedian");
    if ($agsFilePath) {
        $B->addAction("$efiSbDir/merge_data.pl $mergeArgsShared -merged-protein $proteinMergedGenomeNormalized -merged-cluster $clusterMergedGenomeNormalized -protein-name $proteinFileNameGenomeNormalizedMedian -cluster-name $clusterFileNameGenomeNormalizedMedian");
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
$B->addAction("zip -j $outputDir/quantify_ssn.zip $outputSsn");
$B->addAction("$efiSbDir/lock_merged_output.pl $lockFile unlock")           if $mergeAllRuns;
$B->addAction("$efiSbDir/create_quantify_metadata.pl -protein-abundance $proteinMerged -metadata $metadataFile");
$B->addAction("rm -rf $outputDir/quantify-temp-*");
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
    my $jobDir = shift || "";

    my $baseOutputDir = ($jobDir and -d $jobDir) ? $jobDir : $ENV{PWD};
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



sub mgSortFn {
    my $ag = $metagenomeInfo->{$a}->{gender};
    my $bg = $metagenomeInfo->{$b}->{gender};
    my $ab = $metagenomeInfo->{$a}->{bodysite};
    my $bb = $metagenomeInfo->{$b}->{bodysite};
    my $ao = exists $mgMetadata->{$ab} ? $mgMetadata->{$ab}->{order} : 0;
    my $bo = exists $mgMetadata->{$bb} ? $mgMetadata->{$bb}->{order} : 0;

    # Compare by order
    my $ocmp = $ao cmp $bo;
    return $ocmp if $ocmp;

    # Compare by body site
    my $bscmp = $ab cmp $bb;
    return $bscmp if $bscmp;

    # Compare by gender
    my $gcmp = $ag cmp $bg;
    return $gcmp if $gcmp;

    return $a cmp $b;
}

