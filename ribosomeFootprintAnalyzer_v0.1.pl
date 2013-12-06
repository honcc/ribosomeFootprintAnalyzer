#!/usr/bin/perl -w

#====================================================================================================================================================#
#<use>
$|++; #---turn on the auto flush for the progress bar
use strict;
use File::Path;
use Time::HiRes qw( time );
use Storable;
use Getopt::Long;
use File::Basename;
use File::Spec::Functions qw(rel2abs);
use List::Util qw (sum shuffle min max);
use threads;
use threads::shared;
use Statistics::Descriptive;
use URI::Escape;
use Cwd 'abs_path';
#<\use>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<doc>
#	Description
#		This is a perl script to read the perl pileup storables generated from BAMToReadEndPerlStorable and to charaterize the distribution of read 5’end 
#	distribution along the codons.
#
#	Input
#		--rd5EndStorableIndexPath=		file path[compulsory]; path of the index file of the rd5End pileupStorable to be counted;
#		--gffPath=						file path[compulsory]; path of the reference GFF for gene annotation;
#		--fastaPath=					file path [compulsory]; the path fasta file contains the genome sequence, for generating blank perl storables;
#		--maxThread=					integer [4]; max number of threads to be used;
#		--outDir=						directory path ['./BAMToReadEndPerlStorable/']; output directory;
#
#	Usage
#		/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pileupAndCounting/ribosomeFootprintAnalyzer/v0.1/ribosomeFootprintAnalyzer_v0.1.pl --rd5EndStorableIndexPath=/Volumes/F_Analysis/NGS/results/nicolaiTrypRiboFootprint/N26_footprint/BAMToReadEndPerlStorable/countMode.5.offset.0.baseComp.no/cntgCovPls/index.hsh.pls --gffPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/trypanosome/inUse/927/TbruceiTreu927_TriTrypDB-4.2.gff --fastaPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/trypanosome/inUse/927/TbruceiTreu927Genomic_TriTrypDB-4.2.fasta
#
#	v0.1
#		[Wed 21 Aug 2013 20:43:41 CEST] debut;
#
#<\doc>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<lastCmdCalled>
#
#	[2013-09-28 18:24]	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pileupAndCounting/ribosomeFootprintAnalyzer/v0.1/ribosomeFootprintAnalyzer_v0.1.pl --gffPath=/Volumes/F_Analysis/NGS/results/nicolaiTrypRiboFootprint/NS026_27_29_30_pooled/RNASeqORFFinder_final/VC1.CN2.CP70.UL9.NL.30.KO_no.MG.20/GFF/novelORF.transcribed_yes.longest_yes.nonAnno.gff --fastaPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/trypanosome/inUse/927/v4.2/TbruceiTreu927Genomic_TriTrypDB-4.2.fasta --rd5EndStorableIndexPath=/Volumes/F_Analysis/NGS/results/nicolaiTrypRiboFootprint/NS025_28_pooled/BAMToReadEndPerlStorable/countMode.5.offset.0.baseComp.no/cntgCovPls/index.hsh.pls --outDir=/Volumes/F_Analysis/NGS/results/nicolaiTrypRiboFootprint/NS025_28_pooled/ribosomeFootprintAnalyzer/nonAnno/
#
#	/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pileupAndCounting/ribosomeFootprintAnalyzer/v0.1/ribosomeFootprintAnalyzer_v0.1.pl
#	--gffPath=/Volumes/F_Analysis/NGS/results/nicolaiTrypRiboFootprint/NS026_27_29_30_pooled/RNASeqORFFinder_final/VC1.CN2.CP70.UL9.NL.30.KO_no.MG.20/GFF/novelORF.transcribed_yes.longest_yes.nonAnno.gff
#	--fastaPath=/Volumes/A_MPro2TB/softwareForNGS/resources/genome/trypanosome/inUse/927/v4.2/TbruceiTreu927Genomic_TriTrypDB-4.2.fasta
#	--rd5EndStorableIndexPath=/Volumes/F_Analysis/NGS/results/nicolaiTrypRiboFootprint/NS025_28_pooled/BAMToReadEndPerlStorable/countMode.5.offset.0.baseComp.no/cntgCovPls/index.hsh.pls
#	--outDir=/Volumes/F_Analysis/NGS/results/nicolaiTrypRiboFootprint/NS025_28_pooled/ribosomeFootprintAnalyzer/nonAnno/
#
#<\lastCmdCalled>
#====================================================================================================================================================#

#====================================================================================================================================================#
#<global>
my $global_scriptDir = dirname(rel2abs($0));
open DEBUGLOG, ">", "$global_scriptDir/debug.log.txt";
#<\global>
#====================================================================================================================================================#

#====================================================================================================================================================#
{	#Main sections lexical scope starts
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 0_startingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|746, readParameters|932
#	secondaryDependOnSub: currentTime|277
#
#<section ID="startingTasks" num="0">
########################################################################## 
&printCMDLogOrFinishMessage("CMDLog");#->746

my ($rd5EndStorableIndexPath, $gffPath, $fastaPath, $maxThread, $outDir) = &readParameters();#->932

#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 1_defineHardCodedParam
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineHardCodedParam" num="1">
my $margin = 30;
my $plotMinPos = -$margin;
my $plotMaxPos = $margin;

my $minTotalCov = 10;
my $maxPosCov = 999999;
my $paramTag = "range$plotMinPos"."to"."$plotMaxPos.minTotalCov$minTotalCov.maxPosCov$maxPosCov";
my $resultDir = "$outDir/$paramTag/";
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 2_defineOutDirPath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutDirPath" num="2">
my @mkDirAry;
my $resultLogDir = "$resultDir/log/"; push @mkDirAry, $resultLogDir;
my $ggplotGeneralDirHsh_ref = {};
my $ggplotIndivGeneDirHsh_ref = {};
my @ggplotFileTypeAry = qw /dat pdf R log/;
foreach my $fileType (@ggplotFileTypeAry) {
	$ggplotGeneralDirHsh_ref->{$fileType} = "$resultDir/ggplotGeneral/$fileType"; push @mkDirAry, $ggplotGeneralDirHsh_ref->{$fileType};
	$ggplotIndivGeneDirHsh_ref->{$fileType} = "$resultDir/ggplotIndivGene/$fileType"; push @mkDirAry, $ggplotIndivGeneDirHsh_ref->{$fileType};
}
foreach my $dir (@mkDirAry) {system ("mkdir -pm 777 $dir");}

#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 3_defineOutFilePath
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="defineOutFilePath" num="3">
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 4_processInputData
#	primaryDependOnSub: checkGeneInfo|249, getCtgryGeneInfo|390, getIndivCntgCovPlsPath|427, readGFF_oneRNAPerGene|779, readMultiFasta|878, zipUnzipCntgCovInPlsPathHsh|989
#	secondaryDependOnSub: currentTime|277, reportStatus|968
#
#<section ID="processInputData" num="4">
########################################################################## 

my ($fastaHsh_ref) = &readMultiFasta($fastaPath);#->878

my ($geneInfoHsh_ref) = &readGFF_oneRNAPerGene($gffPath);#->779
&checkGeneInfo($geneInfoHsh_ref);#->249
#----------Get bfmRNA and mRNA ranges
my @mRNAAry = qw/mRNA/;
my ($mRNAInfoHsh_ref, $mRNAByCntgHsh_ref)= &getCtgryGeneInfo($geneInfoHsh_ref, \@mRNAAry);#->390

my ($pileupStorablePathHsh_ref) = &getIndivCntgCovPlsPath($rd5EndStorableIndexPath);#->427
&zipUnzipCntgCovInPlsPathHsh('unzip', $pileupStorablePathHsh_ref);#->989

#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 5_getAndPlotCDSCoverage
#	primaryDependOnSub: getCDSCoverageOnCntg|322, plotFrameFraction|563, plotPooledCDSCoverage|615
#	secondaryDependOnSub: generateThreadHshWithRandomCntg|295, getIndividualGeneCoverage|460, ggplotBarChart|493, ggplotXYLinesMultipleSamples|527, reportStatus|968
#
#<section ID="getAndPlotCDSCoverage" num="5">
my ($CDSCovHsh_ref) = &getCDSCoverageOnCntg($pileupStorablePathHsh_ref, $mRNAInfoHsh_ref, $mRNAByCntgHsh_ref, $maxThread, $margin);#->322
&plotPooledCDSCoverage($CDSCovHsh_ref, $mRNAInfoHsh_ref, $margin, $plotMinPos, $plotMaxPos, $ggplotGeneralDirHsh_ref, $minTotalCov, $maxPosCov, $resultLogDir);#->615
&plotFrameFraction($CDSCovHsh_ref, $margin, $ggplotGeneralDirHsh_ref);#->563
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 6_printStatisticsLog
#	primaryDependOnSub: >none
#	secondaryDependOnSub: >none
#
#<section ID="printStatisticsLog" num="6">
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
#	section 7_finishingTasks
#	primaryDependOnSub: printCMDLogOrFinishMessage|746
#	secondaryDependOnSub: currentTime|277
#
#<section ID="finishingTasks" num="7">
&printCMDLogOrFinishMessage("finishMessage");#->746
close DEBUGLOG;
#<\section>
#====================================================================================================================================================#

#====================================================================================================================================================#
}	#Main sections lexical scope ends
#====================================================================================================================================================#

#====================================================================================================================================================#
#List of subroutines by category
#
#	coverage [n=1]:
#		getIndividualGeneCoverage
#
#	fasta [n=1]:
#		readMultiFasta
#
#	general [n=8]:
#		checkGeneInfo, currentTime, getCtgryGeneInfo
#		printCMDLogOrFinishMessage, readGFF_oneRNAPerGene, readMultiFasta
#		readParameters, reportStatus
#
#	gff [n=3]:
#		checkGeneInfo, getCtgryGeneInfo, readGFF_oneRNAPerGene
#
#	ggplot [n=2]:
#		ggplotBarChart, ggplotXYLinesMultipleSamples
#
#	multithread [n=1]:
#		generateThreadHshWithRandomCntg
#
#	plotInR [n=2]:
#		ggplotBarChart, ggplotXYLinesMultipleSamples
#
#	reporting [n=1]:
#		currentTime
#
#	storable [n=2]:
#		getIndivCntgCovPlsPath, zipUnzipCntgCovInPlsPathHsh
#
#	unassigned [n=3]:
#		getCDSCoverageOnCntg, plotFrameFraction, plotPooledCDSCoverage
#
#====================================================================================================================================================#

sub checkGeneInfo {
#....................................................................................................................................................#
#	subroutineCategory: general, gff
#	dependOnSub: reportStatus|968
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|154
#	secondaryAppearInSection: >none
#	input: $geneInfoHsh_ref
#	output: 
#	toCall: &checkGeneInfo($geneInfoHsh_ref);
#	calledInLine: 164
#....................................................................................................................................................#
	
	my ($geneInfoHsh_ref) = @_;
	
	&reportStatus("Checking gene categories", 0, "\n");#->968
	my $ctrgyCountHsh_ref = {};
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {
		my $ctgry = $geneInfoHsh_ref->{$geneID}{'ctgry'};
		$ctrgyCountHsh_ref->{$ctgry}++;
	}
	
	foreach my $ctgry (sort keys %{$ctrgyCountHsh_ref}) {
		&reportStatus("Item in $ctgry = $ctrgyCountHsh_ref->{$ctgry}", 0, "\n");#->968
	}
	
	return ();
}
sub currentTime {
#....................................................................................................................................................#
#	subroutineCategory: general, reporting
#	dependOnSub: >none
#	appearInSub: getCtgryGeneInfo|390, printCMDLogOrFinishMessage|746, readGFF_oneRNAPerGene|779, reportStatus|968
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 0_startingTasks|94, 4_processInputData|154, 7_finishingTasks|197
#	input: none
#	output: $runTime
#	toCall: my ($runTime) = &currentTime();
#	calledInLine: 409, 422, 766, 769, 774, 799, 984
#....................................................................................................................................................#
	
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
	
	return $runTime;
}
sub generateThreadHshWithRandomCntg {
#....................................................................................................................................................#
#	subroutineCategory: multithread
#	dependOnSub: >none
#	appearInSub: getCDSCoverageOnCntg|322
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 5_getAndPlotCDSCoverage|176
#	input: $cntgAry_ref, $threadToSpawn
#	output: $randCntgInThreadHsh_ref
#	toCall: my ($randCntgInThreadHsh_ref) = &generateThreadHshWithRandomCntg($threadToSpawn, $cntgAry_ref);
#	calledInLine: 336
#....................................................................................................................................................#

	my ($threadToSpawn, $cntgAry_ref) = @_;

	my @shuffleCntgAry = shuffle(@{$cntgAry_ref});
	my $threadNum = 1;
	my $randCntgInThreadHsh_ref = {};
	foreach my $cntg (@{$cntgAry_ref}) {
		$threadNum = 1 if $threadNum > $threadToSpawn;
		push @{$randCntgInThreadHsh_ref->{$threadNum}}, $cntg;
		$threadNum++;
	}
	
	return $randCntgInThreadHsh_ref;

}
sub getCDSCoverageOnCntg {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: generateThreadHshWithRandomCntg|295, getIndividualGeneCoverage|460, reportStatus|968
#	appearInSub: >none
#	primaryAppearInSection: 5_getAndPlotCDSCoverage|176
#	secondaryAppearInSection: >none
#	input: $mRNAByCntgHsh_ref, $mRNAInfoHsh_ref, $margin, $maxThread, $pileupStorablePathHsh_ref
#	output: $CDSCovHsh_ref
#	toCall: my ($CDSCovHsh_ref) = &getCDSCoverageOnCntg($pileupStorablePathHsh_ref, $mRNAInfoHsh_ref, $mRNAByCntgHsh_ref, $maxThread, $margin);
#	calledInLine: 181
#....................................................................................................................................................#
	my ($pileupStorablePathHsh_ref, $mRNAInfoHsh_ref, $mRNAByCntgHsh_ref, $maxThread, $margin) = @_;
	
	my @cntgAry = (keys %{$pileupStorablePathHsh_ref});
	my $randCntgInThreadHsh_ref = &generateThreadHshWithRandomCntg($maxThread, \@cntgAry);#->295
	my $cntgProc :shared = 0;
	my %threadHsh = ();

	foreach my $threadNum (sort {$a <=> $b} keys %{$randCntgInThreadHsh_ref}) {
		my $cntgAry_ref = \@{$randCntgInThreadHsh_ref->{$threadNum}};
		my $cntgNum = @{$randCntgInThreadHsh_ref->{$threadNum}};
		&reportStatus("$cntgNum cntgs spawned to thread $threadNum", 20, "\n");#->968

		#---spawn a new thread
		($threadHsh{$threadNum}) = threads->new(#---refer to http://www.perlmonks.org/?node_id=966781
		
			sub {
				my ($cntgAry_ref) = @_;
				my $CDSCovInThrHsh_ref = {};
				foreach my $cntg (@{$cntgAry_ref}) {
					$cntgProc++;
					my $cntgCovAry_ref = retrieve($pileupStorablePathHsh_ref->{$cntg});
					next if not exists $mRNAByCntgHsh_ref->{$cntg};
					foreach my $geneID (keys %{$mRNAByCntgHsh_ref->{$cntg}}) {
						next if not $mRNAInfoHsh_ref->{$geneID}{'CDSRng'};
						my ($posCovAry_ref) = &getIndividualGeneCoverage($cntgCovAry_ref, $mRNAInfoHsh_ref, $geneID, $margin);#->460
						$CDSCovInThrHsh_ref->{$geneID} = $posCovAry_ref;
					}
					&reportStatus("Finished counting features on $cntgProc cntg", 20, "\r");#->968
				}
				return ($CDSCovInThrHsh_ref);
			}
			,($cntgAry_ref)
		);
	}
	
	my $CDSCovHsh_ref = {};
	
	while (keys %threadHsh) {
		foreach my $threadNum (keys %threadHsh) {
			my $thr = $threadHsh{$threadNum};
			if (not $thr->is_running()) {
				my ($CDSCovInThrHsh_ref) = $thr->join;
				foreach my $geneID (keys %{$CDSCovInThrHsh_ref}) {
					$CDSCovHsh_ref->{$geneID} = $CDSCovInThrHsh_ref->{$geneID};
				}
				delete $threadHsh{$threadNum};
			}
		}
		sleep 1;
	}

	my $numCDS = scalar(keys %{$CDSCovHsh_ref});
	&reportStatus("CDS coverage of $numCDS coding genes were stored", 10, "\n");#->968
	
	return ($CDSCovHsh_ref);
}
sub getCtgryGeneInfo {
#....................................................................................................................................................#
#	subroutineCategory: general, gff
#	dependOnSub: currentTime|277
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|154
#	secondaryAppearInSection: >none
#	input: $ctgryAry_ref, $geneInfoHsh_ref
#	output: $geneCtgryByCntgHsh_ref, $geneCtgryInfoHsh_ref
#	toCall: my ($geneCtgryInfoHsh_ref, $geneCtgryByCntgHsh_ref) = &getCtgryGeneInfo($geneInfoHsh_ref, $ctgryAry_ref);
#	calledInLine: 167
#....................................................................................................................................................#

	my ($geneInfoHsh_ref, $ctgryAry_ref) = @_;
	
	my $geneCtgryInfoHsh_ref = {};
	my $geneCtgryByCntgHsh_ref = {};
	
	my $ctgryStr = join ",", @{$ctgryAry_ref};

	print "[".&currentTime()."] Filtering GFF on cgtry $ctgryStr.\n";#->277
	
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {
		my $ctgry = $geneInfoHsh_ref->{$geneID}{'ctgry'};
		my $cntg = $geneInfoHsh_ref->{$geneID}{'cntg'};
		if (grep /^$ctgry$/, @{$ctgryAry_ref}) {
			%{$geneCtgryInfoHsh_ref->{$geneID}} = %{$geneInfoHsh_ref->{$geneID}};
			$geneCtgryByCntgHsh_ref->{$cntg}{$geneID}++;
		}
	}
	
	my $numGene = keys %{$geneCtgryInfoHsh_ref};
	
	print "[".&currentTime()."] $numGene gene filtered on cgtry $ctgryStr.\n";#->277
	
	return $geneCtgryInfoHsh_ref, $geneCtgryByCntgHsh_ref;
}
sub getIndivCntgCovPlsPath {
#....................................................................................................................................................#
#	subroutineCategory: storable
#	dependOnSub: reportStatus|968
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|154
#	secondaryAppearInSection: >none
#	input: $cntgCovPlsIndexPath
#	output: $cntgCovInPlsPathHsh_ref
#	toCall: my ($cntgCovInPlsPathHsh_ref) = &getIndivCntgCovPlsPath($cntgCovPlsIndexPath);
#	calledInLine: 169
#....................................................................................................................................................#
	
	my ($cntgCovPlsIndexPath) = @_;

	$cntgCovPlsIndexPath =~ s/\.gz$//;

	my $cntgCovInPlsPathHsh_ref = {};
	
	system ("gzip -fd $cntgCovPlsIndexPath.gz") if -s "$cntgCovPlsIndexPath.gz";
	my %plsIndexHsh = %{retrieve($cntgCovPlsIndexPath)};
	
	my (undef, $cntgCovStroableDir, undef) = fileparse($cntgCovPlsIndexPath, qr/\.[^.]*/);
	foreach my $cntg (keys %plsIndexHsh) {
		my $cntgCovPlsPath = "$cntgCovStroableDir/$plsIndexHsh{$cntg}";
		die "cntgCovPlsPath $cntgCovPlsPath is invalid\n" if ((not -s $cntgCovPlsPath) and (not -s $cntgCovPlsPath.".gz"));
		$cntgCovInPlsPathHsh_ref->{$cntg} = $cntgCovPlsPath;
	}
	my $numCntg = keys %{$cntgCovInPlsPathHsh_ref};
	&reportStatus("pls path of $numCntg contig stored.", 0, "\n");#->968
	
	return $cntgCovInPlsPathHsh_ref;
}
sub getIndividualGeneCoverage {
#....................................................................................................................................................#
#	subroutineCategory: coverage
#	dependOnSub: >none
#	appearInSub: getCDSCoverageOnCntg|322
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 5_getAndPlotCDSCoverage|176
#	input: $cntgCovAry_ref, $geneID, $mRNAInfoHsh_ref, $margin
#	output: $posCovAry_ref
#	toCall: my ($posCovAry_ref) = &getIndividualGeneCoverage($cntgCovAry_ref, $mRNAInfoHsh_ref, $geneID, $margin);
#	calledInLine: 357
#....................................................................................................................................................#
	my ($cntgCovAry_ref, $mRNAInfoHsh_ref, $geneID, $margin) = @_;
	
	#---create the margin
	my @CDSRng = sort {$a <=> $b} @{$mRNAInfoHsh_ref->{$geneID}{'CDSRng'}};
	unshift @CDSRng, ($CDSRng[0]-1-$margin, $CDSRng[0]-1);
	push @CDSRng, ($CDSRng[-1]+1, $CDSRng[-1]+1+$margin);
	
	my $posCovAry_ref = ();
	for (my $i=0; $i < $#CDSRng; $i += 2) {
		foreach my $j ($CDSRng[$i]-1..$CDSRng[$i+1]-1) {
			my %tmpCovHsh = ('+'=>0, '-'=>0);
			($tmpCovHsh{'+'}, $tmpCovHsh{'-'}) = split /,/, $cntgCovAry_ref->[$j] if ($cntgCovAry_ref->[$j]);
			#die "$mRNAInfoHsh_ref->{$geneID}{'strnd'}\t$geneID\n" if not defined $tmpCovHsh{$mRNAInfoHsh_ref->{$geneID}{'strnd'}};
			push @{$posCovAry_ref}, $tmpCovHsh{$mRNAInfoHsh_ref->{$geneID}{'strnd'}};
		}
	}
	
	@{$posCovAry_ref} = reverse @{$posCovAry_ref} if $mRNAInfoHsh_ref->{$geneID}{'strnd'} eq '-';
	
	return ($posCovAry_ref);
}
sub ggplotBarChart {
#....................................................................................................................................................#
#	subroutineCategory: ggplot, plotInR
#	dependOnSub: >none
#	appearInSub: plotFrameFraction|563
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 5_getAndPlotCDSCoverage|176
#	input: $RScriptPath, $XAXis, $YAxis, $dataPath, $extraArg, $extraStatment, $height, $logPath, $pdfPath, $plotDataHsh_ref, $width
#	output: 
#	toCall: &ggplotBarChart($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $XAXis, $YAxis, $extraArg, $extraStatment, $height, $width);
#	calledInLine: 608
#....................................................................................................................................................#
	my ($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $XAXis, $YAxis, $extraArg, $extraStatment, $height, $width) = @_;
	
	open (PLOTDATA, ">", $dataPath);
	print PLOTDATA join "", (join "\t", ($XAXis, $YAxis)), "\n";
	foreach my $XVal (sort keys %{$plotDataHsh_ref}) {
		my $YVal = $plotDataHsh_ref->{$XVal};
		print PLOTDATA join "", (join "\t", ($XVal, $YVal)), "\n";
	}
	close PLOTDATA;

	open (R, ">", $RScriptPath);
	print R "library(ggplot2)"."\n";
	print R "$extraStatment"."\n";
	print R "dataFrame = read.table(file=\"$dataPath\", header=TRUE)"."\n";
	print R "ggplot(dataFrame, aes(x=$XAXis, y=$YAxis, fill=$YAxis)) + geom_bar(stat=\"identity\") $extraArg"."\n";
	print R "ggsave(file=\"$pdfPath\", height=$height, width=$width)\n";
	close R;
	
	system ("R --slave --vanilla --file=$RScriptPath &>$logPath");

	return ();
}
sub ggplotXYLinesMultipleSamples {
#....................................................................................................................................................#
#	subroutineCategory: ggplot, plotInR
#	dependOnSub: >none
#	appearInSub: plotPooledCDSCoverage|615
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 5_getAndPlotCDSCoverage|176
#	input: $RScriptPath, $XAXis, $YAxis, $YVariable, $dataPath, $extraArg, $extraStatment, $height, $logPath, $pdfPath, $plotDataHsh_ref, $width
#	output: none
#	toCall: &ggplotXYLinesMultipleSamples($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $XAXis, $YAxis, $YVariable, $extraArg, $extraStatment, $height, $width);
#	calledInLine: 739
#....................................................................................................................................................#

	my ($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $XAXis, $YAxis, $YVariable, $extraArg, $extraStatment, $height, $width) = @_;
	
	open (PLOTDATA, ">", $dataPath);
	print PLOTDATA join "", (join "\t", ($YVariable, $YAxis, $XAXis)), "\n";
	foreach my $YCategory (sort keys %{$plotDataHsh_ref}) {
		foreach my $XVal (sort {$a <=> $b} keys %{$plotDataHsh_ref->{$YCategory}}) {
			my $YVal = $plotDataHsh_ref->{$YCategory}{$XVal};
			print PLOTDATA join "", (join "\t", ($YCategory, $YVal, $XVal)), "\n";
		}
	}
	close PLOTDATA;

	open (R, ">", $RScriptPath);
	print R "library(ggplot2)"."\n";
	print R "$extraStatment"."\n";
	print R "dataFrame = read.table(file=\"$dataPath\", header=TRUE)"."\n";
	print R "ggplot(dataFrame, aes(x=$XAXis, y=$YAxis, colour=$YVariable)) + geom_line() $extraArg"."\n";
	print R "ggsave(file=\"$pdfPath\", height=$height, width=$width)\n";
	close R;
	
	system ("R --slave --vanilla --file=$RScriptPath &>$logPath");

}
sub plotFrameFraction {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: ggplotBarChart|493, reportStatus|968
#	appearInSub: >none
#	primaryAppearInSection: 5_getAndPlotCDSCoverage|176
#	secondaryAppearInSection: >none
#	input: $CDSCovHsh_ref, $ggplotGeneralDirHsh_ref, $margin
#	output: 
#	toCall: &plotFrameFraction($CDSCovHsh_ref, $margin, $ggplotGeneralDirHsh_ref);
#	calledInLine: 183
#....................................................................................................................................................#
	my ($CDSCovHsh_ref, $margin, $ggplotGeneralDirHsh_ref) = @_;
	
	my $pooledPlotDataHsh_ref = {};
	$pooledPlotDataHsh_ref->{'0'} = 0;
	$pooledPlotDataHsh_ref->{'1'} = 0;
	$pooledPlotDataHsh_ref->{'2'} = 0;

	&reportStatus("Pooling count at different coding frame", 10, "\n");#->968

	foreach my $mRNAID (keys %{$CDSCovHsh_ref}) {
		my @posCovAry = @{$CDSCovHsh_ref->{$mRNAID}}[$margin+1..$#{$CDSCovHsh_ref->{$mRNAID}}-$margin+1];
		my $frame = 0;
		foreach my $cov (@posCovAry) {
			$pooledPlotDataHsh_ref->{$frame} += $cov;
			$frame++;
			$frame = 0 if $frame == 3;
		}
	}
	
	{
		my $plotDataHsh_ref = $pooledPlotDataHsh_ref;
		my $ggplotDirHsh_ref = $ggplotGeneralDirHsh_ref;
		my $nameTag = "codingFrame.readProportion";
		my $dataPath = $ggplotDirHsh_ref->{'dat'}."/$nameTag.dat";
		my $pdfPath = $ggplotDirHsh_ref->{'pdf'}."/$nameTag.pdf";
		my $RScriptPath = $ggplotDirHsh_ref->{'R'}."/$nameTag.R";
		my $logPath = $ggplotDirHsh_ref->{'log'}."/$nameTag.log";
		my $XAXis = 'condingFrame';
		my $YAxis = 'readNumber';
		my $extraStatment = "";
		my $extraArg = " + ggtitle(\"proportion of reads in coding frame\")";
		#my $extraArg = '';
		my $height = 10;
		my $width = 6;
		&ggplotBarChart($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $XAXis, $YAxis, $extraArg, $extraStatment, $height, $width);#->493
	}
	
	return ();
}
sub plotPooledCDSCoverage {
#....................................................................................................................................................#
#	subroutineCategory: unassigned
#	dependOnSub: ggplotXYLinesMultipleSamples|527, reportStatus|968
#	appearInSub: >none
#	primaryAppearInSection: 5_getAndPlotCDSCoverage|176
#	secondaryAppearInSection: >none
#	input: $CDSCovHsh_ref, $ggplotGeneralDirHsh_ref, $mRNAInfoHsh_ref, $margin, $maxPosCov, $minTotalCov, $plotMaxPos, $plotMinPos, $resultLogDir
#	output: 
#	toCall: &plotPooledCDSCoverage($CDSCovHsh_ref, $mRNAInfoHsh_ref, $margin, $plotMinPos, $plotMaxPos, $ggplotGeneralDirHsh_ref, $minTotalCov, $maxPosCov, $resultLogDir);
#	calledInLine: 182
#....................................................................................................................................................#
	my ($CDSCovHsh_ref, $mRNAInfoHsh_ref, $margin, $plotMinPos, $plotMaxPos, $ggplotGeneralDirHsh_ref, $minTotalCov, $maxPosCov, $resultLogDir) = @_;
	
	foreach my $referencePoint (qw/ATG TAA/) {
		my $numGene = 0;
		my $overallMaxPosCov = 0;
		my $totalCountHsh_ref = {};
		my $rltvCountHsh_ref = {};
		my %indivDataHsh = ();
	
		&reportStatus("Pooling coverage around $referencePoint", 10, "\n");#->968
			
		foreach my $mRNAID (keys %{$CDSCovHsh_ref}) {
			my $totalPosCov = 0;
			my $maxmRNAPosCov = 0;
			my @posCovAry = ();
			my $rltvPos;
			if ($referencePoint eq 'ATG') {
				@posCovAry = @{$CDSCovHsh_ref->{$mRNAID}};
				$rltvPos = -1*$margin;
			} elsif ($referencePoint eq 'TAA') {
				@posCovAry = reverse @{$CDSCovHsh_ref->{$mRNAID}};
				$rltvPos = $margin;
			}
			
			my $indivPlotDataHsh_ref = {};
			
			#---check CDS length
			my $CDSLen = @posCovAry - $margin*2;
			next if ($CDSLen < abs($plotMaxPos) and $referencePoint eq 'ATG');
			next if ($CDSLen < abs($plotMinPos) and $referencePoint eq 'TAA');
			
			#---get individual coverage
			foreach my $cov (@posCovAry) {
				if ($rltvPos > $plotMinPos and $rltvPos <= $plotMaxPos) {
					$totalPosCov += $cov;
					$indivPlotDataHsh_ref->{'CDS'}{$rltvPos} = $cov;
					$maxmRNAPosCov = $cov if $cov > $maxmRNAPosCov;
				}
				last if $rltvPos > $plotMaxPos;
	
				if ($referencePoint eq 'ATG') {
					$rltvPos++;
				} elsif ($referencePoint eq 'TAA') {
					$rltvPos--;
				}
			}

			#---pool coverage
			if ($maxmRNAPosCov <= $maxPosCov and $totalPosCov >= $minTotalCov) {
				$numGene++;
				$overallMaxPosCov = $maxmRNAPosCov if $maxmRNAPosCov > $overallMaxPosCov;
				my $totalCountAry_ref = [];
				my $rltvCountAry_ref = [];
				foreach my $rltvPos (sort {$a <=> $b} keys %{$indivPlotDataHsh_ref->{'CDS'}}) {
					my $totalCount = $indivPlotDataHsh_ref->{'CDS'}{$rltvPos};
					$totalCountHsh_ref->{'CDS'}{$rltvPos} = 0 if not $totalCountHsh_ref->{'CDS'}{$rltvPos};
					$totalCountHsh_ref->{'CDS'}{$rltvPos} += $totalCount;
					push @{$totalCountAry_ref}, $totalCount;
					
					my $rltvCount = sprintf "%.5f", $totalCount/$maxmRNAPosCov;
					$rltvCountHsh_ref->{'CDS'}{$rltvPos} = 0 if not $rltvCountHsh_ref->{'CDS'}{$rltvPos};
					$rltvCountHsh_ref->{'CDS'}{$rltvPos} += $rltvCount;
					push @{$rltvCountAry_ref}, $rltvCount;
				}
				$indivDataHsh{'rltvCount'}{$mRNAID} = $rltvCountAry_ref;
				$indivDataHsh{'totalCount'}{$mRNAID} = $totalCountAry_ref;
			}
		}
		
		#---print indiv data log
		my @headerAry = ("mRNAID", (sort {$a <=> $b} keys %{$rltvCountHsh_ref->{'CDS'}}));
		foreach my $countType (sort keys %indivDataHsh) {
			open INDIVDATALOG, ">", "$resultLogDir/individual.mRNA.$referencePoint.$countType.xls";
			print INDIVDATALOG join "", (join "\t", (@headerAry)), "\n";
			foreach my $mRNAID (sort keys %{$indivDataHsh{$countType}}) {
				my @outputAry = ($mRNAID, @{$indivDataHsh{$countType}{$mRNAID}});
				print INDIVDATALOG join "", (join "\t", (@outputAry)), "\n";
			}
			close INDIVDATALOG;
		}
		
		#---divide the rltv count by number of genes 
		$rltvCountHsh_ref->{'CDS'}{$_} = $rltvCountHsh_ref->{'CDS'}{$_}/$numGene foreach (keys %{$rltvCountHsh_ref->{'CDS'}});
		
		my %tmpPlotDataHsh = (
			"totalCount" => $totalCountHsh_ref,
			"rltvCount" => $rltvCountHsh_ref,
		);
		
		foreach my $countType (keys %tmpPlotDataHsh) {

			my $vLineStart = $plotMinPos + 1;
			my $vLineEnd = $plotMaxPos -2;
			if ($referencePoint eq 'TAA') {
				$vLineStart--;
				$vLineEnd--;
			}
			
			my $plotDataHsh_ref = $tmpPlotDataHsh{$countType};
			my $ggplotDirHsh_ref = $ggplotGeneralDirHsh_ref;
			my $nameTag = "$referencePoint.$countType";
			my $dataPath = $ggplotDirHsh_ref->{'dat'}."/$nameTag.dat";
			my $pdfPath = $ggplotDirHsh_ref->{'pdf'}."/$nameTag.pdf";
			my $RScriptPath = $ggplotDirHsh_ref->{'R'}."/$nameTag.R";
			my $logPath = $ggplotDirHsh_ref->{'log'}."/$nameTag.log";
			my $XAXis = 'relativePositon';
			my $YAxis = 'coverage';
			my $YVariable = 'covType';
			my $extraStatment = "vLinePos <- seq($vLineStart,$vLineEnd,by=3); df <- data.frame(vLinePos);"; #---based on http://stackoverflow.com/questions/15165383/error-when-adding-vertical-lines-in-ggplot2
			my $extraArg = " + ggtitle(\"N=$numGene overallMaxPosCov=$overallMaxPosCov\") + scale_x_continuous(breaks=seq($plotMinPos, $plotMaxPos, by=1)) + geom_vline(data=df, aes(xintercept = vLinePos), alpha=0.2, size=4);";
			#my $extraArg = '';
			my $height = 6;
			my $width = 30;
			&ggplotXYLinesMultipleSamples($plotDataHsh_ref, $dataPath, $pdfPath, $RScriptPath, $logPath, $XAXis, $YAxis, $YVariable, $extraArg, $extraStatment, $height, $width);#->527
		}
	}
	
	return ();
}
sub printCMDLogOrFinishMessage {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|277
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|94, 7_finishingTasks|197
#	secondaryAppearInSection: >none
#	input: $CMDLogOrFinishMessage
#	output: none
#	toCall: &printCMDLogOrFinishMessage($CMDLogOrFinishMessage);
#	calledInLine: 100, 202
#....................................................................................................................................................#

	my ($CMDLogOrFinishMessage) = @_;
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $absoluteScriptPath = abs_path($0);
		my $dirPath = dirname(rel2abs($absoluteScriptPath));
		my ($scriptName, $callScriptPath, $scriptSuffix) = fileparse($absoluteScriptPath, qr/\.[^.]*/);
		open (CMDLOG, ">>$dirPath/$scriptName.cmd.log.txt"); #---append the CMD log file
		print CMDLOG "[".&currentTime()."]\t"."$dirPath/$scriptName$scriptSuffix ".(join " ", @ARGV)."\n";#->277
		close CMDLOG;
		print "\n=========================================================================\n";
		print "[".&currentTime()."] starts running ...... \n";#->277
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		print "\n=========================================================================\n";
		print "[".&currentTime()."] finished running .......\n";#->277
		print "=========================================================================\n\n";
	}
}
sub readGFF_oneRNAPerGene {
#....................................................................................................................................................#
#	subroutineCategory: general, gff
#	dependOnSub: currentTime|277
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|154
#	secondaryAppearInSection: >none
#	input: $gffPath
#	output: $geneInfoHsh_ref
#	toCall: my ($geneInfoHsh_ref) = &readGFF_oneRNAPerGene($gffPath);
#	calledInLine: 163
#....................................................................................................................................................#

	my ($gffPath) = @_;

	my $geneInfoHsh_ref = {};
	
	#---read the gff
	my $geneByRNAHsh_ref = {};

	open (GFF, $gffPath);
	print "[".&currentTime()."] Reading: $gffPath\n";#->277
	while (my $theLine = <GFF>) {

		chomp $theLine;
		
		last if $theLine =~ m/^##FASTA/;
		
		if ($theLine !~ m/^\#|^\@/ and $theLine !~ m/\tsupercontig\t/) {

			my ($seq, undef, $geneCategory, $featureStart, $featureEnd, undef, $geneStrd, undef, $dscrptns) = split (/\t/, $theLine);
			
			#----assigne all non -/+ will be treated as plus
			$geneStrd = "+" if (($geneStrd ne "-") and ($geneStrd ne "+"));
			
			my @dscrptnsSplt = split /;/, $dscrptns;
			my ($unqID, $parent);
			my $geneName = "unknown";
			foreach my $theDscptn (@dscrptnsSplt) {
				if ($theDscptn =~ m/^ID=/) {$unqID = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^Parent=/) {$parent = substr ($theDscptn, index ($theDscptn, "=")+1);}
				if ($theDscptn =~ m/^description=/) {$geneName = substr ($theDscptn, index ($theDscptn, "=")+1);}
			}

			if ($geneCategory eq "gene") {#---gene
				
				my $geneID = $unqID;
				
				$geneInfoHsh_ref->{$geneID}{'strnd'} = $geneStrd;
				$geneInfoHsh_ref->{$geneID}{'cntg'} = $seq;
				$geneInfoHsh_ref->{$geneID}{'description'} = uri_unescape($geneName);
				$geneInfoHsh_ref->{$geneID}{'description'} =~ s/\+/ /g;
				@{$geneInfoHsh_ref->{$geneID}{'geneRng'}} = ($featureStart, $featureEnd);

			} elsif ($geneCategory eq "CDS") {#---Only for coding genes
				
				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh_ref->{$RNAID};
				push @{$geneInfoHsh_ref->{$geneID}{'CDSRng'}}, ($featureStart, $featureEnd);
				
			} elsif ($geneCategory eq "exon") {#---exon, may be exons of alternative transcripts, wiull sort out later
				my $RNAID = $parent;
				my $geneID = $geneByRNAHsh_ref->{$RNAID};
				push @{$geneInfoHsh_ref->{$geneID}{'exonRng'}}, ($featureStart, $featureEnd);
				
			} else {#---can be tRNA, rRNA, mRNA, repRNA, ncRNA
				my $RNAID = $unqID;
				my $geneID = $parent;
				$geneByRNAHsh_ref->{$RNAID} = $geneID;
				$geneInfoHsh_ref->{$geneID}{'ctgry'} = $geneCategory;
				@{$geneInfoHsh_ref->{$geneID}{'RNARng'}} = ($featureStart, $featureEnd);
				$geneInfoHsh_ref->{$geneID}{'RNAID'} = $RNAID;
			}
		}#---end of if (($theLine !~ m/^\#|^\@/) and ($theLine !~ m/\tsupercontig\t/)) {
	}#---end of while (my $theLine = <INFILE>)
	close GFF;
	
	#---get the UTR if any
	my $minUTRLength = 10;
	foreach my $geneID (keys %{$geneInfoHsh_ref}) {

		if (exists $geneInfoHsh_ref->{$geneID}{'CDSRng'}) {
			my $exonMin = min(@{$geneInfoHsh_ref->{$geneID}{'exonRng'}});
			my $exonMax = max(@{$geneInfoHsh_ref->{$geneID}{'exonRng'}});
			my $CDSMin = min(@{$geneInfoHsh_ref->{$geneID}{'CDSRng'}});
			my $CDSMax = max(@{$geneInfoHsh_ref->{$geneID}{'CDSRng'}});

			if ($geneInfoHsh_ref->{$geneID}{'strnd'} eq '+') {
				@{$geneInfoHsh_ref->{$geneID}{'UTR5Rng'}} = ($exonMin, $CDSMin-1) if ($CDSMin-$exonMin > $minUTRLength);
				@{$geneInfoHsh_ref->{$geneID}{'UTR3Rng'}} = ($CDSMax+1, $exonMax) if ($exonMax-$CDSMax > $minUTRLength);
			} else {
				@{$geneInfoHsh_ref->{$geneID}{'UTR3Rng'}} = ($exonMin, $CDSMin-1) if ($CDSMin-$exonMin > $minUTRLength);
				@{$geneInfoHsh_ref->{$geneID}{'UTR5Rng'}} = ($CDSMax+1, $exonMax) if ($exonMax-$CDSMax > $minUTRLength);
			}
		}
	}
	
	return ($geneInfoHsh_ref);
}
sub readMultiFasta {
#....................................................................................................................................................#
#	subroutineCategory: fasta, general
#	dependOnSub: reportStatus|968
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|154
#	secondaryAppearInSection: >none
#	input: $fastaPath
#	output: $fastaHsh_ref
#	toCall: my ($fastaHsh_ref) = &readMultiFasta($fastaPath);
#	calledInLine: 161
#....................................................................................................................................................#

	my ($fastaPath) = @_;

	my ($seq, $seqName);
	my $fastaHsh_ref = {};
	my $i = 0;

	&reportStatus("Reading: $fastaPath", 0, "\n");#->968
	
	open (INFILE, $fastaPath);
	chomp (my $curntLine = <INFILE>); #get the first line
	while (my $nextLine = <INFILE>) {
		chomp $nextLine;
		
		#---Only two types of line in current line, the header or seq
		if ($curntLine =~ m/^>/) {#-- header line
			my @theLineSplt = split (/ *\|/, $curntLine);
			$seqName = $theLineSplt[0]; #---get the first tag
			$seqName =~ s/ //g; #---remove space
			$seqName =~ s/>//g; #---remove space
		} else {#--seq line
			$seq = $seq.$curntLine;
		}
		
		#---check if next line has a > or that's the end of file
		if ($nextLine =~ m/^>/) {
			$seq =~ tr/a-z/A-Z/;
			$fastaHsh_ref->{$seqName} = $seq;
			$seq = "";
		} elsif (eof(INFILE)) {#---this is the last line
			$seq =~ tr/a-z/A-Z/;
			$seq = $seq.$nextLine;
			$fastaHsh_ref->{$seqName} = $seq;
		}
		
		#---next line becomes current line
		$curntLine = $nextLine;
	}

	close INFILE;
	return ($fastaHsh_ref);
}
sub readParameters {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: >none
#	appearInSub: >none
#	primaryAppearInSection: 0_startingTasks|94
#	secondaryAppearInSection: >none
#	input: none
#	output: $fastaPath, $gffPath, $maxThread, $outDir, $rd5EndStorableIndexPath
#	toCall: my ($rd5EndStorableIndexPath, $gffPath, $fastaPath, $maxThread, $outDir) = &readParameters();
#	calledInLine: 102
#....................................................................................................................................................#
	
	my ($rd5EndStorableIndexPath, $gffPath, $fastaPath, $maxThread, $outDir);
	
	my $dirPath = dirname(rel2abs($0));
	$outDir = "$dirPath/ribosomeFootprintAnalyzer/";
	$maxThread = 4;
	
	GetOptions 	("rd5EndStorableIndexPath=s"  => \$rd5EndStorableIndexPath,
				 "gffPath=s"  => \$gffPath,
				 "fastaPath=s"  => \$fastaPath,
				 "maxThread:i"  => \$maxThread,
				 "outDir:s"  => \$outDir)

	or die		("Error in command line arguments\n");
	
	#---check file
	foreach my $fileToCheck ($gffPath, $fastaPath, $rd5EndStorableIndexPath) {
		die "Can't read $fileToCheck" if not -s $fileToCheck;
	}

	system "mkdir -p -m 777 $outDir/";
	
	return($rd5EndStorableIndexPath, $gffPath, $fastaPath, $maxThread, $outDir);
}
sub reportStatus {
#....................................................................................................................................................#
#	subroutineCategory: general
#	dependOnSub: currentTime|277
#	appearInSub: checkGeneInfo|249, getCDSCoverageOnCntg|322, getIndivCntgCovPlsPath|427, plotFrameFraction|563, plotPooledCDSCoverage|615, readMultiFasta|878, zipUnzipCntgCovInPlsPathHsh|989
#	primaryAppearInSection: >none
#	secondaryAppearInSection: 4_processInputData|154, 5_getAndPlotCDSCoverage|176
#	input: $lineEnd, $message, $numTrailingSpace
#	output: 
#	toCall: &reportStatus($message, $numTrailingSpace, $lineEnd);
#	calledInLine: 263, 271, 343, 360, 385, 455, 581, 635, 896, 1004
#....................................................................................................................................................#
	my ($message, $numTrailingSpace, $lineEnd) = @_;

	my $trailingSpaces = '';
	$trailingSpaces .= " " for (1..$numTrailingSpace);
	
	print "[".&currentTime()."] ".$message.$trailingSpaces.$lineEnd;#->277

	return ();
}
sub zipUnzipCntgCovInPlsPathHsh {
#....................................................................................................................................................#
#	subroutineCategory: storable
#	dependOnSub: reportStatus|968
#	appearInSub: >none
#	primaryAppearInSection: 4_processInputData|154
#	secondaryAppearInSection: >none
#	input: $cntgCovInPlsPathHsh_ref, $zipUnzip
#	output: none
#	toCall: &zipUnzipCntgCovInPlsPathHsh($zipUnzip, $cntgCovInPlsPathHsh_ref);
#	calledInLine: 170
#....................................................................................................................................................#

	my ($zipUnzip, $cntgCovInPlsPathHsh_ref) = @_;
	
	foreach my $cntg (sort keys %{$cntgCovInPlsPathHsh_ref}) {
		&reportStatus("Trying to $zipUnzip cntg ary", 20, "\r");#->968

		my $cntgCovPlsPath = "$cntgCovInPlsPathHsh_ref->{$cntg}";

		if ($zipUnzip eq 'unzip') {
			system ("gzip -df $cntgCovPlsPath.gz") if (-s "$cntgCovPlsPath.gz");
		} else {
			system ("gzip -f $cntgCovPlsPath") if (-s "$cntgCovPlsPath");
		}
	}
	print "\n";
}

exit;
