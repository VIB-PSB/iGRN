#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use FileHandle;

my ($network, $mode, $exprData, $matrixData, $clusterData, $dhData, $cnsData, $chipData, $genie3Data, $knnData); #input using getopt

GetOptions(	"network=s"			=> \$network,
			"mode=s"			=> \$mode,
			"expr=s"			=> \$exprData,
			"matrix=s"			=> \$matrixData,
			"cluster=s"			=> \$clusterData,
			"DH=s"				=> \$dhData,
			"CNS=s"				=> \$cnsData,
			"ChIP=s"			=> \$chipData,
			"GENIE3=s"			=> \$genie3Data,
			"KNN=s"				=> \$knnData
) || die "There were errors reading the arguments\n";

my ($networkHash, $tf, $tg) = read_network($network);
my $combinedHash = {};
my $exprHash = {};
my $matrixHash = {};
my $clusterHash={};
my $dhHash = {};
my $cnsHash = {};
my $chipHash = {};
my $genie3Hash = {};
my $knnHash = {};

if($mode eq 'train'){
	$combinedHash = combine_interactions($networkHash, $tf, $tg);
	$exprHash = read_feature_network($exprData);
	$matrixHash = read_feature_network($matrixData);
	$clusterHash = read_feature_network($clusterData);
	$dhHash = read_feature_network($dhData);
	$cnsHash = read_feature_network($cnsData);
	$chipHash = read_feature_network($chipData);
	$genie3Hash = read_feature_network($genie3Data);
	$knnHash = read_feature_network($knnData);
	&make_train_matrix($combinedHash, $exprHash, $matrixHash, $clusterHash, $dhHash, $cnsHash, $chipHash, $genie3Hash, $knnHash);
}elsif($mode eq 'predict'){
	$exprHash = read_feature_network($exprData);
	$matrixHash = read_feature_network($matrixData);
	$clusterHash = read_feature_network($clusterData);
	$dhHash = read_feature_network($dhData);
	$cnsHash = read_feature_network($cnsData);
	$chipHash = read_feature_network($chipData);
	$genie3Hash = read_feature_network($genie3Data);
	$knnHash = read_feature_network($knnData);
	&make_test_matrix($networkHash, $exprHash, $matrixHash, $clusterHash, $dhHash, $cnsHash, $chipHash, $genie3Hash, $knnHash);
}

sub read_network {
	my $input = shift;
	my $network_hash={};
	my $tf={};
	my $tg={};
	open(R, $input);
	while(<R>){
		chomp;
		my @r = split(/\t/, $_);
		$network_hash->{$r[0]}->{$r[1]}++;
		$tf->{$r[0]}++;
		$tg->{$r[1]}++;
	}
	return($network_hash, $tf, $tg);
}

sub combine_interactions {
	my $network_hash=shift;
	my $tf=shift;
	my $tg=shift;
	my $combinedHash = {};
	foreach my $trans (keys %{$tf}){
		foreach my $target (keys %{$tg}){
			if(exists $network_hash->{$trans}->{$target}){
				$combinedHash->{$trans}->{$target}=1;
			}else{
				$combinedHash->{$trans}->{$target}=0;
			}
		}
	}
	return($combinedHash);
}

sub read_feature_network {
	my $input = shift;
	my $network_hash={};
	open(R, $input);
	while(<R>){
		chomp;
		my @r = split(/\t/, $_);
		$network_hash->{$r[0]}->{$r[1]}=$r[2];
	}
	return($network_hash);
}

sub make_train_matrix{
	my $combinedHash = shift;
	my $exprHash = shift; 
	my $matrixHash = shift;
	my $clusterHash=shift;
	my $dhHash = shift;
	my $cnsHash = shift;
	my $chipHash = shift;
	my $genie3Hash = shift;
	my $knnHash = shift;
	open(D, ">learning_matrix.txt");
	print D "Interaction\tCOE\tPWM\tCLUSTER\tDH\tCNS\tChIP\tGENIE3\tKNN\tClass\n";
	foreach my $tf (keys %{$combinedHash}){
		foreach my $tg (keys %{$combinedHash->{$tf}}){
			print D $tf."_".$tg."\t";
			if (exists $exprHash->{$tf}->{$tg}){
				print D $exprHash->{$tf}->{$tg}."\t";
			}else{
				print D "0\t";
			}
			if (exists $matrixHash->{$tf}->{$tg}){
				print D $matrixHash->{$tf}->{$tg}."\t";
			}else{
				print D "0\t";
			}
			if (exists $clusterHash->{$tf}->{$tg}){
				print D $clusterHash->{$tf}->{$tg}."\t";
			}else{
				print D "0\t";
			}
			if (exists $dhHash->{$tf}->{$tg}){
				print D $dhHash->{$tf}->{$tg}."\t";
			}else{
				print D "0\t";
			}
			if (exists $cnsHash->{$tf}->{$tg}){
				print D $cnsHash->{$tf}->{$tg}."\t";
			}else{
				print D "0\t";
			}
			if (exists $chipHash->{$tf}->{$tg}){
				print D $chipHash->{$tf}->{$tg}."\t";
			}else{
				print D "0\t";
			}
			if (exists $genie3Hash->{$tf}->{$tg}){
				print D $genie3Hash->{$tf}->{$tg}."\t";
			}else{
				print D "0\t";
			}
			if (exists $knnHash->{$tf}->{$tg}){
				print D $knnHash->{$tf}->{$tg}."\t";
			}else{
				print D "0\t";
			}
			print D $combinedHash->{$tf}->{$tg}."\n";
		}
	}
	
}

sub make_test_matrix{
	my $combinedHash = shift;
	my $exprHash = shift; 
	my $matrixHash = shift;
	my $clusterHash=shift;
	my $dhHash = shift;
	my $cnsHash = shift;
	my $chipHash = shift;
	my $genie3Hash = shift;
	my $knnHash = shift;
	open(D, ">test_matrix.txt");
	print D "Interaction\tCOE\tPWM\tCLUSTER\tDH\tCNS\tChIP\tGENIE3\tKNN\n";
	foreach my $tf (keys %{$combinedHash}){
		foreach my $tg (keys %{$combinedHash->{$tf}}){
			print D $tf."_".$tg."\t";
			if (exists $exprHash->{$tf}->{$tg}){
				print D $exprHash->{$tf}->{$tg}."\t";
			}else{
				print D "0\t";
			}
			if (exists $matrixHash->{$tf}->{$tg}){
				print D $matrixHash->{$tf}->{$tg}."\t";
			}else{
				print D "0\t";
			}
			if (exists $clusterHash->{$tf}->{$tg}){
				print D $clusterHash->{$tf}->{$tg}."\t";
			}else{
				print D "0\t";
			}
			if (exists $dhHash->{$tf}->{$tg}){
				print D $dhHash->{$tf}->{$tg}."\t";
			}else{
				print D "0\t";
			}
			if (exists $cnsHash->{$tf}->{$tg}){
				print D $cnsHash->{$tf}->{$tg}."\t";
			}else{
				print D "0\t";
			}
			if (exists $chipHash->{$tf}->{$tg}){
				print D $chipHash->{$tf}->{$tg}."\t";
			}else{
				print D "0\t";
			}
			if (exists $genie3Hash->{$tf}->{$tg}){
				print D $genie3Hash->{$tf}->{$tg}."\t";
			}else{
				print D "0\t";
			}
			if (exists $knnHash->{$tf}->{$tg}){
				print D $knnHash->{$tf}->{$tg}."\n";
			}else{
				print D "0\n";
			}
		}
	}
	
}