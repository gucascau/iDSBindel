#!/usr/bin/perl
use strict;
use warnings;
my $version="1.0 version";
use Getopt::Long;

# Copyright (c) 2021 Dr. Kaifu lab    
# Author:Xin Wang 
# email: xin.wang@childrens.harvard.edu
# PI: Kaifu Chen

# Function: This script is to extract the high quality reads that did not contain the large insertion events. In the meantime, we eliminated the three index of PE.

my %opts;
GetOptions(\%opts,"i:s","o:s","f:s","r:s", "h:s","m:s","n:s");
print "*************\n*$version*\n*************\n";
if ( !defined $opts{i}|| !defined $opts{o} ||!defined $opts{r} ||!defined $opts{f} || defined $opts{h}) {
       	die "************************************************************************
       	Usage: $0.pl -i Linsertion -f forward reads -r reverse reads -o Output
				
		Request Parameters:
			-i: file with Read id that contains large insertion events
			-f: forward fastq1
			-r: reverse fastq2
			-o: Output fastq files that eliminated the three customed indexes.
			
		Optional Parameters:
			-m: Trimmed forward index size (Default, 3)
			-n: Trimmed reverese index size (Default, 3)
			-h Help	

************************************************************************\n";
}


my $input=$opts{i};
my $fastq1=$opts{f};
my $fastq2=$opts{r};
my $Ftrimsize=(defined $opts{t})?$opts{t}:3;
my $Rtrimsize=(defined $opts{t})?$opts{t}:3;
my $out=$opts{o};


open FASTQ1,"$fastq1" or die "cannot open file $fastq1";

my $id1; my %string1; my %hash1;
while (<FASTQ1>) {
	chomp;
    # print "$_" ;

    if ($. % 4 == 1)  {
		$id1=(split/\s+/,$_)[0];
		$id1=~s/^@//;
		#print "$id1\n";
		$hash1{$id1}=$_."\n";
	}elsif($.%4 == 2){
		#$hash1{$id1}=$_;
		my $string=substr($_,$Ftrimsize, -$Rtrimsize);
		$hash1{$id1}.=$string."\n";

	}elsif ($.%4 == 0){
		my $string=substr($_,$Ftrimsize, -$Rtrimsize);
		$hash1{$id1}.=$string."\n";	

	}else{
		$hash1{$id1}.=$_."\n";
	}

}
close FASTQ1;

my $id2; my %hash2; 
open FASTQ2, "$fastq2" or die "cannot open file $fastq2";

while (<FASTQ2>) {
	chomp;
    # print "$_" ;

    if ($. % 4 == 1)  {
		$id2=(split/\s+/,$_)[0];
		$id2=~s/^@//;
		$hash2{$id2}=$_."\n";
	}elsif($.%4 == 2){
		my $string=substr($_,$Rtrimsize, -$Ftrimsize);
		$hash2{$id2}.=$string."\n";
		
	}elsif ($.%4==0){
		my $string=substr($_,$Rtrimsize, -$Ftrimsize);
		$hash2{$id2}.=$string."\n";
	}else{
		$hash2{$id2}.=$_."\n";
	}
}
close FASTQ2;

# create a index of large insertion
open IN,"$input" or die $!;
my %insertion;
while (<IN>){
	chomp;
	my $id =(split /\t/,$_)[0];
	$insertion{$id}++;
	#$phix{$_}++;
}
close IN;



open PH1,">$out.R1.fastq" or die $!;
open PH2,">$out.R2.fastq" or die $!;
#open OUT,">$out.assembled.fasta" or die $!;

foreach my $i (keys %hash1){
	next if (exists $insertion{$i});
	next unless (exists $hash1{$i} && exists $hash2{$i});
	print PH1 "$hash1{$i}";
	print PH2 "$hash2{$i}";
	
}

close PH1;
close PH2;


