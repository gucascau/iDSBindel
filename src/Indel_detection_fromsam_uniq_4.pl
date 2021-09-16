#!/usr/bin/perl
use strict;
use warnings;
my $version="1.0 version";
use Getopt::Long;

# Copyright (c) 2021 Dr. Kaifu lab    
# Author:Xin Wang 
# email: xin.wang@childrens.harvard.edu
# PI: Kaifu Chen

# Function: This scripts is to devided the mapped sequence into small insertions (more than 1bp) and small deletions (more than 1bp), considering the highly error read mapping

my %opts;
GetOptions(\%opts,"i:s","o:s","m:s","n:s","c:s","a:s","b:s","d:s","e:s","h:s","f:s","r:s","l:s","q:s","u:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{o} || !defined $opts{q} || defined $opts{h} ) {
	die "************************************************************************
	Usage: $0.pl -i Mapping results -q Sample ID -o output folder [Option]
	
			Request Parameters:
				-i Mapping results with SAM format (mapped from bwa)
				-o Output files that divided reads into small insertion, small deletion and no others
				-q Sample ID
			
			Optional Parameters:
				-m Mata size (default 84)
				-n Short insertion size limit (default 10)
				-c Minimum read counts to support the insertion or deletion(default 5)
				-a Mata mapping chromosome (default ChrIII)
				-b Mata mapping start (default 294300)
				-e Mata mapping end (default 294500)
				-f The size of left primer, extended 5bp if there are deletion on the primer (default 25)
				-r The size of right primer, extended 5bp if there are deletion on the primer (default 22)
				-l The size of left side of Mata region, here we allowed 4 nucleotide shift (default 45)
				-d The size of right side of Mata region, here we allowed 4 nucleotide shift (default 51)
				-u The collect upstream and downstream for the unique of deletion or insertions from the raw reads (default 5bp)
				-h Help
				
************************************************************************\n";
}

my $output=$opts{o};
my $input=$opts{i};
my $Smata=(defined $opts{m})?$opts{m}:84;
my $Sindel=(defined $opts{n})?$opts{n}:10;
my $MinRc=(defined $opts{c})?$opts{c}:5;
my $MataChr=(defined $opts{a})?$opts{a}:"ChrIII";
my $MataStart=(defined $opts{b})?$opts{b}:294300;
my $MataEnd=(defined $opts{e})?$opts{e}:294500;

my $Lprimer=(defined $opts{f})?$opts{f}:25;
my $Rprimer=(defined $opts{r})?$opts{r}:22;


#### Here we should allow 2bp frameshift of Mapping.
my $MataSleft =(defined $opts{l})?$opts{l}:45;
my $MataSright=(defined $opts{d})?$opts{d}:51;

my $Mcutedge=(defined $opts{u})?$opts{u}:5;

my $SampleID=$opts{q};

my $Maxsize=$Smata+$Sindel;
#### input of sam file ####

my %deletion; my %insertion; my %normal; my %hash; my %st; my %num;

open I,"$input" or die "cannot open file $input";

mkdir $output;
open OUTP,">$output/$output.noinsertion.txt" or die $!;
open UN,">$output/$output.undetermined.txt" or die $!;

open CONFUSE,">$output/$output.confusingindel.txt" or die $!;

my $Cnoinsertion=0;
my $Totalnum=0;

open T,">$output/$output.test.txt" or die $!;
while (<I>){
	chomp;
	s/\r//g;
	next if (/^\@/);

	my ($i,$chr,$start,$match,$string)=(split /\t/,$_)[0,2,3,5,9];
	my $length=length $string;
	
	
	##### Reads filter ####
	#### remove the events that have the longer insertion, this might contain some cases that forward read and reverese reads are not consistent.
	if ($length>$Maxsize){
		print UN "$_\n";
		next;
	}

	next unless ($string =~ m/^GCA/ && $string =~m /TTG$/);
	
	### remove reads have showed up this is to remove the duplicates
	next if (exists $hash{$i});
	
	##### remove reads that cannot map to this region
	next unless ($chr eq $MataChr && $start>$MataStart && $start <=$MataEnd);
	$hash{$i}++;
	$Totalnum++;

	##### Matches to identify insertion and deletion #####
	
	#### Change the matches that (\d+)M(\d+)S to ($1+$2)M, and (\d+)S(\d+)M to ($1+$2)M
	
	#$match=~s/(\d+)M(\d+)(I|D)(\d+)M(\d+)S/($1)M($2)($3)($4+$5)M/g;
	#$match=~s/(\d+)S(\d+)M(\d+)(I|D)(\d+)M/($1+$2)M($3)($4)($5)M/g;
	
	$match=~s{(\d+)M(\d+)S}{($1+$2).'M'}eg;
	
	$match=~s{(\d+)S(\d+)M}{($1+$2).'M'}eg;
	
	#$match=~/(\d+)H//g;
	
	### In the following we ignore the mismatches
	
	if ($match =~m/^(\d+)M$/ || $match =~m/^(\d+)M(\d+)S$/ || $match =~m/^(\d+)S(\d+)M$/){
		
		### print out perfect matches that did not show any indel events:
		print OUTP "$_\n";
		$Cnoinsertion++;
		
	}elsif($match =~m/^(\d+)M(\d+)(I|D)(\d+)M$/){
		
		### the deletion cannot exceed to the primer: set up the primer size as cut of of Matches. 
		if ($1>=$Lprimer && $1<=$MataSleft && $4 >=$Rprimer && $4<=$MataSright ){
			
			my $part1=$1-5;
			### short Deletion events:
			if($3 eq "D"){
				$deletion{$2}.="$_\n";
				### part 2 show the insertion/deletion size.
				my $part2=$2; 
				
				### part3 is the right side starts, deletion will be not added
				my $part3=$part1+$Mcutedge;
				
				### the left side 5bp string
				my $string1=substr($string,$part1,$Mcutedge);
				
				#### the right side 5bp string
				my $string2=substr($string,$part3,$Mcutedge);
				
				### join the left 5bp string, "D", deletion size, "D", and the right 5bp string
				my $str=join '',($string1,"D",$part2,"D",$string2);
				
				### Check the read counts of this deletion.
				$num{$str}++;
				
				### identify the best of the string;
				$st{$str}=$_;
				print T "$str\t$_\n";
				
			### Short insertion events:	
			}else{
				$insertion{$2}.="$_\n";
				my $part2=$1; my $part3=$part1+$Mcutedge+$2;
			
				my $string1=substr($string,$part1,$Mcutedge);
				my $string3=substr($string,$part3,$Mcutedge);
				
				### the inner insertion sequences
				my $string2=substr($string, $part2,$2);
				my $str=join '',($string1,$2,"I",$string2,"I",$string3);
				$num{$str}++;
				
				$st{$str}=$_ ;		
				print T "$str\t$_\n";
			}

		}else{

			print OUTP "$_\n";
			$Cnoinsertion++;
		}
	
	}elsif($match =~m/^(\d+)M(\d+)(I|D)(\d+)M(\d+)(I|D)(\d+)M$/){
		
		###### because of certain indels in the upstream sequences
		###### here is the format :15M1I25M8I35M
		
		## 38M2D28M1I16M 		
		### sequence errors in the right part:
 	#$1>=$Lprimer && $1<=$MataSleft && $4 >=$Rprimer && $4<=$MataSright
		if ($1>=$Lprimer && $1<=$MataSleft && ($4+$7) >=$Rprimer && ($4+$7)<=$MataSright){
			
			### set up the start site
			my $part1=$1-$Mcutedge; 
			
			if ($3 eq "I"){
				$insertion{$2}.="$_\n" ;
				my $part2=$part1+$Mcutedge; my $part3=$part1+$2+$Mcutedge;
				
				my $string1=substr($string,$part1,$Mcutedge);
				my $string3=substr($string,$part3,$Mcutedge);
				my $string2=substr($string, $part2,$2);
				
				my $str=join '',($string1,$2,"I",$string2,"I",$string3);
				$num{$str}++;
				
				$st{$str}=$_ if (!exists $st{$str});
				
				#$st{$str}=$_ ;
			
				print T "$str\t$_\n";
				
			###### here is the format : 3S37M2D19M1D27M3S	
			}else{
				$deletion{$2}.="$_\n";
				my $part2=$2; my $part3=$part1+$Mcutedge;
			
				my $string1=substr($string,$part1,$Mcutedge);
				my $string2=substr($string,$part3,$Mcutedge);
				#my $str=join "$part2"."D", ($string1,$string2);
				my $str=join '',($string1,"D",$part2,"D",$string2);
				$num{$str}++;
				$st{$str}=$_ if (!exists $st{$str});
				print T "$str\t$_\n";		
			}
			
		### sequence errors in the left part:			
		}elsif(($1+$4)>=$Lprimer && ($1+$4)<= $MataSleft && $7 >=$Rprimer && $7 <=$MataSright){
			
			### set up the start site
			my $part1=($3 eq "D")? ($1+$4-$Mcutedge):($1+$2+$4-$Mcutedge);
			
			if ($6 eq "I"){
				$insertion{$5}.="$_\n" ;

				my $part2=$part1+$Mcutedge;
				my $part3=$part1+$5+$Mcutedge;
				
				my $string1=substr($string,$part1,$Mcutedge);
				my $string3=substr($string,$part3,$Mcutedge);
				my $string2=substr($string, $part2,$5);
			
				my $str=join '',($string1,$5,"I",$string2,"I",$string3);
			
				$num{$str}++;	
				$st{$str}=$_ if (!exists $st{$str});

				print T "$str\t$_\n";
			###### because of certain indels in the upstream sequences
			###### here is the format : 3S15M1D25M8D35M3S
				
			}else{
			 	$deletion{$5}.="$_\n";
				
				my $part2=$5;
				my $part3=$part1+$Mcutedge;
				

				my $string1=substr($string,$part1,$Mcutedge);
				my $string2=substr($string,$part3,$Mcutedge);
				
				#my $str=join "$6"."D", ($string1,$string2);
				my $str=join '',($string1,"D",$part2,"D",$string2);

				$num{$str}++;
				#$dellength{$str}=$3 if (!exists $st{$str});
				$st{$str}=$_ if (!exists $st{$str});
				print T "$str\t$_\n";			
			 }
			 
			 ### the others we considered them as no insertion or deletion events, but the indel caused by sequence errors.			 	
		}else{
			print OUTP "$_\n";
			$Cnoinsertion++;
		}
		
	### This is mainly due to the low quality of the reads, so we also try to fix this problem caused by sequence errors.		
	}elsif($match =~m/^(\d+)M(\d+)(I|D)(\d+)M(\d+)(I|D)(\d+)M(\d+)(I|D)(\d+)M$/){

		### sequence errors in the left part:	
		if($1>=$Lprimer && $1 <=$MataSleft && ($4+$7+$10)>=$Rprimer  && ($4+$7+$10) <=$MataSright){
			
			### sequence error on the right part, we just ignore them. 
			### set up the start site
			my $part1=$1-5;
			
			if($3 eq "I"){
				$insertion{$2}.="$_\n";
			
				my $part2=$part1+$Mcutedge; my $part3=$part1+$2 +$Mcutedge;
			
				my $string1=substr($string,$part1,$Mcutedge);
				my $string3=substr($string,$part3,$Mcutedge);
				my $string2=substr($string, $part2,$2);
			
				my $str=join '',($string1,$2,"I",$string2,"I",$string3);
			
				$num{$str}++;
				#$inslength{$str}=$3 ;
				#$inssequence{$str}=$string2;
			
				#### here collect the best alignment of the sequence
				$st{$str}=$_ if (!exists $st{$str});	
				print T "$str\t$_\n";				
			}else{
				$deletion{$2}.="$_\n";
				my $part2=$2; my $part3=$part1+$Mcutedge;
			
				my $string1=substr($string,$part1,$Mcutedge);
				my $string2=substr($string,$part3,$Mcutedge);
				#my $str=join "$part2"."D", ($string1,$string2);
				my $str=join '',($string1,"D",$part2,"D",$string2);
				$num{$str}++;
				$st{$str}=$_ if (!exists $st{$str});
				print T "$str\t$_\n";				
			}
			
		## sequence error in the right part 
		}elsif(($1+$4+$7)>=$Lprimer && ($1+$4+$7)<=$MataSleft && $10>=$Rprimer  && $10 <=$MataSright){
			
			#### Set up the start site considering the sequence errors.
			my $part1;	
			if($3 eq "I" && $6 eq "D"){
				$part1=$1+$4+$7+$2-$Mcutedge;
			}elsif($3 eq "D" && $6 eq "I"){
				$part1=$1+$4+$7+$5-$Mcutedge;
								
			}elsif($3 eq "I" && $6 eq "I"){
				$part1=$1+$4+$7+$2+$5-$Mcutedge;
			}else{	
				$part1=($1+$4+$7-$Mcutedge);		
			}
			
			if($9 eq "I"){
				$insertion{$8}.="$_\n";
			

				my $part2=$part1+$Mcutedge; my $part3=$part1 +$8+$Mcutedge;
			
				my $string1=substr($string,$part1,$Mcutedge);
				my $string3=substr($string,$part3,$Mcutedge);
				my $string2=substr($string, $part2,$8);
			
				my $str=join '',($string1,$2,"I",$string2,"I",$string3);
			
				$num{$str}++;
				#$inslength{$str}=$3 ;
				#$inssequence{$str}=$string2;
			
				#### here collect the best alignment of the sequence
				$st{$str}=$_ if (!exists $st{$str});
				print T "$str\t$_\n";				
	
			}else{
				$deletion{$8}.="$_\n";
				my $part2=$8; my $part3=$part1 +$Mcutedge;
			
				my $string1=substr($string,$part1,$Mcutedge);
				my $string2=substr($string,$part3,$Mcutedge);
				#my $str=join "$part2"."D", ($string1,$string2);
				my $str=join '',($string1,"D",$part2,"D",$string2);
				$num{$str}++;
				$st{$str}=$_ if (!exists $st{$str});
				print T "$str\t$_\n";				

			}
			## sequence error in the both sides of MAT
		}elsif(($1+$4)>=$Lprimer && ($1+$4)<=$MataSleft && ($7+$10)>=$Rprimer  && ($7+$10) <=$MataSright){
			
			
			# 3M1D34M2I36M1D9M
			#### Set up the start site considering the sequence errors.
			my $part1=($3 eq "D")? ($1+$4-$Mcutedge):($1+$2+$4-$Mcutedge);
			
			if($6 eq "I"){
				$insertion{$5}.="$_\n";
			

				my $part2=$part1+$Mcutedge; my $part3=$part1 +$5+$Mcutedge;
			
				my $string1=substr($string,$part1,$Mcutedge);
				my $string3=substr($string,$part3,$Mcutedge);
				my $string2=substr($string, $part2,$5);
			
				my $str=join '',($string1,$2,"I",$string2,"I",$string3);
			
				$num{$str}++;
				#$inslength{$str}=$3 ;
				#$inssequence{$str}=$string2;
			
				#### here collect the best alignment of the sequence
				$st{$str}=$_ if (!exists $st{$str});
				print T "$str\t$_\n";				
	
			}else{
				$deletion{$5}.="$_\n";
				my $part2=$5; my $part3=$part1 +$Mcutedge;
			
				my $string1=substr($string,$part1,$Mcutedge);
				my $string2=substr($string,$part3,$Mcutedge);
				#my $str=join "$part2"."D", ($string1,$string2);
				my $str=join '',($string1,"D",$part2,"D",$string2);
				$num{$str}++;
				$st{$str}=$_ if (!exists $st{$str});
				print T "$str\t$_\n";				

			}
			
		}else{
			## sequence error in the Middle part: this part would be really had to define whether it is the sequence errors or the changes in the junction region,
			## due to the lower number of this events, we will directly print them out.
			print CONFUSE "$_\n";
			
		}
	}else{
		print CONFUSE "$_\n";
	}
	
	# elsif($match =~m/^(\d+)M(\d+)(I|D)(\d+)M(\d+)S$/){
#
# 		my $part1=$1-$Mcutedge;
#
# 		if($1>=$Lprimer && $1 <=$MataSleft && ($4+$5)>=$Rprimer  && ($4+$5) <=$MataSright){
# 			if ($3 eq "I"){
#
# 			$insertion{$2}.="$_\n";
# 				my $part2=$1; my $part3=$part1+$Mcutedge+$2;
#
# 				my $string1=substr($string,$part1,$Mcutedge);
# 				my $string3=substr($string,$part3,$Mcutedge);
#
# 			### the inner insertion sequences
# 				my $string2=substr($string, $part2,$2);
# 				my $str=join '',($string1,$2,"I",$string2,"I",$string3);
# 				$num{$str}++;
#
# 				$st{$str}=$_  if (!exists $st{$str});
# 				print T "$str\t$_\n";
#
#
# 			}else{
#
# 				$deletion{$2}.="$_\n";
# 			### part 2 show the insertion/deletion size.
# 				my $part2=$2;
#
# 			### part3 is the right side starts, deletion will be not added
# 				my $part3=$part1+$Mcutedge;
#
# 			### the left side 5bp string
# 				my $string1=substr($string,$part1,$Mcutedge);
#
# 			#### the right side 5bp string
# 				my $string2=substr($string,$part3,$Mcutedge);
#
# 			### join the left 5bp string, "D", deletion size, "D", and the right 5bp string
# 				my $str=join '',($string1,"D",$part2,"D",$string2);
#
# 			### Check the read counts of this deletion.
# 				$num{$str}++;
#
# 			### identify the best of the string;
# 				$st{$str}=$_ if (!exists $st{$str});
# 				print T "$str\t$_\n";
#
# 			}
# 		}else{
# 			print OUTP "$_\n";
# 		}
#
# 	}# elsif($match =~m/^(\d+)S(\d+)M(\d+)(I|D)(\d+)M$/){
# 		my $part1=$1+$2-$Mcutedge;
# 		if(($1+$2)>=$Lprimer && ($1+$2) <=$MataSleft && $5>=$Rprimer  && $5 <=$MataSright){
# 			if ($4 eq "I"){
#
# 				$insertion{$3}.="$_\n";
#
# 				my $part2=$3; my $part3=$part1+$Mcutedge+$3;
#
# 				my $string1=substr($string,$part1,$Mcutedge);
# 				my $string3=substr($string,$part3,$Mcutedge);
#
# 			### the inner insertion sequences
# 				my $string2=substr($string, $part2,$3);
# 				my $str=join '',($string1,$3,"I",$string2,"I",$string3);
# 				$num{$str}++;
#
# 				$st{$str}=$_  if (!exists $st{$str});
# 				print T "$str\t$_\n";
#
# 			}else{
#
# 				$deletion{$3}.="$_\n";
# 			### part 2 show the insertion/deletion size.
# 				my $part2=$3;
#
# 			### part3 is the right side starts, deletion will be not added
# 				my $part3=$part1+$Mcutedge;
#
# 			### the left side 5bp string
# 				my $string1=substr($string,$part1,$Mcutedge);
#
# 			#### the right side 5bp string
# 				my $string2=substr($string,$part3,$Mcutedge);
#
# 			### join the left 5bp string, "D", deletion size, "D", and the right 5bp string
# 				my $str=join '',($string1,"D",$part2,"D",$string2);
#
# 			### Check the read counts of this deletion.
# 				$num{$str}++;
#
# 			### identify the best of the string;
# 				$st{$str}=$_ if (!exists $st{$str});
# 				print T "$str\t$_\n";
#
# 			}
# 		}else{
# 			print OUTP "$_\n";
# 		}
# 	}elsif($match =~m/^(\d+)S(\d+)M(\d+)(I|D)(\d+)M(\d+)S$/){
#
# 	}else{
# 		print CONFUSE "$_\n";
# 	}
}		
close ;


open INS,">$output/$output.insertion.stat" or die $!;
open DEL,">$output/$output.deletion.stat" or die $!;
open INSF,">$output/$output.insertion.fasta" or die $!;
open DELF,">$output/$output.deletion.fasta" or die $!;
open STAT,">$output/$output.finalstatitic.txt" or die $!;



my $Cinsertion=0; my $Cdeletion=0;
foreach my $i (sort {$num{$b} <=> $num{$a}} keys %num){
	
	my ($id,$index,$sequence)=(split/\t/,$st{$i})[0,5,9];
	
	if ($i =~/([ATGC]+)(\d+)I/){
		print INS "$i\t$num{$i}\t$2\t$st{$i}\n";
		
		print INSF ">$id\t$index\t$i\t$num{$i}\t$2\n$sequence\n";
		
		$Cinsertion +=$num{$i};
		
		
	}elsif($i=~/([ATGC]+)D(\d+)D/){
		
		print DEL "$i\t$num{$i}\t$2\t$st{$i}\n";
		print DELF ">$id\t$index\t$i\t$num{$i}\t$2\n$sequence\n";
		$Cdeletion +=$num{$i};
		
	}
	
}
my $Cundefined=$Totalnum-$Cnoinsertion-$Cinsertion-$Cdeletion;

print STAT "ReadInfo:\t$SampleID\nTotalReads:\t$Totalnum\nNoextra:\t$Cnoinsertion\nShortInsertion:\t$Cinsertion\nShortDeletion:\t$Cdeletion\nUnclear:\t$Cundefined\n";

close INS;
close DEL;
close STAT;
close INSF;
close DELF;

