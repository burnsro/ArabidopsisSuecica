#!/usr/bin/perl
#parameter: <inputfile>

use warnings;
use strict;
local $| = 1;

use POSIX qw/floor/;


if($#ARGV < 0)
{
	print "input parameter missing"."\n";
	exit;
}


my $input_filename = $ARGV[0];
my $input_file;
my $input_line;
my $input_linesn;
my @input_line_split;
my $input_lines_splitn;


my @splitdata;
my $splitdatan;

my $joineddatas;

if($input_filename !~ m/\.sam$/)
{
	print "wrong input format";
	exit;
}

my $temps = $input_filename;
$temps =~ s/\.bam\.sam$//;
$temps = $temps.".fastq";


my $output_filename = $temps;
my $output_file;

my $i;
my $j;

my $state = 0;
my $check = 0;
my $lastId = "";
my (@firstRead, @secondRead);
$firstRead[0] = "";
$secondRead[0] = "";
my ($firstReadFlag, $secondReadFlag);

open($input_file,'<'.$input_filename)  || die "Can't open".$input_filename;
open($output_file,'>'.$output_filename)  || die "Can't open".$output_filename;
$j=0;
while(<$input_file>)
{
	$input_line = $_;

	if($input_line =~ /^@/)
	{	next;	}

	if(length($input_line) > 5)
	{
		chomp($input_line);
		@splitdata = split("\t", $input_line);
		$splitdatan = @splitdata;
	
		@input_line_split = @splitdata;


		if(floor($input_line_split[1] / 0x40) % 2 > 0)
		{	$firstReadFlag = 1;	}
		else
		{	$firstReadFlag = 0;	}
		
		if(floor($input_line_split[1] / 0x80) % 2 > 0)
		{	$secondReadFlag = 1;	}
		else
		{	$secondReadFlag = 0;	}	



		if($state == 1 && not($input_line_split[0] eq $lastId))
		{
			print "no second read at line ".$input_line."\n";

			$state = 0;
		}
		
		if($state == 0)
		{
			if(!$firstReadFlag && !$secondReadFlag)
			{
				print "error at ".$input_line."\n";
				next;
			}
			
			if($firstReadFlag)
			{	@firstRead = @input_line_split;	}
			else
			{	@secondRead = @input_line_split;	}

			#print $firstReadFlag."\n\n";


			$state = 1;
		}
		else
		{
			if(!$firstReadFlag && !$secondReadFlag)
			{
				print "error at ".$input_line."\n";
				next;
			}
			
			if($firstReadFlag && length($firstRead[0]) > 1 || $secondReadFlag && length($secondRead[0]) > 1)
			{
				print "wrong flag at ".$input_line."\n";	
			}
			if($firstReadFlag)
			{	@firstRead = @input_line_split;	}
			else
			{	@secondRead = @input_line_split;	}

			#ready to output, plus checking
			#$check = 1;
			#print $firstRead[0]."\n";
			#print $secondRead[0]."\n";
			#if(length($firstRead[9]) != 100 || length($firstRead[10]) != 100 ||
			#	length($secondRead[9]) != 100 || length($secondRead[10]) != 100)
			#{
			#	print "sequence length error at ".$input_line."\n";
			#	$check = 0;
			#}
			
			$temps = "@".$firstRead[0]."\n";
			print($output_file $temps);
			
			$temps = $firstRead[9].$secondRead[9]."\n";
			print($output_file $temps);
			
			print($output_file "+\n");
			
			$temps = $firstRead[10].$secondRead[10]."\n";
			print($output_file $temps);
			
			$firstRead[0] = "";
			$secondRead[0] = "";
			$state = 0;
		}

		$lastId = $input_line_split[0];
		
	}

	$j++;
}
$input_lines_splitn = $j;
close($input_file);
close($output_file);

if($state != 0)
{	print "unfinished file\n";	}


