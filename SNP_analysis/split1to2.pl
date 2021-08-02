#!/usr/bin/perl
#parameter: <inputfilename> <outputfilename1> <outputfilename2> <length>

use warnings;
use strict;
local $| = 1;

use Switch;

if(@ARGV<3){
       print "split one fastq file into two files with paired-end reads.
              Usage: <in> <out1> <out2> <rd_len>\n";
      exit(1);
}

open (IN, "<$ARGV[0]") || die "Can't open input file\n";
open (Rd1,">$ARGV[1]");
open (Rd2,">$ARGV[2]");

my $len=$ARGV[3];

#in which line we are in the 4 recurring..
my $state = 0;
my $i = 0;


while(my $line=<IN>)
{
	chomp $line;
	
	
	switch($state)
	{
		case 0
			{	if($line=~/^@/ && length($line)!=$len*2)
				{	#starting line, but the quality string can also start with that
					#if length == len * 2 it is likely no ID..
					print Rd1 $line,"/1\n";
					print Rd2 $line,"/2\n";
					$state = 1;
				}			
			}
		case 1
			{	if(length($line) == $len*2)
				{	print Rd1 substr($line,0,$len),"\n";
					print Rd2 substr($line,$len,$len),"\n";
					$state = 2;
				}
				else
				{	print "Rd_length wrong! (data line $i) \n";

					print Rd1 "\n";
					print Rd2 "\n";
					
					$state = 0;
				}
			}
		case 2
			{	if($line=~/^\+/)
				{	print Rd1 "\+\n";
					print Rd2 "\+\n";
					$state = 3;
				}
				else
				{	print "Format problem (data line $i)\n";
				
					print Rd1 "\n";
					print Rd2 "\n";
				
					$state = 0;
				}			
			}
		case 3
			{	if(length($line)==$len*2)
				{	print Rd1 substr($line,0,$len),"\n";
					print Rd2 substr($line,$len,$len),"\n";					
				}
				else
				{	print Rd1 "\n";
					print Rd2 "\n";			
					
					print "Q_length wrong! (data line $i)\n"; 
				}
				
				$state = 0;
			}	
	}
	
	$i++;
}

close(Rd2);
close(Rd1);
close(IN);


