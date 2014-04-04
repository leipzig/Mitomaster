#!/usr/bin/perl
use Array::Diff;
use Array::Utils qw(:all);
use Data::ArrayList;


use Bio::Mitomaster::SpeciesRef;
use Bio::Mitomaster::AASeq;
use Bio::Mitomaster::RNASeq;
my $ref = Bio::Mitomaster::SpeciesRef->new(species=>'human', reference=>'rCRS');

use Mitoprocess::Parsing;
use Mitoprocess::FormHandling;
use Mitoprocess::Results;
use Mitoprocess::MitoprocessConfig;
use Mitoprocess::Web;
#use strict;
use warnings;
use Bio::Tools::Run::StandAloneBlast;
use IO::String;
use DBI;
use File::Basename;
use List::Util;

use Benchmark;

# BLAST
$ENV{BLASTDIR}='/usr/local/bin/blast-2.2.25/bin/';

#######
print $ARGV[0]."\n" if (DEBUG);

my $start = new Benchmark;

my $seqfasta = loadFasta($ARGV[0]);

my $skipto;
if($ARGV[1]){
    $skipto=$ARGV[1];
}

print_report_header();
#http://www.webquills.net/web-development/perl/perl-5-hash-slices-can-replace.html
#let's iterate over the hashref
foreach $def(sort keys %{$seqfasta}){
	my $hashslice={$def=>$seqfasta->{$def}};
	#>gi|383403118|gb|JF824858.2| Homo sapiens isolate C136 mitochondrion, complete genome
	#gi|320339343|gb|HQ713447.1| Homo sapiens isolate GC-5 mitochondrion, complete genome
	if($def =~ m/gi\|\d+\|\S+\|(\S+)\|/){
	    if($skipto){
		if($skipto eq $1){$skipto=0;}
		else{next;}
	    }
	    if(genbank_record_exists($1)){
		print STDERR "seen $1\n";
		next;
	    }
	}else{
	    print STDERR "unknown format $def\n";
	}
	print STDERR $def."\n";
	my $run_id=process_seqs( $hashslice , 'file' , 'sequences', 'cli', undef  );
	my @query_ids=Mitoprocess::Results::get_query_ids($run_id );
    
    foreach my $query_id(@query_ids){   
	print "here is query $query_id\n" if (DEBUG);
	print_var_report_txt($query_id);
    }
}
   
my $end = new Benchmark;

my $elapsed = timediff ($end, $start);
print STDERR "Elapsed time: ", timestr ($elapsed), "\n";
