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

my $gb = Mitoprocess::Genbank->new( query => $ARGV[0] );
my $sequences = $gb->get_fastas();

my $tmpFileName = UPLOADDIR . "/tmp_genbank_ids" . time;

#saving the file
open( UPLOADFILE, ">$tmpFileName" ) or die "cannot find $tmpFileName";
binmode UPLOADFILE;
print UPLOADFILE $sequences;
close UPLOADFILE;

#read the fasta file in a hash
my $entries = loadFasta($tmpFileName);
my $run_id=process_seqs( $entries , 'file' , 'sequences', 'cli', undef  );
my @query_ids=Mitoprocess::Results::get_query_ids($run_id );
print_report_header();
foreach my $query_id(@query_ids){
   print "here is query $query_id\n" if (DEBUG);
   print_var_report_txt($query_id);
}
   
my $end = new Benchmark;

my $elapsed = timediff ($end, $start);
print STDERR "Elapsed time: ", timestr ($elapsed), "\n";