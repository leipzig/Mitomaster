#!/usr/bin/perl
use Mitoprocess::Parsing;
use Mitoprocess::FormHandling;
use Mitoprocess::Results;
use Mitoprocess::MitoprocessConfig;

use Benchmark;

# BLAST
$ENV{BLASTDIR}='/usr/local/bin/blast-2.2.25/bin/';

#######
print $ARGV[0]."\n" if (DEBUG);

my $start = new Benchmark;

my $snvptr = loadSNVs($ARGV[0]);
my $run_id=process_seqs( $snvptr , 'file' , 'snvlist' , 'cli', undef );
my @query_ids=get_query_ids($run_id );
print_report_header();
foreach my $query_id(@query_ids){
   print "here is query $query_id\n" if (DEBUG);
   print_var_report_txt($query_id);
}
   
my $end = new Benchmark;

my $elapsed = timediff ($end, $start);
print STDERR "Elapsed time: ", timestr ($elapsed), "\n";
