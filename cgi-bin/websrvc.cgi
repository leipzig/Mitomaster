#!/usr/bin/perl
use CGI;
use CGI::Carp qw(fatalsToBrowser);

#use strict;
use warnings;
use Bio::Tools::Run::StandAloneBlast;
use IO::String;

use File::Basename;
use List::Util;

use Mitoprocess::Parsing;
use Mitoprocess::FormHandling;
use Mitoprocess::Results;
use Mitoprocess::MitoprocessConfig;
use Mitoprocess::Web;
# Create the CGI object
my $query = new CGI;

# Output the HTTP header
#print $query->header();

# BLAST
$ENV{BLASTDIR}='/usr/local/bin/blast-2.2.25/bin/';

#######

# Process form if submitted; otherwise display it


if ( $query->param("file"))
{
    my $filetype = $query->param("fileType") || "sequences";
	my $seqfasta=process_post($query);
	my $run_id=process_seqs( $seqfasta , 'file' , $filetype , 'websrvc', undef);
    my @query_ids=get_query_ids($run_id );
    print_report_header();
    foreach my $query_id(@query_ids){
        print_var_report_txt($query_id);
    }
} else
{
    print $query->header( -status => 200 );
    print "<html>missing file or file argument<br/>\n</html>";
}
