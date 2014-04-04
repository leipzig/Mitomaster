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
print $query->header();

# BLAST
$ENV{BLASTDIR}='/usr/local/bin/blast-2.2.25/bin/';

#######

# Process form if submitted; otherwise display it

if ( defined $query->param("runid") )
{
    if ( defined $query->param("queryid") ){
        print_details_html($query->param("refp"),$query->param("queryid"),$query->param('runid') );
    }else{
        print_details_html($query->param("refp"),undef,$query->param('runid') );
    }
}else
{
	display_sequence_form($query);
}