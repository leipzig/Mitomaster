#!/usr/bin/perl
use CGI;
use CGI::Carp qw(fatalsToBrowser);

#use strict;
use warnings;

use Mitoprocess::Web;

# Create the CGI object
my $query = new CGI;

# Output the HTTP header
print $query->header();

if ( defined $query->param("locus") && defined $query->param("pos"))
{
    print_conservation_popup($query->param("rcrs"),$query->param("locus"),$query->param("pos"),$query->param("res"),$query->param("query"));
}else{
    print "missing parameters<br/>\n";
}