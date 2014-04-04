#!/usr/bin/perl
use CGI;
use CGI::Session;
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

#create session and save id as a cookie
my $sid = $query->cookie("CGISESSID") || undef;
my $session = new CGI::Session(undef, $sid, {Directory=>'/tmp'});
my $cookie = $query->cookie(CGISESSID => $session->id);
print $query->header( -cookie=>$cookie );
#print $query->header(-type => "application/json", -charset => "utf-8");

if ( $query->param("species") )
{
	open (LOGFILE, "> /tmp/cons.txt");
	print LOGFILE $query->param("species");
	close (LOGFILE);
    my @specptr=$query->param('species');
    $session->param(-name=>'species', -value=>\@specptr);
    $session->flush();
    # create a JSON string according to the database result
    #, "species" : scalar(@specptr)
    my $json =  qq{{"success" : "species saved"}};

    # return JSON string
    
    print $json;
    
}else{
    print_checkbox_popup($query);
}