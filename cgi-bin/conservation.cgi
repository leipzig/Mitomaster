#!/usr/bin/perl
use CGI;
use CGI::Session;
use CGI::Carp qw(fatalsToBrowser);

#use strict;
use warnings;

use Mitoprocess::Conservation;

# Create the CGI object
my $query = new CGI;

#create session and save id as a cookie
my $sid = $query->cookie("CGISESSID") || undef;
my $session = new CGI::Session(undef, $sid, {Directory=>'/tmp'});
my $cookie = $query->cookie(CGISESSID => $session->id);
print $query->header( -cookie=>$cookie );

# Process form if submitted; otherwise display it

if ( defined $query->param("locus") && defined $query->param("pos"))
{
    if(defined $query->cookie("CGISESSID")){
        my $sid = $query->cookie("CGISESSID");
        my $session = new CGI::Session(undef, $sid, {Directory=>'/tmp'});
        if($session->param(-name=>'species')){
            my $specptr = $session->param(-name=>'species');
            print_conservation_window_as_javascript($query->param("rcrs"),$query->param("locus"),$query->param("pos"),$query->param("res"),$specptr);
        }else{
            print_conservation_window_as_javascript($query->param("rcrs"),$query->param("locus"),$query->param("pos"),$query->param("res"),undef);
        }
    }else{
        print_conservation_window_as_javascript($query->param("rcrs"),$query->param("locus"),$query->param("pos"),$query->param("res"),undef);
    }
}