#!/usr/bin/perl
use CGI;
use CGI::Session;
use CGI::Carp qw(fatalsToBrowser);

#use strict;
use warnings;
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

# BLAST
$ENV{BLASTDIR}='/usr/local/bin/blast-2.2.25/bin/';

#######

# Process form if submitted; otherwise display it

if ( $query->param("submit") )
{    
    my $snvs=process_form($query);
	if(keys %{$snvs}){
	    print %{$snvs}."<br/>\n" if (DEBUG);
	    my $species_list;
	    if($session->param(-name=>'species')){
            $species_list='\''.join('\',\'',@{$session->param(-name=>'species')}).'\'';
        }
	    
	    
	    my @happtr=$query->param('haplotype');
        $session->param(-name=>'haplotype', -value=>\@happtr);
        $session->flush();

        
        
	    print_summary_preload_html('snvs','');

        my $run_id=process_seqs( $snvs , 'file' , 'snvlist', 'web', $species_list );
        
        if ( $query->param(-name => 'error_message') ) {
            my @errors = $query->param( -name => 'error_message');
            print "<div class=\"alert alert-error\">
                    <button type=\"button\" class=\"close\" data-dismiss=\"alert\">&times;</button>
                    <strong>Error:</strong> "
              . join('<br/>',@errors) . "</div>";
        }
        
        print_summary_html('snvs',$run_id);
    }else{
        $query->append(-name=>'error_message',-values=>["No valid SNV entries found.<br/>\n"]);
        
        display_snvs_form($query);
    }
} elsif($query->param("runid")){
        print_summary_preload_html('snvs',$query->param("runid"));
        print_summary_html('snvs',$query->param("runid"));
    }else{
	display_snvs_form($query);
}