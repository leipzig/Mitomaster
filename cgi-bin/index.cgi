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

    
# BLAST
$ENV{BLASTDIR}='/usr/local/bin/blast-2.2.25/bin/';

#######

# Process form if submitted; otherwise display it

if ( $query->param("submit") )
{  
    my $seqfasta=process_form($query);
	if(keys %{$seqfasta}){
	    print %{$seqfasta}."<br/>\n" if (DEBUG);
	    my $species_list;
        if($session->param(-name=>'species')){
            $species_list='\''.join('\',\'',@{$session->param(-name=>'species')}).'\'';
        }
        
        #assume haplogrouping on
	    my @happtr=('true');
        $session->param(-name=>'haplotype', -value=>\@happtr);
        $session->flush();
        
        print_summary_preload_html('sequences','');
        my $run_id=process_seqs( $seqfasta , 'file' , 'sequences', 'web', $species_list);
        print_summary_html('sequences',$run_id);
    }else{
        $query->append(-name=>'error_message',-values=>"No valid sequences found<br/>\n");
        display_sequence_form($query);
    }
} elsif($query->param("runid")){
    print_summary_preload_html('sequences',$query->param("runid"));
    print_summary_html('sequences',$query->param("runid"));
}else
{
	display_sequence_form($query);
}