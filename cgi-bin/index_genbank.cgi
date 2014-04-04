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
    my $ids=process_form($query);
	if(keys %{$ids}){
	    print %{$ids}."<br/>\n" if (DEBUG);
	    my $species_list;
	    if($session->param(-name=>'species')){
            $species_list='\''.join('\',\'',@{$session->param(-name=>'species')}).'\'';
        }
	    
	    #assume haplogrouping on
	    my @happtr=('true');
        $session->param(-name=>'haplotype', -value=>\@happtr);
        $session->flush();
        
        
	    print_summary_preload_html('genbank','');
        my $run_id=process_seqs( $ids , 'file' , 'sequences', 'web', $species_list );
        print_summary_html('genbank',$run_id);
    }else{
        $query->append(-name=>'error_message',-values=>"No valid genbank entries found<br/>\n");
        display_genbank_form($query);
    }
} else
{
	display_genbank_form($query);
}