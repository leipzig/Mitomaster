package Mitoprocess::Parsing;
require Exporter;

use Array::Diff;
use Array::Utils qw(:all);
use Data::ArrayList;
use Mitoprocess::MitoprocessConfig;
use Carp;
use strict;
use warnings;
use Mitoprocess::SNV;

our @ISA = qw(Exporter);
our @EXPORT =
  qw(loadFasta loadSNVs acceptHSP arraytoString generate_random_string);

sub loadFasta {

    my ($filepath) = @_;
    open( my $IN, "$filepath" );
    my %seqs;
    my $defline;
    while (<$IN>) {
        if ( /^#/ || /^\n/ ) {

            #comment line or empty line do nothing
        }
        elsif (/>/) {
            chomp($_);
            s/>//g;

            #s/\s//g;
            s/[\/\\]//g;    #these will mess up file paths
            $defline = $_;
            $seqs{$defline} = "";
        }
        elsif ($defline) {
            chomp($_);
            s/[^a-zA-Z]//g;
            $_ = uc $_;
            s/U/T/g;
            $seqs{$defline} .= $_;
        }
        else {
            carp('Empty or invalid FASTA file');
            return ();
        }
    }

    #jerm 11/8 return only seqs>MIN_SEQ_LENGTH
    my @keepers = grep { length( $seqs{$_} ) >= MIN_SEQ_LENGTH } keys %seqs;
    my %keepSeqs;
    @keepSeqs{@keepers} = @seqs{@keepers};
    return \%keepSeqs;
}

sub loadSNVs {
    my ($filepath, $query) = @_;
    open( my $IN, "$filepath" );
    my $snvs = Mitoprocess::SNV->new($query);
    while (<$IN>) {
        chomp;
        if ( /^\#/ || /^\n/ || /sample\s+pos\s+ref\s+var/ ) {

            #comment line or empty line or header do nothing
        }
        else{
            $snvs->add($_);
        }
    }
    return $snvs->{snvs};
}

# this subroutin checks if the new HSP location is within any other previously stored HSP
# it takes the following parameteres: newstart,newend,oldHSP which is a hash containin already stored start and ends

sub acceptHSP {

    my $newStart = $_[0];
    my $newEnd   = $_[1];
    my $oldHSP   = $_[2];

    my $hashSize = scalar keys %$oldHSP;

    #print "hash size $hashSize <br>";
    if ( $hashSize < 1 ) {
        return 1;
    }

    #print "newstart $newStart, newEnd $newEnd <br>";

    foreach my $start ( keys %$oldHSP ) {
        my $end = $oldHSP->{$start};

        #print "stored hsp start: $start, end:$end<br>";
        if ( $newStart > $start and $newEnd < $oldHSP->{$start} ) {
            return 0;
        }
    }
    return 1;
}

sub generate_random_string {
    my $length_of_randomstring = shift;    # the length of
                                           # the random string to generate

    my @chars = ( 'a' .. 'z', 'A' .. 'Z', '0' .. '9', '_' );
    my $random_string;
    foreach ( 1 .. $length_of_randomstring ) {

        # rand @chars will generate a random
        # number between 0 and scalar @chars
        $random_string .= $chars[ rand @chars ];
    }
    return $random_string;
}

1;
