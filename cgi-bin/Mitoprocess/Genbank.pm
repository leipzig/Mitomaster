package Mitoprocess::Genbank;

use Carp;
use strict;

use Bio::DB::EUtilities;
use Mitoprocess::MitoprocessConfig;

require Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw(FULL_LENGTH ANY_LENGTH print_javascript get_fastas);

use constant FULL_LENGTH => { 'min' => '014000', 'max' => '017000' };
use constant ANY_LENGTH  => { 'min' => '000100', 'max' => '017000' };

# Synopsis:
# use Mitoprocess::Genbank;
# use warnings;
# use strict;
#
#
# my $gb = Mitoprocess::Genbank->new(length_range_href => FULL_LENGTH, numresults => 10);
#
# foreach my $gi(keys %{$gb->{'mitoseqs'}}){
#     print $gi."\t".$gb->{'mitoseqs'}{$gi}{'accession'}."\n";
# }

#to test from cli
#perl -e 'use Mitoprocess::Genbank; $gb=Mitoprocess::Genbank->new(numresults=>20);$gb->print_javascript()

#initialize with a range and numresults (FULL_LENGTH and 20000 by default )
sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self  = {
        length_range_href => FULL_LENGTH,
        numresults        => 20000,
        @_,    #override previous arguments
    };

    bless( $self, $class );

    if ( !defined $self->{query} ) {
        $self->{query} = _get_mt_query($self);
    }
    $self->{mitoseqs} = {};
    $self->_load_mitoseqs($self);

    return $self;
}

#return a query string to be used as a 'term'
sub _get_mt_query {
    my $self = shift;

#the original genbank query from mitomap
#(014000[SLEN] : 017000[SLEN]) AND "Homo"[Organism] AND mitochondrion[FILT] NOT (("Homo sp. Altai"[Organism] OR "Denisova hominin"[ALL]) OR "neanderthalensis"[ALL])

    #let's build this in a more structured manner
    my @good_orgs = ('Homo');
    my @bad_orgs  = ("Homo sp. Altai");
    my @bad_alls  = ( "Denisova hominin", "neanderthalensis" );
    my @filts     = ('mitochondrion');

    my $mt_query = '('
      . ( $self->{length_range_href}{'min'} )
      . '[SLEN] : '
      . ( $self->{length_range_href}{'max'} )
      . '[SLEN])';
    foreach my $good_org (@good_orgs) {
        $mt_query .= " AND \"$good_org\"[ORGANISM]";
    }

    foreach my $bad_org (@bad_orgs) {
        $mt_query .= " NOT \"$bad_org\"[ORGANISM]";
    }
    foreach my $bad_all (@bad_alls) { $mt_query .= " NOT \"$bad_all\"[ALL]" }
    foreach my $filt    (@filts)    { $mt_query .= " AND \"$filt\"[FILT]" }
    return $mt_query;
}

#given genbank query term, return array ptr of eutil ids (GIs)
sub _load_mitoseqs {
    my $self       = shift;
    my $batch_size = 200;                        #GET limit problems over 300
    my $eutil      = Bio::DB::EUtilities->new(
        -eutil  => 'esearch',
        -term   => $self->{query},
        -db     => 'nucleotide',
        -retmax => $self->{numresults},
        -email  => CONTACT_EMAIL
    );                                           # please use your real email

    # eutil => any of esearch, esummary, elink
    my @ids = $eutil->get_ids();    # returns array or array ref of IDs

    #we have to do these in batches otherwise it is too slow
    for ( my $id_cnt = 0 ; $id_cnt < scalar(@ids) ; $id_cnt += $batch_size ) {

        #don't slice beyond max index
        my $upper_limit =
          $id_cnt + $batch_size > scalar(@ids)
          ? scalar(@ids) - 1
          : $id_cnt + $batch_size - 1;
        my @gis_slice = @ids[ $id_cnt .. $upper_limit ];
        my $accs_ptr = _get_mt_acc( $self, \@gis_slice );
        for (
            my $gis_cnt = $id_cnt, my $acc_num = 0 ;
            $acc_num < scalar( @{$accs_ptr} ) ;
            $acc_num++, $gis_cnt++
          )
        {
            $self->{mitoseqs}{ $ids[$gis_cnt] }{'accession'} =
              $accs_ptr->[$acc_num];
        }
    }

}

#given GI, return accession
#straight from
#http://www.bioperl.org/wiki/HOWTO:EUtilities_Cookbook#Get_accessions_.28actually_accession.versions.29_for_a_list_of_GenBank_IDs_.28GIs.29
sub _get_mt_acc {
    my $self    = shift;
    my $gis     = shift;
    my $factory = Bio::DB::EUtilities->new(
        -eutil   => 'efetch',
        -db      => 'nucleotide',
        -id      => $gis,
        -email   => CONTACT_EMAIL,
        -rettype => 'acc'
    );
    my @accs = split( m{\n}, $factory->get_Response->content );
    return \@accs;
}

#return javascript formatted arrays for autocomplete section of mitomaster
sub print_javascript {
    my $self = shift;
    print "mitoids=[";
    foreach my $gi ( keys( %{ $self->{mitoseqs} } ) ) {
        print '"' . $gi . '","' . $self->{mitoseqs}{$gi}{'accession'} . '",';
    }
    print '""]' . "\n";
}

sub get_fastas {
    my $self    = shift;
    my @gis     = keys %{ $self->{mitoseqs} };
    if(DEBUG){print join ',',@gis; print "\n";}
    my $factory = Bio::DB::EUtilities->new(
        -eutil   => 'efetch',
        -db      => 'nucleotide',
        -id      => \@gis,
        -email   => CONTACT_EMAIL,
        -rettype => 'fasta'
    );
    if(DEBUG){print $factory->get_Response->content."\n";}
    return $factory->get_Response->content;
}
1;
