package Mitoprocess::SNV;

use strict;
use Carp;
use Bio::Mitomaster::SpeciesRef;
use Bio::Mitomaster::Seq;

# Synopsis:
# my @muts = ('Sample1  73  A   G','Sample1 146 T   C','Sample2 263 A   :','Sample2 709 :   G','9028T','291d','291del','573insC','573.1C','573.2C','C9028T','A291d',
# 'A291del','A291-','A291:','C573CC','C573CCC','Sample1 A9028T','Sample2    291d',' Sample 3    291del','573insC','573.1C','573.2C','C9028T','A291d','A291del','A291-','A291:','C573CC','C573CCC');
# my %snvs;
# 
# my $snvtest = Mitoprocess::SNV->new();
# 
# foreach my $mut(@muts){
#     $snvtest->add($mut);
# }
# $snvtest->print_report();

my $rCRS_ref = Bio::Mitomaster::SpeciesRef->new(species=>'human',reference=>'rCRS');
my $rCRS_seq =  Bio::Mitomaster::Seq->new(species_ref=>$rCRS_ref,  variants=>{});
my $rCRS = $rCRS_seq->seq();

sub new {
    my ($proto,$query) = @_;
    my $class = ref($proto) || $proto;
    my $self  = {
        query => $query
    };

    bless( $self, $class );

    $self->{snvs} = {};

    return $self;
}


sub add {
    my ($self,$mut) = @_;
    $mut =~ s/^\s+//g;
    $mut =~ s/\s+$//g;
    my ($cleanRef,$cleanPos,$cleanSNV,$cleanSample,$rstate,$pstate,$sstate);
    if(my ($sample, $pos, $ref, $snv ) = $mut =~ /^(.+)\t(\d+)\t(\S)\t(\S+)$/){
             #print "typeI $mut\tsample $sample\tref $ref\tpos $pos\tsnv $snv\n";
             ($cleanRef,$rstate) = $self->_formatRef($pos,$ref);
             ($cleanPos,$pstate) = $self->_formatPos($pos);
             ($cleanSNV,$sstate) = $self->_formatSNV($snv);
             $cleanSample = $sample;
             if($self->_add($rstate,$pstate,$sstate,$cleanSample,$cleanRef,$cleanPos,$cleanSNV)){
                 return 1;
             }
    }
    elsif($mut =~ /^(.+)\t(\S+)$/){
        if(my ($sample, $ref, $pos, $snv ) = $mut =~ /^(.+)\t([ACGT:])?([0-9.]+)((ins)?[ACGT]+|-|:|d|del)$/){
            $cleanSample = $sample || "query";
            #print "typeII $mut\tsample $sample\tref $ref\tpos $pos\tsnv $snv\n";
             ($cleanRef,$rstate) = $self->_formatRef($pos,$ref);
             ($cleanPos,$pstate) = $self->_formatPos($pos);
             ($cleanSNV,$sstate) = $self->_formatSNV($snv);
             if($self->_add($rstate,$pstate,$sstate,$cleanSample,$cleanRef,$cleanPos,$cleanSNV)){
                 return 1;
             }
        }
    }elsif($mut =~ /^(\S+)$/){
        if (my ($ref, $pos, $snv ) = $mut =~ /^([ACGT:])?([0-9.]+)((ins)?[ACGT]+|-|:|d|del)$/){
            #print "typeIII $mut\tref $ref\tpos $pos\tsnv $snv\n";
            ($cleanRef,$rstate) = $self->_formatRef($pos,$ref);
            ($cleanPos,$pstate) = $self->_formatPos($pos);
            ($cleanSNV,$sstate) = $self->_formatSNV($snv);
            $cleanSample = "query";
            if($self->_add($rstate,$pstate,$sstate,$cleanSample,$cleanRef,$cleanPos,$cleanSNV)){
                return 1;
            }
        }
    }
    
    carp "$mut not parsable\n";
    if($self->{query}){
        $self->{query}->append(
            -name   => 'error_message',
            -values => ["There was an error with the input."]
        );
    }
    
    
    
}

sub _add {
    my ($self,$rstate,$pstate,$sstate,$cleanSample,$cleanRef,$cleanPos,$cleanSNV) = @_;
    
    if(!defined $cleanRef || !defined $cleanPos || !defined $cleanSNV){
        return 0;
    }
    
    #aggregate
    if($rstate eq 'insertion' || $pstate eq 'insertion' || $sstate eq 'insertion'){
        if(defined $self->{snvs}{$cleanSample}{$cleanPos}){
            $cleanSNV = $self->{snvs}{$cleanSample}{$cleanPos}[1].$cleanSNV;
        }
        #C573CC -> :573C
        if($cleanRef =~ /[ACTG]/ && length($cleanSNV)>1){$cleanSNV = substr($cleanSNV,1)}
        $cleanRef = ':';
        #bump position to be compatible with what sequence runs return
        #:573C becomes :574C
        $cleanPos++;
    }
    #print "{$cleanSample}{$cleanPos} = [ $cleanRef, $cleanSNV ]\n";
    
    
    $self->{snvs}{$cleanSample}{$cleanPos} = [ $cleanRef, $cleanSNV ];
    return 1;
}

sub print_report {
    my $self = shift;
    foreach my $snvsamp(keys %{$self->{snvs}}){
        foreach my $pos(keys %{$self->{snvs}{$snvsamp}}){
            print $snvsamp."\t".$pos."\t".$self->{snvs}{$snvsamp}{$pos}[0]."\t".$self->{snvs}{$snvsamp}{$pos}[1]."\n";
        }
    }
}
sub _formatRef{
    my ($self,$pos,$ref) = @_;
    if(!defined $ref){
        return (uc(substr($rCRS,$pos-1,1)),undef);
    }
    if($ref eq '-' || $ref eq ':'){
        return (':','insertion');
    }
    unless ($ref =~ /^[ACGTacgt]$/){
        carp("Reference base should be [ACGT] of length 1");
        if($self->{query}){
            $self->{query}->append(
                -name   => 'error_message',
                -values => ["Reference base should be [ACGT] of length 1"]
            );
        }
        #print $ref."\n".length($ref)."\n";
        return (undef,undef);
    }
    
    #print $rCRS."\n";
    if(uc(substr($rCRS,$pos-1,1)) ne uc($ref) && substr($rCRS,$pos-1,1)){
        carp("rCRS position ".$pos." should be ".uc(substr($rCRS,$pos-1,1)));
        if($self->{query}){
            $self->{query}->append(
                -name   => 'error_message',
                -values => ["rCRS position ".$pos." should be ".uc(substr($rCRS,$pos-1,1))]
            );
        }
        return (undef,undef)
    }
    
    return (uc($ref),undef);
}

sub _formatPos{
    my ($self,$posarg) = @_;
    my ($pos) = $posarg =~ m/^(\d+)/;
    if ($posarg =~ /\./){
        return ($pos,'insertion');
    }
    
    if ($pos<1 || $pos > length($rCRS)){
        carp("Position should be between 1 and ".length($rCRS));
        if($self->{query}){
            $self->{query}->append(
                -name   => 'error_message',
                -values => ["Position should be between 1 and ".length($rCRS)]
            );
        }
        return (undef,undef);
    }
    
    return ($pos,undef);
}

sub _formatSNV{
    my ($self,$snvarg) = @_;
    if($snvarg eq '-' || $snvarg eq ':' || $snvarg eq 'd' || $snvarg eq 'del'){
        return (':','deletion');
    }
    if ($snvarg =~ m/[Ii]/){
        my ($snv) = $snvarg =~ m/([ACGTacgt])$/; #insC
        return (uc($snv), 'insertion');
    }
    if (my ($snv) = $snvarg =~ m/^([ACTGactg]{2,})$/){
        return (uc($snv), 'insertion');
    }
    return (uc($snvarg),undef);
}