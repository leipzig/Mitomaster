package Mitoprocess::Blast;
use Bio::Tools::Run::StandAloneBlast;
use Bio::AlignIO;
use Mitoprocess::MitoprocessConfig;
use LWP::UserAgent;
use HTTP::Request::Common;
use HTML::Entities;

require Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw(get_blast_report_local get_blast_report_remote);

#use local blast (old blast, not blast+) server
sub get_blast_report_local {
    my ( $key, $seqfastaptr ) = @_;

    #ref stuff
    my $refseqDataHR = getRefseqData();
    my $stringfhref  = new IO::String( $refseqDataHR->{refseqseq} );
    my $format       = "fasta";
    my $seqio_inref =
      Bio::SeqIO->new( -fh => $stringfhref, -format => $format );
    my $refseqblast = $seqio_inref->next_seq;

    #creating a sequence object for the query
    my $fastaseq = ">$key\n" . $seqfastaptr->{$key};
    my $stringfh = new IO::String($fastaseq);
    my $seqio_in = Bio::SeqIO->new( -fh => $stringfh, -format => $format );
    my $seqio    = $seqio_in->next_seq;

    #running bl2seq
    my @params  = getParams();
    my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);
    my $report  = $factory->bl2seq( $refseqblast, $seqio );
    return $report;
}

#use remote server running

sub get_blast_report_remote {
    my ( $key, $seqfastaptr ) = @_;
    my $userAgent = LWP::UserAgent->new(timeout => 30); 
    #        Content_Type => 'multipart/form-data',
    my $request = POST BLAST_SERVER,
        Content => [ name => encode_entities($key), seq => encode_entities($seqfastaptr->{$key})];

    my $response = $userAgent->request($request);
    print $response->error_as_HTML . "
    " if $response->is_error;

    if ($response->is_success) {
         my $tmpFileName = UPLOADDIR . "/tmp_alignments" . time;
         #saving the file
         open( UPLOADFILE, ">$tmpFileName" ) or die "$tmpFileName $!";
         binmode UPLOADFILE;
         print UPLOADFILE $response->decoded_content;
         close UPLOADFILE;
         #$str = Bio::AlignIO->new(-file=> $tmpFileName,-format => 'bl2seq');
         my $str = new Bio::SearchIO(-format => 'blast', 
                                    -file   => $tmpFileName);
         return $str;
    } else {
         die $response->status_line;
    }
}
