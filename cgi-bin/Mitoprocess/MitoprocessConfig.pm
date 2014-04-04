package Mitoprocess::MitoprocessConfig;
require Exporter;
use DBI;
use strict;
use warnings;
use Carp;
use Bio::Mitomaster::SpeciesRef;
use Bio::Mitomaster::Seq;
use Mitoprocess::LocalSettings;
our @ISA = qw(Exporter);
our @EXPORT =
  qw(HUMAN_NAME CONS_FLANK HAPLOGROUP_PROGRAM HAPLOGREP_SERVER USE_REMOTE_HAPLOGREP MIN_HAPLOTYPING_COVERAGE GET_LOCLIST_FROM_SNV_RANGE PHYLOTREE_VERSION CONTACT_EMAIL BLAST_SERVER USE_REMOTE_BLAST getParams getRefseqData getDSN $rCRS_SR $rCRS_seq $rCRS SITE_NAME GENOMETOOLS_HOME TEMP_PICS TEMP_PICS_PATH REFERENCE_SEQ DB_PATH DB_USER DB_PASS MIN_SEQ_LENGTH UPLOADDIR GFFDIR DATALIFE  DEBUG USE_RESULTMV WEB_HOME MITOMAP_HOME CGI_SCRIPTALIAS USE_MARIE_STYLE_INSERTIONS SHOW_3107 MAX_NONSENSE ACCORDION_TOGGLE);

#CONTACTS
use constant CONTACT_EMAIL => 'leipzigj@email.chop.edu';

# CONSTANTS DEFINITION
use constant REFERENCE_SEQ  => "rCRS";
use constant MIN_SEQ_LENGTH => 100;
#use constant UPLOADDIR      => "../uploads";

use constant DATALIFE =>
  10;    # Minimum lifetime of stored data in database (minutes)

    #"/var/www/mitomasterbeta/h-mito/h-mito.py"
use constant DEBUG => 0;

use constant HAPLOGROUP_PROGRAM => 'HaploGrep';
use constant PHYLOTREE_VERSION => 16;
# Updating the resultmv which is a dynamic table for improving performance of data retrieval
# it can also result in DBD::Pg::st execute failed: ERROR:  cache lookup failed for relation 27474
#CONTEXT:  SQL statement "DROP TABLE IF EXISTS mitomasterb.resultmv"
#PL/pgSQL function "updateresults" line 4 at SQL statement at Mitoprocess/FormHandling.pm line 176.
use constant USE_RESULTMV => 0;

use constant USE_MARIE_STYLE_INSERTIONS => 1; #C->CT instead of :->T
use constant SHOW_3107 => 0;
use constant MAX_NONSENSE               => 7; #roughly how many amino acids we are willing to show
our $rCRS_SR =
  Bio::Mitomaster::SpeciesRef->new( species => 'human', reference => 'rCRS' );
our $rCRS_seq =
  Bio::Mitomaster::Seq->new( species_ref => $rCRS_SR, variants => {} );
our $rCRS = $rCRS_seq->seq();

#does the species selection opening of an element in turn close the others?
use constant ACCORDION_TOGGLE =>0;

#conservation window flank
use constant CONS_FLANK => 5;

use constant HUMAN_NAME => 'Homo_sapiens_(rCRS)';

use constant BLAST_SERVER => 'http://rescommap02.research.chop.edu:8080/';
use constant USE_REMOTE_BLAST => 1;

use constant HAPLOGREP_SERVER => 'http://rescommap02.research.chop.edu:8081/';
use constant USE_REMOTE_HAPLOGREP => 0;
use constant MIN_HAPLOTYPING_COVERAGE => 100;

#do we use the list of snvs in a genotype or just assume its covered?
use constant GET_LOCLIST_FROM_SNV_RANGE => 1;

sub getDSN {

    #connect to postgres database
    our $dsn = DBI->connect( DB_PATH, DB_USER, DB_PASS );
    if ( !defined $dsn ) {

        croak "Cannot connect to database!</br>";

    }
    return $dsn;
}

sub getParams {

    #
    # DEFINING BLAST PARAMETERS USED LATER BY BL2SEQ
    #
    our @params = (
        'program' => 'blastn',
        '-F'      => 'F',
        '-e'      => '1',
        '-r'      => '1',
        '-E'      => '0',
        '-G'      => '0',
        '-X'      => '50',
        '-W'      => '28',
        '-m'      => 'T',
        '-q'      => '-2',

    );

    return @params;
}

sub getRefseqData {

    #Querying the reference sequence from database
    our $refseqName = REFERENCE_SEQ;
    our $refseqseq  = '';
    our $refseqid   = '';
    our $sequence   = '';
    my $selquery = "SELECT id,sequence FROM mitomasterb.refseq WHERE title=?";
    my $dsn      = getDSN();
    my $query_handle = $dsn->prepare($selquery);
    $query_handle->execute($refseqName);

    $query_handle->bind_columns( undef, \$refseqid, \$sequence );

    # LOOP THROUGH RESULTS
    if ( $query_handle->fetch() ) {
        $refseqseq = ">refSeq\n" . $sequence;
    }
    else {
        print "Reference sequence not found in the database <br />";
        die;
    }
    $query_handle->finish();
    $dsn->disconnect();
    my %hash1 = (
        "refseqid"  => $refseqid,
        "refseqseq" => $refseqseq
    );
    return \%hash1;
}
1;
