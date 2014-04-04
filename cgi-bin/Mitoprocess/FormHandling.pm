package Mitoprocess::FormHandling;
use Mitoprocess::MitoprocessConfig;
use Mitoprocess::Parsing;
use Mitoprocess::Haplotype;
use Mitoprocess::Blast;
use Bio::Mitomaster::SpeciesRef;
use Bio::Mitomaster::AASeq;
use Bio::Mitomaster::RNASeq;

use strict;
no strict "refs";
use warnings;


use IO::String;
use File::Basename;
use List::Util;
use CGI qw(:standard);

require Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw(process_seqs);

#process all seqs or snvs
#return db-issued run_id
sub process_seqs {
    my ( $seqfastaptr, $containerType, $fileType, $caller, $species_list ) = @_;
    my %seqfasta      = %{$seqfastaptr};
    my $totalseqcount = keys %seqfasta;

    my $dsn = getDSN();

    my @params = getParams();
    print "params:" . join( ' ', @params ) . "\n" if (DEBUG);
    my $blastparams = join( '.', @params );
    print "blastparams:" . $blastparams . "\n" if (DEBUG);

    my $run_id = _get_run_id( $seqfastaptr, $blastparams, $containerType );

    my $seqProcessedCounter = 0;
    my $progressRatio       = 0;
    my $sequenceReported    = 0;
    my $progress            = 0;

    #looping through uploaded sequences
    print "i have " . scalar( keys %seqfasta ) . " keys<br/>\n" if (DEBUG);
    foreach my $key ( keys %seqfasta ) {
        if ( $fileType eq 'sequences' ) {
            _process_one_seq( $key, $seqfastaptr, $run_id, $species_list );
        }
        elsif ( $fileType eq 'snvlist' ) {
            _process_one_snv_sample( $key, $seqfastaptr, $run_id,
                $species_list );
        }
        $seqProcessedCounter++;
        $progress =  ($seqProcessedCounter * 100 / $totalseqcount);

        #if (   ( ( $seqProcessedCounter - $sequenceReported ) > 9 )
        #    || ( ( $progress - $progressRatio ) > 29 ) )
        #{
            $sequenceReported = $seqProcessedCounter;
            $progressRatio    = $progress;
            #print $seqProcessedCounter. " sequences processed<br/>\n";
            #print "$progress%\n";
            if($caller eq 'web'){
                print "<script type=\"text/javascript\">\$(\"#progressbar\").progressbar('value', $progress);</script>\n";
        }
            
        #}
    }    # end of query sequences loop

    #update run datemodified command to when analysis finished
    my $updaterunfinish =
"update mitomasterb.run set datemodified = now() where mitomasterb.run.id=?";
    my $update_runobj = $dsn->prepare($updaterunfinish);
    $update_runobj->execute($run_id);
    $update_runobj->finish();

# Updating the resultmv which is a dynamic table for improving performance of data retrieval
    if (USE_RESULTMV) {
        my $updateresultcomm = "select * from mitomasterb.updateresults();";
        my $update_objv      = $dsn->prepare($updateresultcomm);
        $update_objv->execute();
        $update_objv->finish();
    }

    $dsn->disconnect();
    return $run_id;
}

#get a run_id issued from database
sub _get_run_id {
    my ( $seqfastaptr, $blastparams, $containerType ) = @_;
    my %seqfasta = %{$seqfastaptr};

    # Database connection is established
    my $dsn = getDSN();

    # cleaning the stored data from old runs
    my $lifetime = DATALIFE;
    my $cleancommand =
"delete from mitomasterb.run where (EXTRACT(day from now() - mitomasterb.run.datecreated) * 1440 + EXTRACT(hour from now() - mitomasterb.run.datecreated) * 60 + EXTRACT(minute from now() - mitomasterb.run.datecreated)) >? AND ((select count(*) from mitomasterb.variant)>20000) ";
    my $clean_obj = $dsn->prepare($cleancommand);
    $clean_obj->execute($lifetime);

#creating a sequence object for the reference from database to be used for bl2seq
    my $refseqDataHR = getRefseqData();
    my $stringfhref  = new IO::String( $refseqDataHR->{refseqseq} );

    my $run_id = '';
    if ( ( scalar( keys %seqfasta ) ) > 0 ) {

        #inserting run information to the database
        my $insertcommand =
"insert into mitomasterb.run(datecreated,note,\"refseqID\",\"type\",\"blastparams\") values (now(),'no description',$refseqDataHR->{refseqid},'$containerType','$blastparams') returning id;";

        print $insertcommand. "<br/>\n" if (DEBUG);
        my $insert_obj = $dsn->prepare($insertcommand);
        $insert_obj->execute();
        my $foo_rec = $insert_obj->fetchrow_hashref();
        $run_id = $$foo_rec{"id"};
        if(!$run_id){
            my $cgi = new CGI;
            print $cgi->header( -status => 204 );
            print "Unable to create run id. Please contact administrator.<br/>\n";
            exit(0);
        }
        print "Assigned run ID:$run_id<br/>\n" if (DEBUG);
        print "Inserted run id: $run_id<br/>\n" if (DEBUG);

    }
    else {
        my $cgi = new CGI;
        print $cgi->header( -status => 204 );
        print "No sequences.<br/>\n";
        exit(0);
    }

    $dsn->disconnect();
    return $run_id;
}

sub _new_query {
    my ( $key, $seqfastaptr, $run_id ) = @_;
    my $queryseq = $seqfastaptr->{$key};

    #trimming sequence label
    $key =~ m/^(\S+)/;
    my $subquerylabel = $1;

    my $dsn = getDSN();

    my $insertcommand =
"insert into mitomasterb.query(label,\"length\",pstart,pend,\"runID\",sequence) values ('$subquerylabel',"
      . ( length $queryseq )
      . ",0,0,$run_id, '') returning id;";
    print $insertcommand. "<br/>\n" if (DEBUG);
    my $insert_objq = $dsn->prepare($insertcommand);
    $insert_objq->execute();
    my $foo_recq = $insert_objq->fetchrow_hashref();
    my $query_id = $$foo_recq{"id"};
    $insert_objq->finish();

    print "inserted query details ID: $query_id <br/>" if (DEBUG);
    $dsn->disconnect();
    if(!$query_id){
        my $cgi = new CGI;
        print $cgi->header( -status => 204 );
        print "Unable to create query id. Please contact administrator.<br/>\n";
        exit(0);
    }
    return $query_id;
}



sub _process_one_snv_sample {
    my ( $key, $snvsptr, $run_id, $species_list ) = @_;
    my %samp            = %{ $snvsptr->{$key} };
    my $query_id        = _new_query( $key, $snvsptr, $run_id );
    my $delcount        = 0;
    my $inscount        = 0;
    my $variants        = '';
    my $min_var_pos     = 100000;
    my $max_var_pos     = 0;
    my $seq_region_type = "Full_genome";
    foreach my $varpos ( keys %samp ) {
        $min_var_pos = ( $varpos < $min_var_pos ) ? $varpos : $min_var_pos;
        $max_var_pos = ( $varpos > $max_var_pos ) ? $varpos : $max_var_pos;
        my $hit_ch = $samp{$varpos}->[1];
        my $ref_ch = $samp{$varpos}->[0];

        $variants .= $ref_ch.$varpos.$hit_ch.";";

        _insert_variation( $ref_ch, $hit_ch, \$delcount, \$inscount, $varpos, 0,
            $query_id, $species_list );

    }
    my $loclist="";
    if(GET_LOCLIST_FROM_SNV_RANGE){
        $loclist = $min_var_pos . "-" . $max_var_pos . ";";
    }else{
        $loclist = "1-16569"
    }
    $variants =~ s/;$//;
    _haplogroup( $variants, $query_id, $seq_region_type, $loclist );
}

sub _process_one_seq {
    my ( $key, $seqfastaptr, $run_id, $species_list ) = @_;

    my $query_id = _new_query( $key, $seqfastaptr, $run_id );
    
    my $seq_region_type = "Full_genome";

#$variants stores all variants by reference position in a string separated by ; to be used later for haplotyping
    my ( $variants, $loclist ) =
      _get_sequence_based_variants( $key, $seqfastaptr, $query_id,
        $seq_region_type, $species_list );
    print "variants: $variants loclist: $loclist\n" if (DEBUG);
    _haplogroup( $variants, $query_id, $seq_region_type, $loclist );
}

sub _get_sequence_based_variants {
    my ( $key, $seqfastaptr, $query_id, $seq_region_type, $species_list ) = @_;
    my $queryseq = $seqfastaptr->{$key};
    my $dsn      = getDSN();
    my $report = (USE_REMOTE_BLAST ? get_blast_report_remote( $key, $seqfastaptr ) : get_blast_report_local( $key, $seqfastaptr ));

# refalignpos is a hash that stores reference start and end positions of aligned segments returned from blast and stored in database
    my %refalignpos;

#qalignpos is a hash that stores reference start and end positions of aligned segments returned from blast and stored in database
    my %qalignpos;

    #iterating through the blast results
    my $result_count = 0;
    my $variants     = "";
    my $loclist;
    if ( my $result = $report->next_result ) {
        $result_count++;

        # print "result: $result_count <br/>" if (DEBUG);
        my $hit_count = 0;
        if ( my $hit = $result->next_hit ) {
            $hit_count++;
            my $hsp_count = 0;

            # print "hit: $hit_count <br/>\n"
            while ( my $hsp = $hit->next_hsp ) {
                $hsp_count++;

                # print "hsp: $hsp_count <br/>\n" if(DEBUG);
                # $queryname = $result->query_name;
                # print "Query: " . $queryname . "<br/>\n" if(DEBUG);

                #start and end position of current HSP on the reference sequence
                my $refstart = $hsp->start('query');
                my $refend   = $hsp->end('query');

                #start and end position of current HSP on the query sequence
                my $qstart = $hsp->start('hit');
                my $qend   = $hsp->end('hit');

# Here we check if current HSP is a within another HSP that has been seen before. If so, we just ignore it.
# Since the previous HSP has had a better score.
#
                if ( !acceptHSP( $refstart, $refend, \%refalignpos ) ) {
                    next;
                }
                if ( !acceptHSP( $qstart, $qend, \%qalignpos ) ) {
                    next;
                }
                ######

                $refalignpos{$refstart} = $refend;

                print "stored in hash refend: $refstart =>" if (DEBUG);
                print $refalignpos{$refstart}               if (DEBUG);
                print "<br/>\n"                             if (DEBUG);

                $qalignpos{$qstart} = $qend;

                print "stored in hash qend: $qstart => $qend <br/>\n"
                  if (DEBUG);

                printf( "%7s: %s<br/>\n", "start", $refstart ) if (DEBUG);
                printf( "%7s: %s<br/>\n", "end",   $refend )   if (DEBUG);
                printf( "%7s: %s<br/>\n",
                    "query", length( $hsp->query_string ) )
                  if (DEBUG);
                printf( "%7s: %s<br/>\n", "hom", $hsp->homology_string )
                  if (DEBUG);
                printf( "%7s: %s<br/>\n", "hit", $hsp->length('total') )
                  if (DEBUG);
                printf( "%7s: %s<br/>\n",
                    "percent_identity", $hsp->percent_identity )
                  if (DEBUG);

                #blast e-value
                my $evalue = $hsp->evalue();

                print "evalue: $evalue <br/>\n" if (DEBUG);

                #inserting HSP information to the query table
                $loclist .= $refstart . "-" . $refend .";";
                #not sure I want this.";";
                my $insertcommand =
"insert into mitomasterb.aligndetails(\"queryID\",\"refstart\",\"refend\",\"qstart\",\"qend\",\"evalue\") values ($query_id,\'$refstart\', \'$refend\',\'$qstart\', \'$qend\',\'$evalue\')  returning id;";

                print "$insertcommand <br/>\n" if (DEBUG);

                my $insert_objv = $dsn->prepare($insertcommand);
                $insert_objv->execute();
                my $foo_recv = $insert_objv->fetchrow_hashref();
                $insert_objv->finish();
                my $aligndetail_id = $$foo_recv{"id"};

                if ( length($queryseq) < 6000 ) {

# INCOMPLETE: This should become smarter to detect "HVS1", "Control_region", "Full_genome".
# This information should be passed to h-mito for haplogrouping. Check h-mito.py for more the reasoning behind this.
                    if ( $refend < 600 ) {
                        $seq_region_type = "Control_region";
                    }
                }

                my $delcount           = 0;
                my $inscount           = 0;
                my $refpolymer         = '';    #multideletion
                my $hitpolymer         = '';    #mutliinsertion
                my $insintohomopolymer = ''
                  ; #homopolymeric reference area indel (e.g. C->CC), credit last
                    #looping through HSP sequence position by position
                for (
                    my $key = 0 ;
                    $key < length( $hsp->query_string ) ;
                    $key++
                  )
                {
                    my $ref_ch =
                      uc( $refpolymer . substr( $hsp->query_string, $key, 1 ) );

                    my $hit_ch =
                      uc( $hitpolymer . substr( $hsp->hit_string, $key, 1 ) );

                    if ($insintohomopolymer) {

#assignthe reference a -, assign the query that reference residue (since it is a homopolymer we know that is also the query residue)
                        $ref_ch = '-';
                        $hit_ch = $refpolymer
                          . uc( substr( $hsp->query_string, $key, 1 ) );
                    }
                    
                    $refpolymer         = '';
                    $hitpolymer         = '';
                    $insintohomopolymer = '';
                    if ( $ref_ch ne $hit_ch ) {

                        # A variant detected!
                        #notice the inscount is supposed to shift varpos to the proper position
                        my $varpos = $key + $hsp->start('query') - $inscount;
                        
                        

                       #finding the corresponding position on the query sequence
                        my $queryvarpos = $hsp->start('hit') + $key - $delcount;

                        print "ref position: ", $varpos, " query position:",
                          $queryvarpos, " query:", $hit_ch, " ref:",
                          $ref_ch, "next ref",
                          uc( substr( $hsp->query_string, ( $key + 1 ), 1 ) ),
                          " <br/>\n" if (DEBUG);
                          print $inscount."\n" if (DEBUG);

                        if ( $ref_ch eq '-' || $hit_ch eq '-' ) {

                            #if a polymer indel hold off
                            unless ( $key == length( $hsp->query_string ) + 1 )
                            {
                                if (
                                    $ref_ch eq '-'
                                    && uc(
                                        substr(
                                            $hsp->query_string, ( $key + 1 ),
                                            1
                                        )
                                    ) eq '-'
                                  )
                                {
                                    print "insertion\n" if (DEBUG);
                                    $inscount++;
                         #hold off and assign this query del to next position
                                    $hitpolymer = $hit_ch;
                                }
                                elsif (
                                    $hit_ch eq '-'
                                    && uc(
                                        substr(
                                            $hsp->hit_string, ( $key + 1 ),
                                            1
                                        )
                                    ) eq '-'
                                  )
                                {
                         #hold off and assign this query insert to next position
                            print " deletion\n" if (DEBUG);
                                        $refpolymer = $ref_ch;
                                }
                                elsif (
                                    $ref_ch eq '-'
                                    && $hit_ch eq uc(
                                        substr(
                                            $hsp->query_string, ( $key + 1 ),
                                            1
                                        )
                                    )
                                  )
                                {
                                    print "homopolymer\n"  if (DEBUG);
                                    #homopolymer - shift this guy one over
                                    $insintohomopolymer = 1;
                                }
                            }#not a polymer
                        }
                        unless ( $hitpolymer
                            || $refpolymer
                            || $insintohomopolymer )
                        {
                            _insert_variation(
                                $ref_ch,    $hit_ch, \$delcount,
                                \$inscount, $varpos, $queryvarpos,
                                $query_id,  $species_list
                            );
                            #Adding the variant position to $variants list
                           
                             $variants .= $ref_ch.$varpos.$hit_ch.";";
                            
                        }

                    }    #end of variant detected if statement
                } # end variantz loop i.e. looping through HSP string position by position
            }    # next hsp
        }
    }

    print "inputvariants"
      . $variants
      . "query id $query_id seqregion $seq_region_type\n"
      if (DEBUG);
    $dsn->disconnect();
    $variants =~ s/;$//;
    return ( $variants, $loclist );
}



sub _insert_variation {
    my ( $ref_ch, $hit_ch, $delcountref, $inscountref, $varpos, $queryvarpos,
        $query_id, $species_list )
      = @_;
    my $dsn = getDSN();


    #mitomaster locus types are t,m,n,r
    #trna, protein, non-coding, ribosomalRNA
    my $locus_sql =
"SELECT locus_id, common_name,type FROM mitomasterb.locus WHERE (starting <= ending AND starting<=? AND ending>=?) OR (ending < starting AND (? <=ending OR ? >=starting)) ORDER BY locus_id;";
    print $locus_sql. "\n" if (DEBUG);
    my $query_handle = $dsn->prepare($locus_sql);
    $query_handle->execute($varpos,$varpos,$varpos,$varpos);
    my (
        $common_name,   $locus_type, $assigned_locus,
        $coding_effect, $conservation, $conservation_pos,
        $conservation_res, $variant_res, $locus_id
    );
    $query_handle->bind_columns( undef, \$locus_id, \$common_name, \$locus_type );

    if ( $query_handle->fetch() ) {
        $assigned_locus = $common_name;

        if ( $locus_type eq 'n' ) {
            $coding_effect = "non-coding";
            $conservation  = undef;
            $conservation_pos = undef;
            $conservation_res = undef;
            $variant_res = undef;
        }
        else {
            ( $coding_effect, $conservation , $conservation_pos, $conservation_res, $variant_res ) =
              _get_coding_effect( $locus_type, $varpos, $ref_ch, $hit_ch,
                $species_list );
        }
        #if this variation falls into more than one region, call it
        while ( $query_handle->fetch() ) {
            $assigned_locus.="/$common_name";
        }
    }
    else {
        $assigned_locus = "unknown";
        $coding_effect  = "non-coding";
        $conservation   = undef;
        $conservation_pos = undef;
        $conservation_res = undef;
        $variant_res = undef;
    }
    $query_handle->finish();

    #inserting sequence information to the query table

    my ( $ref_ch_disp, $hit_ch_disp, $varpos_disp, $queryvarpos_disp, $note ) =
      get_insertion_format( $ref_ch, $hit_ch, $delcountref, $inscountref,
        $varpos, $queryvarpos );
        
        #if($rCRS_seq->seq($varpos_disp,$varpos_disp) ne $ref_ch_disp){print "$varpos should be ".$rCRS_seq->seq($varpos_disp,$varpos_disp)." not ".$ref_ch_disp."\n";}
    
    unless($note =~ m/ignored/){
        my $variantinsertcommand =
    "insert into mitomasterb.variant(position,\"refNA\",\"queryNA\",\"queryID\",\"note\",\"qposition\",\"locus\",\"aa_change\",\"conservation\",\"conservation_pos\",\"variant_res\",\"locus_id\") values (?,?,?,?,?,?,?,?,?,?,?,?) returning id;";

        print "$variantinsertcommand <br/>\n" if (DEBUG);
        my $insert_objv = $dsn->prepare($variantinsertcommand);
        $insert_objv->execute(
            $varpos_disp,         $ref_ch_disp,   $hit_ch_disp,
            $query_id,       $note,          $queryvarpos_disp,
            $assigned_locus, $coding_effect, $conservation,
            $conservation_pos, $variant_res, $locus_id
        );
        my $foo_recv     = $insert_objv->fetchrow_hashref();
        my $variation_id = $$foo_recv{"id"};
        $insert_objv->finish();

        print "inserted query variation ID: $variation_id <br/>\n" if (DEBUG);
    }
}

sub get_insertion_format {
    my ( $ref_ch, $hit_ch, $delcountref, $inscountref, $varpos, $queryvarpos ) =
      @_;

    # $note stores variant type
    my $note = "";

    if ( $ref_ch eq 'N') {
        if($varpos == 3107){
            $note.="rCRS placeholder";
            unless(SHOW_3107){
                $note.="(ignored)";
            }
        }else{
            $note .= "(ignored)";
        }
    }
    if($hit_ch eq 'N'){
        $note .= "(ignored)";
    }
    
    if ( $hit_ch eq '-' || $hit_ch eq ':' ) {
        $note .= "deletion";
        $$delcountref++; #this will affect the deletion count
        $hit_ch = ':';
    }
    elsif ( $ref_ch eq '-' || $ref_ch eq ':' ) {
        $$inscountref++; #this will affect the insertion count
        $note   .= "insertion";
        $ref_ch = ':';

        #Marie's insertion style
        if (USE_MARIE_STYLE_INSERTIONS) {

            #e.g. 43 : -> A is equivalent to 42 T->TA
            $varpos--;
            $queryvarpos--;
            my $ref_seq = $rCRS_seq->seq( $varpos, $varpos );

            print $ref_seq."->".$ref_seq.$hit_ch."\n" if (DEBUG);
            $ref_ch = $ref_seq;
            $hit_ch = $ref_seq . $hit_ch;
        }
    }
    elsif ($ref_ch eq 'C'
        && $hit_ch eq 'T' )
    {
        $note .= "transition";
    }
    elsif ($ref_ch eq 'T'
        && $hit_ch eq 'C' )
    {
        $note .= "transition";
    }
    elsif ($ref_ch eq 'A'
        && $hit_ch eq 'G' )
    {
        $note .= "transition";
    }
    elsif ($ref_ch eq 'G'
        && $hit_ch eq 'A' )
    {
        $note .= "transition";
    }
    else {
        if($ref_ch =~ m/[ACGT]/ && $hit_ch =~ m/[ACGT]/){
            $note .= "transversion";
        }else{
            $note .= "ambiguous (ignored)";
        }
    }

    return ( $ref_ch, $hit_ch, $varpos, $queryvarpos, $note );
}

sub _get_coding_effect {
    my ( $locus_type, $varpos, $ref_ch, $hit_ch, $species_list ) = @_;
    my $dsn = getDSN();
        
    # 0    id  int4    0   0   0   0       0
    #    0  locus_id    int4    0   0   0   1       0
    #    0  strand  varchar 0   0   0   1       0
    #    0  pos int4    0   0   0   1       0
    #    0  locus_pos   int4    0   0   0   1       0
    #    0  codon_pos   int4    0   0   0   1       0
    #    0  codon_trio  int4    0   0   0   1       0
    #    0  residue varchar 0   0   0   1       0
    #    0  nt  varchar 0   0   0   1       0
    #    1  33  L   14673   0   1   1   A   T
    #    2  33  L   14672   1   1   2   U   A
    #    3  33  L   14671   2   1   3   G   C
    #    4  33  L   14670   3   2   1   A   T

    my $rna_sql =
"SELECT rna.locus_id as locus,rna.locus_pos,rna.codon_pos,rna.codon_trio,rna.strand,rna.nt,coding.codon,coding.aa
            FROM mitomasterb.rna
            LEFT JOIN mitomasterb.coding ON rna.locus_id=coding.locus_id and rna.codon_pos=coding.codon_pos
            WHERE rna.pos=?;";

    
    print $rna_sql. "<br/>\n" if (DEBUG);
    my $query_handle = $dsn->prepare($rna_sql);
    $query_handle->execute($varpos);
    my (
        $locus,        $locus_pos, $codon_pos, $codon_trio,
        $strand,       $nt,        $codon,     $coding_aa,
        $conservation, $conservation_pos, $conservation_residue, $coding_effect,
        $variant_residue
    );
    $query_handle->bind_columns(
        undef,        \$locus,  \$locus_pos, \$codon_pos,
        \$codon_trio, \$strand, \$nt,        \$codon,
        \$coding_aa
    );

    if ( $query_handle->fetch() ) {
        print "locus: $locus<br/>\n" if (DEBUG);
        my $conservation_sql;
        my $conservation_handle;
        
        my $conservation_count=0;
        my $species_count;
        if ($coding_aa) {
            $conservation_pos     = $codon_pos;
            $conservation_residue = $coding_aa;
        }
        else {
            $conservation_pos     = $locus_pos;
            $conservation_residue = $nt;
        }
        if ($species_list) {
            my $conservation_sql =
"SELECT species,residue FROM mitomasterb.conservation WHERE locus = ? AND pos = ? AND species IN ($species_list)";
            $conservation_handle = $dsn->prepare($conservation_sql);
            $conservation_handle->execute( $locus, $conservation_pos );

            #    $species_list );
            # i don't know how to do binding with IN, doesn't seem to like it
        }
        else {
            my $conservation_sql =
"SELECT species,residue FROM mitomasterb.conservation WHERE locus = ? AND pos = ? ";
            $conservation_handle = $dsn->prepare($conservation_sql);
            $conservation_handle->execute( $locus, $conservation_pos );
        }

        my ( $species, $residue );
        $conservation_handle->bind_columns( undef, \$species, \$residue );

        print "homo residue:" . $conservation_residue . "<br/>\n" if (DEBUG);
        while ( $conservation_handle->fetch() ) {
            print $species, "\t", $residue, "<br/>\n" if (DEBUG);
            if ( $residue eq $conservation_residue ) {
                $conservation_count++;
            }
            $species_count++;
        }
        print
"conservation_count $conservation_count  species_count $species_count<br/>\n"
          if (DEBUG);
        if ($species_count) {
            $conservation = sprintf( "%.2f",
                ( $conservation_count * 100.00 / $species_count ) );
        }

        #are we talking protein coding?
        if ($coding_aa) {

            #hit_ch is the user sequence
            #ref_ch is the ref sequence
            if (   $ref_ch eq ':'
                || $hit_ch eq ':'
                || $ref_ch eq '-'
                || $hit_ch eq '-' )
            {
                $coding_effect = _get_indel_effect(
                    $locus,     $ref_ch,    $hit_ch, $locus_pos,
                    $codon_pos, $coding_aa, $strand
                );
            }    #indel
            elsif ( !( $hit_ch =~ /[ACGTU]/ ) || $ref_ch eq 'N' ) {

                #this includes ambiguous nucleotides
                $coding_effect = "unknown";
            }
            else {
                my $rna_query_nt = $hit_ch;
                $rna_query_nt =~ s/T/U/;

                #comp if other strand, this is ND6 or loci 33
                if ( $strand eq 'L' ) { $rna_query_nt =~ tr/ACGU/UGCA/; }

                my @split_codon = split //, $codon;
                $split_codon[ $codon_trio - 1 ] = $rna_query_nt;
                my $new_codon = join '', @split_codon;

                print "new codon: $new_codon <br/>" if (DEBUG);
		my $variant_codon='';
                if($new_codon =~ /^[ACGTU]{3}$/){
		    $variant_codon = $rCRS_SR->codon_code($new_codon);
		}else{
		    $variant_codon = '?';
		}
                if ( $coding_aa eq $variant_codon ) {
                    $coding_effect = "syn:";
                }
                else { $coding_effect = "non-syn:" }
                $coding_effect .= "$coding_aa=>$variant_codon";
                $variant_residue=$variant_codon;
            }
        }
        else {
             $variant_residue = $hit_ch;
              #$variant_residue =~ s/T/U/;
                 #comp if other strand, several rnas are this way
                 if ( $strand eq 'L' ) { $variant_residue =~ tr/ACGT/TGCA/; }
            if ( $locus_type eq 'r' ) {
                $coding_effect = 'rRNA';
            }
            elsif ( $locus_type eq 't' ) {
                $coding_effect = 'tRNA';
            }
        }
    }    #no query
    else {
        $coding_effect = "no coding locus found";
    }

    $query_handle->finish();
    $dsn->disconnect();
    return ( $coding_effect, $conservation, $conservation_pos, $conservation_residue, $variant_residue );
}

sub _get_indel_effect {
    my ( $locus, $ref_ch, $hit_ch, $locus_pos, $codon_pos, $coding_aa, $strand )
      = @_;
    my $coding_effect;
    my $rna_seq = Bio::Mitomaster::RNASeq->new(
        species_ref => $rCRS_SR,
        locus_id    => $locus,
        variants    => {}
    );
    my $aa_seq = Bio::Mitomaster::AASeq->new(
        species_ref => $rCRS_SR,
        locus_id    => $locus,
        variants    => {}
    );
    my $tx = $rna_seq->seq();

    if (   ( ( $ref_ch eq ':' || $ref_ch eq '-' ) && length($hit_ch) % 3 != 0 )
        || ( ( $hit_ch eq ':' || $hit_ch eq '-' ) && length($ref_ch) % 3 != 0 )
      )
    {
        $coding_effect = "frmshft:";
    }

    $coding_effect .= $coding_aa . "=>";

        $hit_ch =~ s/T/U/;
    if ( $strand eq 'L' ) {
        $hit_ch = reverse $hit_ch;
        $hit_ch =~ tr/ACGU/UGCA/;
    }

    #insertion
    if ( $ref_ch eq ':' || $ref_ch eq '-' ) {

        #print "orig:$tx\n";
        substr( $tx, $locus_pos, 0 ) = $hit_ch;

        #print "now: $tx\n";
        #my $prot=$aa_seq->seq();
        #print "prot:$prot\n";
        #take the codon that got screwed
        #33	8	3	2	L	T	UAU	Y
        #so if i am in locus_pos 8 codon_pos 2 of 1,2,3
    }
    else {

        #deletion
        substr( $tx, $locus_pos, length($ref_ch) ) = '';
    }
    my $remainder = substr( $tx, $locus_pos - $codon_pos - 1 );

    #print "remainder $remainder\n";#  if (DEBUG);
    my $offset_bp = 0;
    my $new_aa;
    for (
        ;	 
        $offset_bp < length($remainder) - 2
	 && substr($remainder,$offset_bp,3) =~ /^[ACGTU]{3}$/
	 && $rCRS_SR->codon_code( substr( $remainder, $offset_bp, 3 ) ) ne
        'TERM' ;
        $offset_bp = $offset_bp + 3
      )
    {
        #print "$offset_bp translating:".substr( $remainder, $offset_bp, 3 )."\n";
        $new_aa .= $rCRS_SR->codon_code( substr( $remainder, $offset_bp, 3 ) );
        #print $new_aa."\n";
    }

        #if we have 12 amino acids and MAX_NONSENSE is 10, just show all 12
        #if we have 15 or more, show 10 and (+5AA)
        if ( $offset_bp >= ( MAX_NONSENSE * 3 + 15 ) ) {
            $new_aa =
              substr( $new_aa, 0, MAX_NONSENSE ) . "(+"
              . ( length($new_aa) - (MAX_NONSENSE) ) . "AA)";
        }
        
        #not confident of this code, TODO: re-examine
        if(length($remainder)-$offset_bp==3){
            if ( substr($remainder,$offset_bp,3) =~ /^[ACGTU]{3}$/ && $rCRS_SR->codon_code( substr( $remainder, $offset_bp, 3 ) ) eq 'TERM' ) {
                $new_aa .= '(TERM)';
            }else{
                        $new_aa.="(RT)";
                    }
        }else{
                $new_aa.="(RT)";
        }
        
        #$new_aa .= length($remainder)-$offset_bp;
        
    $coding_effect .= $new_aa;
    return $coding_effect;
}
sub _haplogroup {
    my ( $variants, $query_id, $seq_region_type, $loclist ) = @_;
    my $dsn = getDSN();

    print $variants."\n".$query_id."\n".$seq_region_type."\n".$loclist."\n" if (DEBUG);
    
    my %actions = ( HaploGrep => \&parseHaploGrepOutput,
                    hmito => \&parseHmitoOutput );
                    
    #my ( $hapgroup, $hapdetid ) =
    #  parseHmitoOutput( $variants, $query_id, $seq_region_type );
    
    my $hp = HAPLOGROUP_PROGRAM;
    my ( $hapgroup, $hapdetid ) = $actions{$hp}->($variants, $query_id, $seq_region_type, $loclist );
    
    print $hapgroup. "<br/>\n" if ( $hapgroup && DEBUG );

    if ( !($hapdetid) ) {
        $hapdetid = undef;
    }
    if ( !($hapgroup) ) {
        $hapgroup = undef;
    }

    #preparing update statement
    my $updatecommand =
"UPDATE mitomasterb.query SET poslist=?, haplogroup=?, haplodetid=? where id=?";
    print $updatecommand. "<br/>\n" if (DEBUG);
    my $update_objv = $dsn->prepare($updatecommand);
    $update_objv->execute( $loclist, $hapgroup, $hapdetid, $query_id );
    $update_objv->finish();

    $dsn->disconnect();
}
1;
