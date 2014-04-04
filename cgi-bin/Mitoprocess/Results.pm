package Mitoprocess::Results;
use Mitoprocess::MitoprocessConfig;
use Mitoprocess::Parsing;
use Bio::Mitomaster::SpeciesRef;
use Bio::Mitomaster::AASeq;
use Bio::Mitomaster::RNASeq;
use Mitoprocess::MitoResults; #for gb freq
require HTTP::Headers;
use strict;
use warnings;

use Bio::Tools::Run::StandAloneBlast;
use IO::String;
use File::Basename;
use List::Util;
use CGI qw(:standard);

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT =
  qw(get_refs print_report_header get_query_name get_query_haplogroup get_var_report print_var_report_txt get_query_ids get_summary rCRS_gff3 get_summary_gff3 get_variant_list print_summary_txt genbank_record_exists);

#Views:
#varrtmut
#SELECT variant.id, variant."refNA" AS refna, variant."queryNA", variant."position", variant.qposition, rtmutation.allele, rtmutation.rna, rtmutation.dz, rtmutation.locus, rtmutation.id AS rtmutid, count(*) AS refcount, COALESCE(true) AS rtmut
#FROM variant, mitomap.rtmutation, mitomap.rtmutation_reference
#WHERE ((((variant."position" = rtmutation."position") AND ((variant."queryNA")::text = (rtmutation.na)::text)) AND ((variant."refNA")::text = "substring"((rtmutation.allele)::text, 1, 1))) AND (rtmutation.id = rtmutation_reference.rtmutation_id))
#GROUP BY variant.id, variant."refNA", variant."queryNA", variant."position", variant.qposition, rtmutation.allele, rtmutation.rna, rtmutation.dz, rtmutation.locus, rtmutation.id
#ORDER BY variant."position"

#varpolymorphism
#SELECT variant.id, variant."refNA" AS refna, variant."queryNA", variant."position", variant.qposition, polymorphism.aachange, polymorphism.id AS polyid, count(*) AS refcount, COALESCE(true) AS polymorphism
#FROM variant, mitomap.polymorphism, mitomap.polymorphism_reference
#WHERE ((((variant."position" = polymorphism."position") AND ((variant."queryNA")::text = (polymorphism.regna)::text)) AND ((variant."refNA")::text = (polymorphism.refna)::text)) AND (polymorphism.id = polymorphism_reference.polymorphism_id))
#GROUP BY polymorphism.id, variant.id, variant."refNA", variant."queryNA", variant."position", variant.qposition, polymorphism.aachange
#ORDER BY variant."position"

#varmmut
#SELECT variant.id, variant."refNA" AS refna, variant."queryNA", variant."position", variant.qposition, mmutation.allele, mmutation.aa, mmutation.dz, mmutation.locus, mmutation.id AS mmutid, count(*) AS refcount, COALESCE(true) AS mmut
#FROM variant, mitomap.mmutation, mitomap.mmutation_reference
#WHERE ((((variant."position" = mmutation."position") AND ((variant."refNA")::text = (mmutation.refna)::text)) AND ((variant."queryNA")::text = (mmutation.regna)::text)) AND (mmutation.id = mmutation_reference.mmutation_id))
#GROUP BY variant.id, variant."refNA", variant."queryNA", variant."position", variant.qposition, mmutation.allele, mmutation.aa, mmutation.dz, mmutation.locus, mmutation.id
#ORDER BY variant."position"

sub _get_report_fields {
    my @fields = (
        "query",         "tpos",
        "qpos",          "tnt",
        "qnt",           "ntchange",
        "allele",        "calc_locus",
        "calc_aachange", "conservation",
        "haplogroup",    "patientphenotype",
      
        "mmutid",        "rtmutid",
        "polyid",        "is_polymorphism",
        "is_mmut",       "is_rtmut"
    );
    return @fields;
}

sub print_report_header {
    my $cgi = new CGI;
    print $cgi->header( -status => 200 );
    my @fields = _get_report_fields();
    print join( "\t", @fields ), "\n";
}

#return name of query sequencegiven queryid
sub get_query_name {
    my $queryID = shift;
    my $sql =
      "SELECT label AS query FROM mitomasterb.query WHERE id=?";
    print $sql. "<br/>\n" if (DEBUG);
    my $dsn = getDSN();
    my $sth = $dsn->prepare($sql);
    $sth->execute($queryID);
    my $hash_ref = $sth->fetchrow_hashref;
    return $hash_ref->{query};
}

#return name of query sequencegiven queryid
sub get_query_haplogroup {
    my $queryID = shift;
    my $sql =
      "SELECT haplogroup FROM mitomasterb.query WHERE id=?";
    print $sql. "<br/>\n" if (DEBUG);
    my $dsn = getDSN();
    my $sth = $dsn->prepare($sql);
    $sth->execute($queryID);
    my $hash_ref = $sth->fetchrow_hashref;
    return $hash_ref->{haplogroup};
}

#get a resultset reference from a queryID
sub get_var_report {
    my $queryID = shift;

    my $primaryCols = "SELECT
            variant.id,
            variant.\"position\"  AS tpos,
            variant.\"refNA\" AS tnt,
            variant.qposition AS qpos,
            variant.\"queryNA\"  AS qnt,
            variant.note  AS ntchange,
            variant.\"queryID\",
            variant.\"locus\"  as calc_locus,
            variant.\"aa_change\" as calc_aachange,
            round(conservation::numeric,2) || '%' as conservation,
            variant.conservation_pos,
            variant.variant_res,
            variant.locus_id,
            query.\"runID\",
            query.label  AS query,
            query.haplogroup";
    my $annoCols = ",COALESCE(mmutation.allele, '' ::character varying) ::text || COALESCE(
            rtmutation.allele, '' ::character varying) ::text AS allele,
            CASE
              WHEN rtmutation.dz IS NOT NULL AND mmutation.dz IS NOT NULL THEN ((
              COALESCE(mmutation.dz, '' ::character varying) ::text || ',' ::text) ||
              COALESCE(rtmutation.dz, '' ::character varying) ::text) ::character
              varying
              ELSE COALESCE(mmutation.dz, rtmutation.dz, '' ::character varying)
            END AS patientphenotype,
            
            COALESCE(rtmutation.id, 0) AS rtmutid,
            COALESCE(rtmutation.rna, '' ::character varying) AS rna,
            CASE WHEN COALESCE(rtmutation.id, 0) = 0 THEN 'false' ELSE 'true' END AS is_rtmut,
            
            COALESCE(mmutation.id, 0) AS mmutid,
            CASE WHEN COALESCE(mmutation.id, 0) = 0 THEN 'false' ELSE 'true' END AS is_mmut,
            
            COALESCE(polymorphism.id, 0) AS polyid,
            COALESCE(polymorphism.id) AS polymorphism,
            CASE WHEN COALESCE(polymorphism.id, 0) = 0 THEN 'false' ELSE 'true' END AS is_polymorphism
            
            ";
            #this might be preferable to the freq queries below:
            #polym_genbank_cnt.cnt as genbank_frequency
            #LEFT JOIN mitomap.polym_genbank_cnt ON polymorphism.id=polym_genbank_cnt.id
    my $freqCols = ",(SELECT count (DISTINCT genbank.genbank_id) FROM mitomap.genbank WHERE genbank.tpos=variant.\"position\" and genbank.tnt=variant.\"refNA\" and genbank.qnt=variant.\"queryNA\") as gb_frac
                    ,(SELECT count (DISTINCT genbank.genbank_id) FROM mitomap.genbank WHERE genbank.tpos=variant.\"position\" and genbank.tnt=variant.\"refNA\" and genbank.qnt=variant.\"queryNA\" and genbank.haplogroup=query.haplogroup) as hap_frac";
    
    my $primaryRels=" FROM mitomasterb.variant JOIN mitomasterb.query ON variant.\"queryID\" = query.id
              JOIN mitomasterb.run ON query.\"runID\" = run.id";
    my $annoRels = " 
          LEFT JOIN mitomap.mmutation ON variant.position = mmutation.position
                    AND variant.\"refNA\"::text = mmutation.refna::text
                    AND variant.\"queryNA\"::text = mmutation.regna::text
                    
          LEFT JOIN mitomap.polymorphism ON variant.position = polymorphism.position
                    AND variant.\"queryNA\"::text= polymorphism.regna::text
                    AND variant.\"refNA\"::text = polymorphism.refna::text
          
          LEFT JOIN mitomap.rtmutation ON variant.position = rtmutation.position
                    AND variant.\"queryNA\"::text = rtmutation.regna::text
                    AND variant.\"refNA\"::text = rtmutation.refna::text";
    
    my $primaryConditions=" WHERE \"queryID\"=? ORDER BY variant.id;";
    my $anno =1;
    my $sql = ($anno ? $primaryCols.$annoCols.$freqCols.$primaryRels.$annoRels.$primaryConditions : $primaryCols.$primaryRels.$primaryConditions);
    
   #print $sql. "\n"; exit(1);# if (DEBUG);
    my $dsn = getDSN();
    my $sth = $dsn->prepare($sql);
    $sth->execute($queryID);
    return \$sth;
}

sub get_refs {
    my ($rtmutid,$mmutid,$polyid) = @_;
    my $refagg;
    my $string_agg_sql_old="array_to_string(array_agg(reference.id::text), ',')";
    my $string_agg_sql_new="string_agg(reference.id::text,',')";
    if($polyid){
        my $polyRefs = "SELECT $string_agg_sql_old as refagg FROM mitomap.reference
              JOIN mitomap.polymorphism_reference ON (reference.id = polymorphism_reference.reference_id)
              WHERE polymorphism_reference.polymorphism_id=?";
              #SELECT array_to_string(array_agg(reference.id::text), ',') as refagg FROM mitomap.reference JOIN mitomap.polymorphism_reference ON (reference.id = polymorphism_reference.reference_id)
              #              WHERE polymorphism_reference.polymorphism_id= 1
              my $dsn = getDSN();
              my $sth = $dsn->prepare($polyRefs);
              $sth->execute($polyid);
              while ( my $hash_ref = $sth->fetchrow_hashref ) {
                  $refagg.=$hash_ref->{"refagg"};
              }
    }
    
    if($rtmutid){
        my $rtRefs = "SELECT $string_agg_sql_old as refagg FROM mitomap.reference
              JOIN mitomap.rtmutation_reference ON (reference.id = rtmutation_reference.reference_id)
              WHERE rtmutation_reference.rtmutation_id=?";
              my $dsn = getDSN();
              my $sth = $dsn->prepare($rtRefs);
              $sth->execute($rtmutid);
              while ( my $hash_ref = $sth->fetchrow_hashref ) {
                    $refagg.=$hash_ref->{"refagg"};
                }
    }
    
    if($mmutid){
        my $mmRefs = "SELECT $string_agg_sql_old as refagg FROM mitomap.reference
              JOIN mitomap.mmutation_reference ON (reference.id = mmutation_reference.reference_id)
              WHERE mmutation_reference.mmutation_id=?";
              my $dsn = getDSN();
              my $sth = $dsn->prepare($mmRefs);
              $sth->execute($mmutid);
              while ( my $hash_ref = $sth->fetchrow_hashref ) {
                    $refagg.=$hash_ref->{"refagg"};
                }
    }
    return $refagg;
}
#for CLI and webservice
sub print_var_report_txt {
    my $queryID = shift;
    my $sth_ptr = get_var_report($queryID);
    my @fields  = _get_report_fields();
    
    while ( my $hash_ref = $$sth_ptr->fetchrow_hashref ) {
        #my $refs = get_refs($hash_ref->{rtmutid},$hash_ref->{mmutid},$hash_ref->{polyid});
        #print "my refs:".$refs."\n";
        foreach ( @$hash_ref{@fields} ) { $_ = '' unless defined }
        print join( "\t", @$hash_ref{@fields} ), "\n";
    }
}

#get an array of query_ids (sequence results) from one run id
sub get_query_ids {
    my $runID = shift;
    my $sql =
"select query.id from \"mitomasterb\".query WHERE \"query\".\"runID\"=$runID"
      . " ORDER BY query.id";
    my $dsn = getDSN();
    my $sth = $dsn->prepare($sql);
    $sth->execute;
    my @query_ids;
    while ( my $hash_ref = $sth->fetchrow_hashref ) {
        push( @query_ids, $hash_ref->{id} );
    }
    return @query_ids;
}

#summary of a run
#uses views: countpolymorphvariants, countmmutvariants, countrtmuvariants, countallvariants
#TODO: evaluate usefulness of views
# SELECT variant."queryID", count(*) AS count FROM mitomap.mmutation, variant WHERE (((variant."position" = mmutation."position") AND ((variant."queryNA")::text = (mmutation.regna)::text)) AND ((variant."refNA")::text = (mmutation.refna)::text)) GROUP BY variant."queryID";
# AS SELECT variant."queryID", count(*) AS count FROM mitomap.polymorphism, variant WHERE (((variant."position" = polymorphism."position") AND ((variant."queryNA")::text = (polymorphism.regna)::text)) AND ((variant."refNA")::text = (polymorphism.refna)::text)) GROUP BY variant."queryID";
# SELECT variant."queryID", count(*) AS count FROM mitomap.rtmutation, variant WHERE (((variant."position" = rtmutation."position") AND ((variant."queryNA")::text = (rtmutation.na)::text)) AND ((variant."refNA")::text = "substring"((rtmutation.allele)::text, 1, 1))) GROUP BY variant."queryID";
sub get_summary {
    my $runID = shift;
    ##  *** put a primary key on the first place
#  $sql="select query.id, query.label, query.pstart, query.pend from mitomasterb.query where query.\"runID\"=$runID";
my $noviewsql =
"SELECT query.id, query.label, query.pstart, query.pend, query.poslist,query.haplogroup, count(mitomasterb.variant.id) as varcount, count(mitomap.polymorphism.id) as polycount, count(mitomap.rtmutation.id) as rtcount, count(mitomap.mmutation.id) as mmcount, count(mitomap.rtmutation.id)+count(mitomap.mmutation.id) as patvariants
FROM 
\"mitomasterb\".query
JOIN mitomasterb.variant on variant.\"queryID\" = query.id
LEFT JOIN mitomap.polymorphism on (((variant.\"position\" = polymorphism.\"position\") AND ((variant.\"queryNA\")::text = (polymorphism.regna)::text)) AND ((variant.\"refNA\")::text = (polymorphism.refna)::text))
LEFT JOIN mitomap.rtmutation on (((variant.\"position\" = rtmutation.\"position\") AND ((variant.\"queryNA\")::text = (rtmutation.regna)::text)) AND ((variant.\"refNA\")::text = (rtmutation.refna)::text))
LEFT JOIN mitomap.mmutation on (((variant.\"position\" = mmutation.\"position\") AND ((variant.\"queryNA\")::text = (mmutation.regna)::text)) AND ((variant.\"refNA\")::text = (mmutation.refna)::text))
WHERE \"query\".\"runID\"=? GROUP BY query.id, query.label, query.pstart, query.pend, query.poslist,query.haplogroup";


    my $viewsql =
"select query.id, query.label, query.pstart, query.pend, query.poslist,query.haplogroup, COALESCE(\"mitomasterb\".\"countallvariants\".\"count\", 0) as varcount, 
    COALESCE(\"mitomasterb\".\"countpolymorphvariants\".count, 0) as polycount,
    COALESCE(\"mitomasterb\".\"countmmutvariants\".\"count\", 0) + COALESCE(\"mitomasterb\".\"countrtmuvariants\".\"count\", 0) as patvariants
    FROM 
    \"mitomasterb\".query
    LEFT OUTER JOIN \"mitomasterb\".countallvariants
    ON \"mitomasterb\".query.id = \"mitomasterb\".\"countallvariants\".\"queryID\"
    LEFT OUTER JOIN \"mitomasterb\".countpolymorphvariants
    ON \"mitomasterb\".query.id = \"mitomasterb\".\"countpolymorphvariants\".\"queryID\"
    LEFT OUTER JOIN \"mitomasterb\".\"countmmutvariants\"
    ON \"mitomasterb\".query.id = \"mitomasterb\".\"countmmutvariants\".\"queryID\"
    LEFT OUTER JOIN \"mitomasterb\".countrtmuvariants
    ON \"mitomasterb\".query.id = \"mitomasterb\".\"countrtmuvariants\".\"queryID\" WHERE \"query\".\"runID\"=?";
    my $dsn = getDSN();
    my $sth = $dsn->prepare($noviewsql);
    $sth->execute($runID);
    return \$sth;
}

#gff3 of alignments
sub get_summary_gff3 {
    my $runID = shift;
    my $sth   = get_summary($runID);
    my @gff3;
    while ( my $hash_ref = $$sth->fetchrow_hashref ) {
        my @ranges = split( ";", $hash_ref->{poslist} );
        foreach my $range (@ranges) {
            my ( $start, $end ) = split( '-', $range );
            push @gff3,
                "rCRS\tUser\tYour_sequences\t" 
              . $start . "\t" 
              . $end
              . "\t.\t.\t.\tName="
              . $hash_ref->{label};
        }
    }
    return join( "\n", @gff3 ) . "\n";
}

sub rCRS_gff3 {
    return <<END
##gff-version 3
##sequence-region rCRS 1 16569
rCRS	Mitomap	rCRS_genes	14149	14673	.	-	.	Name=ND6;ID=33
rCRS	Mitomap	rCRS_genes	12337	14148	.	+	.	Name=ND5;ID=32
rCRS	Mitomap	rCRS_genes	8366	8572	.	+	.	Name=ATPase8;ID=21
rCRS	Mitomap	rCRS_genes	5904	7445	.	+	.	Name=COI;ID=16
rCRS	Mitomap	rCRS_genes	10470	10766	.	+	.	Name=ND4L;ID=27
rCRS	Mitomap	rCRS_genes	4470	5511	.	+	.	Name=ND2;ID=10
rCRS	Mitomap	rCRS_genes	14747	15887	.	+	.	Name=Cytb;ID=35
rCRS	Mitomap	rCRS_genes	7586	8269	.	+	.	Name=COII;ID=19
rCRS	Mitomap	rCRS_genes	10059	10404	.	+	.	Name=ND3;ID=25
rCRS	Mitomap	rCRS_genes	10760	12137	.	+	.	Name=ND4;ID=28
rCRS	Mitomap	rCRS_genes	8527	9207	.	+	.	Name=ATPase6;ID=22
rCRS	Mitomap	rCRS_genes	9207	9990	.	+	.	Name=COIII;ID=23
rCRS	Mitomap	rCRS_genes	3307	4262	.	+	.	Name=ND1;ID=6
rCRS	Mitomap	rCRS_rRNA	1671	3230	.	+	.	Name=16S;ID=4
rCRS	Mitomap	rCRS_rRNA	648	1601	.	+	.	Name=12S;ID=2
rCRS	Mitomap	rCRS_tRNA	4263	4331	.	+	.	Name=I;ID=7
rCRS	Mitomap	rCRS_tRNA	10405	10469	.	+	.	Name=R;ID=26
rCRS	Mitomap	rCRS_tRNA	7518	7585	.	+	.	Name=D;ID=18
rCRS	Mitomap	rCRS_tRNA	8295	8364	.	+	.	Name=K;ID=20
rCRS	Mitomap	rCRS_tRNA	12266	12336	.	+	.	Name=L(CUN);ID=31
rCRS	Mitomap	rCRS_tRNA	5512	5579	.	+	.	Name=W;ID=11
rCRS	Mitomap	rCRS_tRNA	12138	12206	.	+	.	Name=H;ID=29
rCRS	Mitomap	rCRS_tRNA	5587	5655	.	-	.	Name=A;ID=12
rCRS	Mitomap	rCRS_tRNA	5826	5891	.	-	.	Name=Y;ID=15
rCRS	Mitomap	rCRS_tRNA	7446	7516	.	-	.	Name=S(UCN);ID=17
rCRS	Mitomap	rCRS_tRNA	577	647	.	+	.	Name=F;ID=1
rCRS	Mitomap	rCRS_tRNA	12207	12265	.	+	.	Name=S(AGY);ID=30
rCRS	Mitomap	rCRS_tRNA	5761	5826	.	-	.	Name=C;ID=14
rCRS	Mitomap	rCRS_tRNA	9991	10058	.	+	.	Name=G;ID=24
rCRS	Mitomap	rCRS_tRNA	5657	5729	.	-	.	Name=N;ID=13
rCRS	Mitomap	rCRS_tRNA	1602	1670	.	+	.	Name=V;ID=3
rCRS	Mitomap	rCRS_tRNA	15888	15953	.	+	.	Name=T;ID=36
rCRS	Mitomap	rCRS_tRNA	4402	4469	.	+	.	Name=M;ID=9
rCRS	Mitomap	rCRS_tRNA	4329	4400	.	-	.	Name=Q;ID=8
rCRS	Mitomap	rCRS_tRNA	14674	14742	.	-	.	Name=E;ID=34
rCRS	Mitomap	rCRS_tRNA	15955	16023	.	-	.	Name=P;ID=37
rCRS	Mitomap	rCRS_tRNA	3230	3304	.	+	.	Name=L(UUA/G);ID=5
END
}

#this is not used anywhere but it might be handy
sub print_summary_txt {
    my $runID = shift;
    my $sth   = get_summary($runID);
    print
"seq\tLabel\tlocations\thaplogroup\tvariants\tpopvariants\tpatientvariants\tvariants\n";
    while ( my $hash_ref = $$sth->fetchrow_hashref ) {
        
        print $hash_ref->{id} . "\t"
          . $hash_ref->{label} . "\t"
          . $hash_ref->{pstart} . "\t"
          . $hash_ref->{pend} . "\t"
          . $hash_ref->{poslist} . "\t"
          . $hash_ref->{haplogroup} . "\t"
          . $hash_ref->{varcount} . "\t"
          # . $hash_ref->{polycount} . "\t"
          #           . $hash_ref->{patvariants} . "\t"
          . get_variant_list($hash_ref->{id})
    }
}

sub get_variant_list {
    my $queryid = shift;
    my $sql = "select variant.\"position\",variant.\"refNA\",variant.\"queryNA\" from \"mitomasterb\".variant where variant.\"queryID\" = ? ORDER BY position";
    my $dsn = getDSN();
    my $sth = $dsn->prepare($sql);
    $sth->execute($queryid);
    my @results;
    while ( my $hash_ref = $sth->fetchrow_hashref ) {
       my $pretty_query = $hash_ref->{'queryNA'};
       $pretty_query =~ s/:/d/;
       my $res=$hash_ref->{'refNA'}.$hash_ref->{'position'}.$pretty_query;
       push @results,$res;
     }
     return join ', ',@results;
}

#has this genbank sequence already been seen?
sub genbank_record_exists {
    my $genbank_id = shift;
    my $sql = "select * from mitomap.genbank where genbank_id = ?";
    my $dsn = getDSN();
    my $sth = $dsn->prepare($sql);
    $sth->execute($genbank_id);
    if($sth->fetchrow_hashref){
        return 1;
    }
    return 0;
}
