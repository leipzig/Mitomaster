package Mitoprocess::MitoResults;
use Mitoprocess::MitoprocessConfig;
use Mitoprocess::Parsing;
use Mitoprocess::Results;

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
  qw(get_genbank_report get_genbank get_haplo_cnt get_haplo_total get_genbank_mut_haplofrq get_genbank_mut_allfrq get_genbank_haplo get_haplo_id get_genbank_seq_cnt);

#get a resultset reference from a queryID
sub get_genbank_report {
    my $queryID = shift;
    my $sql;

    #resultmv was a caching scheme from Medhi, not sure what it did
        $sql = "SELECT 
  mitomap.genbank.genbank_id,
  mitomap.genbank.tpos,
  mitomap.genbank.qpos,
  mitomap.genbank.tnt,
  mitomap.genbank.qnt,
  mitomap.genbank.ntchange,
  mitomap.genbank.allele,
  mitomap.genbank.calc_locus,
  mitomap.genbank.cal_aachange,
  mitomap.genbank.conservation,
  mitomap.genbank.haplogroup,
  
  mitomap.genbank.disea as patientphenotype,

  mitomap.genbank.mmutid,
  mitomap.genbank.rtmutid,
  mitomap.genbank.poly_id,
  
  mitomap.genbank.is_polym as is_polymorphism,
  mitomap.genbank.ismmut,
  mitomap.genbank.isrtmut
FROM
  mitomap.genbank
 WHERE genbank_id =? ORDER BY tpos";
 
    print $sql. "\n" if (DEBUG);
    my $dsn = getDSN();
    my $sth = $dsn->prepare($sql);
    $sth->execute($queryID);
    return \$sth;
}

#summary of a run
sub get_genbank {
    my ($pos,$ref,$alt)=@_;
    ##  *** put a primary key on the first place
#  $sql="select query.id, query.label, query.pstart, query.pend from mitomasterb.query where query.\"runID\"=$runID";
    my $sql =
    "SELECT distinct polym_genbank_full.genbank_id as genbank_id, genbank_reference.reference as reference, genbank_haplogroup.haplogroup as haplogroup
	FROM mitomap.polym_genbank_full
	LEFT JOIN mitomap.genbank_haplogroup ON polym_genbank_full.genbank_id = genbank_haplogroup.genbank_id
	LEFT JOIN mitomap.genbank_reference ON genbank_reference.genbank_id = genbank_haplogroup.genbank_id
	
	WHERE polym_genbank_full.pos=? and polym_genbank_full.ref=? and polym_genbank_full.alt = ?";
	#I clean up the haplogroup names beforehand
	#LEFT JOIN mitomap.haplogroup ON haplogroup.haplogroup_full = genbank_haplogroup.haplogroup
    my $dsn = getDSN();
    my $sth = $dsn->prepare($sql);
    $sth->execute($pos,$ref,$alt);
    return \$sth;
}

sub get_haplo_cnt{
	my ($pos,$ref,$alt,$haploid) = @_;
my $cnt_sql = "
    SELECT genbank_haplogroup.haplogroup, count (DISTINCT polym_genbank_full.genbank_id) as cnt
	FROM mitomap.polym_genbank_full JOIN mitomap.genbank_haplogroup
	ON polym_genbank_full.genbank_id=genbank_haplogroup.genbank_id
	WHERE polym_genbank_full.pos=? and polym_genbank_full.ref=? and polym_genbank_full.alt = ?
	and genbank_haplogroup.haplogroup=?
	GROUP BY genbank_haplogroup.haplogroup";

    my $dsn = getDSN();
    my $sth = $dsn->prepare($cnt_sql);
    $sth->execute($pos,$ref,$alt,$haploid);
    return \$sth;
}

sub get_haplo_total{
    my $haplo = shift;
    my $cnt_sql_all = 
    "SELECT count (DISTINCT genbank_haplogroup.genbank_id) as cnt
	FROM mitomap.genbank_haplogroup
	WHERE genbank_haplogroup.haplogroup=?";
	
    my $dsn = getDSN();
    my $sth = $dsn->prepare($cnt_sql_all);
    $sth->execute($haplo);
    my $hap_cnt=$sth->fetchrow_hashref;
    return $hap_cnt->{'cnt'};
}

sub get_genbank_mut_haplofrq{
    my ($pos, $ref, $alt, $haplo) = @_;
    my $query_genbank_haplocnt = 
    "SELECT count (DISTINCT genbank.genbank_id) as cnt
	FROM mitomap.genbank
	WHERE genbank.tpos=? and genbank.tnt=? and genbank.qnt=? and genbank.haplogroup=?";
    my $dsn = getDSN();
    my $sth = $dsn->prepare($query_genbank_haplocnt);
    $sth->execute($pos, $ref, $alt, $haplo);
    return \$sth;
}

sub get_genbank_mut_allfrq{
    my ($pos,$ref,$alt)=@_;
    my $query_genbank = 
    "SELECT count (DISTINCT genbank.genbank_id) as cnt
	FROM mitomap.genbank
	WHERE genbank.tpos=? and genbank.tnt=? and genbank.qnt=?";
	
    my $dsn = getDSN();
    my $sth = $dsn->prepare($query_genbank);
    $sth->execute($pos,$ref,$alt);
    return \$sth;
}

sub get_genbank_seq_cnt{
    my $query="SELECT count(distinct genbank_id) as cnt FROM mitomap.genbank";
    my $dsn = getDSN();
    my $sth = $dsn->prepare($query);
    $sth->execute;
    my $genbank_cnt=$sth->fetchrow_hashref;
    return $genbank_cnt->{'cnt'};
}

sub get_genbank_haplo{
    my $haplo = shift;
    my $query_genbank_haplo = 	
	"SELECT distinct genbank_reference.genbank_id as genbank_id, genbank_reference.reference as reference
	FROM mitomap.genbank_haplogroup, mitomap.genbank_reference
	WHERE genbank_reference.genbank_id = genbank_haplogroup.genbank_id and genbank_haplogroup.haplogroup=?";
	
    my $dsn = getDSN();
    my $sth = $dsn->prepare($query_genbank_haplo);
    $sth->execute($haplo);
    return \$sth;
}

sub get_haplo_id{
    my $genbank_id = shift;
    my $query_haplo_id = 	
	"SELECT distinct haplogroup as haplogroup_id
	FROM mitomap.genbank_haplogroup
	WHERE genbank_haplogroup.haplogroup and genbank_haplogroup.haplogroup=?";
	
    my $dsn = getDSN();
    my $sth = $dsn->prepare($query_haplo_id);
    $sth->execute($genbank_id);
    return \$sth;
}
