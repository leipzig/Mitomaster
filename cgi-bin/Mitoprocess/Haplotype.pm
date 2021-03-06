package Mitoprocess::Haplotype;
require Exporter;

use Array::Diff;
use Array::Utils qw(:all);
use Data::ArrayList;
use Mitoprocess::MitoprocessConfig;
use Mitoprocess::LocalSettings;
use Carp;
use strict;
use warnings;
use LWP::UserAgent;
use HTTP::Request::Common;

our @ISA    = qw(Exporter);
our @EXPORT = qw(parseHmitoOutput parseHaploGrepOutput);

sub parseHaploGrepOutput {
    my ( $variants_delimited, $query_id, $seq_region_type, $loclist ) = @_;

#input
# SampleId   Range   Haplogroup  Polymorphisms (delimited with tabs)
#  TestSample1  "16024-576;"    ?   16051G  16182C  16183C  16189C  16519C  73G 152C    263G    315.1C  356.1C  8281d   8282d   8283d   8284d   8285d   8286d   8287d   8288d   8289d   11914A
#  TestSample2  "16024-16569 ;1-576;"   ?   16224C  16311C  16519C  73G 146C    152C    263G    309.1C  315.1C  573.1C  573.2C  573.3C
#  TestSample3  "16024-576;2092;3552;4071;4491;4833;4883;8414;8473;9090;9824;10397;10400;11959;11969;12372;12771;13563;14502;14569;15487"   G2a3    16183C  16189C  16193.1C    16223T  16278T  16362C  73G 260A    263G    309.1C  315.1C  489C    4833G   10400T  13563G  14569A  523d    524d

#output
# SampleID        Range   Haplogroup      Quality Polymorphisms
# TestSample1     16024-16569;1-576;      U2+152  63.8    73G (yes)       152C (yes)      263G (yes)      16051G (yes)    315.1C (hotspot)        356.1C (globalPrivateMut)       8281d (outOfRange)      8282d (outOfRange)      8283d (outOfRange)      8284d (outOfRange)      8285d (outOfRange)      8286d (outOfRange)      8287d (outOfRange)      8288d (outOfRange)      8289d (outOfRange)      11914A (outOfRange)     16182C (hotspot)        16183C (hotspot)        16189C (localPrivateMut)        16519C (hotspot)
# TestSample2     16024-16569;1-576;      K2a     100.0   73G (yes)       146C (yes)      152C (yes)      263G (yes)      16224C (yes)    16311C (yes)    309.1C (hotspot)        315.1C (hotspot)        573.1CCC (globalPrivateMut)     16519C (hotspot)
# TestSample3     16024-16569;1-576;2092;3552;4071;4491;4833;4883;8414;8473;9090;9824;10397;10400;11959;11969;12372;12771;13563;14502;14569;15487;        G2a1d   100.0   73G (yes)       260A (yes)      263G (yes)      489C (yes)      4833G (yes)     10400T (yes)    13563G (yes)    14569A (yes)    16189C (yes)    16223T (yes)    16278T (yes)    16362C (yes)    309.1C (hotspot)        315.1C (hotspot)
#         523d (hotspot)  524d (hotspot)  16183C (hotspot)        16193.1C (hotspot)

    #print "haplogrep varlist".$variants_delimited."\n";

    
    my @variants = split /;/, $variants_delimited;
    my $polymorphisms='';
    if( !$loclist ){
        return ( "n/a no homology", undef );
    }else{
        #if there isn't enough coverage don't even haplotype it
        my $total_coverage=0;
        my @locpairs=split ';',$loclist;
        foreach my $locpair(@locpairs){
            my ($begin,$end) = $locpair =~ /(\d+)-(\d+)/;
            $total_coverage += abs($end-$begin);
        }
        if ($total_coverage < MIN_HAPLOTYPING_COVERAGE){
            return ("n/a low homology", undef);
        }
    }
    if ( scalar(@variants) == 0 ) {

        #this can happen with reference genomes
        $polymorphisms .= "\t1G";
    }
    else {
        my $variants_printed=0;
        foreach my $variant (@variants) {
            $variant =~ s/:/-/g;
            print "a variant:" . $variant . "\n" if (DEBUG);

#;A7-;-10A;C114T;C150T;T195C;G583A;A750G;T850C;T961-;AT14523-;A14527C;-14531GAAA;T15018G;T15019C;T15021-
            if ( my ( $ref, $pos, $mut ) =
                $variant =~ /([ACTGN]+|-)(\d+)([ACTGN]+|-)/ )
            {

                print "ref $ref pos $pos mut $mut\n" if (DEBUG);
                unless ( $pos == 3107 ) {
                    if ( $ref eq '-' ) {
                        my @ins_nts = split //, $mut;

       #don't think this is true:begin at 1 since the insertion format is A->ACC
       #if -316C
       #-14531GAAA
                        for (
                            my $intPos = 0 ;
                            $intPos < scalar(@ins_nts) ;
                            $intPos++
                          )
                        {
                            $variants_printed++;
                            $polymorphisms .=  "\t" 
                              . $pos . "."
                              . ( $intPos + 1 )
                              . $ins_nts[$intPos];
                            print "\t" 
                              . $pos . "."
                              . ( $intPos + 1 )
                              . $ins_nts[$intPos] . "\n"
                              if (DEBUG);
                        }
                    }
                    elsif ( $mut eq '-' ) {

                        #1348	1351	GG	:	deletion
                        my @deletion_nts = split //, $ref;
                        for (
                            my $delPos = $pos, my $nt = 0 ;
                            $nt < scalar(@deletion_nts) ;
                            $nt++, $delPos++
                          )
                        {
                            $variants_printed++;
                            $polymorphisms .=  "\t" . $delPos . "d";
                        }
                    }
                    else {
                        unless ( $ref eq $mut ) {
                            $variants_printed++;
                            $polymorphisms .=  "\t" . $pos . $mut;
                        }
                    }
                }
            }
        }
        if($variants_printed==0){$polymorphisms .=  "\t1G";}
    }


    my $tmpOutFile = '';
    my $success = 0;
    #saving the file
    if(USE_REMOTE_HAPLOGREP){
        my $userAgent = LWP::UserAgent->new(timeout => 30); 
        my $request = POST HAPLOGREP_SERVER,
            Content => [ name => $query_id, range => $loclist, polys => $polymorphisms ];

        my $response = $userAgent->request($request);

        if ($response->is_success) {
             my $tmpOutFile  = TEMP_HAPLO_PATH . $query_id . "_" . time . "hsd.out";
             #saving the file
             open( HAP_OUTPUT_FILE, ">$tmpOutFile" ) or die "$!";
             binmode HAP_OUTPUT_FILE;
             print HAP_OUTPUT_FILE $response->decoded_content;
             close HAP_OUTPUT_FILE;
             $success = 1;
        } else {
             print $response->status_line if DEBUG;
             $success = 0;
        }
    }else{
        my $tmpInFile  = TEMP_HAPLO_PATH . $query_id . "_" . time . "hsd";
        $tmpOutFile = $tmpInFile . ".out";
            open( HAP_INPUT_FILE, ">$tmpInFile" ) or die "$!";
            print HAP_INPUT_FILE "SampleId\tRange\tHaplogroup\tPolymorphisms\n";
            print "$query_id\t\"$loclist\"\t$polymorphisms\n" if (DEBUG);
            print HAP_INPUT_FILE "$query_id\t\"$loclist\"\t$polymorphisms";
            print HAP_INPUT_FILE "\n";
            close(HAP_INPUT_FILE);
            my $hapHome    = HAPLOGREP_HOME;
            my $phyVersion = PHYLOTREE_VERSION;
            my $hap_command =
        "java -jar $hapHome --in $tmpInFile --out $tmpOutFile --phylotree $phyVersion  > /dev/null 2>&1";
        if ( system($hap_command) == 0 ) { $success = 1; }else{$success = 0;}
    }
    if ( $success) {

#SampleID	Range	Haplogroup	Quality	Polymorphisms
#30135	1-9941;14421-16569;13861-14281;	H8+(114)	42.0	114T (yes)	195C (yes)	750G (yes)	146C (no)	263G (no)	709A (no)	1438G (no)	4769G (no)	8860G (no)	15326G (no)	16288C (no)	16362C (no)	7d (globalPrivateMut)	150T (localPrivateMut)	583A (globalPrivateMut)	850C (localPrivateMut)	961d (globalPrivateMut)	14523d (globalPrivateMut)	14524d (globalPrivateMut)	14527C (globalPrivateMut)	14531.1AAA (globalPrivateMut)	15018G (globalPrivateMut)	15019C (globalPrivateMut)	15021d (globalPrivateMut)
        open( HAP_OUTPUT_FILE, "<$tmpOutFile" ) or die "$!";
        my $header = <HAP_OUTPUT_FILE>;
        $header =~ /^SampleID/ or carp('Malformed Haplogrep Output');
        my $report = <HAP_OUTPUT_FILE>;
        chomp $report;
        my @report_fields = split /\t/, $report;
        my $haplogroup;

        if ( $report_fields[2] =~ /^([A-Za-z]+[0-9]+[A-Za-z]?)/ ) {
            ($haplogroup) = $report_fields[2] =~ /^([A-Za-z]+[0-9]+[A-Za-z]?)/;
        }
        else {
            $haplogroup = $report_fields[2];
        }
        my $quality = $report_fields[3];
        my @polymorphisms =
          @report_fields[ 4 .. ( scalar(@report_fields) - 1 ) ];
        my $polymorphisms_str = join ',', @polymorphisms;
        my $hapinsert =
"insert into mitomasterb.haplogrepdetails(\"queryID\",haplogroup,quality,polymorphisms) values (?,?,?,?) returning id;";

        # 1 id  int4    0   0   0   0       0
        # 0 queryID int4    0   0   0   0       0
        # 0 haplogroup  varchar 0   0   0   1       0
        # 0 quality float4  0   0   0   1       0
        # 0 expected_found  varchar 0   0   0   1       0
        # 0 expected_notfound   varchar 0   0   0   1       0
        # 0 global_priv_mutations   varchar 0   0   0   1       0
        # 0 local_priv_mutations    varchar 0   0   0   1       0
        #print $hapinsert. "<br/>";
        my $dsn         = getDSN();
        my $insert_obj2 = $dsn->prepare($hapinsert);
        $insert_obj2->execute( $query_id, $haplogroup, $quality,
            $polymorphisms_str );
        my $foo_rec = $insert_obj2->fetchrow_hashref();
        $insert_obj2->finish();
        $dsn->disconnect();
        my $hapdetail_id = $$foo_rec{"id"};
        return ( $haplogroup, $hapdetail_id );
    }
    else {
        return ( undef, undef );
    }
}

# this function takes a list of variant positions and does haplogrouping by running h-mito
sub parseHmitoOutput {

    my $varlist = $_[0];
    my $queryid = $_[1];
    my $runType = $_[2];

    print $varlist. "\n";

    my $hapdetail_id = undef;
    my $haplotype    = undef;

    #print $varlist,"<br>";
    my @vararray = split ';', $varlist;

# this is a flag to check if there is any variant used in haplogrouping. If not, we will not store any haplogrouping info
    my $enoughVars = 0;

  # checking if any variants have been called! if not, no need for haplogrouping
    my $hmitout    = "";
    my $hmito_home = HMITO_HOME;
    if ( length($varlist) > 0 ) {
        my $find    = ";";
        my $replace = " ";
        $varlist =~ s/$find/$replace/g;

        $hmitout =
          `python $hmito_home -i \"$varlist\" -n \"FULL\"  -r \"$runType\"`;

    }
    else {

        $hmitout = `python $hmito_home  -i \"\" -n \"FULL\"  -r \"$runType\"`;
    }

#print "python $hmito_home -i \"$varlist\" -n \"FULL\"  -r \"$runType\"<br/>\n";
#print $hmitout. "<br>";

#		my $hmitout = `python c:\\workspace\\mitomaster\\trunk\\h-mito\\h-mito.py -i \"$varlist\" -n \"FULL\"  -r \"Control_region\"`;

    # PARSING h-mito RESULTS

    my $hapmindist      = 1000;
    my $hapminlevel     = 1000;
    my $hapmindetail_id = 0;
    my $currentscore    = 2000;
    my $maxmutationused = 0;
    my $isBestScore     = 0;
    my @lines           = split /[\n\r\l]+/, $hmitout;

    # Looping through reported haplogroups and store them in database
    foreach my $line (@lines) {

        my @fields = split /\t/, $line;

        my $haplogroup = $fields[0];

#print $haplogroup."<br/>";
#  The predicted haplogroup depth in the global human mtDNA phylogenetic tree (http://www.phylotree.org/tree/main.htm)
        my $level = $fields[1];

# The total number of extra variants (not used in haplogrouping) and variants missing but required for the predicted haplogroup in the sequence tested.
        my $distance = $fields[2];

        # Number of variants that are missing for the designated haplotype
        my $missing = $fields[3];

        # Variants that are not required for the designated haplotype
        my $extra   = $fields[4];
        my $find    = ";";
        my $replace = " - ";
        $isBestScore = 0;

        if ( $haplogroup eq "haplogroup" ) {
            next;
        }
##################################
        # taking care or the rCRS sequence:
        if ( $haplogroup eq "H2a2a" ) {
            if ( $level == 15 ) {
                if ( $distance == 0 ) {
                    if ( $missing eq "-" ) {
                        if ( $extra eq "-" ) {
                            my $goodmuts    = "";
                            my $goodmutSize = 0;
                            my $newscore =
                              $distance - $goodmutSize * 1.5 + ( $level / 5 );

                   # my $newscore = ($goodmutSize/($distance)) + ( $level / 5 );

                            #my $haplevel     = $fields[1];

                            if ( $newscore < $currentscore ) {
                                $hapminlevel = $level;
                                $hapmindist  = $distance;
                                $haplotype   = $haplogroup;

                                #$hapmindetail_id = $hapdetail_id;
                                $maxmutationused = $goodmutSize;
                                $currentscore    = $newscore;
                                $isBestScore     = 1;
                            }
                            elsif ( $newscore == $currentscore ) {
                                if ( $goodmutSize > $maxmutationused ) {
                                    $hapminlevel = $level;
                                    $hapmindist  = $distance;
                                    $haplotype   = $haplogroup;

                                    #$hapmindetail_id = $hapdetail_id;
                                    $maxmutationused = $goodmutSize;
                                    $currentscore    = $newscore;
                                    $isBestScore     = 1;
                                }
                                elsif ( ( length($haplogroup) ) <
                                    ( length($haplotype) ) )
                                {
                                    $hapminlevel = $level;
                                    $hapmindist  = $distance;
                                    $haplotype   = $haplogroup;

                                    #$hapmindetail_id = $hapdetail_id;
                                    $maxmutationused = $goodmutSize;
                                    $currentscore    = $newscore;
                                    $isBestScore     = 1;
                                }

                            }
                            $missing = "";
                            $extra   = "";
                            my $hapinsert =
"insert into mitomasterb.haplodetails(\"queryID\",haplogroup,haplolevel,mutationdistance,missingmuts,extramuts,goodmuts,score,mutationsused) values ($queryid,E'$haplogroup',$level,$distance,'$missing','$extra','$goodmuts','$newscore',$goodmutSize) returning id;";

                            #print $hapinsert. "<br/>";# if (DEBUG);
                            my $dsn         = getDSN();
                            my $insert_obj2 = $dsn->prepare($hapinsert);
                            $insert_obj2->execute();
                            my $foo_rec = $insert_obj2->fetchrow_hashref();
                            $hapdetail_id = $$foo_rec{"id"};
                            $insert_obj2->finish();
                            $dsn->disconnect();

            # print "rCRS haplodetail_id ".$hapdetail_id."<br>"." $isBestScore";
                            if ($isBestScore) {
                                $hapmindetail_id = $hapdetail_id;

                            }
                            $enoughVars = 1;
                            next;

                        }
                    }

                }

            }
        }    # END OF rCRS handling
####################################

        my @extrarray    = split ';', $extra;
        my @missingarray = split ';', $missing;

        my $extramutSize   = scalar(@extrarray);
        my $missingmutSize = scalar(@missingarray);

   # removing extra mutations from all mutations to be stored in goodmuts column

        my $dal = Data::ArrayList->new( my $initialCapacity = 800 );

        for my $variant (@vararray) {
            if ( length($variant) > 0 ) {
                $dal->add($variant);
            }
            else {

            }
        }

        for my $position (@extrarray) {
            if ( $dal->contains( sub { $_ eq $position } ) ) {

                $dal->remove( $dal->indexOf( sub { $_ eq $position } ) );
            }
        }
        my @minus       = $dal->toArray();
        my $goodmutSize = scalar(@minus);

        if ( $goodmutSize lt 1 ) {
            next;
        }
        else {
            $enoughVars = 1;
        }

        # goodmuts are variants that are used for haplogrouping
        my $goodmuts = join( " - ", @minus );

        $missing =~ s/$find/$replace/g;
        $extra   =~ s/$find/$replace/g;
        $find    = " ";
        $replace = " ";
        $goodmuts =~ s/$find/$replace/g;

        $find    = "'";
        $replace = "\\'";
        $haplogroup =~ s/$find/$replace/g;

        #my $currentscore = $hapmindist + $hapminlevel;

# Scoring formula. If changed, make sure you update it in the # if ( $haplogroup eq "H2a2a" )# statement which handles rCRS.
# Lower score is better. We want to maximize $goodmutSize and minimize $distance and $level. But here we give more weight to $goodmutSize and less weight to $level
        my $newscore = $distance - $goodmutSize * 1.5 + ( $level / 5 );

        #my $haplevel     = $fields[1];

        # Does current haplotype have the best score?
        if ( $newscore < $currentscore ) {
            $hapminlevel = $level;
            $hapmindist  = $distance;
            $haplotype   = $haplogroup;

            #$hapmindetail_id = $hapdetail_id;
            $maxmutationused = $goodmutSize;
            $currentscore    = $newscore;
            $isBestScore     = 1;
        }
        elsif ( $newscore == $currentscore ) {
            if ( $goodmutSize > $maxmutationused ) {
                $hapminlevel = $level;
                $hapmindist  = $distance;
                $haplotype   = $haplogroup;

                #$hapmindetail_id = $hapdetail_id;
                $maxmutationused = $goodmutSize;
                $currentscore    = $newscore;
                $isBestScore     = 1;
            }
            elsif ( ( length($haplogroup) ) < ( length($haplotype) ) ) {
                $hapminlevel = $level;
                $hapmindist  = $distance;
                $haplotype   = $haplogroup;

                #$hapmindetail_id = $hapdetail_id;
                $maxmutationused = $goodmutSize;
                $currentscore    = $newscore;
                $isBestScore     = 1;
            }

        }

        # storing haplogroup information in the database
        my $hapinsert =
"insert into mitomasterb.haplodetails(\"queryID\",haplogroup,haplolevel,mutationdistance,missingmuts,extramuts,goodmuts,score,mutationsused) values ($queryid,E'$haplogroup',$level,$distance,'$missing','$extra','$goodmuts','$newscore',$goodmutSize) returning id;";

        #print $hapinsert. "<br/>";
        my $dsn         = getDSN();
        my $insert_obj2 = $dsn->prepare($hapinsert);
        $insert_obj2->execute();
        my $foo_rec = $insert_obj2->fetchrow_hashref();
        $insert_obj2->finish();
        $dsn->disconnect();
        $hapdetail_id = $$foo_rec{"id"};
        if ($isBestScore) {
            $hapmindetail_id = $hapdetail_id;

        }
    }

    if ( !$enoughVars ) {
        $haplotype       = undef;
        $hapmindetail_id = undef;
    }

    # return best scored haplogroup and its id in the haplodetails database
    return ( $haplotype, $hapmindetail_id );
}
