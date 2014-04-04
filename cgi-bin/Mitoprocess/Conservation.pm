package Mitoprocess::Conservation;
require Exporter;

use Mitoprocess::MitoprocessConfig;
use Mitoprocess::LocalSettings;
use Carp;
use strict;
use warnings;

our @ISA = qw(Exporter);
our @EXPORT =
  qw(get_conservation_window print_conservation_window_as_javascript get_species);

#return sort list of all species as array ptr
sub get_species {
    my $dsn = getDSN();    
    my @species_list;
    my @mammys = qw(Mammals Non-Mammals);
    foreach my $mammy (@mammys) {
        my $mammalsSQL =
"SELECT subheading,replace(subheading,' ','_') as subheading_link, MIN(species_order) as suborder from mitomasterb.species WHERE mainheading=? GROUP BY subheading ORDER BY suborder, subheading;";
        my $sth = $dsn->prepare($mammalsSQL);
        $sth->execute($mammy);

        while ( my $hash_ref = $sth->fetchrow_hashref ) {
            my $speciesSql =
"SELECT  species from mitomasterb.species WHERE subheading=? ORDER BY species_order";
            my $spec_handle = $dsn->prepare($speciesSql);
            $spec_handle->execute( $hash_ref->{subheading} );
            my ($species);
            $spec_handle->bind_columns( undef, \$species );

            while ( $spec_handle->fetch() ) {
                push @species_list, $species;
            }
        }    #species
    }    #mammys
    $dsn->disconnect();
    return \@species_list;
}    #get_species



#return residues for all species in a window of width CONS_FLANK around a position
#fill in missing positions, if any with dashes
#return results as hash ptr of positions, each pointing to a hash of species with residue values
sub get_conservation_window {
    my ( $locus, $pos ) = @_;
    my $cons_flank = CONS_FLANK;

    my $cons_query =
"SELECT  pos, conservation.species,residue FROM mitomasterb.conservation  WHERE locus = ? AND pos >= ? and pos <= ? ORDER BY pos";
    my $dsn                 = getDSN();
    my $conservation_handle = $dsn->prepare($cons_query);
    $conservation_handle->execute(
        $locus,
        $pos - $cons_flank,
        $pos + $cons_flank
    );

    my ( $qpos, $species, $residue );
    $conservation_handle->bind_columns( undef, \$qpos, \$species, \$residue );

    my %res;
    my %species_list;
    while ( $conservation_handle->fetch() ) {
        $res{$qpos}{$species} = $residue;
    }
    $dsn->disconnect();
    my %results;
    my $species_list_ptr = get_species();

#we iterate explicity because a window might be off the edge of a locus, just fill missing cols with dashes
#we can also transpose here
    for (
        my $qpos_iter = $pos - $cons_flank ;
        $qpos_iter <= $pos + $cons_flank ;
        $qpos_iter++
      )
    {

        # print $qpos_iter;
        if ( defined $res{$qpos_iter} ) {
            foreach my $species_member ( @{$species_list_ptr} ) {
                if ( defined $res{$qpos_iter}{$species_member} ) {
                    $results{$species_member}{$qpos_iter} =
                      $res{$qpos_iter}{$species_member};
                }
                else {
                    $results{$species_member}{$qpos_iter} = '-';
                }
            }
        }
        else {
            foreach my $species_member ( @{$species_list_ptr} ) {
                $results{$species_member}{$qpos_iter} = '-';
            }
        }

        # foreach my $species_member(@{$species_list_ptr}){
        #     print "\t".$results{$species_member}{$qpos_iter};
        # }
        # print "\n";
    }
    return \%results;
}

# { "aaData": [
#   [
#       "Trident",
#       "Internet Explorer 4.0",
#       "Win 95+",
#       {
#           "version": "4",
#           "grade": "X"
#       }
#   ],

sub print_conservation_window_as_javascript {
    my ( $rcrs,$locus, $pos, $res, $species_list_ptr ) = @_;

    my $cons_hr          = get_conservation_window( $locus, $pos );
    my $human_name = HUMAN_NAME;
    if(!defined $species_list_ptr || scalar(@{$species_list_ptr})==0){
        $species_list_ptr = get_species();
    }
    my $cons_flank       = CONS_FLANK;
    print "{ \"aaData\": [";
    my $jgrid = '';
    
    #here we print the mutation residue flanked by rCRS
    my $locus_hr=$rCRS_SR->locus($locus);
    $jgrid .= "[\"Mutation at rCRS pos $rcrs (locus pos $pos) of ".$locus_hr->{common_name};
    if($locus_hr->{strand} eq 'L'){$jgrid .="(light strand)"}
    $jgrid .="\"";
    for (
        my $qpos_iter = $pos - $cons_flank ;
        $qpos_iter <= $pos + $cons_flank ;
        $qpos_iter++
      )  {
            if($qpos_iter==$pos){$jgrid .= ",\"<strong>" . $res . "</strong>\""; }else{
            $jgrid .= ",\"<span class='muted'>" . $cons_hr->{$human_name}{$qpos_iter} . "</span>\"";
        }
        }
    $jgrid .= "],";
    
    foreach my $species_member ( @{$species_list_ptr} ) {
        my $pretty_species = $species_member;
        $pretty_species =~ s/\_/ /g;
        $jgrid .= "[\"$pretty_species\"";
        for (
            my $qpos_iter = $pos - $cons_flank ;
            $qpos_iter <= $pos + $cons_flank ;
            $qpos_iter++
          )
        {
            $jgrid .= ",\"" . $cons_hr->{$species_member}{$qpos_iter} . "\"";
        }
        $jgrid .= "],";
    }
    $jgrid =~ s/,$//;
    print $jgrid;
    print "]";

 #print ",\n";
 #print "\t\"aoColumns\": [\n\t{ \"sTitle\": \"Species\" }";
 #for(my $qpos_iter=$pos-$cons_flank;$qpos_iter<=$pos+$cons_flank;$qpos_iter++){
 #    print ",\n\t{\"sTitle\": \"$qpos_iter\"}";
 #}
 #print "\n\t]";
    print "}\n";
}
