use Test::More tests => 3;

BEGIN { use_ok( 'Mitoprocess::Genbank' ); }
  require_ok( 'Mitoprocess::Genbank' );

use Mitoprocess::Genbank;
use warnings;
use strict;

my $numresults=250;
my $gb = Mitoprocess::Genbank->new(length_range_href => FULL_LENGTH, numresults => $numresults);

cmp_ok(scalar(keys %{$gb->{'mitoseqs'}}),'==',$numresults);

foreach my $gi(keys %{$gb->{'mitoseqs'}}){
    print $gi."\t".$gb->{'mitoseqs'}{$gi}{'accession'}."\n";
}
