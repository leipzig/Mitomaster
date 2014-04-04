use Test::More tests => 1;

BEGIN { use_ok( 'Mitoprocess::Conservation' ); }
  require_ok( 'Mitoprocess::Conservation' );
  
my $hash_ptr=get_conservation_window(32,4);

cmp_ok(scalar(keys %{$hash_ptr}),'==',45); #45 species
cmp_ok(scalar(keys %{$hash_ptr->{'Alligator_mississippiensis'}}),'==',11); #1 flanked by 5 each side

my $hash_ptr=get_conservation_window(32,4);
print_conservation_window_as_javascript(32,10);
