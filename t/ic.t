use Test::More tests => 137;
#use Test::More 'no_plan';

use Chemistry::Mol;
use Math::VectorReal;
use strict;
use warnings;

my @files;
my $TOL = 0.0001;

BEGIN { 
    @files = glob "*.ic";
    use_ok('Chemistry::InternalCoords');
}

# Test for conversion of internal to cartesian
for my $in_file (@files) {
    my $out_file = $in_file;
    $out_file =~ s/ic$/out/;
    open IN, "<", $in_file or die "couldn't open $out_file: $!\n";
    open OUT, "<", $out_file or die "couldn't open $out_file: $!\n";
    my @rows_in = map { [split] } <IN>;
    my @rows_out = map { [split] } <OUT>;

    my $mol = Chemistry::Mol->new;
    my $n = 0;
    for my $row_in (@rows_in) {
        my $row_out = shift @rows_out;
        my $atom = $mol->new_atom;
        my $ic = Chemistry::InternalCoords->new($atom, @$row_in);
        my @calc_coords = $ic->add_cartesians->array;
        $n++;
        for my $axis (qw(x y z)) {
            my $got = shift @calc_coords;
            my $expected = shift @$row_out;
            ok(abs($got-$expected) < $TOL, 
                "$in_file: $axis($n); expected $expected, got $got");
        }
    }
}

# check new update and setters

# set up simple molecule with the following geometry (atom 1 at origin)
#
#                      4-3
#                        |
#                      1-2

my $mol = Chemistry::Mol->new;
$mol->new_atom for 1 .. 4;
$mol->atoms(1)->internal_coords(0,0,0,0,0,0);
$mol->atoms(2)->internal_coords(1,1,0,0,0,0);
$mol->atoms(3)->internal_coords(2,1,1,90,0,0);
$mol->atoms(4)->internal_coords(3,1,2,90,1,0);
$_->internal_coords->add_cartesians for $mol->atoms;

my $atom = $mol->atoms(4);
ok($atom->distance(vector(0,1,0))*1 < $TOL, 'check coord before');
my ($ref, $val) = $atom->internal_coords->dihedral;
is($ref, $mol->atoms(1), 'get ref before');
is($val, 0, 'get val before');

# change dihedral to 180
$atom->internal_coords->dihedral(180);
($ref, $val) = $atom->internal_coords->dihedral;
is($val, 180, 'get val after');
$atom->internal_coords->add_cartesians;
ok($atom->distance(vector(2,1,0))*1 < $TOL, 'check coord after dihedral change')
    or print $atom->coords;

# now move atom 4 back in cartesian space and update
$atom->coords(0,1,0);
ok($atom->distance(vector(0,1,0))*1 < $TOL, 'check coord after atom move');
$atom->internal_coords->update;
($ref, $val) = $atom->internal_coords->dihedral;
ok(abs($val - 0) < $TOL, 'check dihedral after atom move and update');


