#!/usr/bin/env perl
use warnings;
use strict;

use Gff::Parser;
use Gff::Feature;

use Data::Dumper;
$Data::Dumper::Sortkeys = 1;

my $gp = Gff::Parser->new(
    @ARGV ? (file => $ARGV[0]) : (),
    # region => {
    #     id => "MaV-gen-2.0",
    #     from => 6000,
    #     to => 8000,
    # }
);

$gp->add_condition(sub{$_[0]->type eq 'gene'});

my $fh = $gp->fh;

while (my $g = $gp->next_groups) {
    print Dumper([keys $g->{groups}]);
}

