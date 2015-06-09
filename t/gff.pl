#!/usr/bin/env perl
use warnings;
use strict;

use Gff::Parser;
use Gff::Feature;

use Data::Dumper;
$Data::Dumper::Sortkeys = 1;

my $gp = Gff::Parser->new();
#$gp->is(sub{
#             my $feat = $_[0];
#             return $feat->type eq 'gene' && $feat->strand eq '-';
#         });

while (my $f = $gp->next_feature) {
    print Dumper $f; last;
}
