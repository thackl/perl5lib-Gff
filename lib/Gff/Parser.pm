package Gff::Parser;

use warnings;
use strict;

# preference libs in same folder over @INC
use lib '../';

#use Gff::Feature;

our $VERSION = '0.1.0';

=head1 NAME

Gff::Parser.pm

=head1 DESCRIPTION

Parser module for SAM format files.

=cut

=head1 SYNOPSIS

  use Sam::Parser;
  use Sam::Alignment ':flags';

  my $sp = Sam::Parser->new(
      # parse from a file
    file => '/some/sam/file.sam',
      #or from a file handle
    fh => \*SAM,
      # or even from a STRINGREF
    file => \$sam_string,
  );

  # print names of all reference sequences from header
  while(%h = $sp->next_header_line('@SQ')){
  	print $h{'SN'}."\n";
  }

  # parser for file handle and with customized is routine
  # read starts with 'C'
  my $sp = Sam::Parser->new(
    fh => \*SAM,
    is => sub{ substr($_[0]->seq, 0, 1) eq 'C' }
  );

  # print read ids of all reads with bad quality
  while( my $aln = $sp->next_aln() ){
    print $aln->qname() if $aln->is_bad_quality();
  }

  # seek the begin of the alignments for reparsing
  $sp->seek_alignment_section();

  # reset the 'is' routine
  $sp->is(MAPPED_BOTH);

  # print sequences of read pairs with both reads mapped
  while( my ($aln1, $aln2) = $sp->next_pair() ){
    print $aln1->seq().", ".$aln2->seq()."\n";
  }

=cut

=head1 Constructor METHOD

=head2 new

Initialize a gff parser object. Takes parameters in key => value format.

  fh => \*STDIN,
  file => undef,
  is => undef,
  mode => '<',   # read,
                 # '+>': read+write (clobber file first)
                 # '+<': read+write (append)
                 # '>' : write (clobber file first)
                 # '>>': write (append)
=back

=cut

sub new{
    my $class = shift;

    my $self = {
        # defaults
        fh => \*STDIN,
        file => undef,
        is => undef,
        mode => '<',
        # overwrite defaults
        @_,
    };

    # open file in read/write mode
    if ($self->{file}) {
        my $fh;
        open ( $fh , $self->{mode}, $self->{file}) or die sprintf("%s: %s, %s",(caller 0)[3],$self->{file}, $!);
        $self->{fh} = $fh;
    }

    bless $self, $class;

    # prepare "is" test routine
    $self->is($self->{is});

    return $self;

}

sub DESTROY{
    # just to be sure :D
    my $self = shift;
    close $self->fh if $self->fh;
}








############################################################################


=head1 Object METHODS

=cut

=head2 next_feature

Loop through gff file and return next 'Gff::Feature' object (meeting the
 'is' criteria if specified).

=cut

sub next_feature{
    my ($self) = @_;
    my $fh = $self->{fh};

    while (<$fh>) {
        if (/^#/) {
            if (/^##FASTA/){
                $self->next_segment() || last; #eof
            }
            next;
        }

        # return gff feat object
        my $feat = Gff::Feature->new($_);
        return $feat if !defined($self->{is}) || &{$self->{is}}($feat);
    }
    return;
}


=head2 next_segment

Loop through gff file until next segment (##-directive).

=cut

sub next_segment{
    my ($self) = @_;
    my $fh = $self->{fh};

    while (<$fh>) {
        return $_ if /^##/
    }
    #eof
    return;
}


=head2 append_feature

Append an alignment to the file, provided as object or string. Returns the
 byte offset position in the file.

NOTE: In case a string is provided, make sure it contains trailing newline
 since no further test is performed.

=cut

sub append_feature{
    my ($self, $feat) = @_;
    my $pos = $self->tell;
    print {$self->{fh}} "$feat";
    return $pos;
}


=head2 tell

Return the byte offset of the current append filehandle position

=cut

sub tell{
    return tell($_[0]->{fh});
}


############################################################################

=head1 Accessor METHODS

=head2 fh

Get/Set the file handle.

=cut

sub fh{
    my ($self, $fh) = @_;
    $self->{fh} = $fh if $fh;
    return $self->{fh};
}


=head2 is

Only return features from parser satisfying custom criteria using a predefined
function. The function is called with the feature object as first
parameter. Only features that evaluate to TRUE are returned by the parser.

  # customize parser to only return 'gene' features from '-' strand.
  $gp->is(sub{
             my $feat = $_[0];
             return $feat->type eq 'gene' && $feat->strand eq '-';
         });


  # deactivate testing
  $gp->is(0);

=cut

sub is{
    my ($self, $is) = @_;

    if (@_== 2) {
        unless($is){
            $self->{is} = undef;
        } elsif (ref($is) eq 'CODE') {
            $self->{is} = $is;
        } else {
            die (((caller 0)[3])." requires CODE reference!\n");
        }
    }
    return $self->{is};
}


##----------------------------------------------------------------------------##

=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;
