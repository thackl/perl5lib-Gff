package Gff::Feature;

use warnings;
use strict;

use overload '""' => \&string;

# preference libs in same folder over @INC
use lib '../';

our $VERSION = '0.1.0';

=head1 NAME

Gff::Feature.pm

=head1 DESCRIPTION

Class for handling gff features.

=cut

=head1 SYNOPSIS

  use Gff::Feature;

=cut


##----------------------------------------------------------------------------##

=head1 Class ATTRIBUTES

=cut

# attributes is composite field and treated specially
our @FIELDS = qw(seqid source type start end score strand phase);

our @ATTRIBUTES = qw(ID Name Alias Parent Target Gap Derives_from Note Dbxref Ontology_term Is_circular);

my %ATTRIBUTES;
@ATTRIBUTES{@ATTRIBUTES} = (1) x @ATTRIBUTES;

##----------------------------------------------------------------------------##

=head1 Class METHODS

=head1 Constructor METHOD

=head2 new

Create a sam alignment object. Takes either a sam entry as as string (one
 line of a sam file) or a key => value representation of the sam fields
 C<qname flag rname pos mapq cigar rnext pnext tlen seq qual opt>.
 While the first eleven fields are regular, C<opt> contains a string of all
 the optional fields added to the line.

Returns a sam alignment object. For more informations on the sam format see
 L<http://samtools.sourceforge.net/SAM1.pdf>.


=cut

sub new{
    my $class = shift;
    my $self;

    if (@_ == 1) {              # input is string to split
        my $gff = $_[0];
        chomp($gff);
        my %gff;
        @gff{@FIELDS, 'attributes'} = split("\t",$gff, 9);
        $self = \%gff;
    } else {                    # input is key -> hash structure
        $self = {
            seqid => undef,
            source => undef,
            type => undef,
            start => undef,
            end => undef,
            score => undef,
            strand => undef,
            phase => undef,
            attributes => {},
            @_,
        };
    }

    bless $self, $class;
    
    # process attributes
    $self->attributes($self->{attributes});
    
    return $self;

}

=head1 Object METHODS

=cut

=head1 Accessor METHODS

Get/Set the field values.

  my $seqid = $feat->seqid(); # get
  $feat->seqid("Some_ID"); # set
  $feat->seqid(undef, 1); # reset

=head2 seqid

Get/set seqid.

=cut

sub seqid{
    my ($self, $seqid, $force) = @_;
    if (defined $seqid || $force) {
        $self->{seqid} = $seqid;
    }
    return $self->{seqid};
}

=head2 source

Get/set source.

=cut

sub source{
    my ($self, $source, $force) = @_;
    if (defined $source || $force) {
        $self->{source} = $source;
    }
    return $self->{source};
}

=head2 type

Get/set type.

=cut

sub type{
    my ($self, $type, $force) = @_;
    if (defined $type || $force) {
        $self->{type} = $type;
    }
    return $self->{type};
}

=head2 start

Get/set start.

=cut

sub start{
    my ($self, $start, $force) = @_;
    if (defined $start || $force) {
        $self->{start} = $start;
    }
    return $self->{start};
}

=head2 end

Get/set end.

=cut

sub end{
    my ($self, $end, $force) = @_;
    if (defined $end || $force) {
        $self->{end} = $end;
    }
    return $self->{end};
}

=head2 score

Get/set score.

=cut

sub score{
    my ($self, $score, $force) = @_;
    if (defined $score || $force) {
        $self->{score} = $score;
    }
    return $self->{score};
}

=head2 strand

Get/set strand.

=cut

sub strand{
    my ($self, $strand, $force) = @_;
    if (defined $strand || $force) {
        $self->{strand} = $strand;
    }
    return $self->{strand};
}

=head2 phase

Get/set phase.

=cut

sub phase{
    my ($self, $phase, $force) = @_;
    if (defined $phase || $force) {
        $self->{phase} = $phase;
    }
    return $self->{phase};
}

=head2 attributes

Get/set attributes.

=cut

sub attributes{
    my ($self, @attr) = @_;
    if (! @attr%2 ) {
        $self->{attributes} = {@attr};
    }elsif (@attr == 1) {
        $self->_hash_attributes(@attr)
    }else {
        die (((caller 0)[3])." requires either single string or HASH!\n");
    }
    
    return $self->{attributes};
}



=head2 string

Get stringified alignment. Overload for "".

=cut

# alias for backward comp.
*raw = \&string;

sub string{
    my ($self) = @_;
    #print Dumper($self); use Data::Dumper; exit;
    my $s = join("\t", @$self{@FIELDS});
    $s.= "\t".$self->_string_attributes;
    return $s."\n";
}


=head2 length

Convenience function, equals C<$feat->end - $feat->start>.

=cut

sub length{
    my ($self) = @_;
    return $self->end - $self->start;
}


=head2 att/attribute

=cut

# alias for lazyness
*attr = \&attribute;

sub attribute{
    my ($self, $tag) = @_;
    my $attr = $self->{attributes}{$tag};
    return defined $attr ? @$attr : ();
}


=head2 _hash_attributes

=cut

sub _hash_attributes{
    my ($self, $attr) = @_;

    unless ($attr) {
        $self->{attributes} = {};
        return {};
    }

    my %attr;
    while ( $attr =~ /([^,=;]+)=([^;]+)(?:;|$)/g ){
        $attr{$1} = [split(",", $2)];
    }

    return $self->{attributes} = \%attr;
}


=head2 _string_attributes

=cut

sub _string_attributes{
    my ($self, $attr) = @_;

    my $as = "";
    foreach (@ATTRIBUTES) {
        if (my @attr = $self->attr($_)){
            $as .= join(",", @attr) . ";";
        }
    }

    my @custom_attr = grep{! $ATTRIBUTES{$_}}keys %{$self->{attributes}};
    if (@custom_attr) {
        foreach (sort @custom_attr) {
            if (my @attr = $self->attr($_)){
                $as .= join(",", @attr) . ";";
            }
        }
    }
    
    return $as;
}


##----------------------------------------------------------------------------##

=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;
