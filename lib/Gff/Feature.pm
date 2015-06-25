package Gff::Feature;

use warnings;
use strict;

use overload '""' => \&string;

# preference libs in same folder over @INC
use lib '../';

our $VERSION = '0.2.1';

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

Create a gff feature object. Takes either a gff line or a key => value
 representation of the gff fields.

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
    $self->attributes(ref $self->{attributes} ? %{$self->{attributes}} : $self->{attributes});

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

=cut

sub seqid{
    my ($self, $seqid, $force) = @_;
    if (defined $seqid || $force) {
        $self->{seqid} = $seqid;
    }
    return $self->{seqid};
}

=head2 source

=cut

sub source{
    my ($self, $source, $force) = @_;
    if (defined $source || $force) {
        $self->{source} = $source;
    }
    return $self->{source};
}

=head2 type

=cut

sub type{
    my ($self, $type, $force) = @_;
    if (defined $type || $force) {
        $self->{type} = $type;
    }
    return $self->{type};
}

=head2 start

=cut

sub start{
    my ($self, $start, $force) = @_;
    if (defined $start || $force) {
        $self->{start} = $start;
    }
    return $self->{start};
}

=head2 end

=cut

sub end{
    my ($self, $end, $force) = @_;
    if (defined $end || $force) {
        $self->{end} = $end;
    }
    return $self->{end};
}

=head2 score

=cut

sub score{
    my ($self, $score, $force) = @_;
    if (defined $score || $force) {
        $self->{score} = $score;
    }
    return $self->{score};
}

=head2 strand

=cut

sub strand{
    my ($self, $strand, $force) = @_;
    if (defined $strand || $force) {
        $self->{strand} = $strand;
    }
    return $self->{strand};
}

=head2 phase

=cut

sub phase{
    my ($self, $phase, $force) = @_;
    if (defined $phase || $force) {
        $self->{phase} = $phase;
    }
    return $self->{phase};
}

=head2 attributes

Get/set attributes. Takes STRING or HASH structure. Returns HASHREF structure.

  $feat->attributes("ID=gene1;NAME=Gene1");
  $feat->attributes(ID=>"gene1", NAME=>"Gene1");

=cut

sub attributes{
    my ($self, @attr) = @_;

    if ( @attr%2 == 0) {
        $self->{attributes} = {@attr};
        die (((caller 0)[3]).": attribute values need to be ARRAYREFs\n")
            if grep{ref $_ ne 'ARRAY'} values %{$self->{attributes}}; # require ARRAYREFS

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


=head2 attr/attribute

Get/set attributes.

NOTE: Always returns LIST of values as some attributes can have multiple values,
e.g. "Parent".

  my ($id) = $feat->attr('ID');
  my @parent_ids = $feat->attr('Parent');

  $feat->attr('ID', ["new-ID"]);

=cut

# alias for lazyness
*attr = \&attribute;

sub attribute{
    my ($self, $tag, $value) = @_;
    die (((caller 0)[3]).": $tag required") unless $tag;
    if (@_ > 2) {
        die (((caller 0)[3])."($tag): Value needs to be ARRAYREF\n") if defined $value && ref($value) ne "ARRAY";
        if (defined $value && @$value) {
            $self->{attributes}{$tag}=$value;
        }else {
            delete $self->{attributes}{$tag};
        }
    }
    my $attr = $self->{attributes}{$tag};
    return defined $attr ? @$attr : ();
}


=head2 id, parents

Convenience functions to attributes.

=cut

sub id{
    my $self = shift;
    my ($id) = $self->attr('ID', @_>0 ? [@_] : ());
    return $id;
}

sub parents{
    my $self = shift;
    return $self->attr('Parent', @_);
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

    my @as = ();
    foreach (@ATTRIBUTES) {
        if (my @attr = $self->attr($_)){
            push @as, "$_=". join(",", @attr);
        }
    }

    my @custom_attr = grep{! $ATTRIBUTES{$_}}keys %{$self->{attributes}};
    if (@custom_attr) {
        foreach (sort @custom_attr) {
            if (my @attr = $self->attr($_)){
                push @as, "$_=". join(",", @attr);
            }
        }
    }

    return join(";", @as);
}


##----------------------------------------------------------------------------##

=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;
