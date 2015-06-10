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
        region => {},
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

    # prepare region conditions
    if ($self->{region}) {
        if ($self->{region}{id}) {
            $self->add_condition( sub{ $_[0]->seqid eq $_[1]->{region}{id} });
        }
        if (defined $self->{region}{to} && defined $self->{region}{from}) {
            $self->add_condition(
                sub{
                    ($_[0]->start >= $_[1]->{region}{from} && $_[0]->start <= $_[1]->{region}{to}) || # ovl start
                        ($_[0]->end >= $_[1]->{region}{from} && $_[0]->end <= $_[1]->{region}{to}) || # ovl end
                            ($_[0]->start < $_[1]->{region}{from} && $_[0]->end > $_[1]->{region}{to}); # containing
                });
        }elsif (defined $self->{region}{from}) {
            $self->add_condition(sub{ $_[0]->end >= $_[1]->{region}{from} });
        }
    }

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

    while ( <$fh> ) {
        next if /^\s*$/;
        if (/^#/) {
            if (/^##FASTA/){
                $self->next_segment() || return; #eof
            }
            next;
        }

        # return gff feat object
        my $feat = Gff::Feature->new($_);
        $self->eval_feature($feat) || next;
        return $feat;
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
        next if /^\s*$/;
        return $_ if /^##/
    }
    #eof
    return;
}


=head2 next_groups

  groups =  {
    groups => {
      type => [{
        primary => feat(parent)
        children =>
      },...],
    }
    groups_by_id => {
      group_id => group, ...
    },
    features_by_id => {
      feat_id => feat, ...
    }
    children => [feat, ...]
  }


primary => {types}[{
     %feat_obj
     children => {types}[],
  },
  children => {type}[feat_objs],
  by_id => {id}{feat_obj}

=cut

sub next_groups{
    my ($self) = @_;
    my %g;
    my $fh = $self->{fh};

    while ( <$fh> ) {
        next if /^\s*$/;
        if (/^#/) {
            if (/^###/) {
                %g && last; # structure segment features
                next; # empty segment
            }
            if (/^##FASTA/){
                $self->next_segment() || last; #eof
            }
            next;
        }

        # return gff feat object
        my $feat = Gff::Feature->new($_);
        $self->eval_feature($feat) || next;


        if ($feat->parents) { # child
            # In case of multiple CDS per mRNA, maker does not create unique
            # ids. Need to make them unique...
            my $id = my $oid = $feat->id;
            my $x = 1;
            while (exists $g{features_by_id}{$id}){
                print "Non-uniq: $id\n";
                $x++;
                $id = "$oid:$x";
            }
            print "X:$x\nUniq: $id\n";
            $feat->id($id) if $x > 1;

            $g{features_by_id}{$id} = $feat;

            push @{$g{children}}, $feat;
        }else {
            die ($feat->id)." is not unique in GFF. Unique IDs are required for primary features." if exists $g{features_by_id}{$feat->id};
            $g{features_by_id}{$feat->id} = $feat;

            my $g = {primary => $feat, children => {}};
            $g{groups_by_id}{$feat->id} = $g;
            push @{$g{groups}{$feat->type}}, $g;
        }

    }

    if (%g) {
        # sort the children
        foreach my $feat ( @{$g{children}} ) {
            my $p = $feat;
            while (my ($id) = $p->parents) { # climb up children to primary parent
                die "$id" unless exists $g{features_by_id}{$id};
                $p = $g{features_by_id}{$id};
            }
            push @{$g{groups_by_id}{$p->id}{children}{$feat->type}}, $feat;
        }
    }

    return %g ? \%g : undef;
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


=head2 seek

Set the filehandle to the specified byte offset. Takes two
optional arguments POSITION (default 0), WHENCE (default 0), see perl "seek" for more.
Returns 'true' on success.

NOTE: this operation only works on real files, not on STDIN.

=cut

sub seek{
	my ($self, $offset, $whence) = (@_, 0, 0);
	return seek($self->fh, $offset, $whence);
}


=head2 add_condition/reset_conditions

Only return features from parser satisfying custom condition using a predefined
function. The function is called with the feature object as first
parameter. Only features that evaluate to TRUE are returned by the parser.

  # customize parser to only return 'gene' features from '-' strand.
  $gp->add_condition(sub{
             my $feat = $_[0];
             return $feat->type eq 'gene' && $feat->strand eq '-';
         });


  # deactivate conditions
  $gp->reset_conditions();

=cut

sub add_condition{
    my ($self, $cond) = @_;

    if ($cond && ref($cond) eq 'CODE') {
        $self->{cond} ||= [];
        push @{$self->{cond}}, $cond;
    } else {
        die (((caller 0)[3])." requires condition as CODE reference!\n");
    }
    return $self->{cond};
}

sub reset_conditions{
    my ($self, $cond) = @_;
    $self->{cond} = [];
}

=head2 eval_feature

Returns TRUE if feature matches "conditions" set for parser.

  $gp->eval_feature($feat)

=cut

sub eval_feature{
    my ($self, $feat) = @_;
    if ($self->{cond}) {
        foreach ( @{$self->{cond}} ){ $_->($feat, $self) || return; }
    }
    return 1;
}

##----------------------------------------------------------------------------##

=head1 AUTHOR

Thomas Hackl S<thomas.hackl@uni-wuerzburg.de>

=cut



1;
