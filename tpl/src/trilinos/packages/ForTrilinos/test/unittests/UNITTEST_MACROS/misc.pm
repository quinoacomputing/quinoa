package misc;

use strict;
use warnings;

######################################################################
# format_block: given a string containing one or more lines, indent
#               them all with nindent spaces
######################################################################
sub format_block {
  my ($block, $nindent) = @_;

  my $ind = ' ' x $nindent;
  $block =~ s/^/${ind}/mg;

  return $block;
}

######################################################################
# get_piece: given an array of terms, return as many as will fit
#            within a string of a specified length; if no terms will
#            fit, break off part of the first term and return it
######################################################################
sub get_piece {
  my ($terms, $len_middle, $extra_last) = @_;

  my $piece = "";
  my $all_done = 0;
  my @remains;

  my $tcnt = scalar @$terms;

  my $tind = 0;
  my $filled = 0;
  while (($tind < $tcnt) && !$filled) {
    my $t = $$terms[$tind];
    # the last line can be longer since no continuation needed
    my $len = ($tind == $tcnt-1 ? $len_middle + $extra_last : $len_middle);
    if (length($piece . $t) > $len) {
      # this term didn't fit, so exit the loop now
      $filled = 1;
    } else {
      # this fit, so get ready to try the next term
      $piece .= $t;
      $tind++;
    }
  }

  if ($piece eq "") {
    # if nothing fit, add part of the first term
    my $t = $$terms[0];
    $piece = substr($t, 0, $len_middle);
    @remains = @$terms;
    $remains[0] = substr($t, $len_middle);
  } elsif ($filled) {
    # let the caller know what we didn't include
    for (my $i = $tind; $i < $tcnt; $i++) {
      push (@remains, $$terms[$i]);
    }
  } else {
    # great, everything fit!
    $all_done = 1;
  }

  return ($piece, $all_done, \@remains);    
}

######################################################################
# wrap_line: wrap each line in the block, adding indentation and
#            continuation characters, so as not to exceed the maximum
#            allowable length
######################################################################
sub wrap_line {
  my ($block, $nindent, $max_len) = @_;

  if (length($block) <= $max_len) {
    return $block;
  }

  my $ind = ' ' x $nindent;
  my $pre_brk = "&";
  my $post_brk = $ind . "&";

  my $out_block = "";
  my @lines = split(/(?<=\n)/, $block);
  foreach my $newline (@lines) {
    chomp($newline);
    # split terms by spaces, commas, and left parens
    my @terms = split(/(?<=[\s(,])/, $newline);

    my $out_line = "";
    my $done = 0;
    while (!$done) {
      # depending if there's already indentation in outline, decide how much space we have
      my $this_mid_len = $max_len-length($out_line)-length($pre_brk);
      (my $piece, $done, my $new_terms) = get_piece(\@terms, $this_mid_len, length($pre_brk));
      @terms = @$new_terms;
      $out_line .= $piece;
      if (!$done) {
        # add this line to the block and get ready for next line
        $out_block .= $out_line . $pre_brk . "\n";
        $out_line = $post_brk;
      } else {
        # add this line to the block, and then we're finished
        $out_block .= $out_line . "\n";
        $out_line = "";
      }
    }
  }

  return $out_block;
}

# module load successful (required)
1;
