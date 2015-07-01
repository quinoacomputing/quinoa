package parser;

use strict;
use warnings;

use Text::Balanced qw(extract_bracketed);

use macroexp;

######################################################################
# distinct_parens: separate each outer paren pair, including the text
#                  between pairs with the pair that follows it; for
#                  text with no following pair, it stands alone
######################################################################
sub distinct_parens {
  my ($line) = @_;

  my $tmp = $line;
  my $lpar = ($tmp =~ tr/\(/X/);
  my $rpar = ($tmp =~ tr/\)/X/);

  die "parentheses are unbalanced" if ($lpar != $rpar);

  my @groups;
  my $search_pos = 0;
  my $begin_group = 0;
  my $rank = 0;
  my $lp = index($line, '(', $search_pos);
  my $rp = index($line, ')', $search_pos);
  while (($lp >= 0) or ($rp >= 0)) {
    if (($lp >= 0) and (($lp < $rp) or ($rp < 0))) {
      $search_pos = $lp+1;
      $rank++;
    } elsif (($rp >= 0) and (($rp < $lp) or ($lp < 0))) {
      $search_pos = $rp+1;
      $rank--;
      if ($rank == 0) {
        my $ss = substr($line, $begin_group, $search_pos-$begin_group);
        push (@groups, $ss);
        $begin_group = $rp+1;
      }
    } else {
      die "unknown error when trying to isolate parentheses ";
    }
    die "incorrect parenthesis order" if ($rank < 0);
    $lp = index($line, '(', $search_pos);
    $rp = index($line, ')', $search_pos);
  }
  die "parentheses are unbalanced" if ($rank != 0);

  if ($search_pos != length($line)) {
    my $last_grp = "";
    if (scalar @groups > 0) {
      $last_grp = pop(@groups);
    }
    $last_grp .= substr($line, $begin_group);
    push(@groups, $last_grp);
  }

  return @groups;
}

######################################################################
# parse_parens: for string containing at most one pair of parens,
#               separate what comes before/within/after parens; if no
#               parens, return everything as "before" variable
######################################################################
sub parse_parens {
  my ($line) = @_;

  my $tmp = $line;
  my $lpar = ($tmp =~ tr/\(/X/);
  my $rpar = ($tmp =~ tr/\)/X/);

  die "parentheses are unbalanced" if ($lpar != $rpar);

  my ($before, $within, $after);

  if ($lpar == 0) {
    $before = $line;
    $within = undef;
    $after = "";
  } elsif ($line =~ /^([^\(]*)\((.*)\)([^\)]*)$/) {
    $before = $1;
    $within = $2;
    $after = $3;

    my $lp = index($within, "(");
    if ($lp >= 0) {
      my $rp = index($within, ")");
      die "incorrect parenthesis order" if ($rp < $lp);
    }
  } else {
    die "unable to parse parentheses";
  }

  return ($before, $within, $after);
}

######################################################################
# separate_key: when given all text before a pair of parens, separate
#               the macro name from anything else that may be there
######################################################################
sub separate_key {
  my ($pre_paren) = @_;

  my ($pre, $key);

  if ($pre_paren =~ /^(.*\W)(\w+)$/) {
    $pre = $1;
    $key = $2;
  } else {
    $pre = "";
    $key = $pre_paren;
  }

  return ($pre, $key);
}

######################################################################
# recursive_replace: recursively process nested and/or adjacent paren
#                    pairs, expanding any macros found
######################################################################
sub recursive_replace {
  my ($line) = @_;

  my $output = "";
  my @groups = distinct_parens($line);
  foreach my $grp (@groups) {
    my ($before, $within, $after) = parse_parens($grp);
    if (defined($within)) {
      my ($pre, $key) = separate_key($before);
      my $replaced = recursive_replace($within);
      if ($key ne "") {
        my $expanded = macroexp::expand_macro($key, $replaced);
        if (defined($expanded)) {
          $output .= $pre . $expanded . $after;
        } else {
          $output .= $grp;
        }
      } else {
        $output .= $grp;
      }
    } else {
      $output .= $grp;
    }
  }
  return $output;
}

######################################################################
# break_args: split up each term of an argument list, carefully
#             watching for commas that are protected by parentheses
######################################################################
sub break_args {
  my ($argstr) = @_;

  my $string = $argstr;
  my @params;
  while ($string ne "") {
    if ($string =~ /^([^(]*?),/) {
      push(@params, $1);
      $string =~ s/^\Q$1\E\s*,?\s*//;
    } else {
      my $tmp = $string;
      my $lpar = ($tmp =~ tr/\(/X/);
      if ($lpar > 0) {
        my ($ext, $new_string, $pre) = extract_bracketed($string,'()','[^()]*');
        $string = $new_string;
        push(@params, "$pre$ext");
        $string =~ s/^\s*,\s*//;
      } else {
        push(@params, "$string");
        $string = "";
      }
    }
  }
  return @params;
}

######################################################################
# expand_line: initiate the processing of non-blank lines
######################################################################
sub expand_line {
  my ($line) = @_;
  $line =~ s/^(\s+)//;
  $line =~ s/!(.*)$//;

  if ($line !~ /^\s*$/) {
    return recursive_replace($line);
  } else {
    return undef;
  }
}

# module load successful (required)
1;
