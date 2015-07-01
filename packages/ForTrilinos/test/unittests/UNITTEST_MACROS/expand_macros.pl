#!/usr/bin/perl

use strict;
use warnings;

use macroexp;
use misc;
use parser;

######################################################################
# process_line: perform macro expansion and output formatting on a
#               single complete line (pre-process lines to concatenate
#               continued lines and remove the continuation
#               characters first)
######################################################################
sub process_line {
  my ($input, $nindent) = @_;
  my $result = parser::expand_line($input);
  if (defined($result)) {
    my $wrapped = misc::wrap_line($result, 4, 100);
    my $output = misc::format_block($wrapped, $nindent);
    return $output;
  }
  return "";
}

######################################################################
# main: take fortran input and expand certain macros as an
#       alternative to using the pre-processor
######################################################################

macroexp::load_macros();

my $nindent = 0;

my $prev_line = "";
while (<>) {
  my $new_line = $_;
  if (($prev_line =~ /\&\s*$/) || ($new_line =~ /^\s*\&/)) {
    # this is part of a continued line, so remove the continuation
    # characters and concatenate the lines
    $prev_line =~ s/\s*\&?\s*$//;
    $new_line =~ s/^\s*\&?\s*//;
    $prev_line .= " " . $new_line;
  } else {
    # now that we have a whole line (either standalone or processed
    # to create a full concatenated line), process it
    my $output = process_line($prev_line, $nindent);
    print $output;
    $prev_line = $new_line;
  }
}

# don't forget about the last one!!
my $output = process_line($prev_line, $nindent);
print $output;

print "\n";
