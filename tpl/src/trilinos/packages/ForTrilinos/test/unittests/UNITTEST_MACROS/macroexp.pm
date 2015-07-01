package macroexp;

#use strict;
use warnings;

use parser;

######################################################################
# stfy: stringify the argument
######################################################################
sub strfy { return "\"" . $_[0] . "\""; }

######################################################################
# printany: print a term
######################################################################
sub printany { return "print *, " . $_[0]; }

######################################################################
# printlit: print a stringified term
######################################################################
sub printlit { return "print *, " . strfy($_[0]); }

######################################################################
# echo: print and execute statement
######################################################################
sub echo { return printlit($_[0]) . "\n" . doit($_[0]); }

######################################################################
# printany_list: print a list of terms
######################################################################
sub printany_list {
  my @terms = @_;
  my $str = "print *";
  foreach my $t (@terms) {
    $str .= ", " . $t;
  }
  return $str;
}

######################################################################
# printlit_list: print a list of stringified terms
######################################################################
sub printlit_list {
  my @terms = @_;
  my $str = "print *";
  foreach my $t (@terms) {
    $str .= ", " . strfy($t);
  }
  return $str;
}

######################################################################
# srnd: surround a term with parentheses
######################################################################
sub srnd { return "(" . $_[0] . ")"; }

######################################################################
# doit: execute statement
######################################################################
sub doit { return "$_[0]"; }

######################################################################
# assert_fail: set failure flag and print warning message
######################################################################
sub assert_fail { return $_[0] . " = .FALSE.\n" .
                         printlit_list("Assertion failed on line ", $_[1]); }

######################################################################
# print_test: verbalize the conditional we are testing
######################################################################
sub print_test {
  return printany_list(strfy("TEST: $_[0] = "), $_[0],
                    strfy(" ?$_[1]? $_[2] = "), $_[2]); }

######################################################################
# binary_test: the meat of the relational assertions below
######################################################################
sub binary_test {
  # pass string with comma-separated args, comparison op string, fail operator
  my @args = parser::break_args($_[0]);
  my $ret = print_test($args[0], $_[1], $args[1]);
  $ret .= "\nif (" . srnd($args[0]) . " " . $_[2] . " " . srnd($args[1]) . ") then";
  $ret .= "\n" . assert_fail("success", "?");
  $ret .= "\nendif";
  return $ret;
}

######################################################################
# relational assertions: can trigger unit test failure
######################################################################
sub test_equality { return binary_test($_[0], "==", ".NE."); }
sub test_inequality { return binary_test($_[0], "!=", ".EQ."); }
sub test_equiv { return binary_test($_[0], ".eqv.", ".NEQV."); }
sub test_lessequal { return binary_test($_[0], "<=", ".GT."); }


# store macro map local to this file
my %macro_map;

######################################################################
# load_macros: initialize the macro processor
######################################################################
sub load_macros {
  $macro_map{'ECHO'} = "macroexp::echo";
  $macro_map{'STRINGIFY'} = "macroexp::strfy";
  $macro_map{'PRINTANY'} = "macroexp::printany";
  $macro_map{'PRINTLIT'} = "macroexp::printlit";
  $macro_map{'TEST_EQUALITY'} = "macroexp::test_equality";
  $macro_map{'TEST_INEQUALITY'} = "macroexp::test_inequality";
  $macro_map{'TEST_EQUIV'} = "macroexp::test_equiv";
  $macro_map{'TEST_LESSEQUAL'} = "macroexp::test_lessequal";
}

######################################################################
# expand_macro: given a macro name (key) and arguments, expand it
######################################################################
sub expand_macro {
  my ($key, $args) = @_;

  if (defined($macro_map{$key})) {
    my $fn = $macro_map{$key};
    return &{$fn}($args);
  }

  return undef;
}

# module load successful (required)
1;
