#!/usr/bin/perl -w
#
# This perl script generates an options file
# with or without comments.
#
use strict;
use strict 'refs';
use FileHandle;
use File::Copy;
use File::Path;
use File::Basename;
use Cwd;
#
my $g_use_msg =
  "\nUse: generate-opt-file.pl [-h] [-s] [-m,-i]\n".
  "  -h : Prints this help message\n".
  "  -s : Ommits comments form the generated Moocho.opt file\n".
  "  -m : Include options for MamaJama (active-set) config (default)\n".
  "  -i : Include options for IP config\n".
  "\nPerl scipt that generates a MOOCHO options file named \"Moocho.opt\" in the\n".
  "current working directory."
  ;
#
my $print_comments = 1;
my $mama_jama = 1;
for (my $i = 0; $i < @ARGV; ++$i) {
  $_ = $ARGV[$i];
  if(/^-h$/) {
    print STDERR
      $g_use_msg;
    exit(0);
  }
  elsif(/^-s$/) {
    $print_comments = 0;
  }
  elsif(/^-m$/) {
    $mama_jama = 1;
  }
  elsif(/^-i$/) {
    $mama_jama = 0;
  }
  else {
    die "The option $_ is not recognised!\n${g_use_msg}Try again!\n";
  }
}
#
open FILE_OUT, ">Moocho.opt" || die "The file Moocho.opt could not be opended for output\n";
#
my ( $this_script_name, $moocho_conf_base_dir, $this_script_suffix)
  = fileparse($0);
#
print FILE_OUT "*** Automatically generated options file\n\nbegin_options\n";
#
output_options( "$moocho_conf_base_dir/Moocho.opt.MoochoSolver", $print_comments );
output_options( "$moocho_conf_base_dir/Shared/Moocho.opt.DecompositionSystemStateStepBuilderStd", $print_comments );
if($mama_jama) {
  output_options( "$moocho_conf_base_dir/MamaJama/Moocho.opt.NLPAlgoConfigMamaJama", $print_comments );
}
else {
  output_options( "$moocho_conf_base_dir/IpConfig/Moocho.opt.NLPAlgoConfigIP", $print_comments );
}
#
print FILE_OUT "\nend_options\n";
#
close FILE_OUT;
#
print "\nGenerated output file \"Moocho.opt\"\n";

################
# Subroutines

sub output_options {
  #
  my $file_name      = shift;
  my $print_comments = shift;
  #print "\nfile_name = $file_name\n";
  #
  open FILE_IN, "<$file_name" || die "Error, could not open \'$file_name\' for input: $!";
  #
  my @file_in_array = <FILE_IN>;
  #
  if($print_comments) {
	foreach(@file_in_array) {
	  print FILE_OUT $_;
	}
  }
  else {
	my $i = 0;
	# Print the first line
	print FILE_OUT "\n", $file_in_array[$i];
	++$i;
	# Print the option groups only along with their (commented out) options
	for( ; $i < scalar(@file_in_array); ++$i ) {
	  $_ = $file_in_array[$i];
	  if(/^options_group/) {
		print FILE_OUT "\n", $_;
		++$i;
		# Print out all of the (commented out) options
		for( ; $i < scalar(@file_in_array); ++$i ) {
		  $_ = $file_in_array[$i];
		  if(/^\*/) {
			print FILE_OUT $_;
		  }
		  elsif(/^\}/) {
			print FILE_OUT $_;
			++$i;
			last;
		  }
		}
	  }
	}
  }
}
