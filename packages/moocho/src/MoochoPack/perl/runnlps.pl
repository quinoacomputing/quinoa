#!/usr/bin/perl -w
#
use strict;
use lib "$ENV{MOOCHO_BASE_DIR}/Moocho/core/MoochoPack/perl";
#
my $g_use_msg =
  "Use: runnlps.pl [-h] [-i] [-k0,-k1] [-e executable] [-a \"arguments\"]\n".
  "                [-s setup_script] [-p post_script] [-b base_opt_file]\n".
  "                [-v vary_opt_file] [-o output_file]\n";
#
my $g_use_msg_opts =
  $g_use_msg.
  "options:\n".
  "  -h  : Prints this help message and the program quits\n".
  "  -i  : Print iteration output (default is no iteration output)\n".
  "  -k0 : Don't keep any of the output files or directories (default)\n".
  "  -k1 : Keep the output files and directories\n".
  "  -e executable\n".
  "     : Gives the absolute path of the executable to run (with arguments).\n".
  "      (default is to use the executable with the name \'solve_nlp\').\n".
  "  -a \"arguments\"\n".
  "      : Arguments to pass on to the executable. (default is no arguments).\n".
  "  -s setup_script\n".
  "      : Script to be run in each directory just before the executable is\n".
  "        called (default is to run no script before hand).\n".
  "  -p post_script\n".
  "      : Script to be run just after the executable is called in each directory\n".
  "        (default is to run no script after at all).\n".
  "  -b base_opt_file\n".
  "      : Gives the path of the options file with the base options to be used.\n".
  "        (default is to use the \'Moocho.opt\' file in the current directory. If\n".
  "          this file does not exist then no base options will be used).\n".
  "  -v vary_opt_file\n".
  "      : Gives the path of the file containing the options to be varied (default\n".
  "        is to use the file \'vary.opt\' in the current directory.  If this file\n".
  "        does not exist then no options are varied)\n".
  "  -o output_file\n".
  "      : This is the file to which output is printed (default is to STDOUT)\n";
#
# This is a perl program for running an NLP while
# varying a set of options for MOOCHO.  This program will print
# lots of nice output and give lots of nice statistics.  The program will return 0 if
# none of the NLP runs throws an exception.
#
# Now the format of the above input files will be defined:
#
#   base_opt_file : Options will be read in for all the options between begin_options and end_options:
#
#       ...
#       begin_options
#       ...
#       end_options
#       ...
#
#   vary_opt_file : Options to vary will be read in between begin_vary_options and end_vary_options:
#
#       ...
#       begin_vary_options
#       ...
#       begin_options_group:
#       NLPAlgoConfigMamaJama
#           begin_option:
#           qp_solver
#               QPKWIK
#                   QPKWIK
#      *         QPOPT
#      *             QPOPT
#               QPSCHUR
#                   QPSCHUR
#           end_option
#       ...
#       end_options_group
#       ...
#       end_vary_options
#       ...
#

use strict 'refs';
use FileHandle;
use File::Copy;
use File::Path;
use Cwd;

use MoochoVaryOptions;

use constant TRUE  => 1;
use constant FALSE => 0;

##########################
#
# NLP runner class
#
package NLPRunner;

use strict;
use strict 'refs';

# Object member names
use constant SETUP_SCRIPT       => "setup_script";
use constant POST_SCRIPT        => "post_script";
use constant EXECUTABLE         => "executable";
use constant EXEC_ARGUMENTS     => "exec_arguments";
use constant PRINT_ITER_SUMMARY => "print_iter_summary";

##########################
# Class methods

# Constructor
sub new {
  # Arguments
  my $class           = shift; # Class name (required)
  my $setup_script    = shift; # Script to run before running executable (required)
  my $post_script     = shift; # Script to run after running executable (required)
  my $executable      = shift; # Executable to be run (required)
  my $exec_arguments  = shift; # Arguments to pass on to the executable
  my $print_iter_sum  = shift; # Print the iteration summary to STDOUT or not (optional)
  #
  my $self = {
			  SETUP_SCRIPT()        => $setup_script
			  ,POST_SCRIPT()        => $post_script
			  ,EXECUTABLE()         => $executable
			  ,EXEC_ARGUMENTS()     => $exec_arguments
			  ,PRINT_ITER_SUMMARY() => $print_iter_sum
			 };
  bless $self, $class;
  return $self;
}

#########################
# Object methods

#
# Setup the current directory and run MOOCHO on this NLP
#
sub run_nlp {
  # Arguments
  my $self = shift;
  #
  my $setup_script    = $self->{SETUP_SCRIPT()};
  my $post_script     = $self->{POST_SCRIPT()};
  my $executable      = $self->{EXECUTABLE()};
  my $exec_arguments  = $self->{EXEC_ARGUMENTS()};
  my $print_iter_sum  = $self->{PRINT_ITER_SUMMARY()};
  # Run the setup script
  my $run_cmnd;
  if(defined($setup_script)) {
	$run_cmnd = $setup_script;
	system($run_cmnd);
  }
  # Run the NLP
  $run_cmnd = $executable;
  if(defined($exec_arguments)) {
	$run_cmnd .= (" " . $exec_arguments);
  }
  if( !$print_iter_sum ) {
	$run_cmnd .= " 1> console.out 2> console.out";
  }
  system($run_cmnd);
  # Run the post script
  if(defined($post_script)) {
	$run_cmnd = $post_script;
	system($run_cmnd);
  }
}

package main; # end class NLPRunner

###################
# Executable code #
###################

#
# Directory locations
#

my $g_basedir = Cwd::cwd();

#
# Get the options
#

if( scalar(@ARGV) == 0 ) {
  print STDERR
	$g_use_msg;
  exit(0);
}

# Print iteration summary or not
my $g_print_iter_summ = FALSE();
# Remove output files or not
my $g_remove_output_files = TRUE();
# Fully qualified name of the executable with arguments
# (if UNDEF then no file was specified)
my $g_executable_str = $g_basedir . "/solve_nlp";
# Arguments to pass on to executable
my $g_exec_arguments_str;
# Fully qualified name of the setup script with arguments
# (if UNDEF then no scrip was specified)
my $g_setup_script_str;
# Fully qualified name of the post script with arguments
# (if UNDEF then no script was specified)
my $g_post_script_str;
# Fully qualified name of the base options file
# (if UNDEF then no file was specified)
my $g_default_options_file_name;
# Fully qualified name of the vary options file
# (if UNDEF then no file was specified)
my $g_vary_options_file_name;
# Output file
# (If UNDEF then output is sent to STDOUT )
my $g_output_file_name;

for (my $i = 0; $i < @ARGV; ++$i) {
  $_ = $ARGV[$i];
  if(/^-h$/) {
	print STDERR
	  $g_use_msg_opts;
	exit(0);
  }
  elsif(/^-i$/) {
	$g_print_iter_summ = TRUE();
  }
  elsif(/^-k0$/) {
	$g_remove_output_files       = TRUE();
  }
  elsif(/^-k1$/) {
	$g_remove_output_files       = FALSE();
  }
  elsif(/^-e$/) {
	$g_executable_str = $ARGV[++$i];
  }
  elsif(/^-a$/) {
	$g_exec_arguments_str = $ARGV[++$i];
  }
  elsif(/^-s$/) {
	$g_setup_script_str = $ARGV[++$i];
  }
  elsif(/^-p$/) {
	$g_post_script_str = $ARGV[++$i];
  }
  elsif(/^-b$/) {
	$g_default_options_file_name = $ARGV[++$i];
  }
  elsif(/^-v$/) {
	$g_vary_options_file_name = $ARGV[++$i];
  }
  elsif(/^-o$/) {
	$g_output_file_name = $ARGV[++$i];
  }
  else {
	die "The option $_ is not recognised!  Try again! : $!";
  }
}

#
# Open the files and get the info you need
#

# Get the default options
my $g_default_options_fh;  # Will be a reference to a file handle!
if( $g_default_options_file_name ) {
  if(!($g_default_options_fh = new FileHandle "<$g_default_options_file_name") )
	{
	  die "The default options file $g_default_options_file_name could not be opened!: $!";
	}
}
else {
  if(!($g_default_options_fh = new FileHandle "<Moocho.opt") )
	{
	  print STDERR "Warning, the file Moocho.opt could not be found!  Not using any default options!\n";
	}
}

# Get a file handle to the vary_opt file
my $g_vary_options_fh; # Will be a FileHandle to the vary options file
if( $g_default_options_file_name ) {
  if(!($g_vary_options_fh = new FileHandle "<$g_vary_options_file_name") )
	{
	  die "Error, the vary options file $g_vary_options_file_name could not be opened: $!";
	}
}
else {
  if(!($g_vary_options_fh = new FileHandle "<vary.opt") )
	{
	  print STDERR "Warning, the file vary.opt could not be found!  Will not vary any options!\n";
	}
}

# Get a handle to the output file
my $g_output_fh; # A FileHandle object or a simple typeglob to *STDOUT
if( $g_output_file_name ) {
  if(!($g_output_fh = new FileHandle ">$g_output_file_name") )
	{
	  die "Error, the output file $g_output_file_name could not be opened: $!";
	}
}
else {
  $g_output_fh = *STDOUT;
}

#
# Print output header
#
print $g_output_fh
  "\n*****************************************",
  "\n*** Output summary for solving an NLP ***",
  "\n*** for a various set of options      ***",
  "\n*****************************************\n";

#
# Readin and print the default options
#
my @g_default_options; # NULL by default
if( defined $g_default_options_fh ) {
  @g_default_options
	= readin_options_lines( $g_default_options_fh, "begin_options", "end_options" );
}
print $g_output_fh
  "\n*********************************",
  "\n*** Default options to be used\n",
  "\n", join("\n",@g_default_options), "\n";

#
# Create the RSQPppVaryOptions object and initialize it
#
my $g_vary_options
  = MoochoVaryOptions->new(
						   "Moocho.opt"          # Options file name to create
						   ,\@g_default_options  # Default options
						   ,$g_vary_options_fh   # Filehandle to get options to vary
						 );
# Print the legend for the options that we will vary
$g_vary_options->print_vary_opts_abbr($g_output_fh);

#
# Create the runs directory
#
mkdir("runs",0777) || die "Error, Could not make directory \'runs\': $!";
chdir("runs") || die "What: $!";

# Run with the various options!
print $g_output_fh
  "\n**************************************",
  "\n*** Summary of running MOOCHO on NLP\n";
$g_vary_options->run_nlps(
						  NLPRunner->new(
										 $g_setup_script_str
										 ,$g_post_script_str
										 ,$g_executable_str
										 ,$g_exec_arguments_str
										 ,$g_print_iter_summ )
						  ,$g_output_fh
						 );
#
# Move out and remove the runs directory (if asked to)
#
chdir("..") || die "What: $!";
if( $g_remove_output_files ) {
  rmtree( "runs", 0, 0 ) || die "Error, Could not remove the directory \'runs\': $!";
}

#
# Print the final statistics for the runs
#
$g_vary_options->print_all_opt_values_stats($g_output_fh);

exit ($g_vary_options->total_num_except() ? -1 : 0);
