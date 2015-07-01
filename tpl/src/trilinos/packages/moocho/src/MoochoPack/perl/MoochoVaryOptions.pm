package      MoochoVaryOptions;
require      Exporter;
@ISA       = qw(Exporter);
@EXPORT    = qw(readin_options_lines);
@EXPORT_OK = qw();

use strict;

##################################################
#
# Class for varying MOOCHO options for an NLP
#
# ToDo: Finish documentation!
#

########################
# External subroutines #
########################

#
# Read in the options by line, taking off the newline
# and ignoring lines that begin with an '*'.
#
sub readin_options_lines {
  # Input arguments
  my $fh          = shift; # FileHandle reference
  my $begin_token = shift; # Beginning token name
  my $end_token   = shift; # Ending token name
  # 
  my @options_lines;
  MAINLOOP: while( <$fh> ) {
	chomp;
	if(/$begin_token/) {
	  # Read in the default options
	  while( <$fh> ) {
		chomp;
		if(/$end_token/) {
		  last MAINLOOP;
		}
		elsif( !/^\s*\*/ && !/^\s?$/ ) {
		  @options_lines = (@options_lines,$_);
		}
	  }	
	}
	next;
  }
  return @options_lines;
}

#####################################################
# Class data

#
# Define data structures for summary file outputing
#

# Define formating for output

my $stat_len = 6;
my $var_opt_len = 30;
my $flt_len = 12;
my $int_len = 5;
my $num_spaces = 3;

my $spc				= " " x $num_spaces;
my $tstat_len 		= $stat_len		+ $num_spaces;
my $tvar_opt_len	= $var_opt_len	+ $num_spaces;
my $tflt_len		= $flt_len		+ $num_spaces;
my $tint_len		= $int_len		+ $num_spaces;
	
my $fprec = 5;
	
my $g_summary_cols_header_format =
	"%-${stat_len}s".		# status
	"%-${tvar_opt_len}s".	# varied options
	"%-${tint_len}s".		# #iter
	"%-${tint_len}s".		# #func
	"%-${tint_len}s".		# #grad
	"%-${tflt_len}s".		# CPU(sec)
	"%-${tflt_len}s".		# f(x)
	"%-${tflt_len}s".		# ||c(x)||s
	"%-${tflt_len}s".		# ||rGL||s
	"%-${tint_len}s".		# nact
	"%-${tint_len}s".		# #BC
	"%-${tint_len}s\n";		# #QN

my $g_summary_cols_lines_format =
	"%-${stat_len}s".			# status
	$spc."%-${var_opt_len}s".	# varied options
	"%${tint_len}s".			# #iter
	"%${tint_len}s".			# #func
	"%${tint_len}s".			# #grad
	"%${tflt_len}e".			# CPU(sec)
	"%${tflt_len}e".			# f(x)
	"%${tflt_len}e".			# ||c(x)||s
	"%${tflt_len}e".			# kkt_error
	"%${tint_len}s".			# nact
	"%${tint_len}s".			# #BC
	"%${tint_len}s\n";			# #QN
		
my $g_summary_cols_header_labels =
	sprintf( $g_summary_cols_header_format,
		"status",
		$spc."Varied Options",
		$spc."#iter",
		$spc."#func",
		$spc."#grad",
		$spc."CPU (sec)",
		$spc."f(x)",
		$spc."||c(x)||s",
		$spc."||rGL||s",
		$spc."nact",
		$spc."#BC",
		$spc."#QN"
	);

my $g_summary_cols_lines_sep =
	sprintf( $g_summary_cols_header_format,
		"-"x$stat_len,				# status
		$spc."-"x$var_opt_len,		# varied options
		$spc."-"x$int_len,			# #iter
		$spc."-"x$int_len,			# #func
		$spc."-"x$int_len,			# #grad
		$spc."-"x$flt_len,			# CPU(sec)
		$spc."-"x$flt_len,			# f(x*)
		$spc."-"x$flt_len,			# ||c(x*)||
		$spc."-"x$flt_len,			# kkt_error
		$spc."-"x$int_len,			# nact*
		$spc."-"x$int_len,			# #BC
		$spc."-"x$int_len			# #QN
	);
	
# Define array indexing for storing and manipulating summary line
# for each NLP

my $g_status_i			= 0;
my $g_varied_options_i	= 1;
my $g_niter_i			= 2;
my $g_nfunc_i			= 3;	
my $g_ngrad_i			= 4;
my $g_CPU_i				= 5;
my $g_fx_i				= 6;
my $g_nrm_cx_i			= 7;
my $g_kkt_error_i		= 8;
my $g_nact_i			= 9;
my $g_nBC_i				= 10;
my $g_nQN_i				= 11;

my $g_real_big   = 1e+10;
my $g_real_small = 1e-50;

# Used to index into arrays for summaries of each run for an NLP

use constant SUM_STATUS     => 0;
use constant SUM_NITER      => 1;
use constant SUM_NFUNC      => 2;
use constant SUM_CPU        => 3;
use constant SUM_OPT_VALUES => 4;

################################################################
#
# Object instance data (hash)
#
# $self (ref to hash)
#  |
#  |- {OPT_FILE_NAME} (string)
#  |
#  |- {DEFAULT_OPTS}  (ref to array)
#  |
#  |- {VARY_OPTS}     (ref to hash)
#  |  |
#  |  |- {V_OPT_GRP}          (ref to array)
#  |  |
#  |  |- {V_FIRST_OPT_IN_GRP} (ref to array)
#  |  |
#  |  |- {V_OPT}              (ref to array)
#  |  |
#  |  |- {V_FIRST_OPT_VALUE}  (ref to array)
#  |  |
#  |  |- {V_OPT_VALUE}        (ref to array)
#  |  |
#  |   - {V_OPT_VALUE_ABBR}   (ref to array)
#  |
#  |- {TOTAL_NLP_RUNS}    (integer)
#  |
#  |- {TOTAL_NLPS_SOLVED} (integer)
#  |
#  |- {TOTAL_NLPS_EXCEPT} (integer)
#  |
#   - {OPT_VAL_STATS} (ref to hash)
#     |
#     |- {STATS_MIN_IT}       (ref to array)
#     |
#     |- {STATS_MIN_FN}       (ref to array)
#     |
#     |- {STATS_MIN_T}        (ref to array)
#     |
#     |- {STATS_MAX_IT}       (ref to array)
#     |
#      - {STATS_EXCEPT}       (ref to array)
#
use constant OPT_FILE_NAME => "opt_file_name";
use constant DEFAULT_OPTS  => "default_opts";
use constant VARY_OPTS     => "vary_opts";

#
# v_opt_grp[i] is the name of the ith (0-based) option group
# 0 <= i <= scalar(@v_opt_grp)
#
#my @v_opt_grp;
use constant V_OPT_GRP => "v_opt_grp";

#
# $g_v_opt[g_v_first_opt_in_grp[i]+j] is the jth (0-based) option in the
# option group v_opt_grp[i].
# 0<= j <= g_v_first_opt_in_grp[i] - g_v_first_opt_in_grp[i+1]
# scalar(@g_v_first_opt_in_grp) == scalar(@v_opt_grp)+1
# scalar(@g_v_opt) == scalar(g_v_first_opt_in_grp(scalar(@v_opt_grp))
#
#my @v_first_opt_in_grp;
use constant V_FIRST_OPT_IN_GRP => "v_first_opt_in_grp";
#my @v_opt;
use constant V_OPT => "v_opt";

#
# v_opt_value[g_v_first_opt_value[v_first_opt_in_group[i]+j]+k] and
# v_opt_value_abbr[g_v_first_opt_value[v_first_opt_in_group[i]+j]+k] are
# the kth (0-based) option value and abbreviation respectively
# of the jth option v_opt[j] of the ith options group v_opt_grp[i].
# 0 <= k <= g_v_first_opt_value[p+j] - g_v_first_opt_value[p+j+1]
#		for: p = v_first_opt_in_group[i]
# scalar(@g_v_first_opt_value) == scalar(@v_opt_value)+1
#
#my @v_first_opt_value;
use constant V_FIRST_OPT_VALUE => "v_first_opt_value";
#my @v_opt_value;
use constant V_OPT_VALUE       => "v_opt_value";
#my @v_opt_value_abbr;
use constant V_OPT_VALUE_ABBR  => "v_opt_value_abbr";

# The following example shows how these things are stored for the
# given example input:
#
#begin_vary_options
#	begin_options_group:
#	NLPAlgoConfigMamaJama
#		begin_option:
#		qp_solver
#			QPKWIK
#				QPKWIK
#			QPSCPD
#				QPSCPD
#		end_option
#		begin_option:
#		line_search
#			DIRECT
#
#			2ND_ORDER_CORRECT
#				2nd
#			WATCHDOG
#				watch
#		end_option
#	end_options_group
#	begin_options_group:
#	DirectLineSearchArmQuad
#		begin_option:
#		min_frac_step
#			0.1
#				minf=0.1
#			0.2
#				minf=0.2
#		end_option
#	end_options_group
#end_vary_options
#
#		v_opt_grp					v_first_opt_in_grp
#	0	NLPAlgoConfigMamaJama		0
#	1	DirectLineSearchArmQuad		2
#	2								3
#
#		v_opt						v_first_opt_value
#	0	qp_solver					0
#	1	line_search					2
#	2	min_frac_step				5
#	3								7
#
#		v_opt_value					v_opt_value_abbr
#	0	QPKWIK						QPKWIK
#	1	QPSCPD						QPSCPD
#	2	DIRECT						
#	3	2ND_ORDER_CORRECT			2nd
#	4	WATCHDOG					watch
#	5	0.1							minf=0.1
#	6	0.2							minf=0.2

#
# Total number of problems attempted
#
use constant TOTAL_NLP_RUNS => "total_nlp_runs";

#
# Total number of problems solved
#
use constant TOTAL_NLPS_SOLVED => "total_nlps_solved";

#
# Total number of problems resulting in exceptions
#
use constant TOTAL_NLPS_EXCEPT => "total_nlps_except";

#
# Declare a hash of storage arrays for the cummutaltive statistics
# for each opion value.
#
use constant OPT_VAL_STATS => "opt_val_stats";
#
# Statistics are keep track for each option value the minimum number of
# iterations (min_it), minimum number of function evaluations (min_fn)
# and minimum CPU time (min_t).
#
# The members of each one of these arrays belongs to a specific option value abs_j.
#
# min_it :	Each element in this array min_it[abs_j] is is a reference to
#           an array for statistis for problems solved.  In this array
#			min_it[abs_j][0] is the number of runs involving this option that
#			where the winner or equal to the winner (niter == min_itr).
#			min_it[abs_j][i] (for i = 1...10) is the number of runs that were
#			within 10*i percent of the winner
#			(i.e. niter <= min_it + (0.1 * i) * min_it).
#			min_it[abs_j][11] is the number of runs that were more than 100% worst than
#			the winner (i.e. niter > min_it + (0.1 * i) * min_it).
#
use constant STATS_MIN_IT => "min_it";
#
# min_fn :	Same array of statistics as min_it except for minumum number
#			of function evaluations.
#
use constant STATS_MIN_FN => "min_fn";
#
# min_t :	Same array of statistics as min_it except for minumum CPU time.
#
use constant STATS_MIN_T => "min_t";
#
# max_it :	max_it[abs_j] is the number of runs where the maximum number
#			of iterations was exceeded.
#
use constant STATS_MAX_IT => "max_it";
#
# except :  except[abs_j] is the number of runs where an exception was thrown.
#
use constant STATS_EXCEPT => "except";

################################################################
# Class methods

# Constructor
sub new {
  # Arguments
  my $class         = shift; # Class name (required)
  my $opt_file_name = shift; # name of the options file to create (required)
  my $default_opts  = shift; # Reference to array of default options (required)
  my $vary_opts_fh  = shift; # Filehandle for reading in what options to vary (optional)
  #
  defined($class)
	|| die "Error, class name required: $!";
  defined($opt_file_name)
	|| die "Error, options file name required: $!";
  (ref($default_opts) eq "ARRAY")
	|| die "Error, reference to array of default options required: $!";
  #
  my $self = {
			  OPT_FILE_NAME()  => $opt_file_name,
			  DEFAULT_OPTS()   => [],
			  VARY_OPTS()      => {
								   V_OPT_GRP()          => [],
								   V_FIRST_OPT_IN_GRP() => [],
								   V_OPT()              => [],
								   V_FIRST_OPT_VALUE()  => [],
								   V_OPT_VALUE()        => [],
								   V_OPT_VALUE_ABBR()   => []
								  },
			  TOTAL_NLP_RUNS()   => 0,
			  TOTAL_NLPS_SOLVED()=> 0,
			  TOTAL_NLPS_EXCEPT()=> 0,
			  OPT_VAL_STATS()  => {
								  STATS_MIN_IT()       => [],
								  STATS_MIN_FN()       => [],
								  STATS_MIN_T()        => [],
								  STATS_MAX_IT()       => [],
								  STATS_EXCEPT()       => []
								 }
			 };
  @{$self->{DEFAULT_OPTS()}} = @$default_opts;
  bless $self, $class;
  if( $vary_opts_fh ) { 
	$self->readin_opts_to_vary($vary_opts_fh) 
  };
  return $self;
}

###############################################################
# Public object methods

#
# Read in options to vary from an filehandle
#
sub readin_opts_to_vary {
  # Arguments
  my $self = shift;
  my $in   = shift;  # File handle to extract options to vary from
  #
  my $v_opt_grp          = $self->{VARY_OPTS()}->{V_OPT_GRP()};         # ref to array
  my $v_first_opt_in_grp = $self->{VARY_OPTS()}->{V_FIRST_OPT_IN_GRP()}; # ref to array
  my $v_opt              = $self->{VARY_OPTS()}->{V_OPT()};             # ref to array
  my $v_first_opt_value  = $self->{VARY_OPTS()}->{V_FIRST_OPT_VALUE()}; # ref to array
  my $v_opt_value        = $self->{VARY_OPTS()}->{V_OPT_VALUE()};       # ref to array
  my $v_opt_value_abbr   = $self->{VARY_OPTS()}->{V_OPT_VALUE_ABBR()};  # ref to array
  #
  my $num_opt_grp    = 0;
  my $num_opt        = 0;
  my $num_opt_values = 0;
  my $end_opt_group  = 0;
  # Read in the options values
  while( <$in> ) {
	chomp;
#	print "find begin_vary_options ->", $_, "\n";
	if( /begin_vary_options/ && !/^\*/ ) {
	  # Read in the options groups
	  while( <$in> ) {
		chomp;
#		print "find begin_options_group: ->", $_, "\n";
		if( /begin_options_group:/ && !/^\*/ ) {
		  # Read in the options group name
		  $v_opt_grp->[$num_opt_grp] = r_ws(scalar(<$in>));
#		  print $v_opt_grp->[$num_opt_grp], "\n";
		  # Set the location of the first option in this group
		  $v_first_opt_in_grp->[$num_opt_grp++] = $num_opt;
		  # Read in the options and values
		  while( <$in> ) {
			chomp;
#			print "find begin_option: ->", $_, "\n";
			if( /begin_option:/ && !/^\*/ ) {
			  # Read in the option name
			  $v_opt->[$num_opt] = r_ws(scalar(<$in>));
#			  print $v_opt->[$num_opt], "\n";
			  # Set the location of the first value for this option
			  $v_first_opt_value->[$num_opt++] = $num_opt_values;
			  # Read in the options values
			  my $loc_num_values = 0;
			  while( <$in> ) {
				chomp;
#				print "read option value ->", $_, "\n";
				if( /end_option/ ) {
				  last;
				}
				elsif( /^\*/ ) {
								# Skip since this is a comment line
				}
				else {
								# Read in the option value
				  $v_opt_value->[$num_opt_values] = r_ws($_);
								# Read in the option value abbreviation.
				  $v_opt_value_abbr->[$num_opt_values++]
					= r_ws(scalar(<$in>));
				  $loc_num_values++;
				}
			  }					# end while looking for option values
			  if( $loc_num_values == 0 ) {
#				print "Error, there are no listed values for the option $v_opt->[$num_opt-1]\n";
				exit(-1);
			  }
			}
			elsif( /end_options_group/  && !/^\*/ ) {
			  last;
			}
		  }						# end while looking for begin_options:
		}						# end if begin_options:
		elsif( /end_vary_options/  && !/^\*/ ) {
		  last;
		}
	  }							# end while looking for begin_options_group:
	  last;
	}
  }								# end while looking for begin_vary_options:
  # Cap off the ends
  $v_first_opt_in_grp->[$num_opt_grp] = $num_opt;
  $v_first_opt_value->[$num_opt] = $num_opt_values;
  # Reset the statistics
  $self->reset_stats();
}

#
# Print out summary header for options to vary
#
sub print_vary_opts_abbr {
  # Arguments
  my $self = shift;
  my $out  = shift;  # File handle to print options to
  #
  my $v_opt_grp          = $self->{VARY_OPTS()}->{V_OPT_GRP()};         # ref to array
  my $v_first_opt_in_grp = $self->{VARY_OPTS()}->{V_FIRST_OPT_IN_GRP()}; # ref to array
  my $v_opt              = $self->{VARY_OPTS()}->{V_OPT()};             # ref to array
  my $v_first_opt_value  = $self->{VARY_OPTS()}->{V_FIRST_OPT_VALUE()}; # ref to array
  my $v_opt_value        = $self->{VARY_OPTS()}->{V_OPT_VALUE()};       # ref to array
  my $v_opt_value_abbr   = $self->{VARY_OPTS()}->{V_OPT_VALUE_ABBR()};  # ref to array
  #
  my $abbr_len = 12;
  my $opt_len  = 52;
  #
  print $out
	"\n*********************************************",
	"\n*** Abbreviations for options to be varied\n",
  	
  	"\nabbreviation",
  	"   option\n",
  	"-" x $abbr_len,
  	"   " . "-" x $opt_len . "\n";
    
    for( my $i = 0; $i < scalar(@$v_opt_grp); $i++ ) {
	  my $opt_grp = $v_opt_grp->[$i];
	  my $abs_j;
	  for( $abs_j = $v_first_opt_in_grp->[$i];
		   $abs_j < $v_first_opt_in_grp->[$i+1];
		   $abs_j++								)
		{
		  my $opt = $v_opt->[$abs_j];
		  my $abs_k;
		  for( $abs_k = $v_first_opt_value->[$abs_j];
			   $abs_k < $v_first_opt_value->[$abs_j+1];
			   $abs_k++								)
			{
			  my $val		= $v_opt_value->[$abs_k];
			  my $val_abbr	= $v_opt_value_abbr->[$abs_k];
			  printf $out
				"%${abbr_len}s", $val_abbr;
			  printf $out
				"   ";
			  printf $out
				"%-${opt_len}s", $opt_grp."::".$opt." = ".$val;
			  printf $out
				"\n";
			}
  	  }
    }
}

#
# Reset the statistics
#
sub reset_stats {
  # Arguments
  my $self = shift;
  #
  my $v_opt_value        = $self->{VARY_OPTS()}->{V_OPT_VALUE()};       # ref to array
  #
  my $min_it             = $self->{OPT_VAL_STATS()}->{STATS_MIN_IT()};  # ref to array
  my $min_fn             = $self->{OPT_VAL_STATS()}->{STATS_MIN_FN()};  # ref to array
  my $min_t              = $self->{OPT_VAL_STATS()}->{STATS_MIN_T()};   # ref to array
  my $max_it             = $self->{OPT_VAL_STATS()}->{STATS_MAX_IT()};  # ref to array
  my $except             = $self->{OPT_VAL_STATS()}->{STATS_EXCEPT()};  # ref to array
  #
  $self->{TOTAL_NLP_RUNS()}    = 0;
  $self->{TOTAL_NLPS_SOLVED()} = 0;
  $self->{TOTAL_NLPS_EXCEPT()} = 0;
  my @zero_array;
  for( my $i = 0; $i <= 11; ++$i ) { 
	$zero_array[$i] = 0;
  }
  for( my $j = 0; $j < @$v_opt_value; ++$j ) {
	$min_it->[$j] = [@zero_array];
	$min_fn->[$j] = [@zero_array];
	$min_t->[$j]  = [@zero_array];
	$max_it->[$j] = 0;
	$except->[$j] = 0;
  }
}

#
# Called by the client to vary the options and run the NLPs
#
sub run_nlps {
  # Arguments
  my $self       = shift;
  my $nlp_runner = shift; # Reference to an object that will perform
                          # additional setup and run the NLP.  This object must
                          # support the method $nlp_runner->run_nlp().  This
                          # method should perform any file setup needed to run the
                          # NLP and then run it, creating the output directories.
  my $out        = shift; # Output file handle
  #
  my @v_set_option_value; # Will be used to hold the current options value
  my %best_run_stats = (  # Keep track of the runs with the best statistics
						STATS_MIN_IT()       => [],
						STATS_MIN_FN()       => [],
						STATS_MIN_T()        => []
					   );
  $best_run_stats{STATS_MIN_IT()}->[$g_niter_i]	= $g_real_big; # Reinitialize memory for the best of the set of runs
  $best_run_stats{STATS_MIN_FN()}->[$g_nfunc_i]	= $g_real_big;
  $best_run_stats{STATS_MIN_T()}->[$g_CPU_i]	= $g_real_big;
  my @runs_summaries; # Keep track of statistics for all of the runs
  @runs_summaries[SUM_STATUS()]      = []; # Will have one array element for each run
  @runs_summaries[SUM_NITER()]       = [];
  @runs_summaries[SUM_NFUNC()]       = [];
  @runs_summaries[SUM_CPU()]         = [];
  @runs_summaries[SUM_OPT_VALUES()]  = [];
  my $new_nlp = 1;  # Set to false once the first NLP is solved
  # Run all the NLPs and vary the options
  $self->_vary_options(0,$nlp_runner,\@v_set_option_value
					   ,\%best_run_stats,\@runs_summaries,\$new_nlp,$out);
  # Print the bestof lines
  $self->_add_best_runs_lines(\@v_set_option_value,\%best_run_stats
							  ,\@runs_summaries,$out);
};

#
# Print the performance statistics for the option values.
#
sub print_all_opt_values_stats {
  # Arguments
  my $self       = shift;
  my $out        = shift; # Output file handle
  #
  my $v_opt_grp          = $self->{VARY_OPTS()}->{V_OPT_GRP()};         # ref to array
  my $v_first_opt_in_grp = $self->{VARY_OPTS()}->{V_FIRST_OPT_IN_GRP()}; # ref to array
  my $v_opt              = $self->{VARY_OPTS()}->{V_OPT()};             # ref to array
  my $v_first_opt_value  = $self->{VARY_OPTS()}->{V_FIRST_OPT_VALUE()}; # ref to array
  my $v_opt_value        = $self->{VARY_OPTS()}->{V_OPT_VALUE()};       # ref to array
  my $v_opt_value_abbr   = $self->{VARY_OPTS()}->{V_OPT_VALUE_ABBR()};  # ref to array
  #
  print $out
	"\n\n********************************************************",
	"\n*** Summary of performance of varied options\n",
	"\nTotal number of NLP runs                 = ", $self->{TOTAL_NLP_RUNS()}, "\n",	
	"\nTotal number of NLP runs solved          = ", $self->{TOTAL_NLPS_SOLVED()}, "\n",	
	"\nTotal number of NLP runs with exceptions = ", $self->{TOTAL_NLPS_EXCEPT()}, "\n";	
  my ($opt_grp_i, $opt_abs_j, $opt_val_abs_k);
  for( $opt_grp_i = 0; $opt_grp_i < @$v_opt_grp; $opt_grp_i++ ) {
	my $opt_grp_name = $v_opt_grp->[$opt_grp_i];
	for(	$opt_abs_j = $v_first_opt_in_grp->[$opt_grp_i];
			$opt_abs_j < $v_first_opt_in_grp->[$opt_grp_i+1];
			$opt_abs_j++ )
	  {
		my $opt_name = $v_opt->[$opt_abs_j];
		print $out
		  "\n\nOption : ", $opt_grp_name , "::" , $opt_name, "\n";
		$self->_print_opt_values_stats_max_it_except($opt_abs_j,$out);
		$self->_print_opt_values_stats($opt_abs_j,STATS_MIN_FN(),$out);
		$self->_print_opt_values_stats($opt_abs_j,STATS_MIN_IT(),$out);
		$self->_print_opt_values_stats($opt_abs_j,STATS_MIN_T(),$out);
	  }
  }
}

#
# Return the total number of runs
#
sub total_num_runs {
  my $self       = shift;
  return $self->{TOTAL_NLP_RUNS()};
}

#
# Return the number of runs solved
#
sub total_num_solved {
  my $self       = shift;
  return $self->{TOTAL_NLPS_SOLVED()};
}

#
# Return the number of runs with exceptions
#
sub total_num_except {
  my $self       = shift;
  return $self->{TOTAL_NLPS_EXCEPT()};
}

##########################################
# Private object methods

#
# Recusively vary the options and run the tests
#
# Returns true if there are more options to be varied, false otherwise
#
sub _vary_options {
  # Arguments
  my $self               = shift;
  my $abs_opt_j          = shift; # The absolution option number
  my $nlp_runner         = shift; # Reference to an object that will perform
                                  # additional setup and run the NLP.
  my $v_set_option_value = shift; # reference to array for current options values
  my $best_run_stats     = shift; # Reference to a hash of arrays (see run_nlps(...))
  my $runs_summaries     = shift; # Refererence to array of arrays (see run_nlps(...))
  my $new_nlp            = shift; # Reference to scalar set to true once the first NLP is solved
  my $out                = shift; # Output file handle
  #
  my $default_opts       = $self->{DEFAULT_OPTS()};                     # ref to array
  my $v_opt_grp          = $self->{VARY_OPTS()}->{V_OPT_GRP()};         # ref to array
  my $v_first_opt_in_grp = $self->{VARY_OPTS()}->{V_FIRST_OPT_IN_GRP()};# ref to array
  my $v_opt              = $self->{VARY_OPTS()}->{V_OPT()};             # ref to array
  my $v_first_opt_value  = $self->{VARY_OPTS()}->{V_FIRST_OPT_VALUE()}; # ref to array
  my $v_opt_value        = $self->{VARY_OPTS()}->{V_OPT_VALUE()};       # ref to array
  my $v_opt_value_abbr   = $self->{VARY_OPTS()}->{V_OPT_VALUE_ABBR()};  # ref to array
  #
  if( $abs_opt_j < scalar(@$v_opt) ) {
	# Vary the values for this option
	my $num_values = $v_first_opt_value->[$abs_opt_j+1]
	  - $v_first_opt_value->[$abs_opt_j];
	my $k;
	for( $k = 0; $k < $num_values; $k++ ) {
	  $v_set_option_value->[$abs_opt_j]
		= $v_first_opt_value->[$abs_opt_j] + $k;
	  # Create the directory for this option and move in
	  my $optvaldir = $v_opt->[$abs_opt_j]
		. "=" . $v_opt_value->[$v_first_opt_value->[$abs_opt_j]+$k];
	  mkdir( $optvaldir ,0777);
	  chdir( $optvaldir );
	  # Vary the rest of the options
	  my $v_r = $self->_vary_options( $abs_opt_j + 1, $nlp_runner, $v_set_option_value
									  , $best_run_stats, $runs_summaries
									  , $new_nlp, $out );
	  chdir( ".." );
	  if( !$v_r ) {
		return 0; # No more options to vary
	  }
	}
  }
  else {
	# This is the end of the recusion so just create the
	# options file and run the test.
	#
	# Create the varied options lines for the options file
	my @varied_options;
	my ( $i, $absj );
	for( $i = 0; $i < scalar(@$v_opt_grp); $i++ ) {
	  @varied_options = ( @varied_options
						  , "options_group " . $v_opt_grp->[$i] . " {" );
	  for( $absj = $v_first_opt_in_grp->[$i];
		   $absj < $v_first_opt_in_grp->[$i+1];
		   $absj++ 								)
		{
		  @varied_options = ( @varied_options
							  , "    " . $v_opt->[$absj] . " = "
							  . $v_opt_value->[$v_set_option_value->[$absj]]
							  . ";" );
		}
	  @varied_options = ( @varied_options, "}" );
	}
	# Create the options file
	$self->_create_options_file(
								"begin_options\n",
								"\n*******************",
								"*** Default options\n",
								@$default_opts,
								"\n*******************",
								"*** Varied options\n",
								@varied_options,
								"\nend_options"			);
	# Print summary of these options
	my $varied_options = "";
	foreach (@$v_set_option_value) {
	  $varied_options = $varied_options . $v_opt_value_abbr->[ $_ ] . ",";
	}
	print STDERR
	  "Running MOOCHHO for options : $varied_options ... ";
	# Run MOOCHO on this NLP with these options
	$nlp_runner->run_nlp();
	# Create the header for this NLP if this is first run
	if( $$new_nlp ) {
	  $self->_create_nlp_header($best_run_stats,$out);
	  $$new_nlp = 0;
	}
	# Read the output and fill in the summary line
	my $status
	  = $self->_add_summary_line($v_set_option_value,$best_run_stats,$runs_summaries,$out);
	print STDERR
	  " -> $status\n";
  }
  return 1; # There are more options to vary
}

#
# Create an options file in the current directory from the
# list of line given in the arguments.
#
sub _create_options_file {
  #
  my $self   = shift;
  #
  my $opt_file_name = $self->{OPT_FILE_NAME()};
  #
  my $opt_fh;
  ($opt_fh= FileHandle->new( ">$opt_file_name" ) )
	|| die "Error, could not open file \'$opt_file_name\': $!";
  foreach (@_) {
	print $opt_fh
	  $_, "\n";
  }
}

#
# Create initial header for th NLP.
#
sub _create_nlp_header {
  # Arguments
  my $self            = shift;
  my $best_run_stats  = shift; # Reference to a hash of arrays (see run_nlps(...))
  my $out             = shift; # Output file handle
  #
  my $opt_val_stats  = $self->{OPT_VAL_STATS()}; # ref to hash of arrays
  #
  print $out
	$g_summary_cols_header_labels,
	$g_summary_cols_lines_sep;
}

#
# Adds a summary line for a run of MOOCHO.
#
# This sub reads the file MoochoStats.out.
#
# Here the runs with the minimum time (min_t), minimum iterations (min_it)
# and minimum function evaluations (min_fn) are remembered for later
# output.
#
# This sub returns the solution status of "sol", "max_it", "max_t"
# or "except".
#
sub _add_summary_line {
  # Arguments
  my $self               = shift;
  my $v_set_option_value = shift; # reference to array for current options values
  my $best_run_stats     = shift; # Reference to a hash of arrays (see run_nlps(...))
  my $runs_summaries     = shift; # Refererence to array of arrays (see run_nlps(...))
  my $out                = shift; # Output file handle
  #
  my $v_opt_value_abbr   = $self->{VARY_OPTS()}->{V_OPT_VALUE_ABBR()};  # ref to array
  #
  my $status = "except";
  my $varied_options = "";
  foreach (@$v_set_option_value) {
	my $abbr = $v_opt_value_abbr->[$_];
	if( length($abbr) ) {
	  $varied_options = $varied_options . $abbr . ",";
	}
  }
  my $niter     = '-';
  my $nfunc     = '-';	
  my $ngrad     = '-';
  my $CPU       = 0;
  my $fx        = 0;
  my $nrm_cx    = -1;
  my $kkt_error = -1;
  my $nact      = '-';
  my $nBC       = '-';
  my $nQN       = '-';
  # Read MoochoStats.out and fill this stuff in
  my $stats_fh;
  ($stats_fh = FileHandle->new("<MoochoStats.out"));
# || die "Error, the file MoochoStats.out could not be opened: $!";
#  if (!$stats_fh)
#    {
#	  ($stats_fh = FileHandle->new("</home/cdlaird/research/rsqpdev/rSQPppApplications/CUTErSQPpp/rSQPppStats_FailedCompile.out"))
#		|| die "Error, MoochoStats.out & rSQPppStats_FailedCompile.out cannot be opened";
#	}
  while(<$stats_fh>) { # This may not read anything if an exception is thrown early!
	my ($stat, $val) = read_stats_val($_);
	if( $stat eq "status" ) {
	  if( $val eq "solved" ) {
		$status = "sol";
	  }
	  elsif( $val eq "except" ) {
		$status = "except";
	  }
	  elsif( $val eq "max_iter" ) {
		$status = "max_it";
	  }
	  elsif( $val eq "max_run_time" ) {
		$status = "max_t";
	  }
	  else {
		die "Error, wrong value of $stat = $val from MoochoStats.out file";
	  }
	}
	elsif( $stat eq "niter" ) {
	  $niter = $val unless $val eq "-";
	}
	elsif( $stat eq "nfunc" ) {
	  $nfunc = $val unless $val eq "-";
	}
	elsif( $stat eq "ngrad" ) {
	  $ngrad = $val unless $val eq "-";
	}
	elsif( $stat eq "CPU" ) {
	  $CPU = $val unless $val eq "-";
	}
	elsif( $stat eq "obj_func" ) {
	  $fx = $val unless $val eq "-";
	}
	elsif( $stat eq "feas_kkt_err" ) {
	  $nrm_cx = $val unless $val eq "-";
	}
	elsif( $stat eq "opt_kkt_err" ) {
	  $kkt_error = $val unless $val eq "-";
	}
	elsif( $stat eq "nact" ) {
	  $nact = $val unless $val eq "-";
	}
	elsif( $stat eq "nbasis_change" ) {
	  $nBC = $val unless $val eq "-";
	}
	elsif( $stat eq "nquasi_newton" ) {
	  $nQN = $val unless $val eq "-";
	}
	else {
	  die "Error, wrong value of $stat = $val from MoochoStats.out file";
	}
  }
  $stats_fh->close() if defined($stats_fh);	
  # Create the array to print the summary line
  my @sline;
  $sline[$g_status_i]			= $status;
  $sline[$g_varied_options_i]	= $varied_options;
  $sline[$g_niter_i]			= $niter;
  $sline[$g_nfunc_i]			= $nfunc;	
  $sline[$g_ngrad_i]			= $ngrad;
  $sline[$g_CPU_i]				= $CPU;
  $sline[$g_fx_i]				= $fx;
  $sline[$g_nrm_cx_i]			= $nrm_cx;
  $sline[$g_kkt_error_i]		= $kkt_error;
  $sline[$g_nact_i]				= $nact;
  $sline[$g_nBC_i]				= $nBC;
  $sline[$g_nQN_i]				= $nQN;
  print_summary_line( $out, \@sline );
  # Remember this summary line for later
  push @{$runs_summaries->[SUM_STATUS()]},        $status;
  push @{$runs_summaries->[SUM_NITER()]},         $niter;
  push @{$runs_summaries->[SUM_NFUNC()]},         $nfunc;
  push @{$runs_summaries->[SUM_CPU()]},           $CPU;
  push @{$runs_summaries->[SUM_OPT_VALUES()]},    [ @$v_set_option_value ];
  # Keep track of who is the best
  my $line_min_it  = $best_run_stats->{STATS_MIN_IT()}; # Ref to array
  my $line_min_fn  = $best_run_stats->{STATS_MIN_FN()}; # Ref to array
  my $line_min_t   = $best_run_stats->{STATS_MIN_T()};  # Ref to array
  #
  if( $status eq "sol" ) {
	if( $niter < $line_min_it->[$g_niter_i] ) {
	  @$line_min_it = @sline;
	  $line_min_it->[$g_status_i] = STATS_MIN_IT();
	}
	if( $nfunc < $line_min_fn->[$g_nfunc_i] ) {
	  @$line_min_fn = @sline;
	  $line_min_fn->[$g_status_i] = STATS_MIN_FN();
	}
	if( $CPU < $line_min_t->[$g_CPU_i] ) {
	  @$line_min_t = @sline;
	  $line_min_t->[$g_status_i] = STATS_MIN_T();
	}
  }
 
  return $status
}

#
# After all the options for an NLP have been run this sub
# outputs the best runs.
#
# This sub reads the file MoochoSummary.out.
# Here the runs with the minimum time (min_t), minimum iterations (min_it)
# and minimum function evaluations (min_fn) are output.
#
sub _add_best_runs_lines {
  # Arguments
  my $self               = shift;
  my $v_set_option_value = shift; # reference to array for current options values
  my $best_run_stats     = shift; # Reference to a hash of arrays (see run_nlps(...))
  my $runs_summaries     = shift; # Refererence to array of arrays (see run_nlps(...))
  my $out                = shift; # Output file handle
  #
  my $line_min_it  = $best_run_stats->{STATS_MIN_IT()}; # Ref to array
  my $line_min_fn  = $best_run_stats->{STATS_MIN_FN()}; # Ref to array
  my $line_min_t   = $best_run_stats->{STATS_MIN_T()};  # Ref to array
  #
  my $any_solved = ( $line_min_it->[$g_niter_i] != $g_real_big );
  my ( $min_it, $min_fn, $min_t );
  if( $any_solved ) {
	$min_it =  $line_min_it->[$g_niter_i];
	$min_fn =  $line_min_fn->[$g_nfunc_i];
	$min_t  =  $line_min_t->[$g_CPU_i];
  }
  # Record the statistics for each option value
  my $run_status     = $runs_summaries->[SUM_STATUS()];     # ref to array
  my $run_niter      = $runs_summaries->[SUM_NITER()];      # ref to array
  my $run_nfunc      = $runs_summaries->[SUM_NFUNC()];      # ref to array
  my $run_CPU        = $runs_summaries->[SUM_CPU()];        # ref to array
  my $run_opt_values = $runs_summaries->[SUM_OPT_VALUES()]; # ref to array of string arrays
  my $stats_min_it   = $self->{OPT_VAL_STATS()}->{STATS_MIN_IT()}; # ref to array
  my $stats_min_fn   = $self->{OPT_VAL_STATS()}->{STATS_MIN_FN()}; # ref to array
  my $stats_min_t    = $self->{OPT_VAL_STATS()}->{STATS_MIN_T()};  # ref to array
  my $stats_max_it   = $self->{OPT_VAL_STATS()}->{STATS_MAX_IT()}; # ref to array
  my $stats_except   = $self->{OPT_VAL_STATS()}->{STATS_EXCEPT()}; # ref to array
  for( my $k = 0; $k < @$run_status; ++$k ) {
	if( $run_status->[$k] eq "sol" ) {
	  my $num_it_ratio = ($run_niter->[$k] - $min_it) / ($min_it + $g_real_small);
	  $self->_add_sol_stats( $num_it_ratio, $run_opt_values->[$k], $stats_min_it );
	  my $num_fn_ratio = ($run_nfunc->[$k] - $min_fn) / ($min_fn + $g_real_small);
	  $self->_add_sol_stats( $num_fn_ratio, $run_opt_values->[$k], $stats_min_fn );
	  my $num_t_ratio = ($run_CPU->[$k] - $min_t) / ($min_t + $g_real_small);
	  $self->_add_sol_stats( $num_t_ratio, $run_opt_values->[$k], $stats_min_t );
	  $self->{TOTAL_NLPS_SOLVED()}++
	}
	elsif( $run_status->[$k] eq "max_it" ) {
	  my $abs_opt_val_k;
	  foreach $abs_opt_val_k (@{$run_opt_values->[$k]}) {
		$stats_max_it->[$abs_opt_val_k]++;
	  }
	}
	else { # assume this is an exception.
	  my $abs_opt_val_k;
	  foreach $abs_opt_val_k (@{$run_opt_values->[$k]}) {
		$stats_except->[$abs_opt_val_k]++;
		$self->{TOTAL_NLPS_EXCEPT()}++
	  }
	}
  }
  # Increment the total number of runs
  $self->{TOTAL_NLP_RUNS()} += scalar(@$run_status);
  # Print the "best of" summary lines
  if( !$any_solved ) {
	return;
  }
  print $out
	$g_summary_cols_lines_sep;
  print_summary_line( $out, $line_min_fn );
  print_summary_line( $out, $line_min_it );
  print_summary_line( $out, $line_min_t );
}

#
# Print the summary statistics to all the option values for
# the given "best of" measure.
#
sub _add_sol_stats {
  # Arguments
  my $self            = shift;
  my $percent_of_best = shift; # scalar
  my $opt_values      = shift; # ref to array
  my $stat_opt_array  = shift; # ref to array
  # Convert the precentage of the best to an integer to index into stats.
  my $sol_i = $percent_of_best * 10;
  {
	use integer;
	$sol_i  = 0 + $sol_i; # This should round down?
  }
  if( $sol_i < $percent_of_best * 10 ) {
	$sol_i++;
  }
  if( $sol_i > 11 ) {
	$sol_i = 11; # It is more than 100% worst than the winner.
  }
  # Loop through and update the statistics.
  if( $sol_i == 11 ) {
	my $opt_val_abs_k;
	foreach $opt_val_abs_k (@{$opt_values}) {
	  $stat_opt_array->[$opt_val_abs_k]->[11]++;
	}
  }
  else {	# $sol_i < 11
	my $opt_val_abs_k;
	foreach $opt_val_abs_k (@{$opt_values}) {
	  my $sol_stats = $stat_opt_array->[$opt_val_abs_k];
	  for( my $i = $sol_i; $i <= 10; $i++) {
		$sol_stats->[$i]++;
	  }
	}
  }
}

#
# Print stats for an option for max_it and except.
#
sub _print_opt_values_stats_max_it_except {
  # Arguments
  my $self      = shift;
  my $opt_abs_j	= shift;
  my $out       = shift;
  #
  my $v_first_opt_value  = $self->{VARY_OPTS()}->{V_FIRST_OPT_VALUE()}; # ref to array
  my $v_opt_value_abbr   = $self->{VARY_OPTS()}->{V_OPT_VALUE_ABBR()};  # ref to array
  #
  my $w1 = 12;
  my $w2 = 10;
  # Print Header
  print $out
	"\nnot solved:\n";
  printf $out
	"%${w1}s", "";
  my $opt_val_abs_k;
  for(	$opt_val_abs_k = $v_first_opt_value->[$opt_abs_j];
		$opt_val_abs_k < $v_first_opt_value->[$opt_abs_j+1];
		$opt_val_abs_k++ )
	{
	  printf $out
		"%${w2}s", $v_opt_value_abbr->[$opt_val_abs_k];
	}
  print $out
	"\n";
  printf $out
	"%${w1}s", "-"x($w1-2);
  for(	$opt_val_abs_k = $v_first_opt_value->[$opt_abs_j];
		$opt_val_abs_k < $v_first_opt_value->[$opt_abs_j+1];
		$opt_val_abs_k++ )
	{
	  printf $out
		"%${w2}s", "-"x($w2-2);
	}
  print $out
	"\n";
  # Print the lines
  $self->_print_opt_values_stats_max_it_failure_line($opt_abs_j,STATS_MAX_IT(),STATS_MAX_IT(),$out);
  $self->_print_opt_values_stats_max_it_failure_line($opt_abs_j,STATS_EXCEPT(),STATS_EXCEPT(),$out);
}

#
# Print a line for the statistics for all the varied values
# for an option for max iter or failed varaibles
#
sub _print_opt_values_stats_max_it_failure_line {
  # Arguments
  my $self      = shift;
  my $opt_abs_j	= shift;
  my $field		= shift;
  my $label		= shift;
  my $out       = shift;
  #
  my $v_first_opt_value  = $self->{VARY_OPTS()}->{V_FIRST_OPT_VALUE()}; # ref to array
  my $stat_stat          = $self->{OPT_VAL_STATS()}->{$field};          # ref to array
  #
  my $w1 = 12;
  my $w2 = 10;
  printf $out
	"%${w1}s", substr($label,0,$w1-2);
  my $opt_val_abs_k;
  for(	$opt_val_abs_k = $v_first_opt_value->[$opt_abs_j];
		$opt_val_abs_k < $v_first_opt_value->[$opt_abs_j+1];
		$opt_val_abs_k++ )
	{
	  my $val = $stat_stat->[$opt_val_abs_k]; # This is a scalar for max_it and execpt!
	  printf $out
		"%${w2}s", $val;
	}
  print $out
	"\n";
}

#
# Print stats for an option and a "best of" measure.
#
sub _print_opt_values_stats {
  # Arguments
  my $self        = shift;
  my $opt_abs_j   = shift;
  my $best_at_cat = shift;
  my $out         = shift;
  #
  my $v_first_opt_value  = $self->{VARY_OPTS()}->{V_FIRST_OPT_VALUE()}; # ref to array
  my $v_opt_value_abbr   = $self->{VARY_OPTS()}->{V_OPT_VALUE_ABBR()};  # ref to array
  #
  my $w1 = 12;
  my $w2 = 10;
  # Print Header
  print $out
	"\n", $best_at_cat, ":\n";
  printf $out
	"%${w1}s", "";
  my $opt_val_abs_k;
  for(	$opt_val_abs_k = $v_first_opt_value->[$opt_abs_j];
		$opt_val_abs_k < $v_first_opt_value->[$opt_abs_j+1];
		$opt_val_abs_k++ )
	{
	  printf $out
		"%${w2}s", $v_opt_value_abbr->[$opt_val_abs_k];
	}
  print $out
	"\n";
  
  printf $out
	"%${w1}s", "-"x($w1-2);
  for(	$opt_val_abs_k = $v_first_opt_value->[$opt_abs_j];
		$opt_val_abs_k < $v_first_opt_value->[$opt_abs_j+1];
		$opt_val_abs_k++ )
	{
	  printf $out
		"%${w2}s", "-"x($w2-2);
	}
  print $out
	"\n";
  # Print the lines
  $self->_print_opt_values_stats_sol_line($opt_abs_j,$best_at_cat,0,"winners",$out);
  # <= 10% ... <= 100%
	my $sol_i;
  for( $sol_i = 1; $sol_i <= 10; $sol_i++ ) {
	$self->_print_opt_values_stats_sol_line($opt_abs_j,$best_at_cat,$sol_i
											,"<= ".$sol_i."0%",$out);
  }
  printf $out
	"%${w1}s", "-"x($w1-2);
  for(	$opt_val_abs_k = $v_first_opt_value->[$opt_abs_j];
		$opt_val_abs_k < $v_first_opt_value->[$opt_abs_j+1];
		$opt_val_abs_k++ )
	{
	  printf $out
		"%${w2}s", "-"x($w2-2);
	}
	print $out
	  "\n";
  $self->_print_opt_values_stats_sol_line($opt_abs_j,$best_at_cat,11,"> 100%",$out);
}

#
# Print a line for the statistics for all the varied values
# for an option
#
sub _print_opt_values_stats_sol_line {
  # Arguments
  my $self        = shift;
  my $opt_abs_j	  = shift;
  my $best_at_cat = shift;
  my $sol_i		  = shift;
  my $label		  = shift;
  my $out         = shift;
  #
  my $v_first_opt_value  = $self->{VARY_OPTS()}->{V_FIRST_OPT_VALUE()}; # ref to array
  my $stat_stat          = $self->{OPT_VAL_STATS()}->{$best_at_cat};    # ref to array
  #
  my $w1 = 12;
  my $w2 = 10;
  printf $out
	"%${w1}s", substr($label,0,$w1-2);
  my $opt_val_abs_k;
  for(	$opt_val_abs_k = $v_first_opt_value->[$opt_abs_j];
		$opt_val_abs_k < $v_first_opt_value->[$opt_abs_j+1];
		$opt_val_abs_k++ )
	{
	  my $val = $stat_stat->[$opt_val_abs_k]->[$sol_i];
	  printf $out
		"%${w2}s", $val ? $val : 0;
	}
  print $out
	"\n";
}

#######################
# Utility subroutines #
#######################

#
# Remove leading and trailing whitespace
#
sub r_ws {
  # Arguments
  my $str = shift;
  #
  $str=~s/^\s*//g;
  $str=~s/\s*$//g;
  #
  return $str;
}

#
# Print out a summary line from an input array.
#
sub print_summary_line {
  # Arguments
  my $out  = shift; # Filehandle
  my $line = shift; # Ref to array of strings
  #
  printf $out
	$g_summary_cols_lines_format,
	@$line; #
}

#
# Read an option from the MoochoStats.out file
#
sub read_stats_val {
  #
  my $str = shift;
  #statistic = value; # Comments
  my $equal_poss      = index($str,"=");
  my $semi_colon_poss = index($str,";");
  return ( r_ws(substr($str,0,$equal_poss))
		   ,r_ws(substr($str,$equal_poss+1,$semi_colon_poss-$equal_poss-1)) );
}

# End with true
1
