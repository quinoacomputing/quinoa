#!/usr/bin/perl 

# 
# Wrapper for compilers having troubles with preprocessor.
#
# Assumes the compiler is in your PATH
# 
# Written for a version of NAG that could not preprocess
# lines longer than 132 chars.
# 
$compiler=$ENV{TRILINOS_NAG_COMPILER} || "nagfor";
$forcepre="-fpp";
#
# Scan through the argument list separating:
#   1. input files needing preprocessing:
#      .Fnn .ffnn  where nn are two digits (covers
#              f77 f90 f95 f03)
#      .F   .ff   
#   2. input files not (necessarily) needing preprocessing:
#      .fnn .f   SHOULD WE ADD .FOR?? 
#
#   3. Options affecting the preprocessor:
#      -Idir -Dmacro
#      
#   4. Output clause: -o outfile or -ooutfile
#   
#   5. Option to force preprocessing even for files in 2.
#
#   6. -c to have compile but not link (affects output file name)
#   
#   7. Any other option, irrelevant to the preprocessor
#      and intended for the base fortran compiler
#
#   8. Script-specific options: currently only -debug
#      They must be removed from @ARGV, lest they get
#      passed to the base compiler. 
#

MAIN: for ($i=0; (($i<=$#ARGV), $_=$ARGV[$i]) ; $i++) {
#    print "debug: $_ \n";
    SWITCH: {
	if (/.*\.F[0-9][0-9]$|.*\.ff[0-9][0-9]$|.*\.F$|.*\.ff$/) {
	    push (@infiles,$_);
	    $preproc=1;
	    last SWITCH;	    
	}
	
	if (/.*\.f[0-9][0-9]$|.*\.f$/) {
	    push (@infiles,$_);
	    last SWITCH;	   
	}
	if (/^-I|^-D/) {
	    push (@ccopt,$_);
	    last SWITCH;   
	}
	
	if (/^-debug$/) { 
	    $debug=1;
	    splice(@ARGV,$i,1); 
	    $i--;
	    last SWITCH;   
	}
		
	if (/^-o$/) {
	    $outfile=$ARGV[++$i];
	    last SWITCH;   
	}
	
	if (/^-o(.*)$/) {
	    $outfile=$1;
	    last SWITCH;   
	}
	
	if (/^$forcepre$/) {
	    $preproc=1;
	    last SWITCH;   
	}
	if (/^-c$/) {
	    $componly=1;
	    push (@fopt,$_);
	    last SWITCH;   
	}
	
        push (@fopt,$_);
	
    }
}

#
# Ok, so: we have to preprocess if 
#  1. the flag is set and
#  2. there is an input file
#   remember we might just have been called to link
#   an existing .o
#

if ((@infiles)&&($preproc)) {
    #
    #
    # Build a name for the temporary file
    # holding preprocessed source. Should be the same basename
    # with single lowercase "f" plus digits, if any.
    # Take out prefix to write in current directory, just
    # in case source dir is not writable, and make sure not
    # to conflict with existing files.
    # 
    #
    if ($debug) { print "Preproc branch @infiles\n";}
    foreach $infile (@infiles) {
	$pfile=$infile;
	$pfile =~ s/.*\///;
	$pfile =~ s/\.[Ff][Ff]*([0-9]*)$/\.f\1/;
	while ( -e $pfile ) {
	    $pfile="XtmpX_$pfile";
	}
	
	#
	# Do the preprocessing; adjust return code if needed. 
	#
	
	push (@cpp,"gcc","-undef","-x", "c","-E","-P");

	if ($pfile =~ /\.[Ff]$/) {
	    # This is really needed for .f files
	    push (@cpp,"-traditional-cpp");
	} else {
	    push (@cpp,"-ansi");
	}
	push (@cpp,@ccopt,$infile,"-o",$pfile);
	if ($debug) { print "@cpp \n";}
	$rc=system(@cpp);
	
	if ($rc != 0) { 
	    foreach $tfile (@pfiles) {
		unlink $tfile;
	    }
	    exit $rc>>8;
	}
	@cpp=();
	push (@pfiles,$pfile);
    }
    
    if ($debug) { print "Compiling: @pfiles\n";}
    #
    # Now for compilation: if outfile was not forced, and if
    # we are supposed to produce a .o file, we have to
    # build the correct name.
    # Need to push the -I options, if any; they are often used
    # to specify search path for module files. 
    #
    #    push (@fcomp,$compiler,@ccopt,@fopt,$pfile);
    if ((!$outfile)&&($componly)) {
	#
	# Treat @pfiles one at a time
	#
	if ($debug) { print "Trating infiles: @infiles $#infiles @pfiles $#pfiles\n";}
      CLOOP: for ($i=0; $i<=$#infiles; $i++) {
	  $infile=$infiles[$i];
	  $pfile=$pfiles[$i];
	  $outfile=$infile;
  	  $outfile =~ s/.*\///;
	  $outfile =~ s/\.[Ff][Ff]*([0-9]*)$/\.o/;
	  push (@fcomp,$compiler,@ccopt,@fopt,$pfile,"-o",$outfile);
	  if ($debug) { print "Invoking: @fcomp\n";}
	  $rc=system(@fcomp);
	  if ($rc != 0) { 
	      last CLOOP; 
	  }
	  @fcomp=();
      }
    } else {
	push (@fcomp,$compiler,@ccopt,@fopt,@pfiles);
	if ($outfile) {
	    push (@fcomp,"-o",$outfile);
	}
	if ($debug) { print "Invoking: @fcomp\n";}
	$rc=system(@fcomp);
    } 
    if ($debug) { print "cleaning up\n";}
    foreach $tfile (@pfiles) {
	unlink $tfile;
    }
    exit $rc>>8;
    
} else {
    #
    # No input file was recognized, and/or no preprocessing is necessary.
    # Perhaps  we were just liking?
    # In any case, invoke the base compiler
    # with the arglist (minus script specific options if any)
    # 
    push (@args,$compiler,@ARGV);
    if ($debug) { print "Execing @args \n";}
    exec @args;
}
