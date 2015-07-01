#!/usr/bin/perl

use Term::ANSIColor;


# ---------------------------------------------------------------------------------------
# Usage: checkTests.pl [list of files]
#
# checkTests.pl reads test result files and classifies each as PASS, FAIL, or CRASHED.
# ---------------------------------------------------------------------------------------

# list of failed tests
@failures = ();
# list of crashed tests
@crashes = ();
# list of missing tests
@missing = ();
# list of inactive tests
@inactive = ();

$terminal = $ENV{"TERM"} || "dumb";

$nocolor{"dumb"} = "true";

# If it's in the list of terminal types not to color, or if
# we're writing to something that's not a tty, don't do color.
if (! -t STDOUT || $nocolor{$terminal})
{
  $skipColor = "true";
}

# Input arguments ARGV are filenames. Loop over filenames, scanning each for status.
foreach $i (0 .. $#ARGV)
  {
    # get filename
    $file = $ARGV[$i];
    if (-f $file)
      {
        # scan file for status. This will return PASSED, FAILED, or CRASHED
        $state = &checkFile($file);
      }
    else
      {
        $state = "MISSING";
      }
    $msgcolor = color("bold blue");
    # if FAILED, add filename to failure list
    if ($state =~ /FAILED/) {
      push(@failures, $file);
      $msgcolor = color("bold yellow");
    }
    # if MISSING, add filename to missing list
    if ($state =~ /MISSING/) {
      push(@missing, $file);
      $msgcolor = color("green");
    }
    # if CRASHED, add filename to crash list
    if ($state =~ /CRASHED/) {
      push(@crashes, $file);
      $msgcolor = color("bold red");
    }
    # if INACTIVE, add filename to crash list
    if ($state =~ /INACTIVE/) {
      push(@inactive, $file);
      $msgcolor = color("green");
    }
    # print status of this file
    printf "%-60s", ${file};
if (! $skipColor) {print $msgcolor;}
printf "%s\n", ${state};
if (! $skipColor) {print color("reset");}
}

# print messages if crashes or failures have occurred

$numFailures = @failures;
$numCrashes = @crashes;

if ($numFailures != 0)
  {
    print "-------------------------------------------------------------------------\n";
    print "\n";
    if (! $skipColor) {print color("bold yellow");}
    print "                   FAILURES DETECTED!!!!                                 \n";
    if (! $skipColor) {print color("reset");}
    print "\n";
    print "The following tests failed: ";
    if (! $skipColor) {print color("cyan");}
    print "@{failures}";
    if (! $skipColor) {print color("reset");}
    print "\n";
    print "-------------------------------------------------------------------------\n";
  }
if ($numCrashes != 0)
  {
    print "-------------------------------------------------------------------------\n";
    print "\n";
    if (! $skipColor) {print color("bold red");}
    print "                   CRASHES DETECTED!!!!                                 \n";
    if (! $skipColor) {print color("reset");}
    print "\n";

    print "The following tests crashed: ";
    if (! $skipColor) {print color("cyan");}
    print "@{crashes}";
    if (! $skipColor) {print color("reset");}
    print "\n";
    print "-------------------------------------------------------------------------\n";
  }

# otherwise, print a happy message

if ($numCrashes == 0 && $numFailures == 0)
  {
    print "-------------------------------------------------------------------------\n";
    print "\n";
    if (! $skipColor) {print color("bold blue");}
    print "                   All tests PASSED                                      \n";
    if (! $skipColor) {print color("reset");}
    print "\n";
    print "-------------------------------------------------------------------------\n";
  }

# ------- end of script ------
# Function definitions follow
# ----------------------------


# ---------------------------------------------------------------------------------------
# Function: Scan a file for PASS or FAIL. If neither is found, the state is CRASHED.
# ---------------------------------------------------------------------------------------
sub checkFile
  {
    open(FILE, pop(@_));

    $passed = 0;
    $failed = 0;
    $crashed = 0;
    $inactive = 0;
    $state = "CRASHED";
    while ($line = <FILE>) {
      $passed = $passed || $line =~ /PASS/;
      $failed = $failed || $line =~ /FAIL/;
      $inactive = $inactive || $line =~ /INACTIVE/;
    }
    
    if ($failed!=0)
      {
        $state = "FAILED";
      }
    else 
      {
        if($passed!=0)
        {
        	$state = "PASSED";
        }
        if($inactive!=0)
          {
            $state = "INACTIVE";
          }
      }

    return $state;
  }
