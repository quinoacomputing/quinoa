#!/usr/bin/perl

require "shellwords.pl";

while(<>) {
  chop;
  if ($_ ne "")
    {
      ($azkey, $azval) = &shellwords($_);
      print "paramMap()[\"",$azval,"\"]=", $azkey, "; \n";
    }
}
