#!/usr/bin/env python

import commands
import os
import re
import sys

if len(sys.argv) != 4:
  raise SystemExit("Syntax:\n  %s oldrev newrev refname" % sys.argv[0])

oldrev=sys.argv[1]
newrev=sys.argv[2]
refname=sys.argv[3]
index=refname.rfind("/")
shortrefname=refname[index+1:].lower()

output=commands.getoutput("git diff --name-only %s %s" % (oldrev, newrev))
dirschanged = [os.path.dirname(filename) for filename in output.splitlines()]
defaultEmail = commands.getoutput("git config --get hooks.mailinglist").strip()
#print "defaultEmail =", defaultEmail

dirs = {}.fromkeys(dirschanged, 1)
#print "dirs =", dirs

# Create a list of (regex, email) in the order they are in the file.  That
# way, the first regex that matches will be used to select the email address.

#print "sys.argv =", sys.argv
dirs_to_email_file = os.path.join(os.path.dirname(sys.argv[0]),'dirs_to_emails')
#print "dirs_to_email_file =", dirs_to_email_file
f = open(dirs_to_email_file, 'r')

emails = []

for raw_line in f:
  line = raw_line.strip()
  #print "\nline = '"+line+"'"
  if line.startswith('#') or line == "":
    continue
  (regex, email) = line.split()
  emails.append((regex,email))

f.close()

#for regex_email in emails:
#  print "%s: %s" % regex_email

recipients = {}

found = False

#if this is a branch for a package then only send the email to that package.
#otherwise send it to the list of packages that were modified.
#this is to cut down on the number of "extra" emails that get sent out when
#someone merges master onto a package branch and then pushes.
if not re.match("^master", shortrefname) and not re.match("^trilinos-release", shortrefname):
  for pair in emails:
    index = pair[1].find("-")
    packageRE = pair[1][:index]
    if re.match(packageRE, shortrefname):
      recipients[pair[1]] = 1
      found = True
      break
else:
  for dir in dirs:
    #print "\ndir =", dir
    for regex_email in emails:
      (regex, email) = regex_email
      #print "\nregex =", regex
      #print "email =", email
      if re.match(regex, dir):
        recipients[email] = 1
        found = True
        #print "FOUND IT!"
        break

if not found:
  recipients[defaultEmail] = 1
  #print "\nNOT FOUND: Use default email!"

print ",".join(sorted(recipients.keys()))
