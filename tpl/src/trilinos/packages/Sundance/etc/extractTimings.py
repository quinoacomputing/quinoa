#!/usr/bin/python

import sys, string

lines = sys.stdin.readlines()

for line in lines:
    words = string.split(line)
    if len(words) < 3:
        continue
    if words[0]=='assembly':
        print words[2]

    
