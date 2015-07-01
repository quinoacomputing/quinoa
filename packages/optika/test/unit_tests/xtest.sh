#!/bin/bash
if [[ ! $(xset -q) ]] ; then
  echo "Couldn't reach XServer"
  exit -1 
else
  echo "XServer reached"
  exit 0
fi
