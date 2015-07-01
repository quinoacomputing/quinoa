#!/usr/bin/env python


import os


# -----------------------------------------------------------------------------------
# Run this from the top-level Sundance source directory.
# -----------------------------------------------------------------------------------

optDir = 'LINUX-PY'
debugDir = 'LINUX-PY-DEBUG'

pyLibDir = 'python/src/build/lib.linux-i686-2.3'
versionString = '2.1.2-RH-07032006'

rebuild = 0

if rebuild:
    print 'building optimized libraries'
    os.system('(cd %s; make)' % optDir)
    print 'building debug libraries'
    os.system('(cd %s; make)' % debugDir)


os.system('cp etc/SundanceConfig.xml %s/%s/PySundance' % (optDir, pyLibDir))
os.system('cp etc/SundanceConfig.xml %s/%s/PySundance' % (debugDir, pyLibDir))

optTarball = 'PySundance-opt-%s.tgz' % versionString
debugTarball = 'PySundance-debug-%s.tgz' % versionString

print 'building optimized tarball...'
os.system('tar zcf %s --directory %s/%s PySundance'
          % (optTarball, optDir, pyLibDir))

print 'building debug tarball...'
os.system('tar zcf %s --directory %s/%s PySundance'
          % (debugTarball, debugDir, pyLibDir))


user = 'krlong'
host = 'software.sandia.gov'
hostdir = '/var/www/html/sundance'

print 'uploading opt tarball...'
os.system('scp %s %s@%s:%s' % (optTarball, user, host, hostdir))

print 'uploading debug tarball...'
os.system('scp %s %s@%s:%s' % (debugTarball, user, host, hostdir))  

