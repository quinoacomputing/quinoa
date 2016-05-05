[![License](https://img.shields.io/github/license/quinoacomputing/ndiff.svg)](https://github.com/quinoacomputing/ndiff/blob/master/LICENSE)

## ndiff

_ndiff_ is an efficient and flexible tool designed to compare unformatted text files with numerical content. It is well suited for regression testing, for validating data versus templates, and for filtering data following templates.

## Author
Laurent Deniau (laurent.deniau@cern.ch) CERN - BE/ABP

## Motivation
The _ndiff_ tool is a program developed for the test system of the [MAD-X application](http://cern.ch/mad). The purpose of this tool is to compare line by line unformatted text files with numerical content in a portable way. The portability of the comparison is resulting from the acceptance of (user-defined) small numerical differences and varying number representation that occur across runs, compilers and platforms. For example, _ndiff_ will consider to be equal the numbers 0.001, 1e-3, 100.00e-5 and 0.0001e+001, while _diff_-like Unix tools will report differences.

To the knowledge of the author, the _ndiff_ tool does not seem to have any equivalent freely available on the world wide web, despite the obvious interest and need for the scientific community. Other similar open source tools exist like [numdiff](http://www.nongnu.org/numdiff) or [ndiff](http://www.math.utah.edu/~beebe/software/ndiff) (same name but not same tool), but they do not provide the flexibility required by the test system of _MAD-X_. Many options and commands have been added to _ndiff_ in order to solve problems related to testing _MAD-X_, hence examples will often refer to the _MAD-X_ test system.

The _ndiff_ tool is written entirely in the C programming language with no external dependencies to ensure maximum portability and performance. It is very efficient to deal with both large input of data and large number of rules (i.e., user-defined constraints). _ndiff_ can process more than 100 MB of data per second (i.e., nearly at the rate of disk I/O) and compare few millions numbers per second on recent computers.

## Documentation
https://github.com/quinoacomputing/ndiff/blob/master/doc/CERN-ACC-NOTE-2013-0005.pdf

## Fork
This forks http://svn.cern.ch/guest/madx/trunk/madX/tools/numdiff.
