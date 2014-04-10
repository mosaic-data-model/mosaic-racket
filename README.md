About MOSAIC
------------

MOSAIC is a data model for molecular simulations. It defines how
various quantities of importance in molecular simulations can be
represented in terms of simple data structures (lists, trees,
etc.). MOSAIC also defines concrete file formats for storing molecular
simulation data. There are currently such file formats, one based on
XML and the other one based on HDF5. See the
[MOSAIC Web site](http://mosaic-data-model.github.io/) and in
particular the
[MOSAIC specification](http://mosaic-data-model.github.io/mosaic-specification/)
for more information about MOSAIC.

About the MOSAIC Racket library
-------------------------------

The MOSAIC Racket library implements an in-memory representation of
MOSAIC data items in the [Racket language](http://racket-lang.org/),
and provides I/O in the MOSAIC XML format. The HDF5 format is
currently not supported, because there is no HDF5 interface for Racket
at this time.

I wrote the MOSAIC Racket library for two reasons. First, I wanted to
learn Racket. Second, I wanted to experiment with MOSAIC in the
context of functional programming with immutable values. This is my
first non-trivial Racket code, so don't expect it to be perfect.  In
particular, it is probably overengineered and uses more Racket
features than it would have to. Remember that one of my goals was to
learn Racket.

Since this is experimental code, I do not make any promise at this
time to support it. In particular, I do not expect future versions
of this library to be 100% compatible with the current one.

Installation
------------

The package source for the Mosaic library is

    git://github.com/mosaic-data-model/mosaic-racket?path=mosaic

Use this with DrRacket's package manager, or with the command-line
tool `raco`. Watch out for the question mark in the package source. In
a command-line shell, you have to put a backslash before it to prevent
the shell from interpreting it inappropriately.

