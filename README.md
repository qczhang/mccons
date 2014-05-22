MC-Cons II: Attack of the look-alike
====================================

## Description
This is a multiobjective similarity search algorithm
based on the nondominated sorting genetic algorithm (NSGA-II).

## Quick start
Examples are available in the *examples* folder.

## Documentation
Documentation is available in the *doc* folder.

## Use case
### input
- *n* sets of RNA sequences
- candidate structures for each of them
- *m* distance functions (at least 2, but too much will make results crappy)

### output
- *s* sets of n structures (*s*>=1) that are non-dominated according to distance functions.


## Limitations
- distance functions must be symmetric
- output is not deterministic
- all the limitations of NSGA-II


## Dependency
- [julia 0.3](https://github.com/JuliaLang/julia)


## Authors
Gabriel C-Parent did most of the software and design.

Stefanie Schirmer extracted test data and helped with validation and documentation.


## License
All software is distributed under MIT licence.

## Original MC-Cons
Succint description in this [supplementary material](http://www.nature.com/nature/journal/v452/n7183/extref/nature06684-s1.pdf).
