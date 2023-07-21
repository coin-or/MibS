MibS Change Log

## 1.2.0

 * Generalizing the Benders binary cut based on analyzing signs of matrix A2
 * Generalizing Benders interdiction cut based on analyzing signs of
matrix G2
 * Generalizing ICs to allow separation of fractional solutions
 * Substantial improvements to performance with default parameter settings
 * Remove need to specify application-specific prefix for parameters
 * Renaming of cuts to match new draft of this paper that details
performance of MibS 1.2.
 * Lots cleanup of the repo, fixes for usability, refining of output
refactoring.
	
## 1.1.3

 * Adding parametric bound cut
 * Improving numerics
 * Debugging cuts
 * Improvements to user interface
	
## 1.1.2

 * Fix the fix
 * Version used for version of this
[paper](http://coral.ie.lehigh.edu/~ted/files/papers/MIBLP16.pdf)
published in MPC.
	'
## 1.1.1

 * Fix bug in configure script related to line endings

## 1.1.0

 * Improve automatic instance structure detection
 * Improve automatic parameter setting based on instance structure
 * Generalize increasing objective cut
 * Add watermelon ICs
 * Add new type of Benders cut
 * General improvements to cuts
 * Debugging and improving heuristics
 * Efficiency improvements
 * Many other miscellaneous improvements
 * Bug fixes
	
	
## 1.0.0

 * Added new name-based file format
 * Fixed small bugs and memory leaks
 * Improve robustness	
 * Improved default parameter settings	
	
## 0.95.1

 * Updates for new experiments in updated draft.
 * Changed parameter names to match paper
 * Fixed bug in linking solution pool

## 0.95.0

 * Version used for initial draft of
[paper](http://coral.ie.lehigh.edu/~ted/files/papers/MIBLP16.pdf)
 * Added linking solution pool
 * Refactoring	

## 0.9.0

 * First public release
 * Separated MibS out into a library with a rudimentary API
 * Added autotools build harness
