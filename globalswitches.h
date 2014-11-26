/*
 * globalswitches.h
 *
 * This file contains global #defines used to control behaviour such as
 * debugging, numerical, etc.
 *
 *  Created on: 19 Nov 2011
 *      Author: dan
 */

#ifndef GLOBALSWITCHES_H_
#define GLOBALSWITCHES_H_

////////////////////////////////////////////////////////////////////////////
// Allow for working in single or double precision
////////////////////////////////////////////////////////////////////////////
#define USE_DOUBLE_PRECISION

////////////////////////////////////////////////////////////////////////////
// Whether or not to use reduced amounts of RAM, depending on the data this
// may reduce the precision and either speed up or slow down the computation
////////////////////////////////////////////////////////////////////////////
#define STORE_FLOAT_MASK_AS_CHAR // reduces precision
#define COMPUTED_BINARY_MASK
#define STORE_DATA_AS_SHORT // reduces precision in time variable // TODO: not yet implemented
#define COMPUTED_CORRECTION_TERM // reduces precision in principle, but seems not to matter, can slow down computation?

// Format for printing data
#define SCALARFMT "%f"

// Switch this definition off to remove bounds checking in array accesses,
// when constructing a SafeArray, bounds checking is still performed when
// constructing from a vector.
//#define SAFEARRAY_BOUNDSCHECKING

///////////////// CACHE OPTIMISATION PARAMETERS ////////////////////////////
// TODO: this value probably not optimal for larger numbers of dimensions, need
// to profile on a larger example. Somewhere around 90 might be optimal for
// larger data sets, but not for the one I'm testing on (nDims=97).
#define COVARIANCE_BLOCKSIZE 128

#endif /* GLOBALSWITCHES_H_ */
