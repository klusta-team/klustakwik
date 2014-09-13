/*
 * numerics.h
 *
 * This header file is used for controlling the core numerical aspects of the
 * program, such as the basic numerical data type (float/double) and the
 * SafeArray class for making array accesses a bit safer while remaining
 * efficient.
 *
 *  Created on: 11 Nov 2011
 *      Author: dan
 */

#ifndef NUMERICS_H_
#define NUMERICS_H_

#include "globalswitches.h"
#include<vector>
#include<stdlib.h>
#include<iostream>
#include<stdint.h>

using namespace std;

//#include "array.h"
//#include <blitz/array.h>

typedef intptr_t integer;
typedef uintptr_t uinteger;

#ifdef USE_DOUBLE_PRECISION
typedef double scalar;
#else
typedef float scalar;
#endif

/*
 * SafeArray class is essentially just a proxy for an array pointer T* but
 * is initialised with a name and a length. If it's initialised from a vector<T>
 * then the length is checked against the size of the vector at initialisation
 * time, so if it is too big it raises an error. If bounds checking is activated
 * via the #define statement below, then array accesses will be bounds checked
 * too, which is slower but can be used for debugging purposes.
 *
 * The advantage of this class, rather than just using stl::vector is that
 * bounds checking can be switched on and off with a compiler switch, whereas
 * with stl::vector you use .at() for bounds checked access, and [] for
 * unchecked access. It can also be used to construct views of an array, i.e.
 * a view into the data of an stl::vector but with an offset. A more
 * sophisticated solution for this would be to use Boost.MultiArray, blitz++ or
 * John Bowman's array.h, but these are all much more complex and not necessary
 * for this particular case.
 *
 * The advantage over the previous Array class is that we primarily use
 * stl::vector which has computational advantages over Array, and we can easily
 * construct views as well.
 *
 * For safer programming, it might be advisable to define 2D/etc. arrays, in
 * which case looking at one of the three solutions above is probably more
 * sensible than writing our own.
 */

#include "log.h"

template <class T>
class SafeArray
{
private:
    T* base;
    size_t length;
    const char *name;
public:
    // Initialise from a pointer - not safe, the char *name is used for debugging
    //SafeArray(T* base, size_t length, char *name) : base(base), length(length), name(name) {};
    // Initialise from a vector, with an optional offset
    SafeArray(vector<T> &vec, const char *vecname);
    SafeArray(vector<T> &vec, size_t offset, const char *vecname);
    inline T& operator[](const integer i) const;
    // TODO: add an .at() method which always does bounds checking
};

template <class T>
SafeArray<T>::SafeArray(vector<T> &vec, const char *vecname)
{
    name = vecname;
    base = &vec[0];
    length = vec.size();
}

template <class T>
SafeArray<T>::SafeArray(vector<T> &vec, size_t offset, const char *vecname)
{
    name = vecname;
    if(offset<0 || offset>vec.size())
    {
        Error("Cannot create SafeArray with offset %d from vector of "
              "size %d (%s)\n", (int)offset, (int)vec.size(), name);
        abort();
    }
    base = &vec[offset];
    length = vec.size()-offset;
}

template <class T>
inline T& SafeArray<T>::operator[](const integer i) const
{
#ifdef SAFEARRAY_BOUNDSCHECKING
    if(i<0 || i>=(integer)length)
    {
        cerr << "Array index " << i << " out of bounds (" << name << ")" << endl;
        abort();
    }
#endif
    return base[i];
}

#endif /* NUMERICS_H_ */
