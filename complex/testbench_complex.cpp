/*************************  testbench_complex.cpp   ***************************
* Author:        Agner Fog
* Date created:  2019-07-10
* Last modified: 2022-07-20
* Version:       2.02.00
* Project:       Testbench for complex number vectors in complexvec1 using vector class library
* Description:
* Compile and run this program to test operators and functions in complex vector library
* This file contains test cases for general operators and functions.
* Each function or operator is tested with many different combinations of input data.
*
* Instructions:
* The following parameters must be defined on the command line or added in the
* top of this file:
*
* vtype: Vector type to test
* testcase: A number defining a function or operator to test. See the cases in this file.
* seed:  Seed for random number generator. May be any integer
*
* Compile with any compiler supported by VCL.
* Specify the desired instruction set and optimization options as parameters
* to the compiler.
*
* (c) Copyright 2019-2022 Agner Fog.
* Apache license 2.0
******************************************************************************

Test cases:
1:   operator +
2:   operator -
3:   operator *     (expect loss of precision)
4:   operator /     (expect loss of precision)
5:   vector + real
6:   vector - real
7:   real - vector
8:   vector * real
9:   vector / real
10:  real / vector
11:  operator ~: complex conjugate
20:  constructor from real scalar
21:  constructor from single real+imag pair
22:  constructor from multiple real+imag pairs
23:  constructor from complex scalar
24:  constructor from two half size complex vectors
25:  constructor from four complex scalars
26:  constructor from eight complex scalars
27:  constructor from 16 complex scalars
30:  real(): Get real part of complex scalar
31:  imag(): Get imaginary part of complex scalar
32:  real(): Get real parts of complex vector
33:  imag(): Get imaginary parts of complex scalar
34:  interleave_c: interleave real and imag parts
35:  interleave_c2: interleave real and imag parts
36:  interleave_c4: interleave real and imag parts
37:  get_low(): Get first half size vector
38:  get_high(): Get second half size vector
39:  extract: Get one complex element
40:  operator ==
41:  operator !=
42:  abs_greater
43:  abs_less
50:  select
51:  to_float: different precision
52:  to_double: higher precision
55:  abs
56:  csqrt
60:  chorizontal_add
103: mul_accurate
104: div_accurate
203: mul with accurate reference function
204: div with accurate reference function
500: cexp
501: clog

*****************************************************************************/

#include <stdio.h>
#include <cmath>
#if defined (__linux__) && !defined(__LP64__)
#include <fpu_control.h>      // set floating point control word
#endif

#define MAX_VECTOR_SIZE 512
//#define MAX_VECTOR_SIZE 256
//#define MAX_VECTOR_SIZE 128
//# define __F16C__

#ifndef INSTRSET
#define INSTRSET  10
#endif
#include <vectorclass.h>

#include <vectormath_exp.h>                 // exp and log functions
#include <vectormath_trig.h>                // trigonometric functions

#include "complexvec1.h"                    // complex number vectors
#include "complexvecfp16.h"                    // half precision complex number vectors


#ifndef testcase 
// ---------------------------------------------------------------------------
//            Specify input parameters here if running from an IDE
// ----------------------------------------------------------------------------

#define testcase 1

#define vtype Complex4f

#define seed 1


#endif  // testcase

#define TESTNAN 
// ----------------------------------------------------------------------------
//             Declarations
// ----------------------------------------------------------------------------

// dummy vectors used for getting vector types and element type
vtype dummyc;                                    // complex vector type
typedef decltype(dummyc.to_vector()) wtype;      // corresponding normal vector type
wtype dummyv;
typedef decltype(dummyv.extract(0)) ST;          // scalar type
ST a0, a1, b0;                                   // scalar operands
const int maxvectorsize = 32;                    // max number of elements in a vector
ST oplist[maxvectorsize];                        // operand vector
int k0;                                          // copy of vector index

struct Cpair {                                   // complex number structure for reference function
    ST re;                                       // real part
    ST im;                                       // imaginary part
    Cpair(){}                                    // default constructor
    Cpair(ST r, ST i) {                          // constructor from elements
        re = r;  im = i;
    }
    Cpair & load(ST * p) {                       // load data from array
        re = *p;  im = *(p+1);
        return *this;
    }
    void store(ST * p) {                         // save data to array
        *p = re;  *(p+1) = im;
    }
};

#ifdef COMPLEXVECFP16_H
// type Float16 emulates _Float16 if _Float16 not available
Float16 signbit_(Float16 x) {
    union { Float16 f; uint16_t i; } u;
    u.f = x; 
    return u.i >> 15 != 0 ? -1.f : 1.f;
}
#endif
float signbit_(float x) {
    union { float f; uint32_t i; } u;
    u.f = x; 
    return u.i >> 31 != 0 ? -1.f : 1.f;
}
double signbit_(double x) {
    union { double f; uint64_t i; } u;
    u.f = x; 
    return u.i >> 63 != 0 ? -1. : 1.;
}

#if defined(__FMA__) || INSTRSET >= 8            // FMA3 or AVX2
// define FMA for reference functions
float fmadd_scalar(float a, float b, float c) {
    return _mm_cvtss_f32(_mm_fmadd_ss(_mm_load_ss(&a), _mm_load_ss(&b), _mm_load_ss(&c)));
}

double fmadd_scalar(double a, double b, double c) {
    return _mm_cvtsd_f64(_mm_fmadd_sd(_mm_load_sd(&a), _mm_load_sd(&b), _mm_load_sd(&c)));
}
#endif


/************************************************************************
*
*                          Test cases
*
************************************************************************/

#if   testcase == 1    // +
vtype testFunction(vtype const& a, vtype const& b) { 
    return a + b; 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair (a.re + b.re, a.im + b.im);
}

#elif testcase == 2    // - 
inline vtype testFunction(vtype const& a, vtype const& b) { 
    return a - b; 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair (a.re - b.re, a.im - b.im);
}

#elif testcase == 3    // * 
inline vtype testFunction(vtype const& a, vtype const& b) {
    return a * b;
    //return b * a;  // note: a*b and b*a give different precision
}
Cpair referenceFunction(Cpair a, Cpair b) {
    ST x, y;
#if 0 // defined(__FMA__) || INSTRSET >= 8            // FMA3 or AVX2
    x = fmadd_scalar(a.re, b.re, - a.im * b.im); // calculate with same precision as test function
    y = fmadd_scalar(a.re, b.im,   b.re * a.im);
#else
    x = a.re * b.re - a.im * b.im;
    y = a.re * b.im + b.re * a.im;
#endif
    return Cpair (x, y);
}
#define FACCURACY 1000   // loss of precision and difference with reference function

#elif testcase == 4    // /
inline vtype testFunction(vtype const& a, vtype const& b) {
    return a / b;
}
Cpair referenceFunction(Cpair a, Cpair b) {
    ST x, y;
#if 0// defined(__FMA__) || INSTRSET >= 8              // FMA3 or AVX2
    x = fmadd_scalar(a.im, b.im ,   a.re * b.re);  // calculate with same precision as test function
    y = fmadd_scalar(a.im, b.re , - a.re * b.im);
#else
    x = a.im * b.im + a.re * b.re;
    y = a.im * b.re - a.re * b.im;
#endif
    ST ssum = b.re * b.re + b.im * b.im;
    return Cpair (x / ssum, y / ssum);
}
#define FACCURACY 10000   // loss of precision and difference with reference function

#elif testcase == 5    // vector + real

inline vtype testFunction(vtype const& a, vtype const& b) { 
    a0 = a.to_vector()[0];
    return a0 + b; 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair (a0 + b.re, b.im);
}

#elif testcase == 6    // vector - real
inline vtype testFunction(vtype const& a, vtype const& b) { 
    a0 = a.to_vector()[0];
    return b - a0; 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair (b.re - a0, b.im);
} 

#elif testcase == 7    // real - vector
inline vtype testFunction(vtype const& a, vtype const& b) { 
    a0 = a.to_vector()[0];
    return a0 - b; 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair (a0 - b.re, - b.im);
} 

#elif testcase == 8    // vector * real
inline vtype testFunction(vtype const& a, vtype const& b) { 
    a0 = a.to_vector()[0];
    return a0 * b; 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair (a0 * b.re, a0 * b.im);
}

#elif testcase == 9    // vector / real
inline vtype testFunction(vtype const& a, vtype const& b) { 
    a0 = a.to_vector()[0];
    return b / a0; 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair (b.re / a0, b.im / a0);
}

#elif testcase == 10   // real / vector
inline vtype testFunction(vtype const& a, vtype const& b) { 
    a0 = a.to_vector()[0];
    return a0 / b; 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    ST c = a0 / (b.re*b.re + b.im*b.im);
    return Cpair (c * b.re, -c * b.im);
}
#define FACCURACY 200   // loss of precision and difference with reference function

#elif testcase == 11    // ~ complex conjugate
inline vtype testFunction(vtype const& a, vtype const& b) { 
    return ~a; 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair (a.re, - a.im);
}

#elif testcase == 20    // constructor from real scalar
inline vtype testFunction(vtype const& a, vtype const& b) { 
    a0 = a.to_vector()[0];
    return vtype(a0); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair (a0, 0);
}

#elif testcase == 21    // constructor from single real+imag pair
inline vtype testFunction(vtype const& a, vtype const& b) { 
    a0 = a.to_vector()[0];
    a1 = a.to_vector()[1];
    return vtype(a0, a1); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair (a0, a1);
}

#elif testcase == 22    // constructor from multiple real+imag pairs
template <typename V>
inline V testFunction(V const& a, V const& b) { 
    a.store(oplist);
    if constexpr (V::size() == 1) {
        return V(oplist[0], oplist[1]);
    }
    else if constexpr (V::size() == 2) {
        return V(oplist[0], oplist[1], oplist[2], oplist[3]);
    }
    else if constexpr (V::size() == 4) {
        return V(oplist[0], oplist[1], oplist[2], oplist[3], oplist[4], oplist[5], oplist[6], oplist[7]);
    }
    else if constexpr (V::size() == 8) {
        return V(oplist[0], oplist[1], oplist[2], oplist[3], oplist[4], oplist[5], oplist[6], oplist[7],
            oplist[8], oplist[9], oplist[10], oplist[11], oplist[12], oplist[13], oplist[14], oplist[15]);
    }
    else if constexpr (V::size() == 16) {
        return V(oplist[0], oplist[1], oplist[2], oplist[3], oplist[4], oplist[5], oplist[6], oplist[7],
            oplist[8], oplist[9], oplist[10], oplist[11], oplist[12], oplist[13], oplist[14], oplist[15],
            oplist[16], oplist[17], oplist[18], oplist[19], oplist[20], oplist[21], oplist[22], oplist[23],
            oplist[24], oplist[25], oplist[26], oplist[27], oplist[28], oplist[29], oplist[30], oplist[31]);
    }
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return a;
}

#elif testcase == 23    // constructor from complex scalar
inline vtype testFunction(vtype const& a, vtype const& b) { 
    a0 = a.to_vector()[0];
    a1 = a.to_vector()[1];
    auto s = a.extract(0);   // complex scalar
    return vtype(s); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair (a0, a1);
}

#elif testcase == 24    // constructor from two half size complex vectors
inline vtype testFunction(vtype const& a, vtype const& b) { 
    a.store(oplist);
    decltype(a.get_low()) lo, hi; // half size vectors
    lo.load(oplist);
    hi.load(oplist + 2*lo.size());
    return vtype(lo, hi); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return a;
}

#elif testcase == 25   //  constructor from four complex scalars
inline vtype testFunction(vtype const& a, vtype const& b) {
    return vtype(a.extract(0), a.extract(1), a.extract(2), a.extract(3));
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return a;
}

#elif testcase == 26   //  constructor from eight complex scalars
inline vtype testFunction(vtype const& a, vtype const& b) {
    return vtype(a.extract(0), a.extract(1), a.extract(2), a.extract(3), a.extract(4), a.extract(5), a.extract(6), a.extract(7));
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return a;
}

#elif testcase == 27   //  constructor from 16 complex scalars
inline vtype testFunction(vtype const& a, vtype const& b) {
    return vtype(a.extract(0), a.extract(1), a.extract(2), a.extract(3), a.extract(4), a.extract(5), a.extract(6), a.extract(7),
        a.extract(8), a.extract(9), a.extract(10), a.extract(11), a.extract(12), a.extract(13), a.extract(14), a.extract(15));
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return a;
}

#elif testcase == 30    // real(): Get real part of complex scalar
inline vtype testFunction(vtype const& a, vtype const& b) { 
    ST r = a.real();
    return vtype(r, 0); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair(a.re, 0);
}

#elif testcase == 31    // imag(): Get imaginary part of complex scalar
inline vtype testFunction(vtype const& a, vtype const& b) { 
    ST i = a.imag();
    return vtype(i, 0); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return Cpair(a.im, 0);
}

#elif testcase == 32    // real(): Get real parts of complex vector
vtype aa0;

inline vtype testFunction(vtype const& a, vtype const& b) { 
    aa0 = a;      // save whole vector
    ST rr[vtype::size()*2 + 8] = {0};
    auto r = a.real();
    r.store(rr);
    return vtype().load(rr);
}
Cpair referenceFunction(Cpair a, Cpair b) {
    if (k0 < vtype::size()) {
        auto aa1 = aa0.to_vector();
        Cpair rr(aa1[k0*2], aa1[k0*2+2]); 
        return rr;
    }
    else return Cpair(0, 0);
}

#elif testcase == 33    // imag(): Get imag parts of complex vector
vtype aa0;

inline vtype testFunction(vtype const& a, vtype const& b) { 
    aa0 = a;      // save whole vector
    ST rr[vtype::size()*2 + 8] = {0};
    auto r = a.imag();
    r.store(rr);
    return vtype().load(rr);
}
Cpair referenceFunction(Cpair a, Cpair b) {
    if (k0 < vtype::size()) {
        auto aa1 = aa0.to_vector();
        Cpair rr(aa1[k0*2+1], aa1[k0*2+3]); 
        return rr;
    }
    else return Cpair(0, 0);
}

#elif testcase == 34    // interleave_c: interleave real and imag parts
inline vtype testFunction(vtype const& a, vtype const& b) {
    auto re = a.real();
    auto im = a.imag();
    vtype c = interleave_c(re,im);
    return c;
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return a;
}

#elif testcase == 35    // interleave_c2: interleave real and imag parts
inline vtype testFunction(vtype const& a, vtype const& b) {
    auto re = a.real();
    auto im = a.imag();
    vtype c = interleave_c2(re,im);
    return c;
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return a;
}

#elif testcase == 36    // interleave_c4: interleave real and imag parts
inline vtype testFunction(vtype const& a, vtype const& b) {
    auto re = a.real();
    auto im = a.imag();
    vtype c = interleave_c4(re,im);
    return c;
}
Cpair referenceFunction(Cpair a, Cpair b) {
    return a;
}

#elif testcase == 37    // get_low(): Get first half-size vector
inline vtype testFunction(vtype const& a, vtype const& b) {
    a.store(oplist);
    auto l = a.get_low();
    return vtype(l, l); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    int k2 = k0 % vtype::size();             // get lower half twice
    return Cpair(oplist[k2], oplist[k2+1]);
}

#elif testcase == 38    // get_high(): Get second half-size vector
inline vtype testFunction(vtype const& a, vtype const& b) {
    a.store(oplist);
    auto h = a.get_high();
    return vtype(h, h); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    int k2 = (k0 % vtype::size()) + vtype::size();// get upper half twice
    return Cpair(oplist[k2], oplist[k2+1]);
}

#elif testcase == 39    // extract: Get one complex element
inline vtype testFunction(vtype const& a, vtype const& b) {
    a.store(oplist);
    a0 = b.to_vector()[0];
    uint32_t i = uint32_t(a0) % vtype::size();
    auto e = a.extract(i);
    return vtype(e); 
} 
Cpair referenceFunction(Cpair a, Cpair b) {
    uint32_t i = uint32_t(a0) % vtype::size();
    return Cpair(oplist[i*2], oplist[i*2+1]);
}

#elif testcase == 49    // insert: Insert one complex element
inline vtype testFunction(vtype & a, vtype const& b) {
    a.store(oplist);
    a0 = b.to_vector()[0];    // value to insert
    a1 = b.to_vector()[1];
    b0 = b.to_vector()[vtype::size() > 1 ? 2 : 1]; // index
    uint32_t i = uint32_t(b0) % vtype::size();
    typedef decltype(a.extract(0)) sctype; // type of complex scalar
    a.insert(i, sctype(a0,a1));
    return a; 
} 
Cpair referenceFunction(Cpair a, Cpair b) {
    uint32_t i = uint32_t(b0) % vtype::size();
    if (k0/2 == i) {
        return Cpair(a0, a1);
    }
    else {    
        return Cpair(oplist[k0], oplist[k0+1]);
    }
}


#elif testcase == 40    // operator ==
inline vtype testFunction(vtype const& a, vtype const& b) {
    auto eq = a == b;   // boolean vector
    return select(eq, vtype(1.f), vtype(0.f));
}
Cpair referenceFunction(Cpair a, Cpair b) {
    if (a.re == b.re && a.im == b.im) {
        return Cpair(1,0);
    }
    else {
        return Cpair(0,0);
    }
}

#elif testcase == 41 || testcase == 50   // operator !=, select
inline vtype testFunction(vtype const& a, vtype const& b) {
    auto eq = a != b;   // boolean vector
    return select(eq, vtype(1.f), vtype(0.f));
}
Cpair referenceFunction(Cpair a, Cpair b) {
    if (a.re != b.re || a.im != b.im) {
        return Cpair(1,0);
    }
    else {
        return Cpair(0,0);
    }
}

#elif testcase == 42    // abs_greater
inline vtype testFunction(vtype const& a, vtype const& b) {
    auto g = abs_greater(a,b);   // boolean vector
    return select(g, vtype(1.f), vtype(0.f));
}
Cpair referenceFunction(Cpair a, Cpair b) {
    double absa2 = a.re * a.re + a.im * a.im;
    double absb2 = b.re * b.re + b.im * b.im; 
    if (absa2 > absb2) {
        return Cpair(1,0);
    }
    else {
        return Cpair(0,0);
    }
}

#elif testcase == 43    // abs_less
inline vtype testFunction(vtype const& a, vtype const& b) {
    auto g = abs_less(a,b);   // boolean vector
    return select(g, vtype(1.f), vtype(0.f));
}
Cpair referenceFunction(Cpair a, Cpair b) {
    double absa2 = a.re * a.re + a.im * a.im;
    double absb2 = b.re * b.re + b.im * b.im; 
    if (absa2 < absb2) {
        return Cpair(1,0);
    }
    else {
        return Cpair(0,0);
    }
}



#elif testcase == 51   // to_float: lower precision
float fllist[maxvectorsize];

inline vtype testFunction(vtype const& a, vtype const& b) {
    to_float(a).store(fllist);
    return a;
}
Cpair referenceFunction(Cpair a, Cpair b) {
    if (float(a.re) != float(a.re) || float(a.im) != float(a.im)) return a; // ignore NAN
    if (float(a.re) != fllist[k0] || float(a.im) != fllist[k0+1]) return Cpair(0,0);
    else return a;
}

#elif testcase == 52   // to_double: higher precision
double dbllist[maxvectorsize];

inline vtype testFunction(vtype const& a, vtype const& b) {
    to_double(a).store(dbllist);
    return a;
}
Cpair referenceFunction(Cpair a, Cpair b) {
    if (double(a.re) != double(a.re) || double(a.im) != double(a.im)) return a; // ignore NAN
    if (double(a.re) != dbllist[k0] || double(a.im) != dbllist[k0+1]) return Cpair(0,0);
    else return a;
}

#elif testcase == 55    // abs
inline vtype testFunction(vtype const& a, vtype const& b) { 
    return abs(a); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    ST ssum = a.re * a.re + a.im * a.im;
    if (sizeof(ST) > 4) return Cpair (std::sqrt(double(ssum)), 0);
    else return Cpair (std::sqrt(float(ssum)), 0);
}

#elif testcase == 56    // csqrt
// To do: This test does not work. fix it!
inline vtype testFunction(vtype const& a, vtype const& b) {
    vtype r = csqrt(b);
    // accuracy problem when input is extremely small
    ST rr[vtype::size()*2] = {0};
    r.store(rr);
    float limit = 1.E-20f;
    for (int i=0; i<vtype::size(); i++) {
        if (abs(float(b.extract(i).real())) < limit 
        &&  abs(float(b.extract(i).imag())) < limit) {
            rr[2*i] = ST(0.f);  rr[2*i+1] = ST(0.f);        
        }    
    }
    return r;
}
Cpair referenceFunction(Cpair a, Cpair b) {
    a0 = b.re; a1 = b.im;
    // accuracy problem when input is extremely small
    float limit = 1.E-20f;
    if (std::abs(float(a0)) < limit && std::abs(float(a1)) < limit) {
        return Cpair(0.f, 0.f);
    }
#ifdef __AVX512FP16__
    double n = std::sqrt(double(b.re * b.re + b.im * b.im));
#else   // half precision emulation calculates intermediates with higher precision
    double n = std::sqrt(double(b.re) * double(b.re) + double(b.im) * double(b.im));
#endif
    Cpair r = Cpair ((ST)(std::sqrt((n + double(b.re))*0.5)), ST((float)signbit_(b.im) * std::sqrt((n - double(b.re))*0.5)));
    return r;
}
#define FACCURACY 100000    // loss of precision and difference with reference function
#define IGNORE_SUBNORMAL

#elif testcase == 60    //  chorizontal_add
inline vtype testFunction(vtype const& a, vtype const& b) {
    a.store(oplist);
    vtype r = chorizontal_add(a); 
    return vtype(r);
}
Cpair referenceFunction(Cpair a, Cpair b) {
    double resum = 0, imsum = 0;
    for (int i = 0; i < vtype::size(); i++) {
        resum += oplist[2*i];
        imsum += oplist[2*i + 1];
    }
    Cpair r = Cpair(ST(resum), ST(imsum));
    return r;
}
#define FACCURACY 10    // loss of precision and difference with reference function

#elif testcase == 103   // mul_accurate. Complex multiplication without loss of precision
// Note: The reference function may not be not sufficiently precise
inline vtype testFunction(vtype const& a, vtype const& b) {
    return mul_accurate(a, b);
}
Cpair referenceFunction(Cpair a, Cpair b) {
    long double are = double(a.re), aim = double(a.im), bre = double(b.re), bim = double(b.im);
    long double x = are * bre - aim * bim;
    long double y = are * bim + bre * aim;
    return Cpair (ST(double(x)), ST(double(y)));
}
#define FACCURACY 1  // reference function is not precise if long double not supported by compiler

#elif testcase == 203   // normal multiply with precise reference function
inline vtype testFunction(vtype const& a, vtype const& b) {
    return a * b;
}
Cpair referenceFunction(Cpair a, Cpair b) {
    //a = Cpair(1.01, 1.02);
   // b = Cpair(1.01, 1.0); 

    long double are = a.re, aim = a.im, bre = b.re, bim = b.im;
    long double x = are * bre - aim * bim;
    long double y = are * bim + bre * aim;
    return Cpair (ST(x), ST(y));
}
#define FACCURACY 1  // reference function is not precise if compiler does not support long double

#elif testcase == 104   // div_accurate. Complex division without loss of precision
inline vtype testFunction(vtype const& a, vtype const& b) {
    return div_accurate(a, b);
}
Cpair referenceFunction(Cpair a, Cpair b) {
    long double are = a.re, aim = a.im, bre = b.re, bim = b.im;
    long double x = are * bre + aim * bim;
    long double y = bre * aim - are * bim;
    long double ssum = bre * bre + bim * bim;
    return Cpair (ST(x/ssum), ST(y/ssum));
}
#define FACCURACY 10  // reference function is not precise?

#elif testcase == 204   // normal division with precise reference function
inline vtype testFunction(vtype const& a, vtype const& b) {
    return a / b;
}
Cpair referenceFunction(Cpair a, Cpair b) {
    long double are = a.re, aim = a.im, bre = b.re, bim = b.im;
    long double x = are * bre + aim * bim;
    long double y = bre * aim - are * bim;
    long double ssum = bre * bre + bim * bim;
    return Cpair (ST(x/ssum), ST(y/ssum));
}
#define FACCURACY 0  // reference function is not precise if compiler does not support long double

#elif testcase == 500    // cexp. Complex exponential function
inline vtype testFunction(vtype const& a, vtype const& b) { 
    return cexp(a); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    long double re = a.re;
    long double im = a.im;
    long double e = expl(re);
    long double s = sinl(im);
    long double c = cosl(im);
    return Cpair(ST(double(e*c)), ST(double(e*s)));
}
#define FACCURACY 4     // test absolute error here, not relative error. 
#define ABS_ERR_EXP     // The relative error on sin and cos is unlimited when close to 0

// This test case generates link error in Intel compiler ICPX version 2022.1.0
// for missing external __extendhfxf2 (probably extend half precision to long double precision)
// It works with all other compilers

#elif testcase == 501   // clog. Complex logarithm function

inline vtype testFunction(vtype const& a, vtype const& b) { 
    return clog(a); 
}
Cpair referenceFunction(Cpair a, Cpair b) {
    long double re = a.re;
    long double im = a.im;
    long double x = 0.5 * log(re*re + im*im);
    long double y = atan2l(im, re);
    return Cpair(ST(x), ST(y));
}
#define FACCURACY 20
#define ABS_ERR_EXP


#else
// End of test cases
#error unknown test case
#endif


// ----------------------------------------------------------------------------
//                           Overhead functions
// ----------------------------------------------------------------------------

const int maxerrors = 10;      // maximum errors to report
int numerr = 0;                // count errors

// type-specific load function
template <typename T, typename E>
void loadData(T & x, E const* p) {
    x.load(p);
}

template <typename T>
inline void loadData(T & x, bool const* p) {
    for (int i = 0; i < x.size(); i++) {
        x.insert(i, p[i]);  // bool vectors have no load function
    }
}


// type-specific printing functions

void printVal(float x) {
    printf("%14.10G", x);
}

void printVal(double x) {
    printf("%14.10G", x);
}

void printVal(bool x) {
    printf("%i", (int)x);
}

// Random number generator
class ranGen {
    // parameters for multiply-with-carry generator
    uint64_t x, carry;
public:
    ranGen(int Seed) {  // constructor
        x = Seed;  carry = 1765;  //initialize with seed
        next();  next();
    }
    uint32_t next() {  // get next random number, using multiply-with-carry method
        const uint32_t fac = 3947008974u;
        x = x * fac + carry;
        carry = x >> 32;
        x = uint32_t(x);
        return uint32_t(x);
    }
};

template <typename T>  // get random number of type T
T get_random(ranGen & rangen) {
    if constexpr (sizeof(T) > 4) {
        return (T)(double)rangen.next();
    }
    else {
        return (T)(float)rangen.next();
    }
}

template <>  // special case uint64_t
uint64_t get_random<uint64_t>(ranGen & rangen) {
    uint64_t xx;
    xx = (uint64_t)rangen.next() << 32;
    xx |= rangen.next();
    return xx;
}

template <>  // special case int64_t
int64_t get_random<int64_t>(ranGen & rangen) {
    return (int64_t)get_random<uint64_t>(rangen);
}

template <>  // special case float
float get_random<float>(ranGen & rangen) {
    union Uif {
        uint32_t i;
        float f;
    };
    Uif u1, u2;
    uint32_t r = rangen.next();                  // get 32 random bits
    // Insert exponent and random mantissa to get random number in the interval 1 <= x < 2
    // Subtract 1.0 if next bit is 0, or 1.0 - 2^-24 = 0.99999994f if next bit is 1
    u1.i = 0x3F800000 - ((r >> 8) & 1);          // bit 8
    u2.i = (r >> 9) | 0x3F800000;                // bit 9 - 31
    return u2.f - u1.f;
}

template <>  // special case float
double get_random<double>(ranGen & rangen) {
    union Uqd {
        uint64_t q;
        double d;
    };
    Uqd u1;
    uint64_t r = get_random<uint64_t>(rangen);   // get 64 random bits
    // Insert exponent and random mantissa to get random number in the interval 1 <= x < 2,
    // then subtract 1.0 to get the interval 0 <= x < 1.
    u1.q = (r >> 12) | 0x3FF0000000000000;       // bit 12 - 63
    return u1.d - 1.0;
}
template <>  // special case bool
bool get_random<bool>(ranGen & rangen) {
    return (rangen.next() & 1) != 0;
}


// make random number generator instance
ranGen ran(seed);

#ifdef VECTORFP16_H
// bit_cast function to make special values
Float16 bit_casth(uint16_t x) {  // uint64_t -> double
    union {
        uint16_t i;
        Float16 f;
    } u;
    u.i = x;
    return u.f;
}
#endif

float bit_castf(uint32_t x) {  // uint64_t -> double
    union {
        uint32_t i;
        float f;
    } u;
    u.i = x;
    return u.f;
}

double bit_castd(uint64_t x) {  // uint32_t -> float
    union {
        uint64_t i;
        double f;
    } u;
    u.i = x;
    return u.f;
}


// template to generate list of testdata
template <typename T>
class TestData {
public:
    enum LS {
        // define array size. Must be a multiple of vector size:
        listsize = 1024
    };
    TestData(int ab) {                          // constructor.ab = 1 for a, 2 for b
        int i = 0;                              // loop counter
        if (T(1.1f) != T(1)) {
            // floating point type
#ifdef TESTNAN   // test also with NAN, INF, and other special data
#ifdef VECTORFP16_H
            // additional special values, float16:
            if constexpr (sizeof(ST) == 2) {
                if (ab == 1) {  // adata
                    list[i++] = (T)1.01f;      // provoke high loss of precision in multiplication
                    list[i++] = (T)(1.02f);
                    list[i++] = (T)1.01f;
                    list[i++] = (T)1.02f;
                }
                else {   // bdata
                    list[i++] = (T)1.01f;
                    list[i++] = (T)(1.0f);
                    list[i++] = (T)1.0f; 
                    list[i++] = (T)(-1.01f);
                }
                list[i++] = (T)bit_castf(0x8000);   // -0
#if testcase != 4 && testcase != 104 && testcase != 204 // avoid division overflow
                list[i++] = (T)bit_casth(0x0080);   // smallest positive normal number
                list[i++] = (T)bit_casth(0x8080);   // largest negative normal number
#endif
                list[i++] = (T)bit_casth(0x3BFF);   // next before 1.0, 0
                list[i++] = (T)bit_casth(0x3C01);   // nextafter 1.0, 2
                list[i++] = (T)bit_casth(0x7C00);   // inf
                list[i++] = (T)bit_casth(0xFC00);   // -inf
                list[i++] = (T)bit_casth(0x7FF0);   // nan
            }
#endif
            if constexpr (sizeof(ST) == 4) {
                // additional special values, float:
                if (ab == 1) {  // adata
                    list[i++] = (T)1.0001f;      // provoke high loss of precision in multiplication
                    list[i++] = (T)(-1.0002f);
                    list[i++] = (T)1.0000f;
                    list[i++] = (T)1.0001f;
                }
                else {   // bdata
                    list[i++] = (T)1.0000f;
                    list[i++] = (T)1.0001f;
                    list[i++] = (T)1.0001f; 
                    list[i++] = (T)(-1.0002f);
                }                
                list[i++] = (T)bit_castf(0x80000000);   // -0
#if testcase != 4 && testcase != 104 && testcase != 204 // avoid division overflow
                list[i++] = (T)bit_castf(0x00800000);   // smallest positive normal number
                list[i++] = (T)bit_castf(0x80800000);   // largest negative normal number
#endif
                list[i++] = (T)bit_castf(0x3F7FFFFF);   // nextafter 1.0, 0
                list[i++] = (T)bit_castf(0x3F800001);   // nextafter 1.0, 2
                list[i++] = (T)bit_castf(0x7F800000);   // inf
                list[i++] = (T)bit_castf(0xFF800000);   // -inf
                list[i++] = (T)bit_castf(0x7FF00000);   // nan
            }
            else { // double
                if (ab == 1) {  // adata
                    list[i++] = (T)1.000001;      // provoke high loss of precision in multiplication
                    list[i++] = (T)(-1.000002);
                    list[i++] = (T)1.000000;
                    list[i++] = (T)1.000001;
                }
                else {   // bdata
                    list[i++] = (T)1.000000;
                    list[i++] = (T)1.000001;
                    list[i++] = (T)1.000001; 
                    list[i++] = (T)(-1.000002);
                } 
                list[i++] = (T)1.000001;      // provoke high loss of precision
                list[i++] = (T)(-1.000002);
                list[i++] = (T)1.000000;
                list[i++] = (T)1.000001;
                list[i++] = (T)bit_castd(0x8000000000000000);   // -0
#if testcase != 4 && testcase != 104 && testcase != 204 // avoid division overflow
                list[i++] = (T)bit_castd(0x0010000000000000);   // smallest positive normal number
                list[i++] = (T)bit_castd(0x8010000000000000);   // largest negative normal number
#endif
                list[i++] = (T)bit_castd(0x3FEFFFFFFFFFFFFF);   // nextafter 1.0, 0
                list[i++] = (T)bit_castd(0x3FF0000000000001);   // nextafter 1.0, 2
                list[i++] = (T)bit_castd(0x7FF0000000000000);   // inf
                list[i++] = (T)bit_castd(0xFFF0000000000000);   // -inf
                list[i++] = (T)bit_castd(0x7FFC000000000000);   // nan
            }
#endif
            // fill boundary data into array
            for (; i < 20; i++) {
                list[i] = T(T(i - 4.f) * T(0.25));
            }
            // fill random data into rest of array
            for (; i < listsize; i++) {
                if (vtype::elementtype() > 3) {                
                    list[i] = get_random<T>(ran) * (T)100;
                }
                else {  // bool
                    list[i] = get_random<T>(ran);
                }
            }
        }
        else {
            // integer type
            // fill boundary data into array
            for (i = 0; i < 6; i++) {
                list[i] = T(i - 2.f);
            }
            // data near mid-point of unsigned integers, or overflow point of signed integers:
            uint64_t m = (uint64_t(1) << (sizeof(T) * 8 - 1)) - 2;
            for (; i < 11; i++) {
                list[i] = T(float(m++));
            }
            // fill random data into rest of array
            for (; i < listsize; i++) {
                list[i] = get_random<T>(ran);
            }
        }
    } 
    T list[listsize];                       // array of test data
    int size() {                            // get list size
        return listsize;
    }
};


// get value of least significant bit
Float16 delta_unit(Float16 x) {
    union {
        Float16 f;
        uint32_t i;
    } u;
    x = std::fabs(float(x));
    Vec8h xv = Vec8h(x);
    if (!(is_finite(xv)[0])) return 1.f;
    if (float(x) == 0.f || is_subnormal(xv)[0]) {
        u.i = 0x0400;  // smallest positive normal number
        return u.f;
    }
    Float16 x1 = x;
    u.f = x;
    u.i++;
    return u.f - x1;
}

float delta_unit(float x) {
    union {
        float f;
        uint32_t i;
    } u;
    x = std::fabs(x);
    Vec4f xv = Vec4f(x);
    if (!(is_finite(xv)[0])) return 1.f;
    if (x == 0.f || is_subnormal(xv)[0]) {
        u.i = 0x00800000;  // smallest positive normal number
        return u.f;
    }
    float x1 = x;
    u.f = x;
    u.i++;
    return u.f - x1;
}

double delta_unit(double x) {
    union {
        double f;
        uint64_t i;
    } u;
    x = std::fabs(x);
    Vec2d xv = Vec2d(x);
    if (!(is_finite(xv)[0])) return 1.;
    if (x == 0. || is_subnormal(xv)[0]) {
        u.i = 0x0010000000000000;  // smallest positive normal number
        return u.f;
    }
    double x1 = x;
    u.f = x;
    u.i++;
    return u.f - x1;
}


// compare two scalars. return true if different
template <typename T>
inline bool compare_scalars(T const a, T const b) {
    return a == b;
}

// special cases for float and double:
#ifdef VECTORFP16_H

union fp16toint16 {
    Float16 f;
    int16_t i;
};

template <>
inline bool compare_scalars<Float16>(Float16 const a, Float16 const b) {
    if (a == b || (a != a && b != b)) return true; // return true if equal or both are NAN
#ifdef IGNORE_SUBNORMAL
    fp16toint16 aa, bb;
    aa.f = a; bb.f = b;
    if ((aa.i & 0x7C00) == 0 && (bb.i & 0x7C00) == 0) return true;  // both are zero or subnormal
#endif

#ifdef FACCURACY     // accept minor difference
    float dif = std::fabs(float(a - b)) / (float)delta_unit(a);
    if (dif <= FACCURACY) return true;
    printf("\n%.0f ULP ", dif);
#endif
    return false;
}
#endif

template <>
inline bool compare_scalars<float>(float const a, float const b) {
    if (a == b || (a != a && b != b)) return true; // return false if equal or both are NAN
#ifdef FACCURACY     // accept minor difference
    float dif = std::fabs(a - b) / delta_unit(a);
    if (dif <= FACCURACY) return true;
    printf("\n%.0f ULP ", dif);
#endif
    return false;
}

template <>
inline bool compare_scalars<double>(double const a, double const b) {
    if (a == b || (a != a && b != b)) return true; // return false if equal or both are NAN
#ifdef FACCURACY     // accept minor difference
    double dif = std::fabs(a - b) / delta_unit(a);
    if (dif <= FACCURACY) return true;
    printf("\n%.0f ULP ", dif);
#endif
    return false;
}

// compare two vectors. return true if different
template <typename T>
inline bool compare_vectors(T const& a, T const& b) {
    {
        for (int i = 0; i < a.size(); i++) {
            if (!compare_scalars(a[i], b[i])) return false;
        }
    }
    return true;
}

#ifndef FACCURACY
#define FACCURACY 1
#endif

// compare two complex numbers. return true if different
inline ST compare_complex(Cpair const& a, Cpair const& b) {
    ST dif1 = 0, dif2 = 0;
    if (a.re != b.re) {
        dif1 = ST(std::fabs(float(a.re - b.re)) / (float)delta_unit(b.re));
    }
    if (a.im != b.im) {
        dif2 = ST(std::fabs(float(a.im - b.im)) / (float)delta_unit(b.im));
    }
    if (std::isinf((double)dif1) || std::isinf((double)dif2)) { // check for overflow
        if (sizeof(ST) == 2) {  // float16
            if (std::fabs(float(a.re)) > 65000.f && std::fabs(float(b.re)) > 65000.f) dif1 = 0;
            if (std::fabs(float(a.im)) > 65000.f && std::fabs(float(b.im)) > 65000.f) dif2 = 0;
        }
        else if (sizeof(ST) == 4) {  // float
            if (std::fabs(float(a.re)) > 1.E37f && std::fabs(float(b.re)) > 1.E37f) dif1 = 0;
            if (std::fabs(float(a.im)) > 1.E37f && std::fabs(float(b.im)) > 1.E37f) dif2 = 0;
        }
        else {  // double
            if (std::fabs(double(a.re)) > 1.E300 && std::fabs(double(b.re)) > 1.E300) dif1 = 0;
            if (std::fabs(double(a.im)) > 1.E300 && std::fabs(double(b.im)) > 1.E300) dif2 = 0;
        }
    }    
#ifdef ABS_ERR_EXP  // the sin and cos must be checked with absolute errors, not relative errors
    double n = std::fabs(double(b.re)) + std::fabs(double(b.im));
    dif1 = dif1 / (ST)n;  dif2 = dif2 / (ST)n;
#endif
    if (double(dif1) > FACCURACY || double(dif2) > FACCURACY) {
        return dif1 > dif2 ? dif1 : dif2;
    }
    return 0;
}

// program entry
int main() {

#if defined (__linux__) && !defined(__LP64__)
    // Some 32-bit compilers use x87 calculations with long double precision for 
    // the reference function. This may give slightly different results because
    // the value is rounded twice. To get exactly the same value in the test function
    // and the reference function, we change the precision of x87 calculations.
    // (the fpu control function is different in Windows, but the precision is already
    // reduced in Windows anyway)
    fpu_control_t fpcw = 0x27f;
    _FPU_SETCW(fpcw);
#endif

    vtype a(ST(0)), b(ST(0));  // complex vectors for operands
    vtype result;  // complex vector for result

    // make lists of test data
    TestData<ST> adata(1), bdata(2);

    // make list of results
    ST resultlist[maxvectorsize] = {0};

    int i, j, k = 0;   // loop counters
    ST maxdif = 0;    // maximum observed error

    for (i = 0; i < adata.size(); i += wtype::size()) {
        //a.load(adata.list + i);
        loadData(a, adata.list + i);

        for (j = 0; j < bdata.size(); j += wtype::size()) {
            loadData(b, bdata.list + j);

            // function under test:
            result = testFunction(a, b);
            result.store(resultlist);
            
            // result is vector
            for (k = 0; k < vtype::size()*2; k += 2) {
                Cpair aa, bb, rr, xx;
                aa.load(adata.list + i + k);
                bb.load(bdata.list + j + k);
                k0 = k;
                rr = referenceFunction(aa, bb);
                xx.load(resultlist + k);
                ST dif = compare_complex(rr, xx);
                if (dif != ST(0) && numerr < maxerrors) {
                    // values are different. report error 
                    if (++numerr == 1) {
                        printf("\ntest case %i:", testcase);  // print test case first time
                    }
                    printf("\nError at %i, %i, element %i, dif = %.2G:", i, j, k / 2, double(dif));
                    printf("\n(%10.7G, %10.7G) op (%10.7G, %10.7G) -> (%10.7G, %10.7G), expected (%10.7G, %10.7G)",
                        (double)aa.re, (double)aa.im, (double)bb.re, (double)bb.im, (double)xx.re, (double)xx.im, (double)rr.re, (double)rr.im);
                }
                if (abs((double)dif) > (double)maxdif) maxdif = ST(abs((float)dif));
                //if (numerr > maxerrors) exit(1);      // stop after maxerrors
            }
        }
    }

    if (numerr == 0) {
        printf("\nsuccess\n");
    }
    else {
        printf("\nmax error %12.4G\n", double(maxdif));
    }
    printf("\n");

    return numerr;
}
