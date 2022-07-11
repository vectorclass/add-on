/****************************  matrixv.h   ************************************
* Author:        Agner Fog
* Date created:  2022-07-08
* Last modified: 2022-07-11
* Version:       2.02.00
* Description:
* Header file for matrix container class template.
* The MatrixV container can store a matrix as VCL vectors
*
* For further instructions, see containers_manual.pdf
*
* (c) Copyright 2022 Agner Fog.
* Apache License version 2.0 or later.
******************************************************************************/

#ifndef MATRIXV_H
#define MATRIXV_H 20200

#ifndef VECTORCLASS_H
#include "vectorclass.h"
#endif

#ifdef VCL_NAMESPACE         // optional namespace
namespace VCL_NAMESPACE {
#endif


// Container class template for matrix, contained in VCL vectors
// This template is optimized for large matrixes that are accessed mainly by rows.
// Each row is contained in one or more vectors. The template will not pack multiple rows into the same vector.
// V is the vector type to use for storage or a vector size of maximum size. The template
// will automatically use the smallest vector size that can contain a whole row.
template <typename V, int rows, int columns>
class MatrixV {
protected:
    // round up number of columns to nearest power of 2
    static int constexpr rowlength_pow2() {
        if constexpr ((columns & columns - 1) == 0) return columns;  // columns is a power of 2 
        else return 1 << (bit_scan_reverse_const(columns) + 1);      // round up to nearest power of 2
    }    
    // get smallest vector type that can contain a whole row, or full V type if at least one full size vector is needed per row.
    static auto constexpr get_vector_type() {
        static_assert((V::size() & V::size()-1) == 0, "vector size must be a power of 2"); // check that vector size is a power of 2
        static_assert((V::elementtype() & 0x3F) > 2, "booleans not supported");
        constexpr int rl = rowlength_pow2();
        if constexpr (rl * 4 <= V::size() && sizeof(V) >= 64) return decltype(V().get_low().get_low())();  // quarter size vector
        else if constexpr (rl * 2 <= V::size() && sizeof(V) >= 32) return decltype(V().get_low())();  // half size vector
        else return V();       // full size vector    
    }
public:
    // default constructor
    MatrixV() {
        errorfunction = 0;
    }

    // define vector type used for rows
    typedef decltype(get_vector_type()) row_vector_type;

    // define element type
    typedef decltype(V()[0]) etype;

    // get number of rows
    static int constexpr nrows() {
        return rows;
    }

    // get number of columns
    static int constexpr ncolumns() {
        return columns;
    }

    // get number of vectors of type row_vector_type per row. This includes partially used vectors
    static int constexpr vectors_per_row() {
        if constexpr (columns <= row_vector_type::size()) return 1;  // one vector is enough
        else return (columns + row_vector_type::size() - 1) / row_vector_type::size();  // round up to whole number of vectors per row
    }

    // get number of full vectors of type row_vector_type per row
    static int constexpr full_vectors_per_row() {
        return columns / row_vector_type::size();
    }

    // get number of remaining vector elements if columns is not a multiple of vector size
    static int constexpr partial_vector_elements() {
        return columns & row_vector_type::size() - 1;        
    }

protected:
    // storage matrix
    row_vector_type mat[rows][vectors_per_row()];

    // pointer to error handling function
    void (*errorfunction)(void);

    // tell if elements are a floating point type
    static constexpr bool is_fp() {
        return (V::elementtype() & 0x3F) >= 15;
    }
public:
    // set function pointer to error handler
    void set_error_handler(void (*err)(void)) {
        errorfunction = err;
    }
    // get error handler function
    void (*get_error_handler())(void) {
        return errorfunction;
    }

    // read row as one or more vectors
    row_vector_type get_row(int r, int i = 0) const {
        if (uint32_t(r) < rows && uint32_t(i) < uint32_t(vectors_per_row())) {
            return mat[r][i];
        }
        else {
            (*errorfunction)();                     // call error handler
            if constexpr(is_fp()) {
                return nan_vec<row_vector_type>(2); // floating point type. return NAN
            }
            else {
                return row_vector_type(etype(0));   // integer type. return 0
            }
        }
    }
    // write row as one or more vectors
    void set_row(row_vector_type x, int r, int i = 0) {
        if (uint32_t(r) < rows && uint32_t(i) < uint32_t(vectors_per_row())) {
            int limit = columns - i * row_vector_type::size();
            if (limit < row_vector_type::size()) x.cutoff(limit);  // set unused vector elements to zero
            mat[r][i] = x;
        }
        else {
            (*errorfunction)();                     // call error handler
        }
    }
    // get single array element
    etype get_element(int row, int column) const {
        constexpr uint32_t rvsize = row_vector_type::size();
        if (uint32_t(row) < uint32_t(rows) && uint32_t(column) < uint32_t(columns)) {
            return mat[row][uint32_t(column) / rvsize][column & (rvsize - 1)];
        }
        else {
            (*errorfunction)();                     // call error handler
            if constexpr(is_fp()) {
                return nan_vec<row_vector_type>(2)[0]; // floating point type. return NAN
            }
            else {
                return etype(0);   // integer type. return 0
            }        
        }
    }
    // change single array element
    void set_element(etype x, int row, int column) {
        constexpr uint32_t rvsize = row_vector_type::size();
        if (uint32_t(row) < uint32_t(rows) && uint32_t(column) < uint32_t(columns)) {        
            mat[row][uint32_t(column) / rvsize].insert(column & rvsize - 1, x);
        }
        else {
            (*errorfunction)();                     // call error handler
        }
    }
    // load entire matrix from an array or C-style matrix (row major order)
    void load(void * p) {
        int r, c;   // row and column index
        for (r = 0; r < rows; r++) {   // row loop     
            for (c = 0; c < full_vectors_per_row(); c++) {  // column loop               
                mat[r][c].load((etype*)p + r * columns + c * row_vector_type::size());  // load full vector
            }
            if constexpr (partial_vector_elements()) {
                mat[r][c].load_partial(partial_vector_elements(), (etype*)p + r * columns + c * row_vector_type::size()); // load remaining partial vector
            } 
        }    
    }
    // store entire matrix to an array or C-style matrix (row major order)
    void store(void * p) const {
        int constexpr partial = columns & row_vector_type::size() - 1;
        int r, c;   // row and column index
        for (r = 0; r < rows; r++) {   // row loop     
            for (c = 0; c < full_vectors_per_row(); c++) {  // column loop               
                mat[r][c].store((etype*)p + r * columns + c * row_vector_type::size());  // store full vector
            }
            if constexpr (partial_vector_elements()) {
                mat[r][c].store_partial(partial_vector_elements(), (etype*)p + r * columns + c * row_vector_type::size()); // store remaining partial vector
            } 
        }    
    }
    // set the matrix to all zeroes
    void zero() {
        int r, c;   // row and column index
        for (r = 0; r < rows; r++) {   // row loop     
            for (c = 0; c < vectors_per_row(); c++) {  // column loop               
                mat[r][c] = row_vector_type(etype(0));
            }
        }    
    }
};

// used internally: find a vector type big enough to contain n elements
template <typename V, int n>
auto constexpr match_vector_type() {
    if constexpr (V::size() >= n) return V();  // V is big enough
    else if constexpr (sizeof(V) < 64) {
        typedef decltype(extend_z(V())) double_size;  // double size of V
        if constexpr (double_size::size() >= n) return double_size();
        else if constexpr (sizeof(double_size) < 64) {
            auto constexpr quad_size = decltype(extend_z(double_size()))();
            if constexpr (quad_size.size() >= n) return quad_size;
            else return V();  // failure
        }
        else return V();  // failure
    }
    else return V();  // failure
}


// pack two consecutive rows of a MatrixV matrix into a single vector,
// provided that a vector type of sufficient size exists
template <typename M>
auto pack2rows(M & matrix, int first_row) {
    // number of columns
    int constexpr ncol = M::ncolumns();
    // number of elements in result
    int constexpr npack = ncol * 2;
    // vectortype for row
    typedef decltype(matrix.get_row(0, 0)) row_type;
    // vectortype for 2 rows
    typedef decltype(match_vector_type<row_type, npack>()) pack_type;
    // check if sufficiently big vector type exists
    static_assert(pack_type::size() >= npack, "no vector type is big enough");
    // check if row index out of range
    if (uint32_t(first_row + 1) >= uint32_t(matrix.nrows())) {
        (*(matrix.get_error_handler()))();                   // row index out of range. call error handler
        return pack_type(decltype(pack_type()[0])(0));     // return vector of zeroes
    }

    // read two rows
    auto r1 = matrix.get_row(first_row, 0);
    auto r2 = matrix.get_row(first_row+1, 0);

    // facilitate optimal permutation method. AVX512VL: use compress instruction with zeroing. otherwise: use an unused element as zero
    int constexpr u = INSTRSET >= 10 ? -1 : pack_type::size() - 1; 

    if constexpr (pack_type::size() == row_type::size()) {
        // same vector size for result
        if constexpr (row_type::size() == 4 && ncol == 2) {
            return blend4<0, 1, 4, 5>(r1, r2);
        }
        else if constexpr (row_type::size() == 8 && ncol == 2) return blend8<0,1,8,9,u,u,u,u>(r1, r2);        
        else if constexpr (row_type::size() == 8 && ncol == 3) return blend8<0,1,2,8,9,10,u,u>(r1, r2);        
        else if constexpr (row_type::size() == 8 && ncol == 4) return blend8<0,1,2,3,8,9,10,11>(r1, r2);
        else static_assert(ncol == 0, "unsupported vector combination"); // dummy false, template dependent
    }
    else {
        // result has double vector size
        auto r12 = concatenate2(r1, r2);  // concatenate rows

        if constexpr (row_type::size() == 2 && ncol == 2) return r12;
        else if constexpr (row_type::size() == 4 && ncol == 3) return permute8<0,1,2,4,5,6,u,u>(r12);
        else if constexpr (row_type::size() == 4 && ncol == 4) return r12;
        else if constexpr (row_type::size() == 8 && ncol == 5) return permute16<0,1,2,3,4,8,9,10,11,12,u,u,u,u,u,u>(r12);
        else if constexpr (row_type::size() == 8 && ncol == 6) return permute16<0,1,2,3,4,5,8,9,10,11,12,13,u,u,u,u>(r12);
        else if constexpr (row_type::size() == 8 && ncol == 7) return permute16<0,1,2,3,4,5,6,8,9,10,11,12,13,14,u,u>(r12);
        else if constexpr (row_type::size() == 8 && ncol == 8) return r12;
        else static_assert(ncol == 0, "unsupported vector combination"); // dummy false, template dependent
    }
}


// pack three consecutive rows of a MatrixV matrix into a single vector,
// provided that a vector type of sufficient size exists
template <typename M>
auto pack3rows(M & matrix, int first_row) {
    // number of columns
    int constexpr ncol = M::ncolumns();
    // number of elements in result
    int constexpr npack = ncol * 3;
    // vectortype for row
    typedef decltype(matrix.get_row(0, 0)) row_type;
    // vectortype for 3 rows
    typedef decltype(match_vector_type<row_type, npack>()) pack_type;
    // check if sufficiently big vector type exists
    static_assert(pack_type::size() >= npack, "no vector type is big enough");
    // check if row index out of range
    if (uint32_t(first_row + 2) >= uint32_t(matrix.nrows())) {
        (*(matrix.get_error_handler()))();                   // row index out of range. call error handler
        return pack_type(decltype(pack_type()[0])(0));       // return vector of zeroes
    }
    // read three rows
    auto r12 = pack2rows(matrix, first_row);         // row 1+2
    auto r3 = matrix.get_row(first_row+2);           // row 3
    decltype(r12) r3e;                               // make r3e same size as r12
    if constexpr(r12.size() > r3.size()) {
        r3e = extend_z(r3);                          // zero-extend r3 to r3e
    }
    else {
        r3e = r3;                                    // r3 already has same size as r12
    }
    // facilitate optimal permutation method. AVX512VL: use compress instruction with zeroing. otherwise: use an unused element as zero
    int constexpr u = INSTRSET >= 10 ? -1 : pack_type::size() - 1; 

    if constexpr (r12.size() == 8 && pack_type::size() == 8 && ncol == 2) {
        return blend8<0, 1, 2, 3, 8, 9, -1, -1>(r12, r3e);    
    }
    else if constexpr (r12.size() == 16 && pack_type::size() == 16) {
        if constexpr (ncol == 5) return blend16<0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 16, 17, 18, 19, 20, -1>(r12, r3e);    
        else static_assert(ncol == 0, "unsupported vector combination"); // dummy false, template dependent
    }
    else if constexpr (r12.size() < pack_type::size()) {
        // the result vector has double size of r12 and r3e       
        auto r123 = concatenate2(r12, r3e);  // concatenate rows

        if constexpr (r123.size() == 8 && ncol == 2) return r123;
        else if constexpr (r123.size() == 16 && ncol == 3) return permute16<0, 1, 2, 3, 4, 5, 8, 9, 10, u, u, u, u, u, u, u>(r123);
        else if constexpr (r123.size() == 16 && ncol == 4) return r123;
        else static_assert(ncol == 0, "unsupported vector combination"); // dummy false, template dependent
    }
    else static_assert(ncol == 0, "unsupported vector combination"); // dummy false, template dependent
}


// pack four consecutive rows of a MatrixV matrix into a single vector,
// provided that a vector type of sufficient size exists
template <typename M>
auto pack4rows(M & matrix, int first_row) {
    // number of columns
    int constexpr ncol = M::ncolumns();
    // number of elements in result
    int constexpr npack = ncol * 4;
    // vectortype for row
    typedef decltype(matrix.get_row(0, 0)) row_type;
    // vectortype for 4 rows
    typedef decltype(match_vector_type<row_type, npack>()) pack_type;
    // check if sufficiently big vector type exists
    static_assert(pack_type::size() >= npack, "no vector type is big enough");
    // check if row index out of range
    if (uint32_t(first_row + 3) >= uint32_t(matrix.nrows())) {
        (*(matrix.get_error_handler()))();                   // row index out of range. call error handler
        return pack_type(decltype(pack_type()[0])(0));       // return vector of zeroes
    }
    // read four rows
    auto r12 = pack2rows(matrix, first_row);         // row 1+2
    auto r34 = pack2rows(matrix, first_row+2);       // row 3+4
    if constexpr (r12.size() == 8 && ncol == 2) {
        return blend8<0, 1, 2, 3, 8, 9, 10, 11>(r12, r34);
    }
    else {  // row size > r12.size()
        auto r1234 = concatenate2(r12, r34);         // concatenate rows
        // facilitate optimal permutation method. AVX512VL: use compress instruction with zeroing. otherwise: use an unused element as zero
        int constexpr u = INSTRSET >= 10 ? -1 : pack_type::size() - 1;
        if constexpr (r1234.size() == 8 && ncol == 2) return r1234;
        else if constexpr (r1234.size() == 16 && ncol == 3) return permute16<0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 13, u, u, u, u>(r1234);
        else if constexpr (r1234.size() == 16 && ncol == 4) return r1234;
        else if constexpr (r1234.size() == 32 && ncol == 5) return permute32<0,1,2,3,4,5,6,7,8,9,16,17,18,19,20,21,22,23,24,25,u,u,u,u,u,u,u,u,u,u,u,u>(r1234);
        else static_assert(ncol == 0, "unsupported vector combination"); // dummy false, template dependent
    }
}


// pack five consecutive rows of a MatrixV matrix into a single vector,
// provided that a vector type of sufficient size exists
template <typename M>
auto pack5rows(M & matrix, int first_row) {
    // number of columns
    int constexpr ncol = M::ncolumns();
    // number of elements in result
    int constexpr npack = ncol * 5;
    // vectortype for row
    typedef decltype(matrix.get_row(0, 0)) row_type;
    // vectortype for 5 rows
    typedef decltype(match_vector_type<row_type, npack>()) pack_type;
    // check if sufficiently big vector type exists
    static_assert(pack_type::size() >= npack, "no vector type is big enough");
    // check if row index out of range
    if (uint32_t(first_row + 4) >= uint32_t(matrix.nrows())) {
        (*(matrix.get_error_handler()))();                   // row index out of range. call error handler
        return pack_type(decltype(pack_type()[0])(0));       // return vector of zeroes
    }
    // read five rows
    auto r123 = pack3rows(matrix, first_row);         // row 1+2+3
    auto r45  = pack2rows(matrix, first_row+3);       // row 4+5
    decltype(r123) r45e;                              // make r45e same size as r123
    if constexpr(r123.size() > r45.size()) {
        r45e = extend_z(r45);                         // zero-extend r3 to r3e
    }
    else {
        r45e = r45;                                   // r3 already has same size as r12
    }
    if constexpr (r123.size() == 8 && ncol == 2) {
        auto r12345 = concatenate2(r123, r45e);           // concatenate rows
        return permute16<0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 15, 15, 15, 15, 15, 15>(r12345);
    }
    else if constexpr (r123.size() == 16 && ncol == 3) return blend16<0, 1, 2, 3, 4, 5, 6, 7, 8, 16, 17, 18, 19, 20, 21, 22>(r123, r45e);

    else if constexpr (r123.size() == 16 && ncol == 4) {
        auto r12345 = concatenate2(r123, r45e);           // concatenate rows
        int constexpr d = -1;
        return permute32<0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d>(r12345);
    }
    else static_assert(ncol == 0, "unsupported vector combination"); // dummy false, template dependent
}


// unpack a vector into two consecutive rows of a MatrixV.
// The length of the vector rr must be sufficient to hold two rows of matrix mm
template <typename V, typename M>
void unpack2rows(V rr, M & matrix, int first_row) {
    // number of columns
    int constexpr ncol = M::ncolumns();
    // number of elements in input
    int constexpr npack = ncol * 2;
    // vectortype for row
    typedef decltype(matrix.get_row(0, 0)) row_type;
    // vectortype for input
    typedef decltype(rr) pack_type;
    // check if sufficiently big vector type exists
    static_assert(pack_type::size() >= npack, "input vector type is too small");
    static_assert(pack_type::size() < npack * 2, "input vector type is too large");
    static_assert(pack_type::elementtype() == row_type::elementtype(), "wrong vector type");

    // check if row index out of range
    if (uint32_t(first_row + 1) >= uint32_t(matrix.nrows())) {
        (*(matrix.get_error_handler()))();                   // row index out of range. call error handler
    }
    int constexpr d = V_DC;  // don't care in permutation

    row_type r1, r2;                                       // row vectors
    if constexpr (row_type::size() == pack_type::size()) { // same vector size
        r1 = rr;
        if constexpr (row_type::size() == 4) r2 = permute4<2, 3, d, d>(rr);
        else if constexpr (row_type::size() == 8) r2 = permute8<ncol, ncol+1, ncol+2, ncol+3, d, d, d, d>(rr);
        else static_assert(ncol == 0, "unsupported vector combination"); // dummy false, template dependent
    }
    else {
        // (row_type::size() < pack_type::size()
        r1 = rr.get_low();        // first row. extra elements ignored
        pack_type r2e;            // second row, in bigger vector to match rr
        if constexpr (pack_type::size() == 4) r2e = permute4<2, 3, d, d>(rr);
        else if constexpr (pack_type::size() == 8) r2e = permute8<ncol, ncol+1, ncol+2, ncol+3, d, d, d, d>(rr);
        else if constexpr (pack_type::size() == 16) r2e = permute16<ncol, ncol+1, ncol+2, ncol+3, ncol+4, ncol+5, ncol+6, ncol+7, d, d, d, d, d, d, d, d>(rr);
        else static_assert(ncol == 0, "unsupported vector combination"); // dummy false, template dependent
        r2 = r2e.get_low();       // second row, smaller vector
    }
    // insert rows
    matrix.set_row(r1, first_row);
    matrix.set_row(r2, first_row+1);
}


// unpack a vector into three consecutive rows of a MatrixV.
// The length of the vector rr must be sufficient to hold three rows of matrix mm
template <typename V, typename M>
void unpack3rows(V rr, M & matrix, int first_row) {
    // number of columns
    int constexpr ncol = M::ncolumns();
    // number of elements in input
    int constexpr npack = ncol * 3;
    // vectortype for row
    typedef decltype(matrix.get_row(0, 0)) row_type;
    // vectortype for input
    typedef decltype(rr) pack_type;
    // check if sufficiently big vector type exists
    static_assert(pack_type::size() >= npack, "input vector type is too small");
    static_assert(pack_type::size() < npack * 2, "input vector type is too large");
    static_assert(pack_type::elementtype() == row_type::elementtype(), "wrong vector type");

    // check if row index out of range
    if (uint32_t(first_row + 2) >= uint32_t(matrix.nrows())) {
        (*(matrix.get_error_handler()))();                   // row index out of range. call error handler
    }
    row_type r1, r2, r3;                                       // row vectors
    int constexpr d = V_DC;  // don't care in permutation
    if constexpr (row_type::size() * 4 == pack_type::size()) { // quarter vector size
        if constexpr (ncol == 2 || ncol == 4) {        // 2 or 4 columns
            auto r12 = rr.get_low();
            auto r3_ = rr.get_high();
            r1 = r12.get_low();
            r2 = r12.get_high();
            r3 = r3_.get_low();
        }
        else {
            static_assert (ncol == 3, "unsupported vector combination");     // 3 columns
            r1 = rr.get_low().get_low();
            auto r23 = permute16<3, 4, 5, d, 6, 7, 8, d, d, d, d, d, d, d, d, d>(rr);
            auto r23e = r23.get_low();
            r2 = r23e.get_low();
            r3 = r23e.get_high();
        }
    }
    else if constexpr (row_type::size() * 2 == pack_type::size()) { // half vector size
        if constexpr (ncol == 2 && row_type::size() == 4) {        // 2 columns
            r1 = rr.get_low();
            r2 = permute4<2, 3, d, d>(r1);
            r3 = rr.get_high();
        }
        else if constexpr (ncol == 3 && row_type::size() == 8) {      // 3 columns
            r1 = rr.get_low();
            auto r23 = permute16<3, 4, 5, d, d, d, d, d, 6, 7, 8, d, d, d, d, d>(rr);
            r2 = r23.get_low();
            r3 = r23.get_high();
        }
        else if constexpr (ncol == 4 && row_type::size() == 8) {      // 3 columns
            r1 = rr.get_low();
            auto r23 = permute16<4, 5, 6, 7, d, d, d, d, 8, 6, 10, 11, d, d, d, d>(rr);
            r2 = r23.get_low();
            r3 = r23.get_high();
        }
        else if constexpr (ncol == 5 && row_type::size() == 8) {      // 5 columns
            r1 = rr.get_low();
            auto r23 = permute16<5, 6, 7, 8, 9, d, d, d, 10, 11, 12, 13, 14, d, d, d>(rr);
            r2 = r23.get_low();
            r3 = r23.get_high();
        }
        else static_assert(ncol == 0, "unsupported vector combination"); // dummy false, template dependent
    }
    else { // same vector size
        static_assert(row_type::size() == 8 && ncol == 2, "unsupported vector combination");
        r1 = rr;
        r2 = permute8<2, 3, d, d, d, d, d, d>(rr);
        r3 = permute8<4, 5, d, d, d, d, d, d>(rr);     
    }
    // insert rows
    matrix.set_row(r1, first_row);
    matrix.set_row(r2, first_row+1);
    matrix.set_row(r3, first_row+2);
}


// unpack a vector into four consecutive rows of a MatrixV.
// The length of the vector rr must be sufficient to hold four rows of matrix mm
template <typename V, typename M>
void unpack4rows(V rr, M & matrix, int first_row) {
    // number of columns
    int constexpr ncol = M::ncolumns();
    // number of elements in input
    int constexpr npack = ncol * 4;
    // vectortype for row
    typedef decltype(matrix.get_row(0, 0)) row_type;
    // vectortype for input
    typedef decltype(rr) pack_type;
    // check if sufficiently big vector type exists
    static_assert(pack_type::size() >= npack, "input vector type is too small");
    static_assert(pack_type::size() < npack * 2, "input vector type is too large");
    static_assert(pack_type::elementtype() == row_type::elementtype(), "wrong vector type");

    // check if row index out of range
    if (uint32_t(first_row + 3) >= uint32_t(matrix.nrows())) {
        (*(matrix.get_error_handler()))();                   // row index out of range. call error handler
    }
    row_type r1, r2, r3, r4;                                       // row vectors
    int constexpr d = V_DC;  // don't care in permutation

    if constexpr (row_type::size() * 4 == pack_type::size()) { // quarter vector size
        pack_type ru;
        if constexpr (ncol == 3 && row_type::size() == 4) {   // 3 columns
            ru = permute16<0, 1, 2, d, 3, 4, 5, d, 6, 7, 8, d, 9, 10, 11, d>(rr); // unpack into quarters
        }
        else if constexpr (ncol == 5 && row_type::size() == 8) {   // 5 columns
            ru = permute32<0, 1, 2, 3, 4, d, d, d, 5, 6, 7, 8, 9, d, d, d, 10, 11, 12, 13, 14, d, d, d, 15, 16, 17, 18, 19, d, d, d>(rr); // unpack into quarters
        }
        else {
            ru = rr;
            static_assert((row_type::size() == 2 && ncol == 2) || (row_type::size() == 4 && ncol == 4), "unsupported vector combination");
        }
        auto r12 = ru.get_low();
        auto r34 = ru.get_high();
        r1 = r12.get_low();
        r2 = r12.get_high();
        r3 = r34.get_low();
        r4 = r34.get_high();
    }
    else if constexpr(row_type::size() == 4 && ncol == 2) {    // half vector size
        //static_assert(, "unsupported vector combination");
        r1 = rr.get_low();
        r2 = permute4<2, 3, d, d>(r1);
        r3 = rr.get_high();
        r4 = permute4<2, 3, d, d>(r3);
    }
    else if constexpr (row_type::size() == 8 && ncol == 2) {
        r1 = rr;
        r2 = permute8<2, 3, d, d, d, d, d, d>(rr);
        r3 = permute8<4, 5, d, d, d, d, d, d>(rr);
        r4 = permute8<6, 7, d, d, d, d, d, d>(rr);
    }
    else if constexpr (row_type::size() == 8 && ncol == 3) {
        r1 = rr.get_low();
        r2 = permute8<3, 4, 5, d, d, d, d, d>(r1);
        auto r34 = permute16<6, 7, 8, d, d, d, d, d, 9, 10, 11, d, d, d, d, d>(rr);
        r3 = r34.get_low();
        r4 = r34.get_high();
    }
    else if constexpr (row_type::size() == 8 && ncol == 4) {
        r1 = rr.get_low();
        r2 = permute8<4, 5, 6, 7, d, d, d, d>(r1);
        r3 = rr.get_high();
        r4 = permute8<4, 5, 6, 7, d, d, d, d>(r3);
    }
    else static_assert(ncol == 0, "unsupported vector combination"); // dummy false, template dependent
    // insert rows
    matrix.set_row(r1, first_row);
    matrix.set_row(r2, first_row+1);
    matrix.set_row(r3, first_row+2);
    matrix.set_row(r4, first_row+3);
}


// unpack a vector into five consecutive rows of a MatrixV.
// The length of the vector rr must be sufficient to hold five rows of matrix mm
template <typename V, typename M>
void unpack5rows(V rr, M & matrix, int first_row) {
    // number of columns
    int constexpr ncol = M::ncolumns();
    // number of elements in input
    int constexpr npack = ncol * 5;
    // vectortype for row
    typedef decltype(matrix.get_row(0, 0)) row_type;
    // vectortype for input
    typedef decltype(rr) pack_type;
    // check if sufficiently big vector type exists
    static_assert(pack_type::size() >= npack, "input vector type is too small");
    static_assert(pack_type::elementtype() == row_type::elementtype(), "wrong vector type");

    // check if row index out of range
    if (uint32_t(first_row + 4) >= uint32_t(matrix.nrows())) {
        (*(matrix.get_error_handler()))();                   // row index out of range. call error handler
    }
    row_type r1, r2, r3, r4, r5;                                       // row vectors
    int constexpr d = V_DC;  // don't care in permutation

    if constexpr (row_type::size() * 4 == pack_type::size()) {  // quarter vector size
        if constexpr (ncol == 2 && row_type::size() == 4) {
            auto r1234 = rr.get_low();
            r1 = r1234.get_low();
            r2 = permute4<2, 3, d, d>(r1);
            r3 = r1234.get_high();
            r4 = permute4<2, 3, d, d>(r3);
            r5 = rr.get_high().get_low();
        }
        else if constexpr (ncol == 3 && row_type::size() == 4) {
            auto ru = permute16<0, 1, 2, d, 3, 4, 5, d, 6, 7, 8, d, 9, 10, 11, d>(rr); // unpack first 4 rows into quarters
            auto r12 = ru.get_low();
            auto r34 = ru.get_high();
            r1 = r12.get_low();
            r2 = r12.get_high();
            r3 = r34.get_low();
            r4 = r34.get_high();
            r5 = rr.get_high().get_high();
        }
        else if constexpr (ncol == 4 && row_type::size() == 8) {
            auto r1243 = rr.get_low();
            r1 = r1243.get_low();
            r2 = permute8<4, 5, 6, 7, d, d, d, d>(r1);
            r3 = r1243.get_high();
            r4 = permute8<4, 5, 6, 7, d, d, d, d>(r3);
            r5 = rr.get_high().get_low();
        }
        else {
            static_assert(ncol == 0, "unsupported vector combination"); // dummy false, template dependent
        }
    }
    else {    
        static_assert(row_type::size() == 8, "unsupported vector combination");
        if constexpr (ncol == 2) {
            r1 = rr.get_low();
            r2 = permute8<2, 3, d, d, d, d, d, d>(r1);
            r3 = permute8<4, 5, d, d, d, d, d, d>(r1);
            r4 = permute8<6, 7, d, d, d, d, d, d>(r1);
            r5 = rr.get_high();
        } 
        else if constexpr (ncol == 3) {
            r1 = rr.get_low();
            r2 = permute8<3, 4, 5, d, d, d, d, d>(r1);
            auto r345 = permute16<6, 7, 8, d, 9, 10, 11, d, 12, 13, 14, d, d, d, d, d>(rr);            
            r3 = r345.get_low();
            r4 = permute8<4, 5, 6, d, d, d, d, d>(r3);
            r5 = r345.get_high();
        }
        else {
            static_assert(ncol == 0, "unsupported vector combination"); // dummy false, template dependent
        }
    }
    // insert rows
    matrix.set_row(first_row,   r1);
    matrix.set_row(first_row+1, r2);
    matrix.set_row(first_row+2, r3);
    matrix.set_row(first_row+3, r4);
    matrix.set_row(first_row+4, r5);
}

#ifdef VCL_NAMESPACE
}
#endif

#endif // MATRIXV_H
