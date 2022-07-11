/************************  vector_containers.h   ******************************
* Author:        Agner Fog
* Date created:  2022-07-04
* Last modified: 2022-07-11
* Version:       2.02.00
* Project:       vector class library
* Description:
* Header file for container classes
* These containers can contain vector class objects and matrixes
*
* For instructions, see containers_manual.pdf
*
* (c) Copyright 2022 Agner Fog.
* Apache License version 2.0 or later.
******************************************************************************/

#ifndef VECTOR_CONTAINERS_H
#define VECTOR_CONTAINERS_H 20200

#ifdef VCL_NAMESPACE
namespace VCL_NAMESPACE {
#endif


// Container class to store n vector class objects of type V
template <typename V, int n>
class ContainerV {
protected:
    V buf[n];                                    // array of vectors
    int s_count() const {                        // used internally
        constexpr int s = V::size();             // vector size
        static_assert((s & s-1) == 0, "vector size must be power of 2"); // check that vector size is a power of 2
        return bit_scan_reverse_const(s);        // shift count for fast division by vector size
    }
    void (*errorfunction)(void) = 0;             // pointer to error handling function
public:
    ContainerV() = default;                      // default constructor
    void set_error_handler(void (*e)(void)) {    // set function pointer to error handler
        errorfunction = e;
    }
    typedef decltype (buf[0][0]) etype;          // type of vector elements
    static constexpr int n_vectors() {           // get number of vectors
        return n;
    }
    static constexpr int n_elements() {          // get number of vector elements
        return n * V::size();
    }
    static constexpr int elementtype() {         // info about vector element type and container type
        return  V::elementtype() | 0x1000;
    }
    static constexpr bool is_fp() {              // true if elements are a floating point type
        return (V::elementtype() & 0x3F) >= 15;
    }
    V get_vector(int index) const {              // extract one vector
        if (uint32_t(index) < n) {
            return buf[index];                   // get vector
        }
        else {                                   // index out of range
            (*errorfunction)();                  // call error handler
            if constexpr(is_fp()) {
                return nan_vec<V>(2);            // floating point type. return NAN
            }
            else {
                return V(etype(0));              // integer type. return 0
            }
        }
    }
    void set_vector(V x, int index) {            // insert one vector
        if (uint32_t(index) < n) {
            buf[index] = x;                      // set vector
        }
        else {                                   // error
            (*errorfunction)();                  // call error handler
        }
    }
    etype get_element(uint32_t index) const {    // extract one vector element
        if (index < (uint32_t)n_elements()) {
            return buf[index >> s_count()][index & (V::size() - 1)];        
        }
        else {                                   // index out of range
            (*errorfunction)();                  // call error handler
            if constexpr(is_fp()) {
                return nan_vec<V>(2)[0];         // floating point type. return NAN
            }
            else {
                return 0;                        // integer type. return 0
            }
        }
    }
    void set_element(etype x, uint32_t index) {  // insert one vector element
        if (index < (uint32_t)n_elements()) {
            buf[index >> s_count()].insert(index & (V::size()-1), x);
        }
        else {                                   // error
            (*errorfunction)();                  // call error handler
        }
    }
    void load(int nn, void * p) {                // load nn elements from array
        if (nn <= 0) return;                     // nothing to do
        if (nn > n_elements()) nn = n_elements();// max size
        int m = (uint32_t)nn >> s_count();       // number of full vectors to load
        int i;                                   // loop counter
        for (i = 0; i < m; i++) {
            buf[i].load((etype*)p + i * V::size()); // store one vector
        }
        int partial = nn & (V::size() - 1);      // any partial store needed
        if (partial) {                           // nn is not divisible by vector size
            // load partial vector in the end
            buf[i].load_partial(partial, (etype*)p + i * V::size()); // load part of last vector       
        }
    }
    void store(int nn, void * p) {               // store nn elements to array
        if (nn <= 0) return;                     // nothing to do
        if (nn > n_elements()) nn = n_elements();// max size
        int m = (uint32_t)nn >> s_count();       // number of full vectors to store
        int i;                                   // loop counter
        for (i = 0; i < m; i++) {
            buf[i].store((etype*)p + i * V::size()); // store one vector
        }
        int partial = nn & (V::size() - 1);      // any partial store needed
        if (partial) {                           // nn is not divisible by vector size
            // store partial vector in the end
            buf[i].store_partial(partial, (etype*)p + i * V::size()); // store part of last vector       
        }
    }
    V * get_buf() {                              // get address of internal buffer
        return buf;
    }
    void zero() {                                // set all contents to zero
        for (int i = 0; i < n; i++) {
            buf[i] = V(etype(0));
        }
    }
};



// Container class to store a variable number of vector class objects of type V
template <typename V>
class ContainerV <V, 0> {
protected:
    V * buf;                                     // allocated memory buffer containing array of vectors
    uint32_t allocatedSize;                      // size of allocated buffer
    uint32_t nvectors;                           // number of vectors currently used (includes partially used)
    uint32_t nelements;                          // number of vector elements currently used
    void (*errorfunction)(void);                 // pointer to error handling function
    int s_count() const {                        // used internally
        constexpr int s = V::size();             // vector size
        static_assert((s & s-1) == 0, "vector size must be a power of 2"); // check that vector size is a power of 2
        return bit_scan_reverse_const(s);        // shift count for fast division by vector size
    }
public:
    ContainerV() {                               // constructor
        buf = 0;  allocatedSize = 0;  nvectors = 0;  nelements = 0;  errorfunction = 0;
    }
    ~ContainerV() {                              // destructor
        if (buf) delete[] buf;                   // free allocated memory
    }
    ContainerV(ContainerV&) = delete;            // prevent copying entire container (a copy constructor would have to allocate a new buffer)
    ContainerV operator = (ContainerV&) = delete;// prevent copying entire container
    void set_error_handler(void (*e)(void)) {    // set function pointer to error handler
        errorfunction = e;
    }
    typedef decltype (buf[0][0]) etype;          // type of vector elements
    static constexpr int elementtype() {         // info about vector element type and container type
        return V::elementtype() | 0x1000;
    }
    static constexpr bool is_fp() {              // true if elements are a floating point type
        return (V::elementtype() & 0x3F) >= 15;
    }
    int n_vectors() const {                      // get number of vectors
        return nvectors;
    }
    int n_elements() const {                     // get number of vector elements
        return nelements;
    }
    int allocated_size() const {                 // maximum size that can be set without reallocation
        return allocatedSize;
    }
    void set_nvectors(int size) {
        // Allocate, reallocate or deallocate buffer of specified size. size is the number of vectors.
        // Setting size > allocated_size will allocate more buffer and fill it with zeroes
        // Setting size < allocated_size will decrease size so that some of the data are inaccessible
        // Setting size = 0 will discard all data and de-allocate the buffer.
        if (size <= 0) {                         // discard everything         
            if (buf) delete[] buf;               // de-allocate buffer
            buf = 0;  allocatedSize = 0;  nvectors = 0;  nelements = 0;
        }
        else if (uint32_t(size) <= allocatedSize) { // grow or shrink within allocated size
            nvectors = size;  nelements = size * V::size();
        }
        else {                                   // increase allocated size
            uint32_t newallocsize;               // new size to allocate
            if (uint32_t(size) >= allocatedSize + allocatedSize/2) {
                newallocsize = size;             // first time or big increase. allocate only the specified size
            }
            else {
                newallocsize = size*2;           // small increase. allocate more than requested to avoid frequent reallocations
            }
            V * buf2 = 0;                        // pointer to new buffer
            buf2 = new V[newallocsize]();        // allocate new buffer. () means initialize to zero
            uint32_t i = 0;                      // loop counter
            if (buf) {                           // previously allocated buffer exists
                for (i = 0; i < allocatedSize; i++) {
                    buf2[i] = buf[i];            // copy from old to new buffer
                }
                delete [] buf;                   // deallocate old buffer         
            }
            // store pointer to new buffer
            buf = buf2;  allocatedSize = newallocsize;
            nvectors = size;  nelements = size * V::size(); // new used size        
        }
    }
    void set_nelements(int n) {
        // Allocate, reallocate or deallocate buffer of specified size, not necessarily a multiple of the vector size
        int nv = uint32_t(n + V::size() - 1) >> s_count(); // round up to nearest multiple of the vector size
        set_nvectors(nv);
        nelements = n;
    }
    V get_vector(int index) const {              // extract one vector
        if (uint32_t(index) < nvectors) {
            return buf[index];                   // get vector
        }
        else {                                   // index out of range
            (*errorfunction)();                  // call error handler
            if constexpr(is_fp()) {
                return nan_vec<V>(2);            // floating point type. return NAN
            }
            else {
                return V(etype(0));              // integer type. return 0
            }
        }
    }
    void set_vector(V x, int index) {            // insert one vector
        if (uint32_t(index) < nvectors) {
            buf[index] = x;                      // set vector
        }
        else {                                   // error
            (*errorfunction)();                  // call error handler
        }
    }
    etype get_element(uint32_t index) const {    // extract one vector element
        if (index <  uint32_t(nelements)) {
            return buf[index >> s_count()][index & (V::size() - 1)];        
        }
        else {                                   // index out of range
            (*errorfunction)();                  // call error handler
            if constexpr(is_fp()) {
                return nan_vec<V>(2)[0];         // floating point type. return NAN
            }
            else {
                return 0;                        // integer type. return 0
            }
        }
    }
    void set_element(etype x, uint32_t index) {  // insert one vector element
        if (index <  uint32_t(nelements)) {
            buf[index >> s_count()].insert(index & (V::size()-1), x);
        }
        else {                                   // error
            (*errorfunction)();                  // call error handler
        }
    }
    void load(int n, void * p) {                 // load n elements from array
        if (n <= 0) return;                      // nothing to do
        if (uint32_t(n) > nelements) n = nelements;// max size
        int m = (uint32_t)n >> s_count();        // number of full vectors to load
        int i;                                   // loop counter
        for (i = 0; i < m; i++) {
            buf[i].load((etype*)p + i * V::size()); // store one vector
        }
        int partial = n & (V::size() - 1);       // any partial store needed
        if (partial) {                           // n is not divisible by vector size
            // load partial vector in the end
            buf[i].load_partial(partial, (etype*)p + i * V::size()); // load part of last vector       
        }
    }
    void store(int n, void * p) {                // store n elements to array
        if (n <= 0) return;                      // nothing to do
        if (uint32_t(n) > nelements) n = nelements;// max size
        int m = (uint32_t)n >> s_count();        // number of full vectors to store
        int i;                                   // loop counter
        for (i = 0; i < m; i++) {
            buf[i].store((etype*)p + i * V::size()); // store one vector
        }
        int partial = n & (V::size() - 1);       // any partial store needed
        if (partial) {                           // n is not divisible by vector size
            // store partial vector in the end
            buf[i].store_partial(partial, (etype*)p + i * V::size()); // store part of last vector       
        }
    }
    V * get_buf() {                              // get address of internal buffer. warning: address may change
        return buf;
    }
    void zero() {                                // set all contents to zero
        for (uint32_t i = 0; i < nvectors; i++) {
            buf[i] = V(etype(0));
        }
    }
};

#ifdef VCL_NAMESPACE
}
#endif

#endif // VECTOR_CONTAINERS_H
