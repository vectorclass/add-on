/************************  general_containers.h   *****************************
* Author:        Agner Fog
* Date created:  2022-07-05
* Last modified: 2022-07-11
* Version:       2.02.00
* Description:
* Header file for general container classes
* These containers are independent of the vector class library and intended 
* for objects that are not VCL vectors.
* It may not be suitable for objects that have non-standard constructors, 
* copy constructors, move constructors, or destructors.
*
* Example:

  ContainerG<double> c;      // make container for type double
  c.set_size(10);            // allocate space for 10 objects
  c[2] = 88.8;               // change value of one object
  // print out all objects
  for (int i = 0; i < c.size(); i++) printf(" %.1f", c[i]); 
 
* For further instructions, see containers_manual.pdf
*
* (c) Copyright 2022 Agner Fog.
* Apache License version 2.0 or later.
******************************************************************************/

#ifndef GENERAL_CONTAINERS_H
#define GENERAL_CONTAINERS_H 20200

// Container class to store a variable number of objects of any type.
// This container does not rely on the vector class library
template <typename T>
class ContainerG {
protected:
    T * buf;                                     // allocated memory buffer containing array
    unsigned int allocatedSize;                  // size of allocated buffer
    unsigned int nobjects;                       // number of objects currently used
    void (*errorfunction)(void);                 // pointer to error handling function
public:
    ContainerG() {                               // constructor
        buf = 0;  allocatedSize = 0;  nobjects = 0;  errorfunction = 0;
    }
    ~ContainerG() {                              // destructor
        if (buf) delete[] buf;                   // free allocated memory
    }
    ContainerG(ContainerG&) = delete;            // prevent copying entire container (a copy constructor would have to allocate a new buffer)
    ContainerG operator = (ContainerG&) = delete;// prevent copying entire container
    int size() const {                           // get size as number of objects
        return nobjects;
    }
    int allocated_size() const {                 // maximum size that can be set without reallocation
        return allocatedSize;
    }
    void set_error_handler(void (*e)(void)) {    // set function pointer to error handler
        errorfunction = e;
    }
    void set_size(int size) {
        // Allocate, reallocate or deallocate buffer of specified size. size is the number of objects.
        // Setting size > allocated_size will allocate more buffer and fill it with zeroes
        // Setting size < allocated_size will decrease size so that some of the data are inaccessible
        // Setting size = 0 will discard all data and de-allocate the buffer.
        if (size <= 0) {                         // discard everything         
            if (buf) delete[] buf;               // de-allocate buffer
            buf = 0;  allocatedSize = 0;  nobjects = 0;
        }
        else if ((unsigned int)size <= allocatedSize) { // grow or shrink within allocated size
            nobjects = size;
        }
        else {                                   // increase allocated size
            unsigned int newallocsize;           // new size to allocate
            if ((unsigned int)size >= allocatedSize + allocatedSize/2) {
                newallocsize = size;             // first time or big increase. allocate only the specified size
            }
            else {
                newallocsize = size*2;           // small increase. allocate more than requested to avoid frequent reallocations
            }
            T * buf2 = 0;                        // pointer to new buffer
            buf2 = new T[newallocsize]();        // allocate new buffer. () means initialize to zero
            if (buf) {                           // previously allocated buffer exists
                for (unsigned int i = 0; i < allocatedSize; i++) {
                    buf2[i] = buf[i];            // copy from old to new buffer
                }
                delete [] buf;                   // deallocate old buffer         
            }
            // store pointer to new buffer
            buf = buf2;  allocatedSize = newallocsize;
            nobjects = size;                     // new used size        
        }
    }
    T & operator [] (int index)  {               // access one object
        if ((unsigned int)index < nobjects) {
            return buf[index];                   // get reference to object
        }
        else {                                   // index out of range
            (*errorfunction)();                  // report error
            return buf[0];
        }
    }
    void load(int n, T * p) {                    // load n objects from array
        if (n <= 0) return;                      // nothing to do
        if ((unsigned int)n > nobjects) n = nobjects;// max size
        for (int i = 0; i < n; i++) {
            buf[i] = p[i];                       // load n objects
        }
    }
    void store(int n, T * p) {                   // store n elements to array
        if (n <= 0) return;                      // nothing to do
        if (uint32_t(n) > nobjects) n = nobjects;// max size
        for (int i = 0; i < n; i++) {
            p[i] = buf[i];                       // store n objects
        }
    }
    T * get_buf() {                              // get address of internal buffer. warning: address may change
        return buf;
    }
};


#endif // GENERAL_CONTAINERS_H
