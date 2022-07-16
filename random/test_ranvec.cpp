/*************************  test_ranvec.cpp   *********************************
* Author:        Agner Fog
* Date created:  2019-07-08
* Last modified: 2022-07-16
* Version:       2.02
* Project:       add-on package for vector class library
* Description:
* Test program for ranvec1.cpp
*
******************************************************************************/


#include <stdio.h>


#ifndef INSTRSET
#define INSTRSET 10                    // instruction set
#endif

#define MAX_VECTOR_SIZE 512

#include "vectorclass.h"               // vector class library
#include "ranvec1.cpp"                 // random number generator
#include "physseed.cpp"


int main() {
    // Make instance of random number generator class, type 3.
    Ranvec1 ran(3);
    //Ranvec1 ran(3, 0);   // constructor with seed

#if true           // initialize with single seed
    ran.init(0);
#else              // initialize with array of seeds
    const int numseeds = 5;
    const int seeds[numseeds] = {5,4,3,2,1};
    ran.initByArray(seeds, numseeds);
#endif

    Vec16i ri = ran.random16i(0, 99);            // random integers in interval 0 - 99
    Vec16f rf = ran.random16f();                 // random floats in interval 0 - 1

    for (int i=0; i<ri.size(); i++) {            // print random integers
        printf("%3i  ", ri[i]);
    }
    printf("\n\n");

    for (int i=0; i<rf.size(); i++) {            // print random floats
        printf("%8.4f  ", rf[i]);
    }
    printf("\n\n");    

    for (int i = 0; i < 1000; i++) {             // call 1000 times
        ri += ran.random512b();
    }

    printf("%X\n", ri[7]);                       // print the sum

    // test physical seed generator
    printf("\nSeed type %i\n", physicalSeedType());

    // print two physical seeds. Must be independent if seed type > 1
    printf("\nSeed = %08X %08X\n", physicalSeed(), physicalSeed());

    return 0;
}
