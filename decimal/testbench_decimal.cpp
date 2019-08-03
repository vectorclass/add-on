/*************************  testbench_decimal.cpp   **************************
* Author:        Agner Fog
* Date created:  2019-07-14
* Last modified: 2019-07-14
* Version:       2.00
* Project:       Vector class library add-on package 'decimal'
* Description:   Testbench for decimal.cpp using vector class library
* Compile and run this program to test functions in decimal.h package
*
* Instructions:
* The following parameters must be defined on the command line or added in the
* top of this file:
*
* testcase: A number defining a function or operator to test. See the cases in this file.
*
* Compile with any compiler supported by VCL.
* Specify the desired instruction set and optimization options as parameters
* to the compiler.
*
* (c) Copyright 2019 Agner Fog.
* Apache license 2.0
******************************************************************************

Test cases:
1:  bin2ascii
2:  ascii2bin

*****************************************************************************/

#include <stdio.h>
#include <cmath>

#include <string.h>

#define MAX_VECTOR_SIZE 512
#ifndef INSTRSET
#define INSTRSET  8
#endif

//#define __AVX512VBMI2__

#include <vectorclass.h>
#include "decimal.cpp"


// ---------------------------------------------------------------------------
//            Specify input parameters here if running from an IDE
// ----------------------------------------------------------------------------

#ifndef testcase

#define testcase 1

#endif  // testcase 

// ----------------------------------------------------------------------------
//             Declarations
// ----------------------------------------------------------------------------
int globalError = 0;       // any error indicated in program return


/************************************************************************
*
*                          Test cases
*
************************************************************************/

#if testcase == 1   // test bin2ascii

// check results of ascii2bin
void checkb2a (int len, const char * res, const char * expected) {
    int slen = (int)strlen(res);
    if (strcmp(res, expected) != 0) {
        printf("\nbin2ascii error. Result:\n  >%s<\nExpected:\n  >%s<",
            res, expected);    
        globalError++;
    }
    else if (len != slen) {
        printf("\nbin2ascii length error. Actual length: %i, Reported length: %i\n  (%s)",
            slen, len, res);
        globalError++;
    }
}


int main() {
    char text[1024];
    int r = 0;

#if 0  // debugging
    Vec4i a0 (-87654321,-200000,3000000,40000000);
    r = bin2ascii(a0, text, 10, 4, '*', ',', true, true);
    printf("\nr=%i, \ntext=%s", r, text);

#else

    //static int bin2ascii (
    //   Vec4i const & a, char * string, int fieldlen = 8, int numdat = 4, bool signd = true, char ovfl = '*', char separator = ',', bool term = true) {
    Vec4i a1 (101,- 202,30303,-4040404);
    r = bin2ascii(a1, text, 10, 4, '*', ',', true, true);
    checkb2a(43, text, "       101,      -202,     30303,  -4040404");
    r = bin2ascii(a1, text, 5, 4, 0, ';', true, true);
    checkb2a(26, text, "  101; -202;30303;-4040404");
    r = bin2ascii(a1, text, 5, 4, '*', ',', true, true);
    checkb2a(23, text, "  101, -202,30303,*****");

    Vec4i a2 (101,-20202,-30303030,404040404);
    r = bin2ascii(a2, text, 10, 4, '*', ',', true, true);
    checkb2a(43, text, "       101,    -20202, -30303030, 404040404");
    r = bin2ascii(a2, text, 9, 4, '*', ',', true, true);
    checkb2a(39, text, "      101,   -20202,-30303030,404040404");
    r = bin2ascii(a2, text, 8, 4, '*', ',', true, true);
    checkb2a(35, text, "     101,  -20202,********,********");
    r = bin2ascii(a2, text, 6, 4, '*', ',', true, true);
    checkb2a(27, text, "   101,-20202,******,******");

    Vec4i a3 (-1,-100,10000,-10000);
    r = bin2ascii(a3, text, 6, 4, '*', ',', true, true);
    checkb2a(27, text, "    -1,  -100, 10000,-10000");
    r = bin2ascii(a3, text, 2, 4, '*', ',', true, true);
    checkb2a(11, text, "-1,**,**,**");
    r = bin2ascii(a3, text, 1, 4, '*', ',', true, true);
    checkb2a(7, text, "*,*,*,*");
    r = bin2ascii(a3, text, 0, 4, '*', ',', true, true);
    checkb2a(0, text, "");
    r = bin2ascii(a3, text, 5, 3, '*', ',', true, true);
    checkb2a(17, text, "   -1, -100,10000");
    r = bin2ascii(a3, text, 5, 3, '*', 0, true, true);
    checkb2a(15, text, "   -1 -10010000");
    r = bin2ascii(a3, text, 1, 3, 0, ',', true, true);
    checkb2a(13, text, "-1,-100,10000");
    r = bin2ascii(a3, text, 1, 2, 0, ',', true, true);
    checkb2a(7, text, "-1,-100");
    
    Vec4i a4 (-100000,-1000000,10000000,-100000000);
    r = bin2ascii(a4, text, 6, 4, 0, ',', true, true);
    checkb2a(36, text, "-100000,-1000000,10000000,-100000000");
    r = bin2ascii(a4, text, 7, 4, '*', ',', true, true);
    checkb2a(31, text, "-100000,*******,*******,*******");
    r = bin2ascii(a4, text, 7, 4, '*', ',', false, true);
    checkb2a(31, text, "*******,*******,*******,*******");
    r = bin2ascii(a4, text, 7, 4, 0, ',', false, true);
    checkb2a(41, text, "4294867296,4293967296,10000000,4194967296");
    
    Vec4i a5 (10000000,1000000000,2000000000,3000000000u);
    r = bin2ascii(a5, text, 8, 4, 0, ',', true, true);
    checkb2a(42, text, "10000000,1000000000,2000000000,-1294967296");
    r = bin2ascii(a5, text, 8, 4, 0, ',', false, true);
    checkb2a(41, text, "10000000,1000000000,2000000000,3000000000");
    
    Vec4i a6 (1,2,3,4);
    r = bin2ascii(a6, text, 2, 4, '*', ',', true, false); // no terminator. the rest of the previous result remains
    checkb2a(41, text, " 1, 2, 3, 400000000,2000000000,3000000000");
    r = bin2ascii(a6, text, 2, 4, 0, ',', true, true);
    checkb2a(11, text, " 1, 2, 3, 4");
    r = bin2ascii(a6, text, 2, 4, 0, 0, true, true);
    checkb2a(8, text, " 1 2 3 4");
    r = bin2ascii(a6, text, 1, 4, 0, 0, true, true);
    checkb2a(4, text, "1234");


//static int bin2ascii (
//  Vec8i const & a, char * string, int fieldlen = 8, int numdat = 8, bool signd = true, char ovfl = '*', char separator = ',', bool term = true) {
    Vec8i b1 (1,-22,333,-4321,55555,-666,7,8000);
    r = bin2ascii(b1, text, 10, 8, '*', ',', true, true);
    checkb2a(87, text, "         1,       -22,       333,     -4321,     55555,      -666,         7,      8000");
    r = bin2ascii(b1, text, 5, 8, 0, ';', true, true);
    checkb2a(47, text, "    1;  -22;  333;-4321;55555; -666;    7; 8000");
    r = bin2ascii(b1, text, 4, 8, '*', '|', true, true);
    checkb2a(39, text, "   1| -22| 333|****|****|-666|   7|8000");
    r = bin2ascii(b1, text, 6, 7, '*', ',', true, true);
    checkb2a(48, text, "     1,   -22,   333, -4321, 55555,  -666,     7");
    r = bin2ascii(b1, text, 6, 6, '*', ',', true, true);
    checkb2a(41, text, "     1,   -22,   333, -4321, 55555,  -666");
    r = bin2ascii(b1, text, 6, 5, '*', ',', true, true);
    checkb2a(34, text, "     1,   -22,   333, -4321, 55555");    
    r = bin2ascii(b1, text, 6, 4, '*', ',', true, true);
    checkb2a(27, text, "     1,   -22,   333, -4321");    
    r = bin2ascii(b1, text, 6, 3, '*', ',', true, true);
    checkb2a(20, text, "     1,   -22,   333");
    r = bin2ascii(b1, text, 6, 2, '*', ',', true, true);
    checkb2a(13, text, "     1,   -22");
    r = bin2ascii(b1, text, 6, 1, '*', ',', true, true);
    checkb2a(6, text, "     1");
    r = bin2ascii(b1, text, 6, 0, '*', ',', true, true);
    checkb2a(0, text, "");
    
    Vec8i b2 (1,-20,300,4000,50000,654321,7000000,87654321);
    r = bin2ascii(b2, text, 10, 8, '*', ',', true, true);
    checkb2a(87, text, "         1,       -20,       300,      4000,     50000,    654321,   7000000,  87654321");
    r = bin2ascii(b2, text, 9, 8, '*', ',', true, true);
    checkb2a(79, text, "        1,      -20,      300,     4000,    50000,   654321,  7000000, 87654321");
    r = bin2ascii(b2, text, 8, 8, '*', ',', true, true);
    checkb2a(71, text, "       1,     -20,     300,    4000,   50000,  654321, 7000000,87654321");
    r = bin2ascii(b2, text, 7, 8, '*', ',', true, true);
    checkb2a(63, text, "      1,    -20,    300,   4000,  50000, 654321,7000000,*******");
    r = bin2ascii(b2, text, 6, 8, '*', ',', true, true);
    checkb2a(55, text, "     1,   -20,   300,  4000, 50000,654321,******,******");
    r = bin2ascii(b2, text, 5, 8, '*', ',', true, true);
    checkb2a(47, text, "    1,  -20,  300, 4000,50000,*****,*****,*****");
    r = bin2ascii(b2, text, 4, 8, '*', ',', true, true);
    checkb2a(39, text, "   1, -20, 300,4000,****,****,****,****");
    r = bin2ascii(b2, text, 3, 8, '*', ',', true, true);
    checkb2a(31, text, "  1,-20,300,***,***,***,***,***");
    r = bin2ascii(b2, text, 2, 8, '*', ',', true, true);
    checkb2a(23, text, " 1,**,**,**,**,**,**,**");
    r = bin2ascii(b2, text, 1, 8, '*', ',', true, true);
    checkb2a(15, text, "1,*,*,*,*,*,*,*");
    r = bin2ascii(b2, text, 0, 8, '*', ',', true, true);
    checkb2a(0, text, "");

    r = bin2ascii(b2, text, 9, 8, '*', 0, true, true);
    checkb2a(72, text, "        1      -20      300     4000    50000   654321  7000000 87654321");
    r = bin2ascii(b2, text, 8, 8, '*', 0, true, true);
    checkb2a(64, text, "       1     -20     300    4000   50000  654321 700000087654321");
    r = bin2ascii(b2, text, 7, 8, '*', 0, true, true);
    checkb2a(56, text, "      1    -20    300   4000  50000 6543217000000*******");
    r = bin2ascii(b2, text, 6, 8, '*', 0, true, true);
    checkb2a(48, text, "     1   -20   300  4000 50000654321************");
    r = bin2ascii(b2, text, 5, 8, '*', 0, true, true);
    checkb2a(40, text, "    1  -20  300 400050000***************");
    r = bin2ascii(b2, text, 4, 8, '*', 0, true, true);
    checkb2a(32, text, "   1 -20 3004000****************");
    r = bin2ascii(b2, text, 3, 8, '*', 0, true, true);
    checkb2a(24, text, "  1-20300***************");
    r = bin2ascii(b2, text, 2, 8, '*', 0, true, true);
    checkb2a(16, text, " 1**************");
    r = bin2ascii(b2, text, 1, 8, '*', 0, true, true);
    checkb2a(8, text, "1*******");
    r = bin2ascii(b2, text, 0, 8, '*', 0, true, true);
    checkb2a(0, text, "");

    // fields too long
    Vec8i b3 (1000,-200000,3000000,40000000,205050505,3060606060u,-432100000,-87654321);
    r = bin2ascii(b3, text, 8, 8, '*', ',', true, true);
    checkb2a(71, text, "    1000, -200000, 3000000,40000000,********,********,********,********");
    r = bin2ascii(b3, text, 8, 8, 0, ',', true, true);
    checkb2a(78, text, "    1000, -200000, 3000000,40000000,205050505,-1234361236,-432100000,-87654321");
    r = bin2ascii(b3, text, 9, 8, '*', ',', true, true);
    checkb2a(79, text, "     1000,  -200000,  3000000, 40000000,205050505,*********,*********,-87654321");
    r = bin2ascii(b3, text, 10, 8, '*', ',', true, true);
    checkb2a(87, text, "      1000,   -200000,   3000000,  40000000, 205050505,**********,-432100000, -87654321");
    r = bin2ascii(b3, text, 11, 8, '*', ',', true, true);
    checkb2a(95, text, "       1000,    -200000,    3000000,   40000000,  205050505,-1234361236, -432100000,  -87654321");
    r = bin2ascii(b3, text, 12, 8, '*', ',', true, true);
    checkb2a(103, text, "        1000,     -200000,     3000000,    40000000,   205050505, -1234361236,  -432100000,   -87654321");
    r = bin2ascii(b3, text, 10, 8, '*', ',', false, true);
    checkb2a(87, text, "      1000,4294767296,   3000000,  40000000, 205050505,3060606060,3862867296,4207312975");
    r = bin2ascii(b3, text, 5, 1, '*', ',', true, false); // no terminator. overwrite previous string
    checkb2a(87, text, " 1000 1000,4294767296,   3000000,  40000000, 205050505,3060606060,3862867296,4207312975");
    
    if (!globalError) printf("\nsuccess\n");
#endif
    return globalError;
}


#else           // test ascii2bin

// check results of ascii2bin
void checka2b (Vec8i res, Vec8i expected, int length, int lengthExp, int err, int errExp) {
    bool dataerr = horizontal_or(res != expected);
    bool lengtherr = length != lengthExp;
    bool errorerr  = err != errExp;

    if (dataerr || lengtherr || errorerr) {
        printf("\nascii2bin error:");    
    } 
    if (dataerr) {
        globalError |= 1;
        printf("\n  data error:\n    found:  expected:");
        for (int i = 0; i < res.size(); i++) {
            printf("\n%10i %10i", res[i], expected[i]);
        }
    }
    if (lengtherr) {
        globalError |= 2;
        printf("\n  length error: found: %i, expected: %i", length, lengthExp);
    }
    if (errorerr) {
        globalError |= 4;
        printf("\n  error code: found: 0x%X, expected: 0x%X", err, errExp);
    }
    if ((lengtherr || errorerr) && !dataerr) { // print numbers to help identify the case
        printf("\n(");
        for (int i=0; i < res.size(); i++) {
            printf("%i ", res[i]);
        }
        printf(")\n");
    }
}

int main() {

    int error = 0;
    int n_read = 0;
    Vec8i dat;

#if 0  // for debugging only:
    //                      10        20        30        40        50        60
    //                       v         v         v         v         v         v
    //             01234567890123456789012345678901234567890123456789012345678901234567890
    char num0[] = " 1, +2 ,-3, -4321,, 007777, 88888,98765432";

    dat = ascii2bin(num0, &n_read, &error, 1000, 7, ',');
    checka2b (dat, Vec8i(1,2,3,-4321,0,7777,88888,98765432), n_read, 41, error, 0);

    printf ("\nnread %i, error 0x%X\n", n_read, error);
    for (int i=0; i<8; i++) printf("%i ", dat[i]);
    return 1;

#endif
    char num1[] = " 1, +21  , 321, -4321, 55, 7777, 88888,98765432";
    dat = ascii2bin(num1, &n_read, &error, 64, 8, ',');
    checka2b (dat, Vec8i(1,21,321,-4321,55,7777,88888,98765432),n_read, 47, error, 0);

    // test no numbers
    dat = ascii2bin(num1, &n_read, &error, 64, 0, ',');
    checka2b (dat, Vec8i(0), n_read, 0, error, 0);
    // test fewer numbers
    dat = ascii2bin(num1, &n_read, &error, 64, 3, ',');
    checka2b (dat, Vec8i(1,21,321,0,0,0,0,0), n_read, 15, error, 0);
    // test fewer numbers
    dat = ascii2bin(num1, &n_read, &error, 64, 7, ',');
    checka2b (dat, Vec8i(1,21,321,-4321,55,7777,88888,0), n_read, 39, error, 0);

    // test short string
    dat = ascii2bin(num1, &n_read, &error, 40, 7, ',');
    checka2b (dat, Vec8i(1,21,321,-4321,55,7777,88888,0), n_read, 39, error, 0);
    // test short string
    dat = ascii2bin(num1, &n_read, &error, 26, 7, ',');
    checka2b (dat, Vec8i(1,21,321,-4321,55,0,0,0), n_read, 26, error, 8);

    // test string 64 bytes long
    char num2[] = "1     , +22  ,    300, - 4444,   55555, 666666, 7777777,88888888";
    dat = ascii2bin(num2, &n_read, &error, 64, 8, ',');
    checka2b (dat, Vec8i(1,22,300,-4444,55555,666666,7777777,88888888),n_read, 64, error, 0);

    // test string > 64 bytes long
    char num3[] = "1   , +22  ,    300, - 4444,   55555, 666666, 7777777,88888888,999,101010,111111";
    dat = ascii2bin(num3, &n_read, &error, 100, 8, ',');
    checka2b (dat, Vec8i(1,22,300,-4444,55555,666666,7777777,88888888), n_read, 63, error, 0);

    // test missing numbers
    char num4[] = ",- 321,+,-9876543";
    dat = ascii2bin(num4, &n_read, &error, 17, 8, ',');
    checka2b (dat, Vec8i(0,-321,0,-9876543,0,0,0,0), n_read, 17, error, 8);
    // test unfinished numbers
    dat = ascii2bin(num4, &n_read, &error, 14, 8, ',');
    checka2b (dat, Vec8i(0,-321,0,-9876,0,0,0,0), n_read, 14, error, 8);
    dat = ascii2bin(num4, &n_read, &error, 10, 8, ',');
    checka2b (dat, Vec8i(0,-321,0,0,0,0,0,0), n_read, 10, error, 8);

    // test misplaced character and illegal character
    char num5[] = "111 ,  -222 , 333-, 444., 555E6, 6666";
    dat = ascii2bin(num5, &n_read, &error, 37, 6, ',');
    checka2b (dat, Vec8i(111,-222,0,444,0,6666,0,0), n_read, 37, error, 2+4);

    dat = ascii2bin(num5, &n_read, &error, 80, 4, ',');
    checka2b (dat, Vec8i(111,-222,0,444,0,0,0,0), n_read, 25, error, 2+4);

    // test field too long
    char num6[] = "111 ,1234567890, -1234567890  , 4444, 55555, 666666";
    dat = ascii2bin(num6, &n_read, &error, 51, 6, ',');
    checka2b (dat, Vec8i(111,1234567890,-1234567890,4444,55555,666666,0,0), n_read, 51, error, 0);

    // test overflow
    char num7[] = "111 ,12345678901, -1234567890  , 4444, 55555, 666666";
    dat = ascii2bin(num7, &n_read, &error, 64, 6, ',');
    checka2b (dat, Vec8i(111,2147483647,-1234567890,4444,55555,666666,0,0), n_read, 52, error, 16);

    // test chain
    char num8[] = "-111, 222 , -333 , +4444, -55555,+ 666666, -777, 888, -999, 1010, -1111";
    dat = ascii2bin(num8, &n_read, &error, 53, 8, ',');
    checka2b (dat, Vec8i(-111,222,-333,4444,-55555,666666,-777,888), n_read, 53, error, 0);
    dat = ascii2bin(num8 + n_read, &n_read, &error, 64, 3, ',');
    checka2b (dat, Vec8i(-999, 1010, -1111,0,0,0,0,0), n_read, 18, error, 0);

    // test garbage after string. multiple signs, tab as separator
    char num9[] = "111\t+-2\t---3\t4444\t55555\t666666\t-7\t8\t 1.2E3\ttext\t'''\t\t%&/()";
    dat = ascii2bin(num9, &n_read, &error, 64, 8, '\t');
    checka2b (dat, Vec8i(111,0,0,4444,55555,666666,-7,8), n_read, 36, error, 4);

    // test newline as end of string
    char num9a[] = "111,+2,-3,4444,55555,666666,-7\n8, 1.2E3";
    dat = ascii2bin(num9a, &n_read, &error, 64, 8, ',');
    checka2b (dat, Vec8i(111,2,-3,4444,55555,666666,-7,0), n_read, 30, error, 8);

    // test error in first field
    char num10[]= "1 1 1,22,333,4444,55555,666666,7777777,";
    dat = ascii2bin(num10, &n_read, &error, 64, 8, ',');
    checka2b (dat, Vec8i(0,22,333,4444,55555,666666,7777777,0), n_read, 39, error, 4+8);

    char num11[]= "+-+-0,22,333,4444,55555,666666,7777777,";
    dat = ascii2bin(num11, &n_read, &error, 64, 8, ',');
    checka2b (dat, Vec8i(0,22,333,4444,55555,666666,7777777,0), n_read, 39, error, 4+8);

    // test error in last field
    char num12[]= "1 ,22,333,4444,55555,666666,7777777,+-8";
    dat = ascii2bin(num12, &n_read, &error, 64, 8, ',');
    checka2b (dat, Vec8i(1,22,333,4444,55555,666666,7777777,0), n_read, 39, error, 4);

    char num13[]= "1 ,22,333,4444,55555,666666,7777777,8.8";
    dat = ascii2bin(num13, &n_read, &error, 64, 8, ',');
    checka2b (dat, Vec8i(1,22,333,4444,55555,666666,7777777,0), n_read, 39, error, 2+4);

    char num14[]= "1 ,22,333,4444,55555,666666,7777777,...garbage 1 more garbage   ";
    dat = ascii2bin(num14, &n_read, &error, 64, 8, ',');
    checka2b (dat, Vec8i(1,22,333,4444,55555,666666,7777777,1), n_read, 48, error, 2);

    if (!globalError) printf("\nsuccess\n");
    return globalError;
}

#endif

