/***************************  decimal.h   *************************************
* Author:        Agner Fog
* Date created:  2012-07-08
* Last modified: 2019-07-20
* Version:       2.00
* Project:       Extension to vector class library
* Description:
* Functions for conversion between binary number vectors and comma-separated
* decimal ASCII lists.
*
* Please see decimal_manual.pdf for instructions
*
* (c) Copyright 2012-2019 Agner Fog. Apache License version 2.0 or later.
******************************************************************************/

#pragma once
#include "vectorclass.h"

#ifdef VCL_NAMESPACE
namespace VCL_NAMESPACE {
#endif


/*****************************************************************************
*
*               Conversion from binary to decimal ASCII string
*
*****************************************************************************/

// Convert binary numbers to decimal ASCII string.
// The numbers will be written to the string as decimal numbers in human-readable format.
// Each number will be right-justified with leading spaces in a field of the specified length.
int bin2ascii (
    Vec4i const a,           // vector of integers to convert
    char * string,           // string to receive the decimal ascii numbers
    int fieldlen = 8,        // length of each field
    int numdat = 4,          // number of data
    char ovfl = '*',         // overflow indicated by this character. 
                             // ovfl = 0 will make the field wide enough to contain the number
    char separator = ',',    // character to separate fields. 0 for no separator
    bool signd = true,       // data are interpreted as signed integers
    bool term  = true);      // write a zero-terminated string

int bin2ascii (
    Vec8i const a,           // vector of integers to convert 
    char * string,           // string to receive the decimal ascii numbers
    int fieldlen = 8,        // length of each field
    int numdat = 4,          // number of data
    char ovfl = '*',         // overflow indicated by this character. 
                             // ovfl = 0 will make the field wide enough to contain the number
    char separator = ',',    // character to separate fields. 0 for no separator
    bool signd = true,       // data are interpreted as signed integers
    bool term  = true);      // write a zero-terminated string



/*****************************************************************************
*
* Conversion from comma-separated decimal ASCII string to binary number vector
*
*****************************************************************************/

/*
The function ascii2bin shows how it is possible to parse a string of 
variable-length fields without looping through the characters of the sting. 
It is quite a challenge, though. There are many special cases to take care
of and to test. Whether it is worth the effort depends on whether string 
parsing is a bottleneck. In many cases, data transfer is the bottleneck 
that limits the speed, not data parsing. 
This code may serve as a source of inspiration anyway.
*/

Vec8i ascii2bin(
    const char * string,     // ASCII string containing numdat comma-separated integers
    int * chars_read,        // Number of characters read
    int * error,             // Errors will be indicated here
    int max_stringlen = 64,  // Maximum length of string
    int numdat = 8,          // Expected number of data in string. Max 8
    char separator = ',');

// Error codes returned in *error:
// 1:  parameters out of range
// 2:  illegal character.    value will be interpreted as if this was a space
// 4:  misplaced character.  value will be zero
// 8:  too few separators.   value will be zero
// 16: overflow.             value will be INT_MAX or INT_MIN


#ifdef VCL_NAMESPACE
}
#endif
