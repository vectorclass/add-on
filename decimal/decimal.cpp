/***************************  decimal.h   *************************************
* Author:        Agner Fog
* Date created:  2012-07-08
* Last modified: 2019-07-20
* Version:       2.00
* Project:       Add-on package to vector class library
* Description:
* Functions for conversion between binary number vectors and comma-separated
* decimal ASCII lists.
*
* Please see decimal_manual.pdf for instructions
*
* (c) Copyright 2012-2019 Agner Fog. Apache License version 2.0 or later.
******************************************************************************/


#include "decimal.h"

#ifdef VCL_NAMESPACE
namespace VCL_NAMESPACE {
#endif

#if INSTRSET >= 10 && VECTORCLASS_H >= 20000
#define D_COMPACT_BOOLEANS true           // all boolean vectors are compact
#else
#define D_COMPACT_BOOLEANS false          // only the largest boolean vectors are compact
#endif


/*****************************************************************************
*
*               Conversion from binary to BCD
*               These functions are used by bin2ascii
*
*****************************************************************************/

// convert from binary to BCD, for values 0 - 99999999
// no check for overflow
static Vec4ui bin2bcd (Vec4ui const a) {
    // split into halves by dividing by 10000
    Vec4ui thi = a / const_uint(10000);            // a / 10000
    Vec4ui tlo = a - thi * 10000;                  // a % 10000
    Vec4ui split = (thi << 16) | tlo;
    // split into quarters by dividing by 100
    Vec8us thi2 = Vec8us(split) / const_uint(100); // a / 100
    Vec8us tlo2 = Vec8us(split) - thi2 * 100;      // a % 100
    Vec8us split2 = (thi2 << 8) | tlo2;
    // convert each quarter to BCD, using double-dabble algorithm    
    Vec16uc mask3(0x30), mask8(0x80);
    Vec16uc x = Vec16uc(split2);
    for (int p = 4; p > 0; p--) {
        // find nibbles that are bigger than 4
        // (adding 3 will not overflow into next nibble because BCD nibbles are < 10)
        Vec16uc bigger4 = (x + mask3) & mask8;
        // generate value 3 for nibbles that are > 4
        Vec16uc add3 = (bigger4 >> 2) + (bigger4 >> 3);
        // add 3 to generate carry when shifting left
        x += add3;
        // shift masks right instead of shifting x left
        mask3 = mask3 >> 1;  mask8 = mask8 >> 1;
    }
    return Vec4ui(x);
}

// convert from binary to BCD, for values 0 - 99999999
// no check for overflow
static Vec8ui bin2bcd (Vec8ui const a) {
    // split into halves by dividing by 10000
    Vec8ui thi = a / const_uint(10000);               // a / 10000
    Vec8ui tlo = a - thi * 10000;                     // a % 10000
    Vec8ui split = (thi << 16) | tlo;
    // split into quarters by dividing by 100
    Vec16us thi2 = Vec16us(split) / const_uint(100);  // a / 100
    Vec16us tlo2 = Vec16us(split) - thi2 * 100;       // a % 100
    Vec16us split2 = (thi2 << 8) | tlo2;
    // convert each quarter to BCD, using double-dabble algorithm    
    Vec32uc mask3(0x30), mask8(0x80);
    Vec32uc x = Vec32uc(split2);
    for (int p = 4; p > 0; p--) {
        // find nibbles that are bigger than 4
        // (adding 3 will not overflow into next nibble because BCD nibbles are < 10)
        Vec32uc bigger4 = (x + mask3) & mask8;
        // generate value 3 for nibbles that are > 4
        Vec32uc add3 = (bigger4 >> 2) + (bigger4 >> 3);
        // add 3 to generate carry when shifting left
        x += add3;
        // shift masks right instead of shifting x left
        mask3 = mask3 >> 1;  mask8 = mask8 >> 1;
    }
    return Vec8ui(x);
}


/*****************************************************************************
*
*               Conversion from binary to decimal ASCII string
*
*****************************************************************************/

// Convert binary numbers to decimal ASCII string.
// The numbers will be written to the string as decimal numbers in human-readable format.
// Each number will be right-justified with leading spaces in a field of the specified length.
int bin2ascii (Vec4i const a, char * string, int fieldlen, int numdat, 
    char ovfl, char separator, bool signd, bool term) {
    // Parameters:
    // a         vector of 4 numbers to convert
    // string    buffer to receive ASCII string
    // fieldlen  string length of each converted number
    // numdat    number of data elements to output
    // signd     each number will be interpreted as signed if signd is true, unsigned if false. 
    //           Negative numbers will be indicated by a preceding '-'
    // ovfl      Output string will be filled with this character if the number is too big to fit in the field.
    //           The size of a field will be extended in case of overflow if ovfl = 0.
    // separator This character is inserted between data fields, but not after the last field.
    //           The separator character is not included in fieldlen. Separator = 0 for no separator.
    // term      The output string will have a terminating zero ('\0') if term is true.
    // Return value: The length of the written string is returned. The terminating zero is not included in the count.

    if (fieldlen <= 0 || (uint32_t)(numdat - 1) > 4 - 1) {
        // write nothing but terminator
        if (term) *string = 0;
        return 0;
    }
    Vec4ui aa = Vec4ui(a);   // a or abs(a)
    Vec16c space16(' ');     // 16 spaces
#if D_COMPACT_BOOLEANS       // compact boolean vectors
    Vec4ib signa;            // sign of a
    Vec4ib ovfla;            // overflow of a
    int    numwrit;          // number of bytes written to string, not including terminating zero

    // limits depending on fieldlength
    const int limits[9] = { 0,9,99,999,9999,99999,999999,9999999,99999999 };
    int flen = fieldlen < 8 ? fieldlen : 8;                // max fieldlength for vector processing
    if (signd) {                                           // signed
        aa = abs(aa);                                      // abs a
        signa = a < 0;                                     // sign
        ovfla = (a > limits[flen]) | (a < -limits[flen - 1]);// overflow
    }
    else {                                                 // unsigned
        ovfla = aa > limits[flen];                         // overflow
        signa = false;
    }
    if (!(horizontal_or(ovfla) && (fieldlen > 8 || ovfl == 0))) {
        // normal case
        Vec16uc bcd = Vec16uc(bin2bcd(aa));                // converted to BCD
        Vec16uc bcd0246 = bcd & 0x0F;                      // low  nibbles of BCD code
        Vec16uc bcd1357 = bcd >> 4;                        // high nibbles of BCD code
        // interleave nibbles and reverse byte order
        Vec16c  declo = blend16<19, 3, 18, 2, 17, 1, 16, 0, 23, 7, 22, 6, 21, 5, 20, 4>(bcd0246, bcd1357);
        Vec16c  dechi = blend16<27, 11, 26, 10, 25, 9, 24, 8, 31, 15, 30, 14, 29, 13, 28, 12>(bcd0246, bcd1357);
        Vec32c  dec = Vec32c(declo, dechi);                // all digits, big endian digit order
        Vec32c  ascii = dec + 0x30;                        // add '0' to get ascii digits
        // find most significant nonzero digit, or digit 0 if all zero
        Vec32cb decnz = (dec != Vec32c(0));
        uint32_t scan = to_bits(decnz) | 0x80808080;
        // broadcast nonzero indicator bit to the lower bits
        scan |= (scan & 0x7F7F7F7F) << 1;
        scan |= (scan & 0x3F3F3F3F) << 2;
        scan |= (scan & 0x0F0F0F0F) << 4;
        // insert spaces to the left of most significant nonzero digit
        ascii = select(Vec32cb(scan), ascii, Vec32c(' '));
        if (signd) {
            uint32_t minuspos = (scan >> 1) & ~scan;  // position of minus sign
            Vec4uq minussign = select(signa, Vec4uq(Vec32c('-')), Vec4uq(ascii));
            // insert minus sign
            ascii = select(Vec32cb(minuspos), Vec32c(minussign), ascii);
        }
        // insert overflow indicator
        ascii = Vec32c(select(ovfla, Vec4q(Vec32c(ovfl)), Vec4q(ascii)));
#else
    Vec4i  signa;            // sign of a
    Vec4q  signe;            // sign of a, extended
    Vec4ib ovfla;            // overflow of a
    Vec4q  ovfle;            // overflow of a, extended
    int    numwrit;          // number of bytes written to string, not including terminating zero

    // limits depending on fieldlength
    const int limits[9] = { 0,9,99,999,9999,99999,999999,9999999,99999999 };
    int flen = fieldlen < 8 ? fieldlen : 8;
    if (signd) {                                           // signed
        aa = abs(aa);                                      // abs a
        signa = a >> 31;                                   // sign
        ovfla = (a > limits[flen]) | (a < -limits[flen - 1]);// overflow
    }
    else {                                                 // unsigned
        ovfla = aa > limits[flen];                         // overflow
        signa = 0;
    }
    if (!(horizontal_or(ovfla) && (fieldlen > 8 || ovfl == 0))) {
        // normal case
        ovfle = Vec4q(extend_low(Vec4i(ovfla)), extend_high(Vec4i(ovfla)));// overflow, extended
        Vec16uc bcd = Vec16uc(bin2bcd(aa));                // converted to BCD
        Vec16uc bcd0246 = bcd & 0x0F;                      // low  nibbles of BCD code
        Vec16uc bcd1357 = bcd >> 4;                        // high nibbles of BCD code
        // interleave nibbles and reverse byte order
        Vec16c  declo = blend16<19, 3, 18, 2, 17, 1, 16, 0, 23, 7, 22, 6, 21, 5, 20, 4>(bcd0246, bcd1357);
        Vec16c  dechi = blend16<27, 11, 26, 10, 25, 9, 24, 8, 31, 15, 30, 14, 29, 13, 28, 12>(bcd0246, bcd1357);
        Vec32c  dec = Vec32c(declo, dechi);                // all digits, big endian digit order
        Vec32c  ascii = dec + 0x30;                        // add '0' to get ascii digits
        signe = Vec4q(extend_low(signa), extend_high(signa));// sign, extended
        // find most significant nonzero digit, or digit 0 if all zero
        Vec32c decnz = Vec32c(dec != 0) | Vec32c(0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1);
        Vec4q  scan = Vec4q(decnz);
        scan |= scan << 8; scan |= scan << 16; scan |= scan << 32;
        // insert spaces to the left of most significant nonzero digit
        ascii = select(Vec32cb(Vec32c(scan)), ascii, Vec32c(' '));
        if (signd) {
            Vec32c minuspos = Vec32c(andnot(scan >> 8, scan)) & Vec32c(signe);// position of minus sign
            ascii = select(Vec32cb(minuspos), Vec32c('-'), ascii); // insert minus sign
        }
        // insert overflow indicator
        ascii = select(Vec32cb(Vec32c(ovfle)), Vec32c(ovfl), ascii);
#endif
        const int d = V_DC;  // V_DC means don't care in permute functions
        if (separator) {
            numwrit = (fieldlen + 1) * numdat - 1;
            // write output fields with separator
            Vec32c sep(separator);
            if (fieldlen <= 7) {
#if INSTRSET >= 10 && defined (__AVX512VBMI2__)  // more efficient with future AVX512VBMI2 instruction set
                uint32_t fieldbits = (1 << fieldlen) - 1;                         // sequence of fieldlen 1-bits
                uint32_t fieldmask = (fieldbits << (8-fieldlen)) * 0x01010101;    // left-justify and broadcast
                Vec32c   compact = _mm256_maskz_compress_epi8 (fieldmask, ascii); // pick fieldlen chars from each 8-bytes field
                uint32_t expandmask = fieldbits | (fieldbits << (fieldlen + 1));  // indicate positions without comma
                expandmask |= expandmask << (fieldlen + 1) * 2;
                ascii = _mm256_mask_expand_epi8(sep, expandmask, compact);        // distribute fieldlen chars into each field
#else
                switch (fieldlen) {
                case 1:
                    ascii = blend32<7,32,15,32,23,32,31,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d>(ascii, sep);
                    break;
                case 2:
                    ascii = blend32<6,7,32,14,15,32,22,23,32,30,31,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d>(ascii, sep);
                    break;
                case 3:
                    ascii = blend32<5,6,7,32,13,14,15,32,21,22,23,32,29,30,31,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d>(ascii, sep);
                    break;
                case 4:
                    ascii = blend32<4,5,6,7,32,12,13,14,15,32,20,21,22,23,32,28,29,30,31,d,d,d,d,d,d,d,d,d,d,d,d,d>(ascii, sep);
                    break;
                case 5:
                    ascii = blend32<3,4,5,6,7,32,11,12,13,14,15,32,19,20,21,22,23,32,27,28,29,30,31,d,d,d,d,d,d,d,d,d>(ascii, sep);
                    break;
                case 6:
                    ascii = blend32<2,3,4,5,6,7,32,10,11,12,13,14,15,32,18,19,20,21,22,23,32,26,27,28,29,30,31,d,d,d,d,d>(ascii, sep);
                    break;
                case 7:
                    ascii = blend32<1,2,3,4,5,6,7,32,9,10,11,12,13,14,15,32,17,18,19,20,21,22,23,32,25,26,27,28,29,30,31,d>(ascii, sep);
                    break;
                }
#endif
                // store to string
                ascii.store_partial(numwrit, string);
                if (term) string[numwrit] = 0;
            }
            else {
                // fieldlen > 7
                int f;  // field counter
                int c;  // space counter
                // loop for each field
                for (f = 0; f < numdat; f++) {
                    // loop for multiples of 16 spaces
                    for (c = fieldlen - 8; c >= 16; c -= 16) {
                        space16.store(string);  string += 16;
                    }
                    // remaining < 16 spaces
                    space16.store_partial(c, string);  string += c;
                    // insert number (8 digits)
                    ascii.store(string);  string += 8;
                    if (f < numdat-1) {
                        // not last field. insert separator
                        *(string++) = separator;
                        // get next number into position 0
                        ascii = Vec32c(permute4<1,2,3,d>(Vec4q(ascii)));
                    }
                }
                if (term) *string = 0;
            }
        }
        else {
            // write output fields without separator
            numwrit = fieldlen * numdat;
            if (fieldlen <= 8) {
#if INSTRSET >= 10 && defined (__AVX512VBMI2__)
                // more efficient with the future AVX512VBMI2 instruction set
                uint32_t fieldbits = (1 << fieldlen) - 1;                      // sequence of fieldlen 1-bits
                uint32_t fieldmask = (fieldbits << (8-fieldlen)) * 0x01010101; // left-justify and broadcast
                ascii = _mm256_maskz_compress_epi8 (fieldmask, ascii);         // pick fieldlen chars from each 8-bytes field
#else
                switch (fieldlen) {
                case 1:
                    ascii = permute32<7,15,23,31,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d>(ascii);
                    break;
                case 2:
                    ascii = permute32<6,7,14,15,22,23,30,31,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d>(ascii);
                    break;
                case 3:
                    ascii = permute32<5,6,7,13,14,15,21,22,23,29,30,31,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d>(ascii);
                    break;
                case 4:
                    ascii = permute32<4,5,6,7,12,13,14,15,20,21,22,23,28,29,30,31,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d>(ascii);
                    break;
                case 5:
                    ascii = permute32<3,4,5,6,7,11,12,13,14,15,19,20,21,22,23,27,28,29,30,31,d,d,d,d,d,d,d,d,d,d,d,d>(ascii);
                    break;
                case 6:
                    ascii = permute32<2,3,4,5,6,7,10,11,12,13,14,15,18,19,20,21,22,23,26,27,28,29,30,31,d,d,d,d,d,d,d,d>(ascii);
                    break;
                case 7:
                    ascii = permute32<1,2,3,4,5,6,7,9,10,11,12,13,14,15,17,18,19,20,21,22,23,25,26,27,28,29,30,31,d,d,d,d>(ascii);
                    break;
                }
#endif
                // store to string
                ascii.store_partial(numwrit, string);
                if (term) string[numwrit] = 0;
            }
            else {
                // fieldlen > 8
                int f;                    // field counter
                int c;                    // space counter
                // loop for each field
                for (f = 0; f < numdat; f++) {
                    // loop for multiples of 16 spaces
                    for (c = fieldlen - 8; c >= 16; c -= 16) {
                        space16.store(string);  string += 16;
                    }
                    // remaining < 16 spaces
                    space16.store_partial(c, string);  string += c;
                    // insert number (8 digits)
                    ascii.store_partial(8, string);  string += 8;
                    // get next number into position 0
                    ascii = Vec32c(permute4<1,2,3,d>(Vec4q(ascii)));
                }
                if (term) *string = 0;
            }
        }
    }
    else {
        // two special cases are handled here by making one number at a time:
        // (1) more than 8 characters needed
        // (2) one or more fields need to be extended to handle overflow
        Vec4i x;
        int i;
        numwrit = 0;
#if INSTRSET < 4  // SSSE3
        char temp[32];                                     // needed in loop below
#endif
        // loop for numdat data
        // (to do: do this without a loop using future AVX512VBMI2 instruction set)
        for (i = 0; i < numdat; i++) {
            x = permute4<0,-1,-1,-1>(aa);                  // zero-extend abs(a)
            bool big = aa[0] > 99999999u;                  // needs two 64-bit fields
            if (big) {
                // extend into next field
                x = Vec4i(uint32_t(aa[0]) % 100000000u, uint32_t(aa[0]) / 100000000u, 0, 0);
            }
            Vec4ui  bcd = bin2bcd(Vec4ui(x));              // converted to BCD
            Vec16uc bcd0246    = Vec16uc(bcd) & 0x0F;      // low  nibbles of BCD code
            Vec16uc bcd1357    = Vec16uc(bcd) >> 4;        // high nibbles of BCD code
            // interleave nibbles and reverse byte order
            Vec16uc dec        = blend16<-1,-1,-1,-1,-1,-1,20,4,19,3,18,2,17,1,16,0 >(bcd0246, bcd1357);
            Vec16c  ascii      = dec + 0x30;               // add '0' to get ascii digits
            // find most significant nonzero digit, or digit 0 if all zero
            uint16_t decnz = to_bits(dec != Vec16c(0));
            if (big) {
                decnz |= 0xFF80;
            }
            else {
                decnz |= 0x8000;
            } 
            uint16_t scan = decnz;
            // broadcast nonzero indicator bits to the lower bits
            scan |= scan << 1; 
            scan |= scan << 2;
            scan |= scan << 4;

            // insert spaces to the left of most significant nonzero digit
            ascii = select(Vec16cb().load_bits(scan), ascii, space16);
            // count digits
            int bs = bit_scan_reverse(uint32_t(uint16_t(~scan)));
            int numchars = 15 - bs;
            if (signa[i]) {
                uint16_t mpos = 1u << bs;
                ascii = select(Vec16cb().load_bits(mpos), Vec16c('-'), ascii);// insert minus sign
                numchars++;
            }
            int flen2 = fieldlen;
            if (numchars > flen2) {
                // number too big for field
                if (ovfl) {
                    ascii = Vec16c(ovfl);                  // fill whole field with ovfl character
                }
                else {
                    flen2 = numchars;                      // make field longer
                }
            }
            numwrit += flen2;
            // write field
            int c = flen2 - 16;                            // extra spaces needed
            if (c > 0) {                  
                flen2 -= c;
                // loop for multiples of 16 spaces
                for (; c >= 16; c -= 16) {
                    space16.store(string);  string += 16;
                }
                // remaining < 16 spaces
                space16.store_partial(c, string);  string += c;
            }
            if (flen2 < 16) {
                // shift string (16-flen2) characters down.
                // To do: use AVX512VBMI2 if available
#if INSTRSET >= 4  // SSSE3
                const char shufindex[32] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
                Vec16c shiftindex = Vec16c().load(shufindex + (16-flen2));
                ascii = _mm_shuffle_epi8(ascii, shiftindex);
#else 
                ascii.store(temp);
                ascii.load(temp + 16-flen2);
#endif
            }
            // insert number
            ascii.store_partial(flen2, string);  string += flen2;
            if (i < numdat-1) {
                // not last field. insert separator
                if (separator) {
                    *(string++) = separator;
                    numwrit++;
                }
                // get next number into position 0
                aa = permute4<1,2,3,V_DC>(aa);  // shift down next value
            }
        }
        if (term) *string = 0;
    }
    return numwrit;
}


// Convert binary numbers to decimal ASCII string.
// The numbers will be written to the string as decimal numbers in human-readable format.
// Each number will be right-justified with leading spaces in a field of the specified length.
int bin2ascii (Vec8i const a, char * string, int fieldlen, int numdat, 
    char ovfl, char separator, bool signd, bool term) {
    // Parameters:
    // a         vector of 8 numbers to convert
    // string    buffer to receive ASCII string
    // fieldlen  string length of each converted number
    // numdat    number of data elements to output
    // signd     each number will be interpreted as signed if signd is true, unsigned if false. 
    //           Negative numbers will be indicated by a preceding '-'
    // ovfl      Output string will be filled with this character if the number is too big to fit in the field.
    //           The size of a field will be extended in case of overflow if ovfl = 0.
    // separator This character is inserted between data fields, but not after the last field.
    //           The separator character is not included in fieldlen. Separator = 0 for no separator.
    // term      The output string will have a terminating zero ('\0') if term is true.
    // Return value: The length of the written string is returned. The terminating zero is not included in the count.

    if (fieldlen <= 0 || (uint32_t)(numdat - 1) > 8 - 1) {
        // write nothing but terminator
        if (term) *string = 0;
        return 0;
    }
    if (numdat <= 4) {
        return bin2ascii (a.get_low(), string, fieldlen, numdat, ovfl, separator, signd, term);
    }
    Vec8ui aa = Vec8ui(a);   // a or abs(a)
    Vec16c space16(' ');     // 16 spaces
#if D_COMPACT_BOOLEANS       // compact boolean vectors
    Vec8ib signa;            // sign of a
    Vec8ib ovfla;            // overflow of a
    int    numwrit;          // number of bytes written to string, not including terminating zero

    // limits depending on fieldlength
    const int limits[9] = { 0,9,99,999,9999,99999,999999,9999999,99999999 };
    int flen = fieldlen < 8 ? fieldlen : 8;                // max fieldlength for vector processing
    if (signd) {                                           // signed
        aa = abs(aa);                                      // abs a
        signa = a < 0;                                     // sign
        ovfla = (a > limits[flen]) | (a < -limits[flen - 1]);// overflow
    }
    else {                                                 // unsigned
        ovfla = aa > limits[flen];                         // overflow
        signa = false;
    }
    if (!(horizontal_or(ovfla) && (fieldlen > 8 || ovfl == 0))) {
        // normal case
        Vec32uc bcd = Vec32uc(bin2bcd(aa));                // converted to BCD
        Vec32uc bcd0246 = bcd & 0x0F;                      // low  nibbles of BCD code
        Vec32uc bcd1357 = bcd >> 4;                        // high nibbles of BCD code
        // interleave nibbles and reverse byte order
        Vec32c  declo = blend32<35,3,34,2,33,1,32,0, 39,7,38,6,37,5,36,4, 43,11,42,10,41,9,40,8, 47,15,46,14,45,13,44,12>(bcd0246, bcd1357);
        Vec32c  dechi = blend32<51,19,50,18,49,17,48,16, 55,23,54,22,53,21,52,20, 59,27,58,26,57,25,56,24, 63,31,62,30,61,29,60,28>(bcd0246, bcd1357);
        Vec64c  dec = Vec64c(declo, dechi);                // all digits, big endian digit order
        Vec64c  ascii = dec + 0x30;                        // add '0' to get ascii digits
        // find most significant nonzero digit, or digit 0 if all zero
        Vec64cb decnz = (dec != Vec64c(0));
        uint64_t scan = to_bits(decnz) | 0x8080808080808080;
        // broadcast nonzero indicator bit to the lower bits
        scan |= (scan & 0x7F7F7F7F7F7F7F7F) << 1;
        scan |= (scan & 0x3F3F3F3F3F3F3F3F) << 2;
        scan |= (scan & 0x0F0F0F0F0F0F0F0F) << 4;
        // insert spaces to the left of most significant nonzero digit
        ascii = select(Vec64cb().load_bits(scan), ascii, Vec64c(' '));
        if (signd) {
            uint64_t minuspos = (scan >> 1) & ~scan;  // position of minus sign
            Vec8uq minussign = select(signa, Vec8uq(Vec64c('-')), Vec8uq(ascii));
            // insert minus sign
            ascii = select(Vec64cb().load_bits(minuspos), Vec64c(minussign), ascii);
        }
        // insert overflow indicator
        ascii = Vec64c(select(ovfla, Vec8q(Vec64c(ovfl)), Vec8q(ascii)));
#else  // broad boolean vectors
    Vec8i  signa;            // sign of a
    Vec8q  signe;            // sign of a, extended
    Vec8ib ovfla;            // overflow of a
    Vec8q  ovfle;            // overflow of a, extended
    int    numwrit;          // number of bytes written to string, not including terminating zero

    // limits depending on fieldlength
    const int limits[9] = { 0,9,99,999,9999,99999,999999,9999999,99999999 };
    int flen = fieldlen < 8 ? fieldlen : 8;
    if (signd) {                                           // signed
        aa = abs(aa);                                      // abs a
        signa = a >> 31;                                   // sign
        ovfla = (a > limits[flen]) | (a < -limits[flen - 1]);// overflow
    }
    else {                                                 // unsigned
        ovfla = aa > limits[flen];                         // overflow
        signa = 0;
    }
    if (!(horizontal_or(ovfla) && (fieldlen > 8 || ovfl == 0))) {
        // normal case
        ovfle = Vec8q(extend_low(Vec8i(ovfla)), extend_high(Vec8i(ovfla)));// overflow, extended
        Vec32uc bcd = Vec32uc(bin2bcd(aa));                // converted to BCD
        Vec32uc bcd0246 = bcd & 0x0F;                      // low  nibbles of BCD code
        Vec32uc bcd1357 = bcd >> 4;                        // high nibbles of BCD code
        // interleave nibbles and reverse byte order
        Vec32c  declo = blend32<35,3,34,2,33,1,32,0,39,7,38,6,37,5,36,4,43,11,42,10,41,9,40,8,47,15,46,14,45,13,44,12>(bcd0246, bcd1357);
        Vec32c  dechi = blend32<51,19,50,18,49,17,48,16,55,23,54,22,53,21,52,20,59,27,58,26,57,25,56,24, 63,31,62,30,61,29,60,28>(bcd0246, bcd1357);
        Vec64c  dec = Vec64c(declo, dechi);                // all digits, big endian digit order
        Vec64c  ascii = dec + 0x30;                        // add '0' to get ascii digits
        signe = Vec8q(extend_low(signa), extend_high(signa));// sign, extended
        // find most significant nonzero digit, or digit 0 if all zero
        Vec64c decnz = Vec64c(dec != 0) | Vec64c(Vec8uq(0xFF00000000000000u));
        Vec8q  scan = Vec8q(decnz);
        scan |= scan << 8; scan |= scan << 16; scan |= scan << 32;
        // insert spaces to the left of most significant nonzero digit
        ascii = select(Vec64cb(Vec64c(scan)), ascii, Vec64c(' '));
        if (signd) {
            Vec64c minuspos = Vec64c(andnot(scan >> 8, scan)) & Vec64c(signe);// position of minus sign
            ascii = select(Vec64cb(minuspos), Vec64c('-'), ascii); // insert minus sign
        }
        // insert overflow indicator
        ascii = select(Vec64cb(Vec64c(ovfle)), Vec64c(ovfl), ascii);
#endif
        const int d = V_DC;  // V_DC means don't care in permute functions
        if (separator) {
            numwrit = (fieldlen + 1) * numdat - 1;
            // write output fields with separator
            Vec64c sep(separator);
            if (fieldlen <= 7) {
#if INSTRSET >= 10 && defined (__AVX512VBMI2__)
                // more efficient with the future AVX512VBMI2 instruction set
                uint64_t fieldbits = (1 << fieldlen) - 1;   // sequence of fieldlen 1-bits
                uint64_t fieldmask = (fieldbits << (8-fieldlen)) * 0x0101010101010101; // left-justify and broadcast
                Vec64c   compact = _mm512_maskz_compress_epi8 (fieldmask, ascii); // pick fieldlen chars from each 8-bytes field
                uint64_t expandmask = fieldbits | (fieldbits << (fieldlen + 1));  // indicate positions without comma
                expandmask |= expandmask << (fieldlen + 1) * 2;
                expandmask |= expandmask << (fieldlen + 1) * 4;
                ascii = _mm512_mask_expand_epi8(sep, expandmask, compact);        // distribute fieldlen chars into each field
#else
                switch (fieldlen) {
                case 1:
                    ascii = blend64<7,64,15,64,23,64,31, 64,39,64,47,64,55,64,63 ,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d>(ascii, sep);
                    break;
                case 2:
                    ascii = blend64<6,7,64,14,15,64,22,23,64,30,31, 64,38,39,64,46,47,64,54,55,64,62,63 ,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d>(ascii, sep);
                    break;
                case 3:
                    ascii = blend64<5,6,7,64,13,14,15,64,21,22,23,64,29,30,31, 64,37,38,39,64,45,46,47,64,53,54,55,64,61,62,63 ,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d>(ascii, sep);
                    break;
                case 4:
                    ascii = blend64<4,5,6,7,64,12,13,14,15,64,20,21,22,23,64,28,29,30,31, 64,36,37,38,39,64,44,45,46,47,64,52,53,54,55,64,60,61,62,63 ,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d>(ascii, sep);
                    break;
                case 5:
                    ascii = blend64<3,4,5,6,7,64,11,12,13,14,15,64,19,20,21,22,23,64,27,28,29,30,31, 64,35,36,37,38,39,64,43,44,45,46,47,64,51,52,53,54,55,64,59,60,61,62,63 ,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d>(ascii, sep);
                    break;
                case 6:
                    ascii = blend64<2,3,4,5,6,7,64,10,11,12,13,14,15,64,18,19,20,21,22,23,64,26,27,28,29,30,31, 64,34,35,36,37,38,39,64,42,43,44,45,46,47,64,50,51,52,53,54,55,64,58,59,60,61,62,63 ,d,d,d,d,d,d,d,d,d>(ascii, sep);
                    break;
                case 7:
                    ascii = blend64<1,2,3,4,5,6,7,64,9,10,11,12,13,14,15,64,17,18,19,20,21,22,23,64,25,26,27,28,29,30,31, 64,33,34,35,36,37,38,39,64,41,42,43,44,45,46,47,64,49,50,51,52,53,54,55,64,57,58,59,60,61,62,63 ,d>(ascii, sep);
                    break;
                }
#endif
                // store to string
                ascii.store_partial(numwrit, string);
                if (term) string[numwrit] = 0;
            }
            else {
                // fieldlen > 7
                int f;  // field counter
                int c;  // space counter
                // loop for each field
                for (f = 0; f < numdat; f++) {
                    // loop for multiples of 16 spaces
                    for (c = fieldlen - 8; c >= 16; c -= 16) {
                        space16.store(string);  string += 16;
                    }
                    // remaining < 16 spaces
                    space16.store_partial(c, string);  string += c;
                    // insert number (8 digits)
                    ascii.store(string);  string += 8;
                    if (f < numdat-1) {
                        // not last field. insert separator
                        *(string++) = separator;
                        // get next number into position 0
                        ascii = Vec64c(permute8<1,2,3,4,5,6,7,d>(Vec8q(ascii)));
                    }
                }
                if (term) *string = 0;
            }
        }
        else {
            // write output fields without separator
            numwrit = fieldlen * numdat;
            if (fieldlen <= 8) {
#if INSTRSET >= 10 && defined (__AVX512VBMI2__)
                // more efficient with the future AVX512VBMI2 instruction set
                uint64_t fieldbits = (1 << fieldlen) - 1;   // sequence of fieldlen 1-bits
                uint64_t fieldmask = (fieldbits << (8-fieldlen)) * 0x0101010101010101; // left-justify and broadcast
                ascii = _mm512_maskz_compress_epi8 (fieldmask, ascii); // pick fieldlen chars from each 8-bytes field
#else
                switch (fieldlen) {
                case 1:
                    ascii = permute64<7,15,23,31,39,47,55,63,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d>(ascii);
                    break;
                case 2:
                    ascii = permute64<6,7,14,15,22,23,30,31,38,39,46,47,54,55,62,63,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d>(ascii);
                    break;
                case 3:
                    ascii = permute64<5,6,7,13,14,15,21,22,23,29,30,31,37,38,39,45,46,47,53,54,55,61,62,63,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d>(ascii);
                    break;
                case 4:
                    ascii = permute64<4,5,6,7,12,13,14,15,20,21,22,23,28,29,30,31,36,37,38,39,44,45,46,47,52,53,54,55,60,61,62,63,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d>(ascii);
                    break;
                case 5:
                    ascii = permute64<3,4,5,6,7,11,12,13,14,15,19,20,21,22,23,27,28,29,30,31,35,36,37,38,39,43,44,45,46,47,51,52,53,54,55,59,60,61,62,63,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d>(ascii);
                    break;
                case 6:
                    ascii = permute64<2,3,4,5,6,7,10,11,12,13,14,15,18,19,20,21,22,23,26,27,28,29,30,31,34,35,36,37,38,39,42,43,44,45,46,47,50,51,52,53,54,55,58,59,60,61,62,63,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d>(ascii);
                    break;
                case 7:
                    ascii = permute64<1,2,3,4,5,6,7,9,10,11,12,13,14,15,17,18,19,20,21,22,23,25,26,27,28,29,30,31,33,34,35,36,37,38,39,41,42,43,44,45,46,47,49,50,51,52,53,54,55,57,58,59,60,61,62,63,d,d,d,d,d,d,d,d>(ascii);
                    break;
                }
#endif
                // store to string
                ascii.store_partial(numwrit, string);
                if (term) string[numwrit] = 0;
            }
            else {
                // fieldlen > 8
                int f;                    // field counter
                int c;                    // space counter
                // loop for each field
                for (f = 0; f < numdat; f++) {
                    // loop for multiples of 16 spaces
                    for (c = fieldlen - 8; c >= 16; c -= 16) {
                        space16.store(string);  string += 16;
                    }
                    // remaining < 16 spaces
                    space16.store_partial(c, string);  string += c;
                    // insert number (8 digits)
                    ascii.store_partial(8, string);  string += 8;
                    // get next number into position 0
                    ascii = Vec64c(permute8<1,2,3,4,5,6,7,d>(Vec8q(ascii)));
                }
                if (term) *string = 0;
            }
        }
    }
    else {
        // two special cases are handled here by making one number at a time:
        // (1) more than 8 characters needed
        // (2) one or more fields need to be extended to handle overflow
#if INSTRSET < 4  // SSSE3
        char temp[32];                                     // needed in loop below
#endif
        Vec4i x;
        int i;
        numwrit = 0;
        // loop for numdat data
        // (to do: do this without a loop using future AVX512VBMI2 instruction set)
        for (i = 0; i < numdat; i++) {
            x = Vec4i(aa[0],0,0,0);                        // zero-extend abs(a)
            bool big = aa[0] > 99999999u;                  // needs two 64-bit fields
            if (big) {
                // extend into next field
                x = Vec4i(uint32_t(aa[0]) % 100000000u, uint32_t(aa[0]) / 100000000u, 0, 0);
            }
            Vec4ui  bcd = bin2bcd(Vec4ui(x));              // converted to BCD
            Vec16uc bcd0246    = Vec16uc(bcd) & 0x0F;      // low  nibbles of BCD code
            Vec16uc bcd1357    = Vec16uc(bcd) >> 4;        // high nibbles of BCD code
            // interleave nibbles and reverse byte order
            Vec16uc dec        = blend16<-1,-1,-1,-1,-1,-1,20,4,19,3,18,2,17,1,16,0 >(bcd0246, bcd1357);
            Vec16c  ascii      = dec + 0x30;               // add '0' to get ascii digits
            // find most significant nonzero digit, or digit 0 if all zero
            uint16_t decnz = to_bits(dec != Vec16uc(0));
            if (big) {
                decnz |= 0xFF80;
            }
            else {
                decnz |= 0x8000;
            } 
            uint16_t scan = decnz;
            // broadcast nonzero indicator bits to the lower bits
            scan |= scan << 1; 
            scan |= scan << 2;
            scan |= scan << 4;
            // insert spaces to the left of most significant nonzero digit
            ascii = select(Vec16cb().load_bits(scan), ascii, space16);
            // count digits
            int bs = bit_scan_reverse(uint32_t(uint16_t(~scan)));
            int numchars = 15 - bs;
            if (signa[i]) {
                uint16_t mpos = 1u << bs;
                ascii = select(Vec16cb().load_bits(mpos), Vec16c('-'), ascii);// insert minus sign
                numchars++;
            }
            int flen2 = fieldlen;
            if (numchars > flen2) {
                // number too big for field
                if (ovfl) {
                    ascii = Vec16c(ovfl);                  // fill whole field with ovfl character
                }
                else {
                    flen2 = numchars;                      // make field longer
                }
            }
            numwrit += flen2;
            // write field
            int c = flen2 - 16;                            // extra spaces needed
            if (c > 0) {                  
                flen2 -= c;
                // loop for multiples of 16 spaces
                for (; c >= 16; c -= 16) {
                    space16.store(string);  string += 16;
                }
                // remaining < 16 spaces
                space16.store_partial(c, string);  string += c;
            }
            if (flen2 < 16) {
                // shift string 16-flen2 characters down.
                // To do: use AVX512VBMI2
#if INSTRSET >= 4  // SSSE3
                const char shufindex[32] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
                Vec16c shiftindex = Vec16c().load(shufindex + (16-flen2));
                ascii = _mm_shuffle_epi8(ascii, shiftindex);
#else 
                ascii.store(temp);
                ascii.load(temp + 16-flen2);
#endif
            }
            // insert number
            ascii.store_partial(flen2, string);  string += flen2;
            if (i < numdat-1) {
                // not last field. insert separator
                if (separator) {
                    *(string++) = separator;
                    numwrit++;
                }
                // get next number into position 0
                aa = permute8<1,2,3,4,5,6,7,V_DC>(aa);  // shift down next value
            }
        }
        if (term) *string = 0;
    }
    return numwrit;
}


/*****************************************************************************
*
* Conversion from comma-separated decimal ASCII string to binary number vector
*
*****************************************************************************/

/*
The function ascii2bin shows how it is possible to parse a string of 
variable-length fields without looping through the characters of the sting. 
The method is to indicate the positions of digits, signs, and commas using
64-bit integers as bitfields. Tricks with propagating carries are then 
used to find where each number begins and ends.

This is quite a challenge, though. There are many special cases to take care
of and to test. Whether it is worth the effort depends on whether string 
parsing is a bottleneck. In many cases, data transfer is the bottleneck 
that limits the speed, not data parsing. 
This code may serve as a source of inspiration anyway.
*/

// fallback function for ascii2bin when parallel processing fails
static Vec8i ascii2bin_fallback(const char * string, int * chars_read, int * error, 
    int max_stringlen, int numdat, char separator);

Vec8i ascii2bin(
    const char * string,       // ASCII string containing numdat comma-separated integers
    int * chars_read,          // Number of characters read
    int * error,               // Errors will be indicated here
    int max_stringlen,         // Maximum length of string
    int numdat,                // Expected number of data in string. Max 8
    char separator) {          // Separator character between data
    // Error codes returned in *error:
    // 1:  parameters out of range
    // 2:  illegal character.   value will be interpreted as if this was a space
    // 4:  misplaced character. value will be zero
    // 8:  too few separators.  value will be zero
    // 16: overflow.            value will be INT_MAX or INT_MIN

    if ((uint32_t)numdat > 8 || (uint32_t)max_stringlen > 10000) {
        if (chars_read) *chars_read = 0;
        *error = 1; return 0;
    }
    if (numdat == 0) {
        if (chars_read) *chars_read = 0;
        return 0;
    }
#if INSTRSET < 5  // no 32-bit multiply instruction. fallback function may be faster
    return ascii2bin_fallback(string, chars_read, error, max_stringlen, numdat, separator);
#endif
    char separator1 = separator;
    int err = 0;                                 // reset error indicator
    Vec64c str;                                  // vector containing string
    str.load_partial(max_stringlen, string);     // read string
    if (separator1 <= ' ') {
        // separator is a control character. Replace by comma
        str = select(str == separator1, Vec64c(','), str);
        separator1 = ',';
    }
    // Replace newline and other control characters by end of string
    str = select(str < Vec64c(' '), Vec64c(0), str);
    // check if string contains a terminating zero
    uint64_t termz = to_bits(str == 0);
    if (termz) {
        int strlength = bit_scan_forward(termz);
        if (max_stringlen > strlength) {
            max_stringlen = strlength;           // reduce length
        }
    }
    // remove any garbage after end of string
    if (max_stringlen < 64) {
        str.cutoff(max_stringlen);
    }

    // find positions of all commas
    uint64_t commas = to_bits(str == separator1);// one bit for each comma
    // check if enough fields
    int num_commas = (int)vml_popcnt(commas);    // count commas
    int numdatr = numdat;
    if (num_commas < numdat) {
        if (max_stringlen > 64) {
            // there are not numdat data fields within 64 chars. cannot use vector calculation        
            return ascii2bin_fallback(string, chars_read, error, max_stringlen, numdat, separator);
        }
        if (num_commas + 1 < numdat) {
            numdatr = num_commas + 1;  // don't attempt to read more than one number after the last comma
            err |= 8;
        }
    }

    // identify tokens in string
    uint64_t digits = to_bits(Vec64uc(str) - Vec64uc('0') <= Vec64uc(9)); // position of digits 0-9
    uint64_t signs  = to_bits(str == '+' | str == '-');                   // position of +/- signs
    uint64_t space  = to_bits(Vec64c(str & Vec64c(~0x20)) == 0);          // position of space/blanks

    // check for illegal characters
    uint64_t illeg = ~(digits | commas | signs | space);                  // illegal characters

#if INSTRSET >= 10 && defined (__AVX512VBMI2__)
    Vec64c const sequence(   // this sequence is used later as well
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,
        32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63);

    // Warning: Bug here in MS VS2019 with optimization on
    if (illeg) {
        // cut off any illegal characters after last comma
        if (num_commas >= numdatr) {
            // find position of comma number numdatr
            Vec64uc commalist = _mm512_maskz_compress_epi8(commas, sequence);  // get bit positions
            int strlength2 = commalist[numdatr-1] + 1;
            // remove any garbage after end of string
            if (strlength2 < max_stringlen) {
                max_stringlen = strlength2;
                str.cutoff(strlength2);
                uint64_t mm = ((uint64_t)1 << strlength2) - 1;
                illeg &= mm;           // cut off illegal characters after the end
            }
            if (illeg) err |= 2;       // still error
        }
        else {
            err |= 2;
        }
        // remove illegal characters
        str = _mm512_maskz_mov_epi8(~illeg, str);
    }
#else
    if (illeg) {
        // cut off any illegal characters after last comma
        if (num_commas >= numdatr) {
            // find position of comma number numdatr
            uint64_t commared = commas;
            for (int i = 1; i < numdatr; i++) { // remove lowest numdatr - 1 bits from commas
                commared &= commared - 1;       // remove one bit
            }
            int strlength2 = bit_scan_forward(commared) + 1;
            // remove any garbage after end of string
            if (strlength2 < max_stringlen) {
                max_stringlen = strlength2;
                str.cutoff(strlength2);
                uint64_t mm = ((uint64_t)1 << strlength2) - 1;
                illeg &= mm;                     // cut off illegal characters after the end
            }
            if (illeg) err |= 2;                 // still error after cutting off superfluous characters
        }
        else {
            err |= 2;                            // report error: illegal characters
        }
        // remove illegal characters
        str = select(Vec64cb().load_bits(illeg), 0, str);
    }
#endif
    // store string with 8 empty chars before string
    union UU {
        char      s[64+8];                       // char string
        long long L[9];                          // intrinsic gather functions use long long *
    } u;
    u.L[0] = 0;                                  // fill first 8 chars with zeroes
    str.store(u.s + 8);                          // put data string next

    // Parse string. Find begin and end of each number
    // Note that decimal ASCII numbers are stored in big endian format:
    // The least significant digit is at the highest address
    uint64_t btw = ((digits | commas) - signs) & ~ commas; // fill any space between sign and digits
    uint64_t num = digits | btw;                 // whole number, including sign, and any space between sign and digits
    uint64_t cn  = commas - num;
    uint64_t begn = cn & num;                    // position of begin of number
    uint64_t emptn = cn & commas;                // find empty fields between commas
    uint64_t endn = (num + begn) >> 1;           // position of end of number
    begn |= emptn;                               // add empty fields to list of begin positions
    endn |= num & (1ull << 63);                  // fix any overflow of endn
    endn |= emptn;                               // add empty fields to list of end positions
    int end_string = 0;                          // find end of string

#if INSTRSET >= 10 && defined (__AVX512VBMI2__)  
    // Make list of begin and end of each field
    // The forthcoming AVX512VBMI2 enables us to avoid the loop through field positions.
    Vec64uc beglist1 = _mm512_maskz_compress_epi8(begn, sequence);  // get bit positions
    Vec64uc endlist1 = _mm512_maskz_compress_epi8(endn, sequence);  // get bit positions
    Vec64uc field_size = endlist1 - beglist1;                       // size of each field
    field_size.cutoff(numdatr);                                     // cut off unused fields
    Vec8i end_index = _mm256_cvtepi8_epi32(_mm512_castsi512_si128(endlist1));
    end_string = endlist1[numdatr-1];           // position of end of last number
    if (end_string == 0 && numdatr > 1) end_string = endlist1[numdatr-2]; // last field missing

    // Certain syntax errors, such as multiple signs or space between digits will cause
    // extra bits in endn. To detect these errors, we will look for more endn positions
    // before any terminating comma
    uint64_t comm2 = commas >> end_string;       // any remaining commas
    int remain_char = 0;                         // any remaining characters
    int end_string2 = max_stringlen;             // space to look for syntax errors in
    if (comm2) {
        remain_char = bit_scan_forward(comm2);   // find terminating comma
        end_string += remain_char;               // make chars_read include terminating comma
        end_string2 = end_string;
    }
    int nextend = endlist1[numdatr];             // next unused position in endlist1
    if (nextend > 0 && nextend <= end_string2) { // there are extra endn positions before terminating comma
        endn = 1;                                // use this as error indicator to call fallback function
    }
    else {
        endn = 0;
    }
#else
    // make list of positions of first and last character in each field, without spaces
    int32_t beglist[8] = {0};
    int32_t endlist[8] = {0};

    for (int f = 0; f < numdatr; f++) {          // loop to find each field
        uint64_t begn2 = begn & begn - 1;        // remove lowest bit
        uint64_t bpos = begn & ~ begn2;
        if (bpos) beglist[f] = bit_scan_forward(bpos);// index to this bit = begin of number
        uint64_t endn2 = endn & endn - 1;        // remove lowest bit
        uint64_t epos = endn & ~ endn2;
        if (epos == 0) break;
        endlist[f] = bit_scan_forward(epos);// index to this bit = end of number
        begn = begn2;                            // remove this bit to find the next bit
        endn = endn2;                            // remove this bit to find the next bit
    }
    Vec8i  beg_index = Vec8i().load_partial(numdatr, beglist); // load list
    Vec8i  end_index = Vec8i().load_partial(numdatr, endlist); // load list
    Vec8ui field_size = Vec8ui(end_index) - Vec8ui(beg_index); // this is actually the field size - 1
    end_string = endlist[numdatr-1];             // position of end of last number
    if (end_string == 0 && numdatr > 1) {        // last field missing
        end_string = endlist[numdatr-2];         // get previous position instead
    }

    uint64_t comm2 = commas >> end_string;       // any remaining commas
    endn >>= end_string;                         // test if there are any unused bits in endn
    int remain_char = 0;                         // any remaining characters
    if (comm2) {
        remain_char = bit_scan_forward(comm2);   // find terminating comma
        endn &= ((uint64_t)1 << remain_char)-1;  // cut off bits beyond terminating comma
        end_string += remain_char;               // make chars_read include terminating comma
    }
#endif

    // Certain syntax errors, such as multiple signs or space between digits will cause
    // extra bits in endn. We are using the fallback function to handle these
    // syntax errors in order to get the right values in the non-error fields.
    // We are also using the fallback function if any number is more than 8 characters
    // because they don't fit into the vectors
    if (horizontal_or(field_size > 7) || endn != 0) {
        // Use fallback function
        return ascii2bin_fallback(string, chars_read, error, max_stringlen, numdat, separator);
    }
    end_string += 1;                             // number of characters read
    if (chars_read) *chars_read = end_string;    // save number of characters read

    // We need to distinguish missing fields from the first field at position zero.
    // Missing fields will be pointed to the initial zero block.
    // If the first field is missing, it will be zero anyway
    Vec8ib missing_fields = (end_index == 0) & Vec8ib(0,1,1,1,1,1,1,1);
    // right-justify numbers with the last significant digit at the end of the block.
    // each read block ends at the end index and begins 8 chars before
    Vec8i  ix = select(missing_fields, -8, end_index - 7);

    // gather 8-character fields, possibly overlapping
#if INSTRSET >= 9      // AVX512. Gather 8 fields in one instruction
    Vec8q fields = _mm512_mask_i32gather_epi64 (Vec8q(0), __mmask8((1 << numdatr) - 1), ix, u.L+1, 1);

#elif INSTRSET >= 8    // AVX2. Gather 2x4 fields
    uint8_t mask = (1 << numdatr) - 1;  // bit mask for read
    Vec4qb m0 = Vec4qb().load_bits(mask);    
    Vec4q fields0 = _mm256_mask_i32gather_epi64 (Vec4q(0), u.L+1, ix.get_low(), m0, 1);
    Vec4q fields1(0);
    if (numdatr > 4) {
        Vec4qb m1 = Vec4qb().load_bits(mask >> 4);   
        fields1 = _mm256_mask_i32gather_epi64 (Vec4q(0), u.L+1, ix.get_high(), m1, 1);
    }
    Vec8q fields(fields0, fields1);              // join fields
#else
    // no gather instruction available, read fields one by one
    uint64_t fieldsa[8];
    int f;
    for (f = 0; f < numdatr; f++) {
        fieldsa[f] = *(uint64_t*)(u.s + 8 + ix[f]);
    }
    for (; f < 8; f++) {
        fieldsa[f] = 0;                          // set remaining fields to zero
    }
    Vec8q fields = Vec8q().load(fieldsa);        // join 64-bit fields into vector
#endif
    // Find commas in fields. 
    // Everyting after a comma belongs to the next field and must be deleted
    Vec64cb ecommas = Vec64c(fields) == separator1;
#if INSTRSET >= 10   // AVX512BW: ecommas is a compact boolean vector
    uint64_t re = to_bits(ecommas);
    re |= (re & 0xFEFEFEFEFEFEFEFE) >> 1;        // remove everything after comma
    re |= (re & 0xFCFCFCFCFCFCFCFC) >> 2;
    re |= (re & 0xF0F0F0F0F0F0F0F0) >> 4;
    fields = _mm512_maskz_mov_epi8(~re, fields);
#else  // ecommas is a broad boolean vector
    Vec8q remove = Vec8q(Vec64c(ecommas));
    remove |= remove >>  8;                      // remove everything after comma
    remove |= remove >> 16;
    remove |= remove >> 32;
    fields = andnot(fields, remove);
#endif

    // Find '+' and '-' signs
    Vec64cb minussigns = Vec64c(fields) == '-';
    Vec64cb plussigns  = Vec64c(fields) == '+';
    // Remove '+' and '-' signs
    fields = Vec8q(select(minussigns | plussigns, Vec64c(0), Vec64c(fields)));

    // Find digits '0' - '9'
    Vec64cb digitpos = Vec64uc(Vec64c(fields)) - Vec64uc(0x30) <= Vec64uc(9);

    // Check that we have nothing to the right of digits or between digits
#if INSTRSET >= 10     // AVX512BW: digitpos is a compact boolean vector
    uint64_t digp = to_bits(digitpos);
    uint64_t misplaced = ((digp & 0x7F7F7F7F7F7F7F7F) << 1) & ~ digp;
    if (misplaced) {   // possibly a misplaced sign after the digits
        // find corresponding field and set value to zero
        Vec64c mispos = _mm512_movm_epi8(misplaced);
        fields = select(Vec8q(mispos) != 0, Vec8q(0), fields);
        err |= 4;
    }
#else   // digitpos is a broad boolean vector
    Vec8uq digp = Vec8uq(digitpos);
    Vec8uq misplaced = Vec8uq(andnot(digp << 8 , digp));
    if (horizontal_or(misplaced != 0)) {  
        // set field to zero if it contains a misplaced character
        fields = select(misplaced != 0, Vec8q(0), fields);
        err |= 4;
    }
#endif
    // convert ASCII digits to binary. Spaces become zeroes
    fields = Vec8q(Vec64c(fields) & Vec64c(~0x30));

    // Convert to binary by multiplying
    // Note that the digit order is big endian within each 8-byte field
    // The code below takes the reversed digit order into account
    Vec64uc string1 = Vec64uc(Vec64c(fields));
    
    // Multiply even-numbered digits by 10
    Vec64uc string2 = (string1 << 1) + (string1 << 3);  // there is no efficient 8-bit multiply. shift and add instead
    
    // Add odd-numbered digits to even-numbered digits * 10
    Vec32us string3 = (Vec32us(string1) >> 8) + (Vec32us(string2) & Vec32us(0x00FF));

    // Now we have 8x4 blocks of 00-99. Multiply every second number by 100
    Vec32us string4 = Vec32us(string3 * Vec32us(Vec16ui(100 + (1 << 16))));
    string4 += Vec32us(Vec16ui(string4) >> 16);    // add
    // zero-extend
#if INSTRSET >= 10
    string4 = _mm512_maskz_mov_epi16 (0x55555555, string4);
#else
    string4 &= Vec32us(Vec16ui(0x0000FFFF));
#endif
    // Now we have 8x2 blocks of 0000-9999. Multiply every second number by 10000
    Vec16ui mu =  Vec16ui(Vec8uq(10000 + (1ull << 32)));
    Vec16ui string5 = Vec16ui(Vec16ui(string4) * mu);
    string5 += Vec16ui(Vec8uq(string5) >> 32); // add

    // apply minus sign, and get even-numbered elements only
#if INSTRSET >= 10  // AV512BW
    Vec64c minus  = _mm512_movm_epi8(minussigns);
    Vec8q string6 = select(Vec8q(minus) != 0, -Vec8q(string5), Vec8q(string5));
    // convert to 32-bit integers. the upper 32 bits are truncated
    Vec8i string7 = _mm512_cvtepi64_epi32(string6);
#elif INSTRSET == 9  // AV512 (512 bit vectors have compact booleans, smaller vectors have broad booleans)
    Vec8q minus = Vec8q(minussigns);
    Vec8qb minus1 = minus != 0;      // detect minus sign anywhere in each 64-bits block
    // apply minussign to 8 64-bit integers, only the lower 32 bits of each are valid
    Vec8q string6 = select(minus1, -Vec8q(string5), Vec8q(string5));
    // convert to 32-bit integers. the upper 32 bits of each 64 bit integer is truncated
    Vec8i string7 = _mm512_cvtepi64_epi32(string6);
#else
    Vec8q minus = Vec8q(minussigns);
    Vec8qb minus1 = minus != 0;      // detect minus sign anywhere in each 64-bits block
    // apply minussign to 8 64-bit integers, only the lower 32 bits of each are valid
    Vec8q string6 = select(minus1, -Vec8q(string5), Vec8q(string5));
    Vec8i s60 = Vec16i(string6).get_low();
    Vec8i s61 = Vec16i(string6).get_high();
    Vec8i string7 = blend8<0,2,4,6,8,10,12,14>(s60,s61); // get lower half of each 64-bit block
#endif

    *error = err;            // error indicator
    return string7;          // return vector
}

// fallback function for ascii2bin if parallel processing fails
static Vec8i ascii2bin_fallback(const char * string, int * chars_read, int * error, 
    int max_stringlen, int numdat, char separator) {
    int err = 0;             // error
    int s = 0;               // string index
    int state = 0;           // 0: begin of number, 
                             // 1: after sign, 
                             // 2: after digit, 
                             // 3: trailing space, 
                             // 4: error
    int n = 0;               // index to number being interpreted
    int last_pos = 0;        // end of last number
    bool negative = false;   // number is negative
    int32_t results[8] = {0};// result array
    int64_t val = 0;         // current value

#ifdef  _DEBUG   // for debugging only. You may remove this
    // printf("\nFallback");
#endif

    // loop through string
    for (s = 0; s < max_stringlen; s++) {
        char c = string[s];
        if ((uint8_t)c < 0x20u && c != separator) break;  // end of string
        switch (state) {
        case 0:  // 0: begin of number
            if (c == '+') {
                state = 1;
            }
            else if (c == '-') {
                negative = true; state = 1;
            }
            else if (c == separator) {
                goto STORE_RESULT;
            }
            else if (c >= '0' && c <= '9') {
                goto DIGIT;
            }
            else if (c != ' ') {
                err |= 2;
            }
            break;
        case 1:  // 1: after sign
            if (c == separator) {
                goto STORE_RESULT;
            }
            else if (c >= '0' && c <= '9') {
                goto DIGIT;
            }
            else if (c == '-' || c == '+') {
                err |= 4;
                state = 4;
                // negative = true;
            }
            else if (c != ' ') {
                err |= 2;
            }
            break;
        case 2:  // 2: after digit
            if (c >= '0' && c <= '9') {
            DIGIT:
                val = val * 10 + (c - '0');
                state = 2;
            }
            else if (c == separator || s == max_stringlen - 1) {
            STORE_RESULT:
                if (negative) val = -val;
                if (val > 0x7FFFFFFF || val < int32_t(0x80000000)) {
                    val = 0x7FFFFFFFu + (uint32_t)negative;// overflow
                    err |= 16;
                }
                results[n++] = int32_t(val);
                val = 0;  negative = false;  state = 0;    // prepare for next
            }
            else if (c == ' ') {
                last_pos = s;
                state = 3;
            }
            else if (c == '-' || c == '+') {
                err |= 4;
                state = 4;
            }
            else {
                err |= 2;
                state = 3;  // illegal character. treat as space
            }
            break;
        case 3:  // 3: trailing space. expecting separator
            if (c == separator) {
                goto STORE_RESULT;
            }
            else if (c == '+' || c == '-' || (c >= '0' && c <= '9')) {
                err |= 4;
                state = 4;
            }
            else if (c != ' ') {
                err |= 2;
            }
            break;
        case 4:  // 4: error: misplaced character. set value to zero
            val = 0;
            if (c == separator) {
                goto STORE_RESULT;
            }
        }
        if (n >= numdat) {
            if (numdat) s++;
            break;
        }
    }
    if (state > 0 && n < numdat) {     // pending number without terminating comma
        if (state == 4) val = 0;
        if (state == 3) s = last_pos;  // trailing space is not included in chars_read
        results[n++] = negative ? -int32_t(val) : int32_t(val);
        if (val > 0x7FFFFFFF || val < int32_t(0x80000000)) err |= 16; // overflow
    }
    if (n < numdat) err |= 8;          // fields missing
    *error = err;
    if (chars_read) {                  // get number of characters read
        if (s > 0 && string[s-1] == 0) s--;  // don't point beyond end of string
        *chars_read = s;
    }
    return Vec8i().load(results);
}

#ifdef VCL_NAMESPACE
}
#endif
