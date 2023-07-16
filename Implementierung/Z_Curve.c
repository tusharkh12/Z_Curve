#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <immintrin.h>
#include <emmintrin.h>
#include <unistd.h>
#include <ctype.h>
#include <getopt.h>

static struct option long_options[] = {
        {"V", required_argument, NULL, 'V'},
        {"B", optional_argument, NULL, 'B'},
        {"t", no_argument, NULL, 't'},
        {"c", no_argument, NULL, 'c'},
        {"d", required_argument, NULL, 'd'},
        {"p", no_argument, NULL, 'p'},
        {"i", required_argument, NULL, 'i'},
        {"a", no_argument, NULL, 'a'},
        {"h", no_argument, NULL, 'h'},
        {NULL, 0, NULL, 0}
};

typedef int coord_t;

void z_curve_recursive(unsigned degree, coord_t* x, coord_t* y);
void z_curve_recursive_helper(unsigned degree, coord_t start_x, coord_t start_y, coord_t* x, coord_t* y, unsigned* index);
void z_curve_iterative(unsigned degree, coord_t* x, coord_t* y);
void z_curve_iterative_simd(unsigned degree, coord_t* x, coord_t* y);
void z_curve_magic(unsigned degree, coord_t* x, coord_t* y);
void z_curve_magic_simd(unsigned degree, coord_t* x, coord_t* y);

void z_curve_at(unsigned degree, size_t idx, coord_t* x, coord_t* y);
void z_curve_at_iterative(unsigned degree, size_t idx, coord_t* x, coord_t* y);
void z_curve_at_magic(unsigned degree, size_t idx, coord_t* x, coord_t* y);

size_t z_curve_pos(unsigned degree, coord_t x, coord_t y);
size_t z_curve_pos_iterative(unsigned degree, coord_t x, coord_t y);
size_t z_curve_pos_magic(unsigned degree, coord_t x, coord_t y);

coord_t ReversePart1By1(coord_t z);
coord_t Part1By1(coord_t x);
__m128i ReversePart1By1_simd(__m128i z);

void tester3();

void print_usage();
void testRun_Method2();
void testRun_Method1();
void describe_Vmethod_option();
void rahmenProgram_helper( unsigned implementation ,unsigned repetitions,unsigned degree,coord_t* xCoords,
                           coord_t* yCoords,int measure_runtime);
void rahmen_svg(coord_t* xCoords, coord_t* yCoords,unsigned numberOfPoints);
void addDefaultArgument(char *argv[], int *argc);
void rahmenProgram(int argc, char *argv[]);




int main(int argc, char *argv[]) {
    rahmenProgram(argc, argv);
    return 0;
}



//inline methods
inline coord_t ReversePart1By1(coord_t z) {
    z &= 0x55555555;                  
    z = (z ^ (z >> 1)) & 0x33333333;  
    z = (z ^ (z >> 2)) & 0x0f0f0f0f;  
    z = (z ^ (z >> 4)) & 0x00ff00ff;  
    z = (z ^ (z >> 8)) & 0x0000ffff;  

    return z;
}

inline coord_t Part1By1(coord_t x) {
  x &= 0x0000ffff;                  
  x = (x ^ (x <<  8)) & 0x00ff00ff; 
  x = (x ^ (x <<  4)) & 0x0f0f0f0f; 
  x = (x ^ (x <<  2)) & 0x33333333; 
  x = (x ^ (x <<  1)) & 0x55555555; 
  return x;
}

inline __m128i ReversePart1By1_simd(__m128i z) {
    const __m128i mask1 = _mm_set1_epi32(0x55555555);
    const __m128i mask2 = _mm_set1_epi32(0x33333333);
    const __m128i mask3 = _mm_set1_epi32(0x0f0f0f0f);
    const __m128i mask4 = _mm_set1_epi32(0x00ff00ff);
    const __m128i mask5 = _mm_set1_epi32(0x0000ffff);

    z = _mm_and_si128(z, mask1);
    z = _mm_xor_si128(z, _mm_srli_epi32(z, 1));
    z = _mm_and_si128(z, mask2);
    z = _mm_xor_si128(z, _mm_srli_epi32(z, 2));
    z = _mm_and_si128(z, mask3);
    z = _mm_xor_si128(z, _mm_srli_epi32(z, 4));
    z = _mm_and_si128(z, mask4);
    z = _mm_xor_si128(z, _mm_srli_epi32(z, 8));
    z = _mm_and_si128(z, mask5);

    return z;
}

//1st Methods
void z_curve_recursive(unsigned degree, coord_t* x, coord_t* y) {
    //Calculate the total number of points
    unsigned size = 1 << (2 * degree);
    //index for buffer
    unsigned index = 0;

    //starting coordinates
    coord_t start_x = 0;
    coord_t start_y = 0;

    z_curve_recursive_helper(degree, start_x, start_y, x, y, &index);
}

void z_curve_recursive_helper(unsigned degree, coord_t start_x, coord_t start_y, coord_t* x, coord_t* y, unsigned* index) {
    if (degree == 0) {
        x[*index] = start_x;
        y[*index] = start_y;
        (*index)++;
        return;
    }

    unsigned sub_size = 1 << (degree - 1);
    //unsigned sub_size_half = sub_size >> 1;

    z_curve_recursive_helper(degree - 1, start_x, start_y, x, y, index);
    z_curve_recursive_helper(degree - 1, start_x + sub_size, start_y, x, y, index);
    z_curve_recursive_helper(degree - 1, start_x, start_y + sub_size, x, y, index);
    z_curve_recursive_helper(degree - 1, start_x + sub_size, start_y + sub_size, x, y, index);
}

void z_curve_iterative(unsigned degree, coord_t* x, coord_t* y) {
    //calculate the total number of points
    //2^(2*degree)
    unsigned total_points = 1 << (2 * degree);

    //iterate over each point
    for (unsigned i = 0; i < total_points; i++) {
        //index and starting coordinate values of x and y for the current point
        unsigned n = i;
        unsigned x2 = 0;
        unsigned y2 = 0;

        //calculate coordinates for current point
        for (unsigned s = 0; s < degree; s++) {
            //bitwise-or of x2 with the least significant bit of n
            //which we get via (n & 1)
            //this is then left-shifted by s positions
            //basically the opposite of interleaving
            //similar approach in z_curve_at2
            x2 |= (n & 1) << s;

            //bitwise-or of y2 with the second least significan bit of n
            //which we get via (n & 2)
            //which we then shift once to the right before again
            //left-shifting by s positions
            //basically the opposite of interleaving
            //similar approach in z_curve_at2
            y2 |= ((n & 2) >> 1) << s;

            //discard the two bits used for x2 and y2 by shifting twice to the right
            n >>= 2;
        }

        //store the calculated coordinates
        x[i] = x2;
        y[i] = y2;
    }
}

void z_curve_iterative_simd(unsigned degree, coord_t* x, coord_t* y) {
    //calculate the total number of points
    //2^(2*degree)
    unsigned total_points = 1 << (2 * degree);

    //iterate over each point
    for (unsigned i = 0; i < total_points; i += 4) {
        //index and starting coordinate values of x and y for the current point
        __m128i n = _mm_set_epi32(i + 3, i + 2, i + 1, i);
        __m128i x2 = _mm_setzero_si128();
        __m128i y2 = _mm_setzero_si128();

        //calculate coordinates for current point
        for (unsigned s = 0; s < degree; s++) {
            //bitwise-or of x2 with the least significant bit of n
            //which we get via (n & 1)
            //this is then left-shifted by s positions
            //basically the opposite of interleaving
            //similar approach in z_curve_at2
            __m128i n_and_1 = _mm_and_si128(n, _mm_set1_epi32(1));
            __m128i x2_shifted = _mm_slli_epi32(n_and_1, s);
            x2 = _mm_or_si128(x2, x2_shifted);

            //bitwise-or of y2 with the second least significan bit of n
            //which we get via (n & 2)
            //which we then shift once to the right before again
            //left-shifting by s positions
            //basically the opposite of interleaving
            //similar approach in z_curve_at2
            __m128i n_and_2 = _mm_and_si128(n, _mm_set1_epi32(2));
            __m128i y2_shifted = _mm_slli_epi32(_mm_srli_epi32(n_and_2, 1), s);
            y2 = _mm_or_si128(y2, y2_shifted);

            //discard the two bits used for x2 and y2 by shifting twice to the right
            n = _mm_srli_epi32(n, 2);
        }

        //store the calculated coordinates
        _mm_storeu_si128((__m128i*)(x + i), x2);
        _mm_storeu_si128((__m128i*)(y + i), y2);
    }

}

void z_curve_magic(unsigned degree, coord_t* x, coord_t* y) {
    unsigned total_points = 1 << (2 * degree);

    for (unsigned i = 0; i < total_points; i ++) {
        x[i] = ReversePart1By1(i & 0x55555555);
        y[i] = ReversePart1By1((i >> 1) & 0x55555555);
    }
}

void z_curve_magic_simd(unsigned degree, coord_t* x, coord_t* y) { //not faster it seems
    unsigned total_points = 1 << (2 * degree);

    __m128i indices;
    for (unsigned i = 0; i < total_points; i += 4) {
        indices = _mm_setr_epi32(i, i + 1, i + 2, i + 3);

        __m128i x_values = ReversePart1By1_simd(_mm_and_si128(indices, _mm_set1_epi32(0x55555555)));
        __m128i y_values = ReversePart1By1_simd(_mm_and_si128(_mm_srli_epi32(indices, 1), _mm_set1_epi32(0x55555555)));

        _mm_storeu_si128((__m128i*)&x[i], x_values);
        _mm_storeu_si128((__m128i*)&y[i], y_values);
    }
}

//2nd Methods
void z_curve_at(unsigned degree, size_t idx, coord_t* x, coord_t* y) {
    // Calculate the number of points based on the degree
    unsigned numberOfPoints = 1 << (2 * degree);

//    if ( idx<0 || idx>degree) {
//        printf("Invalid index: %s\n", "ERROR");
//        return;
//    }

    // Check if the given index is valid
    if (idx >= numberOfPoints || idx < 0) {
        printf("Invalid index: %s\n", "ERROR");
        return;
    }

    // Starting coordinates
    coord_t start_x = 0;
    coord_t start_y = 0;

    // Loop through the levels of the Z-curve
    for (unsigned d = 0; d < degree; d++) {
        unsigned sub_size = 1 << (degree - d - 1);

        // Calculate the quadrant and offset for the current level
        unsigned quadrant = idx >> (2 * (degree - d - 1));
        coord_t offset_x = (quadrant & 1) ? sub_size : 0;
        coord_t offset_y = (quadrant & 2) ? sub_size : 0;

        // Update the starting coordinates
        start_x += offset_x;
        start_y += offset_y;

        // Update the index for the next level
        idx %= (1 << (2 * (degree - d - 1)));
    }


    // Store the final coordinates
    *x = start_x;
    *y = start_y;

}

void z_curve_at_iterative(unsigned degree, size_t idx, coord_t* x, coord_t* y) {
    //opposite of interleaving similar approach in the z_curve_iterative method (maybe one is faster than the other?)
    coord_t num1 = 0; // Stores the bits of the first number
    coord_t num2 = 0; // Stores the bits of the second number
    coord_t mask = 1; // Bit mask to extract alternating bits

    while (idx > 0) {
        if (idx & 1) {
            num1 |= mask;
            idx >>= 1;
        } else {
            idx >>= 1;
        }
        if (idx & 1) {
            num2 |= mask;
            idx >>= 1;
        } else {
            idx >>= 1;
        }
        mask <<= 1;
    }

    *x = num1;
    *y = num2;
}

void z_curve_at_magic(unsigned degree, size_t idx, coord_t* x, coord_t* y) {
    //needs testing
    *x = ReversePart1By1(idx & 0x55555555);
    *y = ReversePart1By1((idx >> 1) & 0x55555555);
}

//3rd Methods
size_t z_curve_pos(unsigned degree, coord_t x, coord_t y) {
    // Calculate the number of points based on the degree
    unsigned numberOfPoints =  1 << (2 * degree);
    unsigned biggestCoordinate = (1 << degree) -1;

    // Allocate space for the coordinates
    coord_t* xCoords = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
    coord_t* yCoords = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));

    // Generate the Z-curve coordinates
    z_curve_recursive(degree, xCoords, yCoords);

    // Find the index corresponding to the given (x, y) coordinates
    size_t index = 0;
    for (size_t i = 0; i < numberOfPoints; i++) {
        if (xCoords[i] == x && yCoords[i] == y) {
            index = i;
            break;
        }
    }
    if (biggestCoordinate < x || biggestCoordinate < y  || x<0 || y <0) {
        printf("Invalid coordinates: %s\n","ERROR");

    }

    // Free allocated memory
    free(xCoords);
    free(yCoords);

    return index;
}

size_t z_curve_pos_iterative(unsigned degree, coord_t x, coord_t y) {
    //interleaving
    size_t result = 0;
    coord_t mask = 1;

    while (x > 0 || y > 0) {
        if (y & 1) {
            result |= (mask << 1);
        }
        if (x & 1) {
            result |= mask;
        }
        x >>= 1;
        y >>= 1;
        mask <<= 2;
    }

    return result;
}

size_t z_curve_pos_magic(unsigned degree, coord_t x, coord_t y) {
    //needs testing
    degree = (Part1By1(y) << 1) + Part1By1(x);
    return degree;
}

//private testing
void tester3(){
    
    for (int degree; degree <16; degree ++) {
        unsigned numberOfPoints = 1 << (2 * degree);
        coord_t* x1 = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
        coord_t* y1 = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
        coord_t* x2 = (coord_t*)malloc(1 * sizeof(coord_t));
        coord_t* y2 = (coord_t*)malloc(1 * sizeof(coord_t));
        z_curve_magic_simd(degree, x1, y1);
        int atCounter = 0;
        int posCounter = 0;
        for (size_t i = 0; i < numberOfPoints; i++) {
            z_curve_at_magic(degree, i, x2, y2);
            if (x1[i] != x2[0]){
                atCounter++;
            }
            if (y1[i] != y2[0]){
                atCounter++;
            }

            if (i != z_curve_pos_magic(degree, x1[i], y1[i])) {
                posCounter++;
            }
        }
         printf("Degree: %d\natCounter: %d\nposCounter: %d\n", degree, atCounter, posCounter);
    }
}

//Rahmenprogram etc.
void print_usage() {
    // Output a description of all program options and usage examples
    printf("Usage:\n");
    printf("  ./z_curve_program -V <number> -B <number> -d <number> -i <number> -p -a -t -c -h [--help] <number> <number>\n");
    printf("Options:\n");
    printf("  -V <number>     Specify the implementation to use (default: 0)\n");
    printf("  -B <number>     Measure and output the runtime of the specified implementation\n");
    printf("  -d <number>     Specify the degree of the Z-curve to construct\n");
    printf("  -i <number>     Specify the index\n");
    printf("  -p              Call z_curve_pos function\n");
    printf("  -a              Call z_curve_at function\n");
    printf("  -t              Test Method 1\n");
    printf("  -c              Test Method 2\n");
    printf("  -h              Display help\n");
    printf("  --help          Display help\n");
    printf("Positional Arguments:\n");
    printf("  <number>        x\n");
    printf("  <number>        y\n");
}

////GERMAN VERSION
//void print_usage() {
//    printf("Verwendung:\n");
//    printf("  ./z_curve_program -V <Zahl> -B <Zahl> -d <Zahl> -i <Zahl> -p -a -t -c -h [--hilfe] <Zahl> <Zahl>\n");
//    printf("Optionen:\n");
//    printf("  -V <Zahl>       Verwendete Implementierung angeben (Standard: 0)\n");
//    printf("  -B <Zahl>       Laufzeit der angegebenen Implementierung messen und ausgeben\n");
//    printf("  -d <Zahl>       Grad der zu konstruierenden Z-Kurve angeben\n");
//    printf("  -i <Zahl>       idx angeben\n");
//    printf("  -p              z_curve_pos-Funktion aufrufen\n");
//    printf("  -a              z_curve_at-Funktion aufrufen\n");
//    printf("  -t              Methode 1 testen\n");
//    printf("  -c              Methode 2 testen\n");
//    printf("  -h              Hilfe anzeigen\n");
//    printf("  --hilfe         Hilfe anzeigen\n");
//    printf("Positionale Argumente:\n");
//    printf("  <Zahl>          x\n");
//    printf("  <Zahl>          y\n");
//}

void testRun_Method2() {
    unsigned degree = 9;
    unsigned numberOfPoints = 1 << (2 * degree);
    coord_t x_1, y_1;
    ////METHOD @1 AT
    coord_t* x1 = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
    if (x1 == NULL) {
        return ;
    }

    coord_t* y1 = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
    if (y1 == NULL) {
        return;
    }
    for (unsigned i = 0; i <numberOfPoints; ++i) {
        z_curve_at(degree, i, x1, y1);
    }

    ////METHOD @2
    coord_t x_3, y_3;
    coord_t* x3 = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
    if (x3 == NULL) {
        return ;
    }

    coord_t* y3 = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
    if (y3 == NULL) {
        return;
    }
    for (int i = 0; i <numberOfPoints; ++i) {
        z_curve_at_iterative(degree, i, x3, y3);
    }
    int mismatchCount = 0;
    unsigned mismatchIndices[numberOfPoints];
    for (unsigned i = 0; i < numberOfPoints; i++) {
        if (x1[i] != x3[i] || y1[i] != y3[i]) {
            mismatchIndices[mismatchCount] = i;
            mismatchCount++;
        }
    }

    printf("CHECKING METHOD @1 AT AND METHOD @2 AT_ITERATIVE FOR Correctness\n");
    printf("Number of discrepancies :%d\n",mismatchCount);

    ////CHECKING METHOD @1 AT AND METHOD @2 AT_ITERATIVE
    if (mismatchCount > 0) {
        printf("Coordinates are different at the following indices:\n");
    } else {
        printf("Coordinates are identical for x and y arrays.\n");
    }

    printf("\n");
    printf("\n");


    ////METHOD @3

    coord_t* x4 = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
    if (x4 == NULL) {
        return ;
    }

    coord_t* y4 = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
    if (y4 == NULL) {
        return;
    }
    for (unsigned i = 0; i <numberOfPoints; ++i) {
      z_curve_at_magic(degree, i, x4, y4);
    }
    mismatchCount = 0;
    mismatchIndices[numberOfPoints];
    for (unsigned i = 0; i < numberOfPoints; i++) {
        if (x1[i] != x4[i] || y1[i] != y4[i]) {
            mismatchIndices[mismatchCount] = i;
            mismatchCount++;
        }
    }

    printf("CHECKING METHOD @1 AT AND METHOD @3 AT_MAGIC FOR Correctness\n");
    printf("Number of discrepancies :%d\n",mismatchCount);

    ////CHECKING METHOD @1 AT AND METHOD @3 AT_MAGIC
    if (mismatchCount > 0) {
        printf("Coordinates are different at the following indices:\n");
    } else {
        printf("Coordinates are identical for x and y arrays.\n");
    }

    printf("\n");
    printf("\n");

    free(x1);
    free(y1);
    free(x3);
    free(y3);
    free(x4);
    free(y4);
}

void testRun_Method1(){
    unsigned degree =10;

    unsigned numberOfPoints = 1 << (2 * degree);

    ////METHOD @1 RECURSIV
    coord_t* x1 = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
    if (x1 == NULL) {
        return ;
    }

    coord_t* y1 = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
    if (y1 == NULL) {
        return;
    }
    z_curve_magic(degree,x1,y1);


    ////METHOD @2 ITERATIVE
    coord_t* x2 = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
    if (x2 == NULL) {
        return ;
    }

    coord_t* y2 = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
    if (y2 == NULL) {
        return;
    }
    z_curve_iterative(degree,x2,y2);




    int mismatchCount = 0;
    unsigned mismatchIndices[numberOfPoints];
    for (unsigned i = 0; i < numberOfPoints; i++) {
        if (x1[i] != x2[i] || y1[i] != y2[i]) {
            mismatchIndices[mismatchCount] = i;
            mismatchCount++;
        }
    }

    printf("CHECKING METHOD @1 RECURSIVE AND METHOD @2 ITERATIVE FOR Correctness\n");
    printf("Number of discrepancies :%d\n",mismatchCount);

    ////CHECKING METHOD @1 RECURSIVE AND METHOD @2 ITERATIVE
    if (mismatchCount > 0) {
        printf("Coordinates are different at the following indices:\n");
        for (int i = 0; i < mismatchCount; i++) {


            printf("Index %u: x1[%u] = %d, y1[%u] = %d | x2[%u] = %d, y2[%u] = %d\n",
                   mismatchIndices[i], mismatchIndices[i], x1[mismatchIndices[i]], mismatchIndices[i],
                   y1[mismatchIndices[i]], mismatchIndices[i], x2[mismatchIndices[i]], mismatchIndices[i],
                   y2[mismatchIndices[i]]);
        }

    } else {
        printf("Coordinates are identical for x and y arrays.\n");
    }

    printf("\n");
    printf("\n");


    ////METHOD @3 MAGIC
    coord_t* x3 = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
    if (x3 == NULL) {
        return ;
    }

    coord_t* y3 = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
    if (y3 == NULL) {
        return;
    }
    z_curve_magic(degree,x3,y3);



    int mismatchCount1 = 0;
    for (unsigned i = 0; i < numberOfPoints; i++) {
        mismatchIndices[i] = -1;
    }
    mismatchIndices[numberOfPoints];
    for (unsigned i = 0; i < numberOfPoints; i++) {
        if (x1[i] != x3[i] || y1[i] != y3[i]) {
            mismatchIndices[mismatchCount1] = i;
            mismatchCount1++;
        }
    }

    printf("CHECKING METHOD @1 RECURSIVE AND METHOD @3 MAGIC FOR Correctness\n");
    printf("Number of discrepancies :%d\n",mismatchCount1);

    ////CHECKING METHOD @1 RECURSIVE AND METHOD @3 MAGIC
    if (mismatchCount1 > 0) {
        printf("Coordinates are different at the following indices:\n");
        for (int i = 0; i < mismatchCount1; i++) {


            printf("Index %u: x1[%u] = %d, y1[%u] = %d | x3[%u] = %d, y3[%u] = %d\n",
                   mismatchIndices[i], mismatchIndices[i], x1[mismatchIndices[i]], mismatchIndices[i],
                   y1[mismatchIndices[i]], mismatchIndices[i], x3[mismatchIndices[i]], mismatchIndices[i],
                   y3[mismatchIndices[i]]);
        }

    } else {
        printf("Coordinates are identical for x and y arrays.\n");
    }

    printf("\n");
    printf("\n");
    ////METHOD @4 ITERATIVE SIMD
    coord_t* x4 = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
    if (x4 == NULL) {
        return ;
    }

    coord_t* y4 = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
    if (y4 == NULL) {
        return;
    }
    z_curve_iterative_simd(degree,x4,y4);


    int mismatchCount2 = 0;
    for (unsigned i = 0; i < numberOfPoints; i++) {
        mismatchIndices[i] = -1;
    }
     mismatchIndices[numberOfPoints];
    for (unsigned i = 0; i < numberOfPoints; i++) {
        if (x1[i] != x4[i] || y1[i] != y4[i]) {
            mismatchIndices[mismatchCount2] = i;
            mismatchCount2++;
        }
    }

    printf("CHECKING METHOD @1 RECURSIVE AND METHOD @4 ITERATIVE SIMD FOR Correctness\n");
    printf("Number of discrepancies :%d\n",mismatchCount2);

    ////CHECKING METHOD @1 RECURSIVE AND METHOD @4 ITERATIVE SIMD
    if (mismatchCount2 > 0) {
        printf("Coordinates are different at the following indices:\n");
        for (int i = 0; i < mismatchCount2; i++) {


            printf("Index %u: x1[%u] = %d, y1[%u] = %d | x4[%u] = %d, y4[%u] = %d\n",
                   mismatchIndices[i], mismatchIndices[i], x1[mismatchIndices[i]], mismatchIndices[i],
                   y1[mismatchIndices[i]], mismatchIndices[i], x4[mismatchIndices[i]], mismatchIndices[i],
                   y4[mismatchIndices[i]]);
        }

    } else {
        printf("Coordinates are identical for x and y arrays.\n");
    }

    printf("\n");
    printf("\n");
    ////METHOD @4 MAGIC SIMD
    coord_t* x5 = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
    if (x5 == NULL) {
        return ;
    }

    coord_t* y5 = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
    if (y5 == NULL) {
        return;
    }
      z_curve_magic_simd(degree,x5,y5);


    int mismatchCount3 = 0;
    for (unsigned i = 0; i < numberOfPoints; i++) {
        mismatchIndices[i] = -1;
    }
     mismatchIndices[numberOfPoints];
    for (unsigned i = 0; i < numberOfPoints; i++) {
        if (x1[i] != x5[i] || y1[i] != y5[i]) {
            mismatchIndices[mismatchCount3] = i;
            mismatchCount3++;
        }
    }

    printf("CHECKING METHOD @1 RECURSIVE AND METHOD @5 MAGIC SIMD FOR Correctness\n");
    printf("Number of discrepancies :%d\n",mismatchCount3);

    ////CHECKING METHOD @1 RECURSIVE AND METHOD @5 MAGIC SIMD
    if (mismatchCount3 > 0) {
        printf("Coordinates are different at the following indices:\n");
        for (int i = 0; i < mismatchCount3; i++) {


            printf("Index %u: x1[%u] = %d, y1[%u] = %d | x5[%u] = %d, y5[%u] = %d\n",
                   mismatchIndices[i], mismatchIndices[i], x1[mismatchIndices[i]], mismatchIndices[i],
                   y1[mismatchIndices[i]], mismatchIndices[i], x5[mismatchIndices[i]], mismatchIndices[i],
                   y5[mismatchIndices[i]]);
        }

    } else {
        printf("Coordinates are identical for x and y arrays.\n");
    }


    free(x1);
    free(y1);
    free(x2);
    free(y2);
    free(x3);
    free(y3);
    free(x4);
    free(y4);
    free(x5);
    free(y5);
}

void describe_Vmethod_option() {

            printf("Option -V 0: Uses the z_curve_magic2_simd method\n");


            printf("Option -V 1: Uses the z_curve_iterative_simd method\n");

            printf("Option -V 2: Uses the z_curve_magic2 method\n");

            printf("Option -V 3: Uses the z_curve_iterative method\n");

            printf("Option -V 4: Uses the z_curve method\n");



    }

void rahmenProgram_helper( unsigned implementation ,unsigned repetitions,unsigned degree,coord_t* xCoords,
                           coord_t* yCoords,int measure_runtime){
    for (int i = 0; i <= repetitions; i++) {
        double sum = 0;

            struct timespec start, end;
            float z;
            clock_gettime(CLOCK_MONOTONIC, &start);


        switch (implementation) {
            case 0:
                z_curve_magic_simd(degree, xCoords, yCoords);
                break;
            case 1:
                z_curve_iterative_simd(degree, xCoords, yCoords);
                break;
            case 2:
                z_curve_magic(degree, xCoords, yCoords);
                break;
            case 3:
                z_curve_iterative(degree, xCoords, yCoords);
                break;
            case 4:
                z_curve_recursive(degree, xCoords, yCoords);
                break;
            default:
                fprintf(stderr, "Invalid method argument: \n");
                describe_Vmethod_option();
                break;
        }

            if(measure_runtime) {
                clock_gettime(CLOCK_MONOTONIC, &end);


                double elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
                printf("Runtime: %.9f seconds\n", elapsed_time);


            }

        }

    }

void rahmen_svg(coord_t* xCoords, coord_t* yCoords,unsigned numberOfPoints){
    coord_t scalingFactor = 10;
    coord_t minX = xCoords[0];
    coord_t minY = yCoords[0];
    coord_t maxX = xCoords[0];
    coord_t maxY = yCoords[0];
    for (unsigned i = 1; i < numberOfPoints; i++) {
        if (xCoords[i] < minX) minX = xCoords[i];
        if (xCoords[i] > maxX) maxX = xCoords[i];
        if (yCoords[i] < minY) minY = yCoords[i];
        if (yCoords[i] > maxY) maxY = yCoords[i];
    }
    coord_t width = (maxX - minX + 1) * scalingFactor ;
    coord_t height = (maxY - minY + 1) * scalingFactor;
    FILE* svgFile = fopen("zcurve.svg", "w");
    if (svgFile == NULL) {
        printf("Error opening the SVG file.\n");
        return ;
    }

    // Write the SVG header
    fprintf(svgFile, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\">\n", width, height);

    // Write the lines representing the Z curve
    for (unsigned i = 0; i < numberOfPoints - 1; i++) {
        coord_t adjustedX1 = (xCoords[i] - minX) * scalingFactor;
        coord_t adjustedY1 = (yCoords[i] - minY) * scalingFactor;
        coord_t adjustedX2 = (xCoords[i + 1] - minX) * scalingFactor;
        coord_t adjustedY2 = (yCoords[i + 1] - minY) * scalingFactor;
        fprintf(svgFile, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"black\" />\n",
                adjustedX1, adjustedY1, adjustedX2, adjustedY2);
    }

    // Write the SVG footer
    fprintf(svgFile, "</svg>");

    // Close the SVG file and free the allocated memory
    fclose(svgFile);


}

void addDefaultArgument(char *argv[], int *argc) {
    const char *target = "-B";
    const char *defaultArg = "0";
    for (int i = 0; i < *argc - 1; i++) {
        if (strcmp(argv[i], target) == 0 && argv[i + 1][0] == '-') {
            // Insert default argument after "-B"
            int len = strlen(defaultArg);
            char *newArg = (char *) malloc(len + 1);
            strcpy(newArg, defaultArg);

            // Shift arguments to make room for default argument
            for (int j = *argc; j > i + 1; j--) {
                argv[j] = argv[j - 1];
            }

            // Insert default argument
            argv[i + 1] = newArg;
            (*argc)++;  // Increase argument count
        }
    }
}

void rahmenProgram(int argc, char *argv[]) {
    unsigned implementation = 0;  // Default implementation set to 0
    unsigned repetitions = 0;     // Default repetitions set to 0
    unsigned degree = 2;          // Default degree set to 1
    coord_t x = 0;                // Default x set to 0
    coord_t y = 0;                // Default y set to 0
    int call_z_curve_pos = 0;     // Default call_z_curve_pos set to 0
    unsigned idx = 0;             // Default idx set to 0
    int call_z_curve_at = 0;      // Default call_z_curve_at set to 0
    int measure_runtime = 0;
    int testMethod1 = 0;          // Default testMethod1 set to 0
    int testMethod2 = 0;          // Default testMethod2 set to 0

    addDefaultArgument(argv, &argc);

//    for (int i = 1; i <=argc; i++) {
//        printf("Argument %d: %s\n", i, argv[i]);
//    }
    // Parse command-line arguments
    int opt;
    int option_index = 0;
    char *endptr;  // Pointer to store the conversion endpoint

    while ((opt = getopt_long(argc, argv, "V:B:d:i:patch", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'V':
                implementation = strtoul(optarg, &endptr, 10);
                if (*endptr != '\0') {
                    fprintf(stderr, "Invalid option: %s\n", optarg);
                    print_usage();
                    return;
                }
                break;
            case 'B':
                repetitions = strtoul(optarg, &endptr, 10);
                if (*endptr != '\0') {
                    fprintf(stderr, "Invalid option: %s\n", optarg);
                    print_usage();
                    return;
                }
                measure_runtime = 1; // Set the flag for runtime measurement
                break;
            case 't':
                testMethod1 = 1;
                break;
            case 'c':
                testMethod2 = 1;
                break;
            case 'd':
                degree = strtoul(optarg, &endptr, 10);
                if (*endptr != '\0') {
                    fprintf(stderr, "Invalid option: %s\n", optarg);
                    print_usage();
                    return;
                }
                break;
            case 'p':
                call_z_curve_pos = 1;
                break;
            case 'i':
                idx = strtoul(optarg, &endptr, 10);
                if (*endptr != '\0') {
                    fprintf(stderr, "Invalid option: %s\n", optarg);
                    print_usage();
                    return;
                }
                break;
            case 'a':
                call_z_curve_at = 1;
                break;
            case 'h':
                print_usage();
                return;
            default:
                fprintf(stderr, "Invalid option\n");
                print_usage();
                return;
        }
    }

    // Handle missing options and positional arguments
    // ...

    // Assign positional arguments to x and y
    if (argc - optind >= 1) {
        x = strtol(argv[optind], &endptr, 10);
        if (*endptr != '\0') {
            fprintf(stderr, "Invalid positional argument: %s\n", argv[optind]);
            print_usage();
            return;
        }
    }
        if (argc - optind == 2) {
            y = strtol(argv[optind + 1], &endptr, 10);
            if (*endptr != '\0') {
                fprintf(stderr, "Invalid positional argument: %s\n", argv[optind + 1]);
                print_usage();
                return;
            }
        } else {
            // Assign default value to y
            y = y;
        }
//     else {
//        fprintf(stderr, "Missing positional arguments\n");
//        print_usage();
//        return;
//    }

    if (argc > optind + 2) {
        fprintf(stderr, "Invalid arguments\n");
        print_usage();
        return;
    }

    // unsigned numberOfPoints = pow(2,(2 * degree));
    unsigned numberOfPoints = 1 << (2 * degree);

    ////Implementation Check
    if(implementation>5 || implementation<0){
        fprintf(stderr,"Invalid Implement String: %s\n\n", "ERROR");
        describe_Vmethod_option();
        return;
    }
    ////DEGREE CHECK
    if(degree>16 || degree <0){
        printf("Invalid degree for z_curve method: %s\n", "ERROR");
        fprintf(stderr,"Allowed degrees for z_curve method are 0 to 15: %s\n\n", "ERROR");
        //print_usage();
        return;
    }

     ////Repetitions Check
    if (repetitions < 0) {
        fprintf(stderr, "Invalid repetitions value\n");
        //print_usage();
        return;
    }


    ////2nd method check

    if ((idx >= numberOfPoints && call_z_curve_at == 1) || (idx < 0 && call_z_curve_at == 1)) {
        printf("Invalid index for z_curve_at method: %s\n", "ERROR");
        printf("Index should be between 0 to %u for degree %u\n",numberOfPoints-1,degree);
        fprintf(stderr, "Wrong index arguments\n\n");
        //print_usage();
        return ;
    }

    ////3rd method check
    unsigned checkX = 1 << (degree);
    checkX=checkX-1;
    unsigned checkY = 1 << (degree);
    checkY=checkY-1;
    if(call_z_curve_pos) {
        if (x>checkX || x < 0 || y>checkY || y < 0) {
            printf("Invalid coordinates for z_curve_pos method: %s\n", "ERROR");
            printf("For degree %u :\n x should be between 0 to %u\n and y should be between 0 to %u\n",degree,checkX,checkY);
            fprintf(stderr, "Wrong coordinates arguments\n\n");
            //print_usage();
            return;
        }
    }

    // Call the appropriate functions based on the options


    coord_t* xCoords = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
    coord_t* yCoords = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));


    // Call the appropriate implementation based on the specified option
    rahmenProgram_helper(  implementation , repetitions, degree, xCoords, yCoords,measure_runtime);
    printf("\n");

    if (call_z_curve_pos) {
        ////CAN CHANGE TO ANY POS METHOD
        size_t index = z_curve_pos_magic(degree, x, y);
        printf("Calling z_curve_pos\n");
        printf("Z-curve index: %zu\n", index);
    }
    printf("\n");

    coord_t coordX, coordY;

    if (call_z_curve_at) {
        z_curve_at(degree, idx, &coordX, &coordY);
        printf("Calling z_curve_at\n");
        printf("Z-curve coordinates: (%d, %d)\n", x, y);
    }
     printf("\n");

    if(testMethod1 == 1){
        printf("Test for Method 1\n");
        testRun_Method1();
    }

    if(testMethod2 == 1){
        printf("Test for Method 2\n");
        testRun_Method2();
    }

    rahmen_svg(xCoords, yCoords,numberOfPoints);
    free(xCoords);
    free(yCoords);


}
