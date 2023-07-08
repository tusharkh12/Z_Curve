#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
// #include <immintrin.h>

typedef int coord_t;

void z_curve_iterative(unsigned degree, coord_t* x, coord_t* y);
void z_curve(unsigned degree, coord_t* x, coord_t* y);
void z_curve_recursive(unsigned degree, coord_t start_x, coord_t start_y, coord_t* x, coord_t* y, unsigned* index);

void z_curve_at(unsigned degree, size_t idx, coord_t* x, coord_t* y);
size_t z_curve_pos(unsigned degree, coord_t x, coord_t y);

void z_curve_at2(unsigned degree, size_t idx, coord_t* x, coord_t* y) {
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

size_t z_curve_pos2(unsigned degree, coord_t x, coord_t y) {
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

size_t z_curve_pos3(unsigned degree, coord_t x, coord_t y) {
    //interleaving
    size_t z = 0;
    for (int i = 0; i < sizeof(unsigned int) * 8; i++) {
        z |= (x & (1u << i)) << i | (y & (1u << i)) << (i + 1);
    }
    return z;
}

// void z_curve_recursive_simd(unsigned degree, unsigned start_index, unsigned sub_size, coord_t* x, coord_t* y, unsigned* index) {
//     if (degree == 0) {
//         unsigned end_index = start_index + sub_size * sub_size;
//         for (unsigned i = start_index; i < end_index && *index < (1 << (2 * degree)); i++) {
//             x[*index] = i % sub_size;
//             y[*index] = i / sub_size;
//             (*index)++;
//         }
//         return;
//     }

//     unsigned num_points = sub_size * sub_size;

//     __m128i start_index_vec = _mm_set1_epi32(start_index);
//     __m128i sub_size_vec = _mm_set1_epi32(sub_size);

//     for (unsigned i = 0; i < num_points; i += 4) {
//         if (*index >= (1 << (2 * degree)))
//             break;

//         __m128i indices = _mm_add_epi32(start_index_vec, _mm_set_epi32(i + 3, i + 2, i + 1, i));
//         __m128i x_indices = _mm_and_si128(indices, _mm_set1_epi32(sub_size - 1));
//         __m128i y_indices = _mm_srli_epi32(indices, 1);

//         _mm_storeu_si128((__m128i*) &x[*index], x_indices);
//         _mm_storeu_si128((__m128i*) &y[*index], y_indices);

//         (*index) += 4;
//     }

//     unsigned next_start_index = start_index;
//     unsigned next_sub_size = sub_size >> 1;

//     z_curve_recursive_simd(degree - 1, next_start_index, next_sub_size, x, y, index); // Top-left quadrant
//     z_curve_recursive_simd(degree - 1, next_start_index + next_sub_size, next_sub_size, x, y, index); // Top-right quadrant
//     z_curve_recursive_simd(degree - 1, next_start_index + next_sub_size * sub_size, next_sub_size, x, y, index); // Bottom-left quadrant
//     z_curve_recursive_simd(degree - 1, next_start_index + next_sub_size * sub_size + next_sub_size, next_sub_size, x, y, index); // Bottom-right quadrant
// }

// void z_curve_simd(unsigned degree, coord_t* x, coord_t* y) {
//     unsigned size = 1 << degree;
//     unsigned index = 0;

//     z_curve_recursive_simd(degree, 0, size, x, y, &index);
// }


void z_curve_morton(unsigned degree, coord_t* x, coord_t* y) {
    // calculate the total number of points
    // 2^(2*degree)
    unsigned total_points = 1 << (2 * degree);

    // iterate over each point
    for (unsigned i = 0; i < total_points; i++) {
        // calculate Morton code (Z-curve)
        unsigned x2 = 0;
        unsigned y2 = 0;

        for (unsigned bit = 0; bit < degree; bit++) {
            x2 |= (i & (1 << (2 * bit))) >> bit;
            y2 |= (i & (1 << (2 * bit + 1))) >> (bit + 1);
        }

        // store the calculated coordinates
        x[i] = x2;
        y[i] = y2;
    }
}


int main(int argc, char *argv[]) {
    //for testing purpose will be replaced with inputs later
    unsigned testDegree = 15;
    //coord_t testX = 0;
    //coord_t testY = 1;
    //3rd method test check
    //size_t index = z_curve_pos(testDegree, testX, testY);
    //printf("3rd Method :Index for (%d, %d): %zu\n", testX, testY, index);

    //calculate number of points based on degree
    unsigned numberOfPoints = 1 << (2 * testDegree);

    //allocate space of coordinates
    coord_t* x = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
    if (x == NULL) {
        return -1;
    }

    coord_t* y = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
    if (y == NULL) {
        return -1;
    }

    //create image (SVG for now)

    // FILE* svgFile = fopen("zcurve.svg", "wb");
    // if (svgFile == NULL) {
    //     printf("Error opening the SVG file.\n");
    //     return -1;
    // }

    //clock setup (could prob be done cleaner)
    clock_t begin, end;
    float z;
    begin = clock();

    z_curve(testDegree, x, y);
    //z_curve_iterative(testDegree, x, y);
    //z_curve_simd(testDegree, x, y);
    // z_curve_morton(testDegree, x, y);

    //2nd method test check
    //coord_t coordX, coordY;
    //size_t index =64;

    //z_curve_at(testDegree, index,  &coordX, &coordY);
    //void z_curve_at(unsigned degree, size_t idx, coord_t* x, coord_t* y)

    // Print the resulting coordinates
    //printf("2nd Method :Coordinates for index %zu on the Z-curve:\n", index);
    //printf("x = %d\n", coordX);
    //printf("y = %d\n", coordY);

    end = clock();
    z = end - begin;
    z /= CLOCKS_PER_SEC;
    if (z != 0) {
        printf("Runtime: %f seconds\n", z);
    }


    // //to determine the width and height of the SVG Image we need to find min/max of x and y
    // coord_t scalingFactor = 10;
    // coord_t minX = x[0];
    // coord_t minY = y[0];
    // coord_t maxX = x[0];
    // coord_t maxY = y[0];
    // for (unsigned i = 1; i < numberOfPoints; i++) {
    //     if (x[i] < minX) minX = x[i];
    //     if (x[i] > maxX) maxX = x[i];
    //     if (y[i] < minY) minY = y[i];
    //     if (y[i] > maxY) maxY = y[i];
    // }

    // //dimensions of the SVG image
    // coord_t width = (maxX - minX + 1) * scalingFactor+100;
    // coord_t height = (maxY - minY + 1) * scalingFactor;

    // fprintf(svgFile, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\">\n", width, height); //header

    // for (unsigned i = 0; i < numberOfPoints; i++) {
    //     //adjust the coordinates by subtracting the minimum values
    //     coord_t adjustedX1 = (x[i] - minX) * scalingFactor;
    //     coord_t adjustedY1 = (y[i] - minY) * scalingFactor;

    //     //check if there is a next point to draw a line
    //     if (i + 1 < numberOfPoints) {
    //         coord_t adjustedX2 = (x[i + 1] - minX) * scalingFactor;
    //         coord_t adjustedY2 = (y[i + 1] - minY) * scalingFactor;

    //         fprintf(svgFile, "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" stroke=\"black\" />\n",
    //                 adjustedX1, adjustedY1, adjustedX2, adjustedY2);
    //     }
    // }


    // fprintf(svgFile, "</svg>"); //footer

    //free allocated memory
    free(x);
    free(y);
    // fclose(svgFile);

    return 0;
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

void z_curve(unsigned degree, coord_t* x, coord_t* y) {
    //Calculate the total number of points
    unsigned size = 1 << (2 * degree);
    //index for buffer
    unsigned index = 0;

    //starting coordinates
    coord_t start_x = 0;
    coord_t start_y = 0;

    z_curve_recursive(degree, start_x, start_y, x, y, &index);
}

void z_curve_recursive(unsigned degree, coord_t start_x, coord_t start_y, coord_t* x, coord_t* y, unsigned* index) {
    if (degree == 0) {
        x[*index] = start_x;
        y[*index] = start_y;
        (*index)++;
        return;
    }

    unsigned sub_size = 1 << (degree - 1);
    //unsigned sub_size_half = sub_size >> 1;

    z_curve_recursive(degree - 1, start_x, start_y, x, y, index);
    z_curve_recursive(degree - 1, start_x + sub_size, start_y, x, y, index);
    z_curve_recursive(degree - 1, start_x, start_y + sub_size, x, y, index);
    z_curve_recursive(degree - 1, start_x + sub_size, start_y + sub_size, x, y, index);
}


//2nd METHOD
void z_curve_at(unsigned degree, size_t idx, coord_t* x, coord_t* y) {
    // Calculate the number of points based on the degree
    unsigned numberOfPoints = pow(2,(2 * degree));

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


//3rd METHOD

size_t z_curve_pos(unsigned degree, coord_t x, coord_t y) {
    // Calculate the number of points based on the degree
    unsigned numberOfPoints = pow(2,(2 * degree));

    // Allocate space for the coordinates
    coord_t* xCoords = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));
    coord_t* yCoords = (coord_t*)malloc(numberOfPoints * sizeof(coord_t));

    // Generate the Z-curve coordinates
    z_curve(degree, xCoords, yCoords);

    // Find the index corresponding to the given (x, y) coordinates
    size_t index = 0;
    for (size_t i = 0; i < numberOfPoints; i++) {
        if (xCoords[i] == x && yCoords[i] == y) {
            index = i;
            break;
        }
    }
    if (degree < x || degree < y  || x<0 || y <0) {
        printf("Invalid coordinates: %s\n","ERROR");

    }

    // Free allocated memory
    free(xCoords);
    free(yCoords);

    return index;
}








