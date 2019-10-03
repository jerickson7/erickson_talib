#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>

#include "entropy.h"

// Maximum likelihood entropy estimator
int entropy(int start_idx, int end_idx, int window, int word_length, int encoding_states,
            const int input[], double out[]){

    // Window must be larger than word length. It should be much larger
    if (window <= word_length){
        printf("Entropy - Window must be larger than word length. It should be much larger to get an accurate estimate.\n");
        return 1;
    }


    // word_length ^ encoding states = total number of possible words.
    // dictionary stores the number of occurrences of each word.
    int total_words = pow(encoding_states, word_length);
    int* dictionary = (int*) calloc(total_words, sizeof(int));


    int i, j, n, word_index;

    // 1 / Total occurrences. Multiply the counts in dictionary to get the pmf
    double pmf_mult = 1.0 / ((double) (window - word_length));
    double pmf;


    double entropy;
    for (i = start_idx + window; i <= end_idx; i++){

        // Add words to dictionary
        for (j = 0; j < window - word_length; j++){
            word_index = 0;
            //printf("Word: %i ", j);
            for (n = 0; n < word_length; n++){
                word_index += input[i - j - n] * pow(encoding_states, n);
                //printf("%i ", input[i - j - n]);
            }
            //printf(" Index: %i \n", word_index);
            dictionary[word_index]++; // Increment that word
        }


        // Sum the prob * log2(prob)
        entropy = 0.0;
        for (j = 0; j < total_words; j++) {
            pmf = dictionary[j] * pmf_mult;
            if (pmf > 0) entropy += pmf * log2f(pmf);
        }

        // Multiply by -1 / w
        entropy *= -1;
        entropy /= word_length;
        out[i] = entropy;
        memset(dictionary, 0, total_words * sizeof(*dictionary));
    }

    free(dictionary);
    return 0;
}
