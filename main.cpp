#include <iostream>
#include "data.h"
#include "lda_ftree.h"
#include "lda_combine.h"

Corpus corpus;

void lda_main() {
    srand(time(0));
    string path = "data/nytimes.data";
    corpus.read_corpus(path);

    int K = 1000;
    double alpha = 50.0 / K;
    double beta = 0.01;

    LDA_FTree lda(K, alpha, beta, corpus);
    lda.experiment();
}

int main() {
    lda_main();
}
