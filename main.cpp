#include <fstream>
#include <sstream>
#include "maps.h"
#include "FTreeLDA.h"
using namespace std;

const unsigned MAX_LINE_LEN = 1 << 20;  //1MB

LDA_Map maps;
vector<document> docs;

void read(char* path) {

    char* buf = new char[MAX_LINE_LEN];

    ifstream fin(path);
    while (fin.getline(buf, MAX_LINE_LEN)) {
        stringstream ss(buf);

        doc_t doc;
        ss >> doc.id1 >> doc.id2;
        document d = new vector<word_id_t>;
        docs.push_back(d);
        doc_id_t doc_id = docs.size() - 1;
        maps.set_doc_id(doc, doc_id);

        string word;

        while (ss >> word) {
            word_id_t word_id = maps.get_word_id(word);
            d->push_back(word_id);
        }
    }
}

int main(int argc, char* argv[]) {
    char* path = argv[1];
    read(path);
    int M = (int) docs.size();
    int V = maps.max_word_id;
    int K = atoi(argv[2]);
    float alpha = 50.0F / K;
    float beta = 0.01F;
    int I = atoi(argv[3]);
    int taskNum = atoi(argv[4]);

    FTreeLDA* lda = new FTreeLDA(M, V, K, alpha, beta, &docs, taskNum);
    lda->initDk();
    lda->initParallel();

    for (int i = 0; i < I; i++) {
        lda->iteration(i);
    }
    return 0;
}