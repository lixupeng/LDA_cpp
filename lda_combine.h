//
// Created by Xupeng Li on 13/06/2017.
//

#ifndef TOPIC_MODEL_LDA_COMBINE_H
#define TOPIC_MODEL_LDA_COMBINE_H

#include "data.h"
#include "random.h"
#include "sparse_counter.h"
#include "atomic_int.h"

class LDA_Combine {
public:
    int n_topics;
    double alpha, beta;
    Corpus& corpus;

    LDA_Combine(int topic_num, double alpha, double beta, Corpus& corpus) :
            n_topics(topic_num), alpha(alpha), beta(beta), corpus(corpus) {

        count_topic.resize(topic_num);
        topics.resize(corpus.n_tokens);
        count_doc_topic.resize(corpus.n_docs * n_topics);
        count_word_topic.resize(corpus.n_words);
        for (int i = 0; i < corpus.n_words; ++i) {
            count_word_topic[i].resize(n_topics);
        }

        is_train.resize(corpus.n_docs);
        ftree_method.resize(corpus.n_words);
        int p = rand[0].next_prime(corpus.n_docs + 1), q = p;
        for (int i = 0; i < corpus.n_docs * 0.8; ++i) {
            q = (q + p) % corpus.n_docs;
            is_train[q] = true;
        }
    }

    void init();
    void iterate(bool test);
    double loglikelihood(bool test);

    void start_test();
    void finish_test();
    void experiment();

private:
    vector<int> topics;
    vector<CountTopic> count_word_topic;
    vector<int> count_doc_topic;
    vector <AtomicInt> count_topic;
    vector<bool> ftree_method, is_train;
    Random rand[32];
    AtomicInt ft, mh;

    double cal_doc_topic(int doc, int topic) {
        return (count_doc_topic[doc * n_topics + topic] + alpha) / (count_topic[topic] + corpus.n_words * beta);
    }
};

#endif //TOPIC_MODEL_LDA_COMBINE_H
