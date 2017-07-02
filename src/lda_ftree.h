//
// Created by lixupeng on 17-6-9.
//

#ifndef TOPIC_MODEL_LDA_FTREE_H
#define TOPIC_MODEL_LDA_FTREE_H

#include "data.h"
#include "random.h"
#include "sparse_counter.h"
#include "atomic_int.h"

class LDA_FTree {
public:
    int n_topics;
    double alpha, beta;
    Corpus& corpus;

    LDA_FTree(int topic_num, double alpha, double beta, Corpus& corpus) :
            n_topics(topic_num), alpha(alpha), beta(beta), corpus(corpus) {

        count_topic.resize(topic_num);
        topics.resize(corpus.n_tokens);
        count_word_topic.resize(corpus.n_words * n_topics);
        count_doc_topic.resize(corpus.n_docs);
        for (int i = 0; i < corpus.n_docs; ++i) {
            count_doc_topic[i].resize(n_topics);
        }
/*
        is_train.resize(corpus.n_docs);
        int p = rand[0].next_prime(corpus.n_docs + 1), q = p;
        for (int i = 0; i < corpus.n_docs * 0.9; ++i) {
            q = (q + p) % corpus.n_docs;
            is_train[q] = true;
        }
	*/
    }

    void init();
    void iterate(bool test);
    double loglikelihood(bool test);

    void start_test();
    void finish_test();
    void experiment();

private:
    vector<int> topics;
    vector<CountTopic> count_doc_topic;
    vector<int> count_word_topic;
    vector <AtomicInt> count_topic;
    //vector<bool> is_train;
    Random rand[64];

    double cal_word_topic(int word, int topic) {
        return (count_word_topic[word * n_topics + topic] + beta) / (count_topic[topic] + corpus.n_words * beta);
    }
};

#endif //TOPIC_MODEL_LDA_FTREE_H
