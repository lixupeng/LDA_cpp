//
// Created by lixupeng on 17-6-9.
//

#include <queue>
#include "lda_ftree.h"
#include "ftree.h"
#include <omp.h>
#include <time.h>
#include <sys/time.h>

void LDA_FTree::init() {
    for (int i = 0; i < corpus.n_words; ++i) {
	for (int j = corpus.word_offset[i]; j < corpus.word_offset[i + 1]; ++j) {
	    int doc = corpus.words[j];
	    int topic = rand[0].rand_int(n_topics);
	    topics[j] = topic;
	    count_doc_topic[doc].inc(topic);
	   // if (is_train[doc]) {
		count_word_topic[i * n_topics + topic]++;
		count_topic[topic]++;
	  //  }
	}
    }
}

void LDA_FTree::start_test() {
    for (int i = 0; i < corpus.n_words; ++i) {
	for (int j = corpus.word_offset[i]; j < corpus.word_offset[i + 1]; ++j) {
	    int doc = corpus.words[j];
	   // if (is_train[doc]) continue;
	    int topic = topics[j];
	    count_word_topic[i * n_topics + topic]++;
	    count_topic[topic]++;
	}
    }
}

void LDA_FTree::finish_test() {
    for (int i = 0; i < corpus.n_words; ++i) {
	for (int j = corpus.word_offset[i]; j < corpus.word_offset[i + 1]; ++j) {
	    int doc = corpus.words[j];
	  //  if (is_train[doc]) continue;
	    int topic = topics[j];
	    count_word_topic[i * n_topics + topic]--;
	    count_topic[topic]--;
	}
    }
}

void LDA_FTree::iterate(bool test) {
    vector<double> psum;
    FTree tree(n_topics);

#pragma omp parallel for schedule(dynamic) private(psum), firstprivate(tree)
    for (int word = 0; word < corpus.n_words; ++word) {
	int thread = omp_get_thread_num();
	psum.resize(n_topics);
	for (int i = 0; i < n_topics; ++i) {
	    tree.set(i, cal_word_topic(word, i));
	}
	tree.build();
	for (int i = corpus.word_offset[word]; i < corpus.word_offset[word + 1]; ++i) {
	    int doc = corpus.words[i];
	  //  if (is_train[doc] == test) continue;
	    int topic = topics[i];
	    CountTopic &count_doc = count_doc_topic[doc];
	    count_doc.mtx.lock();

	    count_word_topic[word * n_topics + topic]--;
	    count_doc.dec(topic);
	    count_topic[topic]--;
	    tree.update(topic, cal_word_topic(word, topic));

	    double prob_left = tree.sum() * alpha;
	    double prob_all = prob_left;
	    for (int t = 0, s = count_doc.item.size(); t < s; ++t) {
		CountItem &item = count_doc.item[t];
		double p = item.count * tree.get(item.item);
		prob_all += p;
		psum[t] = p;
		if (t > 0) psum[t] += psum[t - 1];
	    }

	    double prob = rand[thread].rand_double(prob_all);
	    int new_topic;
	    if (prob < prob_left) {
		new_topic = tree.sample(prob / alpha);
	    }
	    else {
		prob -= prob_left;
		int p = lower_bound(psum.begin(), psum.begin() + count_doc.item.size(), prob) - psum.begin();
		new_topic = count_doc.item[p].item;
	    }

	    count_doc.inc(new_topic);
	    count_doc.mtx.unlock();

	    topics[i] = new_topic;
	    count_word_topic[word * n_topics + new_topic]++;
	    count_topic[new_topic]++;
	    tree.update(new_topic, cal_word_topic(word, new_topic));
	}
    }
}

double LDA_FTree::loglikelihood(bool test) {
    double llh = 0;

    for (int doc = 0; doc < corpus.n_docs; ++doc) {
//	if (is_train[doc] == test) continue;
	for (int i = 0, s = count_doc_topic[doc].item.size(); i < s; ++i) {
	    llh += lgamma(count_doc_topic[doc].item[i].count + alpha);
	}
	llh += (n_topics - count_doc_topic[doc].item.size()) * lgamma(alpha);
    }

    if (test) {
	vector<int> count_topic_test(n_topics);
	for (int word = 0; word < corpus.n_words; ++word) {
	    vector<int> count_word(n_topics);
	    for (int j = corpus.word_offset[word]; j < corpus.word_offset[word + 1]; ++j) {
		//if (is_train[corpus.words[j]]) continue;
		count_topic_test[topics[j]]++;
		count_word[topics[j]]++;
	    }
	    for (int t = 0; t < n_topics; ++t) {
		llh += lgamma(count_word[t] + beta);
	    }
	}
	for (int topic = 0; topic < n_topics; ++topic) {
	    llh -= lgamma(count_topic_test[topic] + corpus.n_words * beta);
	}
    }
    else {
	for (int word = 0; word < corpus.n_words; ++word) {
	    for (int t = 0; t < n_topics; ++t) {
		llh += lgamma(count_word_topic[word * n_topics + t] + beta);
	    }
	}
	for (int topic = 0; topic < n_topics; ++topic) {
	    llh -= lgamma(count_topic[topic].load() + corpus.n_words * beta);
	}
    }
    return llh;
}

void LDA_FTree::experiment() {
    char logfile[32];
    sprintf(logfile, "log_LDA_FTree_%d", n_topics);
    std::ofstream fout(logfile);

    struct timeval time_start, time_end;
    double sample_time, total_time = 0;

    init();
    int iter = 0;
    fout << "iter\tsample_time\ttotal_time\ttrain_llh" << std::endl;
    while (true) {
	++iter;
	gettimeofday(&time_start, NULL);
	iterate(false);
	gettimeofday(&time_end, NULL);
	sample_time = (double) ((time_end.tv_sec - time_start.tv_sec) * 1000000 + time_end.tv_usec - time_start.tv_usec) / 1e6;
	total_time += sample_time;
	double train_llh = loglikelihood(false);
/*
	start_test();
	for (int i = 0; i < 5; ++i) {
	    iterate(true);
	}
	double test_llh = loglikelihood(true);
	finish_test();
*/
	std::cout << iter << "\t" << sample_time << "s\t" << total_time << "s\t" << train_llh << std::endl;
	fout << iter << "\t" << sample_time << "\t" << total_time << "\t" << train_llh << std::endl;
	fout.flush();
    }
}
