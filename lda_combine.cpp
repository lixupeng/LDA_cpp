//
// Created by Xupeng Li on 13/06/2017.
//

#include <queue>
#include "lda_combine.h"
#include "htree.h"
#include <omp.h>
#include <time.h>

void LDA_Combine::init() {
    for (int i = 0; i < corpus.n_docs; ++i) {
        if (is_train[i]) {
            for (int j = corpus.doc_offset[i]; j < corpus.doc_offset[i + 1]; ++j) {
                int word = corpus.docs[j];
                int topic = rand[0].rand_int(n_topics);
                topics[j] = topic;
                count_doc_topic[i * n_topics + topic]++;
                count_word_topic[word].inc(topic);
                count_topic[topic]++;
            }
        }
        else {
            for (int j = corpus.doc_offset[i]; j < corpus.doc_offset[i + 1]; ++j) {
                int topic = rand[0].rand_int(n_topics);
                topics[j] = topic;
                count_doc_topic[i * n_topics + topic]++;
            }
        }
    }
}

void LDA_Combine::start_test() {
    for (int i = 0; i < corpus.n_docs; ++i) {
        if (is_train[i]) continue;
        for (int j = corpus.doc_offset[i]; j < corpus.doc_offset[i + 1]; ++j) {
            int word = corpus.docs[j];
            int topic = topics[j];
            count_word_topic[word].inc(topic);
            count_topic[topic]++;
        }
    }
}

void LDA_Combine::finish_test() {
    for (int i = 0; i < corpus.n_docs; ++i) {
        if (is_train[i]) continue;
        for (int j = corpus.doc_offset[i]; j < corpus.doc_offset[i + 1]; ++j) {
            int word = corpus.docs[j];
            int topic = topics[j];
            count_word_topic[word].dec(topic);
            count_topic[topic]--;
        }
    }
}

void LDA_Combine::iterate(bool test) {
    static HTree tree;
    static vector<double> psum;

#pragma omp parallel for schedule(dynamic) private(psum, tree)
    for (int doc = 0; doc < corpus.n_docs; ++doc) {
        if (is_train[doc] == test) continue;
        int thread = omp_get_thread_num();
        psum.resize(n_topics);
        for (int i = 0; i < n_topics; ++i) {
            psum[i] = cal_doc_topic(doc, i);
        }
        tree.init(psum);
        for (int i = corpus.doc_offset[doc]; i < corpus.doc_offset[doc + 1]; ++i) {
            int word = corpus.docs[i];
            int topic = topics[i];
            int new_topic;

            CountTopic &count_word = count_word_topic[word];
            count_word.mtx.lock();

            count_doc_topic[doc * n_topics + topic]--;
            count_word_topic[word].dec(topic);
            count_topic[topic]--;

            tree.update(topic, cal_doc_topic(doc, topic));

            if (count_word.item.size() < n_topics + 1) {
                double prob_left = tree.sum() * beta;
                double prob_all = prob_left;
                for (int t = 0, s = (int) count_word.item.size(); t < s; ++t) {
                    CountItem &item = count_word.item[t];
                    double p = item.count * tree.get(item.item);
                    prob_all += p;
                    psum[t] = p;
                    if (t > 0) psum[t] += psum[t - 1];
                }

                double prob = rand[thread].rand_double(prob_all);
                if (prob < prob_left) {
                    new_topic = tree.sample(prob / beta);
                } else {
                    prob -= prob_left;
                    int p = lower_bound(psum.begin(), psum.begin() + count_word.item.size(), prob) - psum.begin();
                    new_topic = count_word.item[p].item;
                }
            }
			else {
				for (int c = 0; c < 2; ++c) {
					double sample = rand[thread].rand_double(tree.sum());
					double accept = rand[thread].rand_double();
					new_topic = tree.sample(sample);
					if (accept > (count_word.count(new_topic) + beta) / (count_word.count(topic) + beta)) {
						new_topic = topic;
					}
				}
			}
             
            count_word_topic[word].inc(new_topic);
            count_word.mtx.unlock();

            topics[i] = new_topic;
            count_doc_topic[doc * n_topics + new_topic]++;
            count_topic[new_topic]++;
            tree.update(new_topic, cal_doc_topic(doc, new_topic));
        }
		//cout << double(tree.total_depth) / tree.sample_count << endl;
    }
}

double LDA_Combine::loglikelihood(bool test) {
    double llh = 0;

    for (int doc = 0; doc < corpus.n_docs; ++doc) {
        if (is_train[doc] == test) continue;
        for (int topic = 0; topic < n_topics; ++topic) {
            llh += lgamma(count_doc_topic[doc * n_topics + topic] + alpha);
        }
    }

    for (int doc = 0; doc < corpus.n_docs; ++doc) {
        if (test && is_train[doc]) {
            for (int i = corpus.doc_offset[doc]; i < corpus.doc_offset[doc + 1]; ++i) {
                int word = corpus.docs[i];
                int topic = topics[i];
                count_word_topic[word].dec(topic);
                count_topic[topic]--;
            }
        }
    }

    for (int word = 0; word < corpus.n_words; ++word) {
        for (int t = 0; t < count_word_topic[word].item.size(); ++t) {
            llh += lgamma(count_word_topic[word].item[t].count + beta);
        }
        llh += lgamma(beta) * (n_topics - count_word_topic[word].item.size());
    }
    for (int topic = 0; topic < n_topics; ++topic) {
        llh -= lgamma(count_topic[topic].load() + corpus.n_words * beta);
    }

    for (int doc = 0; doc < corpus.n_docs; ++doc) {
        if (test && is_train[doc]) {
            for (int i = corpus.doc_offset[doc]; i < corpus.doc_offset[doc + 1]; ++i) {
                int word = corpus.docs[i];
                int topic = topics[i];
                count_word_topic[word].inc(topic);
                count_topic[topic]++;
            }
        }
    }
    return llh;
}

void LDA_Combine::experiment() {
    char logfile[32];
    sprintf(logfile, "log_LDA_Combine_%d", n_topics);
    std::ofstream fout(logfile);

	clock_t time_start, time_end;
    double sample_time, total_time = 0;

    init();
    int iter = 0;
    fout << "iter\tsample_time\ttotal_time\ttrain_llh\ttest_llh\ttoken_ft\ttoken_mh" << std::endl;
    while (true) {
        ++iter;
        ft.store(0);
        mh.store(0);

		time_start = clock();
		iterate(false);
		time_end = clock();
		sample_time = (double)(time_end - time_start) / (double)CLOCKS_PER_SEC;
		total_time += sample_time;
		
		double train_llh = 0;// loglikelihood(false);
		/*
        start_test();
        for (int i = 0; i < 5; ++i) {
            iterate(true);
        }
        double test_llh = loglikelihood(true);
        finish_test();
		*/
        std::cout << iter << "\t" << sample_time << "s\t" << total_time << "s\t" << train_llh << "\t" << 0 << "\t" << ft.load() << "\t" << mh.load() << std::endl;
        fout << iter << "\t" << sample_time << "\t" << total_time << "\t" << 0 << "\t" << 0 << "\t" << ft.load() << "\t" << mh.load() << std::endl;
        fout.flush();
    }
}
