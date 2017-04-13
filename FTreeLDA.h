//
// Created by lixupeng on 4/3/17.
//

#ifndef LDA_FTREELDA_H
#define LDA_FTREELDA_H

#include <boost/math/special_functions/gamma.hpp>
#include <tbb/concurrent_queue.h>
#include <atomic>
#include <cstdlib>
#include "TraverseHashMap.h"
#include "BinarySearch.h"
#include "FTree.h"
#include "types.h"
#include <iostream>
using namespace std;

class FTreeLDA {

public:
    class SampleTask;
    class TestTask;
    class InitTask;

    int D, V, K, N;
    double alpha, beta, vbeta;

    int *ws, *cols, *nk, *docLens;
    short* topics;
    double* results;

    int taskNum;

    vector<TraverseHashMap*> ndk;

    vector<SampleTask*> sampleTasks;
    vector<TestTask*> testTasks;
    vector<InitTask*> initTasks;

    double lgammaBeta;
    double lgammaAlpha;
    double lgammaAlphaSum;

    int step = 5;

    bool test = false;

    void buildMat(vector<document>* docs) {
        document doc;
        int* wcnt = new int[V];
        ws = new int[V + 1];
        memset(wcnt, 0, sizeof(int) * V);

        N = 0;
        // count word
        for (unsigned d = 0; d < D; d ++) {
            doc = (*docs)[d];
            for (int w = 0; w < doc->size(); w++)
                wcnt[(*doc)[w]] ++;
            docLens[d] = (int) doc->size();
            N += doc->size();
        }

        cols = new int[N];

        // build word start idx
        ws[0] = 0;
        for (int i = 0; i < V; i ++) {
            ws[i + 1] = ws[i] + wcnt[i];
        }

        for (int d = D - 1; d >= 0; d --) {
            doc = (*docs)[d];
            int wid;
            for (int w = 0; w < doc->size(); w++) {
                wid = (*doc)[w];
                int pos = ws[wid] + (--wcnt[wid]);
                cols[pos] = d;
            }
        }
        delete[] wcnt;
    }

    FTreeLDA(int D, int V, int K, float alpha, float beta, vector<document>* docs, int taskNum) {
        this->D = D;
        this->V = V;
        this->K = K;
        this->alpha = alpha;
        this->beta  = beta;
        this->vbeta = V * beta;

        docLens = new int[D];

        buildMat(docs);

        topics = new short[N];
        nk = new int[K];
        ndk.resize((unsigned long) D);
        memset(nk, 0, sizeof(int) * K);

        this->taskNum = taskNum;
        sampleTasks.resize((unsigned long) taskNum);
        testTasks.resize((unsigned long) taskNum);
        initTasks.resize((unsigned long) taskNum);
        for (int i = 0; i < taskNum; i ++) {
            sampleTasks[i] = new SampleTask(this);
            testTasks[i] = new TestTask(this);
        }
        results = new double[taskNum];

        lgammaBeta = lgamma(beta);
        lgammaAlpha = lgamma(alpha);
        lgammaAlphaSum = lgamma(alpha * K);
    }

    void initDk() {
        for (int d = 0; d < D; d ++) {
            if (docLens[d] < 127) {
                ndk[d] = new S2BTraverseMap(docLens[d]);
            } else if (docLens[d] < (K / 2)) {
                ndk[d] = new S2STraverseMap(min(K, docLens[d]));
            } else {
                ndk[d] = new S2STraverseArray(K);
            }
        }
    }


    void buildFTree(int* wk, FTree* tree, float* p, int* nk) {
        for (int topic = 0; topic < K; topic ++) {
            p[topic] = (float) ((wk[topic] + beta) / (nk[topic] + vbeta));
        }
        tree->build(p, K);
    }

    float build(S2STraverseMap* map, float* p, short* tidx, FTree* tree) {
        float psum = 0.0F;
        short topic;
        short count;
        int idx = 0;

        for (int i = 0; i < map->size; i++) {
            topic = map->key[map->idx[i]];
            count = map->value[map->idx[i]];
            psum += count * tree->get(topic);
            p[idx] = psum;
            tidx[idx ++] = topic;
        }
        return psum;
    }

    float build(S2BTraverseMap* map, float* p, short* tidx, FTree* tree) {
        float psum = 0.0F;
        short topic;
        short count;
        int idx = 0;

        for (int i = 0; i < map->size; i ++) {
            topic = map->key[map->idx[i]];
            count = map->value[map->idx[i]];
            psum += count * tree->get(topic);
            p[idx] = psum;
            tidx[idx ++] = topic;
        }
        return psum;
    }

    float build(S2STraverseArray* map, float* p, short* tidx, FTree* tree) {
        float psum = 0.0F;
        short topic;
        short count;
        int idx = 0;

        for (int i = 0; i < map->size; i ++) {
            topic = map->idx[i];
            count = map->value[topic];
            psum += count * tree->get(topic);
            p[idx] = psum;
            tidx[idx++] = topic;
        }
        return psum;
    }

    float buildDocDist(int did, float* p, short* tidx, FTree* tree) {
        TraverseHashMap* map = ndk[did];

        if (typeid(*map) == typeid(S2STraverseMap))
            return build((S2STraverseMap*) map, p, tidx, tree);
        if (typeid(*map) == typeid(S2BTraverseMap))
            return build((S2BTraverseMap*) map, p, tidx, tree);
        if (typeid(*map) == typeid(S2STraverseArray))
            return build((S2STraverseArray*) map, p, tidx, tree);
        return 0;
    }

    void init(int wid, int* nk) {
        int si = ws[wid];
        int ei = ws[wid + 1];
        for (int wi = si; wi < ei; wi ++) {
            short kk = (short) (rand() % K);
            topics[wi] = kk;
            nk[kk] ++;
            int d = cols[wi];
            ndk[d]->mtx.lock();
            ndk[d]->inc(kk);
            ndk[d]->mtx.unlock();
        }
    }

    void init() {
        for (int w = 0; w < V; w ++) {
            init(w, NULL);
        }
    }

    void initParallel() {
        tbb::concurrent_queue<int>* queue = new tbb::concurrent_queue<int>();
        for (int w = 0; w < V; w ++)
            queue->push(w);

        //openmp
        for (int i = 0; i < taskNum; i ++) {
            initTasks[i] = new InitTask(this, queue);
            initTasks[i]->call();
        }
        reduceNkInit();
    }

    void reduceNkInit() {
        for (int i = 0; i < taskNum; i ++) {
            int* localNk = initTasks[i]->nk;
            for (int k = 0; k < K; k ++)
                nk[k] += localNk[k];
        }
    }

    long sample(int wid, FTree* tree, float* p, short* tidx, int* nk, int* wk) {
        float value;
        int si, ei;
        si = ws[wid];
        ei = ws[wid + 1];
        int* tokens = cols;
        int d, size, idx;
        short kk;
        float psum, u;
        TraverseHashMap* dk;
        long K_d = 0;

        for (int wi = si; wi < ei; wi ++) {
            d = tokens[wi];
            dk = ndk[d];
            kk = topics[wi];

            wk[kk] --;
            nk[kk] --;
            value = (float) ((wk[kk] + beta) / (nk[kk] + vbeta));
            tree->update(kk, value);

            dk->mtx.lock();
            dk->dec(kk);
            size = dk->size;
            K_d += size;
            psum = buildDocDist(d, p, tidx, tree);

            u = (float) (rand() / float(RAND_MAX) * (psum + alpha * tree->first()));

            if (u < psum) {
                u = rand() / float(RAND_MAX) * psum;
                idx = BinarySearch::binarySearch(p, u, 0, size - 1);
                kk = tidx[idx];
            } else {
                kk = (short) tree->sample(rand() / float(RAND_MAX) * tree->first());
            }

            dk->inc(kk);
            dk->mtx.unlock();

            wk[kk] ++;
            nk[kk] ++;
            value = (float) ((wk[kk] + beta) / (nk[kk] + vbeta));
            tree->update(kk, value);
            topics[wi] = kk;
        }

        return K_d;
    }

    double trainParallel(int it) {

        atomic_int wids(0);
        double ll = 0;

        for (int t = 0; t < taskNum; t ++) {
            sampleTasks[t]->setInteger(&wids, nk);
            results[t] = sampleTasks[t]->call();
        }
        for (int t = 0; t < taskNum; t++) {
            ll += results[t];
        }

        reduceNk();
        long K_d_sum = 0;
        for (int t = 0; t < taskNum; t ++)
            K_d_sum += sampleTasks[t]->K_d;

        double K_d = K_d_sum * 1.0 / N;

        //System.out.format("iteration=%d k_d_sum=%d k_d=%f\n", it, K_d_sum, K_d);

        return ll;
    }

    void reduceNk() {
        for (int i = 0; i < taskNum; i ++) {
            int* localNk = sampleTasks[i]->nk;
            for (int k = 0; k < K; k ++) {
                localNk[k] -= nk[k];
            }
        }

        for (int i = 0; i < taskNum; i ++) {
            int* localNk = sampleTasks[i]->nk;
            for (int k = 0; k < K; k ++)
                nk[k] += localNk[k];
        }
    }

    double loglikelihood() {
        double ll = computeDocLLH();
        ll += computeWordLLHSummary();
        return ll;
    }

    void iteration(int it) {

        test = it % step == 0;
        long start;
        long train_tt, eval_tt;
        //start = System.currentTimeMillis();
        double ll = trainParallel(it);
        //train_tt = System.currentTimeMillis() - start;

        //start = System.currentTimeMillis();
        if (test) {
            ll += loglikelihood();
            cout << ll << endl;
        }

        //eval_tt = System.currentTimeMillis() - start;
        //System.out.format("it=%d ll=%f train_tt=%d eval_tt=%d\n",
        //                 it, ll, train_tt, eval_tt);
    }

    double computeDocLLH() {

        atomic_int dids(0);
        double ll = 0;

        //openmp
        for (int i = 0; i < taskNum; i ++) {
            testTasks[i]->setDids(&dids);
            results[i] = testTasks[i]->call();
        }
        for (int i = 0; i < taskNum; i++) {
            ll += results[i];
        }

        return ll;
    }

    double computeDocLLH(int d, TraverseHashMap* dk) {
        double ll = 0;
        for (int j = 0; j < dk->size; j++) {
            short count = dk->get(j);
            ll += lgamma(alpha + count) - lgammaAlpha;
        }
        ll -= lgamma(alpha * K + docLens[d]) - lgammaAlphaSum;
        return ll;
    }

    double computeWordLLHSummary() {
        double ll = 0.0;
        ll += K * lgamma(beta * V);
        for (int k = 0; k < K; k ++) {
            ll -= lgamma(nk[k] + beta * V);
        }
        return ll;
    }

    class SampleTask {
    public:
        tbb::concurrent_queue<int>* taskQueue;
        float* p;
        short* tidx;
        FTree* tree;
        int* nk;
        int* wk;
        long K_d;
        FTreeLDA* lda;

        atomic_int* wids;

        SampleTask(FTreeLDA* lda) {
            this->lda = lda;
            p = new float[lda->K];
            tidx = new short[lda->K];
            tree = new FTree(lda->K);
            wk = new int[lda->K];
            nk = new int[lda->K];
        }

        void setTaskQueue(tbb::concurrent_queue<int>* taskQueue, int* nk) {
            this->taskQueue = taskQueue;
            memcpy(this->nk, nk, sizeof(int) * lda->K);
            K_d = 0;
        }

        void setInteger(atomic_int* wids, int* nk) {
            this->wids = wids;
            memcpy(this->nk, nk, sizeof(int) * lda->K);
            K_d = 0;
        }

        double call() {
            double ll = 0;

            while (true) {
                int wid = (*wids)++;
                if (wid >= lda->V)
                    break;

                memset(wk, 0, sizeof(int) * lda->K);
                for (int wi = lda->ws[wid]; wi < lda->ws[wid + 1]; wi++) {
                    wk[lda->topics[wi]]++;
                }

                lda->buildFTree(wk, tree, p, nk);
                K_d += lda->sample(wid, tree, p, tidx, nk, wk);

                if (lda->test) {
                    for (int k = 0; k < lda->K; k++)
                        if (wk[k] > 0) {
                            ll += lgamma(wk[k] + lda->beta) - lda->lgammaBeta;
                        }
                }
            }
            return ll;
        }

        int* getNk() {
            return nk;
        }
    };

    class InitTask {
    public:
        tbb::concurrent_queue<int>* taskQueue;
        int* nk;
        FTreeLDA* lda;

        InitTask(FTreeLDA* lda, tbb::concurrent_queue<int>* queue) {
            this->lda = lda;
            taskQueue = queue;
            nk = new int[lda->K];
            memset(nk, 0, sizeof(int) * lda->K);
        }

        double call() {
            int wid;
            while (taskQueue->try_pop(wid)) {
                lda->init(wid, nk);
            }
            return 0.0;
        }
    };

    class TestTask {
    public:
        atomic_int* dids;
        FTreeLDA* lda;

        TestTask(FTreeLDA* lda) {
            this->lda = lda;
        }

        void setDids(atomic_int* dids) {
            this->dids = dids;
        }

        double call() {
            double ll = 0;

            while (true) {
                int did = (*dids)++;
                if (did >= lda->D)
                    break;
                ll += lda->computeDocLLH(did, lda->ndk[did]);
            }
            return ll;
        }
    };
};

#endif //LDA_FTREELDA_H
