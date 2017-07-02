//
// Created by lixupeng on 17-5-7.
//

#ifndef TOPIC_MODEL_FTREE_H
#define TOPIC_MODEL_FTREE_H

#include <vector>
using std::vector;
using std::min;

class FTree {
private:
    vector<double> value;
    int n;
    int topic_num;

public:
    FTree() {}

    FTree(int size): topic_num(size) {
       resize(size);
    }

    void resize(int size) {
        n = 1;
        while (n < size) n <<= 1;
        value.resize(n * 2);
    }

    void set(int i, double val) {
        value[n + i] = val;
    }

    void build() {
        // build must be called after values set
        for (int i = n - 1; i >= 1; --i) {
            value[i] = value[i + i] + value[i + i + 1];
        }
    }

    double get(int i) {
        return value[n + i];
    }

    double sum() {
        return value[1];
    }

    void update(int i, double val) {
        i += n;
        value[i] = val;
        while (i > 1) {
            i >>= 1;
            value[i] = value[i + i] + value[i + i + 1];
        }
    }

    int sample(double prob) {
        int i = 1;
        while (i < n) {
            if (prob < value[i + i]) {
                i = i + i;
            }
            else {
                prob -= value[i + i];
                i = i + i + 1;
            }
        }
        return min(i - n, topic_num - 1);
    }
};

#endif //TOPIC_MODEL_FTREE_H
