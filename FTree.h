//
// Created by lixupeng on 4/3/17.
//

#ifndef LDA_FTREE_H
#define LDA_FTREE_H

#include <vector>
#include <algorithm>
#include <cstring>
using namespace std;

class FTree {
public:
    float* tree;
    int length;
    int K;

    FTree(int length) {
        int len = nextPowerOfTwo(length);
        tree = new float[2 * len];
        this->length = len;
        this->K = length;
        memset(tree, 0, sizeof(float) * 2 * len);
    }

    FTree(float* p, int size, int length) : FTree(length) {
        build(p, size);
    }

    void build(float* p, int size) {
        int start = min(2 * length - 1, length + size - 1);
        for (int i = start; i > 0; i --)  {
            if (i >= length) {
                tree[i] = p[i - length];
            }
            else {
                tree[i] = tree[i << 1] + tree[(i << 1) + 1];
            }
        }
    }

    void update(int index, float value) {
        int i = index + length;
        float delta = value - tree[i];
        while (i > 0) {
            tree[i] += delta;
            i >>= 1;
        }
    }

    int sample(float u) {
        int i = 1;
        while (i < length) {
            if (u < tree[i << 1]) {
                i <<= 1;
            } else {
                u = u - tree[i << 1];
                i = i * 2 + 1;
            }
        }
        return min(i - length, K - 1);
    }

    static int nextPowerOfTwo(int x) {
        if(x == 0) {
            return 1;
        } else {
            --x;
            x |= x >> 1;
            x |= x >> 2;
            x |= x >> 4;
            x |= x >> 8;
            return (x | x >> 16) + 1;
        }
    }

    float first() {
        return tree[1];
    }

    float get(int index) {
        return tree[index + length];
    }
};

#endif //LDA_FTREE_H
