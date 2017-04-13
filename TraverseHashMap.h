//
// Created by lixupeng on 4/3/17.
//

#ifndef LDA_TRAVERSEHASHMAP_H
#define LDA_TRAVERSEHASHMAP_H

#include <vector>
#include <cstring>
#include <mutex>
using std::vector;
using std::mutex;

class TraverseHashMap {
public:
    short size;
    short* key;
    short n;
    int mask;
    mutex mtx;

    short keySize(int expected) {
        short n = 1;
        while (n < expected) {
            n <<= 1;
        }
        return n;
    }

    TraverseHashMap(int expected) {
        this->n = keySize(expected + 1);
        this->key = new short[n];
        this->size = 0;
        this->mask = n - 1;
        memset(this->key, 0xff, sizeof(short) * n);
    }

    TraverseHashMap() {}

    virtual short get(const short k) { return 0; }
    virtual short get(const int k) { return 0; }
    virtual void put(const short k, const short v) {}
    virtual void rehash() {}
    virtual short dec(short k) { return 0; }
    virtual short dec(int k) { return 0; }
    virtual short inc(short k) { return 0; }
    virtual short inc(int k) { return 0; }
};

class S2BTraverseMap : public TraverseHashMap {
public:
    char* value;
    bool* used;
    char* idx;
    char* poss;

    S2BTraverseMap(char expected) : TraverseHashMap(expected) {
        value = new char[n];
        used = new bool[n];
        idx = new char[n];
        poss = new char[n];
        memset(used, 0, sizeof(bool) * n);
    }

    S2BTraverseMap(int expected) : S2BTraverseMap((char)expected){
    }

    short get(const short k) {
        // The starting point
        int pos = k & mask;

        // There's always an unused entry.
        int cnt = 0;
        while (used[pos]) {
            if (key[pos] == k) {
                return value[pos];
            }
            pos = (pos + 1) & mask;
            cnt ++;

            if (cnt > n) {
                rehash();
                return get(k);
            }
        }
        return 0;
    }

    short get(const int k) {
        return get((short) k);
    }

    void put(const short k, const short v) {
        put(k, (char) v);
    }

    void put(const short k, const char v) {
        if (v == 0)
            return;

        // The starting point
        int pos = k & mask;

        // There's always an unused entry.
        while (used[pos]) {
            if (key[pos] == k) {
                value[pos] = v;
                return;
            }
            pos = (pos + 1) & mask;
        }

        used[pos] = true;
        key[pos] = k;
        value[pos] = v;
        idx[size] = (char) pos;
        poss[(char) pos] = (char) size;
        size ++;
    }

    void rehash() {
        short* kkey = key;
        char* vvalue = value;

        key = new short[n];
        value = new char[n];

        memset(used, 0, sizeof(bool) * n);

        int temp = size;
        size = 0;

        for (int i = 0; i < temp; i ++) {
            short k = kkey[idx[i]];
            char v = vvalue[idx[i]];
            put(k, v);
        }
    }

    short dec(short k) {
        int pos = k & mask;

        while (used[pos]) {
            if (key[pos] == k) {
                value[pos] --;
                if (value[pos] == 0) {
                    size --;
                    idx[poss[pos]] = idx[size];
                    poss[idx[size]] = poss[pos];
                }
                return value[pos];
            }
            pos = (pos + 1) & mask;
        }
        return 0;
    }

    short dec(int k) {
        return dec((short) k);
    }

    short inc(short k) {
        int pos = k & mask;

        int cnt = 0;
        while (used[pos]) {
            if (key[pos] == k) {
                value[pos] ++;
                if (value[pos] == 1) {
                    idx[size] = (char) pos;
                    poss[(char) pos] = (char) size;
                    size ++;
                }

                return value[pos];
            }

            cnt ++;
            if (cnt > n) {
                rehash();
                return inc(k);
            }
            pos = (pos + 1) & mask;
        }

        key[pos] = k;
        value[pos] = 1;
        used[pos] = true;
        idx[size] = (char) pos;
        poss[(char) pos] = (char) size;
        size ++;
        return 1;
    }

    short inc(int k) {
        return inc((short) k);
    }

    short getKey(int idx) {
        return key[this->idx[idx]];
    }

    short getVal(int idx) {
        return value[this->idx[idx]];
    }
};

class S2STraverseMap : public TraverseHashMap {
public:
    short* value;
    bool* used;
    short* idx;
    short* poss;

    S2STraverseMap(short expected) : TraverseHashMap(expected) {
        value = new short[n];
        used  = new bool[n];
        idx   = new short[n];
        poss  = new short[n];
        memset(used, 0, sizeof(bool) * n);
    }

    S2STraverseMap(int expected) : S2STraverseMap((short)expected) {
    }

    short get(const short k) {
        // The starting point
        int pos = k & mask;

        // There's always an unused entry.
        int cnt = 0;
        while (used[pos]) {
            if (key[pos] == k) {
                return value[pos];
            }
            pos = (pos + 1) & mask;
            cnt ++;

            if (cnt > n) {
                rehash();
                return get(k);
            }
        }
        return 0;
    }

    short get(const int k) {
        return get((short) k);
    }

    void put(const int k, const int v) {
        put((short) k, (short) v);
    }

    void put(const short k, const short v) {
        if (v == 0)
            return;

        // The starting point
        int pos = k & mask;

        // There's always an unused entry.
        while (used[pos]) {
            if (key[pos] == k) {
                value[pos] = v;
                return;
            }
            pos = (pos + 1) & mask;
        }

        used[pos] = true;
        key[pos] = k;
        value[pos] = v;
        idx[size] = (short) pos;
        poss[(short) pos] = size;
        size ++;
    }

    void rehash() {

        short* kkey = key;
        short* vvalue = value;

//    print();

        key = new short[n];
        value = new short[n];

        memset(used, 0, sizeof(bool) * n);

        int temp = size;
        size = 0;

        for (int i = 0; i < temp; i ++) {
            short k = kkey[idx[i]];
            short v = vvalue[idx[i]];
            put(k, v);
        }

    }

    short dec(int k) {
        return dec((short) k);
    }

    short dec(short k) {
        int pos = k & mask;

        while (used[pos]) {
            if (key[pos] == k) {
                value[pos] --;
                if (value[pos] == 0) {
                    size --;
                    idx[poss[pos]] = idx[size];
                    poss[idx[size]] = poss[pos];
                }
                return value[pos];
            }

            pos = (pos + 1) & mask;
        }
        return 0;
    }

    short inc(int k) {
        return inc((short) k);
    }

    short inc(short k) {
        int pos = k & mask;

        int cnt = 0;
        while (used[pos]) {
            if (key[pos] == k) {
                value[pos] ++;
                if (value[pos] == 1) {
                    idx[size] = (short) pos;
                    poss[pos] = size;
                    size ++;
                }

                return value[pos];
            }

            cnt ++;
            if (cnt > n) {
                rehash();
                return inc(k);
            }
            pos = (pos + 1) & mask;
        }

        key[pos] = k;
        value[pos] = 1;
        used[pos] = true;
        idx[size] = (short) pos;
        poss[(short) pos] = size;
        size ++;
        return 1;
    }

    short getKey(int idx) {
        return key[this->idx[idx]];
    }

    short getVal(int idx) {
        return value[this->idx[idx]];
    }
};

class S2STraverseArray : public TraverseHashMap {
public:
    short* value;
    short* idx;
    short* poss;

    S2STraverseArray(int expected) : TraverseHashMap(expected) {
        this->n = (short) expected;
        this->value = new short[n];
        this->idx   = new short[n];
        this->poss  = new short[n];
        this->size  = 0;
    }

    short get(const short k) {
        return value[k];
    }

    short get(const int k) {
        return value[k];
    }

    void put(const short k, const short v) {
        value[k] = v;
    }

    void rehash() {
        // Do nothing.
    }

    short dec(short k) {
        value[k] --;
        if (value[k] == 0) {
            size --;
            idx[poss[k]] = idx[size];
            poss[idx[size]] = poss[k];
        }
        return value[k];
    }

    short dec(int k) {
        return dec((short) k);
    }

    short inc(short k) {
        value[k] ++;
        if (value[k] == 1) {
            poss[k] = size;
            idx[size ++] = k;
        }
        return value[k];
    }

    short inc(int k) {
        return inc((short) k);
    }

    short getKey(int idx) {
        return this->idx[idx];
    }

    short getVal(int idx) {
        return value[this->idx[idx]];
    }
};

#endif //LDA_TRAVERSEHASHMAP_H
