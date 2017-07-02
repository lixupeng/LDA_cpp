//
// Created by lixupeng on 17-5-7.
//

#ifndef TOPIC_MODEL_SPARSE_COUNTER_H
#define TOPIC_MODEL_SPARSE_COUNTER_H

#include <vector>
#include <mutex>
using std::vector;
using std::mutex;

struct CountItem {
    int item;
    int count;

    CountItem(int it, int cnt): item(it), count(cnt) {}
    CountItem() {}
};

struct CountTopic {
    vector<CountItem> item;
    vector<int> index;
    mutex mtx;

    bool dec(int idx) {
        if (item.size() <= index[idx] || item[index[idx]].item != idx) {
            return false;
        }
        int i = index[idx];
        item[i].count--;
        if (item[i].count == 0) {
            item[i] = *item.rbegin();
            index[item[i].item] = i;
            item.pop_back();
        }
        return true;
    }

    int count(int idx) {
        if (item.size() <= index[idx] || item[index[idx]].item != idx) {
            return 0;
        }
        return item[index[idx]].count;
    }

    void inc(int idx) {
        int i = index[idx];
        if (item.size() <= i || item[i].item != idx) {
            add(idx, 1);
        }
        else {
            item[i].count++;
        }
    }

    void add(int idx, int cnt) {
        index[idx] = item.size();
        item.push_back(CountItem(idx, cnt));
    }

    void resize(int size) {
        item.reserve(size);
        index.resize(size);
    }

    void lock() {
        mtx.lock();
    }

    void unlock() {
        mtx.unlock();
    }

    CountTopic() {};
    CountTopic(int size) { item.reserve(size); index.resize(size); };
    CountTopic(CountTopic &&a) { item = a.item; }
};

#endif //TOPIC_MODEL_SPARSE_COUNTER_H
