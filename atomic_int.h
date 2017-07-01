//
// Created by lixupeng on 17-5-23.
//

#ifndef TOPIC_MODEL_ATOMIC_INT_H
#define TOPIC_MODEL_ATOMIC_INT_H

#include <atomic>
using namespace std;

struct AtomicInt : atomic_int {
    AtomicInt() : atomic_int(0) {}
    AtomicInt(AtomicInt &&b) {
        this->store(b);
    }
};

#endif //TOPIC_MODEL_ATOMIC_INT_H
