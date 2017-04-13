//
// Created by lixupeng on 4/2/17.
//

#ifndef LDA_TYPES_H
#define LDA_TYPES_H

#include <vector>
#include <string>
#include <utility>
using namespace std;

typedef string word_t;
typedef unsigned word_id_t;

typedef pair<string, string> doc_t;
typedef unsigned long doc_id_t;

#define id1 first
#define id2 second

typedef vector<word_id_t> *document;

#endif //LDA_TYPES_H
