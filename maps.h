//
// Created by lixupeng on 4/2/17.
//

#ifndef LDA_MAPS_H
#define LDA_MAPS_H

#include "types.h"
#include <map>
using std::map;

class LDA_Map {
public:
    map<word_t, word_id_t> word_to_id;
    map<word_id_t, word_t> id_to_word;
    map<doc_t, doc_id_t> doc_to_id;
    map<doc_id_t, doc_t> id_to_doc;
    word_id_t max_word_id;

    LDA_Map() {}

    word_id_t get_word_id(word_t word) {
        if (word_to_id.count(word) == 0) {
            word_to_id[word] = max_word_id;
            id_to_word[max_word_id] = word;
            max_word_id++;
        }
        return word_to_id[word];
    }

    void set_doc_id(doc_t doc, doc_id_t id) {
        doc_to_id[doc] = id;
        id_to_doc[id] = doc;
    }

    word_t get_word(word_id_t id) {
        return id_to_word[id];
    }

    doc_t get_doc(doc_id_t id) {
        return id_to_doc[id];
    }
};

#endif //LDA_MAPS_H
