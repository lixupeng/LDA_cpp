//
// Created by leleyu on 2017/5/3.
//

#ifndef TOPIC_MODEL_DATA_H
#define TOPIC_MODEL_DATA_H

#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <mutex>
#include <cmath>
using namespace std;

struct Corpus {
    // Number of documents
    int n_docs = 0;
    // Number of vocabulary
    int n_words = 0;
    // Number of tokens
    int n_tokens = 0;
    // Document Ids organized by word
    vector<int> words;
    // Word Ids organized by document
    vector<int> docs;
    // index map from docs to words
    vector<int> doc_to_word;
    // index map from words to docs
    vector<int> word_to_doc;
    // Cumulative counts for the number of words
    vector<int> word_offset;
    // Cumulative counts for the number of documents
    vector<int> doc_offset;
    // Vocabulary map
    map<string, int> word_to_int;
    // Document map
    map<string, int> doc_to_int;
    // Vocabulary list
    vector<string> word_list;
    // Document list
    vector<string> doc_list;

    // get doc id
    int get_doc_id(string &doc) {
        if (doc_to_int.count(doc) == 0) {
            doc_list.push_back(doc);
            doc_to_int[doc] = n_docs;
            tmp_docs.push_back(vector<int>());
            tmp_doc_to_word.push_back(vector<int>());
            return n_docs++;
        }
        else {
            return doc_to_int[doc];
        }
    }

    int get_word_id(string &word) {
        if (word_to_int.count(word) == 0) {
            word_list.push_back(word);
            word_to_int[word] = n_words;
            tmp_words.push_back(vector<int>());
            tmp_word_to_doc.push_back(vector<int>());
            return n_words++;
        }
        else {
            return word_to_int[word];
        }
    }

    // Add a token (doc, word) to the corpus
    void add_token(string &doc, string &word) {
        int doc_id = get_doc_id(doc);
        int word_id = get_word_id(word);
        tmp_doc_to_word[doc_id].push_back(tmp_words[word_id].size());
        tmp_word_to_doc[word_id].push_back(tmp_docs[doc_id].size());
        tmp_docs[doc_id].push_back(word_id);
        tmp_words[word_id].push_back(doc_id);
        n_tokens++;
    }

    // call this method when all tokens have been added
    // set up words, docs
    // calculate word_offset, doc_offset
    void done() {
        word_offset.resize(n_words + 1);
        doc_offset.resize(n_docs + 1);
        doc_to_word.reserve(n_tokens);

        word_offset[0] = 0;
        for (int i = 0; i < n_words; ++i) {
            word_offset[i] = words.size();
            words.insert(words.end(), tmp_words[i].begin(), tmp_words[i].end());
        }
        word_offset[n_words] = words.size();

        doc_offset[0] = 0;
        for (int i = 0; i < n_docs; ++i) {
            doc_offset[i] = docs.size();
            docs.insert(docs.end(), tmp_docs[i].begin(), tmp_docs[i].end());
        }
        doc_offset[n_docs] = docs.size();


        for (int i = 0; i < n_words; ++i) {
            for (int j = 0; j < tmp_word_to_doc[i].size(); ++j) {
                word_to_doc.push_back(doc_offset[tmp_words[i][j]] + tmp_word_to_doc[i][j]);
            }
        }
        for (int i = 0; i < n_docs; ++i) {
            for (int j = 0; j < tmp_doc_to_word[i].size(); ++j) {
                doc_to_word.push_back(word_offset[tmp_docs[i][j]] + tmp_doc_to_word[i][j]);
            }
        }
        tmp_words.clear();
        tmp_docs.clear();
        tmp_doc_to_word.clear();
        tmp_word_to_doc.clear();
    }

    void read_corpus(string& doc_path) {
        ifstream fin(doc_path);
        string doc;

        int read_buffer_size = 1 << 20;
        char read_buffer[1 << 20];

        int n = 0;
        // read documents
        while (fin >> doc) {
            // read words
            doc = to_string(n++);
            fin.getline(read_buffer, read_buffer_size, '\n');
            stringstream doc_s(read_buffer);
            string word;
            while (doc_s >> word) {
                add_token(doc, word);
            }
        }
        fin.close();

        done();

        cout << "total words " << n_words << endl;
        cout << "total documents " << n_docs << endl;
        cout << "total tokens " << n_tokens << endl;
    }

private:
    vector<vector<int>> tmp_docs, tmp_words, tmp_doc_to_word, tmp_word_to_doc;
};

#endif //TOPIC_MODEL_DATA_H
