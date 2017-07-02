//
// Created by lixupeng on 17-5-7.
//

#ifndef TOPIC_MODEL_ALIAS_TABLE_H
#define TOPIC_MODEL_ALIAS_TABLE_H

#include <cmath>
#include <vector>
#include "random.h"
using std::vector;

struct AliasTable {

    struct Column {
        double prob;
        int less;
        int greater;
    };

    vector<Column> table;
    double prob_sum, prob_per_column;
    int n;

    void init(vector<double> &p) {
        n = p.size();
        table.resize(n);
        prob_sum = 0;

        for (int i = 0; i < n; ++i) {
            prob_sum += p[i];
        }
        prob_per_column = prob_sum / p.size();

        for (int i = 0; i < n; ++i) {
            table[i].less = table[i].greater = i;
            table[i].prob = p[i];
        }
        for (int i = 0, j = 0; i < n; ++i) {
            if (fabs(table[i].prob - prob_per_column) < 1e-10) {
                if (i > j) {
                    Column tmp = table[i];
                    table[i] = table[j];
                    table[j] = tmp;
                }
                ++j;
                continue;
            }
            if (table[i].prob > prob_per_column == table[j].prob > prob_per_column)
                continue;
            if (table[i].prob > prob_per_column) {
                table[j].greater = table[i].greater;
                table[i].prob -= prob_per_column - table[j].prob;
            }
            else {
                double new_prob = table[j].prob + table[i].prob - prob_per_column;
                table[j].prob = table[i].prob;
                table[j].less = table[i].less;
                table[i].prob = new_prob;
                table[i].less = table[i].greater = table[j].greater;
            }
            --i, ++j;
        }
    }

    int sample(Random &R) {
        double p = R.rand_double() * prob_sum;
        int n = (int) (p / prob_per_column);
        p -= n * prob_per_column;
        if (p < table[n].prob)
            return table[n].less;
        else
            return table[n].greater;
    }
};

#endif //TOPIC_MODEL_ALIAS_TABLE_H
