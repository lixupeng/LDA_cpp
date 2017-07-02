//
// Created by lixupeng on 17-6-27.
//

#ifndef TOPIC_MODEL_HTREE_H
#define TOPIC_MODEL_HTREE_H

#include <queue>
#include <vector>
#include <algorithm>
using namespace std;

struct HTree {

    struct HNode {
        double psum;
        int left, right, parent;
    };

    vector<HNode> nodes;
	vector<int> item;
    int n, root;

	struct cmp {
		vector<HNode>& nd;
		cmp(vector<HNode>& n) : nd(n) {}
		bool operator() (int a, int b)
		{
			return nd[a].psum > nd[b].psum;
		}
	};

    void init(vector<double>& p) {
        n = p.size();
		root = n + n - 1;
        nodes.resize(n * 2);
        item.resize(n);

        for (int i = 0; i < n; ++i) {
			item[i] = i;
            nodes[i].psum = p[i];
        }
		sort(item.begin(), item.end(), [&p](int i, int j) {return p[i] < p[j]; });

        int i = 0, j = n;

        for (int t = n; t <= root; ++t) {
            nodes[t].psum = 0;
            if (i < n && (j == t || p[item[i]] < nodes[j].psum)) {
                nodes[t].psum += p[item[i]];
                nodes[t].left = item[i];
                nodes[item[i]].parent = t;
                i++;
            }
            else {
                nodes[t].psum += nodes[j].psum;
                nodes[t].left = j;
                nodes[j].parent = t;
                j++;
            }
            if (i < n && (j == t || p[item[i]] < nodes[j].psum)) {
                nodes[t].psum += p[item[i]];
                nodes[t].right = item[i];
                nodes[item[i]].parent = t;
                i++;
            }
            else {
                nodes[t].psum += nodes[j].psum;
                nodes[t].right = j;
                nodes[j].parent = t;
                j++;
            }
        }
    }

    double get(int i) {
        return nodes[i].psum;
    }

    double sum() {
        return nodes[root].psum;
    }

    void update(int i, double val) {
        double delta = val - nodes[i].psum;
        nodes[i].psum = val;
        do {
            i = nodes[i].parent;
            nodes[i].psum += delta;
        }
        while (i != root);
    }
	int total_depth = 0;
	int sample_count = 0;

    int sample(double prob) {
        int i = root;
        while (i >= n) {
            int left = nodes[i].left;
            if (prob < nodes[left].psum) {
                i = left;
            }
            else {
                prob -= nodes[left].psum;
                i = nodes[i].right;
            }
        }
        return i;
    }
};

#endif //TOPIC_MODEL_HTREE_H
