//
// Created by lixupeng on 17-5-7.
//

#ifndef TOPIC_MODEL_RANDOM_H
#define TOPIC_MODEL_RANDOM_H

#include <cmath>
#include <cstdlib>

struct Random {
    unsigned x, y, z;
    double V1, V2, S;
    int phase;

    unsigned next_prime(int n) {
        while (true) {
            bool prime = true;
            for (int i = 2; i * i <= n; ++i) {
                if (n % i == 0) {
                    prime = false;
                    break;
                }
            }
            if (prime) break;
            ++n;
        }
        return n;
    }
    Random() {
        x = next_prime(rand());
        y = next_prime(rand());
        z = next_prime(rand());
		for (int i = 0; i < 10; ++i) {
			next();
		}
    }

    unsigned next() {
        return x = y * x + z;
    }

    double rand_double(double x = 1.0) {
        return x * double(next()) / ~0U;
    }

    int rand_int(int n) {
        return next() % n;
    }

    double rand_norm(double mean = 0, double var = 1) {
        double X;
        if ( phase == 0 ) {
            do {
                double U1 = rand_double();
                double U2 = rand_double();

                V1 = 2 * U1 - 1;
                V2 = 2 * U2 - 1;
                S = V1 * V1 + V2 * V2;
            } while(S >= 1 || S == 0);

            X = V1 * sqrt(-2 * log(S) / S);
        }
        else {
            X = V2 * sqrt(-2 * log(S) / S);
        }
        phase = 1 - phase;
        return X * var + mean;
    }
};

#endif //TOPIC_MODEL_RANDOM_H
