#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <queue>

//#define NUM_PER_TASK 2048

typedef int32_t num;

static num NUM_PER_TASK;

struct node {
    std::vector<num>::iterator iterator;
    std::vector<num>::iterator end;

    bool operator<(const node &other) const {
        return *iterator >= *other.iterator;
    }
};

void primeRecursive(num upper_bound, num min_bound, std::vector<num> *primes);

void primeGenerator(num start, num end, const std::vector<num> &primes, std::vector<num> *foundPrimes);

void mergeInto(std::vector<num> *parent, std::vector<std::vector<num>> *children);

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "Two arguments expected" << std::endl;

        return -1;
    }

    num N = std::stoi(argv[1]);
    NUM_PER_TASK = std::stoi(argv[2]);

    std::vector<num> primes;
    primes.reserve(static_cast<num>(std::sqrt(N)));

    primes.push_back(2);
    primes.push_back(3);
    primeRecursive(N, 10, &primes);
    std::cout << primes.size() << std::endl;
//    for (size_t i = 0; i < primes.size(); i++) std::cout << primes[i] << std::endl;
}

void primeRecursive(num upper_bound, num min_bound, std::vector<num> *primes) {
    auto &p = *primes;
    if (upper_bound < min_bound) {
        std::vector<num> foundPrimes;
        primeGenerator(3, upper_bound, p, &foundPrimes);
        for (const auto &num: foundPrimes) {
            p.push_back(num);
        }
    } else {
        num lower_bound = static_cast<num>(std::sqrt(upper_bound));

        primeRecursive(lower_bound, min_bound, primes);

        num size = std::ceil(static_cast<float>(upper_bound - lower_bound) / NUM_PER_TASK);

        std::vector<std::vector<num>> foundPrimes(size);

#pragma omp parallel for default(none) shared(lower_bound, upper_bound, foundPrimes, p, NUM_PER_TASK)

        for (num i = lower_bound + 1; i < upper_bound; i += NUM_PER_TASK) {
            primeGenerator(i, std::min(i + NUM_PER_TASK - 1, upper_bound), p,
                           &foundPrimes[(i - lower_bound) / NUM_PER_TASK]);
        }

        for (const auto &fprime: foundPrimes) {
            for (const auto &prime: fprime) {
                p.push_back(prime);
            }
        }
        size_t k = primes->size();
    }
}

void mergeInto(std::vector<num> *parent, std::vector<std::vector<num>> *children) {
    std::vector<num>::iterator child;
    auto &p = *parent;
    auto &c = *children;

    std::priority_queue<node> q;
    for (num i = 0; i < c.size(); i++) {
        if (c[i].size() > 0) {
            q.push({c[i].begin(), c[i].end()});
        }
    }

    while (!q.empty()) {
        auto node = q.top();
        q.pop();
        p.push_back(*node.iterator);

        node.iterator++;
        if (node.iterator != node.end) {
            q.push({node.iterator, node.end});
        }
    }

    int k = q.size();
}

void primeGenerator(num start, num end, const std::vector<num> &primes, std::vector<num> *foundPrimes) {
    int j;
    int k;
    int n;
    int quo, rem;
    auto &p = *foundPrimes;

    p.reserve(end - start);

    P1:
    n = start % 2 == 0 ? start + 1 : start;
    goto P5;

    P2:
    p.push_back(n);

    P3:
    n = n + 2;

    if (n > end) {
        return;
    }
    P5:
    k = 1;

    P6:
    if (k >= primes.size()) {
        goto P2;
    }

    quo = n / primes[k];
    rem = n % primes[k];

    if (rem == 0) {
        goto P3;
    }

    P7:
    if (quo <= primes[k]) {
        goto P2;
    }

    P8:
    k = k + 1;
    goto P6;
}

void givenGenerator() {

}