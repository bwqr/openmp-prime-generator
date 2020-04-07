#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>

//#define NUM_PER_TASK 2048

typedef int32_t num;

static num NUM_PER_TASK;

void mergeSortedVectors(std::vector<num> *dest, std::vector<std::vector<num>> &sortedVectors);

void primeRecursive(num upper_bound, num min_bound, std::vector<num> *primes);

void primeGenerator(num start, num end, const std::vector<num> &primes, std::vector<std::vector<num>> *foundPrimes);

num numberOfPrimes(const num &number);

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "An argument expected" << std::endl;

        return -1;
    }

    num N = std::stoi(argv[1]);
    NUM_PER_TASK = 1000;
    std::vector<num> primes;

    if (N < 2) {
        goto exit;
    }

    primes.reserve(static_cast<unsigned long>(numberOfPrimes(N)));
    primes.push_back(2);

    if (N == 2) {
        goto print;
    }
    primes.push_back(3);
    primeRecursive(N, 10, &primes);

    print:
//    std::cout << primes.size() << std::endl;
    for (size_t i = 0; i < primes.size(); i++) std::cout << primes[i] << std::endl;
    exit:
    return 0;
}

void primeRecursive(num upper_bound, num min_bound, std::vector<num> *primes) {

    static int num_threads = omp_get_max_threads();

    if (upper_bound < min_bound) {
        std::vector<std::vector<num>> foundPrimes(num_threads);

        primeGenerator(3, upper_bound, *primes, &foundPrimes);

        mergeSortedVectors(primes, foundPrimes);
    } else {
        num lower_bound = static_cast<num>(std::ceil(std::sqrt(upper_bound)));

        primeRecursive(lower_bound, min_bound, primes);

        std::vector<std::vector<num>> foundPrimes(num_threads);

        primeGenerator(lower_bound + 1, upper_bound, *primes, &foundPrimes);

        mergeSortedVectors(primes, foundPrimes);
    }
}

num numberOfPrimes(const num &number) {
    return static_cast<num>(static_cast<float>(number) / std::log(number - 1));
}

void primeGenerator(num start, num end, const std::vector<num> &primes, std::vector<std::vector<num>> *foundPrimes) {
    int j;
    int n;

    for (auto &f: *foundPrimes) {
        f.reserve(static_cast<unsigned long>(std::max<num>(numberOfPrimes(end) - numberOfPrimes(start), 1)) /
                  foundPrimes->size());
    }

    P1:
    n = start % 2 == 0 ? start + 1 : start;

#pragma omp parallel for default(none) shared(foundPrimes, primes, end, n)
    for (num i = n; i <= end; i += 2) {
        int k;
        int quo, rem;

        auto index = omp_get_thread_num();

        goto P5;

        P2:
        (*foundPrimes)[index].push_back(i);

        P3:
//        n = n + 2;
        continue;
//        if (n > end) {
//            return;
//        }
        P5:
        k = 1;

        P6:
        if (k >= primes.size()) {
            goto P2;
        }

        quo = i / primes[k];
        rem = i % primes[k];

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
}

void mergeSortedVectors(std::vector<num> *dest, std::vector<std::vector<num>> &sortedVectors) {
    size_t currentPointerIndex = 0;
    std::vector<std::vector<num>::iterator> pointers;

    for (size_t i = 0, j = 0; i < sortedVectors.size(); i++) {
        if(sortedVectors[i - j].size() > 0) {
            pointers.push_back(sortedVectors[i - j].begin());
        } else {
            sortedVectors.erase(sortedVectors.begin() + i - j);
            j++;
        }
    }

    while (!pointers.empty()) {
        if (pointers[currentPointerIndex] == sortedVectors[currentPointerIndex].end()) {
            if (pointers.size() <= 1) {
                break;
            }

            sortedVectors.erase(sortedVectors.begin() + currentPointerIndex);
            pointers.erase(pointers.begin() + currentPointerIndex);
            currentPointerIndex %= pointers.size();
            continue;
        }

        num min = *pointers[0];
        std::vector<num>::iterator it = pointers[0];
        for (size_t i = 1; i < pointers.size(); i++) {
            if (*pointers[i] < min) {
                min = *pointers[i];
                it = pointers[i];
                currentPointerIndex = i;
            }
        }

        dest->push_back(*pointers[currentPointerIndex]);
        pointers[currentPointerIndex]++;
    }
}