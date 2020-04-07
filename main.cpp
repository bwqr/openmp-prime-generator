#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>

//#define NUM_PER_TASK 2048

typedef int32_t num;

static num NUM_PER_TASK;

static num pointer = 0;

void primeRecursive(num upper_bound, num min_bound, std::vector<num> *primes);

void primeGenerator(num start, num end, const std::vector<num> &primes, std::vector<num> *foundPrimes);

num numberOfPrimes(const num &number);

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "An argument expected" << std::endl;

        return -1;
    }

    num N = std::stoi(argv[1]);
    NUM_PER_TASK = 100;
    std::vector<num> primes;
    auto size = static_cast<unsigned long>(numberOfPrimes(N)) * 2;

    if (N < 2) {
        goto exit;
    }
    primes.resize(size);
    primes[pointer++] = 2;

    if (N == 2) {
        goto print;
    }
    primes[pointer++] = 3;
    primeRecursive(N, 10, &primes);
    print:
//    std::cout << primes.size() << std::endl;
    for (size_t i = 0; i < pointer; i++) std::cout << primes[i] << std::endl;
    exit:
    return 0;
}

void primeRecursive(num upper_bound, num min_bound, std::vector<num> *primes) {
    auto &p = *primes;
    if (upper_bound < min_bound) {
        std::vector<num> foundPrimes;
        primeGenerator(3, upper_bound, p, &foundPrimes);
        for (const auto &num: foundPrimes) {
            p[pointer++] = num;
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
            memcpy(&p[pointer], fprime.data(), fprime.size() * sizeof(num));
            pointer += fprime.size();
//            for (const auto &prime: fprime) {
//                p[pointer++] = prime;
//            }
        }
    }
}

num numberOfPrimes(const num &number) {
    return static_cast<num>(static_cast<float>(number) / std::log(number - 1));
}

void primeGenerator(num start, num end, const std::vector<num> &primes, std::vector<num> *foundPrimes) {
    int j;
    int k;
    int n;
    int quo, rem;

    auto &p = *foundPrimes;
    p.reserve(static_cast<unsigned long>(std::max<num>(numberOfPrimes(end) - numberOfPrimes(start), 1)));

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
    if (k >= pointer) {
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
