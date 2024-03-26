// compilation by: g++ genetic.cpp -O3 -ffast-math -fopenmp
#define _CRT_SECURE_NO_WARNINGS

#include <cmath>
#include <algorithm>
#include <vector>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <chrono>
#include <windows.h>
#include <functional>
#include <sstream>

using namespace std;

LARGE_INTEGER Frequency;

class Stopwatch final {
private:
    LARGE_INTEGER StartingTime{};

public:
    Stopwatch() {
        Reset();
    }

    void Reset() {
        QueryPerformanceCounter(&StartingTime);
    }

    double Elapsed() const {
        LARGE_INTEGER EndingTime;
        QueryPerformanceCounter(&EndingTime);
        auto elapsed = (double)(EndingTime.QuadPart - StartingTime.QuadPart);
        elapsed /= (double)Frequency.QuadPart;
        return elapsed;
    }
};

std::mt19937_64 gen;

const int GENE_BITS = 16;
const double MAX_GENE = (double)(uint16_t)~0;

struct {
    double a, b;
} Interval = { 1, 10 };

struct chromo {
    uint16_t gene;

    void mutate() {
        int mi = (int)gen() % (GENE_BITS - 1);
        gene ^= 1 << mi;
    }
};

double get_x(const chromo& c) {
    double percent = (double)c.gene / MAX_GENE;
    return Interval.a + (Interval.b - Interval.a) * percent;
}

double fitness(const chromo& c) {
    double x = get_x(c);
    double value = log(x) * cos(3*x - 15);
    return value;
}

void crossovermy(chromo& c1, chromo& c2, const int at) {
    uint16_t lower_mask = (1 << at) - 1;
    uint16_t upper_mask = ~lower_mask;
    uint16_t c1l = c1.gene & lower_mask;
    uint16_t c2l = c2.gene & lower_mask;
    c1.gene = (c1.gene & upper_mask) | c2l;
    c2.gene = (c2.gene & upper_mask) | c1l;
}

bool cmpfun(const std::pair<int, double>& r1, const std::pair<int, double>& r2) {
    return r1.second > r2.second;
}

uint64_t coin(uint64_t crp) { // a cointoss
    //std::cout << "rand: " << rand() % 1000 << "\n";
    if (gen() % 1000 < crp) return 1; // crossover
    else return 0; // mutation
}

void init_population(std::vector<chromo>& pop) {
    for (chromo& i : pop) {
        i.gene = gen();
    }
}


class Logger {
public:
    explicit Logger(const std::string& filename) : filename(filename) {
        file.open(filename,
            std::ios_base::app);
        if (!file.is_open()) {
            std::cerr << "file not open " << filename << std::endl;
        }
    }

    ~Logger() {
        if (file.is_open()) {
            file.close();
        }
    }

    void log(const std::string& message) {
        if (file.is_open()) {
            file << message << std::endl;
        }
        else {
            std::cerr << "file not open " << filename << std::endl;
        }
    }

private:
    std::string filename;
    std::ofstream file;
};

int number_of_tests = 40;
const double MAX_TIME = 25;
const uint64_t GENERATION_LIMIT = 10000; // maximum number of generations

double min_possible_fitness = 100;
double fittest;

typedef struct run_results_t {
    int count = 0;
    double sum_x_error = 0;
    double sum_time = 0;
    int sum_gens = 0;
    int incorrect = 0;
} run_results_t;

run_results_t run_tests(uint64_t crossover_prob, uint64_t mutation_prob, uint64_t population_size, std::function<void(int, double, double)> callback) {
    run_results_t res;
    Logger logger("logfile.txt");
    for (res.count = 0; res.count < number_of_tests; ++res.count) {
        vector<chromo> population(population_size);
        vector<double> fitness_values(population_size);
        Stopwatch clock;

        init_population(population);

        double best_x = 0;
        for (uint64_t p = 0; p < GENERATION_LIMIT; p++, res.sum_gens++) {
            //--------------------mutate--------------------
            for (chromo& i : population) {
                if (coin(mutation_prob) == 1) {
                    i.mutate();
                }
            }

            //--------------------crossover--------------------
            for (int i = 0; i < population_size; i++) {
                if (coin(crossover_prob) == 1) {
                    int ind1 = i; // choosing parents for crossover
                    int ind2 = ind1;
                    while (ind2 == ind1) {
                        ind2 = (int)gen() % population_size;
                    }

                    // choose a crossover strategy here
                    crossovermy(population[ind1], population[ind2], (int)gen() % (GENE_BITS - 1));
                }
            }

            //--------------------fitness recount--------------------
            for (int i = 0; i < population_size; i++) {
                fitness_values[i] = fitness(population[i]);
            }

            //--------------------selection--------------------
            vector<chromo> tmp;
            tmp.reserve(population_size);

            for (int i = 0; i < population_size; i++) {
                int ind1 = gen() % population_size;
                int ind2 = gen() % population_size;
                int winner;
                if (fitness_values[ind1] < fitness_values[ind2]) {
                    winner = ind1;
                }
                else {
                    winner = ind2;
                }
                chromo cp(population[winner]);
                tmp.push_back(cp);
            }
            population = tmp;

            //--------------------fitness recount--------------------
            double best = 0;
            for (int i = 0; i < population_size; i++) {
                fitness_values[i] = fitness(population[i]);
                if (fitness_values[i] < best) {
                    best = fitness_values[i];
                    best_x = get_x(population[i]);
                }
            }
            /*if (callback != nullptr) {
                callback(p, best, best_x);
            }*/
            std::ostringstream ist;
            for (auto gen : population) {
                ist << " " << get_x(gen);
            }
            logger.log(std::to_string(p) + ist.str());
           // std::cout << min_possible_fitness << "\n";
            double error_percent = abs((min_possible_fitness - best) / min_possible_fitness);
            if (p % 1000 == 0 || error_percent < 0.0001) {
                printf("\n#%llu\t", p);
                printf("diff fitness: %f \t", best - min_possible_fitness);
                printf("best fitness: %f \t", best);
                double t = clock.Elapsed();
                if (t > MAX_TIME) {
                    res.incorrect += 1;
                    goto end;
                }
                if (error_percent < 0.0001) goto end; // psst...don't tell anyone
            }
        }

    end:
        res.sum_x_error += abs((fittest - best_x) / fittest);
        double t = clock.Elapsed();
        res.sum_time += t;
        cout << "\nCompletion time: " << t << "s.\n";
    }
    return res;
}

int main() {
    Logger logger("logfile.txt");
    logger.log("//--------------------------------");
    QueryPerformanceFrequency(&Frequency);

    random_device::result_type seed = std::random_device{}();
    //random_device::result_type seed = 1337;
    std::cout << "seed: " << seed << std::endl;
    gen = std::mt19937_64(seed);

    for (int i = 0; i < MAX_GENE + 1; ++i) {
        chromo c;
        c.gene = i;
        double fit = fitness(c);
        if (fit < min_possible_fitness) {
            min_possible_fitness = fit;
            fittest = get_x(c);
        }
    }
    std::cout << min_possible_fitness << " minposs\n";

    std::cout << min_possible_fitness << " " << fittest << " results\n";

    //--------------------generation parameters--------------------
    
    vector<int> mut_probabilities = { 10, 20, 40, 60, 80, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000 };
    vector<int> crossover_probabilities = { 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000 };
    vector<int> population_sizes = { 25, 50, 100, 200, 500, 1000 };

    logger.log("mutation");
    for (auto mut : mut_probabilities) {
        auto res = run_tests(750, mut, 100, nullptr);

        logger.log(
            std::to_string(mut) + " " +
            std::to_string(res.sum_gens / number_of_tests) + " " +
            std::to_string(res.sum_time / number_of_tests) + " " +
            std::to_string(res.sum_x_error / number_of_tests) + " " +
            std::to_string(res.incorrect)
        );
    }

    logger.log("crossover");
    for (auto crossover_probability : crossover_probabilities) {
        auto res = run_tests(crossover_probability, 800, 100, nullptr);

        logger.log(
            std::to_string(crossover_probability) + " " +
            std::to_string(res.sum_gens / number_of_tests) + " " +
            std::to_string(res.sum_time / number_of_tests) + " " +
            std::to_string(res.sum_x_error / number_of_tests) + " " +
            std::to_string(res.incorrect)
        );
    }

    logger.log("population");
    for (auto population_size : population_sizes) {
        auto res = run_tests(750, 800, population_size, nullptr);

        logger.log(
            std::to_string(population_size) + " " +
            std::to_string(res.sum_gens / number_of_tests) + " " +
            std::to_string(res.sum_time / number_of_tests) + " " +
            std::to_string(res.sum_x_error / number_of_tests) + " " +
            std::to_string(res.incorrect)
        );
    }
   
    number_of_tests = 1;
    logger.log("iterations");
    run_tests(750, 800, 30, [&](int iteration, double best, double best_x) {
        logger.log(std::to_string(iteration) + " " + std::to_string(best) + " " + std::to_string(best_x));
        });
    
    logger.log("iterations");
    run_tests(750, 800, 25, [&](int iteration, double best, double best_x) {
        logger.log(std::to_string(iteration) + " " + std::to_string(best) + " " + std::to_string(best_x));
        });
        
    return 0;
}
