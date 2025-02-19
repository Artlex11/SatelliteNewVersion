#include "Generators.h"

namespace RandomGenerators {

    double generateGauss(double mean, double stddev) {
        std::mt19937 engine;
        std::seed_seq seeds{
            static_cast<uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count()),
            static_cast<uint64_t>(std::chrono::system_clock::now().time_since_epoch().count()),
            static_cast<uint64_t>(std::hash<std::thread::id>{}(std::this_thread::get_id())),
        };
        engine.seed(seeds);
        std::normal_distribution<double> distribution(mean, stddev);
        return distribution(engine);
    }

    int generateUniform(int min, int max) {
        std::mt19937 engine;
        std::seed_seq seeds{
            static_cast<uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count()),
            static_cast<uint64_t>(std::chrono::system_clock::now().time_since_epoch().count()),
            static_cast<uint64_t>(std::hash<std::thread::id>{}(std::this_thread::get_id())),
        };
        engine.seed(seeds);
        std::uniform_int_distribution<int> distribution(min, max);
        return distribution(engine);
    }

    double generateUniform(double min, double max) {
        std::mt19937 engine;
        std::seed_seq seeds{
            static_cast<uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count()),
            static_cast<uint64_t>(std::chrono::system_clock::now().time_since_epoch().count()),
            static_cast<uint64_t>(std::hash<std::thread::id>{}(std::this_thread::get_id())),
        };
        engine.seed(seeds);
        std::uniform_real_distribution<double> distribution(min, max);
        return distribution(engine);
    }

    int generateUniformFromVector(const std::vector<int>& values) {
        std::mt19937 engine;
        std::seed_seq seeds{
            static_cast<uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count()),
            static_cast<uint64_t>(std::chrono::system_clock::now().time_since_epoch().count()),
            static_cast<uint64_t>(std::hash<std::thread::id>{}(std::this_thread::get_id())),
        };
        engine.seed(seeds);
        std::uniform_int_distribution<size_t> distribution(0, values.size() - 1);
        return values[distribution(engine)];
    }

    double generateUniformFromVector(const std::vector<double>& values) {
        std::mt19937 engine;
        std::seed_seq seeds{
            static_cast<uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count()),
            static_cast<uint64_t>(std::chrono::system_clock::now().time_since_epoch().count()),
            static_cast<uint64_t>(std::hash<std::thread::id>{}(std::this_thread::get_id())),
        };
        engine.seed(seeds);
        std::uniform_int_distribution<size_t> distribution(0, values.size() - 1);
        return values[distribution(engine)];
    }

}
