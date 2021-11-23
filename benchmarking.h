#ifndef BENCHMARKING_H
#define BENCHMARKING_H

#include <chrono>
#include <string>
#include <utility>
#include <unordered_map>

class benchmarking
{
private:
    struct Category
    {
        std::chrono::time_point<std::chrono::system_clock,std::chrono::duration<long,std::ratio<1,1000000000>>> start_time;
        unsigned long long accumulated_microseconds = 0;
    };

    static std::unordered_map<std::string,Category> categories;

public:
    static void start(std::string category)
    {
        categories[category].start_time = std::chrono::high_resolution_clock::now();
    }

    static void stop(std::string category)
    {
        auto stop_time = std::chrono::high_resolution_clock::now();
        Category &cat = categories.at(category);
        cat.accumulated_microseconds +=
            std::chrono::duration_cast<std::chrono::microseconds>(stop_time-cat.start_time).count();
    }

    static void reset()
    {
        categories.clear();
    }

    static std::vector<std::pair<std::string,unsigned long long>> results()
    {
        std::vector<std::pair<std::string,unsigned long long>> result;
        for(auto const &kvp : categories) result.emplace_back(kvp.first, kvp.second.accumulated_microseconds);
        return result;
    }

};

std::unordered_map<std::string,benchmarking::Category> benchmarking::categories;

#endif // BENCHMARKING_H
