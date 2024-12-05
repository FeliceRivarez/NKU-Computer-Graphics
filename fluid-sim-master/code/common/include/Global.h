#pragma once
#ifndef GLOBAL_H
#define GLOBAL_H
#include <chrono>
#include <random>
#include <unordered_map>
#include <glm/glm.hpp>
#include <string>

namespace Glb {

    class Timer {
    public:
        static Timer& getInstance() {
            static Timer instance;
            return instance;
        }

    private:
        Timer() {

        };

        Timer(const Timer&) = delete;
        Timer& operator=(const Timer&) = delete;

        std::chrono::system_clock::time_point fpsLastTime;
        std::chrono::system_clock::time_point fpsNow;

        std::chrono::system_clock::time_point lastTime;
        std::chrono::system_clock::time_point now;

        std::unordered_map<std::string, unsigned long long int> record;

    public:

        bool empty() {
            return record.empty();
        }

        void clear() {
            record.clear();
        }

        void start() {
            lastTime = std::chrono::system_clock::now();
        } 

        void timeFPS() {
            fpsLastTime = fpsNow;
            fpsNow = std::chrono::system_clock::now();
        }

        std::string getFPS() {
            auto dur = fpsNow - fpsLastTime;
            float dt = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
            if (dt >= 1000) {
                return "0.00";
            }
            else if (dt <= 1) {
                return "1000.00";
            }
            return std::to_string(1000 / dt).substr(0, 5);
        }

        void recordTime(std::string str) {
            now = std::chrono::system_clock::now();
            auto dur = now - lastTime;
            lastTime = now;
            auto it = record.find(str);
            if (it != record.end()) {
                it->second = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
            }
            else {
                record[str] = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
            }
        }

        std::string currentStatus() {
            std::string str;
            int total_time = 0;
            for (const auto& timing : record) {
                total_time += timing.second;
            }

            for (const auto& timing : record) {
                float percentage = static_cast<float>(timing.second) / total_time * 100;
                str += timing.first + ": " + std::to_string(percentage).substr(0, 5) + "%% \n";
            }

            return str;
        }
    };


    class RandomGenerator {
    private:
        std::random_device dev;
    public:
        float GetUniformRandom(float min = 0.0f, float max = 1.0f) {
            std::mt19937 rng(dev());
            std::uniform_real_distribution<float> dist(min, max); // distribution in range [min, max]
            return dist(rng);
        }
    };

}


#endif // !GLOBAL_H




