#pragma once
#include <random>
#include <cassert>


struct PRNG{
    std::mt19937 engine;
};
void initGenerator(PRNG& generator);
// Генерация целого число в диапазоне [minValue, maxValue)
unsigned getRandomInt(PRNG& generator, unsigned minValue, unsigned maxValue);
// Генерация числа с плавающей точкой в диапазоне [minValue, maxValue)
double getRandomDouble(PRNG& generator, double minValue, double maxValue);
// Генерация массива нормированных на единицу чисел
double* getNormRandomDouble(PRNG& generator, int n);
