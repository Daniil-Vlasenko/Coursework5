#include "random.h"


void initGenerator(PRNG& generator){
    // Создаём псевдо-устройство для получения случайного зерна.
    std::random_device device;
    // Получаем случайное зерно последовательности
    generator.engine.seed(device());
}

// Генерация целого число в диапазоне [minValue, maxValue)
unsigned getRandomInt(PRNG& generator, unsigned minValue, unsigned maxValue){
    // Проверяем корректность аргументов
    assert(minValue < maxValue);
    // Создаём распределение
    std::uniform_int_distribution<unsigned> distribution(minValue, maxValue);
    // Вычисляем псевдослучайное число: вызовем распределение как функцию,
    //  передав генератор произвольных целых чисел как аргумент.
    return distribution(generator.engine);
}

// Генерация числа с плавающей точкой в диапазоне [minValue, maxValue)
double getRandomDouble(PRNG& generator, double minValue, double maxValue){
    // Проверяем корректность аргументов
    assert(minValue < maxValue);
    // Создаём распределение
    std::uniform_real_distribution<double> distribution(minValue, maxValue);
    // Вычисляем псевдослучайное число: вызовем распределение как функцию,
    //  передав генератор произвольных целых чисел как аргумент.
    return distribution(generator.engine);
}

// Генерация массива нормированных на единицу чисел
double* getNormRandomDouble(PRNG& generator, int n){
    double *mas = new double[n], sum{};
    for(int i = 0; i < n; ++i){
        mas[i] = getRandomDouble(generator, 0, 1);
        sum += mas[i];
    }
    for(int i = 0; i < n; ++i)
        mas[i] /= sum;
    
    return mas;
}