#pragma once
#include <vector>
#include <cmath>
#include "random.h"
#include <fstream>

class DishonestCasino{

    double inDistribution[2] = {0.9, 0.1};
    double trProbability[2][2] = {{0.9, 0.1}, {0.1, 0.9}};
    double resProbability[2][6] = {{1./6, 1./6, 1./6, 1./6, 1./6, 1./6}, {0.1, 0.1, 0.1, 0.1, 0.1, 0.5}};

public:
    DishonestCasino();
    DishonestCasino(double inDistribution[2], double trProbability[2][2],  double resProbability[2][6]);

    double getInDistribution(int state);
    double getTrProbability(int from, int to);
    double getResProbability(int state, int result);

    /* ------------------------------------------------Далее алгоритмы------------------------------------------------*/

    // Форвард процедура
    double forwardProcedure(std::vector<int> const &obsSequence, int state);
    // Бэкворд процедура
    double backwardProcedure(std::vector<int> const &obsSequence, int state);
    // Вероятность наблюдаемой последовательности (через форвард процедуру)
    double probabilityOfObsSeqF(std::vector<int> const &obsSequence);
    // Вероятность наблюдаемой последовательности (через бэкворд процедуру)
    double probabilityOfObsSeqB(std::vector<int> const &obsSequence);
    // Алгоритм Витерби
    std::vector<int> viterbiAlgorithm(std::vector<int> const &obsSequence);
    // Угадывание параметров HMM через алгоритм Витерби
    DishonestCasino* viterbiGuessing(std::vector<int> &obsSequence);
    // Алгоритм Баума-Велча
    DishonestCasino* baumWelchAlgorithm(PRNG& generator, std::vector<int> &obsSequence);
    // Первые пять шагов из сегментационного алгоритма к-средник
    DishonestCasino* kMediumSegmentationAlgorithm(std::vector<int> &obsSequence, int w = 10);
    // Вывод значений функции правдоподобия в файл
    static void likelihoodFunToFile(std::vector<int> &obsSequence);

/* -------------------------------------------Далее вспомогательные функции-------------------------------------------*/

    // Генерация наблюдаемого символа в зависимости от распределения для кости
    int generateObs(PRNG& generator, int state);
    // Генерация последовательности состояний и последовательности наблюдаемых символов
    std::vector<std::vector<int>> generateSequences(PRNG& generator, unsigned long lengthOfStateSequence);
    // Взятие подвектора
    template<typename T> std::vector<T> subVector(std::vector<T> const &v, unsigned long m, unsigned long n);
};
