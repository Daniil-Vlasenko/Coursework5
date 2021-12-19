#include "HMM.h"
#include <iostream>
#include <cmath>
#include <unistd.h>


DishonestCasino::DishonestCasino() = default;
DishonestCasino::DishonestCasino(double inDistribution[2], double trProbability[2][2],  double resProbability[2][6]){
    this->inDistribution[0] = inDistribution[0];
    this->inDistribution[1] = inDistribution[1];
    for(int i = 0; i < 2; ++i)
        for(int j = 0; j < 2; ++j)
            this->trProbability[i][j] = trProbability[i][j];
    for(int i = 0; i < 2; ++i)
        for(int j = 0; j < 6; ++j)
            this->resProbability[i][j] = resProbability[i][j];
}

double DishonestCasino::getInDistribution(int state){return inDistribution[state];}
double DishonestCasino::getTrProbability(int from, int to){ return trProbability[from][to];}
double DishonestCasino::getResProbability(int state, int result){return resProbability[state][result];}

/* --------------------------------------------------Далее алгоритмы--------------------------------------------------*/

// Форвард процедура
double DishonestCasino::forwardProcedure(std::vector<int> const &obsSequence, int state){
    double probability[2] = {getInDistribution(0) * getResProbability(0, obsSequence[0]),
                             getInDistribution(1) * getResProbability(1, obsSequence[0])},
           probabilityNext[2]{};

    unsigned long lengthOfObsSequence = obsSequence.size();
    for (unsigned long i = 1; i < lengthOfObsSequence; ++i){
        probabilityNext[0] = (probability[0] * getTrProbability(0, 0) +
                              probability[1] * getTrProbability(1, 0)) * getResProbability(0,obsSequence[i]);
        probabilityNext[1] = (probability[0] * getTrProbability(0, 1) +
                              probability[1] * getTrProbability(1, 1)) * getResProbability(1,obsSequence[i]);

        probability[0] = probabilityNext[0];
        probability[1] = probabilityNext[1];
    }
    if (state)
        return probability[1];
    else
        return probability[0];
}

// Бэкворд процедура
double DishonestCasino::backwardProcedure(std::vector<int> const &obsSequence, int state) {
    double probability[2]{}, probabilityPrev[2] = {1, 1};

    long lengthOfObsSequence = obsSequence.size();
    for (long i = lengthOfObsSequence - 2; i >= -1; --i) {
        probability[0] = probabilityPrev[0];
        probability[1] = probabilityPrev[1];

        probabilityPrev[0] = getTrProbability(0, 0) * getResProbability(0, obsSequence[i + 1]) * probability[0]
                             + getTrProbability(0, 1) * getResProbability(1, obsSequence[i + 1]) * probability[1];
        probabilityPrev[1] = getTrProbability(1, 0) * getResProbability(0, obsSequence[i + 1]) * probability[0]
                             + getTrProbability(1, 1) * getResProbability(1, obsSequence[i + 1]) * probability[1];
    }

    if (state)
        return probabilityPrev[1];
    else
        return probabilityPrev[0];
}

// Вероятность наблюдаемой последовательности (через форвард процедуру)
double DishonestCasino::probabilityOfObsSeqF(std::vector<int> const &obsSequence) {
    return forwardProcedure(obsSequence, 0) + forwardProcedure(obsSequence, 1);
}

// Вероятность наблюдаемой последовательности (через бэкворд процедуру)
double DishonestCasino::probabilityOfObsSeqB(std::vector<int> const &obsSequence) {
    return getInDistribution(0) * getResProbability(0, obsSequence[0]) * backwardProcedure(obsSequence, 0) +
            getInDistribution(1) * getResProbability(1, obsSequence[0]) * backwardProcedure(obsSequence, 1);
}

// Алгоритм Витерби
std::vector<int> DishonestCasino::viterbiAlgorithm(std::vector<int> const &obsSequence){
    double probability[2] = {-log(getInDistribution(0)) - log(getResProbability(0, obsSequence[0])),
                             -log(getInDistribution(1)) - log(getResProbability(1, obsSequence[0]))},
           probabilityNext[2]{};
    std::vector<std::vector<int>> stSequence(2);

    unsigned long lengthOfObsSequence = obsSequence.size();
    double a, b;
    for(unsigned long i = 1; i < lengthOfObsSequence; ++i){
        a = probability[0] - log(getTrProbability(0, 0));
        b = probability[1] - log(getTrProbability(1, 0));

        if (a <= b){
            probabilityNext[0] = -log(getResProbability(0, obsSequence[i])) + a;
            stSequence[0].push_back(0);
        } else{
            probabilityNext[0] = -log(getResProbability(0, obsSequence[i])) + b;
            stSequence[0].push_back(1);
        }

        a = probability[0] - log(getTrProbability(0, 1));
        b = probability[1] - log(getTrProbability(1, 1));

        if (a <= b){
            probabilityNext[1] = -log(getResProbability(1, obsSequence[i])) + a;
            stSequence[1].push_back(0);
        } else{
            probabilityNext[1] = -log(getResProbability(1, obsSequence[i])) + b;
            stSequence[1].push_back(1);
        }

        probability[0] = probabilityNext[0];
        probability[1] = probabilityNext[1];
    }

    if (probability[0] <= probability[1]){
        stSequence[0].push_back(0);
        return stSequence[0];
    }else{
        stSequence[1].push_back(1);
        return stSequence[1];
    }
}

// Угадывание параметров HMM с помощью алгоритм Витерби
DishonestCasino* DishonestCasino::viterbiGuessing(std::vector<int> &obsSequence) {
    std::vector<int> vStSequence = viterbiAlgorithm(obsSequence);
    double tmpInDistribution[2]{}, tmpTrProbability[2][2]{}, tmpResProbability[2][6]{};
    unsigned long countTransitions[2][2]{}, countResults[2][6]{}, lengthOfObsSequence = obsSequence.size();;

    for(unsigned long i = 0; i < lengthOfObsSequence - 1; ++i) {
        ++countTransitions[vStSequence[i]][vStSequence[i + 1]];
        ++countResults[vStSequence[i]][obsSequence[i]];
    }
    ++countResults[vStSequence[lengthOfObsSequence - 1]][obsSequence[lengthOfObsSequence - 1]];

    tmpInDistribution[vStSequence[0]] = 1;
    for(int i = 0; i < 2; ++i)
        for(int j = 0; j < 2; ++j)
            if(countTransitions[i][j] != 0)
                tmpTrProbability[i][j] = static_cast<double>(countTransitions[i][j]) /
                        static_cast<double>(countTransitions[i][0] + countTransitions[i][1]);
            else
                tmpTrProbability[i][j] = 0;
    for(int i = 0; i < 2; ++i)
        for(int j = 0; j < 6; ++j)
            if(countResults[i][j] != 0)
                tmpResProbability[i][j] = static_cast<double>(countResults[i][j]) /
                        static_cast<double>(countResults[i][0] + countResults[i][1] + countResults[i][2] +
                        countResults[i][3] + countResults[i][4] + countResults[i][5]);
            else
                tmpResProbability[i][j] = 0;

    auto *newCasino = new DishonestCasino(tmpInDistribution, tmpTrProbability, tmpResProbability);
    return newCasino;
}

// Алгоритм Баума-Велча
DishonestCasino* DishonestCasino::baumWelchAlgorithm(PRNG& generator, std::vector<int> &obsSequence){
    // Задание точности для остановки алгоритма
    double prEps = 1.e-225, trEps = 1.e-225, resEps = 1.e-225;
    // Задание случайных начальных параметров
    unsigned long lengthOfObsSequence = obsSequence.size(), count{};
    for (int i = 0; i < lengthOfObsSequence; ++i)
        if(obsSequence[i] == 5) {++count;}
    double estProb = (double) count / lengthOfObsSequence * 2;
    if(estProb > 0.9) {estProb = 0.9;}
    if(estProb < 0.01) {estProb = 0.01;}

    double rInDistribution[2] = {0.5, 0.5}, rTrProbability[2][2] = {{0.5, 0.5}, {0.5, 0.5} },
        rResProbability[2][6] = {{1./6, 1./6, 1./6, 1./6, 1./6, 1./6, },
                                 {(1 - estProb) / 5, (1 - estProb) / 5, (1 - estProb) / 5, (1 - estProb) / 5, (1 - estProb) / 5, estProb}};

    auto *newCasino = new DishonestCasino(rInDistribution, rTrProbability, rResProbability),
            *casino = new DishonestCasino();

    // Поиск оптимальной HMM
    double prDif, trDif[2][2], resDif[2][6];
    do {
        casino = newCasino;
        double ysum[2]{}, esum[2][2]{}, yksum[2][6]{};
        double probability = casino->probabilityOfObsSeqF(obsSequence);

        // Высчитываем новые параметры
        for(unsigned long t = 0; t < lengthOfObsSequence; ++t) {
            // Вычисление ysum(i) и ysumk(i)
            for(int i = 0; i < 2; ++i) {
                double tmp = (casino->forwardProcedure(subVector(obsSequence, 0, t), i) *
                              casino->backwardProcedure(subVector(obsSequence, t + 1, lengthOfObsSequence - 1), i)) /
                              probability;
                ysum[i] += tmp;
                yksum[i][obsSequence[t]] += tmp;
            }
            // Вычисление esum[i][j]
            if (t != lengthOfObsSequence - 1)
                for(int i = 0; i < 2; ++i)
                    for(int j = 0; j < 2; ++j){
                        esum[i][j] += (casino->forwardProcedure(subVector(obsSequence, 0, t), i) *
                                       casino->getTrProbability(i, j) * casino->getResProbability(j, obsSequence[t + 1]) *
                                       casino->backwardProcedure(subVector(obsSequence, t + 2, lengthOfObsSequence - 1), j)) /
                                       probability;
                    }
            // Вычисление pi(i)
            if(t == 0) {
                for(int i = 0; i < 2; ++i)
                    rInDistribution[i] = ysum[i];
            }
            // Вычисление a(i,j)
            if(t == lengthOfObsSequence - 2) {
                for(int i = 0; i < 2; ++i)
                    for(int j = 0; j < 2; ++j)
                        rTrProbability[i][j] = esum[i][j] / ysum[i];
            }
            // Вычисление b(j,k)
            if(t == lengthOfObsSequence - 1) {
                for(int i = 0; i < 2; ++i)
                    for(int k = 0; k < 6; ++k)
                        rResProbability[i][k] = yksum[i][k] / ysum[i];
            }
        }
        newCasino = new DishonestCasino(rInDistribution, rTrProbability, rResProbability);
        // Высчитываем новые параметры

        // Считаем разницу по параметрам и вероятности
        prDif = fabs(probability - newCasino->probabilityOfObsSeqF(obsSequence));
        for(int i = 0; i < 2; ++i)
            for(int j = 0; j < 2; ++j)
                trDif[i][j] = fabs(newCasino->getTrProbability(i, j) - casino->getTrProbability(i, j));
        for(int i = 0; i < 2; ++i)
            for(int k = 0; k < 6; ++k)
                resDif[i][k] = fabs(newCasino->getResProbability(i, k) - casino->getResProbability(i, k));
        delete casino;
    } while(prDif > prEps && trDif[0][0] > trEps && trDif[0][1] > trEps && trDif[1][0] > trEps && trDif[1][1] > trEps &&
            resDif[0][0] > resEps && resDif[0][1] > resEps && resDif[0][2] > resEps && resDif[0][3] > resEps && resDif[0][4] > resEps && resDif[0][5] > resEps &&
            resDif[1][0] > resEps && resDif[1][1] > resEps && resDif[1][2] > resEps && resDif[1][3] > resEps && resDif[1][4] > resEps && resDif[1][5] > resEps);
    // Поиск оптимальной HMM
    return newCasino;
}

/* -------------------------------------------Далее вспомогательные функции-------------------------------------------*/

// Генерация наблюдаемого символа в зависимости от распределения для кости
int DishonestCasino::generateObs(PRNG& generator, int state) {
    double random = getRandomDouble(generator, 0, 1), sum = 0;
    if (random <= (sum += getResProbability(state, 0))){return 0;}
    if (random <= (sum += getResProbability(state, 1))){return 1;}
    if (random <= (sum += getResProbability(state, 2))){return 2;}
    if (random <= (sum += getResProbability(state, 3))){return 3;}
    if (random <= (sum += getResProbability(state, 4))){return 4;}
    else {return 5;}
}

// Генерация последовательности состояний и последовательности наблюдаемых символов
std::vector<std::vector<int>> DishonestCasino::generateSequences(PRNG& generator, unsigned long lengthOfStateSequence){
    std::vector<int> stSequence, obsSequence;
    stSequence.push_back(getRandomDouble(generator, 0, 1) <= getInDistribution(0) ? 0 : 1);
    obsSequence.push_back(generateObs(generator, stSequence[0]));

    for (unsigned long i = 1; i < lengthOfStateSequence; ++i){
        if (stSequence[i - 1])
            stSequence.push_back(getRandomDouble(generator, 0, 1) <= getTrProbability(1, 1) ? 1 : 0);

        else
            stSequence.push_back(getRandomDouble(generator, 0, 1) <= getTrProbability(0, 0) ? 0 : 1);

        obsSequence.push_back(generateObs(generator, stSequence[i]));
    }
    std::vector<std::vector<int>> resultSq = {stSequence, obsSequence};
    return resultSq;
}

// Взятие подвектора
template<typename T>
std::vector<T> DishonestCasino::subVector(std::vector<T> const &v, unsigned long m, unsigned long n) {
    auto first = v.begin() + m;
    auto last = v.begin() + n + 1;
    std::vector<T> vector(first, last);
    return vector;
}


