#include <algorithm>
#include <cstdlib>
#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <random>
#include <math.h>
#include <limits.h>
#define KNIGHT_DIAGS_NUM 4

// TODO(eldmitro): check if it norm
std::random_device rd;
std::mt19937 gen(rd());


class BoardManager {
public:
    BoardManager() {}

    BoardManager(long long boardSize) {
        bSize = boardSize;
        resetStats();
    }

    BoardManager(long long boardSize, long long firstFreeIdx) {
        bSize = boardSize;
        resetStats();

        freeIdx = firstFreeIdx;
    }

    void copy(const BoardManager &other) {
        bSize = other.bSize;
        qNum = other.qNum;
        freeIdx = other.freeIdx;

        qCoords = other.qCoords;
        nCollisions = other.nCollisions;

        colsCumSum = other.colsCumSum;
        rowsCumSum = other.rowsCumSum;

        majorDiagsCumSum = other.majorDiagsCumSum;
        minorDiagsCumSum = other.minorDiagsCumSum;

        knightDiagsCumSum = other.knightDiagsCumSum;
        knightMoves = other.knightMoves;

    }

    BoardManager(const BoardManager &other) {
        copy(other);
    }

    BoardManager& operator=(const BoardManager &other) {
        copy(other);
    }

    long long getDiagIdx(long long x, long long y) const {
        return y - x + bSize - 1;
    }

    long long getMinorDiagIdx(long long x , long long y) const {
        return getDiagIdx(x, bSize - y - 1);
    }

    long long getKnightDiagIdx(long long x, long long y,
                         const std::pair<long long, long long> move) const {
        long long m_x = move.first;
        long long m_y = move.second;
        if (m_y < 0) {
            m_y = -m_y;
            y = bSize - y - 1;
        }

        long long x_steps = x / m_x;
        long long y_steps = y / m_y;
        long long steps = std::max(std::min(x_steps, y_steps), 0ll);

        x = x - m_x * steps;
        y = y - m_y * steps;

        if (m_x < m_y) {
            std::swap(x, y);
        }

        if (x < 2) {
            return x * bSize + y;
        }

        return 2 * bSize + x - 2;
    }

    void put(long long x, long long y) {
        if ((x < 0) || (x >= bSize)) {
            throw "PUT: out of range";
        }
        if ((y < 0) || (y >= bSize)) {
            throw "PUT: out of range";
        }

        ++qNum;
        qCoords.push_back(std::make_pair(x, y));
        updateStatsPut(x, y);
    }

    void move(long long qIdx, long long x, long long y) {
        std::pair<long long, long long> coord = qCoords[qIdx];
        long long x_old = coord.first;
        long long y_old = coord.second;

        if ((x == x_old) && (y == y_old)) {
            return;
        }

        qCoords[qIdx] = std::make_pair(x, y);
        updateStatsTake(x_old, y_old);
        updateStatsPut(x, y);
    }

    long long getBoardSize() const {
        return bSize;
    }

    long long size() const {
        return qNum;
    }

    std::pair<long long, long long> operator[](size_t i) {
        return qCoords[i];
    }

    std::pair<long long, long long> operator[](size_t i) const {
        return qCoords[i];
    }

    long long getNCollisions() const {
        return nCollisions;
    }

    long long countPosCollisions(long long x, long long y, bool pseudo_put) const {
        long long diff = !pseudo_put;

        long long ncollisions = 0;
        ncollisions += colsCumSum[y] - diff;
        ncollisions += rowsCumSum[x] - diff;

        long long majorDIdx = getDiagIdx(x, y);
        long long minorDIdx = getMinorDiagIdx(x, y);

        ncollisions += majorDiagsCumSum[majorDIdx] - diff;
        ncollisions += minorDiagsCumSum[minorDIdx] - diff;

        long long idx;
        for (long long i = 0; i < KNIGHT_DIAGS_NUM; ++i) {
            idx = getKnightDiagIdx(x, y, knightMoves[i]);
            ncollisions += knightDiagsCumSum[i][idx] - diff;
        }

        return ncollisions;
    }

    void showBoard() {
        std::vector<std::vector<bool>> board(bSize, std::vector<bool>(bSize, false));

        for (auto coord : qCoords) {
            board[coord.first][coord.second] = true;
        }

        for (auto row : board) {
            for (auto el : row) {
                if (el) {
                    std::cout << '*';
                    continue;
                }

                std::cout << '#';
            }

            std::cout << std::endl;
        }
    }

    long long getFreeIdx() const {
        return freeIdx;
    }

    const std::vector<long long>& getColsCumSum() const {
        return colsCumSum;
    }

    const std::vector<long long>& getRowsCumSum() const {
        return rowsCumSum;
    }

    const std::vector<long long>& getMajorDCumSum() const {
        return majorDiagsCumSum;
    }

    const std::vector<long long>& getMinorDCumSum() const {
        return minorDiagsCumSum;
    }

    const std::vector<std::vector<long long>>& getKnightDCumSum() const {
        return knightDiagsCumSum;
    }

private:
    long long bSize;
    long long qNum;
    long long nCollisions;
    long long freeIdx;
    std::vector<std::pair<long long, long long>> qCoords;

    std::vector<long long> colsCumSum;
    std::vector<long long> rowsCumSum;

    std::vector<long long> majorDiagsCumSum;
    std::vector<long long> minorDiagsCumSum;

    std::vector<std::vector<long long>> knightDiagsCumSum;
    std::vector<std::pair<long long, long long>> knightMoves;


    void resetStats() {
        qNum = 0;
        nCollisions = 0;
        freeIdx = 0;
        qCoords.clear();
        colsCumSum.assign(bSize, 0);
        rowsCumSum.assign(bSize, 0);

        majorDiagsCumSum.assign(2 * bSize - 1, 0);
        minorDiagsCumSum.assign(2 * bSize - 1, 0);

        knightDiagsCumSum.assign(KNIGHT_DIAGS_NUM, std::vector<long long>(0));
        for (long long i = 0; i < KNIGHT_DIAGS_NUM; ++i) {
            for (long long j = 0; j < 3 * bSize - 2; ++j) {
                knightDiagsCumSum[i].push_back(0);
            }
        }

        knightMoves.push_back(std::make_pair(2, 1));
        knightMoves.push_back(std::make_pair(1, 2));
        knightMoves.push_back(std::make_pair(2, -1));
        knightMoves.push_back(std::make_pair(1, -2));
    }

    void updateStatsPut(long long x, long long y) {
        nCollisions += colsCumSum[y];
        ++colsCumSum[y];
        nCollisions += rowsCumSum[x];
        ++rowsCumSum[x];

        long long majorDIdx = getDiagIdx(x, y);
        long long minorDIdx = getMinorDiagIdx(x, y);

        nCollisions += majorDiagsCumSum[majorDIdx];
        ++majorDiagsCumSum[majorDIdx];
        nCollisions += minorDiagsCumSum[minorDIdx];
        ++minorDiagsCumSum[minorDIdx];

        long long idx;
        for (long long i = 0; i < KNIGHT_DIAGS_NUM; ++i) {
            idx = getKnightDiagIdx(x, y, knightMoves[i]);
            nCollisions += knightDiagsCumSum[i][idx];
            ++knightDiagsCumSum[i][idx];
        }
    }

    void updateStatsTake(long long x, long long y) {
        nCollisions -= colsCumSum[y] - 1;
        --colsCumSum[y];
        nCollisions -= rowsCumSum[x] - 1;
        --rowsCumSum[x];

        long long majorDIdx = getDiagIdx(x, y);
        long long minorDIdx = getMinorDiagIdx(x, y);

        nCollisions -= majorDiagsCumSum[majorDIdx] - 1;
        --majorDiagsCumSum[majorDIdx];
        nCollisions -= minorDiagsCumSum[minorDIdx] - 1;
        --minorDiagsCumSum[minorDIdx];

        long long idx;
        for (long long i = 0; i < KNIGHT_DIAGS_NUM; ++i) {
            idx = getKnightDiagIdx(x, y, knightMoves[i]);
            nCollisions -= knightDiagsCumSum[i][idx] - 1;
            --knightDiagsCumSum[i][idx];
        }
    }
};


std::istream& operator>>(std::istream &input_stream, BoardManager &obj) {
    long long bSize;
    long long qNum_input;
    input_stream >> bSize >> qNum_input;

    obj = BoardManager(bSize, qNum_input);

    long long x, y;
    for (long long i = 0; i < qNum_input; ++i) {
        input_stream >> x >> y;
        obj.put(x - 1, y - 1);
    }
}


std::ostream& operator<<(std::ostream &output_stream, const BoardManager &obj) {
    std::vector<long long> rows(obj.getBoardSize());

    std::pair<long long, long long> pair;
    for (long long i = 0; i < obj.size(); ++i) {
        pair = obj[i];

        rows[pair.first] = pair.second;
    }

    for (auto value : rows) {
        output_stream << value + 1 << std::endl;
    }

    return output_stream;
}


long long getRandomLineFreeInit(const BoardManager &baseBoard, BoardManager &board) {
    board = baseBoard;
    long long bSize = board.getBoardSize();
    const std::vector<long long> rowsCumSum = board.getRowsCumSum();
    const std::vector<long long> colsCumSum = board.getColsCumSum();

    long long col;
    std::vector<int> freeCols;
    for (long long i = 0; i < bSize; ++i) {
        if (colsCumSum[i] > 0) {
            continue;
        }

        freeCols.push_back(i);
    }

    std::shuffle(freeCols.begin(), freeCols.end(), gen);
    long long k = 0;
    for (long long i = 0; i < bSize; ++i) {
        if (rowsCumSum[i] > 0) {
            continue;
        }

        board.put(i, freeCols[k]);
        ++k;
    }

    return board.size() == bSize;
}


bool isTerminatedState(const BoardManager &state) {
    return state.getNCollisions() == 0;
}


bool isTerminatedState(const BoardManager &state, long long &prevFunc, long long &numberPlateu, long long numberPlateuLimit=100) {
    if (state.getNCollisions() == prevFunc) {
        ++numberPlateu;
    } else {
        numberPlateu = 0;
        prevFunc = state.getNCollisions();
    }

    if (numberPlateu > numberPlateuLimit) {
        return 1;
    }
    return state.getNCollisions() == 0;
}


void getRandomPair(BoardManager &solution, long long &a, long long &b) {
    long long freeIdx = solution.getFreeIdx();
    long long qNum = solution.size();

    std::uniform_int_distribution<int> distribution(0, qNum - 1);

    a = distribution(gen);
    b = distribution(gen);
    while (b == a) {
        b = distribution(gen);
    }
}


long long swapPair(BoardManager &solution, long long i, long long j) {
    long long x, y;
    long long a, b;

    x = solution[i].first;
    y = solution[i].second;

    a = solution[j].first;
    b = solution[j].second;

    solution.move(i, x, b);
    solution.move(j, a, y);
}


long long SwapnScore(BoardManager &solution, long long i, long long j) {
    long long nCollisions = solution.getNCollisions();

    swapPair(solution, i, j);
    return solution.getNCollisions() - nCollisions;
}


long long swapPair(BoardManager &solution, long long i, long long j,
                   std::vector<std::vector<long long>> &tabu, long long iter, long long L=10) {
    long long x, y;
    long long a, b;

    x = solution[i].first;
    y = solution[i].second;

    a = solution[j].first;
    b = solution[j].second;

    if ((tabu[i][b] > iter) || (tabu[j][y] > iter)) {
        return 1;
    }

    solution.move(i, x, b);
    solution.move(j, a, y);

    tabu[i][b] = iter + L;
    tabu[j][y] = iter + L;
    return 0;
}


long long SwapnScore(BoardManager &solution, long long i, long long j,
                     std::vector<std::vector<long long>> &tabu,
                     long long iter, long long L=10) {
    long long nCollisions = solution.getNCollisions();

    long long tabuFlag = swapPair(solution, i, j, tabu, iter);
    return solution.getNCollisions() - nCollisions;
}


long long scoreSwap(BoardManager &solution, long long i, long long j) {
    long long score = SwapnScore(solution, i, j);

    swapPair(solution, i, j);

    return score;
}


long long localSearch(BoardManager &initSolution, long long steps) {
    long long score;
    long long bSize = initSolution.getBoardSize();
    long long freeIdx = initSolution.getFreeIdx();

    long long a, b;
    long long prevFunc;
    long long plateuNumber = 0;
    for (int i = 0; i < steps; ++i) {
        if (isTerminatedState(initSolution, prevFunc, plateuNumber, 0.5 * (bSize - freeIdx))) {
            return initSolution.getNCollisions() == 0;
        }

        getRandomPair(initSolution, a, b); 
        score = SwapnScore(initSolution, a, b);

        if (score >= 0) {
            swapPair(initSolution, a, b);
        }
    }
}


long long tabuLocalSearch(BoardManager &initSolution, long long steps) {
    long long score;
    long long qNum = initSolution.size();
    long long bSize = initSolution.getBoardSize();
    long long freeIdx = initSolution.getFreeIdx();

    long long a, b;
    long long prevFunc;
    long long plateuNumber = 0;
    std::vector<std::vector<long long>> tabu(qNum, std::vector<long long>(bSize, 0));
    for (int i = 0; i < steps; ++i) {
        if (isTerminatedState(initSolution, prevFunc, plateuNumber)) {
            return initSolution.getNCollisions() == 0;
        }

        getRandomPair(initSolution, a, b); 
        score = SwapnScore(initSolution, a, b, tabu, i, 10);

        if (score >= 0) {
            swapPair(initSolution, a, b);
        }
    }
}


void getHeavyPair(BoardManager &solution, long long &a, long long &b) {
    long long maxPosCollision = -1;
    long long freeIdx = solution.getFreeIdx();
    long long qNum = solution.size();

    long long x, y;
    long long currPosCollisions;
    for (int i = freeIdx; i < qNum; ++i) {
        currPosCollisions = solution.countPosCollisions(solution[i].first,
                                                        solution[i].second,
                                                        false);

        if (currPosCollisions > maxPosCollision) {
            maxPosCollision = currPosCollisions;
            a = i;
        }
    }

    long long score;
    long long nCollisions = solution.getNCollisions();
    long long optNCollisions = nCollisions;
    b = a;
    for (int i = freeIdx; i < qNum; ++i) {
        if (i == a) {
            continue;
        }

        score = scoreSwap(solution, a, i);

        if (nCollisions + score < optNCollisions) {
            b = i;
            optNCollisions = nCollisions + score;
        }
    }
}



long long localSearchHeavy(BoardManager &initSolution, long long steps) {
    long long score;
    long long bSize = initSolution.getBoardSize();
    long long freeIdx = initSolution.getFreeIdx();

    long long a, b;
    long long prevFunc;
    long long plateuNumber = 0;
    for (int i = 0; i < steps; ++i) {
        if (isTerminatedState(initSolution, prevFunc, plateuNumber, 0.5 * (bSize - freeIdx))) {
            return initSolution.getNCollisions() == 0;
        }

        getHeavyPair(initSolution, a, b); 
        score = SwapnScore(initSolution, a, b);

        if (score >= 0) {
            swapPair(initSolution, a, b);
        }
    }
}


bool getBerRandom(double prob) {
    std::bernoulli_distribution d(prob);
    return d(gen);
}

long long simulatedAnnealing(BoardManager &initSolution, long long steps, double temp, double alpha=0.9) {
    long long score;
    long long bSize = initSolution.getBoardSize();
    long long freeIdx = initSolution.getFreeIdx();

    long long a, b;
    long long prevFunc;
    long long plateuNumber = 0;
    for (int i = 0; i < steps; ++i) {
        if (isTerminatedState(initSolution, prevFunc, plateuNumber, (bSize - freeIdx))) {
            return initSolution.getNCollisions() == 0;
        }

        getRandomPair(initSolution, a, b); 
        score = SwapnScore(initSolution, a, b);

        if ((score >= 0) && (!getBerRandom(exp(-score * 1. / temp)))) {
            swapPair(initSolution, a, b);
        }
        temp = alpha * temp;
    }
}


int MAX_TIME = 9000;
typedef long long (*RandomInitFunc_t)(const BoardManager&, BoardManager&);
typedef long long (*OptimizerFunc_t)(BoardManager&, long long);
long long RandomMetaSearch(RandomInitFunc_t getRandomInit, OptimizerFunc_t optimizer,
                     const BoardManager& baseBoard, BoardManager &solution,
                     long long maxIter, long long locMaxIter) {
    auto start = std::chrono::high_resolution_clock::now();
    long long success;
    BoardManager initSolution;
    for (long long i = 0; i < maxIter; ++i) {
        auto point = std::chrono::high_resolution_clock::now();
        if (std::chrono::duration_cast<std::chrono::milliseconds>(point - start).count() > MAX_TIME) {
            return 0;
        }
        success = getRandomInit(baseBoard, initSolution);
        if (!success) {
            continue;
        }
        success = optimizer(initSolution, locMaxIter);
        if (success) {
            solution = initSolution;
            return 1;
        }
    }

    return 0;
}

typedef long long (*AnnealOptimizerFunc_t)(BoardManager&, long long, double, double);
long long RandomMetaSearch(RandomInitFunc_t getRandomInit, AnnealOptimizerFunc_t optimizer,
                     const BoardManager& baseBoard, BoardManager &solution,
                     long long maxIter, long long locMaxIter, double temp, double alpha) {
    auto start = std::chrono::high_resolution_clock::now();
    long long success;
    BoardManager initSolution;
    for (long long i = 0; i < maxIter; ++i) {
        auto point = std::chrono::high_resolution_clock::now();
        if (std::chrono::duration_cast<std::chrono::milliseconds>(point - start).count() > MAX_TIME) {
            return 0;
        }
        success = getRandomInit(baseBoard, initSolution);
        if (!success) {
            continue;
        }
        success = optimizer(initSolution, locMaxIter, temp, alpha);
        if (success) {
            solution = initSolution;
            return 1;
        }
    }

    return 0;
}


int main(int argc, char *argv[]) {
    std::srand(time(NULL));

    BoardManager boardManager;
    BoardManager empty(3);
    getRandomLineFreeInit(empty, boardManager);

    int bSize = boardManager.getBoardSize();
    long long a, b;
    std::cout << boardManager.getNCollisions() << std::endl;
    boardManager.showBoard();
    std::cout << std::endl;
    getRandomPair(boardManager, a, b);
    SwapnScore(boardManager, a, b);
    std::cout << boardManager.getNCollisions() << std::endl;

    BoardManager test(bSize);
    for (int i = 0; i < boardManager.size(); ++i) {
        a = boardManager[i].first;
        b = boardManager[i].second;

        test.put(a, b);
    }

    std::cout << test.getNCollisions() << std::endl;

    if (argc > 2) {
        boardManager.showBoard();
    }
}
