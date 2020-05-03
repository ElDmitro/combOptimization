#include <algorithm>
#include <cstdlib>
#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <random>
#include <math.h>
#define KNIGHT_DIAGS_NUM 4

int MAX_ITER = 51000;
double ITER_RATIO = 10;
std::random_device rd;
std::mt19937 gen(rd());


class BoardManager {
public:
    BoardManager() {}

    BoardManager(long long boardSize) {
        bSize = boardSize;
        resetStats();
    }

    BoardManager(const BoardManager &other) {
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

    BoardManager& operator=(const BoardManager &other) {
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

    void readBoard() {
        long long qNum_tmp;
        std::cin >> bSize >> qNum_tmp;

        resetStats();

        long long x, y;
        for (long long i = 0; i < qNum_tmp; ++i) {
            ++freeIdx;
            std::cin >> x >> y;
            put(x - 1, y - 1);
        }
    }

    void readBoard(std::ifstream &input_stream) {
        long long qNum_tmp;
        input_stream >> bSize >> qNum_tmp;

        resetStats();

        long long x, y;
        for (long long i = 0; i < qNum_tmp; ++i) {
            ++freeIdx;
            input_stream >> x >> y;
            put(x - 1, y - 1);
        }
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
        if (qIdx < freeIdx) {
            throw "HEH exception";
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

    long long evaluateCollisions() const {
        long long ncollisions = 0;
        for (auto value : colsCumSum) {
            ncollisions += value * (value - 1) / 2;
        }
        for (auto value : rowsCumSum) {
            ncollisions += value * (value - 1) / 2;
        }
        for (auto value : majorDiagsCumSum) {
            ncollisions += value * (value - 1) / 2;
        }
        for (auto value : minorDiagsCumSum) {
            ncollisions += value * (value - 1) / 2;
        }
        for (auto cumSum : knightDiagsCumSum) {
            for (auto value : cumSum) {
                ncollisions += value * (value - 1) / 2;
            }
        }

        return ncollisions;
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

    void printDebugLog() {
        for (long long i = 0; i < colsCumSum.size(); ++i) {
            std::cout << "Column " << i << " |" << colsCumSum[i] << std::endl;
        }
        std::cout << "#########" << std::endl;
        for (long long i = 0; i < rowsCumSum.size(); ++i) {
            std::cout << "Row " << i << " |" << rowsCumSum[i] << std::endl;
        }
        std::cout << "#########" << std::endl;

        for (long long i = 0; i < majorDiagsCumSum.size(); ++i) {
            std::cout << "Major Diag " << i - (bSize - 1) << " |" << majorDiagsCumSum[i] << std::endl;
        }
        std::cout << "#########" << std::endl;
        for (long long i = 0; i < minorDiagsCumSum.size(); ++i) {
            std::cout << "Minor Diag " << i - (bSize - 1) << " |" << minorDiagsCumSum[i] << std::endl;
        }
        std::cout << "#########" << std::endl;

        for (long long i = 0; i < 4; ++i) {
            for (long long j = 0; j < knightDiagsCumSum[i].size(); ++j) {
                std::cout << j / bSize << " " << j % bSize << std::endl;
                std::cout << "Knight diag type " << i << " |" << knightDiagsCumSum[i][j] << std::endl;
            }
            std::cout << "##########" << std::endl;
        }
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


long long naiveQueenPlacement(BoardManager &solution) {
    long long bSize = solution.getBoardSize();

    long long nCollisions = solution.getNCollisions();
    if (nCollisions > 0) {
        return 0;
    }
    const std::vector<long long> &colsCumSum = solution.getColsCumSum();
    const std::vector<long long> &rowsCumSum = solution.getRowsCumSum();

    long long j = 0;
    for (long long i = 0; i < bSize; ++i) {
        if (colsCumSum[i] > 0) {
            continue;
        }
        while ((j < bSize) && (rowsCumSum[j] > 0)) {
            ++j;
        }
        if (j >= bSize) {
            break;
        }

        solution.put(j, i);
    }

    if (solution.size() != bSize) {
        return 0;
    }

    return 1;
}


long long naiveQueenPlacement(const BoardManager &board, BoardManager &solution) {
    solution = board;
    long long bSize = solution.getBoardSize();

    long long nCollisions = solution.getNCollisions();
    if (nCollisions > 0) {
        return 0;
    }
    const std::vector<long long> &colsCumSum = solution.getColsCumSum();
    const std::vector<long long> &rowsCumSum = solution.getRowsCumSum();

    long long j = 0;
    for (long long i = 0; i < bSize; ++i) {
        if (colsCumSum[i] > 0) {
            continue;
        }
        while ((j < bSize) && (rowsCumSum[j] > 0)) {
            ++j;
        }
        if (j >= bSize) {
            break;
        }

        solution.put(j, i);
    }

    if (solution.size() != bSize) {
        return 0;
    }

    return 1;
}


void getOptimaRowMove(const BoardManager &board, long long qIdx,
                      long long &x_opt, long long &y_opt) {
    long long bSize = board.getBoardSize();
    long long initNCollisions = board.getNCollisions();

    long long row = board[qIdx].first;
    long long col = board[qIdx].second;

    long long optNCollisions = initNCollisions;
    x_opt = row;
    y_opt = col;

    long long posCollisions;
    long long currentCollisions = board.countPosCollisions(row, col, false);
    for (long long j = 0; j < bSize; ++j) {
        posCollisions = board.countPosCollisions(row, j, true) - 1;

        if ((initNCollisions + posCollisions - currentCollisions) < optNCollisions) {
            optNCollisions = initNCollisions + posCollisions - currentCollisions;
            y_opt = j;;
        }
    }
}


long long localSearch(BoardManager &solution, long long maxIter) {
    long long bsize = solution.getBoardSize();
    long long qNum = solution.size();

    long long qIdx;
    long long colIdx, rowIdx;
    for (size_t i = 0; i < maxIter; ++i) {
        qIdx = std::rand() % qNum;

        getOptimaRowMove(solution, qIdx, rowIdx, colIdx);
        solution.move(qIdx, rowIdx, colIdx);
    }

    return solution.getNCollisions() == 0;
}


long long colsPresearchCheck(const BoardManager& board) {
   return board.getNCollisions() == 0; 
}


long long getRandomRowFreeInit(const BoardManager &baseBoard, BoardManager &board) {
    board = baseBoard;
    long long bSize = board.getBoardSize();
    const std::vector<long long> rowsCumSum = board.getRowsCumSum();

    long long col;
    for (long long i = 0; i < bSize; ++i) {
        if (rowsCumSum[i] != 0) {
            continue;
        }
    
        col = std::rand() % bSize;
        board.put(i, col);
    }

    return board.size() == bSize;
}


long long getRandomRowFreeInit(BoardManager &board) {
    long long bSize = board.getBoardSize();
    const std::vector<long long> rowsCumSum = board.getRowsCumSum();

    long long col;
    for (long long i = 0; i < bSize; ++i) {
        if (rowsCumSum[i] != 0) {
            continue;
        }
    
        col = std::rand() % bSize;
        board.put(i, col);
    }

    return board.size() == bSize;
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


long long getRandomLineFreeInit(BoardManager &board) {
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

    std::random_shuffle(freeCols.begin(), freeCols.end());
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


typedef long long (*PrecheckFunc_t)(const BoardManager&);
typedef long long (*RandomInitFunc_t)(const BoardManager&, BoardManager&);
typedef long long (*OptimizerFunc_t)(BoardManager&, long long);
long long RandomMetaSearch(PrecheckFunc_t precheck, RandomInitFunc_t getRandomInit, OptimizerFunc_t optimizer,
                     const BoardManager& baseBoard, BoardManager &solution,
                     long long maxIter, long long locMaxIter) {
    long long success = precheck(baseBoard);
    if (!success) {
        return 0;
    }
    BoardManager initSolution;
    for (long long i = 0; i < maxIter; ++i) {
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


typedef long long (*ConfigOptimizerFunc_t)(PrecheckFunc_t, BoardManager&, long long&, double, double, double);
long long RandomMetaSearchAnn(PrecheckFunc_t precheck, RandomInitFunc_t getRandomInit, ConfigOptimizerFunc_t optimizer,
                     const BoardManager& baseBoard, BoardManager &solution,
                     long long maxIter, long long locMaxIter, double temp, double lower_temp, double alpha) {
    long long success;
    if (precheck) {
        success = precheck(baseBoard);
        if (!success) {
            return 0;
        }
    }

    BoardManager initSolution;
    for (long long i = 0; maxIter > 0; ++i) {
        success = getRandomInit(baseBoard, initSolution);
        if (!success) {
            continue;
        }
        long long locIters = std::min(locMaxIter, maxIter);
        maxIter -= locIters;
        success = optimizer(precheck, initSolution, locIters, temp, lower_temp, alpha);
        if (success) {
            solution = initSolution;
            return 1;
        }

        maxIter += locIters;
    }

    return 0;
}


/*
long long RandomMetaSearchAnnRest(PrecheckFunc_t precheck, RandomInitFunc_t getRandomInit, ConfigOptimizerFunc_t optimizer,
                     const BoardManager& baseBoard, BoardManager &solution,
                     long long maxIter, long long locMaxIter, double temp, double lower_temp, double alpha) {
    long long success;
    if (precheck) {
        success = precheck(baseBoard);
        if (!success) {
            return 0;
        }
    }

    BoardManager initSolution;
    getRandomInit(baseBoard, initSolution);
    for (long long i = 0; i < maxIter; ++i) {
        success = optimizer(precheck, initSolution, locMaxIter, temp, lower_temp, alpha);
        if (success) {
            solution = initSolution;
            return 1;
        }
    }

    return 0;
}
*/


bool getBerRandom(double prob) {
    std::bernoulli_distribution d(prob);
    return d(gen);
}


void metropolisStep(BoardManager &solution, double temp) {
    long long bSize = solution.getBoardSize();
    long long qNum = solution.size();
    long long nCollisions = solution.getNCollisions();
    long long freeIdx = solution.getFreeIdx();

    long long x, y;
    long long idx = std::rand() % (qNum - freeIdx) + freeIdx;
    x = solution[idx].first;
    y = solution[idx].second;

    long long yNew = std::rand() % bSize;
    if (yNew == y) {
        return;
    }
    long long currPosColls = solution.countPosCollisions(x, y, false);
    long long newPosColls = solution.countPosCollisions(x, yNew, true) - 1;

    long long nCollisionsNew = nCollisions - currPosColls + newPosColls;

    if (nCollisionsNew < nCollisions) {
        solution.move(idx, x, yNew);
    } else {
        if (getBerRandom(exp(-(nCollisionsNew - nCollisions) * 1. / temp))) {
            solution.move(idx, x, yNew);
        }
    }
}


void metropolisStepSwap(BoardManager &solution, double temp) {
    long long bSize = solution.getBoardSize();
    long long qNum = solution.size();
    long long nCollisions = solution.getNCollisions();
    long long freeIdx = solution.getFreeIdx();

    if (qNum - freeIdx == 1) {
        return;
    }

    long long x, y;
    long long a, b;
    // long long idx1 = std::rand() % (qNum - freeIdx) + freeIdx;
    // long long idx2 = std::rand() % (qNum - freeIdx) + freeIdx;
    std::uniform_int_distribution<int> distribution(freeIdx, qNum - 1);

    long long idx1 = distribution(gen);

    long long idx2 = distribution(gen);
    while (idx2 == idx1) {
        // idx2 = std::rand() % (qNum - freeIdx) + freeIdx;
        idx2 = distribution(gen);
    }
    
    x = solution[idx1].first;
    y = solution[idx1].second;

    a = solution[idx2].first;
    b = solution[idx2].second;

    solution.move(idx1, x, b);
    solution.move(idx2, a, y);

    long long nCollisionsNew = solution.getNCollisions();
    if (nCollisionsNew < nCollisions) {
        return;
    } else {
        if (getBerRandom(exp(-(nCollisionsNew - nCollisions) * 1. / temp))) {
            return;
        }

        solution.move(idx1, x, y);
        solution.move(idx2, a, b);
    }
}

/*
long long simulatedAnnealing(PrecheckFunc_t precheck, BoardManager &solution, long long maxIter,
                             double temp, double alpha=0.9) {
    long long success;
    if (precheck) {
        success = precheck(solution);
        if (!success) {
            return 0;
        }
    }

    for (long long i = 0; i < maxIter; ++i) {
        metropolisStepSwap(solution, temp);
        temp = temp * alpha; 
    }

    return solution.getNCollisions() == 0;
}
*/


long long simulatedAnnealing(PrecheckFunc_t precheck, BoardManager &solution, long long &maxIter,
                             double temp, double lower_temp, double alpha=0.9) {
    long long success;
    if (precheck) {
        success = precheck(solution);
        if (!success) {
            return 0;
        }
    }

    // long long waitIters = 50;
    long long waitIters = 10. * (solution.getBoardSize() - solution.getFreeIdx()) / 999 + 49850. / 999;
    // long long waitIters = 10 * (solution.getBoardSize() - solution.getFreeIdx());

    long long prev = solution.getNCollisions();
    long long plateuIters = 0;
    for (long long i = 0; maxIter > 0; ++i) {
        --maxIter;
        if (solution.getNCollisions() == 0) {
            return 1;
        }
        if (plateuIters > waitIters) {
            return 0;
        }
        if (temp < lower_temp) {
            return solution.getNCollisions() == 0;
        }
        metropolisStepSwap(solution, temp);
        temp = temp * alpha;

        if (prev <= solution.getNCollisions()) {
            ++plateuIters;
        } else {
            plateuIters = 0;
            prev = solution.getNCollisions();
        }

     }

    return solution.getNCollisions() == 0;
}


int main(int argc, char *argv[]) {
    std::srand(time(NULL));

    BoardManager boardManager;
    if (argc > 1) {
        std::ifstream input_stream;
        input_stream.open(argv[1]);

        boardManager.readBoard(input_stream);
    } else {
        boardManager.readBoard();
    }

    if (argc > 2) {
        std::cout << boardManager.getNCollisions() << std::endl;
    }

    if (argc > 2) {
        boardManager.showBoard();
        std::cout << std::endl;
    }
    int bSize = boardManager.getBoardSize();

    BoardManager solution;
    long long success = RandomMetaSearchAnn(NULL, getRandomLineFreeInit, simulatedAnnealing,
                                   boardManager, solution, 8000000ll, 5000ll, 2.7 * (bSize / 2), 0.1, 0.99);
    if (argc > 2) {
        solution.showBoard();
    }

    if (success and (solution.getNCollisions() == 0)) {
        std::cout << "YES" << std::endl;
        std::cout << solution;
    } else {
        std::cout << "NO" << std::endl;
    }
}
