#include "DBReader.h"
#include "IndexReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "structureto3diseqdist.h"
#include "StructureUtil.h"
#include "TMaligner.h"
#include "Coordinate16.h"
#include "MemoryMapped.h"
#include "algorithm"

#ifdef OPENMP
#include <omp.h>
#endif

struct Chain {
    Chain() {}
    Chain(unsigned int complexId, unsigned int chainKey) : complexId(complexId), chainKey(chainKey) {}
    unsigned int complexId;
    unsigned int chainKey;
    std::vector<float> caVecX;
    std::vector<float> caVecY;
    std::vector<float> caVecZ;
    int startPos;
    unsigned int length;

    void setLength(unsigned int len){
        length = len;
    }

    void setCaData(std::vector<float> caX, std::vector<float> caY, std::vector<float> caZ){
        caVecX = caX;
        caVecY = caY;
        caVecZ = caZ;
    }

    void setStartPos(int sPos){
        startPos = sPos;
    }

    bool operator<(const Chain &o) const {
        if (complexId < o.complexId ) {
            return true;
        }
        if (complexId > o.complexId) {
            return false;
        }
        if (chainKey < o.chainKey) {
            return true;
        }
        if (chainKey > o.chainKey) {
            return false;
        }
        return false;
    }

    bool operator>(const Chain &o) const {
        if (complexId > o.complexId ) {
            return true;
        }
        if (complexId < o.complexId) {
            return false;
        }
        if (chainKey > o.chainKey) {
            return true;
        }
        if (chainKey < o.chainKey) {
            return false;
        }
        return false;
    }
    bool operator==(const Chain &o) const {
        return complexId == o.complexId && chainKey == o.chainKey;
    }
};

struct FeatureVector{
    FeatureVector(){}
    FeatureVector(float u[3][3], float t[3]) {
        features[0]=0;
        float w2 = (u[0][0] + u[1][1] + u[2][2] + 1)/4;
        float x2 = -(u[1][1]+u[2][2])/2;
        float y2 = (1-u[2][2])/2;
        if (w2>0) {
            features[0] = sqrt(w2);
            features[1] = (u[1][2] - u[2][1]) / (4 * features[0]);
            features[2] = (u[2][0] - u[0][2]) / (4 * features[0]);
            features[3] = (u[0][1] - u[1][0]) / (4 * features[0]);
        } else if (x2>0) {
            features[0] = 0.0;
            features[1] = sqrt(x2);
            features[2] = u[0][1] / (2 * features[1]);
            features[3] = u[0][2] / (2 * features[1]);
        } else if (y2>0) {
            features[0] = 0.0;
            features[1] = 0.0;
            features[2] = sqrt(y2);
            features[3] = u[1][2] / (2 * features[2]);
        } else {
            features[0] = 0.0;
            features[1] = 0.0;
            features[2] = 0.0;
            features[3] = 1.0;
        }
        for (size_t i=0; i<3; i++){
            features[4 + i]=t[i];
        }
        label = 0;
    }
    float features[7];
    int label;

    void resetLabel(){
        label=0;
    }

    double getDistance(const FeatureVector &o){
        double dist = 0;
        for (size_t i=0; i<7; i++){
            dist += std::pow(features[i] - o.features[i], 2);
        }
        dist = std::sqrt(dist);
        return dist;
    }
};

struct ChainToChainAln {
    ChainToChainAln() {}
    ChainToChainAln(Chain qChain, Chain dbChain, std::string backtrace): qChain(qChain), dbChain(dbChain), backtrace(backtrace) {}
    ChainToChainAln (Chain qInputChain, Chain dbInputChain, float * queryCaData, float * targetCaData, Matcher::result_t alnResult, TMaligner::TMscoreResult tmResult) {
        int qStartPos = alnResult.qStartPos;
        int dbStartPos = alnResult.dbStartPos;
        unsigned int qLength = alnResult.qLen;
        unsigned int dbLength = alnResult.dbLen;
        std::string  inputBacktrace = alnResult.backtrace;
        std::vector<float> qCaXVec;
        std::vector<float> qCaYVec;
        std::vector<float> qCaZVec;
        std::vector<float> dbCaXVec;
        std::vector<float> dbCaYVec;
        std::vector<float> dbCaZVec;
//        int numMatches = 0;
        std::string newBacktrace;
        int qPos = qStartPos;
        int dbPos = dbStartPos;
        int qXPos =0;
        int qYPos = qLength;
        int qZPos = qLength*2;
        int dbXPos = 0;
        int dbYPos = dbLength;
        int dbZPos = dbLength*2;
        for (size_t j=0; j<inputBacktrace.size(); j++){
            char cigar = inputBacktrace[j];
            switch (cigar) {
                case 'M':
                    newBacktrace += "M";
                    qCaXVec.emplace_back(queryCaData[qXPos + qPos]);
                    qCaYVec.emplace_back(queryCaData[qYPos + qPos]);
                    qCaZVec.emplace_back(queryCaData[qZPos + qPos]);
                    dbCaXVec.emplace_back(targetCaData[dbXPos + dbPos]);
                    dbCaYVec.emplace_back(targetCaData[dbYPos + dbPos]);
                    dbCaZVec.emplace_back(targetCaData[dbZPos + dbPos]);
                    qPos++;
                    dbPos++;
                    break;
                case 'I':
                    qPos++;
                    break;
                case 'D':
                    dbPos++;
                    break;
                default:
                    // TODO
                    // alert errors in backtrace
                    break;
            }
        }
        qInputChain.setStartPos(0);
        qInputChain.setCaData(qCaXVec, qCaYVec, qCaZVec);
        qInputChain.setLength(qLength);
        dbInputChain.setStartPos(0);
        dbInputChain.setCaData(dbCaXVec, dbCaYVec, dbCaZVec);
        dbInputChain.setLength(dbLength);
        qChain = qInputChain;
        dbChain = dbInputChain;
        backtrace = newBacktrace;
        alnLength = inputBacktrace.size();
        featureVector = FeatureVector(tmResult.u, tmResult.t);
    }
    Chain qChain;
    Chain dbChain;
    std::string backtrace;
    unsigned int alnLength;
    FeatureVector featureVector;
};

struct DistAndIndexPair{
    DistAndIndexPair(std::pair<unsigned int, unsigned int> indexPair, double distance) : indexPair(indexPair), distance(distance){}
    std::pair<unsigned int, unsigned int> indexPair;
    double distance;
};

struct Complex {
    Complex() {}
    Complex(unsigned int complexId, std::vector<unsigned int> chainKeys) : complexId(complexId), chainKeys(chainKeys) {
        alnVec = std::vector<ChainToChainAln>();
    }
    Complex(unsigned int complexId, std::vector<unsigned int> chainKeys, std::vector<ChainToChainAln> alnVec) : complexId(complexId), chainKeys(chainKeys), alnVec(alnVec) {}
    unsigned int complexId;
    std::vector<unsigned int> chainKeys;
    std::vector<ChainToChainAln> alnVec;

    void filterAlnVec(float compatibleCheckRatio){
        if (alnVec.empty()) {
            return;
        }
        std::vector<ChainToChainAln> newAlnVec;
        std::vector<ChainToChainAln> currDbComplexAlnVec;

        ChainToChainAln firstAln = alnVec[0];
        unsigned int dbPrevComplexId = firstAln.dbChain.complexId;
        unsigned int qChainCount = 1;
        unsigned int qPrevChainKey = firstAln.qChain.chainKey;

        for (size_t i=0; i < alnVec.size(); i++) {
            ChainToChainAln aln = alnVec[i];
            if (aln.dbChain.complexId != dbPrevComplexId) {
                if (qChainCount >= (float)(chainKeys.size() * compatibleCheckRatio)) {
                    newAlnVec.insert(newAlnVec.end(), currDbComplexAlnVec.begin(), currDbComplexAlnVec.end());
                }
                currDbComplexAlnVec.clear();
                dbPrevComplexId = aln.dbChain.complexId;
                qChainCount = 1;
                qPrevChainKey = aln.qChain.chainKey;
                currDbComplexAlnVec.emplace_back(aln);
            } else if (aln.qChain.chainKey != qPrevChainKey) {
                qChainCount ++;
                qPrevChainKey = aln.qChain.chainKey;
                currDbComplexAlnVec.emplace_back(aln);
            } else {
                currDbComplexAlnVec.emplace_back(aln);
            }
        }
        if (qChainCount >= (float)(chainKeys.size() * compatibleCheckRatio)) {
            newAlnVec.insert(newAlnVec.end(), currDbComplexAlnVec.begin(), currDbComplexAlnVec.end());
        }
        alnVec = newAlnVec;
    }

    void normalize(){
        double length = (double)(alnVec.size());
        double mean[7] = {0,0,0,0,0,0,0};
        double var[7] = {0,0,0,0,0,0,0};
        double sd[7] = {0,0,0,0,0,0,0};
        for (size_t i=0; i < length; i++){
            for (size_t j=0; j<7; j++){
                double val = (double)(alnVec[i].featureVector.features[j]);
                mean[j] += val / length;
            }
        }
        for (size_t i=0; i < length; i++){
            for (size_t j=0; j<7; j++){
                double val = (double)(alnVec[i].featureVector.features[j]);
                val -= mean[j];
                val *= val;
                var[j] += val / length;
            }
        }
        for (size_t j=0; j<7; j++){
            sd[j] = std::sqrt(var[j]);
        }
        for (size_t i=0; i < length; i++){
            for (size_t j=0; j<7; j++){
                alnVec[i].featureVector.features[j] = sd[j] == 0 ? 0 : (alnVec[i].featureVector.features[j] - mean[j]) / sd[j];
            }
        }
        return;
    }
};

struct ComplexToComplexAln{
    ComplexToComplexAln() {}
    ComplexToComplexAln(unsigned int qComplexId, unsigned int dbComplexId): qComplexId(qComplexId), dbComplexId(dbComplexId){}
    ComplexToComplexAln(ChainToChainAln aln) {
        qComplexId = aln.qChain.complexId;
        dbComplexId = aln.dbChain.complexId;
        qChainKeys = {aln.qChain.chainKey};
        dbChainKeys = {aln.dbChain.chainKey};
        backtrace = aln.backtrace;
        qLength = aln.qChain.length;
        dbLength = aln.dbChain.length;
        alnLength = aln.alnLength;
        qCaXVec = aln.qChain.caVecX;
        qCaYVec = aln.qChain.caVecY;
        qCaZVec = aln.qChain.caVecZ;
        dbCaXVec = aln.dbChain.caVecX;
        dbCaYVec = aln.dbChain.caVecY;
        dbCaZVec = aln.dbChain.caVecZ;
    }
    unsigned int qComplexId;
    unsigned int dbComplexId;
    std::vector<unsigned int> qChainKeys;
    std::vector<unsigned int> dbChainKeys;
    unsigned int qLength;
    unsigned int dbLength;
    unsigned int alnLength;
    std::string backtrace;
    std::vector<float> qCaXVec;
    std::vector<float> qCaYVec;
    std::vector<float> qCaZVec;
    std::vector<float> dbCaXVec;
    std::vector<float> dbCaYVec;
    std::vector<float> dbCaZVec;
    double tmScore;

    void setTmScore(double tmscore){
        tmScore = tmscore;
    }

    void appendChainToChainAln(ChainToChainAln aln){
        qChainKeys.emplace_back(aln.qChain.chainKey);
        dbChainKeys.emplace_back(aln.dbChain.chainKey);
        backtrace += aln.backtrace;
        qLength += aln.qChain.length;
        dbLength += aln.dbChain.length;
        alnLength += aln.alnLength;
        qCaXVec.insert(qCaXVec.end(), aln.qChain.caVecX.begin(), aln.qChain.caVecX.end());
        qCaYVec.insert(qCaYVec.end(), aln.qChain.caVecY.begin(), aln.qChain.caVecY.end());
        qCaZVec.insert(qCaZVec.end(), aln.qChain.caVecZ.begin(), aln.qChain.caVecZ.end());
        dbCaXVec.insert(dbCaXVec.end(), aln.dbChain.caVecX.begin(), aln.dbChain.caVecX.end());
        dbCaYVec.insert(dbCaYVec.end(), aln.dbChain.caVecY.begin(), aln.dbChain.caVecY.end());
        dbCaZVec.insert(dbCaZVec.end(), aln.dbChain.caVecZ.begin(), aln.dbChain.caVecZ.end());
    }
    void reset(){
        qChainKeys.clear();
        dbChainKeys.clear();
        qLength = 0;
        dbLength = 0;
        alnLength = 0;
        backtrace = "";
        qCaXVec.clear();
        qCaYVec.clear();
        qCaZVec.clear();
        dbCaXVec.clear();
        dbCaYVec.clear();
        dbCaZVec.clear();
    }

    void reset(unsigned int newDbId){
        dbComplexId = newDbId;
        reset();
    }
};

struct OutputLine{
    OutputLine(unsigned int qComplexId, unsigned int dbComplexId, double tmScore, std::string line) : qComplexId(qComplexId), dbComplexId(dbComplexId), tmScore(tmScore), line(line) {}
    OutputLine(unsigned int qComplexId, unsigned int dbComplexId, double tmScore, std::string qComplexName, std::string dbComplexName, std::string qChainNames, std::string dbChainNames): qComplexId(qComplexId), dbComplexId(dbComplexId), tmScore(tmScore), qComplexName(qComplexName), dbComplexName(dbComplexName), qChainNames(qChainNames), dbChainNames(dbChainNames){}
    unsigned int qComplexId;
    unsigned int dbComplexId;
    std::string qComplexName;
    std::string dbComplexName;
    std::string qChainNames;
    std::string dbChainNames;
    double tmScore;
    std::string line;
    void writeLine(unsigned int idx){
         char sep = '\t';
         char nl = '\n';
         line += qComplexName;
         line += sep;
         line += dbComplexName;
         line += sep;
         line += std::to_string(tmScore);
         line += sep;
         line += qChainNames;
         line += sep;
         line += dbChainNames;
         line += sep;
         line += std::to_string(idx);
         line += nl;
    }
};

bool compareComplex(const Complex &first, const Complex &second){
    if (first.complexId < second.complexId) {
        return true;
    }
    if (first.complexId > second.complexId) {
        return false;
    }
    return false;
}

bool compareChainToChainAlnByComplexIdAndChainKey(const ChainToChainAln &first, const ChainToChainAln &second){
    if (first.qChain.complexId < second.qChain.complexId)
        return true;
    if (first.qChain.complexId > second.qChain.complexId)
        return false;
    if (first.dbChain.complexId < second.dbChain.complexId)
        return true;
    if (first.dbChain.complexId > second.dbChain.complexId)
        return false;
    if (first.qChain.chainKey < second.qChain.chainKey)
        return true;
    if (first.qChain.chainKey > second.qChain.chainKey)
        return false;
    if (first.dbChain.chainKey < second.dbChain.chainKey)
        return true;
    if (first.dbChain.chainKey > second.dbChain.chainKey)
        return false;
    return false;
}

bool compareChainToChainAlnByClusterLabel(const ChainToChainAln &first, const ChainToChainAln &second){
    if (first.featureVector.label < second.featureVector.label)
        return true;
    if (first.featureVector.label > second.featureVector.label)
        return false;
    if (first.qChain.chainKey < second.qChain.chainKey)
        return true;
    if (first.qChain.chainKey > second.qChain.chainKey)
        return false;
    if (first.dbChain.chainKey < second.dbChain.chainKey)
        return true;
    if (first.dbChain.chainKey > second.dbChain.chainKey)
        return false;

    return false;
}

bool compareDistAndIndexPairByDist(const DistAndIndexPair &first, const DistAndIndexPair &second){
    if (first.distance < second.distance){
        return true;
    }
    if (first.distance > second.distance){
        return false;
    }
    return false;
}

bool compareOutputLine(const OutputLine &first, const OutputLine &second){
    if (first.qComplexId < second.qComplexId) {
        return true;
    }
    if (first.qComplexId > second.qComplexId) {
        return false;
    }
    if (first.dbComplexId < second.dbComplexId) {
        return true;
    }
    if (first.dbComplexId > second.dbComplexId) {
        return false;
    }
    if (first.tmScore > second.tmScore) {
        return true;
    }
    if (first.tmScore < second.tmScore) {
        return false;
    }
    return false;
}

class DBSCANCluster {
public:
    DBSCANCluster(){
        recursiveNum = 0;
        minClusterSize = 2;
    }
    DBSCANCluster(unsigned int minSize, float defEps){
        recursiveNum = 0;
        minClusterSize = minSize;
        defaultEps = defEps;
    }
    void  clusterAlns(Complex & qComplex, float eps, unsigned int clusterSize) {
        initializeAlnlabeling(qComplex);
        if (++recursiveNum > MAX_RECURSIVE_NUM) return;
        if (clusterSize < minClusterSize) return;
        if (clusterSize==2) return getClusterWithClusterSize2(qComplex);
        int numCluster = 0;
        clearClusterVectors();

        for (size_t i=0; i<qComplex.alnVec.size(); i++) {
            ChainToChainAln &centerAln = qComplex.alnVec[i];
            if (centerAln.featureVector.label != 0) continue;
            std::vector<unsigned int> neighbors = getNeighbors(qComplex, i, eps);
            if (neighbors.size() < MIN_PTS){
                centerAln.featureVector.label = -1;
                continue;
            }
            centerAln.featureVector.label = ++numCluster;
            std::vector<unsigned int> neighborsOfNeighbors = neighbors;
            unsigned int j = 0;
            neighbors.clear();
            while (j < neighborsOfNeighbors.size()) {
                unsigned int q = neighborsOfNeighbors[j++];
                if (i==q) continue;
                ChainToChainAln &neighborAln = qComplex.alnVec[q];
                switch (neighborAln.featureVector.label) {
                    case 0:
                        neighborAln.featureVector.label = numCluster;
                        break;
                    case -1:
                        neighborAln.featureVector.label = numCluster;
                    default:
                        continue;
                }
                std::vector<unsigned int> neighborsOfCurrNeighbor = getNeighbors(qComplex, q, eps);
                if (neighborsOfCurrNeighbor.size() >= MIN_PTS){
                    for (size_t k=0; k < neighborsOfCurrNeighbor.size(); k++){
                        if (std::find(neighborsOfNeighbors.begin(), neighborsOfNeighbors.end(), neighborsOfCurrNeighbor[k]) == neighborsOfNeighbors.end())
                            neighborsOfNeighbors.emplace_back(neighborsOfCurrNeighbor[k]);
                    }
                }
            }
            std::vector<unsigned int> qFoundChainKeys;
            std::vector<unsigned int> dbFoundChainKeys;
            bool isDefectiveCluster = false;
            for (size_t k=0; k < neighborsOfNeighbors.size(); k++) {
                unsigned int s = neighborsOfNeighbors[k];
                ChainToChainAln currAln = qComplex.alnVec[s];
                unsigned int qChainKey = currAln.qChain.chainKey;
                unsigned int dbChainKey = currAln.dbChain.chainKey;
                bool isNewQchainKey = std::find(qFoundChainKeys.begin(), qFoundChainKeys.end(), qChainKey)==qFoundChainKeys.end();
                bool isNewDbChainKey = std::find(dbFoundChainKeys.begin(), dbFoundChainKeys.end(), dbChainKey)==dbFoundChainKeys.end();
                if (isNewQchainKey && isNewDbChainKey) {
                    qFoundChainKeys.emplace_back(qChainKey);
                    dbFoundChainKeys.emplace_back(dbChainKey);
                } else {
                    isDefectiveCluster = true;
                    break;
                }
            }
            if (isDefectiveCluster)
                defectiveClusters.emplace_back(numCluster);
            else if (neighborsOfNeighbors.size() > clusterSize)
                tooBigClusters.emplace_back(numCluster);
            else if (neighborsOfNeighbors.size() < clusterSize)
                tooSmallClusters.emplace_back(numCluster);
            else
                validClusters.emplace_back(numCluster);
        }

        if (!validClusters.empty()){
            keepValidClustersOnly(qComplex);
            std::sort(qComplex.alnVec.begin(), qComplex.alnVec.end(), compareChainToChainAlnByClusterLabel);
            return;
        }
        else if (numCluster == 0) return clusterAlns(qComplex, eps * (1 + LEARNING_RATE), clusterSize);
        else if (!tooSmallClusters.empty()) return clusterAlns(qComplex, eps * (1 + LEARNING_RATE), clusterSize);
        else if (!tooBigClusters.empty()) return clusterAlns(qComplex, eps * (1 - LEARNING_RATE), clusterSize);
        else if (!defectiveClusters.empty()) return clusterAlns(qComplex, defaultEps, clusterSize - 1);
        else return clusterAlns(qComplex, eps * (1 + LEARNING_RATE), clusterSize);
    }

private:
    const unsigned int MAX_RECURSIVE_NUM = 1000;
    const double LEARNING_RATE = 0.05;
    const unsigned int MIN_PTS = 2;
    float defaultEps;
    unsigned int recursiveNum;
    unsigned int minClusterSize;
    std::vector<unsigned int> validClusters;
    std::vector<unsigned int> tooSmallClusters;
    std::vector<unsigned int> tooBigClusters;
    std::vector<unsigned int> defectiveClusters;

    std::vector<unsigned int> getNeighbors(Complex & qComplex, unsigned int centerIdx, float eps){
        ChainToChainAln centerAln = qComplex.alnVec[centerIdx];
        std::vector<unsigned int> neighbors;
        neighbors.emplace_back(centerIdx);
        for (size_t neighborIdx=0; neighborIdx < qComplex.alnVec.size(); neighborIdx++) {
            if (neighborIdx == centerIdx) continue;
            ChainToChainAln neighborAln = qComplex.alnVec[neighborIdx];
            double dist = centerAln.featureVector.getDistance(neighborAln.featureVector);
            if (dist<eps)
                neighbors.emplace_back(neighborIdx);
        }
        return neighbors;
    }

    void getClusterWithClusterSize2(Complex & qComplex) {
        unsigned int C = 0;
        std::vector<DistAndIndexPair> distAndIndexPairVec;
        for (size_t i=0; i<qComplex.alnVec.size(); i++){
            ChainToChainAln prevAln = qComplex.alnVec[i];
            for (size_t j=i+1; j<qComplex.alnVec.size(); j++){
                ChainToChainAln currAln = qComplex.alnVec[j];
                if (qComplex.alnVec[i].qChain.chainKey==qComplex.alnVec[j].qChain.chainKey || qComplex.alnVec[i].dbChain.chainKey==qComplex.alnVec[j].dbChain.chainKey) continue;
                double dist = prevAln.featureVector.getDistance(currAln.featureVector);
                distAndIndexPairVec.emplace_back(DistAndIndexPair(std::pair<unsigned int, unsigned int>(i, j), dist));
            }
        }
        std::sort(distAndIndexPairVec.begin(), distAndIndexPairVec.end(), compareDistAndIndexPairByDist);
        for (size_t i=0; i < distAndIndexPairVec.size(); i++){
            unsigned int alnIdx1 = distAndIndexPairVec[i].indexPair.first;
            unsigned int alnIdx2 = distAndIndexPairVec[i].indexPair.second;
            if (qComplex.alnVec[alnIdx1].featureVector.label > 0 || qComplex.alnVec[alnIdx2].featureVector.label > 0) continue;
            qComplex.alnVec[alnIdx1].featureVector.label = ++C;
            qComplex.alnVec[alnIdx2].featureVector.label = C;
        }
        std::sort(qComplex.alnVec.begin(), qComplex.alnVec.end(), compareChainToChainAlnByClusterLabel);
        return ;
    }

    void initializeAlnlabeling(Complex & qComplex){
        for (size_t i=0; i<qComplex.alnVec.size(); i++){
            ChainToChainAln &P = qComplex.alnVec[i];
            P.featureVector.resetLabel();
        }
    }

    void clearClusterVectors() {
        validClusters.clear();
        tooSmallClusters.clear();
        tooBigClusters.clear();
        defectiveClusters.clear();
    }

    void keepValidClustersOnly(Complex & qComplex) {
        for (size_t i=0; i<qComplex.alnVec.size(); i++){
            ChainToChainAln &P = qComplex.alnVec[i];
            P.featureVector.label = std::find(validClusters.begin(), validClusters.end(), P.featureVector.label) == validClusters.end() ? -1 : P.featureVector.label;
        }
    }
};


class ComplexScorer {
public:
    ComplexScorer(
            IndexReader *qDbr3Di, IndexReader *tDbr3Di, std::string qLookupFile, std::string tLookupFile,
            DBReader<unsigned int> &alnDbr, IndexReader *qCaDbr, IndexReader *tCaDbr, int thread_idx)
            : alnDbr(alnDbr), qCaDbr(qCaDbr), tCaDbr(tCaDbr), thread_idx(thread_idx) {
         maxSeqLen = std::max(qDbr3Di->sequenceReader->getMaxSeqLen()+1, tDbr3Di->sequenceReader->getMaxSeqLen()+1);
        getMaps(qLookupFile, qLookup, qComplexMap, qChainMap);
        getMaps(tLookupFile, tLookup, dbComplexMap, dbChainMap);
    }

    std::vector<Complex> getQComplexes(){
        std::vector<Complex> qTempComplexes = parseInputData();
        std::vector<Complex> qComplexes;
        for (size_t i=0; i<qTempComplexes.size(); i++) {
            Complex qComplex = qTempComplexes[i];
            qComplex.filterAlnVec(1.0);
            if (!qComplex.alnVec.empty()) {
                qComplex.normalize();
                qComplexes.emplace_back(qComplex);
            }
        }
        qTempComplexes.clear();
        return qComplexes;
    }
    std::vector<ComplexToComplexAln> getComplexAlns(Complex qComplex) {
        tmAligner = new TMaligner((unsigned int)(maxSeqLen*qComplex.chainKeys.size()), false);
        DBSCANCluster dbscanCluster = DBSCANCluster(2, 0.5);
//        DBSCANCluster dbscanCluster = DBSCANCluster((unsigned int)std::ceil(qComplex.chainKeys.size()*0.8));
        dbscanCluster.clusterAlns(qComplex, 0.5, qComplex.chainKeys.size());
        int currLabel = 0;
        std::vector<ComplexToComplexAln> complexAlns;
        ComplexToComplexAln complexAln = ComplexToComplexAln(qComplex.complexId, qComplex.alnVec[0].dbChain.complexId);
        for (size_t chainToChainAlnIdx=0; chainToChainAlnIdx < qComplex.alnVec.size(); chainToChainAlnIdx++) {
            ChainToChainAln currAln = qComplex.alnVec[chainToChainAlnIdx];
            if (currAln.featureVector.label == 0 || currAln.featureVector.label == -1) continue;
            if (currAln.featureVector.label == currLabel){
                complexAln.appendChainToChainAln(currAln);
            } else {
                if (currLabel>0){
                    complexAln.setTmScore(getTmScore(complexAln));
                    complexAlns.emplace_back(complexAln);
                }
                complexAln.reset(currAln.dbChain.complexId);
                complexAln.appendChainToChainAln(currAln);
                currLabel = currAln.featureVector.label;
            }
        }
        if (currLabel>0){
            complexAln.setTmScore(getTmScore(complexAln));
            complexAlns.emplace_back(complexAln);
        }
        return complexAlns;
    }

    std::vector<OutputLine>  getOutputLines(std::vector<ComplexToComplexAln> complexAlns){
        std::vector<OutputLine> outputLines;
        if (complexAlns.size()==0) {
            Debug(Debug::WARNING) << "Assignment is not made. Nothing will be returned." << "\n";
            return outputLines;
        }

        for (size_t i=0; i<complexAlns.size(); i++) {
            ComplexToComplexAln complexAln = complexAlns[i];
            std::string qComplexName = qComplexMap.at(complexAln.qComplexId);
            std::string dbComplexName = dbComplexMap.at(complexAln.dbComplexId);
            std::string qChainNames = getChainNanes(complexAln.qChainKeys, qChainMap);
            std::string dbChainNames = getChainNanes(complexAln.dbChainKeys, dbChainMap);
            OutputLine outputLine = OutputLine(complexAln.qComplexId, complexAln.dbComplexId, complexAln.tmScore, qComplexName, dbComplexName, qChainNames, dbChainNames);
            outputLines.emplace_back(outputLine);
        }

        std::sort(outputLines.begin(), outputLines.end(), compareOutputLine);
        unsigned int prevDbComplexId = outputLines[0].dbComplexId;
        unsigned int resultIdx = 0;
        for (size_t outputIdx=0; outputIdx < outputLines.size(); outputIdx++) {
            if (outputLines[outputIdx].dbComplexId!=prevDbComplexId){
                resultIdx = 0;
                prevDbComplexId = outputLines[outputIdx].dbComplexId;
                outputLines[outputIdx].writeLine(resultIdx++);
            } else {
                outputLines[outputIdx].writeLine(resultIdx++);
            }
        }
        return outputLines;
    }

private:
    TMaligner * tmAligner;
    unsigned int maxSeqLen;
    std::map<unsigned int, unsigned int> qLookup;
    std::map<unsigned int, unsigned int> tLookup;
    std::map<unsigned int, std::string> qComplexMap;
    std::map<unsigned int, std::string> dbComplexMap;
    std::map<unsigned int, std::string> qChainMap;
    std::map<unsigned int, std::string> dbChainMap;
    DBReader<unsigned int> &alnDbr;
    IndexReader *qCaDbr;
    IndexReader *tCaDbr;
    Coordinate16 qCoords;
    Coordinate16 tCoords;
    int thread_idx;

    std::vector<Complex> parseInputData() {
        tmAligner = new TMaligner((unsigned int)(maxSeqLen), false);
        std::vector<Complex> qComplexes;
        unsigned int prevComplexId = qLookup.at(0);
        unsigned int complexId;
        std::vector<unsigned int> qTempChainKeys;
        for (size_t chainKey = 0; chainKey < qLookup.size(); chainKey++) {
            complexId = qLookup.at(chainKey);
            if (complexId != prevComplexId) {
                qComplexes.emplace_back(Complex(prevComplexId, qTempChainKeys));
                qTempChainKeys.clear();
            }
            prevComplexId = complexId;
            qTempChainKeys.emplace_back(chainKey);
        }
        qComplexes.emplace_back(Complex(prevComplexId, qTempChainKeys));
        qTempChainKeys.clear();
        std::sort(qComplexes.begin(), qComplexes.end(), compareComplex);

        for (size_t i = 0; i < alnDbr.getSize(); i++) {
            size_t queryKey = alnDbr.getDbKey(i);
            const unsigned queryComplexId = qLookup.at(queryKey);
            char *data = alnDbr.getData(i, thread_idx);
            if (*data == '\0') continue;

            Matcher::result_t qAlnResult = Matcher::parseAlignmentRecord(data);
            size_t qId = qCaDbr->sequenceReader->getId(queryKey);
            char *qCaData = qCaDbr->sequenceReader->getData(qId, thread_idx);
            size_t qCaLength = qCaDbr->sequenceReader->getEntryLen(qId);
            float* queryCaData = qCoords.read(qCaData, qAlnResult.qLen, qCaLength);
            Chain qChain = Chain(queryComplexId, queryKey);
            tmAligner->initQuery(queryCaData, &queryCaData[qAlnResult.qLen], &queryCaData[qAlnResult.qLen*2], NULL, qAlnResult.qLen);
            while (*data != '\0') {
                char dbKeyBuffer[255 + 1];
                Util::parseKey(data, dbKeyBuffer);
                const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                const unsigned int dbComplexId = tLookup.at(dbKey);
                Matcher::result_t alnResult =  Matcher::parseAlignmentRecord(data);
                size_t tCaId = tCaDbr->sequenceReader->getId(dbKey);
                char *tCaData = tCaDbr->sequenceReader->getData(tCaId, thread_idx);
                size_t tCaLength = tCaDbr->sequenceReader->getEntryLen(tCaId);
                float* targetCaData = tCoords.read(tCaData, alnResult.dbLen, tCaLength);
                Chain dbChain = Chain(dbComplexId, dbKey);
                TMaligner::TMscoreResult tmResult = tmAligner->computeTMscore(targetCaData, &targetCaData[alnResult.dbLen], &targetCaData[alnResult.dbLen + alnResult.dbLen], alnResult.dbLen, alnResult.qStartPos, alnResult.dbStartPos, Matcher::uncompressAlignment(alnResult.backtrace));
                ChainToChainAln chainAln(qChain, dbChain, queryCaData, targetCaData, alnResult, tmResult);
                qComplexes[queryComplexId].alnVec.emplace_back(chainAln);
                data = Util::skipLine(data);
            }
            std::sort(qComplexes[queryComplexId].alnVec.begin(), qComplexes[queryComplexId].alnVec.end(), compareChainToChainAlnByComplexIdAndChainKey);
        }
        return qComplexes;
    }

    double getTmScore(ComplexToComplexAln aln){
        bool chainOverlapAllowed = false;
        std::vector<unsigned int> foundDbKeys;
        unsigned int dbKeyIdx = 0;
        while (!chainOverlapAllowed && dbKeyIdx<aln.dbChainKeys.size()) {
            unsigned int currDbKey = aln.dbChainKeys[dbKeyIdx++];
            bool found = std::count(foundDbKeys.begin(), foundDbKeys.end(), currDbKey) > 0;
            if (found) return 0;
            foundDbKeys.emplace_back(currDbKey);
        }
        tmAligner->initQuery(&aln.qCaXVec[0], &aln.qCaYVec[0], &aln.qCaZVec[0], NULL, aln.qLength);
        TMaligner::TMscoreResult tmResult = tmAligner->computeTMscore2(&aln.dbCaXVec[0], &aln.dbCaYVec[0], &aln.dbCaZVec[0], aln.dbLength, 0, 0, Matcher::uncompressAlignment(aln.backtrace), aln.alnLength);
        return tmResult.tmscore;
    }

    static void getMaps(const std::string& file, std::map<unsigned int, unsigned int> &lookupMap, std::map<unsigned int, std::string> &complexNameMap, std::map<unsigned int, std::string> &chainNameMap){
        if (file.length() == 0) return;
        MemoryMapped lookup(file, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
        char* data = (char *) lookup.getData();
        const char* entry[255];
        while (*data != '\0') {
            const size_t columns = Util::getWordsOfLine(data, entry, 255);
            if (columns < 3) {
                Debug(Debug::WARNING) << "Not enough columns in lookup file " << file << "\n";
                continue;
            }
            lookupMap.emplace(Util::fast_atoi<unsigned int>(entry[0]), Util::fast_atoi<unsigned int>(entry[2]));
            data = Util::skipLine(data);
            std::string name(entry[1], data - entry[1] - 1);
            complexNameMap.emplace(Util::fast_atoi<unsigned int>(entry[2]), name.substr(0,name.find_last_of('_')));
            chainNameMap.emplace(Util::fast_atoi<unsigned int>(entry[0]), name.substr(name.find_last_of('_')+1,name.find('\t')-name.find_last_of('_')-1));
        }
        lookup.close();
    }

    static std::string getChainNanes(std::vector<unsigned int> chainKeys, std::map<unsigned int,std::string>map){
        std::string chainListString;
        unsigned int idx = 0;
        while (true) {
            chainListString += map.at(chainKeys[idx]);
            idx++;
            if (idx < chainKeys.size()){ chainListString += ',';} else {break;}
        }
        return chainListString;
    }
};

int scorecomplex(int argc, const char **argv, const Command& command) {

    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    IndexReader q3DiDbr(StructureUtil::getIndexWithSuffix(par.db1, "_ss"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
    IndexReader *t3DiDbr = NULL;
    IndexReader *qCaDbr = new IndexReader(par.db1, par.threads, IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY), touch ? IndexReader::PRELOAD_INDEX : 0,  DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA, "_ca" );
    IndexReader *tCaDbr = NULL;

    bool sameDB = false;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        t3DiDbr = &q3DiDbr;
        tCaDbr = qCaDbr;
    } else {
        t3DiDbr = new IndexReader(StructureUtil::getIndexWithSuffix(par.db2, "_ss"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
        tCaDbr = new IndexReader(par.db2, par.threads, IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY), touch ? IndexReader::PRELOAD_INDEX : 0, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA, "_ca");
    }

    std::string qLookupFile = par.db1 + ".lookup";
    std::string dbLookupFile = par.db2 + ".lookup";

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    resultWriter.open();
    Debug::Progress progress(alnDbr.getSize());

//#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::vector<Matcher::result_t> alignmentResult;
        std::string backtrace;
        std::string resultBuffer;
        ComplexScorer complexScorer(&q3DiDbr, t3DiDbr, qLookupFile, dbLookupFile, alnDbr, qCaDbr, tCaDbr, thread_idx);
        std::vector<Complex> qComplexes = complexScorer.getQComplexes();
//#pragma omp for schedule(dynamic, 1)

        for (size_t qComplexIdx=0; qComplexIdx < qComplexes.size(); qComplexIdx++) {
//            for (size_t i=0; i<qComplexes[qComplexIdx].alnVec.size(); i++){
//                std::cout << qComplexes[qComplexIdx].alnVec[i].qChain.chainKey << "\t" << qComplexes[qComplexIdx].alnVec[i].dbChain.chainKey << "\t";
//                for (size_t j=0; j<7; j++) {
//                    std::cout << qComplexes[qComplexIdx].alnVec[i].featureVector.features[j] << "\t";
//                }
//                std::cout<< std::endl;
//            }
            std::vector<ComplexToComplexAln> complexAlns = complexScorer.getComplexAlns(qComplexes[qComplexIdx]);
            std::vector<OutputLine> resultLines = complexScorer.getOutputLines(complexAlns);
            progress.updateProgress();
            for (size_t resultIdx=0; resultIdx < resultLines.size(); resultIdx++) {
                resultWriter.writeData(resultLines[resultIdx].line.c_str(), resultLines[resultIdx].line.size(), 0, thread_idx);
            }
        }
    }
    alnDbr.close();
    if (!sameDB) {
        delete t3DiDbr;
    }
    resultWriter.close(true);
    return EXIT_SUCCESS;
}
