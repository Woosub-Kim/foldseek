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

    void setCaData(std::vector<float> &caX, std::vector<float> &caY, std::vector<float> &caZ) {
        caVecX = caX;
        caVecY = caY;
        caVecZ = caZ;
    }
};

struct FeatureVector{
    FeatureVector(){}
    FeatureVector(float u[3][3], float t[3]) {
        for (size_t i=0; i<3; i++) {
            for (size_t j=0; j<3; j++) {
                features[3*i+j] =  u[i][j];
            }
            features[9+i] = t[i];
        }
        label = 0;
    }
    double features[12];
    int label;

    double getDistance(const FeatureVector &o){
        double dist = 0;
        for (size_t i=0; i<12; i++){
            dist += std::pow(features[i] - o.features[i], 2);
        }
        dist = std::sqrt(dist);
        return dist;
    }
};

struct ChainToChainAln {
    ChainToChainAln() {}
    ChainToChainAln(Chain queryChain, Chain targetChain, float * queryCaData, float * targetCaData, Matcher::result_t alnResult, TMaligner::TMscoreResult tmResult)
    : qChain(queryChain), dbChain(targetChain)  {
        std::vector<float> qCaXVec;
        std::vector<float> qCaYVec;
        std::vector<float> qCaZVec;
        std::vector<float> dbCaXVec;
        std::vector<float> dbCaYVec;
        std::vector<float> dbCaZVec;
        std::string newBacktrace;
        unsigned int qPos = alnResult.qStartPos;
        unsigned int dbPos = alnResult.dbStartPos;
        unsigned int qXPos = 0;
        unsigned int qYPos = alnResult.qLen;
        unsigned int qZPos = alnResult.qLen*2;
        unsigned int dbXPos = 0;
        unsigned int dbYPos = alnResult.dbLen;
        unsigned int dbZPos = alnResult.dbLen*2;
        for (char cigar : alnResult.backtrace){
            switch (cigar) {
                case 'M':
                    newBacktrace += "M";
                    qCaXVec.emplace_back(queryCaData[qXPos + qPos]);
                    qCaYVec.emplace_back(queryCaData[qYPos + qPos]);
                    qCaZVec.emplace_back(queryCaData[qZPos + qPos++]);
                    dbCaXVec.emplace_back(targetCaData[dbXPos + dbPos]);
                    dbCaYVec.emplace_back(targetCaData[dbYPos + dbPos]);
                    dbCaZVec.emplace_back(targetCaData[dbZPos + dbPos++]);
                    break;
                case 'I':
                    qPos++;
                    break;
                case 'D':
                    dbPos++;
                    break;
//                default:
//                    Debug(Debug::WARNING) << "backtrace ???" << "\n";
//                    break;
            }
        }
        qChain.setCaData(qCaXVec, qCaYVec, qCaZVec);
        dbChain.setCaData(dbCaXVec, dbCaYVec, dbCaZVec);
        qChain.startPos = 0;
        dbChain.startPos = 0;
//        qChain.length = newBacktrace.length();
//        dbChain.length = newBacktrace.length();
//        alnLength = newBacktrace.length();
        qChain.length = alnResult.qLen;
        dbChain.length = alnResult.dbLen;
        alnLength = alnResult.alnLength;
        backtrace = newBacktrace;
        featureVector = FeatureVector(tmResult.u, tmResult.t);
    }
    Chain qChain;
    Chain dbChain;
    std::string backtrace;
    unsigned int alnLength;
    FeatureVector featureVector;
};

struct IndexPairWithDist{
    IndexPairWithDist(std::pair<unsigned int, unsigned int> indexPair, double distance) : indexPair(indexPair), distance(distance){}
    std::pair<unsigned int, unsigned int> indexPair;
    double distance;
};

struct Complex {
    Complex() {}
    Complex(unsigned int complexId, std::vector<unsigned int> &chainKeys) : complexId(complexId), chainKeys(chainKeys), alnVec({}) {}
    unsigned int complexId;
    std::vector<unsigned int> chainKeys;
    std::vector<ChainToChainAln> alnVec;

    void filterAlnVec(float compatibleCheckRatio){
        if (alnVec.empty())
            return;
        std::vector<ChainToChainAln> newAlnVec;
        std::vector<ChainToChainAln> currDbComplexAlnVec;
        ChainToChainAln firstAln = alnVec[0];
        unsigned int dbPrevComplexId = firstAln.dbChain.complexId;
        unsigned int qChainCount = 1;
        unsigned int qPrevChainKey = firstAln.qChain.chainKey;

        for (auto aln : alnVec) {
            if (aln.dbChain.complexId != dbPrevComplexId) {
                if (qChainCount >= (unsigned int)(chainKeys.size() * compatibleCheckRatio))
                    newAlnVec.insert(newAlnVec.end(), currDbComplexAlnVec.begin(), currDbComplexAlnVec.end());
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
        if (qChainCount >= (unsigned int)(chainKeys.size() * compatibleCheckRatio))
            newAlnVec.insert(newAlnVec.end(), currDbComplexAlnVec.begin(), currDbComplexAlnVec.end());
        alnVec = newAlnVec;
    }

    void normalize(){
        unsigned int length = alnVec.size();
        double mean[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
        double var[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
        double sd[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
        double cv[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
        for (size_t i=0; i < length; i++){
            for (size_t j=0; j<12; j++){
                mean[j] += alnVec[i].featureVector.features[j]/(double)length;
            }
        }
        for (size_t i=0; i < length; i++){
            for (size_t j=0; j<12; j++){
                double value = alnVec[i].featureVector.features[j];
                value -= mean[j];
                value *= value;
                var[j] += value/(double)length;
            }
        }
        for (size_t j=0; j<12; j++){
            sd[j] = std::sqrt(var[j]);
            cv[j] = abs(mean[j]) > 1.0 ? sd[j]/std::abs(mean[j]) : sd[j];
        }
        for (size_t i=0; i < length; i++){
            for (size_t j=0; j<12; j++){
                alnVec[i].featureVector.features[j] = (cv[j] < 0.1) ? 0.0 : (alnVec[i].featureVector.features[j] - mean[j]) / sd[j];
            }
        }
    }
};

struct ComplexToComplexAln{
    ComplexToComplexAln() {}
    ComplexToComplexAln(unsigned int qComplexId, unsigned int dbComplexId): qComplexId(qComplexId), dbComplexId(dbComplexId){}
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
    OutputLine(unsigned int qComplexId, unsigned int dbComplexId, double tmScore, std::string &qComplexName, std::string &dbComplexName, std::string &qChainNames, std::string &dbChainNames): qComplexId(qComplexId), dbComplexId(dbComplexId), tmScore(tmScore), qComplexName(qComplexName), dbComplexName(dbComplexName), qChainNames(qChainNames), dbChainNames(dbChainNames){}
    unsigned int qComplexId;
    unsigned int dbComplexId;
    double tmScore;
    std::string qComplexName;
    std::string dbComplexName;
    std::string qChainNames;
    std::string dbChainNames;
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
    if (first.complexId < second.complexId)
        return true;
    if (first.complexId > second.complexId)
        return false;
    return false;
}

bool compareChainToChainAlnByDbComplexId(const ChainToChainAln &first, const ChainToChainAln &second){
    if (first.dbChain.complexId < second.dbChain.complexId)
        return true;
    if (first.dbChain.complexId > second.dbChain.complexId)
        return false;
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

bool compareIndexPairWithDistByDist(const IndexPairWithDist &first, const IndexPairWithDist &second){
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
    DBSCANCluster(unsigned int minSize, float defEps){
        recursiveNum = 0;
        minClusterSize = minSize;
        defaultEps = defEps;
    }
    void  clusterAlns(Complex & qComplex, float eps, unsigned int clusterSize) {
        initializeAlnLabeling(qComplex);
        if (++recursiveNum > MAX_RECURSIVE_NUM)
            return;
        if (clusterSize < minClusterSize)
            return;
        if (clusterSize==2)
            return getClusterWithClusterSize2(qComplex);
        int cLabel = 0;
        clearClusterVectors();

        for (size_t i=0; i<qComplex.alnVec.size(); i++) {
            ChainToChainAln &centerAln = qComplex.alnVec[i];
            if (centerAln.featureVector.label != 0) continue;
            std::vector<unsigned int> neighbors = getNeighbors(qComplex, i, eps);
            if (neighbors.size() < MIN_PTS){
                centerAln.featureVector.label = -1;
                continue;
            }
            centerAln.featureVector.label = ++cLabel;
            std::vector<unsigned int> neighborsOfNeighbors = neighbors;
            unsigned int j = 0;
            neighbors.clear();
            while (j < neighborsOfNeighbors.size()) {
                unsigned int q = neighborsOfNeighbors[j++];
                if (i==q) continue;
                ChainToChainAln &neighborAln = qComplex.alnVec[q];
                switch (neighborAln.featureVector.label) {
                    case 0:
                        neighborAln.featureVector.label = cLabel;
                        break;
                    case -1:
                        neighborAln.featureVector.label = cLabel;
                    default:
                        continue;
                }
                std::vector<unsigned int> neighborsOfCurrNeighbor = getNeighbors(qComplex, q, eps);
                if (neighborsOfCurrNeighbor.size() >= MIN_PTS){
                    for (unsigned int & neighbor : neighborsOfCurrNeighbor){
                        if (std::find(neighborsOfNeighbors.begin(), neighborsOfNeighbors.end(), neighbor) == neighborsOfNeighbors.end())
                            neighborsOfNeighbors.emplace_back(neighbor);
                    }
                }
            }
            std::vector<unsigned int> qFoundChainKeys;
            std::vector<unsigned int> dbFoundChainKeys;
            bool isDefectiveCluster = false;
            for (auto neighborIdx : neighborsOfNeighbors) {
                ChainToChainAln currAln = qComplex.alnVec[neighborIdx];
                unsigned int qChainKey = currAln.qChain.chainKey;
                unsigned int dbChainKey = currAln.dbChain.chainKey;
                bool isNewQChainKey = std::find(qFoundChainKeys.begin(), qFoundChainKeys.end(), qChainKey) == qFoundChainKeys.end();
                bool isNewDbChainKey = std::find(dbFoundChainKeys.begin(), dbFoundChainKeys.end(), dbChainKey)==dbFoundChainKeys.end();
                if (isNewQChainKey && isNewDbChainKey) {
                    qFoundChainKeys.emplace_back(qChainKey);
                    dbFoundChainKeys.emplace_back(dbChainKey);
                } else {
                    isDefectiveCluster = true;
                    break;
                }
            }
            if (isDefectiveCluster)
                defectiveClusters.emplace_back(cLabel);
            else if (neighborsOfNeighbors.size() > clusterSize)
                bigClusters.emplace_back(cLabel);
            else if (neighborsOfNeighbors.size() < clusterSize)
                smallClusters.emplace_back(cLabel);
            else
                validClusters.emplace_back(cLabel);
        }

        if (!validClusters.empty()){
            keepValidClustersOnly(qComplex);
            std::sort(qComplex.alnVec.begin(), qComplex.alnVec.end(), compareChainToChainAlnByClusterLabel);
            return;
        }
        else if (cLabel==0 || !smallClusters.empty()) return clusterAlns(qComplex, eps * (1 + LEARNING_RATE), clusterSize);
        else if (!bigClusters.empty()) return clusterAlns(qComplex, eps * (1 - LEARNING_RATE), clusterSize);
        else if (!defectiveClusters.empty()) return clusterAlns(qComplex, defaultEps, clusterSize - 1);
        else return clusterAlns(qComplex, eps * (1 + LEARNING_RATE), clusterSize);
    }

private:
    const unsigned int MAX_RECURSIVE_NUM = 1000;
    const float LEARNING_RATE = 0.05;
    const unsigned int MIN_PTS = 2;
    float defaultEps;
    unsigned int recursiveNum;
    unsigned int minClusterSize;
    std::vector<unsigned int> validClusters;
    std::vector<unsigned int> smallClusters;
    std::vector<unsigned int> bigClusters;
    std::vector<unsigned int> defectiveClusters;

    static std::vector<unsigned int> getNeighbors(Complex & qComplex, unsigned int centerIdx, float eps){
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

    static void getClusterWithClusterSize2(Complex & qComplex) {
        int cLabel = 0;
        std::vector<IndexPairWithDist> IndexPairs;
        for (size_t i=0; i<qComplex.alnVec.size(); i++){
            ChainToChainAln prevAln = qComplex.alnVec[i];
            for (size_t j=i+1; j<qComplex.alnVec.size(); j++){
                ChainToChainAln currAln = qComplex.alnVec[j];
                if (qComplex.alnVec[i].qChain.chainKey==qComplex.alnVec[j].qChain.chainKey || qComplex.alnVec[i].dbChain.chainKey==qComplex.alnVec[j].dbChain.chainKey) continue;
                double dist = prevAln.featureVector.getDistance(currAln.featureVector);
                IndexPairs.emplace_back(IndexPairWithDist(std::pair<unsigned int, unsigned int>(i, j), dist));
            }
        }
        std::sort(IndexPairs.begin(), IndexPairs.end(), compareIndexPairWithDistByDist);
        for (auto & IndexPair : IndexPairs){
            unsigned int alnIdx1 = IndexPair.indexPair.first;
            unsigned int alnIdx2 = IndexPair.indexPair.second;
            if (qComplex.alnVec[alnIdx1].featureVector.label > 0 || qComplex.alnVec[alnIdx2].featureVector.label > 0) continue;
            qComplex.alnVec[alnIdx1].featureVector.label = ++cLabel;
            qComplex.alnVec[alnIdx2].featureVector.label = cLabel;
        }
        std::sort(qComplex.alnVec.begin(), qComplex.alnVec.end(), compareChainToChainAlnByClusterLabel);
    }

    static void initializeAlnLabeling(Complex & qComplex){
        for (auto & aln : qComplex.alnVec){
            aln.featureVector.label=0;
        }
    }

    void clearClusterVectors() {
        validClusters.clear();
        smallClusters.clear();
        bigClusters.clear();
        defectiveClusters.clear();
    }

    void keepValidClustersOnly(Complex & qComplex) {
        for (auto & aln : qComplex.alnVec){
            aln.featureVector.label = std::find(validClusters.begin(), validClusters.end(), aln.featureVector.label) == validClusters.end() ? -1 : aln.featureVector.label;
        }
    }
};


class ComplexScorer {
public:
    ComplexScorer(
            IndexReader *qDbr3Di, IndexReader *tDbr3Di, const std::string& qLookupFile, const std::string& tLookupFile,
            DBReader<unsigned int> &alnDbr, IndexReader *qCaDbr, IndexReader *tCaDbr, unsigned int thread_idx, float minAssignedChainsRatio
            ): alnDbr(alnDbr), qCaDbr(qCaDbr), tCaDbr(tCaDbr), thread_idx(thread_idx), minAssignedChainsRatio(minAssignedChainsRatio) {
        maxSeqLen = std::max(qDbr3Di->sequenceReader->getMaxSeqLen()+1, tDbr3Di->sequenceReader->getMaxSeqLen()+1);
        getMaps(qLookupFile, qLookup, qComplexMap, qChainMap);
        getMaps(tLookupFile, tLookup, dbComplexMap, dbChainMap);

    }

    std::vector<Complex> getQComplexes(){
        tmAligner = new TMaligner((unsigned int)(maxSeqLen), false);
        std::vector<Complex> qComplexes;
        std::vector<Complex> qTempComplexes;
        unsigned int prevComplexId = qLookup.at(0);
        unsigned int complexId;
        std::vector<unsigned int> qTempChainKeys;
        for (size_t chainKey = 0; chainKey < qLookup.size(); chainKey++) {
            complexId = qLookup.at(chainKey);
            if (complexId != prevComplexId) {
                qTempComplexes.emplace_back(Complex(prevComplexId, qTempChainKeys));
                qTempChainKeys.clear();
            }
            prevComplexId = complexId;
            qTempChainKeys.emplace_back(chainKey);
        }
        qTempComplexes.emplace_back(Complex(prevComplexId, qTempChainKeys));
        qTempChainKeys.clear();
        std::sort(qTempComplexes.begin(), qTempComplexes.end(), compareComplex);

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
                const auto dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                const unsigned int dbComplexId = tLookup.at(dbKey);
                Matcher::result_t alnResult =  Matcher::parseAlignmentRecord(data);
                size_t tCaId = tCaDbr->sequenceReader->getId(dbKey);
                char *tCaData = tCaDbr->sequenceReader->getData(tCaId, thread_idx);
                size_t tCaLength = tCaDbr->sequenceReader->getEntryLen(tCaId);
                float* targetCaData = tCoords.read(tCaData, alnResult.dbLen, tCaLength);
                Chain dbChain = Chain(dbComplexId, dbKey);
                TMaligner::TMscoreResult tmResult = tmAligner->computeTMscore(targetCaData, &targetCaData[alnResult.dbLen], &targetCaData[alnResult.dbLen + alnResult.dbLen], alnResult.dbLen, alnResult.qStartPos, alnResult.dbStartPos, Matcher::uncompressAlignment(alnResult.backtrace));
                ChainToChainAln chainAln(qChain, dbChain, queryCaData, targetCaData, alnResult, tmResult);
                qTempComplexes[queryComplexId].alnVec.emplace_back(chainAln);
                data = Util::skipLine(data);

            }
            std::sort(qTempComplexes[queryComplexId].alnVec.begin(), qTempComplexes[queryComplexId].alnVec.end(), compareChainToChainAlnByDbComplexId);
        }
        for (size_t queryComplexId=0; queryComplexId < qTempComplexes.size(); queryComplexId++) {
            Complex qComplex = Complex(qTempComplexes[queryComplexId].complexId, qTempComplexes[queryComplexId].chainKeys);
            unsigned int currDbComplexId = qTempComplexes[queryComplexId].alnVec[0].dbChain.complexId;
            for (auto& aln : qTempComplexes[queryComplexId].alnVec) {
                if (aln.dbChain.complexId!=currDbComplexId){
                    qComplex.filterAlnVec(1.0);
                    if (!qComplex.alnVec.empty()) {
                        qComplex.normalize();
                        qComplexes.emplace_back(qComplex);
                    }
                    currDbComplexId = aln.dbChain.complexId;
                    qComplex.alnVec.clear();
                }
                qComplex.alnVec.emplace_back(aln);
            }
            qComplex.filterAlnVec(1.0);
            if (!qComplex.alnVec.empty()) {
                qComplex.normalize();
                qComplexes.emplace_back(qComplex);
            }
            std::sort(qComplexes[queryComplexId].alnVec.begin(), qComplexes[queryComplexId].alnVec.end(), compareChainToChainAlnByComplexIdAndChainKey);
        }
        delete tmAligner;
        return qComplexes;
    }



    std::vector<ComplexToComplexAln> getComplexAlns(Complex qComplex) {
        tmAligner = new TMaligner((unsigned int)(maxSeqLen*qComplex.chainKeys.size()), false);
        DBSCANCluster dbscanCluster = DBSCANCluster((unsigned int)std::ceil(qComplex.chainKeys.size()*minAssignedChainsRatio), 0.5);
        dbscanCluster.clusterAlns(qComplex, 0.5, qComplex.chainKeys.size());
        int currLabel = 0;
        std::vector<ComplexToComplexAln> complexAlns;
        ComplexToComplexAln complexAln = ComplexToComplexAln(qComplex.complexId, qComplex.alnVec[0].dbChain.complexId);
        for (const auto& currAln : qComplex.alnVec) {
            if (currAln.featureVector.label == 0 || currAln.featureVector.label == -1) continue;
            if (currAln.featureVector.label == currLabel){
                complexAln.appendChainToChainAln(currAln);
            } else {
                if (currLabel>0){
                    complexAln.tmScore = getTmScore(complexAln);
                    complexAlns.emplace_back(complexAln);
                }
                complexAln.reset(currAln.dbChain.complexId);
                complexAln.appendChainToChainAln(currAln);
                currLabel = currAln.featureVector.label;
            }
        }
        if (currLabel>0){
            complexAln.tmScore = getTmScore(complexAln);
            complexAlns.emplace_back(complexAln);
        }
        delete tmAligner;
        return complexAlns;
    }

    std::vector<OutputLine>  getOutputLines(const std::vector<ComplexToComplexAln>& complexAlns){
        std::vector<OutputLine> outputLines;
        if (complexAlns.empty()) {
            Debug(Debug::WARNING) << "Any Assignment is not made. Nothing will be returned." << "\n";
            return outputLines;
        }

        for (const auto& aln : complexAlns) {
            std::string qComplexName = qComplexMap.at(aln.qComplexId);
            std::string dbComplexName = dbComplexMap.at(aln.dbComplexId);
            std::string qChainNames = getChainNames(aln.qChainKeys, qChainMap);
            std::string dbChainNames = getChainNames(aln.dbChainKeys, dbChainMap);
            OutputLine outputLine = OutputLine(aln.qComplexId, aln.dbComplexId, aln.tmScore, qComplexName, dbComplexName, qChainNames, dbChainNames);
            outputLines.emplace_back(outputLine);
        }

        std::sort(outputLines.begin(), outputLines.end(), compareOutputLine);
        unsigned int prevDbComplexId = outputLines[0].dbComplexId;
        unsigned int resultIdx = 0;
        for (auto & outputLine : outputLines) {
            if (outputLine.dbComplexId!=prevDbComplexId){
                resultIdx = 0;
                prevDbComplexId = outputLine.dbComplexId;
                outputLine.writeLine(resultIdx++);
            } else {
                outputLine.writeLine(resultIdx++);
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
    unsigned int thread_idx;
    float minAssignedChainsRatio;

    double getTmScore(ComplexToComplexAln aln){
        unsigned int foo = aln.backtrace.length();
        tmAligner->initQuery(&aln.qCaXVec[0], &aln.qCaYVec[0], &aln.qCaZVec[0], NULL, foo);
        TMaligner::TMscoreResult tmResult = tmAligner->computeTMscore2(&aln.dbCaXVec[0], &aln.dbCaYVec[0], &aln.dbCaZVec[0], foo, 0, 0, Matcher::uncompressAlignment(aln.backtrace), foo);
        return tmResult.tmscore*foo/aln.alnLength;
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

    static std::string getChainNames(std::vector<unsigned int> chainKeys, std::map<unsigned int,std::string>map){
        std::string chainListString;
//        char sep = ',';
        unsigned int idx = 0;
        while (true) {
            chainListString += map.at(chainKeys[idx]);
            idx++;
            if (idx < chainKeys.size()){
                chainListString += ',';
            } else {
                break;
            }
        }
//        for (auto chainKey : chainKeys) {
//            std::string chainName = map.at(chainKey);
//        }

        return chainListString;
    }
};

int scorecomplex(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    IndexReader q3DiDbr(StructureUtil::getIndexWithSuffix(par.db1, "_ss"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
    IndexReader *t3DiDbr = NULL;
    auto *qCaDbr = new IndexReader(par.db1, par.threads, IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY), touch ? IndexReader::PRELOAD_INDEX : 0,  DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA, "_ca" );
    IndexReader *tCaDbr = NULL;

    bool sameDB = false;
    if (par.db1 == par.db2) {
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
    float minAssignedChainsRatio = par.minAssignedChainsThreshold/100;
    minAssignedChainsRatio = minAssignedChainsRatio > 1.0 ? 1.0 : minAssignedChainsRatio;


#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::vector<Matcher::result_t> alignmentResult;
        std::string backtrace;
        std::string resultBuffer;

        ComplexScorer complexScorer(&q3DiDbr, t3DiDbr, qLookupFile, dbLookupFile, alnDbr, qCaDbr, tCaDbr, thread_idx, minAssignedChainsRatio);
        std::vector<Complex> qComplexes = complexScorer.getQComplexes();
#pragma omp for schedule(dynamic, 1)
        for (size_t qComplexId=0; qComplexId< qComplexes.size(); qComplexId++) {
            Complex & qComplex = qComplexes[qComplexId];
//            for (size_t i=0; i < qComplex.alnVec.size(); i++){
//                std::cout << qComplex.alnVec[i].qChain.chainKey << "\t" << qComplex.alnVec[i].dbChain.chainKey << "\t";
//                for (size_t j=0; j<12; j++) {
//                    std::cout << qComplex.alnVec[i].featureVector.features[j] << "\t";
//                }
//                std::cout<< std::endl;
//            }
            std::vector<ComplexToComplexAln> complexAlns = complexScorer.getComplexAlns(qComplex);
            std::vector<OutputLine> resultLines = complexScorer.getOutputLines(complexAlns);
            progress.updateProgress();
            for (auto & resultLine : resultLines) {
                resultWriter.writeData(resultLine.line.c_str(), resultLine.line.size(), 0, thread_idx);
            }
        }
    }
    alnDbr.close();
    delete qCaDbr;
    if (!sameDB) {
        delete t3DiDbr;
        delete tCaDbr;
    }
    resultWriter.close(true);
    return EXIT_SUCCESS;
}
