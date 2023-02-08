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
    std::vector<float> caDataX;
    std::vector<float> caDataY;
    std::vector<float> caDataZ;
    int startPos;
    unsigned int length;

    void setLength(unsigned int len){
        length = len;
    }

    void setCaData(std::vector<float> caX, std::vector<float> caY, std::vector<float> caZ){
        caDataX = caX;
        caDataY = caY;
        caDataZ = caZ;
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

struct ChainToChainAln {
    ChainToChainAln() {}
    ChainToChainAln(Chain qChain, Chain dbChain, std::string backtrace): qChain(qChain), dbChain(dbChain), backtrace(backtrace) {}
    ChainToChainAln (Chain qInputChain, Chain dbInputChain, float * queryCaData, float * targetCaData, Matcher::result_t alnResult) {
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
        int numMatches = 0;
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
                    numMatches++;
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

        qInputChain.setStartPos(0); // 0,0
        qInputChain.setCaData(qCaXVec, qCaYVec, qCaZVec);
//        qInputChain.setLength(numMatches);
        qInputChain.setLength(qLength);
        dbInputChain.setStartPos(0); // 0,0
        dbInputChain.setCaData(dbCaXVec, dbCaYVec, dbCaZVec);
//        dbInputChain.setLength(numMatches);
        dbInputChain.setLength(dbLength);
        qChain = qInputChain;
        dbChain = dbInputChain;
        backtrace = newBacktrace;
        bitScore = alnResult.score;
        alnLength = inputBacktrace.size();
    }
    Chain qChain;
    Chain dbChain;
    std::string backtrace;
    unsigned int alnLength;
    int bitScore;
};

struct Complex {
    Complex() {}
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

//    void filterAlnVec(float compatibleCheckRatio, float alnFilterRatio){
//        if (alnVec.empty()) {
//            return;
//        }
//        std::vector<ChainToChainAln> newAlnVec;
//        std::vector<ChainToChainAln> currDbComplexAlnVec;
//        std::vector<ChainToChainAln> tempAlnVec;
//
//        ChainToChainAln firstAln = alnVec[0];
//        unsigned int dbPrevComplexId = firstAln.dbChain.complexId;
//        unsigned int qChainCount = 1;
//        unsigned int qPrevChainKey = firstAln.qChain.chainKey;
//        int maxBitScore = firstAln.bitScore;
//
//        for (size_t i=0; i < alnVec.size(); i++) {
//            ChainToChainAln aln = alnVec[i];
//            if (aln.dbChain.complexId != dbPrevComplexId) {
//                for (size_t j = 0; j < tempAlnVec.size(); j++) {
//                    ChainToChainAln tempAln = tempAlnVec[j];
//                    if (tempAln.bitScore > (int)maxBitScore*alnFilterRatio) {
//                        currDbComplexAlnVec.emplace_back(tempAln);
//                    }
//                }
//                tempAlnVec.clear();
//
//                if (qChainCount >= (float)(chainKeys.size() * compatibleCheckRatio)) {
//                    newAlnVec.insert(newAlnVec.end(), currDbComplexAlnVec.begin(), currDbComplexAlnVec.end());
//                }
//                currDbComplexAlnVec.clear();
//
//                dbPrevComplexId = aln.dbChain.complexId;
//                qChainCount = 1;
//                qPrevChainKey = aln.qChain.chainKey;
//                maxBitScore = aln.bitScore;
//                tempAlnVec.emplace_back(aln);
//            } else if (aln.qChain.chainKey != qPrevChainKey) {
//                for (size_t j = 0; j < tempAlnVec.size(); j++) {
//                    ChainToChainAln tempAln = tempAlnVec[j];
//                    if (tempAln.bitScore > (int)maxBitScore*alnFilterRatio) {
//                        currDbComplexAlnVec.emplace_back(tempAln);
//                    }
//                }
//                tempAlnVec.clear();
//
//                qChainCount ++;
//                qPrevChainKey = aln.qChain.chainKey;
//                maxBitScore = aln.bitScore;
//                tempAlnVec.emplace_back(aln);
//            } else {
//                maxBitScore = std::max(maxBitScore, aln.bitScore);
//                tempAlnVec.emplace_back(aln);
//            }
//        }
//        for (size_t j = 0; j < tempAlnVec.size(); j++) {
//            ChainToChainAln tempAln = tempAlnVec[j];
//            if (tempAln.bitScore > (int)maxBitScore*alnFilterRatio) {
//                currDbComplexAlnVec.emplace_back(tempAln);
//            }
//        }
//        tempAlnVec.clear();
//
//        if (qChainCount >= (float)(chainKeys.size() * compatibleCheckRatio)) {
//            newAlnVec.insert(newAlnVec.end(), currDbComplexAlnVec.begin(), currDbComplexAlnVec.end());
//        }
//        alnVec = newAlnVec;
//    }
};

struct ComplexToComplexAln{
    ComplexToComplexAln() {}
    ComplexToComplexAln(unsigned int qComplexId, unsigned int dbComplexId, double tmScore): qComplexId(qComplexId), dbComplexId(dbComplexId), tmScore(tmScore){}
    ComplexToComplexAln(ChainToChainAln aln) {
        qComplexId = aln.qChain.complexId;
        dbComplexId = aln.dbChain.complexId;
        qChainKeys = {aln.qChain.chainKey};
        dbChainKeys = {aln.dbChain.chainKey};
        backtrace = aln.backtrace;
        qLength = aln.qChain.length;
        dbLength = aln.dbChain.length;
        alnLength = aln.alnLength;
        qCaXVec = aln.qChain.caDataX;
        qCaYVec = aln.qChain.caDataY;
        qCaZVec = aln.qChain.caDataZ;
        dbCaXVec = aln.dbChain.caDataX;
        dbCaYVec = aln.dbChain.caDataY;
        dbCaZVec = aln.dbChain.caDataZ;
    }
    unsigned int qComplexId;
    unsigned int dbComplexId;
    double tmScore;
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

    void setTmScore(double tm){
        tmScore = tm;
    }

    void appendChainToChainAln(ChainToChainAln aln){
        qChainKeys.emplace_back(aln.qChain.chainKey);
        dbChainKeys.emplace_back(aln.dbChain.chainKey);
        backtrace += aln.backtrace;
        qLength += aln.qChain.length;
        dbLength += aln.dbChain.length;
        alnLength += aln.alnLength;
        qCaXVec.insert(qCaXVec.end(), aln.qChain.caDataX.begin(), aln.qChain.caDataX.end());
        qCaYVec.insert(qCaYVec.end(), aln.qChain.caDataY.begin(), aln.qChain.caDataY.end());
        qCaZVec.insert(qCaZVec.end(), aln.qChain.caDataZ.begin(), aln.qChain.caDataZ.end());
        dbCaXVec.insert(dbCaXVec.end(), aln.dbChain.caDataX.begin(), aln.dbChain.caDataX.end());
        dbCaYVec.insert(dbCaYVec.end(), aln.dbChain.caDataY.begin(), aln.dbChain.caDataY.end());
        dbCaZVec.insert(dbCaZVec.end(), aln.dbChain.caDataZ.begin(), aln.dbChain.caDataZ.end());
    }
};

struct OutputLine{
    OutputLine(unsigned int qComplexId, double tmScore, std::string line) : qComplexId(qComplexId), tmScore(tmScore), line(line) {}
    OutputLine(unsigned int qComplexId, double tmScore, std::vector<std::string> fields) : qComplexId(qComplexId), tmScore(tmScore) {
        std::string separatedLine;
        unsigned int i=0;
        while (true) {
            separatedLine.append(fields[i]);
            i++;
            if (i<fields.size()){
                separatedLine.append("\t");
            } else {
                separatedLine.append("\n");
                break;
            }
        }
        line = separatedLine;
    }
    unsigned int qComplexId;
    double tmScore;
    std::string line;
    unsigned int idx;
};

bool compareChainToChainAlnByQChainDbChain(const ChainToChainAln &first, const ChainToChainAln &second){
    if (first.qChain.complexId < second.qChain.complexId)
        return true;
    if (first.qChain.complexId > second.qChain.complexId)
        return false;
    if (first.qChain.chainKey < second.qChain.chainKey)
        return true;
    if (first.qChain.chainKey > second.qChain.chainKey)
        return false;
    if (first.dbChain.complexId < second.dbChain.complexId)
        return true;
    if (first.dbChain.complexId > second.dbChain.complexId)
        return false;
    if (first.dbChain.chainKey < second.dbChain.chainKey)
        return true;
    return false;
}

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

bool compareComplexToComplexAlnByTmScore(const ComplexToComplexAln &first, const ComplexToComplexAln &second){
    if (first.tmScore > second.tmScore) {
        return true;
    }
    if (first.tmScore < second.tmScore) {
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
    if (first.tmScore > second.tmScore) {
        return true;
    }
    if (first.tmScore < second.tmScore) {
        return false;
    }
    return false;
}

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
                qComplexes.emplace_back(qComplex);
            }
        }
        qTempComplexes.clear();
        return qComplexes;
    }

    std::vector<OutputLine>  getOutputLines(const Complex qComplex){
        tmAligner = new TMaligner((unsigned int)(maxSeqLen*qComplex.chainKeys.size()), false);
        std::vector<OutputLine> outputLines;
        std::vector<unsigned int> qChainKeys = qComplex.chainKeys;
        std::vector<ChainToChainAln> tempChainAlnVec;
        std::vector<ComplexToComplexAln> complexAlnVec;
        unsigned int prevDbComplexId = qComplex.alnVec[0].dbChain.complexId;
        unsigned int prevQChainKey = qComplex.alnVec[0].qChain.chainKey;

        for (size_t chainToChainAlnIdx=0; chainToChainAlnIdx < qComplex.alnVec.size(); chainToChainAlnIdx++) {
            ChainToChainAln currAln = qComplex.alnVec[chainToChainAlnIdx];
            unsigned int dbCurrComplexId = currAln.dbChain.complexId;
            unsigned int qCurrChainKey = currAln.qChain.chainKey;

            if (prevDbComplexId != dbCurrComplexId) {
                complexAlnVec = updateComplexAlnVec(complexAlnVec, tempChainAlnVec);
                for (size_t lineIdx=0; lineIdx < complexAlnVec.size(); lineIdx++) {
                    ComplexToComplexAln aln = complexAlnVec[lineIdx];
                    if (aln.tmScore > 0){
                        std::vector<std::string> fields = std::vector<std::string>{
                            qComplexMap.at(aln.qComplexId),
                            dbComplexMap.at(aln.dbComplexId),
                            std::to_string(aln.tmScore),
                            getChainListString(aln.qChainKeys, qChainMap),
                            getChainListString(aln.dbChainKeys, dbChainMap),
                            std::to_string(lineIdx)
                        };
                        outputLines.emplace_back(OutputLine(aln.qComplexId, aln.tmScore, fields));
                    }
                }
                prevDbComplexId = dbCurrComplexId;
                prevQChainKey = qCurrChainKey;
                complexAlnVec.clear();
                tempChainAlnVec.clear();

            } else if (prevQChainKey != qCurrChainKey) {
                complexAlnVec = updateComplexAlnVec(complexAlnVec, tempChainAlnVec);
                prevQChainKey = qCurrChainKey;
                tempChainAlnVec.clear();
            }
            tempChainAlnVec.emplace_back(currAln);
        }
        complexAlnVec = updateComplexAlnVec(complexAlnVec, tempChainAlnVec);
        for (size_t lineIdx=0; lineIdx < complexAlnVec.size(); lineIdx++) {
            ComplexToComplexAln aln = complexAlnVec[lineIdx];
            if (aln.tmScore > 0){
                std::vector<std::string> fields = std::vector<std::string>{
                    qComplexMap.at(aln.qComplexId),
                    dbComplexMap.at(aln.dbComplexId),
                    std::to_string(aln.tmScore),
                    getChainListString(aln.qChainKeys, qChainMap),
                    getChainListString(aln.dbChainKeys, dbChainMap),
                    std::to_string(lineIdx)
                };
                outputLines.emplace_back(OutputLine(aln.qComplexId, aln.tmScore, fields));
            }
        }
        std::sort(outputLines.begin(), outputLines.end(), compareOutputLine);
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
        std::vector<Complex> qComplexes;
        unsigned int prevComplexId = qLookup.at(0);
        unsigned int complexId;
        std::vector<unsigned int> qTempChainKeys;
        for (size_t chainKey = 0; chainKey < qLookup.size(); chainKey++) {
            complexId = qLookup.at(chainKey);
            if (complexId != prevComplexId) {
                qComplexes.emplace_back(Complex(prevComplexId, qTempChainKeys, std::vector<ChainToChainAln>()));
                qTempChainKeys.clear();
            }
            prevComplexId = complexId;
            qTempChainKeys.emplace_back(chainKey);
        }
        qComplexes.emplace_back(Complex(prevComplexId, qTempChainKeys, std::vector<ChainToChainAln>()));
        qTempChainKeys.clear();
        std::sort(qComplexes.begin(), qComplexes.end(), compareComplex);

        for (size_t i = 0; i < alnDbr.getSize(); i++) {
            size_t queryKey = alnDbr.getDbKey(i);
            const unsigned queryComplexId = qLookup.at(queryKey);
            char *data = alnDbr.getData(i, thread_idx);
            if (*data == '\0') {
                continue;
            }
            Matcher::result_t qAlnResult = Matcher::parseAlignmentRecord(data);
            size_t qId = qCaDbr->sequenceReader->getId(queryKey);
            char *qCaData = qCaDbr->sequenceReader->getData(qId, thread_idx);
            size_t qCaLength = qCaDbr->sequenceReader->getEntryLen(qId);
            float* queryCaData = qCoords.read(qCaData, qAlnResult.qLen, qCaLength);
            Chain qChain = Chain(queryComplexId, queryKey);
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
                ChainToChainAln chainToChainAln(qChain, dbChain, queryCaData, targetCaData, alnResult);
                qComplexes[queryComplexId].alnVec.emplace_back(chainToChainAln);
                data = Util::skipLine(data);
            }
            std::sort(qComplexes[queryComplexId].alnVec.begin(), qComplexes[queryComplexId].alnVec.end(), compareChainToChainAlnByComplexIdAndChainKey);
        }
        return qComplexes;
    }

    std::vector<ComplexToComplexAln> updateComplexAlnVec(std::vector<ComplexToComplexAln> inputComplexAlnVec, std::vector<ChainToChainAln> chainAlnVec){
        // k cluster will be here
        std::vector<ComplexToComplexAln> tempComplexAlnVec;
        bool isInputEmpty = inputComplexAlnVec.empty();
        double maxTmScore;

        for (size_t alnVecIdx=0; alnVecIdx < chainAlnVec.size(); alnVecIdx++){
            // TODO Multithreading
            ChainToChainAln aln = chainAlnVec[alnVecIdx];
            std::vector<ChainToChainAln> alnVec;
            if (isInputEmpty) {
                ComplexToComplexAln currComplexAln(aln);
                double tmScore = getTmScore(currComplexAln);
                currComplexAln.setTmScore(tmScore);
                tempComplexAlnVec.emplace_back(currComplexAln);
            } else {
                for (size_t inputMatrixIdx = 0; inputMatrixIdx < inputComplexAlnVec.size(); inputMatrixIdx++) {
                    ComplexToComplexAln currComplexAln = inputComplexAlnVec[inputMatrixIdx];
                    currComplexAln.appendChainToChainAln(aln);
                    double tmScore = getTmScore(currComplexAln);
                    currComplexAln.setTmScore(tmScore);
                    tempComplexAlnVec.emplace_back(currComplexAln);
                }
            }
        }

        std::sort(tempComplexAlnVec.begin(), tempComplexAlnVec.end(), compareComplexToComplexAlnByTmScore);
        maxTmScore = tempComplexAlnVec[0].tmScore;
        unsigned int tempVecIdx=0;
        
        while (tempVecIdx < tempComplexAlnVec.size() && tempComplexAlnVec[tempVecIdx].tmScore>0){//maxTmScore-0.05){
//             TEMP
            ComplexToComplexAln tempCompAln = tempComplexAlnVec[tempVecIdx];
            std::cout<< "q:\t";
            for (size_t j=0; j<tempCompAln.qChainKeys.size(); j++) {
                std:: cout << tempCompAln.qChainKeys[j] << "\t";
            }
            std::cout<< "t:\t";
            for (size_t j=0; j<tempCompAln.dbChainKeys.size(); j++) {
                std:: cout << tempCompAln.dbChainKeys[j] << "\t";
            }
            std::cout << "tm:\t" << tempCompAln.tmScore << std::endl;

            tempVecIdx++;
        }
        std::vector<ComplexToComplexAln> outputComplexAlnVec(tempComplexAlnVec.begin(), tempComplexAlnVec.begin() + tempVecIdx);
        return outputComplexAlnVec;
    }

    double getTmScore(ComplexToComplexAln aln){
        bool chainOverlapAllowed = false;
        std::vector<unsigned int> foundDbKeys;
        unsigned int dbKeyIdx = 0;
        while (!chainOverlapAllowed && dbKeyIdx<aln.dbChainKeys.size()) {
            unsigned int currDbKey = aln.dbChainKeys[dbKeyIdx++];
            bool found = std::count(foundDbKeys.begin(), foundDbKeys.end(), currDbKey) > 0;
            if (found) {
                return 0;
            }
            foundDbKeys.emplace_back(currDbKey);
        }
        tmAligner->initQuery(&aln.qCaXVec[0], &aln.qCaYVec[0], &aln.qCaZVec[0], NULL, aln.qLength);
        TMaligner::TMscoreResult tmResult = tmAligner->computeTMscore2(&aln.dbCaXVec[0], &aln.dbCaYVec[0], &aln.dbCaZVec[0], aln.dbLength, 0, 0, Matcher::uncompressAlignment(aln.backtrace), aln.alnLength);
        return tmResult.tmscore;
    }

    static void getMaps(const std::string& file, std::map<unsigned int, unsigned int> &lookupMap, std::map<unsigned int, std::string> &complexNameMap, std::map<unsigned int, std::string> &chainNameMap){
        if (file.length() == 0) {
            return;
        }
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
            complexNameMap.emplace(Util::fast_atoi<unsigned int>(entry[2]), name.substr(0,4));
            chainNameMap.emplace(Util::fast_atoi<unsigned int>(entry[0]), name.substr(name.find_last_of('_')+1,name.find('\t')-name.find_last_of('_')-1));
        }
        lookup.close();
    }

    static std::string getChainListString(std::vector<unsigned int> chainKeys, std::map<unsigned int,std::string>map){
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
//        Coordinate16 qCoords;
//        Coordinate16 tCoords;
        ComplexScorer complexScorer(&q3DiDbr, t3DiDbr, qLookupFile, dbLookupFile, alnDbr, qCaDbr, tCaDbr, thread_idx);
        std::vector<Complex> qComplexes = complexScorer.getQComplexes();
//#pragma omp for schedule(dynamic, 1)
//        std::cout << "!!!\n";
        for (size_t qComplexIdx=0; qComplexIdx < qComplexes.size(); qComplexIdx++) {
            std::vector<OutputLine> resultLines = complexScorer.getOutputLines(qComplexes[qComplexIdx]);
            progress.updateProgress();
            for (size_t resultIdx=0; resultIdx < resultLines.size(); resultIdx++) {
                OutputLine result = resultLines[resultIdx];
                resultWriter.writeData(result.line.c_str(), result.line.size(), 0, thread_idx);
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