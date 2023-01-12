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
    int endPos;
    unsigned int length;

    void setLength(unsigned int len){
        length = len;
    }

    void setCaData(std::vector<float> caX, std::vector<float> caY, std::vector<float> caZ){
        caDataX = caX;
        caDataY = caY;
        caDataZ = caZ;
    }

    void setStartPosEndPos(int sPos, int ePos){
        startPos = sPos;
        endPos = ePos;
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
        qInputChain.setStartPosEndPos(0, 0);
        qInputChain.setCaData(qCaXVec, qCaYVec, qCaZVec);
        qInputChain.setLength(numMatches);
        dbInputChain.setStartPosEndPos(0, 0);
        dbInputChain.setCaData(dbCaXVec, dbCaYVec, dbCaZVec);
        dbInputChain.setLength(numMatches);
        qChain = qInputChain;
        dbChain = dbInputChain;
        backtrace = newBacktrace;
    }
    Chain qChain;
    Chain dbChain;
    std::string backtrace;
};

struct Complex {
    Complex() {}
    Complex(unsigned int complexId, std::vector<unsigned int> chainKeys, std::vector<ChainToChainAln> alnVec) : complexId(complexId), chainKeys(chainKeys), alnVec(alnVec) {}
    unsigned int complexId;
    std::vector<unsigned int> chainKeys;
    std::vector<ChainToChainAln> alnVec;
};

struct ComplexToComplexAln{
    ComplexToComplexAln() {}
    ComplexToComplexAln(unsigned int qComplexId, unsigned int dbComplexId, double tmScore): qComplexId(qComplexId), dbComplexId(dbComplexId), tmScore(tmScore){}
    unsigned int qComplexId;
    unsigned int dbComplexId;
    double tmScore;
    std::vector<unsigned int> qChainKeys;
    std::vector<unsigned int> dbChainKeys;

    void setTmScore(double tm){
        tmScore = tm;
    }
    void setChainKeys(std::vector<unsigned int> qChainKeyVec, std::vector<unsigned int> dbChainKeyVec){
        qChainKeys = qChainKeyVec;
        dbChainKeys = dbChainKeyVec;
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
//    OutputLine(unsigned int qComplexId, double tmScore, std::vector<std::string> fields, std::string separator, std::string endLine) : qComplexId(qComplexId), tmScore(tmScore) {
//        std::string separatedLine;
//        unsigned int i=0;
//        while (true) {
//            separatedLine.append(fields[i]);
//            i++;
//            if (i<fields.size()){
//                separatedLine.append(separator);
//            } else {
//                separatedLine.append(endLine);
//                break;
//            }
//        }
//        line = separatedLine;
//    }
    unsigned int qComplexId;
    double tmScore;
    std::string line;
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

//void getMaps(const std::string& file, std::map<unsigned int, unsigned int> &lookupMap, std::map<unsigned int, std::string> &complexNameMap, std::map<unsigned int, std::string> &chainNameMap){
//    if (file.length() == 0) {
//        return;
//    }
//    MemoryMapped lookup(file, MemoryMapped::WholeFile, MemoryMapped::SequentialScan);
//    char* data = (char *) lookup.getData();
//    const char* entry[255];
//    while (*data != '\0') {
//        const size_t columns = Util::getWordsOfLine(data, entry, 255);
//        if (columns < 3) {
//            Debug(Debug::WARNING) << "Not enough columns in lookup file " << file << "\n";
//            continue;
//        }
//        lookupMap.emplace(Util::fast_atoi<unsigned int>(entry[0]), Util::fast_atoi<unsigned int>(entry[2]));
//        data = Util::skipLine(data);
//        std::string name(entry[1], data - entry[1] - 1);
//        complexNameMap.emplace(Util::fast_atoi<unsigned int>(entry[2]), name.substr(0,4));
//        chainNameMap.emplace(Util::fast_atoi<unsigned int>(entry[0]), name.substr(name.find('_')+1,name.find('\t')-name.find('_')-1));
//    }
//    lookup.close();
//}

class ComplexScorer {
public:
    ComplexScorer(IndexReader *qDbr3Di, IndexReader *tDbr3Di, std::string qLookupFile, std::string tLookupFile){
        tmAligner = new TMaligner(std::max(qDbr3Di->sequenceReader->getMaxSeqLen()+1, tDbr3Di->sequenceReader->getMaxSeqLen()+1), false);
        getMaps(qLookupFile, qLookup, qComplexMap, qChainMap);
        getMaps(tLookupFile, tLookup, dbComplexMap, dbChainMap);
    }

    std::vector<Complex> getQComplexes(DBReader<unsigned int> &alnDbr, IndexReader *qCaDbr, IndexReader *tCaDbr, Coordinate16 qCoords, Coordinate16 tCoords, int thread_idx){
        std::vector<Complex> qComplexes = getQueryComplexVector();
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
            float *queryCaData = (float *) qCaData;
            if (qCaDbr->getDbtype() == LocalParameters::DBTYPE_CA_ALPHA_F16) {
                qCoords.read(qCaData, qAlnResult.qLen);
                queryCaData = qCoords.getBuffer();
            }
            Chain qChain = Chain(queryComplexId, queryKey);
            // getting chain to chain TM
            tmAligner->initQuery(queryCaData, &queryCaData[qAlnResult.qLen], &queryCaData[qAlnResult.qLen + qAlnResult.qLen], NULL, qAlnResult.qLen);

            while (*data != '\0') {
                char dbKeyBuffer[255 + 1];
                Util::parseKey(data, dbKeyBuffer);
                const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                const unsigned int dbComplexId = tLookup.at(dbKey);
                Matcher::result_t alnResult =  Matcher::parseAlignmentRecord(data);
                // get TM score
                size_t tCaId = tCaDbr->sequenceReader->getId(dbKey);
                char *tCaData = tCaDbr->sequenceReader->getData(tCaId, thread_idx);
                float *targetCaData = (float *) tCaData;
                if (tCaDbr->getDbtype() == LocalParameters::DBTYPE_CA_ALPHA_F16) {
                    tCoords.read(tCaData, alnResult.dbLen);
                    targetCaData = tCoords.getBuffer();
                }
                qChain.setStartPosEndPos(alnResult.qStartPos, alnResult.qEndPos);
                Chain dbChain = Chain(dbComplexId, dbKey);
                // getting chain to chain TM
                TMaligner::TMscoreResult tmResult = tmAligner->computeTMscore(targetCaData, &targetCaData[alnResult.dbLen], &targetCaData[alnResult.dbLen + alnResult.dbLen], alnResult.dbLen, alnResult.qStartPos, alnResult.dbStartPos, Matcher::uncompressAlignment(alnResult.backtrace));
                double tmScore = tmResult.tmscore;

                if (tmScore > 0) {
//                    ChainToChainAln chainToChainAln = getChainToChainAln(qChain, dbChain, queryCaData, targetCaData, alnResult);
                    ChainToChainAln chainToChainAln(qChain, dbChain, queryCaData, targetCaData, alnResult);
                    qComplexes[queryComplexId].alnVec.emplace_back(chainToChainAln);
                }
                data = Util::skipLine(data);
            }
            std::sort(qComplexes[queryComplexId].alnVec.begin(), qComplexes[queryComplexId].alnVec.end(), compareChainToChainAlnByComplexIdAndChainKey);
        }
        qComplexes = filterQComplexes(qComplexes);
        return qComplexes;
    }

    std::vector<OutputLine>  getOutputLines(const Complex qComplex){
        std::vector<OutputLine> outputLines;
        std::vector<unsigned int> qChainKeys = qComplex.chainKeys;
        std::vector<ChainToChainAln> tempChainAlnVec;
        std::vector<std::vector<ChainToChainAln>> chainAlnMatrix;
        unsigned int prevDbComplexId = qComplex.alnVec[0].dbChain.complexId;
        unsigned int prevQChainKey = qComplex.alnVec[0].qChain.chainKey;
        for (size_t chainToChainAlnIdx=0; chainToChainAlnIdx < qComplex.alnVec.size(); chainToChainAlnIdx++) {
            ChainToChainAln currAln = qComplex.alnVec[chainToChainAlnIdx];
            unsigned int dbCurrComplexId = currAln.dbChain.complexId;
            unsigned int qCurrChainKey = currAln.qChain.chainKey;
            if (prevDbComplexId != dbCurrComplexId) {
                chainAlnMatrix = updateChainToChainMatrix(chainAlnMatrix, tempChainAlnVec);
                double bestComplexTM = 0;
                ComplexToComplexAln bestCompAln;
                for (size_t matrixIdx=0; matrixIdx < chainAlnMatrix.size(); matrixIdx++) {
                    std::vector<ChainToChainAln> currChainToChainAlnVec = chainAlnMatrix[matrixIdx];
//                    old compatibility check
//                    if (currChainToChainAlnVec.size()!=qChainKeys.size()) {
//                        continue;
//                    }
                    ComplexToComplexAln currCompToCompAln = getComplexToComplexAlnWithComplexTmScore(
                            currChainToChainAlnVec);
                    if (currCompToCompAln.tmScore > bestComplexTM) {
                        bestComplexTM = currCompToCompAln.tmScore;
                        bestCompAln = currCompToCompAln;
                    }
                }
                if (bestComplexTM > 0){
                    std::vector<std::string> fields = std::vector<std::string>{qComplexMap.at(bestCompAln.qComplexId), dbComplexMap.at(bestCompAln.dbComplexId), std::to_string(bestCompAln.tmScore), getChainListString(bestCompAln.qChainKeys, qChainMap), getChainListString(bestCompAln.dbChainKeys, dbChainMap)};
                    outputLines.emplace_back(OutputLine(bestCompAln.qComplexId, bestCompAln.tmScore, fields));
                }
                prevDbComplexId = dbCurrComplexId;
                prevQChainKey = qCurrChainKey;
                chainAlnMatrix.clear();
                tempChainAlnVec.clear();

            } else if (prevQChainKey != qCurrChainKey) {
                chainAlnMatrix = updateChainToChainMatrix(chainAlnMatrix, tempChainAlnVec);
                prevQChainKey = qCurrChainKey;
                tempChainAlnVec.clear();
            }
            tempChainAlnVec.emplace_back(currAln);
        }
        chainAlnMatrix = updateChainToChainMatrix(chainAlnMatrix, tempChainAlnVec);
        double bestComplexTM = 0;
        ComplexToComplexAln bestCompAln;
        for (size_t matrixIdx=0; matrixIdx < chainAlnMatrix.size(); matrixIdx++) {
            std::vector<ChainToChainAln> currChainAlnVec = chainAlnMatrix[matrixIdx];
//            old compatibility check
//            if (currChainAlnVec.size() != qChainKeys.size()) {
//                continue;
//            }
            ComplexToComplexAln currCompAln = getComplexToComplexAlnWithComplexTmScore(currChainAlnVec);
            if (currCompAln.tmScore > bestComplexTM) {
                bestComplexTM = currCompAln.tmScore;
                bestCompAln = currCompAln;
            }
        }
        if (bestComplexTM > 0){
            std::vector<std::string> fields = std::vector<std::string>{qComplexMap.at(bestCompAln.qComplexId), dbComplexMap.at(bestCompAln.dbComplexId), std::to_string(bestCompAln.tmScore), getChainListString(bestCompAln.qChainKeys, qChainMap), getChainListString(bestCompAln.dbChainKeys, dbChainMap)};
            outputLines.emplace_back(OutputLine(bestCompAln.qComplexId, bestCompAln.tmScore, fields));
        }
        std::sort(outputLines.begin(), outputLines.end(), compareOutputLine);
        return outputLines;
    }
private:
    TMaligner * tmAligner;
    std::map<unsigned int, unsigned int> qLookup;
    std::map<unsigned int, unsigned int> tLookup;
    std::map<unsigned int, std::string> qComplexMap;
    std::map<unsigned int, std::string> dbComplexMap;
    std::map<unsigned int, std::string> qChainMap;
    std::map<unsigned int, std::string> dbChainMap;

    std::vector<Complex> getQueryComplexVector() {
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
        return qComplexes;
    }

//    ChainToChainAln getChainToChainAln(Chain qChain, Chain dbChain, float * queryCaData, float * targetCaData, Matcher::result_t alnResult){
//        int qStartPos = alnResult.qStartPos;
//        int dbStartPos = alnResult.dbStartPos;
//        unsigned int qLength = alnResult.qLen;
//        unsigned int dbLength = alnResult.dbLen;
//        std::string  backtrace = alnResult.backtrace;
//        std::vector<float> qCaXVec;
//        std::vector<float> qCaYVec;
//        std::vector<float> qCaZVec;
//        std::vector<float> dbCaXVec;
//        std::vector<float> dbCaYVec;
//        std::vector<float> dbCaZVec;
//        int numMatches = 0;
//        std::string newBacktrace;
//        int qPos = qStartPos;
//        int dbPos = dbStartPos;
//        int qXPos =0;
//        int qYPos = qLength;
//        int qZPos = qLength*2;
//        int dbXPos = 0;
//        int dbYPos = dbLength;
//        int dbZPos = dbLength*2;
//        for (size_t j=0; j<backtrace.size(); j++){
//            char cigar = backtrace[j];
//            switch (cigar) {
//                case 'M':
//                    numMatches++;
//                    newBacktrace += "M";
//                    qCaXVec.emplace_back(queryCaData[qXPos + qPos]);
//                    qCaYVec.emplace_back(queryCaData[qYPos + qPos]);
//                    qCaZVec.emplace_back(queryCaData[qZPos + qPos]);
//                    dbCaXVec.emplace_back(targetCaData[dbXPos + dbPos]);
//                    dbCaYVec.emplace_back(targetCaData[dbYPos + dbPos]);
//                    dbCaZVec.emplace_back(targetCaData[dbZPos + dbPos]);
//                    qPos++;
//                    dbPos++;
//                    break;
//                case 'I':
//                    qPos++;
//                    break;
//                case 'D':
//                    dbPos++;
//                    break;
//            }
//        }
//        qChain.setStartPosEndPos(0,0);
//        qChain.setCaData(qCaXVec, qCaYVec, qCaZVec);
//        qChain.setLength(numMatches);
//        dbChain.setStartPosEndPos(0,0);
//        dbChain.setCaData(dbCaXVec, dbCaYVec, dbCaZVec);
//        dbChain.setLength(numMatches);
//        return ChainToChainAln(qChain, dbChain, newBacktrace);
//    }

    std::vector<std::vector<ChainToChainAln>> updateChainToChainMatrix(std::vector<std::vector<ChainToChainAln>> inputMatrix, std::vector<ChainToChainAln> currAlnVec){
        std::vector<std::vector<ChainToChainAln>> outputMatrix;
        for (size_t alnVecIdx=0; alnVecIdx < currAlnVec.size(); alnVecIdx++){
            ChainToChainAln aln = currAlnVec[alnVecIdx];
            std::vector<ChainToChainAln> alnVec;
            if (inputMatrix.empty()) {
                alnVec.emplace_back(aln);
                outputMatrix.emplace_back(alnVec);
            } else {
                for (size_t inputMatrixIdx = 0; inputMatrixIdx < inputMatrix.size(); inputMatrixIdx++) {
                    alnVec.clear();
                    alnVec.insert(alnVec.end(), inputMatrix[inputMatrixIdx].begin(), inputMatrix[inputMatrixIdx].end());
                    alnVec.emplace_back(aln);
                    outputMatrix.emplace_back(alnVec);
                }
            }
        }
        return outputMatrix;
    }

    ComplexToComplexAln getComplexToComplexAlnWithComplexTmScore(std::vector<ChainToChainAln> chainToChainAlnVec){
        std::sort(chainToChainAlnVec.begin(), chainToChainAlnVec.end(), compareChainToChainAlnByQChainDbChain);
        std::vector<unsigned int> qChainKeys;
        std::vector<unsigned int> dbChainKeys;
        bool chainOverlapAllowed = false;
        std::vector<float> qCaXVec;
        std::vector<float> qCaYVec;
        std::vector<float> qCaZVec;
        std::vector<float> dbCaXVec;
        std::vector<float> dbCaYVec;
        std::vector<float> dbCaZVec;
        unsigned int numMatches = 0;
        std::string newBacktrace;
        unsigned int qComplexId = chainToChainAlnVec[0].qChain.complexId;
        unsigned int dbComplexId = chainToChainAlnVec[0].dbChain.complexId;
        ComplexToComplexAln outputAln(qComplexId, dbComplexId, 0);
        for (size_t i=0; i < chainToChainAlnVec.size(); i++) {
            Chain qChain = chainToChainAlnVec[i].qChain;
            Chain dbChain = chainToChainAlnVec[i].dbChain;
            if (std::count(dbChainKeys.begin(), dbChainKeys.end(), dbChain.chainKey) > 0 && !chainOverlapAllowed) {
                return outputAln;
            }
            std::string backtrace = chainToChainAlnVec[i].backtrace;
            newBacktrace += backtrace;
            numMatches += chainToChainAlnVec[i].qChain.length;
            qCaXVec.insert(qCaXVec.end(), qChain.caDataX.begin(), qChain.caDataX.end());
            qCaYVec.insert(qCaYVec.end(), qChain.caDataY.begin(), qChain.caDataY.end());
            qCaZVec.insert(qCaZVec.end(), qChain.caDataZ.begin(), qChain.caDataZ.end());
            dbCaXVec.insert(dbCaXVec.end(), dbChain.caDataX.begin(), dbChain.caDataX.end());
            dbCaYVec.insert(dbCaYVec.end(), dbChain.caDataY.begin(), dbChain.caDataY.end());
            dbCaZVec.insert(dbCaZVec.end(), dbChain.caDataZ.begin(), dbChain.caDataZ.end());
            qChainKeys.emplace_back(qChain.chainKey);
            dbChainKeys.emplace_back(dbChain.chainKey);
        }
        tmAligner->initQuery(&qCaXVec[0], &qCaYVec[0], &qCaZVec[0], NULL, numMatches);
        TMaligner::TMscoreResult tmResult = tmAligner->computeTMscore(&dbCaXVec[0], &dbCaYVec[0], &dbCaZVec[0], numMatches, 0, 0, Matcher::uncompressAlignment(newBacktrace));
        outputAln.setTmScore(tmResult.tmscore);
        outputAln.setChainKeys(qChainKeys, dbChainKeys);
        return outputAln;
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
            chainNameMap.emplace(Util::fast_atoi<unsigned int>(entry[0]), name.substr(name.find('_')+1,name.find('\t')-name.find('_')-1));
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
    static std::vector<Complex> filterQComplexes(std::vector<Complex>qComplexes){
        float compatibleCheckRatio = 1.0;
        for (size_t i=0; i<qComplexes.size(); i++) {
            std::vector<unsigned int> qChainKeys = qComplexes[i].chainKeys;
            std::vector<ChainToChainAln> newAlnVec;
            std::vector<ChainToChainAln> currAlnVec;
            unsigned int prevDbComplexId = qComplexes[i].alnVec[0].dbChain.complexId;
            unsigned int qChainCount = 0;
            unsigned int qChainIdx = 0;
            unsigned int currQChain = qChainKeys[qChainIdx];
            for (size_t j=0; j<qComplexes[i].alnVec.size(); j++) {
                ChainToChainAln aln = qComplexes[i].alnVec[j];
                if (aln.dbChain.complexId != prevDbComplexId) {
                    // isCompatible
                    if (qChainCount >= (float)(qChainKeys.size() * compatibleCheckRatio)) {
                        newAlnVec.insert(newAlnVec.end(), currAlnVec.begin(), currAlnVec.end());
                    }
                    prevDbComplexId = aln.dbChain.complexId;
                    qChainCount = 0;
                    qChainIdx = 0;
                    currQChain = qChainKeys[qChainIdx];
                    currAlnVec.clear();
                }
                currAlnVec.emplace_back(aln);
                if (aln.qChain.chainKey == currQChain) {
                    qChainCount++;
                    currQChain = qChainKeys[++qChainIdx];
                }
            }
            qComplexes[i].alnVec = newAlnVec;
        }
        return qComplexes;
    }
};

int scorecomplex(int argc, const char **argv, const Command& command) {
    std::cout << "-1" << std::endl;
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);
    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    IndexReader q3DiDbr(StructureUtil::getIndexWithSuffix(par.db1, "_ss"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
    IndexReader *t3DiDbr = NULL;
    IndexReader *qCaDbr = new IndexReader(par.db1, par.threads, IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY), touch ? IndexReader::PRELOAD_INDEX : 0,  DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA, "_ca" );
    IndexReader *tCaDbr = NULL;
    std::cout << "0" << std::endl;
    bool sameDB = false;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        t3DiDbr = &q3DiDbr;
        tCaDbr = qCaDbr;
    } else {
        t3DiDbr = new IndexReader(StructureUtil::getIndexWithSuffix(par.db2, "_ss"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
        tCaDbr = new IndexReader(par.db2, par.threads, IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY), touch ? IndexReader::PRELOAD_INDEX : 0, DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA, "_ca");
    }
    std::cout << "1" << std::endl;
    std::string qLookupFile = par.db1 + ".lookup";
    std::string dbLookupFile = par.db2 + ".lookup";
    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    DBWriter resultWriter(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, Parameters::DBTYPE_ALIGNMENT_RES);
    resultWriter.open();
    Debug::Progress progress(alnDbr.getSize());
    std::cout << "2" << std::endl;
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::vector<Matcher::result_t> alignmentResult;
        std::string backtrace;
        std::string resultBuffer;
        Coordinate16 qCoords;
        Coordinate16 tCoords;
        // workflow 1,2
        std::cout << "3" << std::endl;
        ComplexScorer complexScorer(&q3DiDbr, t3DiDbr, qLookupFile, dbLookupFile);
        std::vector<Complex> qComplexes = complexScorer.getQComplexes(alnDbr, qCaDbr, tCaDbr, qCoords, tCoords, thread_idx);
        std::cout << "4" << std::endl;
#pragma omp for schedule(dynamic, 1)
        // workflow 3,4
        for (size_t qComplexIdx=0; qComplexIdx < qComplexes.size(); qComplexIdx++) {
            std::vector<OutputLine> resultLines = complexScorer.getOutputLines(qComplexes[qComplexIdx]);
            progress.updateProgress();
            for (size_t resultIdx=0; resultIdx < resultLines.size(); resultIdx++) {
                OutputLine result = resultLines[resultIdx];
                resultWriter.writeData(result.line.c_str(), result.line.size(), 0, thread_idx);
                // TEMP
                std::cout << result.line;
            }
        }
    }
    std::cout << "5" << std::endl;
    alnDbr.close();
    if (!sameDB) {
        delete t3DiDbr;
    }
    std::cout << "6" << std::endl;
    resultWriter.close(true);
    std::cout << "7" << std::endl;
    return EXIT_SUCCESS;
}