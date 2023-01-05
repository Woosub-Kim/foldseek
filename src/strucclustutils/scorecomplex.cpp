#include "DBReader.h"
#include "IndexReader.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "LocalParameters.h"
#include "Matcher.h"
#include "Alignment.h"
#include "structureto3diseqdist.h"
#include "StructureSmithWaterman.h"
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
    Chain(unsigned int complexid, unsigned int chainKey, float *caData, unsigned int length) : complexid(complexid), chainKey(chainKey), caData(caData), length(length) {}
    Chain(unsigned int complexid, unsigned int chainKey, float *caData, int startPos, int endPos, unsigned int length) : complexid(complexid), chainKey(chainKey), caData(caData), startPos(startPos), endPos(endPos), length(length) {}
    unsigned int complexid;
    unsigned int chainKey;
    float *caData;
    int startPos;
    int endPos;
    unsigned int length;

    void setStartPos(int pos){
        startPos = pos;
    }

    void setEndPos(int pos){
        endPos = pos;
    }

    void setStartPosEndPos(int sPos, int ePos){
        startPos = sPos;
        endPos = ePos;
    }

    bool operator<(const Chain &o) const {
        if (complexid < o.complexid ) {
            return true;
        }
        if (complexid > o.complexid) {
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
        if (complexid > o.complexid ) {
            return true;
        }
        if (complexid < o.complexid) {
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
        return complexid==o.complexid && chainKey==o.chainKey;
    }

};

struct ChainToChainAln {
    ChainToChainAln() {}
    ChainToChainAln(Chain chain1, Chain chain2, double tmScore, std::string backtrace): qChain(chain1), dbChain(chain2), tmScore(tmScore), backtrace(backtrace) {}
    Chain qChain;
    Chain dbChain;
    double tmScore;
    std::string backtrace;
};

struct Complex {
    Complex() {}
    Complex(unsigned int complexid, std::vector<unsigned int> chainKeys, std::vector<ChainToChainAln> alns) : complexid(complexid), chainKeys(chainKeys), alns(alns) {}
    unsigned int complexid;
    std::vector<unsigned int> chainKeys;
    std::vector<ChainToChainAln> alns;
};

struct ComplexToComplexAln{
    ComplexToComplexAln() {}
    ComplexToComplexAln(unsigned int qComplexid, unsigned int dbComplexid, double tmScore):qComplexid(qComplexid), dbComplexid(dbComplexid), tmScore(tmScore){}
    unsigned int qComplexid;
    unsigned int dbComplexid;
    double tmScore;
};

//bool compareChain(const Chain &first, const Chain &second){
//    if (first.complexid < second.complexid)
//        return true;
//    if (first.complexid > second.complexid)
//        return false;
//    if (first.chainKey < second.chainKey)
//        return true;
//    if (first.chainKey > second.chainKey)
//        return false;
//    return false;
//}

bool compareChainToChainAlnByQueryChain(const ChainToChainAln &first, const ChainToChainAln &second){
    if (first.qChain.complexid < second.qChain.complexid)
        return true;
    if (first.qChain.complexid > second.qChain.complexid)
        return false;
    if (first.qChain.chainKey < second.qChain.chainKey)
        return true;
    if (first.qChain.chainKey > second.qChain.chainKey)
        return false;
    if (first.dbChain.complexid < second.dbChain.complexid)
        return true;
    if (first.dbChain.complexid > second.dbChain.complexid)
        return false;
    if (first.dbChain.chainKey < second.dbChain.chainKey)
        return true;
    return false;
}

bool compareComplex(const Complex &first, const Complex &second){
    if (first.complexid < second.complexid) {
        return true;
    }
    if (first.complexid > second.complexid) {
        return false;
    }
    return false;
}

bool compareChainToChainAlnByDbChain(const ChainToChainAln &first, const ChainToChainAln &second){
    if (first.dbChain.complexid < second.dbChain.complexid)
        return true;
    if (first.dbChain.complexid > second.dbChain.complexid)
        return false;
//    if (first.dbChain.chainKey < second.dbChain.chainKey)
//        return true;
//    if (first.dbChain.chainKey > second.dbChain.chainKey)
//        return false;
    if (first.qChain.complexid < second.qChain.complexid)
        return true;
    if (first.qChain.complexid > second.qChain.complexid)
        return false;
    if (first.qChain.chainKey < second.qChain.chainKey)
        return true;
    if (first.qChain.chainKey > second.qChain.chainKey)
        return false;

    return false;
}

std::map<unsigned int, unsigned int> getLookupMap(const std::string& file) {
    std::map<unsigned int, unsigned int> mapping;
    if (file.length() == 0) {
        return mapping;
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
        mapping.emplace(Util::fast_atoi<unsigned int>(entry[0]), Util::fast_atoi<unsigned int>(entry[2]));
        data = Util::skipLine(data);
    }
    lookup.close();
    return mapping;
}
class ComplexScorer {
public:
    ComplexScorer(IndexReader *qDbr3Di, IndexReader *tDbr3Di, std::map<unsigned int, unsigned int> qLookup, std::map<unsigned int, unsigned int> tLookup) : qDbr3Di(qDbr3Di), tDbr3Di(tDbr3Di), qLookup(qLookup), tLookup(tLookup) {
        tmaligner = new TMaligner(
                std::max(qDbr3Di->sequenceReader->getMaxSeqLen()+1, tDbr3Di->sequenceReader->getMaxSeqLen()+1),
                false);
    }

    std::vector<Complex> parseAlnDb(
            DBReader<unsigned int> alnDbr,
            IndexReader *qcadbr,
            IndexReader *tcadbr,
            Coordinate16 qcoords,
            Coordinate16 tcoords,
            int thread_idx
    ){
        std::vector<Complex> qComplexes = getComplexVector(qLookup);
        std::vector<ChainToChainAln> resultVector;
        for (size_t i = 0; i < alnDbr.getSize(); i++) {
            size_t queryKey = alnDbr.getDbKey(i);
            const unsigned queryComplexId = qLookup.at(queryKey);
            char *data = alnDbr.getData(i, thread_idx);
            if (*data == '\0') {
                continue;
            }
            Matcher::result_t qAlnResult = Matcher::parseAlignmentRecord(data);
            size_t qId = qcadbr->sequenceReader->getId(queryKey);
            char *qcadata = qcadbr->sequenceReader->getData(qId, thread_idx);
            float *queryCaData = (float *) qcadata;
            if (qcadbr->getDbtype() == LocalParameters::DBTYPE_CA_ALPHA_F16) {
                qcoords.read(qcadata, qAlnResult.qLen);
                queryCaData = qcoords.getBuffer();
            }
            Chain qChain = Chain(queryComplexId, queryKey, queryCaData, qAlnResult.qLen);

            tmaligner->initQuery(queryCaData, &queryCaData[qAlnResult.qLen], &queryCaData[qAlnResult.qLen+qAlnResult.qLen],
                                 NULL, qAlnResult.qLen);

            while (*data != '\0') {
                char dbKeyBuffer[255 + 1];
                Util::parseKey(data, dbKeyBuffer);
                const unsigned int dbKey = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                const unsigned int dbComplexId = tLookup.at(dbKey);

                Matcher::result_t alnResult =  Matcher::parseAlignmentRecord(data);
                // get TM score
                size_t tCaId = tcadbr->sequenceReader->getId(dbKey);
                char *tcadata = tcadbr->sequenceReader->getData(tCaId, thread_idx);
                float *targetCaData = (float *) tcadata;
                if (tcadbr->getDbtype() == LocalParameters::DBTYPE_CA_ALPHA_F16) {
                    tcoords.read(tcadata, alnResult.dbLen);
                    targetCaData = tcoords.getBuffer();
                }
                qChain.setStartPosEndPos(alnResult.qStartPos, alnResult.qEndPos);
                Chain dbChain = Chain(dbComplexId, dbKey, targetCaData, alnResult.dbStartPos, alnResult.dbEndPos, alnResult.dbLen);
                TMaligner::TMscoreResult tmres = tmaligner->computeTMscore(targetCaData, &targetCaData[alnResult.dbLen], &targetCaData[alnResult.dbLen+alnResult.dbLen], alnResult.dbLen,alnResult.qStartPos, alnResult.dbStartPos, Matcher::uncompressAlignment(alnResult.backtrace));
                Coordinates coord = Coordinates(alnResult.dbLen);
                double tmScore = tmres.tmscore;
                ChainToChainAln chainToChainAln = ChainToChainAln(qChain, dbChain, tmScore, alnResult.backtrace);
                if (tmScore > 0) {
                    qComplexes[queryComplexId].alns.emplace_back(chainToChainAln);
                }
                data = Util::skipLine(data);
            }
            std::sort(
                    qComplexes[queryComplexId].alns.begin(),
                    qComplexes[queryComplexId].alns.end(),
                    compareChainToChainAlnByDbChain
                    );
        }
        return qComplexes;
    }

    std::vector<ComplexToComplexAln>  returnOutputResult(std::vector<Complex> resultVector){
        std::vector<ComplexToComplexAln> outputVec;

        for (size_t qComplexIdx=0; qComplexIdx < resultVector.size(); qComplexIdx++) {
            unsigned int qComplexid = resultVector[qComplexIdx].complexid;
            std::vector<unsigned int> qChainKeys = resultVector[qComplexIdx].chainKeys;
            std::vector<ChainToChainAln> chainToChainAlns = resultVector[qComplexIdx].alns;
            std::vector<ChainToChainAln> bestAlnsForCurrDbComplex;
            unsigned int prevDbComplexid = chainToChainAlns[0].dbChain.complexid;

            for (size_t chainToChainAlnIdx=0; chainToChainAlnIdx < chainToChainAlns.size(); chainToChainAlnIdx++) {
                ChainToChainAln currAln = chainToChainAlns[chainToChainAlnIdx];
                unsigned int currDbComplexid = currAln.dbChain.complexid;
                // dbComplex is changed
                if (currDbComplexid != prevDbComplexid) {
                    // isCompatible
                    if (bestAlnsForCurrDbComplex.size() == qChainKeys.size()) {
                        outputVec = calcComplexTmScore(bestAlnsForCurrDbComplex, outputVec);
                    }
                    bestAlnsForCurrDbComplex.clear();
                    prevDbComplexid = currDbComplexid;
                }
                bestAlnsForCurrDbComplex = updateBestAlnsForCurrDbComplex(bestAlnsForCurrDbComplex, currAln);
            }
        }
        return outputVec;
    }
private:
    IndexReader * qDbr3Di;
    IndexReader * tDbr3Di;
    TMaligner * tmaligner;
    std::map<unsigned int, unsigned int> qLookup;
    std::map<unsigned int, unsigned int> tLookup;

    std::vector<Complex> getComplexVector(std::map<unsigned int, unsigned int> lookup) {
        std::vector<Complex> qComplexes;
        unsigned int prevComplexid = lookup.at(0);
        unsigned int complexid;
        std::vector<unsigned int> qTempChainKeys;
        for (size_t chainKey = 0; chainKey < lookup.size(); chainKey++) {
            complexid = lookup.at(chainKey);
            if (complexid != prevComplexid) {
                qComplexes.emplace_back(Complex(prevComplexid, qTempChainKeys, std::vector<ChainToChainAln>()));
                qTempChainKeys.clear();
            }
            prevComplexid = complexid;
            qTempChainKeys.emplace_back(chainKey);
        }
        qComplexes.emplace_back(Complex(prevComplexid, qTempChainKeys, std::vector<ChainToChainAln>()));
        qTempChainKeys.clear();
        std::sort(qComplexes.begin(), qComplexes.end(), compareComplex);
        return qComplexes;
    }

//    std::map<unsigned int, bool> resetFoundChainKeys(std::vector<unsigned int> qChainKeys){
//        std::map<unsigned int, bool> foundChainKeys;
//        for (size_t i=0; i<qChainKeys.size(); i++){
//            foundChainKeys.insert({qChainKeys[i], false});
//        }
//        return foundChainKeys;
//    }
    std::vector<ChainToChainAln> updateBestAlnsForCurrDbComplex(std::vector<ChainToChainAln> bestAlnsForCurrDbComplex, ChainToChainAln currAln){
        if (bestAlnsForCurrDbComplex.size()==0) {
            bestAlnsForCurrDbComplex.emplace_back(currAln);
            return bestAlnsForCurrDbComplex;
        }
        ChainToChainAln prevAln = bestAlnsForCurrDbComplex[bestAlnsForCurrDbComplex.size()-1];
        if (prevAln.qChain.chainKey!=currAln.qChain.chainKey){
            bestAlnsForCurrDbComplex.emplace_back(currAln);
        } else {
            bestAlnsForCurrDbComplex[bestAlnsForCurrDbComplex.size()-1] = (prevAln.tmScore > currAln.tmScore) ? prevAln : currAln;
        }
        return bestAlnsForCurrDbComplex;
    }

    std::vector<ComplexToComplexAln> calcComplexTmScore(std::vector<ChainToChainAln> tempVec, std::vector<ComplexToComplexAln> outputVec){
        std::sort(tempVec.begin(), tempVec.end(), compareChainToChainAlnByQueryChain);
//        std::vector<Chain> qChainVector;
        std::vector<Chain> dbChainVector;
//        std::vector<double> tmScoreVector;
        bool chainOverlapAllowed = false;
        std::vector<float> qCaXVec;
        std::vector<float> qCaYVec;
        std::vector<float> qCaZVec;
        std::vector<float> dbCaXVec;
        std::vector<float> dbCaYVec;
        std::vector<float> dbCaZVec;
        int numMatches = 0;
        std::string newBacktrace = "";
        unsigned int qComplexid = tempVec[0].qChain.complexid;
        unsigned int dbComplexid = tempVec[0].dbChain.complexid;
        for (size_t i=0; i<tempVec.size(); i++) {
            Chain qChain = tempVec[i].qChain;
            Chain dbChain = tempVec[i].dbChain;
            std::string backtrace = tempVec[i].backtrace;
            if (std::count(dbChainVector.begin(), dbChainVector.end(), dbChain)>0 && !chainOverlapAllowed) {
                return outputVec;
            }
            int qPos = qChain.startPos-1;
            int dbPos = dbChain.startPos - 1;
//            qChainVector.emplace_back(qChain);
            dbChainVector.emplace_back(dbChain);
//            tmScoreVector.emplace_back(tempVec[i].tmScore);

            for (size_t j=0; j<backtrace.size(); j++){
                char cigar = backtrace[j];
                switch (cigar) {
                    case 'M':
                        qPos++;
                        dbPos++;
                        numMatches++;
                        newBacktrace += "M";
                        qCaXVec.emplace_back(qChain.caData[0 + qPos]);
                        qCaYVec.emplace_back(qChain.caData[qChain.length + qPos]);
                        qCaZVec.emplace_back(qChain.caData[qChain.length * 2 + qPos]);
                        dbCaXVec.emplace_back(dbChain.caData[0 + dbPos]);
                        dbCaYVec.emplace_back(dbChain.caData[dbChain.length + dbPos]);
                        dbCaZVec.emplace_back(dbChain.caData[dbChain.length * 2 + dbPos]);
                        break;
                    case 'I':
                        qPos++;
                        break;
                    case 'D':
                        dbPos++;
                        break;
                }
            }
        }
        for (size_t i=0; i<tempVec.size(); i++) {
            Chain qChain = tempVec[i].qChain;
            Chain dbChain = tempVec[i].dbChain;
            std::cout << qChain.complexid << "\t" << qChain.chainKey << "\t" << dbChain.complexid << "\t" << dbChain.chainKey << "\t" << tempVec[i].tmScore << std::endl;
        }
        tmaligner->initQuery(&qCaXVec[0], &qCaYVec[0], &qCaZVec[0], NULL, numMatches);
        TMaligner::TMscoreResult tmres = tmaligner->computeTMscore(&dbCaXVec[0], &dbCaYVec[0], &dbCaZVec[0], numMatches, 0, 0, Matcher::uncompressAlignment(newBacktrace));
        double complexTmScore = tmres.tmscore;
        outputVec.emplace_back(qComplexid, dbComplexid, complexTmScore);

        return outputVec;
    }
};

int scorecomplex(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
//    int dbaccessMode = DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA;

    IndexReader qdbr3Di(StructureUtil::getIndexWithSuffix(par.db1, "_ss"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
    IndexReader *t3DiDbr = NULL;
    IndexReader *qcadbr = new IndexReader(
            par.db1,
            par.threads,
            IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY),
            touch ? IndexReader::PRELOAD_INDEX : 0,
            DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
            "_ca"
            );
    IndexReader *tcadbr = NULL;
    bool sameDB = false;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        t3DiDbr = &qdbr3Di;
        tcadbr = qcadbr;
    } else {
        t3DiDbr = new IndexReader(StructureUtil::getIndexWithSuffix(par.db2, "_ss"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
        tcadbr = new IndexReader(
                par.db2,
                par.threads,
                IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY),
                touch ? IndexReader::PRELOAD_INDEX : 0,
                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
                "_ca"
                );
    }
    std::map<unsigned int, unsigned int> qLookup;
    std::map<unsigned int, unsigned int> tLookup;
    std::string file1 = par.db1 + ".lookup";
    std::string file2 = par.db2 + ".lookup";
    qLookup = getLookupMap(file1);
    tLookup = getLookupMap(file2);

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed,  Parameters::DBTYPE_ALIGNMENT_RES);
    dbw.open();
    Debug::Progress progress(alnDbr.getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        std::vector<Matcher::result_t> alignmentResult;
        std::string backtrace;
        char buffer[1024+32768];
        std::string resultBuffer;
        Coordinate16 qcoords;
        Coordinate16 tcoords;
        ComplexScorer complexScorer(&qdbr3Di, t3DiDbr, qLookup, tLookup);
        std::vector<Complex> qComplexes = complexScorer.parseAlnDb(alnDbr, qcadbr, tcadbr, qcoords, tcoords, thread_idx);
        std::vector<ComplexToComplexAln> tempVec3 = complexScorer.returnOutputResult(qComplexes);
#pragma omp for schedule(dynamic, 1)

        for (size_t i=0; i<tempVec3.size(); i++) {
            ComplexToComplexAln complexToComplexAln = tempVec3[i];
            std::cout << complexToComplexAln.qComplexid << "\t" << complexToComplexAln.dbComplexid << "\t" << complexToComplexAln.tmScore << std::endl;
        }
    }


//    free(tinySubMatAA);
//    free(tinySubMat3Di);
    dbw.close();
    alnDbr.close();
    if (sameDB == false) {
        delete t3DiDbr;
//        delete tAADbr;
    }
    return EXIT_SUCCESS;
}