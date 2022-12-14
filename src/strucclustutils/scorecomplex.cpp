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

class ComplexScorer  {

};


struct Chain {
    Chain() {}
    Chain(unsigned int complexid, unsigned int chainKey, float *caData, unsigned int length) : complexid(complexid), chainKey(chainKey), caData(caData), length(length) {}
    unsigned int complexid;
    unsigned int chainKey;
    float *caData;
    unsigned int length;
};

struct ChainToChainAln {
    ChainToChainAln() {}
    ChainToChainAln(Chain chain1, Chain chain2, double tmScore, std::string backtrace): qChain(chain1), dbChain(chain2), tmScore(tmScore), backtrace(backtrace) {}
    Chain qChain;
    Chain dbChain;
    double tmScore;
    std::string backtrace;
};

bool compareChainsWithTmByQueryChain(const ChainToChainAln &first, const ChainToChainAln &second){
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

bool compareChainsWithTmByDbChain(const ChainToChainAln &first, const ChainToChainAln &second){
    if (first.dbChain.complexid < second.dbChain.complexid)
        return true;
    if (first.dbChain.complexid > second.dbChain.complexid)
        return false;
    if (first.dbChain.chainKey < second.dbChain.chainKey)
        return true;
    if (first.dbChain.chainKey > second.dbChain.chainKey)
        return false;
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

std::vector<ChainToChainAln> update(std::vector<ChainToChainAln> tempVector1, std::vector<ChainToChainAln> tempVector2){
    bool NEED_OVERLAP = false;
    if (tempVector2.size()==0){
        return tempVector1;
    }
    std::vector<ChainToChainAln> foo;
    for (size_t i=0; i<tempVector1.size(); i++) {
        ChainToChainAln currAln = tempVector1[i];
        std::vector<ChainToChainAln> egg;
        Chain currQueryChain = currAln.qChain;
        Chain currDbChain = currAln.dbChain;
        unsigned int currDbKey = currDbChain.chainKey;
        bool needKeep = false;
        for (size_t j=0; j<tempVector2.size(); j++){
            ChainToChainAln prevAln = tempVector2[j];
            Chain prevQueryChain = prevAln.qChain;
            Chain prevDbChain = prevAln.dbChain;
            unsigned int prevDbKey = prevDbChain.chainKey;

            if (currQueryChain.complexid == prevQueryChain.complexid and currDbChain.complexid == prevDbChain.complexid){
                if (prevDbKey == currDbKey and not NEED_OVERLAP){
                    needKeep = false;
                    break;
                }
                egg.emplace_back(prevAln);
                needKeep = true;
            }
        }
        if (needKeep){
            std::sort(egg.begin(), egg.end(), compareChainsWithTmByQueryChain);
            for (size_t k=0; k<egg.size(); k++){
                foo.emplace_back(egg[k]);
            }
            egg.clear();
            foo.emplace_back(currAln);
        }
    }
    return foo;
}

void printOut(std::vector<ChainToChainAln> tempVector1, std::vector<ChainToChainAln> tempVector2){
    bool NEED_OVERLAP = false;
    if (tempVector2.size()==0){
        return;
    }
    std::vector<ChainToChainAln> foo;
    for (size_t i=0; i<tempVector1.size(); i++) {
        ChainToChainAln currAln = tempVector1[i];
        std::vector<ChainToChainAln> egg;
        Chain currQueryChain = currAln.qChain;
        Chain currDbChain = currAln.dbChain;
        unsigned int currDbKey = currDbChain.chainKey;
        bool needKeep = false;

        for (size_t j=0; j<tempVector2.size(); j++){
            ChainToChainAln prevAln = tempVector2[j];
            Chain prevQueryChain = prevAln.qChain;
            Chain prevDbChain = prevAln.dbChain;
            unsigned int prevDbKey = prevDbChain.chainKey;

            if (currQueryChain.complexid == prevQueryChain.complexid and currDbChain.complexid == prevDbChain.complexid){
                if (prevDbKey == currDbKey and not NEED_OVERLAP){
                    needKeep = false;
                    break;
                }
                egg.emplace_back(prevAln);
                needKeep = true;
            }
        }
        if (needKeep){
            egg.emplace_back(currAln);
            std::sort(egg.begin(), egg.end(), compareChainsWithTmByQueryChain);
            for (size_t i=0; i<egg.size(); i++) {
                ChainToChainAln aln = egg[i];
                Chain qChain = aln.qChain;
                Chain dbChain = aln.dbChain;
                double tmScore = aln.tmScore;
                std::cout << qChain.complexid << "\t" << qChain.chainKey << "\t" << dbChain.complexid << "\t" << dbChain.chainKey << "\t" << tmScore << std::endl;
            }

            egg.clear();
            foo.emplace_back(currAln);
        }
    }
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


int scorecomplex(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_ALIGN);

    const bool touch = (par.preloadMode != Parameters::PRELOAD_MODE_MMAP);
    IndexReader qdbrAA(par.db1, par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
    IndexReader qdbr3Di(StructureUtil::getIndexWithSuffix(par.db1, "_ss"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);

    IndexReader *t3DiDbr = NULL;
    IndexReader *tAADbr = NULL;
    bool sameDB = false;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        t3DiDbr = &qdbr3Di;
        tAADbr = &qdbrAA;
    } else {
        tAADbr = new IndexReader(par.db2, par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
        t3DiDbr = new IndexReader(StructureUtil::getIndexWithSuffix(par.db2, "_ss"), par.threads, IndexReader::SEQUENCES, touch ? IndexReader::PRELOAD_INDEX : 0);
    }
    bool needLookup = true;
    std::map<unsigned int, unsigned int> qLookup;
    std::map<unsigned int, unsigned int> tLookup;
    if (needLookup) {
        std::string file1 = par.db1 + ".lookup";
        std::string file2 = par.db2 + ".lookup";
        qLookup = getLookupMap(file1);
        tLookup = getLookupMap(file2);
    }

    DBReader<unsigned int> alnDbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
    alnDbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed,  Parameters::DBTYPE_ALIGNMENT_RES);
    dbw.open();

    int dbaccessMode = DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA;
    IndexReader qDbr(par.db1, par.threads,  IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
    IndexReader qDbrHeader(par.db1, par.threads, IndexReader::SRC_HEADERS , (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);

    IndexReader *tDbr;
    IndexReader *tDbrHeader;
    if (sameDB) {
        tDbr = &qDbr;
        tDbrHeader= &qDbrHeader;
    } else {
        tDbr = new IndexReader(par.db2, par.threads, IndexReader::SRC_SEQUENCES, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0, dbaccessMode);
        tDbrHeader = new IndexReader(par.db2, par.threads, IndexReader::SRC_HEADERS, (touch) ? (IndexReader::PRELOAD_INDEX | IndexReader::PRELOAD_DATA) : 0);
    }

    bool needTMaligner = true;//(par.tmScoreThr > 0);
    IndexReader *qcadbr = NULL;
    IndexReader *tcadbr = NULL;
    if(needTMaligner){
        qcadbr = new IndexReader(
                par.db1,
                par.threads,
                IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY),
                touch ? IndexReader::PRELOAD_INDEX : 0,
                DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
                "_ca");
        if (sameDB) {
            tcadbr = qcadbr;
        } else {
            tcadbr = new IndexReader(
                    par.db2,
                    par.threads,
                    IndexReader::makeUserDatabaseType(LocalParameters::INDEX_DB_CA_KEY),
                    touch ? IndexReader::PRELOAD_INDEX : 0,
                    DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA,
                    "_ca"
            );
        }
    }

    SubstitutionMatrix subMat3Di(par.scoringMatrixFile.values.aminoacid().c_str(), 2.1, par.scoreBias);
    std::string blosum;
    for (size_t i = 0; i < par.substitutionMatrices.size(); i++) {
        if (par.substitutionMatrices[i].name == "blosum62.out") {
            std::string matrixData((const char *)par.substitutionMatrices[i].subMatData, par.substitutionMatrices[i].subMatDataLen);
            std::string matrixName = par.substitutionMatrices[i].name;
            char * serializedMatrix = BaseMatrix::serialize(matrixName, matrixData);
            blosum.assign(serializedMatrix);
            free(serializedMatrix);
            break;
        }
    }
    SubstitutionMatrix subMatAA(blosum.c_str(), 1.4, par.scoreBias);
    //temporary output file
    Debug::Progress progress(alnDbr.getSize());

    // sub. mat needed for query profile
    int8_t * tinySubMatAA = (int8_t*) mem_align(ALIGN_INT, subMatAA.alphabetSize * 32);
    int8_t * tinySubMat3Di = (int8_t*) mem_align(ALIGN_INT, subMat3Di.alphabetSize * 32);

    for (int i = 0; i < subMat3Di.alphabetSize; i++) {
        for (int j = 0; j < subMat3Di.alphabetSize; j++) {
            tinySubMat3Di[i * subMat3Di.alphabetSize + j] = subMat3Di.subMatrix[i][j]; // for farrar profile
        }
    }
    for (int i = 0; i < subMatAA.alphabetSize; i++) {
        for (int j = 0; j < subMatAA.alphabetSize; j++) {
            tinySubMatAA[i * subMatAA.alphabetSize + j] = subMatAA.subMatrix[i][j];
        }
    }

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        EvalueNeuralNet evaluer(tAADbr->sequenceReader->getAminoAcidDBSize(), &subMat3Di);
        std::vector<Matcher::result_t> alignmentResult;
        StructureSmithWaterman structureSmithWaterman(par.maxSeqLen, subMat3Di.alphabetSize, par.compBiasCorrection, par.compBiasCorrectionScale);
        StructureSmithWaterman reverseStructureSmithWaterman(par.maxSeqLen, subMat3Di.alphabetSize, par.compBiasCorrection, par.compBiasCorrectionScale);
        TMaligner *tmaligner = NULL;
        if(needTMaligner) {
            tmaligner = new TMaligner(
                    std::max(qdbr3Di.sequenceReader->getMaxSeqLen() + 1, t3DiDbr->sequenceReader->getMaxSeqLen() + 1), false);
        }
        Sequence qSeqAA(par.maxSeqLen, qdbrAA.getDbtype(), (const BaseMatrix *) &subMatAA, 0, false, par.compBiasCorrection);
        Sequence qSeq3Di(par.maxSeqLen, qdbr3Di.getDbtype(), (const BaseMatrix *) &subMat3Di, 0, false, par.compBiasCorrection);
        Sequence tSeqAA(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMatAA, 0, false, par.compBiasCorrection);
        Sequence tSeq3Di(par.maxSeqLen, Parameters::DBTYPE_AMINO_ACIDS, (const BaseMatrix *) &subMat3Di, 0, false, par.compBiasCorrection);
        std::string backtrace;
        char buffer[1024+32768];
        std::string resultBuffer;

        Coordinate16 qcoords;
        Coordinate16 tcoords;

        // write output file

#pragma omp for schedule(dynamic, 1)


        std::vector<ChainToChainAln> resultVector;
        for (size_t i = 0; i < alnDbr.getSize(); i++) {
            size_t queryKey = alnDbr.getDbKey(i);
            const unsigned queryComplexId = qLookup.at(queryKey);

            char *data = alnDbr.getData(i, thread_idx);
            unsigned int queryId = qdbr3Di.sequenceReader->getId(queryKey);
            if (*data == '\0') {
                continue;
            }
            char *querySeqAA = qdbrAA.sequenceReader->getData(queryId, thread_idx);
            char *querySeq3Di = qdbr3Di.sequenceReader->getData(queryId, thread_idx);
            unsigned int querySeqLen = qdbr3Di.sequenceReader->getSeqLen(queryId);
            qSeq3Di.mapSequence(i, queryKey, querySeq3Di, querySeqLen);
            qSeqAA.mapSequence(i, queryKey, querySeqAA, querySeqLen);
            size_t qId = qcadbr->sequenceReader->getId(queryKey);
            char *qcadata = qcadbr->sequenceReader->getData(qId, thread_idx);
            float *queryCaData = (float *) qcadata;
            if (qcadbr->getDbtype() == LocalParameters::DBTYPE_CA_ALPHA_F16) {
                qcoords.read(qcadata, qSeq3Di.L);
                queryCaData = qcoords.getBuffer();
            }
            Chain qChain = Chain(queryComplexId, queryKey, queryCaData, qSeq3Di.L);
            tmaligner->initQuery(queryCaData, &queryCaData[qSeq3Di.L], &queryCaData[qSeq3Di.L + qSeq3Di.L],
                                 NULL, qSeq3Di.L);

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
                Chain dbChain = Chain(dbComplexId, dbKey, targetCaData, alnResult.dbLen);
                TMaligner::TMscoreResult tmres = tmaligner->computeTMscore(targetCaData, &targetCaData[alnResult.dbLen], &targetCaData[alnResult.dbLen+alnResult.dbLen], alnResult.dbLen,alnResult.qStartPos, alnResult.dbStartPos, Matcher::uncompressAlignment(alnResult.backtrace));
                double tmScore = tmres.tmscore;
                ChainToChainAln chainToChainAln = ChainToChainAln(qChain, dbChain, tmScore, alnResult.backtrace);
                if (tmScore > 0) {
                    resultVector.emplace_back(chainToChainAln);
                }
                data = Util::skipLine(data);
            }
        }
        std::sort(resultVector.begin(), resultVector.end(), compareChainsWithTmByQueryChain);
        std::vector<ChainToChainAln> tempVec1;
        std::vector<ChainToChainAln> tempVec2;
        unsigned int prevQueryComplexID = resultVector[0].qChain.complexid;
        unsigned int prevQueryKey = resultVector[0].qChain.chainKey;
        unsigned int prevDbComplexID = resultVector[0].dbChain.complexid;
        unsigned int prevDbKey = resultVector[0].dbChain.chainKey;
        Chain bestDbChain;
        double bestTmScore = 0;
        std::string bestBacktrace;
        for (size_t i=0; i<resultVector.size(); i++) {
            ChainToChainAln chainToChainAln = resultVector[i];
            Chain qChain = chainToChainAln.qChain;
            Chain dbChain = chainToChainAln.dbChain;
            unsigned int currQueryComplexID = qChain.complexid;
            unsigned int currQueryKey = qChain.chainKey;
            unsigned int currDbComplexID = dbChain.complexid;
//            unsigned int currDbKey = dbChain.chainKey;
            double tmScore = chainToChainAln.tmScore;
            std::string backtrace;
            if (currQueryComplexID != prevQueryComplexID) {
                printOut(tempVec1, tempVec2);
                prevQueryComplexID = currQueryComplexID;
                prevQueryKey = currQueryKey;
                prevDbComplexID = currDbComplexID;

                bestDbChain = dbChain;
                bestTmScore = tmScore;
                bestBacktrace = backtrace;
                tempVec1.clear();
                tempVec2.clear();

            } else if (currQueryKey != prevQueryKey){
                tempVec2 = update(tempVec1, tempVec2);
                prevQueryKey = currQueryKey;
                prevDbComplexID = currDbComplexID;

                bestDbChain = dbChain;
                bestTmScore = tmScore;
                bestBacktrace = backtrace;
                tempVec1.clear();

            } else if (currDbComplexID != prevDbComplexID) {
                tempVec1.emplace_back(ChainToChainAln(qChain, bestDbChain, bestTmScore, bestBacktrace));
                prevQueryKey = currQueryKey;
                prevDbComplexID = currDbComplexID;

                bestDbChain = dbChain;
                bestTmScore = tmScore;
                bestBacktrace = backtrace;
            } else if (tmScore > bestTmScore) {
                bestDbChain = dbChain;
                bestTmScore = tmScore;
                bestBacktrace = backtrace;
            }
        }
        printOut(tempVec1, tempVec2);

        if(needTMaligner){
            delete tmaligner;
        }
    }

    free(tinySubMatAA);
    free(tinySubMat3Di);

    dbw.close();
    alnDbr.close();
    if (sameDB == false) {
        delete t3DiDbr;
        delete tAADbr;
    }
    return EXIT_SUCCESS;
}