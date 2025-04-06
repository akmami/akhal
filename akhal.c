#include "validate.h"
#include "analyze.h"
#include "extract.h"
#include "gaf2sam.h"

void printUsage() {
    fprintf(stderr, "[ERROR] Usage: ./akhal PROGRAM ...\n");
    fprintf(stderr, "[OPTION]\n");
    fprintf(stderr, "\tparse\tParse r/GFA and make sure that it is correct.\n");
    fprintf(stderr, "\tstats\tStatistics about r/GFA.\n");
    fprintf(stderr, "\textract\tExtract reference genome.\n");
    fprintf(stderr, "\tgaf2sam\tGAF to SAM conversion.\n");
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printUsage();
        return -1;
    }

    if (strcmp(argv[1], "parse") == 0) {
        validate(argc, argv);
    } else if (strcmp(argv[1], "stats") == 0) {
        analyze(argc, argv);
    } else if (strcmp(argv[1], "extract") == 0) {
        extract(argc, argv);
    } else if (strcmp(argv[1], "gaf2sam") == 0) {
        gaf2sam(argc, argv);
    } else {
        printUsage();
        return -1;
    }

    return 0;
}