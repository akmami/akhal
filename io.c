#include "io.h"

FILE *io_open(const char* file_path, char **line, int cap) {

    FILE* file = fopen(file_path, "r");
    if (!file) { fprintf(stderr, "[ERROR] Failed to open file."); exit(EXIT_FAILURE); }
    
    *line = (char *)malloc(cap);

    if (!(*line)) { fprintf(stderr, "[ERROR] Failed to allocate memory."); exit(EXIT_FAILURE); }

    return file;
}

void io_close(FILE *file, char **str) {
    if (file)
        fclose(file);
    if (*str) {
        free(*str);
    }
    file = NULL;
}

int io_read(FILE *file, char **str, size_t *cap) {

    char *line = *str;
    size_t line_len = *cap;

    if (fgets(line, line_len, file)) {
        size_t len = strlen(line);

        // read the line and fit it to `line`
        while (len == line_len - 1 && line[len - 1] != '\n') {
            line_len *= 2;
            *cap = line_len;
            char *temp_line = (char *)realloc(line, line_len);
            if (!temp_line) { fprintf(stderr, "[ERROR] Memory reallocation failed.\n"); return 0; }
            line = temp_line;
            *str = temp_line;

            if (fgets(line + len, line_len - len, file) == NULL) { break; }
            len = strlen(line);
        }
        if (len && line[len - 1] == '\n') { line[len - 1] = '\0'; len--; }

        return len;
    }

    return 0;
}