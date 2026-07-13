#ifndef __IO_H__
#define __IO_H__

#include "struct_def.h"
#include "utils.h"
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

FILE *io_open(const char* file_path, char **line, int cap);
void io_close(FILE *file, char **str);
int io_read(FILE *file, char **str, size_t *cap);
void write_sam_record(FILE *out_sam, alignment *aln, char *ops, int c_size, const char *rname, int pos, int simplify);

#endif