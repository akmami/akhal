#ifndef __IO_H__
#define __IO_H__

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

FILE *io_open(const char* file_path, char **line, int cap);
void io_close(FILE *file, char **str);
int io_read(FILE *file, char **str, size_t *cap);

#endif