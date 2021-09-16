#ifndef IO_H
#define IO_H

#include "define.h"
#include "structures.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
//*************File management functions*************

bool readTable(double *table, const char *const fileName);
// Read the indicated data file and store it to the given table. Size is not
// controlled, to avoid error the table should be large enough

void saveTable(double *table, const char *const name, const char *const path,
               const int tabSize);
// Save general 1D table containing doubles.

void save2DTable(double **table, const char *const name, const char *const path,
                 const int nRow, const int nCol, const int skipR,
                 const int skipC, const int startC);
// Save general 2D table containing doubles, skipR and skipC allow to save
// respectively only every "skipR" lines and every "skipC" column, startC gives
// the startig column to be saved

void save2DTable_top(double **table, const char *const name,
                     const char *const path, const double *const thickness,
                     const int nrow, const int ncol);

bool copyFile(const char *const inFile, const char *const outFile);

// TO DO
bool readInitFile(model_parameters *const params, const char *const fileName);

void removeSpaces(char *const source);

void removeComments(char *const source);

bool parseLine(char *line, char **arg, char **val);

bool setParameter(model_parameters *const params, const char *const arg,
                  char *const val);

bool checkParameter(model_parameters *const params);

void upper_string(char *s);

bool checkOrCreateDir(char *const path);

bool createOutputDirs(model_parameters *const params);

#endif /* !IO_H */
