#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "define.h"

//*************File management functions*************

void readTable(double* table,char* fileName);
// Read the indicated data file and store it to the given table. Size is not controlled, to avoid error the table should be large enough

void saveTable(double *table, char *name, char* path, int tabSize);
// Save general 1D table containing doubles.

void save2DTable(double **table, char *name, char* path, int nRow, int nCol, int skipR, int skipC, int startC);
// Save general 2D table containing doubles, skipR and skipC allow to save respectively only every "skipR" lines and every "skipC" column, startC gives the startig column to be saved

void save2DTable_top(double **table, char *name, char* path, double* thickness, int nrow, int ncol);
// TO DO

#endif  /* !IO_H */
