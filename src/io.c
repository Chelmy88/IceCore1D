#include "io.h"
#include <dirent.h>
#include <errno.h>
#include <limits.h>
#include <string.h>
#include <unistd.h>
//*************File management functions*************

bool readTable(double *table, const char *const fileName) {
  // Read a table from the file called "filename" and store it into the double
  // table called "table". Table should be 1D. "filename" can contain directory
  // path.
  FILE *fp;
  int li = 0;
  double a = 0;
  if ((fp = fopen(fileName, "r")) == NULL) {
    printf("[E] annot open file: %s\n", fileName);
    return false;
  } else {
    while (fscanf(fp, "%lf", &a) == 1) {
      table[li] = a;
      li++;
    }
    if (fclose(fp) == 0) {
      printf("[I] File imported successfully (%d data) and closed: %s \n", li,
             fileName);
      return true;
    } else {
      printf("[E] Not able to close: %s \n", fileName);
      return false;
    }
  }
}

void saveTable(double *table, const char *const name, const char *const path,
               const int tabSize) {
  // Read the 1D double table called "table" and store it a file called
  // "filename".  Table should be tab delimited. "filename" can contain
  // directory path. The default path is a folder called "export". The table
  // size should be passed as parameter.
  FILE *fp;
  int li = 0;
  char full_path[400] = "";
  sprintf(full_path, "%s/%s", path, name);
  if ((fp = fopen(full_path, "w+")) == NULL) {
    printf("[E] Cannot open file: %s\n", full_path);
  } else {
    printf("[I] File oppened:%s\n[I] ...writing...\n", full_path);
    for (li = 0; li < tabSize; li++) {
      fprintf(fp, "%f \n", table[li]);
    }
    fclose(fp);
    printf("[I] File closed: %s \n", full_path);
  }
}

void save2DTable(double **table, const char *const name, const char *const path,
                 const int nRow, const int nCol, const int skipR,
                 const int skipC, const int startC) {
  // Read the 2D double table called "table" and store it a tab delimited
  // file called "filename". The default path is a folder called "export".
  // The table dimensions should be passed as parameter.
  FILE *fp;
  int li = 0;
  int co = 0;
  char full_path[400] = "";
  sprintf(full_path, "%s/%s", path, name);
  if ((fp = fopen(full_path, "w+")) == NULL) {
    printf("[E] Cannot open file: %s\n", full_path);
  } else {
    printf("[I] File oppened:%s\n[I] ...writing...\n", full_path);
    for (li = 0; li < nRow; li += skipR) {
      for (co = startC; co < nCol; co += skipC) {
        if (co != startC) {
          fprintf(fp, "\t");
        }
        if (table[li][co] > 0) {
          fprintf(fp, "%f", table[li][co]);
        } else {
          fprintf(fp, "NaN");
        }
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
    printf("[I] File closed: %s \n", full_path);
  }
}

// save2DTable_top(density,fileName,path,iceThickness,100,T);

void save2DTable_top(double **table, const char *const name,
                     const char *const path, const double *const thickness,
                     const int nRow, const int nCol) {
  // Read the 2D double table called "table" and store it a tab delimited
  // file called "filename". The default path is a folder called "export".
  // The table dimensions should be passed as parameter.
  FILE *fp;
  int li = 0;
  int co = 0;
  int startC = 0;
  char full_path[400] = "";
  sprintf(full_path, "%s/%s", path, name);
  if ((fp = fopen(full_path, "w+")) == NULL) {
    printf("[I] Cannot open file: %s\n", full_path);
  } else {
    printf("[I] File oppened:%s\n[I] ...writing...\n", full_path);
    for (li = 0; li < nRow; li += 1) {
      for (co = startC; co < nCol; co += 1) {
        int real_li = 0;
        if (co != startC) {
          fprintf(fp, "\t");
          real_li = (int)(thickness[co - 1]) - li;
        } else {
          real_li = (int)(thickness[co]) - li;
        }
        if (table[real_li][co] > 0) {
          fprintf(fp, "%f", table[real_li][co]);
        } else {
          fprintf(fp, "NaN");
        }
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
    printf("[I] File closed: %s \n", full_path);
  }
}

bool copyFile(const char *const inFile, const char *const outFile) {
  FILE *fs, *ft;
  char ch;
  fs = fopen(inFile, "r");
  if (fs == NULL) {
    printf("[E] Error in opening source file: %s", inFile);
    return false;
  }
  ft = fopen(outFile, "w");
  if (ft == NULL) {
    printf("[E] Error in opening destination file: %s", outFile);
    fclose(fs);
    return false;
  }
  while (1) {
    ch = fgetc(fs);
    if (ch == EOF) {
      break;
    } else {
      fputc(ch, ft);
    }
  }
  printf("[I] File copied successfully");
  fclose(fs);
  fclose(ft);
  return true;
}

bool readInitFile(model_parameters *const params, const char *const fileName) {
  FILE *fp;
  if ((fp = fopen(fileName, "r")) == NULL) {
    printf("[E] Cannot open file: %s\n", fileName);
    return false;
  }
  char line[BUFSIZ];
  bool state = true;
  while (fgets(line, BUFSIZ, fp) != NULL) {
    removeComments(line);
    if (line[0] == '\0') {
      continue;
    }
    char *arg = NULL;
    char *val = NULL;
    if (!parseLine(line, &arg, &val)) {
      state = false;
    } else if (!setParameter(params, arg, val)) {
      state = false;
    }
  }
  fclose(fp);
  if (state) {
    if (!checkParameter(params)) {
      state = false;
    }
  }
  if (state) {
    params->values.tot = params->values.mw_n * params->values.QG_n *
                         params->values.TCor_n * params->values.TCor2_n *
                         params->values.PCor_n * params->values.deltaH_n *
                         params->values.len_n * params->values.flat_n;

    printf("[I] Ini file %s read successfully\n", fileName);
    return true;
  } else {
    return false;
  }
}

void removeSpaces(char *source) {
  char *i = source;
  char *j = source;
  do {
    *i = *j;
    if (*i != ' ')
      i++;
  } while (*j++ != 0);
}
void removeComments(char *const source) {
  char *ptr;
  ptr = strchr(source, '#');
  if (ptr != NULL) {
    *ptr = '\0';
  }
  ptr = strchr(source, ';');
  if (ptr != NULL) {
    *ptr = '\0';
  }
  ptr = strchr(source, '\n');
  if (ptr != NULL) {
    *ptr = '\0';
  }
}

bool parseLine(char *line, char **arg, char **val) {
  char *pch;
  pch = strtok(line, ":");
  if (pch != NULL) {
    *arg = pch;
    removeSpaces(*arg);
    upper_string(*arg);
    pch = strtok(NULL, ":");
    if (pch != NULL) {
      *val = pch;
      pch = strtok(NULL, ":");
      if (pch != NULL) {
        printf("[E] Two ':' delimiter found in argument line %s\n", line);
        return false;
      }
    } else {
      printf("[E] No ':' delimiter found in argument line %s\n", line);
      return false;
    }
  } else {
    printf("[E] No ':' delimiter found in argument line %s\n", line);
    return false;
  }
  return true;
}

// Rather long function that test for the various allowed keys
// and strores the corresponding values
bool setParameter(model_parameters *const params, const char *const arg,
                  char *const val) {
  if (strcmp(arg, "Z") == 0) {
    removeSpaces(val);
    if (atoi(val) > 0) {
      params->Z = atoi(val);
    } else {
      printf("[E] Wrong parameter %s for Z. Must be positive integer.\n", val);
      return false;
    }
  }

  else if (strcmp(arg, "T") == 0) {
    removeSpaces(val);
    if (atoi(val) > 0) {
      params->T = atoi(val);
    } else {
      printf("[E] Wrong parameter %s for Z. Must be positive integer.\n", val);
      return false;
    }
  }

  else if (strcmp(arg, "S") == 0) {
    removeSpaces(val);
    if (atoi(val) > 0) {
      params->S = atoi(val);
    } else {
      printf("[E] Wrong parameter %s for Z. Must be positive integer\n", val);
      return false;
    }
  }

  else if (strcmp(arg, "OUTPUT_PATH") == 0) {
    removeSpaces(val);
    params->OUTPUT_PATH = (char *)malloc(sizeof(char) * (strlen(val) + 1));
    strcpy(params->OUTPUT_PATH, val);


  }

  else if (strcmp(arg, "TEMPERATURE_FILE") == 0) {
    removeSpaces(val);
    if (access(val, R_OK) != -1) {
      params->TEMPERATURE_FILE =
          (char *)malloc(sizeof(char) * (strlen(val) + 1));
      strcpy(params->TEMPERATURE_FILE, val);
    } else {
      printf("[E] Wrong parameter %s for TEMPERATURE_FILE, file not exist or "
             "not readable\n",
             val);
      return false;
    }
  }

  else if (strcmp(arg, "ACCUMULATION_FILE") == 0) {
    removeSpaces(val);
    if (access(val, R_OK) != -1) {
      params->ACCUMULATION_FILE =
          (char *)malloc(sizeof(char) * (strlen(val) + 1));
      strcpy(params->ACCUMULATION_FILE, val);
    } else {
      printf("[E] Wrong parameter %s for ACCUMULATION_FILE, file not exist or "
             "not readable\n",
             val);
      return false;
    }
  }

  else if (strcmp(arg, "ICE_THICKNESS_FILE") == 0) {
    removeSpaces(val);
    if (access(val, R_OK) != -1) {
      params->ICE_THICKNESS_FILE =
          (char *)malloc(sizeof(char) * (strlen(val) + 1));
      strcpy(params->ICE_THICKNESS_FILE, val);
    } else {
      printf("[E] Wrong parameter %s for ICE_THICKNESS_FILE, file not exist or "
             "not readable\n",
             val);
      return false;
    }
  }

  else if (strcmp(arg, "AGE_FILE") == 0) {
    removeSpaces(val);
    if (access(val, R_OK) != -1) {
      params->AGE_FILE = (char *)malloc(sizeof(char) * (strlen(val) + 1));
      strcpy(params->AGE_FILE, val);
    } else {
      printf("[E] Wrong parameter %s for AGE_FILE, file not exist or not "
             "readable\n",
             val);
      return false;
    }
  }

  else if (strcmp(arg, "BOREHOLE_TEMPERATURE_FILE") == 0) {
    removeSpaces(val);
    if (access(val, R_OK) != -1) {
      params->BOREHOLE_TEMPERATURE_FILE =
          (char *)malloc(sizeof(char) * (strlen(val) + 1));
      strcpy(params->BOREHOLE_TEMPERATURE_FILE, val);
    } else {
      printf("[E] Wrong parameter %s for BOREHOLE_TEMPERATURE_FILE, file not "
             "exist or not readable\n",
             val);
      return false;
    }
  }

  else if (strcmp(arg, "SAVE_TYPE") == 0) {
    removeSpaces(val);
    if (strcmp(val, "MATRIX") == 0) {
      params->SAVE_TYPE = ST_MATRIX;
    } else if (strcmp(val, "VECTOR") == -0) {
      params->SAVE_TYPE = ST_VECTOR;
    } else {
      printf("[E] Unknown value %s for SAVE_TYPE\n", val);
      return false;
    }
  }

  else if (strcmp(arg, "SCHEME") == 0) {
    removeSpaces(val);
    if (strcmp(val, "CN") == 0) {
      params->SCHEME = SC_CN;
    } else if (strcmp(val, "EXPL") == 0) {
      params->SCHEME = SC_EXPL;
    } else {
      printf("[E] Unknown value %s for SAVE_TYPE\n", val);
      return false;
    }
  }

  else if (strcmp(arg, "RHO_SNOW") == 0) {
    removeSpaces(val);
    if (atoi(val) > 0) {
      params->RHO_SNOW = atoi(val);
    } else {
      printf("[E] Wrong parameter %s for RHO_SNOW. Must be positive integer\n",
             val);
      return false;
    }
  }

  else if (strcmp(arg, "THERMAL_ICE") == 0) {
    removeSpaces(val);
    if (strcmp(val, "CP") == 0) {
      params->THERMAL_ICE = TI_CP;
    } else if (strcmp(val, "GO") == 0) {
      params->THERMAL_ICE = TI_GO;
    } else {
      printf("[E] Unknown value %s for THERMAL\n", val);
      return false;
    }
  }

  else if (strcmp(arg, "THERMAL_FIRN") == 0) {
    removeSpaces(val);
    if (strcmp(val, "SC") == 0) {
      params->THERMAL_FIRN = TF_SC;
    } else if (strcmp(val, "CP") == 0) {
      params->THERMAL_FIRN = TF_CP;
    } else if (strcmp(val, "CP_LIN") == 0) {
      params->THERMAL_FIRN = TF_CP_LIN;
    } else if (strcmp(val, "SC_LIN") == 0) {
      params->THERMAL_FIRN = TF_SC_LIN;
    } else if (strcmp(val, "CP_AL") == 0) {
      params->THERMAL_FIRN = TF_CP_AL;
    } else if (strcmp(val, "SC_AL") == 0) {
      params->THERMAL_FIRN = TF_SC_AL;
    } else if (strcmp(val, "CP_ST") == 0) {
      params->THERMAL_FIRN = TF_CP_ST;
    } else if (strcmp(val, "SC_ST") == 0) {
      params->THERMAL_FIRN = TF_SC_ST;
    } else if (strcmp(val, "CP_WE_ADD") == 0) {
      params->THERMAL_FIRN = TF_CP_WE_ADD;
    } else if (strcmp(val, "CP_WE_LIN") == 0) {
      params->THERMAL_FIRN = TF_CP_WE_LIN;
    } else if (strcmp(val, "SC_WE_ADD") == 0) {
      params->THERMAL_FIRN = TF_SC_WE_ADD;
    } else if (strcmp(val, "SC_WE_LIN") == 0) {
      params->THERMAL_FIRN = TF_SC_WE_LIN;
    } else {
      printf("[E] Unknown value %s for THERMAL_FIRN\n", val);
      return false;
    }
  }

  else if (strcmp(arg, "RHO_FIRN") == 0) {
    removeSpaces(val);
    if (strcmp(val, "HL") == 0) {
      params->RHO_FIRN = RF_HL;
    } else if (strcmp(val, "CONST") == 0) {
      params->RHO_FIRN = RF_CONST;
    } else {
      printf("[E] Unknown value %s for RHO_FIRN\n", val);
      return false;
    }
  }

  else if (strcmp(arg, "HEAT_CAPACITY") == 0) {
    removeSpaces(val);
    if (strcmp(val, "CP") == 0) {
      params->HEAT_CAPACITY = CP_CP;
    } else if (strcmp(val, "CP_AL") == 0) {
      params->HEAT_CAPACITY = CP_CP_AL;
    } else {
      printf("[E] Unknown value %s for HEAT_CAPACITY\n", val);
      return false;
    }
  }

  else if (strcmp(arg, "VERTICAL_PROFILE") == 0) {
    removeSpaces(val);
    if (strcmp(val, "FI") == 0) {
      params->VERTICAL_PROFILE = VP_FI;
    } else if (strcmp(val, "PA") == 0) {
      params->VERTICAL_PROFILE = VP_PA;
    } else {
      printf("[E] Unknown value %s for VERTICAL1\n", val);
      return false;
    }
  }

  else if (strcmp(arg, "INTERNAL_ENERGY") == 0) {
    removeSpaces(val);
    if (strcmp(val, "ON") == 0) {
      params->INTERNAL_ENERGY = IE_ON;
    } else if (strcmp(val, "OFF") == 0) {
      params->INTERNAL_ENERGY = IE_OFF;
    } else {
      printf("[E] Unknown value %s for INTERNAL_ENERGY\n", val);
      return false;
    }
  } else if (strcmp(arg, "MELTING") == 0) {
    removeSpaces(val);
    if (strcmp(val, "FREE_MELT") == 0) {
      params->MELTING = ME_FREE_MELT;
    } else if (strcmp(val, "FREEZING_NO_ICE") == 0) {
      params->MELTING = ME_FREEZING_NO_ICE;
    } else if (strcmp(val, "FREEZING") == 0) {
      params->MELTING = ME_FREEZING;
    } else {
      printf("[E] Unknown value %s for MELTING\n", val);
      return false;
    }
  }

  else if (strcmp(arg, "MW") == 0) {
    char *token;
    // make a copy for the first pass
    char val_copy[40];
    strcpy(val_copy, val);
    // get the first token
    token = strtok(val_copy, " ");
    // walk through other tokens and store the size
    // also validate the input
    size_t n = 0;
    while (token != NULL) {
      if (atof(token) > 0 || strcmp(token, "0") == 0 ||
          strcmp(token, "0.0") == 0) {
        ++n;
      } else {
        printf("[E] Unknown value %s for MW\n", token);
        return false;
      }
      token = strtok(NULL, " ");
    }
    // Now alloc the array
    params->values.mw = (double *)malloc(sizeof(double) * n);
    // and store the data in it
    token = strtok(val, " ");
    n = 0;
    while (token != NULL) {
      params->values.mw[n] = atof(token);
      token = strtok(NULL, " ");
      ++n;
    }
    params->values.mw_n = n;
  }

  else if (strcmp(arg, "QG") == 0) {
    char *token;
    // make a copy for the first pass
    char val_copy[40];
    strcpy(val_copy, val);
    // get the first token
    token = strtok(val_copy, " ");
    // walk through other tokens and store the size
    // also validate the input
    size_t n = 0;
    while (token != NULL) {
      if (atof(token) != 0 || strcmp(token, "0") == 0 ||
          strcmp(token, "0.0") == 0) {
        ++n;
      } else {
        printf("[E] Unknown value %s for QG\n", token);
        return false;
      }
      token = strtok(NULL, " ");
    }
    // Now alloc the array
    params->values.QG = (double *)malloc(sizeof(double) * n);
    // and store the data in it
    token = strtok(val, " ");
    n = 0;
    while (token != NULL) {
      params->values.QG[n] = atof(token);
      token = strtok(NULL, " ");
      ++n;
    }
    params->values.QG_n = n;
  }

  else if (strcmp(arg, "TCOR") == 0) {
    char *token;
    // make a copy for the first pass
    char val_copy[40];
    strcpy(val_copy, val);
    // get the first token
    token = strtok(val_copy, " ");
    // walk through other tokens and store the size
    // also validate the input
    size_t n = 0;
    while (token != NULL) {
      if (atof(token) > 0 || strcmp(token, "0") == 0 ||
          strcmp(token, "0.0") == 0) {
        ++n;
      } else {
        printf("[E] Unknown value %s for TCOR\n", token);
        return false;
      }
      token = strtok(NULL, " ");
    }
    // Now alloc the array
    params->values.TCor = (double *)malloc(sizeof(double) * n);
    // and store the data in it
    token = strtok(val, " ");
    n = 0;
    while (token != NULL) {
      params->values.TCor[n] = atof(token);
      token = strtok(NULL, " ");
      ++n;
    }
    params->values.TCor_n = n;
  }

  else if (strcmp(arg, "TCOR2") == 0) {
    char *token;
    // make a copy for the first pass
    char val_copy[40];
    strcpy(val_copy, val);
    // get the first token
    token = strtok(val_copy, " ");
    // walk through other tokens and store the size
    // also validate the input
    size_t n = 0;
    while (token != NULL) {
      if (atof(token) > 0 || strcmp(token, "0") == 0 ||
          strcmp(token, "0.0") == 0) {
        ++n;
      } else {
        printf("[E] Unknown value %s for TCOR2\n", token);
        return false;
      }
      token = strtok(NULL, " ");
    }
    // Now alloc the array
    params->values.TCor2 = (double *)malloc(sizeof(double) * n);
    // and store the data in it
    token = strtok(val, " ");
    n = 0;
    while (token != NULL) {
      params->values.TCor2[n] = atof(token);
      token = strtok(NULL, " ");
      ++n;
    }
    params->values.TCor2_n = n;
  }

  else if (strcmp(arg, "PCOR") == 0) {
    char *token;
    // make a copy for the first pass
    char val_copy[40];
    strcpy(val_copy, val);
    // get the first token
    token = strtok(val_copy, " ");
    // walk through other tokens and store the size
    // also validate the input
    size_t n = 0;
    while (token != NULL) {
      if (atof(token) > 0 || strcmp(token, "0") == 0 ||
          strcmp(token, "0.0") == 0) {
        ++n;
      } else {
        printf("[E] Unknown value %s for PCOR\n", token);
        return false;
      }
      token = strtok(NULL, " ");
    }
    // Now alloc the array
    params->values.PCor = (double *)malloc(sizeof(double) * n);
    // and store the data in it
    token = strtok(val, " ");
    n = 0;
    while (token != NULL) {
      params->values.PCor[n] = atof(token);
      token = strtok(NULL, " ");
      ++n;
    }
    params->values.PCor_n = n;
  }

  else if (strcmp(arg, "DELTAH") == 0) {
    char *token;
    // make a copy for the first pass
    char val_copy[40];
    strcpy(val_copy, val);
    // get the first token
    token = strtok(val_copy, " ");
    // walk through other tokens and store the size
    // also validate the input
    size_t n = 0;
    while (token != NULL) {
      if (atof(token) > 0 || strcmp(token, "0") == 0 ||
          strcmp(token, "0.0") == 0) {
        ++n;
      } else {
        printf("[E] Unknown value %s for DELTAH\n", token);
        return false;
      }
      token = strtok(NULL, " ");
    }
    // Now alloc the array
    params->values.deltaH = (double *)malloc(sizeof(double) * n);
    // and store the data in it
    token = strtok(val, " ");
    n = 0;
    while (token != NULL) {
      params->values.deltaH[n] = atof(token);
      token = strtok(NULL, " ");
      ++n;
    }
    params->values.deltaH_n = n;
  } else if (strcmp(arg, "LEN") == 0) {
    char *token;
    // make a copy for the first pass
    char val_copy[40];
    strcpy(val_copy, val);
    // get the first token
    token = strtok(val_copy, " ");
    // walk through other tokens and store the size
    // also validate the input
    size_t n = 0;
    while (token != NULL) {
      if (atof(token) > 0 || strcmp(token, "0") == 0 ||
          strcmp(token, "0.0") == 0) {
        ++n;
      } else {
        printf("[E] Unknown value %s for LEN\n", token);
        return false;
      }
      token = strtok(NULL, " ");
    }
    // Now alloc the array
    params->values.len = (double *)malloc(sizeof(double) * n);
    // and store the data in it
    token = strtok(val, " ");
    n = 0;
    while (token != NULL) {
      params->values.len[n] = atof(token);
      token = strtok(NULL, " ");
      ++n;
    }
    params->values.len_n = n;
  }

  else if (strcmp(arg, "FLAT") == 0) {
    char *token;
    // make a copy for the first pass
    char val_copy[40];
    strcpy(val_copy, val);
    // get the first token
    token = strtok(val_copy, " ");
    // walk through other tokens and store the size
    // also validate the input
    size_t n = 0;
    while (token != NULL) {
      if (atof(token) > 0 || strcmp(token, "0") == 0 ||
          strcmp(token, "0.0") == 0) {
        ++n;
      } else {
        printf("[E] Unknown value %s for FLAT\n", token);
        return false;
      }
      token = strtok(NULL, " ");
    }
    // Now alloc the array
    params->values.flat = (double *)malloc(sizeof(double) * n);
    // and store the data in it
    token = strtok(val, " ");
    n = 0;
    while (token != NULL) {
      params->values.flat[n] = atof(token);
      token = strtok(NULL, " ");
      ++n;
    }
    params->values.flat_n = n;
  } else {
    printf("[E] Unkown parameter %s find in ini file\n", arg);
    return false;
  }
  return true;
}

bool checkParameter(model_parameters *const params) {

  bool state = true;
  if (params->Z == -1) {
    printf("[E] Unset parameter Z. This parameter is mandatory\n");
    state = false;
  }
  if (params->T == -1) {
    printf("[E] Unset parameter T. This parameter is mandatory\n");
    state = false;
  }
  if (params->S == -1) {
    params->S = 1500;
    printf("[W] Unset parameter S, the default value of 1500 is used\n");
  }
  if (params->RHO_SNOW == -1) {
    params->RHO_SNOW = 350;
    printf("[W] Unset parameter RHO_SNOW, the default value of 350 is used\n");
  }
  if (!params->OUTPUT_PATH) {
    printf("[E] Unset parameter OUTPUT_PATH. This parameter is mandatory\n");
    state = false;
  }
  if (!params->TEMPERATURE_FILE) {
    printf(
        "[E] Unset parameter TEMPERATURE_FILE. This parameter is mandatory\n");
    state = false;
  }
  if (!params->ACCUMULATION_FILE) {
    printf(
        "[E] Unset parameter ACCUMULATION_FILE. This parameter is mandatory\n");
    state = false;
  }
  if (!params->ICE_THICKNESS_FILE) {
    printf("[E] Unset parameter ICE_THICKNESS_FILE. This parameter is "
           "mandatory\n");
    state = false;
  }
  if (!params->AGE_FILE) {
    printf("[E] Unset parameter AGE_FILE. This parameter is mandatory\n");
    state = false;
  }
  if (!params->BOREHOLE_TEMPERATURE_FILE) {
    printf("[E] Unset parameter BOREHOLE_TEMPERATURE_FILE. This parameter is "
           "mandatory\n");
    state = false;
  }
  if (params->SAVE_TYPE == ST_UNSET) {
    params->SAVE_TYPE = ST_VECTOR;
    printf(
        "[W] Unset parameter SAVE_TYPE, the default value of VECTOR is used\n");
  }
  if (params->SCHEME == SC_UNSET) {
    params->SCHEME = SC_CN;
    printf("[W] Unset parameter SCHEME, the default value of CN is used\n");
  }
  if (params->THERMAL_ICE == TI_UNSET) {
    params->THERMAL_ICE = TI_CP;
    printf("[W] Unset parameter THERMAL, the default value of CP is used\n");
  }
  if (params->THERMAL_FIRN == TF_UNSET) {
    params->THERMAL_FIRN = TF_SC;
    printf("[W] Unset parameter FIRN, the default value of SC is used\n");
  }
  if (params->RHO_FIRN == RF_UNSET) {
    params->RHO_FIRN = RF_HL;
    printf("[W] Unset parameter RHO_FIRN, the default value of RF_HL is used\n");
  }
  if (params->VERTICAL_PROFILE == VP_UNSET) {
    params->VERTICAL_PROFILE = VP_FI;
    printf("[W] Unset parameter VERTICAL_PROFILE, the default value of FI is "
           "used\n");
  }
  if (params->INTERNAL_ENERGY == IE_UNSET) {
    params->INTERNAL_ENERGY = IE_OFF;
    printf("[W] Unset parameter INTERNAL_ENERGY, the default value of OFF is "
           "used\n");
  }
  if (params->MELTING == ME_UNSET) {
    params->MELTING = ME_FREE_MELT;
    printf("[W] Unset parameter MELTING, the default value of FREE_MELT is "
           "used\n");
  }
  if (!params->values.mw) {
    printf("[E] Unset parameter MW. This parameter is mandatory\n");
    state = false;
  }
  if (!params->values.QG) {
    printf("[E] Unset parameter QG. This parameter is mandatory\n");
    state = false;
  }
  if (!params->values.TCor) {
    params->values.TCor = (double *)malloc(sizeof(double));
    params->values.TCor[0] = 0.;
    params->values.TCor_n = 1;
    printf("[W] Unset parameter TCOR, the default value of 0 is used\n");
  }
  if (!params->values.TCor2) {
    params->values.TCor2 = (double *)malloc(sizeof(double));
    params->values.TCor2[0] = 0.;
    params->values.TCor2_n = 1;
    printf("[W] Unset parameter TCOR2, the default value of 0 is used\n");
  }
  if (!params->values.PCor) {
    params->values.PCor = (double *)malloc(sizeof(double));
    params->values.PCor[0] = 0.;
    params->values.PCor_n = 1;
    printf("[W] Unset parameter PCOR, the default value of 0 is used\n");
  }
  if (!params->values.deltaH) {
    params->values.deltaH = (double *)malloc(sizeof(double));
    params->values.deltaH[0] = 0.;
    params->values.deltaH_n = 1;
    printf("[W] Unset parameter DELTAH, the default value of 0 is used\n");
  }
  if (!params->values.len) {
    params->values.len = (double *)malloc(sizeof(double));
    params->values.len[0] = 0.;
    params->values.len_n = 1;
    printf("[W] Unset parameter LEN, the default value of 0 is used\n");
  }
  if (!params->values.flat) {
    params->values.flat = (double *)malloc(sizeof(double));
    params->values.flat[0] = 0.;
    params->values.flat_n = 1;
    printf("[W] Unset parameter FLAT, the default value of 0 is used\n");
  }
  return (state);
}

void upper_string(char *s) {
  size_t c = 0;
  while (s[c] != '\0') {
    if (s[c] >= 'a' && s[c] <= 'z') {
      s[c] = s[c] - 32;
    }
    ++c;
  }
}

bool checkOrCreateDir(char *const path)
{
  struct stat st = {0};
  if (stat(path, &st) == 0 && S_ISDIR(st.st_mode)) {
    if (! S_ISDIR(st.st_mode)) {
      printf("[E] Impossible to create %s directory, path exist but is not a directory\n",
      path);
      return(false);
    }
    return(true);
  }
  if(mkdir(path, 0700)==-1)
  {
    printf("[E] Impossible to create %s directory\n", path);
    return(false);
 }
  return(true);
}




bool createOutputDirs(model_parameters *const params)
{
  if(!checkOrCreateDir(params->OUTPUT_PATH))
  {
    return(false);
  }
  for (size_t mw = 0; mw < params->values.mw_n; mw++) {
    for (size_t QG = 0; QG < params->values.QG_n; QG++) {
      for (size_t Tcor = 0; Tcor < params->values.TCor_n; Tcor++) {
        for (size_t Tcor2 = 0; Tcor2 < params->values.TCor2_n; Tcor2++) {
          for (size_t Pcor = 0; Pcor < params->values.PCor_n; Pcor++) {
            for (size_t deltaH = 0; deltaH < params->values.deltaH_n; deltaH++) {
              for (size_t len = 0; len < params->values.len_n; len++) {
                for (size_t flat = 0; flat < params->values.flat_n; flat++) {
                  char path[400] = "";
                  sprintf(path,
                          "%s/"
                          "m_%.3f_Q_%.2f_Pcor_%.0f_Tcor_%.1f_Tcor2_%.1f_dH_%.0f_len_%.0f_flat_%"
                          ".0f_%s_Rho_Snow_%d_Thermal_Ice_%s_Thermal_Firn_%s_Heat_Capacity_%s_"
                          "Rho_Firn_%s_"
                          "Internal_Energy_%s_Scheme_%s",
                          params->OUTPUT_PATH, params->values.mw[mw], params->values.QG[QG] * 1000,
                          params->values.PCor[Pcor], params->values.TCor[Tcor], params->values.TCor2[Tcor2],
                          params->values.deltaH[deltaH], params->values.len[len], params->values.flat[flat], "EDC",
                          params->RHO_SNOW, params->strings[THERMAL_ICE][params->THERMAL_ICE],
                          params->strings[THERMAL_FIRN][params->THERMAL_FIRN],
                          params->strings[HEAT_CAPACITY][params->HEAT_CAPACITY],
                          params->strings[RHO_FIRN][params->RHO_FIRN],
                          params->strings[INTERNAL_ENERGY][params->INTERNAL_ENERGY],
                          params->strings[SCHEME][params->SCHEME]);
                  // Check if the man export directory and the subdirectory are already existing
                if(!checkOrCreateDir(path))
                {
                  return(false);
                }
  }}}}}}}}
  return(true);
}
