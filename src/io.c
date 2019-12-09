#include "io.h"
#include <dirent.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>
//*************File management functions*************

bool readTable(double* table, const char* const fileName)
{
  // Read a table from the file called "filename" and store it into the double table called "table". Table should be 1D. "filename" can contain directory path.
  FILE *fp;
  int li=0;
  double a=0;
  if((fp=fopen(fileName, "r"))==NULL)
  {
    printf("[E] annot open file: %s\n",fileName);
    return false;
  }
  else
  {
    while(fscanf(fp,"%lf",&a)==1)
    {
      table[li]=a;
      li++;
    }
    if(fclose(fp)==0)
    {
      printf("[I] File imported successfully (%d data) and closed: %s \n\n",
             li,fileName);
      return true;
    }
    else
    {
      printf("[E] Not able to close: %s \n\n",fileName);
      return false;
    }
  }
}


void saveTable(double *table,const char* const name, const char* const path,
               const int tabSize)
{
    // Read the 1D double table called "table" and store it a file called "filename".  Table should be tab delimited. "filename" can contain directory path. The default path is a folder called "export". The table size should be passed as parameter.
  FILE *fp;
  int li=0;
  char full_path[180]="";
  sprintf(full_path,"%s/%s",path,name);
  if((fp=fopen(full_path, "w+"))==NULL)
  {
    printf("Cannot open file: %s\n",full_path);
  }
  else
  {
    printf("File oppened:%s\n...writing...\n",full_path);
    for (li=0; li<tabSize; li++)
    {
      fprintf(fp,"%f \n",table[li]);
    }
    fclose(fp);
    printf("File closed: %s \n\n",full_path);
  }
}

void save2DTable(double **table, const char* const name, const char* const path,
                 const int nRow, const int nCol, const int skipR,
                 const int skipC, const int startC)
{
  // Read the 2D double table called "table" and store it a tab delimited
  // file called "filename". The default path is a folder called "export".
  // The table dimensions should be passed as parameter.
  FILE *fp;
  int li=0;
  int co=0;
  char full_path[180]="";
  sprintf(full_path,"%s/%s",path,name);
  if((fp=fopen(full_path, "w+"))==NULL)
  {
    printf("[E] Cannot open file: %s\n",full_path);
  }
  else
  {
    printf("File oppened:%s\n...writing...\n",full_path);
    for (li=0; li<nRow; li+=skipR)
    {
      for (co=startC; co<nCol; co+=skipC)
      {
        if(co!=startC)
        {
          fprintf(fp,"\t");
        }
        if(table[li][co]>0)
        {
          fprintf(fp,"%f",table[li][co]);
        }
        else
        {
          fprintf(fp,"NaN");
        }
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
    printf("File closed: %s \n\n",full_path);
  }
}

//save2DTable_top(density,fileName,path,iceThickness,100,T);

void save2DTable_top(double **table, const char *const name, const char* const path,
     const double* const thickness, const int nRow, const int nCol)
{
  // Read the 2D double table called "table" and store it a tab delimited
  // file called "filename". The default path is a folder called "export".
  // The table dimensions should be passed as parameter.
  FILE *fp;
  int li=0;
  int co=0;
  int startC=0;
  char full_path[180]="";
  sprintf(full_path,"%s/%s",path,name);
  if((fp=fopen(full_path, "w+"))==NULL)
  {
      printf("Cannot open file: %s\n",full_path);
  }
  else
  {
    printf("File oppened:%s\n...writing...\n",full_path);
    for (li=0; li<nRow; li+=1)
    {
      for (co=startC; co<nCol; co+=1)
      {
        int real_li=0;
        if(co!=startC)
        {
          fprintf(fp,"\t");
          real_li = (int)(thickness[co-1])-li;
        }
        else{
          real_li = (int)(thickness[co])-li;
        }
        if(table[real_li][co]>0)
        {
          fprintf(fp,"%f",table[real_li][co]);
        }
        else
        {
          fprintf(fp,"NaN");
        }
      }
    fprintf(fp,"\n");
    }
    fclose(fp);
    printf("File closed: %s \n\n",full_path);
  }
}

bool readInitFile(model_parameters* const params, const char* const fileName)
{
  FILE *fp;
  if((fp=fopen(fileName, "r"))==NULL)
  {
    printf("[E] Cannot open file: %s\n",fileName);
    return false;
  }
  char *line = NULL;
  size_t len = 0;
  bool state=true;
  while(getline(&line, &len, fp) != -1)
  {
    removeComments(line);
    if(line[0] == '\0')
    {
      continue;
    }
    char *arg=NULL;
    char *val=NULL;
    if(!parseLine(line,&arg,&val))
    {
      state=false;
    }
    else if(!setParameter(params,arg,val))
    {
      state=false;
    }
  }
  fclose(fp);
  if(state)
  {
    if(!checkParameter(params))
    {
      state=false;
    }
  }
  if(state)
  {
    params->values.tot=params->values.mw_n *params->values.QG_n *
                       params->values.TCor_n *params->values.TCor2_n *
                       params->values.PCor_n *params->values.deltaH_n *
                       params->values.len_n *params->values.flat_n;

    printf("[I] Ini file %s read successfully\n",fileName);
    printModelParameters(params);
    return true;
  }
  else{
    return false;
  }
}

void removeSpaces(char* source)
{
  char* i = source;
  char* j = source;
  do {
    *i = *j;
    if(*i != ' ') i++;
  }
  while(*j++ != 0);
}
void removeComments(char* const source)
{
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

bool parseLine(char* line,char** arg, char** val)
{
  char * pch;
  pch = strtok (line,":");
  if(pch != NULL)
  {
    *arg=pch;
    removeSpaces(*arg);
    upper_string(*arg);
    pch = strtok (NULL, ":");
    if(pch != NULL)
    {
      *val=pch;
      pch = strtok (NULL, ":");
      if(pch != NULL)
      {
        printf("[E] Two ':' delimiter found in argument line %s\n",line);
        return false;
      }
    }
    else
    {
      printf("[E] No ':' delimiter found in argument line %s\n",line);
      return false;
    }
  }
  else
  {
    printf("[E] No ':' delimiter found in argument line %s\n",line);
    return false;

  }
  return true;
}

//Rather long function that test for the various allowed keys
//and strores the corresponding values
bool setParameter(model_parameters* const params, const char* const arg,
                  char* const val)
{
  if(strcmp(arg,"Z")==0)
  {
    removeSpaces(val);
    if(atoi(val)>0)
    {
      params->Z1=atoi(val);
    }
    else
    {
      printf("[E] Wrong parameter %s for Z. Must be positive integer.\n",val);
      return false;
    }
  }

  else if(strcmp(arg,"T")==0)
  {
    removeSpaces(val);
    if(atoi(val)>0)
    {
      params->T1=atoi(val);
    }
    else
    {
      printf("[E] Wrong parameter %s for Z. Must be positive integer.\n",val);
      return false;
    }
  }

  else if(strcmp(arg,"S")==0)
  {
    removeSpaces(val);
    if(atoi(val)>0)
    {
      params->S1=atoi(val);
    }
    else
    {
      printf("[E] Wrong parameter %s for Z. Must be positive integer\n",val);
      return false;
    }
  }

  else if(strcmp(arg,"OUTPUT_PATH")==0)
  {
    removeSpaces(val);
    DIR* dir = opendir(val);
    if (dir)
    {
        params->OUTPUT_PATH = (char*)malloc(sizeof(char) * (strlen(val)+1));
        strcpy(params->OUTPUT_PATH,val);
        closedir(dir);
    }
    else if (ENOENT == errno)
    {
      printf("[E] Wrong parameter %s for OUTPUT_PATH, directory not exist\n",val);
      return false;
    }
    else
    {
      printf("[E] Failed to open OUTPUT_PATH %s\n",val);
      return false;
    }
  }

  else if(strcmp(arg,"TEMPERATURE_FILE")==0)
  {
    removeSpaces(val);
    if( access( val, R_OK ) != -1 ) {
      params->TEMPERATURE_FILE = (char*)malloc(sizeof(char) * (strlen(val)+1));
      strcpy(params->TEMPERATURE_FILE,val);
    }
    else
    {
      printf("[E] Wrong parameter %s for TEMPERATURE_FILE, file not exist or not readable\n",val);
      return false;
    }
  }

  else if(strcmp(arg,"ACCUMULATION_FILE")==0)
  {
    removeSpaces(val);
    if( access( val, R_OK ) != -1 ) {
      params->ACCUMULATION_FILE = (char*)malloc(sizeof(char) * (strlen(val)+1));
      strcpy(params->ACCUMULATION_FILE,val);
    }
    else
    {
      printf("[E] Wrong parameter %s for ACCUMULATION_FILE, file not exist or not readable\n",val);
      return false;
    }
  }

  else if(strcmp(arg,"ICE_THICKNESS_FILE")==0)
  {
    removeSpaces(val);
    if( access( val, R_OK ) != -1 ) {
      params->ICE_THICKNESS_FILE = (char*)malloc(sizeof(char) * (strlen(val)+1));
      strcpy(params->ICE_THICKNESS_FILE,val);
    }
    else
    {
      printf("[E] Wrong parameter %s for ICE_THICKNESS_FILE, file not exist or not readable\n",val);
      return false;
    }
  }

  else if(strcmp(arg,"AGE_FILE")==0)
  {
    removeSpaces(val);
    if( access( val, R_OK ) != -1 ) {
      params->AGE_FILE = (char*)malloc(sizeof(char) * (strlen(val)+1));
      strcpy(params->AGE_FILE,val);
    }
    else
    {
      printf("[E] Wrong parameter %s for AGE_FILE, file not exist or not readable\n",val);
      return false;
    }
  }

  else if(strcmp(arg,"BOREHOLE_TEMPERATURE_FILE")==0)
  {
    removeSpaces(val);
    if( access( val, R_OK ) != -1 ) {
      params->BOREHOLE_TEMPERATURE_FILE = (char*)malloc(sizeof(char) * (strlen(val)+1));
      strcpy(params->BOREHOLE_TEMPERATURE_FILE,val);
    }
    else
    {
      printf("[E] Wrong parameter %s for BOREHOLE_TEMPERATURE_FILE, file not exist or not readable\n",val);
      return false;
    }
  }

  else if(strcmp(arg,"SAVE_TYPE")==0)
  {
    removeSpaces(val);
    if (strcmp(val,"MATRIX")==0)
    {
      params->SAVE_TYPE1=ST_MATRIX;
    }
    else if (strcmp(val,"VECTOR")==-0)
    {
      params->SAVE_TYPE1=ST_VECTOR;
    }
    else
    {
      printf("[E] Unknown value %s for SAVE_TYPE\n",val);
      return false;
    }
  }

  else if(strcmp(arg,"SCHEME")==0)
  {
    removeSpaces(val);
    if (strcmp(val,"CN")==0)
    {
      params->SCHEME1=SC_CN;
    }
    else if (strcmp(val,"EXPL")==0)
    {
      params->SCHEME1=SC_EXPL;
    }
    else
    {
      printf("[E] Unknown value %s for SAVE_TYPE\n",val);
      return false;
    }
  }

  else if(strcmp(arg,"RHOSNOW")==0)
  {
    removeSpaces(val);
    if(atoi(val)>0)
    {
      params->rhoSnow1=atoi(val);
    }
    else
    {
      printf("[E] Wrong parameter %s for RHOSNOW. Must be positive integer\n",val);
      return false;
    }
  }

  else if(strcmp(arg,"THERMAL")==0)
  {
    removeSpaces(val);
    if (strcmp(val,"CP")==0)
    {
      params->THERMAL1=TH_CP;
    }
    else if (strcmp(val,"GO")==0)
    {
      params->THERMAL1=TH_GO;
    }
    else
    {
      printf("[E] Unknown value %s for THERMAL\n",val);
      return false;
    }
  }

  else if(strcmp(arg,"FIRN")==0)
  {
    removeSpaces(val);
    if (strcmp(val,"SC")==0)
    {
      params->FIRN1=FI_SC;
    }
    else if (strcmp(val,"CP")==0)
    {
      params->FIRN1=FI_CP;
    }
    else if (strcmp(val,"FI")==0)
    {
      params->FIRN1=FI_FI;
    }
    else
    {
      printf("[E] Unknown value %s for FIRN1\n",val);
      return false;
    }
  }

  else if(strcmp(arg,"RHO")==0)
  {
    removeSpaces(val);
    if (strcmp(val,"FIRN")==0)
    {
      params->RHO1=RHO_FIRN;
    }
    else if (strcmp(val,"CONST")==0)
    {
      params->RHO1=RHO_CONST;
    }
    else
    {
      printf("[E] Unknown value %s for RHO\n",val);
      return false;
    }
  }

  else if(strcmp(arg,"VERTICAL_PROFILE")==0)
  {
    removeSpaces(val);
    if (strcmp(val,"FI")==0)
    {
      params->VERTICAL_PROFILE1=VP_FI;
    }
    else if (strcmp(val,"PA")==0)
    {
      params->VERTICAL_PROFILE1=VP_PA;
    }
    else
    {
      printf("[E] Unknown value %s for VERTICAL1\n",val);
      return false;
    }
  }

  else if(strcmp(arg,"INTERNAL_ENERGY")==0)
  {
    removeSpaces(val);
    if (strcmp(val,"ON")==0)
    {
      params->INTERNAL_ENERGY1=IE_ON;
    }
    else if (strcmp(val,"OFF")==0)
    {
      params->INTERNAL_ENERGY1=IE_OFF;
    }
    else
    {
      printf("[E] Unknown value %s for INTERNAL_ENERGY\n",val);
      return false;
    }
  }
  else if(strcmp(arg,"MELTING")==0)
  {
    removeSpaces(val);
    if (strcmp(val,"FREE_MELT")==0)
    {
      params->MELTING1=ME_FREE_MELT;
    }
    else if (strcmp(val,"NO_ICE")==0)
    {
      params->MELTING1=ME_FREEZING_NO_ICE;
    }
    else if (strcmp(val,"FREEZING")==0)
    {
      params->MELTING1=ME_FREEZING;
    }
    else
    {
      printf("[E] Unknown value %s for MELTING\n",val);
      return false;
    }
  }

  else if(strcmp(arg,"MW")==0)
  {
    char *token;
    // make a copy for the first pass
    char val_copy[40];
    strcpy(val_copy,val);
    // get the first token
    token = strtok(val_copy, " ");
    // walk through other tokens and store the size
    // also validate the input
    size_t n=0;
    while( token != NULL ) {
      if(atof(token)>0 || strcmp(token,"0")==0 || strcmp(token,"0.0")==0)
      {
       ++n;
      }
      else{
       printf("[E] Unknown value %s for MW\n",token);
       return false;
      }
      token = strtok(NULL, " ");
    }
    // Now alloc the array
    params->values.mw = (double *) malloc(sizeof(double) * n);
    // and store the data in it
    token = strtok(val, " ");
    n=0;
    while( token != NULL ) {
      params->values.mw[n]=atof(token);
      token = strtok(NULL, " ");
      ++n;
    }
    params->values.mw_n=n;
  }

  else if(strcmp(arg,"QG")==0)
  {
    char *token;
    // make a copy for the first pass
    char val_copy[40];
    strcpy(val_copy,val);
    // get the first token
    token = strtok(val_copy, " ");
    // walk through other tokens and store the size
    // also validate the input
    size_t n=0;
    while( token != NULL ) {
      if(atof(token)!=0 || strcmp(token,"0")==0 || strcmp(token,"0.0")==0)
      {
       ++n;
      }
      else{
       printf("[E] Unknown value %s for QG\n",token);
       return false;
      }
      token = strtok(NULL, " ");
    }
    // Now alloc the array
    params->values.QG = (double *) malloc(sizeof(double) * n);
    // and store the data in it
    token = strtok(val, " ");
    n=0;
    while( token != NULL ) {
      params->values.QG[n]=atof(token);
      token = strtok(NULL, " ");
      ++n;
    }
    params->values.QG_n=n;
  }

  else if(strcmp(arg,"TCOR")==0)
  {
    char *token;
    // make a copy for the first pass
    char val_copy[40];
    strcpy(val_copy,val);
    // get the first token
    token = strtok(val_copy, " ");
    // walk through other tokens and store the size
    // also validate the input
    size_t n=0;
    while( token != NULL ) {
      if(atof(token)>0 || strcmp(token,"0")==0 || strcmp(token,"0.0")==0)
      {
       ++n;
      }
      else{
       printf("[E] Unknown value %s for TCOR\n",token);
       return false;
      }
      token = strtok(NULL, " ");
    }
    // Now alloc the array
    params->values.TCor = (double *) malloc(sizeof(double) * n);
    // and store the data in it
    token = strtok(val, " ");
    n=0;
    while( token != NULL ) {
      params->values.TCor[n]=atof(token);
      token = strtok(NULL, " ");
      ++n;
    }
    params->values.TCor_n=n;
  }

  else if(strcmp(arg,"TCOR2")==0)
  {
    char *token;
    // make a copy for the first pass
    char val_copy[40];
    strcpy(val_copy,val);
    // get the first token
    token = strtok(val_copy, " ");
    // walk through other tokens and store the size
    // also validate the input
    size_t n=0;
    while( token != NULL ) {
      if(atof(token)>0 || strcmp(token,"0")==0 || strcmp(token,"0.0")==0)
      {
       ++n;
      }
      else{
       printf("[E] Unknown value %s for TCOR2\n",token);
       return false;
      }
      token = strtok(NULL, " ");
    }
    // Now alloc the array
    params->values.TCor2 = (double *) malloc(sizeof(double) * n);
    // and store the data in it
    token = strtok(val, " ");
    n=0;
    while( token != NULL ) {
      params->values.TCor2[n]=atof(token);
      token = strtok(NULL, " ");
      ++n;
    }
    params->values.TCor2_n=n;
  }

  else if(strcmp(arg,"PCOR")==0)
  {
    char *token;
    // make a copy for the first pass
    char val_copy[40];
    strcpy(val_copy,val);
    // get the first token
    token = strtok(val_copy, " ");
    // walk through other tokens and store the size
    // also validate the input
    size_t n=0;
    while( token != NULL ) {
      if(atof(token)>0 || strcmp(token,"0")==0 || strcmp(token,"0.0")==0)
      {
       ++n;
      }
      else{
       printf("[E] Unknown value %s for PCOR\n",token);
       return false;
      }
      token = strtok(NULL, " ");
    }
    // Now alloc the array
    params->values.PCor = (double *) malloc(sizeof(double) * n);
    // and store the data in it
    token = strtok(val, " ");
    n=0;
    while( token != NULL ) {
      params->values.PCor[n]=atof(token);
      token = strtok(NULL, " ");
      ++n;
    }
    params->values.PCor_n=n;
  }

  else if(strcmp(arg,"DELTAH")==0)
  {
    char *token;
    // make a copy for the first pass
    char val_copy[40];
    strcpy(val_copy,val);
    // get the first token
    token = strtok(val_copy, " ");
    // walk through other tokens and store the size
    // also validate the input
    size_t n=0;
    while( token != NULL ) {
      if(atof(token)>0 || strcmp(token,"0")==0 || strcmp(token,"0.0")==0)
      {
       ++n;
      }
      else{
       printf("[E] Unknown value %s for DELTAH\n",token);
       return false;
      }
      token = strtok(NULL, " ");
    }
    // Now alloc the array
    params->values.deltaH = (double *) malloc(sizeof(double) * n);
    // and store the data in it
    token = strtok(val, " ");
    n=0;
    while( token != NULL ) {
      params->values.deltaH[n]=atof(token);
      token = strtok(NULL, " ");
      ++n;
    }
    params->values.deltaH_n=n;
  }
  else if(strcmp(arg,"LEN")==0)
  {
    char *token;
    // make a copy for the first pass
    char val_copy[40];
    strcpy(val_copy,val);
    // get the first token
    token = strtok(val_copy, " ");
    // walk through other tokens and store the size
    // also validate the input
    size_t n=0;
    while( token != NULL ) {
      if(atof(token)>0 || strcmp(token,"0")==0 || strcmp(token,"0.0")==0)
      {
       ++n;
      }
      else{
       printf("[E] Unknown value %s for LEN\n",token);
       return false;
      }
      token = strtok(NULL, " ");
    }
    // Now alloc the array
    params->values.len = (double *) malloc(sizeof(double) * n);
    // and store the data in it
    token = strtok(val, " ");
    n=0;
    while( token != NULL ) {
      params->values.len[n]=atof(token);
      token = strtok(NULL, " ");
      ++n;
    }
    params->values.len_n=n;
  }

  else if(strcmp(arg,"FLAT")==0)
  {
    char *token;
    // make a copy for the first pass
    char val_copy[40];
    strcpy(val_copy,val);
    // get the first token
    token = strtok(val_copy, " ");
    // walk through other tokens and store the size
    // also validate the input
    size_t n=0;
    while( token != NULL ) {
      if(atof(token)>0 || strcmp(token,"0")==0 || strcmp(token,"0.0")==0)
      {
       ++n;
      }
      else{
       printf("[E] Unknown value %s for FLAT\n",token);
       return false;
      }
      token = strtok(NULL, " ");
    }
    // Now alloc the array
    params->values.flat = (double *) malloc(sizeof(double) * n);
    // and store the data in it
    token = strtok(val, " ");
    n=0;
    while( token != NULL ) {
      params->values.flat[n]=atof(token);
      token = strtok(NULL, " ");
      ++n;
    }
    params->values.flat_n=n;
  }
  else
  {
    printf("[E] Unkown parameter %s find in ini file\n",arg);
    return false;
  }
  return true;
}

bool checkParameter(model_parameters* const params){

  bool state=true;
  if(params->Z1==-1)
  {
    printf("[E] Unset parameter Z. This parameter is mandatory\n");
    state=false;
  }
  if(params->T1==-1)
  {
    printf("[E] Unset parameter T. This parameter is mandatory\n");
    state=false;
  }
  if(params->S1==-1)
  {
    params->S1=1500;
    printf("[W] Unset parameter S, the default value of 1500 is used\n");
  }
  if(params->rhoSnow1==-1)
  {
    params->rhoSnow1=350;
    printf("[W] Unset parameter RHOSNOW, the default value of 350 is used\n");
  }
  if(!params->OUTPUT_PATH)
  {
    printf("[E] Unset parameter OUTPUT_PATH. This parameter is mandatory\n");
    state=false;
  }
  if(!params->TEMPERATURE_FILE)
  {
    printf("[E] Unset parameter TEMPERATURE_FILE. This parameter is mandatory\n");
    state=false;
  }
  if(!params->ACCUMULATION_FILE)
  {
    printf("[E] Unset parameter ACCUMULATION_FILE. This parameter is mandatory\n");
    state=false;
  }
  if(!params->ICE_THICKNESS_FILE)
  {
    printf("[E] Unset parameter ICE_THICKNESS_FILE. This parameter is mandatory\n");
    state=false;
  }
  if(!params->AGE_FILE)
  {
    printf("[E] Unset parameter AGE_FILE. This parameter is mandatory\n");
    state=false;
  }
  if(!params->BOREHOLE_TEMPERATURE_FILE)
  {
    printf("[E] Unset parameter BOREHOLE_TEMPERATURE_FILE. This parameter is mandatory\n");
    state=false;
  }
  if(params->SAVE_TYPE1==ST_UNSET)
  {
    params->SAVE_TYPE1=ST_VECTOR;
    printf("[W] Unset parameter SAVE_TYPE, the default value of VECTOR is used\n");
  }
  if(params->SCHEME1==SC_UNSET)
  {
    params->SCHEME1=SC_CN;
    printf("[W] Unset parameter SCHEME, the default value of CN is used\n");
  }
  if(params->THERMAL1==TH_UNSET)
  {
    params->THERMAL1=TH_CP;
    printf("[W] Unset parameter THERMAL, the default value of CP is used\n");
  }
  if(params->FIRN1==FI_UNSET)
  {
    params->FIRN1=FI_SC;
    printf("[W] Unset parameter FIRN, the default value of SC is used\n");
  }
  if(params->RHO1==RHO_UNSET)
  {
    params->RHO1=RHO_FIRN;
    printf("[W] Unset parameter RHO, the default value of FIRN is used\n");
  }
  if(params->VERTICAL_PROFILE1==VP_UNSET)
  {
    params->VERTICAL_PROFILE1=VP_FI;
    printf("[W] Unset parameter VERTICAL_PROFILE, the default value of FI is used\n");
  }
  if(params->INTERNAL_ENERGY1==IE_UNSET)
  {
    params->INTERNAL_ENERGY1=IE_OFF;
    printf("[W] Unset parameter INTERNAL_ENERGY, the default value of OFF is used\n");
  }
  if(params->MELTING1==ME_UNSET)
  {
    params->MELTING1=ME_FREE_MELT;
    printf("[W] Unset parameter MELTING, the default value of FREE_MELT is used\n");
  }
  if(!params->values.mw)
  {
    printf("[E] Unset parameter MW. This parameter is mandatory\n");
    state=false;
  }
  if(!params->values.QG)
  {
    printf("[E] Unset parameter QG. This parameter is mandatory\n");
    state=false;
  }
  if(!params->values.TCor)
  {
    params->values.TCor=(double *) malloc(sizeof(double));
    params->values.TCor[0]=0.;
    params->values.TCor_n=1;
    printf("[W] Unset parameter TCOR, the default value of 0 is used\n");
  }
  if(!params->values.TCor2)
  {
    params->values.TCor2=(double *) malloc(sizeof(double));
    params->values.TCor2[0]=0.;
    params->values.TCor2_n=1;
    printf("[W] Unset parameter TCOR2, the default value of 0 is used\n");
  }
  if(!params->values.PCor)
  {
    params->values.PCor=(double *) malloc(sizeof(double));
    params->values.PCor[0]=0.;
    params->values.PCor_n=1;
    printf("[W] Unset parameter PCOR, the default value of 0 is used\n");
  }
  if(!params->values.deltaH)
  {
    params->values.deltaH=(double *) malloc(sizeof(double));
    params->values.deltaH[0]=0.;
    params->values.deltaH_n=1;
    printf("[W] Unset parameter DELTAH, the default value of 0 is used\n");
  }
  if(!params->values.len)
  {
    params->values.len=(double *) malloc(sizeof(double));
    params->values.len[0]=0.;
    params->values.len_n=1;
    printf("[W] Unset parameter LEN, the default value of 0 is used\n");
  }
  if(!params->values.flat)
  {
    params->values.flat=(double *) malloc(sizeof(double));
    params->values.flat[0]=0.;
    params->values.flat_n=1;
    printf("[W] Unset parameter FLAT, the default value of 0 is used\n");
  }
  return(state);
}


void upper_string(char* s) {
   size_t c = 0;
   while (s[c] != '\0') {
      if (s[c] >= 'a' && s[c] <= 'z') {
         s[c] = s[c] - 32;
      }
      ++c;
   }
}
