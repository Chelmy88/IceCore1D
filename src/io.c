#include "io.h"
#include <dirent.h>
#include <errno.h>
//*************File management functions*************

bool readTable(double* table,char* fileName)
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
      printf("[I] File imported successfully (%d data) and closed: %s \n\n",li,fileName);
      return true;
    }
    else
    {
      printf("[E] Not able to close: %s \n\n",fileName);
      return false;
    }
  }
}


void saveTable(double *table,char *name, char* path,int tabSize)
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

void save2DTable(double **table,char *name, char* path,int nRow,int nCol, int skipR, int skipC,int startC)
{
  // Read the 2D double table called "table" and store it a tab delimited file called "filename". The default path is a folder called "export". The table dimensions should be passed as parameter.
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

void save2DTable_top(double **table,char *name, char* path,double* thickness, int nRow,int nCol)
{
    // Read the 2D double table called "table" and store it a tab delimited file called "filename". The default path is a folder called "export". The table dimensions should be passed as parameter.
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

bool readInitFile(model_parameters *params, char* fileName)
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
    removeSpaces(line);
    removeComments(line);
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
void removeComments(char* source)
{
  char *ptr;
  ptr = strchr(source, '//');
  if (ptr != NULL) {
    *ptr = '\0';
  }
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
    pch = strtok (NULL, ":");
    if(pch != NULL)
    {
      *val=pch;
      pch = strtok (NULL, ":");
      if(pch != NULL)
      {
        printf("[E] Two ':' delimiter found in argument line %s\n",line);
        return false;      }
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

bool setParameter(model_parameters *params, char* arg, char* val)
{
  if(strcmp(arg,"Z")==0)
  {
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
  else if(strcmp(arg,"SAVE_TYPE")==0)
  {
    if (strcmp(val,"MATRIX")==0)
    {
      params->SAVE_TYPE1=ST_MATRIX;
    }
    else if (strcmp(val,"VECTOR")==0)
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
  else if(strcmp(arg,"rhoSnow")==0)
  {
    if(atoi(val)>0)
    {
      params->rhoSnow1=atoi(val);
    }
    else
    {
      printf("[E] Wrong parameter %s for rhoSnow. Must be positive integer\n",val);
      return false;
    }
  }
  else if(strcmp(arg,"THERMAL")==0)
  {
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
  else if(strcmp(arg,"VERTICAL")==0)
  {
    if (strcmp(val,"FI")==0)
    {
      params->VERTICAL1=VP_FI;
    }
    else if (strcmp(val,"PA")==0)
    {
      params->VERTICAL1=VP_PA;
    }
    else
    {
      printf("[E] Unknown value %s for VERTICAL1\n",val);
      return false;
    }
  }
  else if(strcmp(arg,"INTERNAL_ENERGY")==0)
  {
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
  else
  {
    printf("[E] Unkown parameter %s find in ini file\n",arg);
    return false;
  }

  return true;
}
