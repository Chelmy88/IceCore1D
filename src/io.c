#include "io.h"

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
    printf("[I] File: %s opened in reading mode\n",fileName);
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

bool readIOFile(model_parameters *params, char* fileName)
{
  FILE *fp;
  if((fp=fopen(fileName, "r"))==NULL)
  {
    printf("[E] Cannot open file: %s\n",fileName);
    return false;
  }

  printf("[I] File: %s opened in reading mode\n",fileName);
  char *line = NULL;
  size_t len = 0;
  while(getline(&line, &len, fp) != -1)
  {
    removeSpaces(line);
    removeComments(line);
    printf("READ: %s \n", line);
  }
  fclose(fp);
  return true;
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
  ptr = strchr(source, '/');
  if (ptr != NULL) {
    *ptr = '\0';
  }
}
