#include "io.h"

//*************File management functions*************

void readTable(double* table,char* fileName)
{
    // Read a table from the file called "filename" and store it into the double table called "table". Table should be 1D. "filename" can contain directory path.
    FILE *fp;
    int li=0;
    double a=0;
    if((fp=fopen(fileName, "r"))==NULL)
    {
        printf("Cannot open file: %s\n",fileName);
        exit(0);
    }
    else
    {
        printf("File: %s opened in reading mode\n",fileName);
        while(fscanf(fp,"%lf",&a)==1)
        {
            table[li]=a;
            li++;
        }
        if(fclose(fp)==0)
        {
            printf("File imported successfully (%d data) and closed: %s \n\n",li,fileName);
        }
        else
        {
            printf("Not able to close: %s \n\n",fileName);
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
        printf("Cannot open file: %s\n",full_path);
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
