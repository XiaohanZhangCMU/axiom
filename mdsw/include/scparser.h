/*
  scparser.h
  by Wei Cai  caiwei@mit.edu
  Last Modified : Thu Dec  7 01:10:06 2006

  FUNCTION  :  Input scripting parser module for MDFrame
*/

#ifndef _SCPARSER_H
#define _SCPARSER_H

#include <stdio.h>
#include <string.h>

enum {INT=10,LONG=20,
      DOUBLE=30,STRING=100};
class SCParser
{
#define MAXTAG_BUF    10000
#define MAXVAR        2000
#define MAXNAMELENGTH 500
//#define MAXTAG_BUF    1000
//#define MAXVAR        1000
//#define MAXNAMELENGTH 100

#define protect(statement) do{statement}while(0)
#define bindcommand_real(a,b,c) protect(if(strcmp(a,b)==0){c;return 0;})
#define read_buffer(fmt,type) sscanf(buffer,fmt,\
                              ((type *)(varptr[curn])+shift+offset))
#define output_buffer(port,ref,type) port(varname[curn]<<" = "\
                              <<ref((type *)varptr[curn]))


#define bindcommand(a,b,c) bindcommand_real(a,b,c)

public:
    char buffer[MAXTAG_BUF], cchar;
    char varname[MAXVAR][MAXNAMELENGTH];
    int vartype[MAXVAR];
    void *varptr[MAXVAR];
    int varn,curn,shift;
    bool willabort;
    int ncpu, shmsize;
    
    SCParser():cchar(-1),varn(0),willabort(false),ncpu(1),shmsize(0){};
    virtual ~SCParser(){};
    virtual int exec(const char *name);
    void abortparser();
    virtual int assignvar(int offset=0);
    void dumpbuffer();
    bool bufferis(const char *s);
    bool bufferbeginswith(char c);
    bool bufferendswith(char c);
    bool buffercontains(char c);
    
    void bindvar(const char *vn,void *p,int type);
    static FILE *getfilehandler(int argc, char *argv[]);

    void init(); //system use
    int readnextstring(FILE *file); //into buffer
    int identify(const char *name);
    int parse(FILE *file);
    int parse_buffer(FILE *file);
};


#endif // _SCPARSER_H

