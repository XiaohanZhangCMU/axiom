/*
  organizer.h
  by Wei Cai  caiwei@mit.edu
  Last Modified : Wed Jan 21 23:21:38 2009

  FUNCTION  : Organizer package to facilitate MD simulation

  Featuring : class Oragnizer, class AUXFile
*/

#ifndef _ORGANIZER_H
#define _ORGANIZER_H

#include "scparser.h"
#include "filecls.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>

#ifdef _PARALLEL
#include <mpi.h>
#endif

class AUXFile
{
public:
    LFile *f; int type; int count; char fname[1000], zipname[1000], extname[1000];
    virtual int writeentry(void *){return 0;}
    virtual int writeblock(void *){return 0;}
    virtual int readblock(void *){return 0;}
    virtual char *describe(){return 0;}
public:
    enum {ENTRIES,BLOCK,SERIES};
    AUXFile(int t):f(NULL),type(t),count(1){};
    virtual ~AUXFile(){close();}
    int close(int zip=0, bool bg=true);
    int open(const char *name,int format=LFile::O_WriteNew);
    int reopen();
    int write(void *p, int zip=1, bool bg=true);
    int read(void *p);
    static const char *Skip(const char *str, const char *blank)
    {
        ASSERT(str!=NULL);
        ASSERT(blank!=NULL);
        const char *ret=str;
        while(*ret!=0 && strchr(blank, *ret)) ret++;
        return ret;
    }
    static const char *SkipBlanks(const char *str)
        { return Skip(str, " \n\t\r"); }

    void setcount(int n){count=n;}
private:
    static int lastoccur(char *n,char c);        
public:
    static void insertindex(char *ext,char *n,int c);
};


class Organizer : public SCParser
{
public:
    bool overwrite,nolog,renew,diropened;
    int sleepseconds;
    LFile *olog;
public:
    char dirname[1000], logfilename[1000];
public:
    Organizer():overwrite(false),nolog(false),renew(false),diropened(false),
        sleepseconds(600),olog(NULL){init();};
    virtual ~Organizer() {} //{closefiles();}

    void bindvar(const char *vn,void *p,int type)
        {SCParser::bindvar(vn,p,type);}
    int exec(const char *name);
    int parse(FILE *file);
    int parse_line(FILE *file);
    int assignvar(int offset=0);
    void quit();
    
public:
    void setoverwrite(bool b=true){overwrite=b;}
    void setnolog(bool b=true){nolog=b;}
    void getsleep();
    void init(); //system use
    int SysParse(FILE *file); //system use
    char *currentime();
    void printSysInfo();
    void printEndInfo();
    void closefiles(int zip=1,bool bg=true);
    int opendir();

public:
    void welcome();

    //dirname
    inline std::string get_dirname(){ std::string str(dirname); return str; };
    inline void set_dirname(std::string s){ strcpy(dirname, s.c_str()); };
};

#define ah HIC "ASSIGN " NOR
#define eh HIB "EXEC   " NOR

#define INFO_buffer_plain(pre,ref,type) if(shift==0) \
                                            INFO(pre<<varname[curn] \
                                                <<" = "\
                                                <<ref((type *)varptr[curn])); \
                                        else INFO(pre<<varname[curn] \
                                                <<"("<<shift<<") = "\
                                                <<ref(((type *)varptr[curn])+shift))
#define INFO_buffer_offset(pre,ref,type) INFO(pre<<varname[curn] \
                                       <<"("<<offset<<") = "\
                                       <<ref(((type *)varptr[curn])+offset))
#define INFO_buffer(pre,ref,type) if(offset==0) \
                                       INFO_buffer_plain(pre,ref,type);\
                                  else INFO_buffer_offset(pre,ref,type);

#endif // _ORGANIZER_H

