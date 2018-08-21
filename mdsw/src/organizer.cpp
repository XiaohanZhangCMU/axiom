/*
  organizer.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Wed Apr 14 11:33:20 2010

  FUNCTION  : Organizer package to facilitate MD simulation

  Featuring : class Oragnizer, class AUXFile
*/

/* 5/24/2001
   bunzip2 before read
   bzip2   after write */

#include "organizer.h"

int AUXFile::close(int zip, bool bg)
{
    LFile::Close(f);

    if(type==SERIES) DUMP("AUXFile::close  "<<extname);
    else DUMP("AUXFile::close  "<<fname);

    if(zip==1)
    {
        if(type==ENTRIES)
            if(strlen(fname)>0) LFile::GZipFile(fname,bg);
    }
    else if(zip==2)
    {
        if(type==ENTRIES)
            if(strlen(fname)>0) LFile::BZipFile(fname,bg);
    }
    return 0;
}

int AUXFile::open(const char *name,int format)
{
    LFile::SubHomeDir(name,fname);
    
    if(type==SERIES)
    {
        insertindex(extname,fname,count);
        INFO(HIM"FILEOPEN "NOR<<extname<<" "<<describe());
    }
    else INFO(HIM"FILEOPEN "NOR<<fname<<" "<<describe());
    if(f!=NULL)
    {
        ERROR("AUXFile open(): filehandler not NULL");
        return -1;
    }
    else
    {
        if(type==SERIES)
        {
            f=LFile::Open(extname,format);
        }
        else
            f=LFile::Open(fname,format);
        if(f==NULL)
        {
            if(type==BLOCK)
            {
                sprintf(zipname,"%s.bz2",name);
                LFile::SubHomeDir(zipname,fname);
                
                f=LFile::Open(fname,format);
                if(f!=NULL) return 0;
                
                sprintf(zipname,"%s.gz",name);
                LFile::SubHomeDir(zipname,fname);

                f=LFile::Open(fname,format);
                if(f!=NULL) return 0;
            }
            ERROR("AUXFile open(): file open failure");
            return -1;
        }
        return 0;
    }
}

int AUXFile::reopen()
{
    LFile::Close(f);
    INFO(HIM"REOPEN "NOR<<extname<<" "<<describe());
    f=LFile::Open(extname,LFile::O_WriteNew);
    return (f==NULL)?-1:0;
}

int AUXFile::write(void *p, int zip, bool bg)
{/* zip = 0: no zip
          1: gzip
          2: bz2
    bg = 1: zip process in background
 */
    if((type==BLOCK)||(type==SERIES))
    {
        int r;
        if(type==BLOCK) INFO(HIM"WRITEFILE "NOR<<fname);
        else INFO(HIM"WRITEFILE "NOR<<extname);
        if(f==NULL)
        {
            if((type==BLOCK)||((type==SERIES)&&(count==0)))
            {
                ERROR("AUXFile write(): file handler NULL");
                return -1;
            }
            else
                reopen();
        }
        r=writeblock(p); if(r!=0) return -1;
        LFile::Close(f);

        /* Zip file */
        if(zip==1)
        {
            if(type==SERIES) LFile::GZipFile(extname,bg);
            else LFile::GZipFile(fname,bg);
        }
        else if(zip==2)
        {
            if(type==SERIES) LFile::BZipFile(extname,bg);
            else LFile::BZipFile(fname,bg);
        }
        
        
        if(type==SERIES)
        {
            count++;
            insertindex(extname,fname,count);
        }
        return r;
    }
    else
    {
        DUMP(HIM"(entry) WRITEFILE "NOR<<fname);
        if(f==NULL){
            ERROR("AUXFile write(): file handler NULL");
            return -1;
        }
        return writeentry(p);
    }
}

int AUXFile::read(void *p)
{
    INFO(HIM"READFILE "NOR<<fname);
    if(f==NULL){
        FATAL("AUXFile read(): file handler NULL");
        return -1;
    }
    if(type!=BLOCK)
    {
        ERROR("AUXFile read(): read file not BLOCK type");
        return -1;
    }
    else
        readblock(p);
    LFile::Close(f);
    return 0;
}


int AUXFile::lastoccur(char *n,char c)
{
    int len,i; int ret;
    len=strlen(n);
    ret=-1;
    for(i=len-1;i>=0;i--)
        if(n[i]==c) {ret=i;break;}
    return ret;
}


void AUXFile::insertindex(char *ext,char *n,int c)
{
    int i; char tmp[100];
    sprintf(tmp,"%04d",c);
    strcpy(ext,n);
    i=lastoccur(n,'.');
    if(i<0) strcat(ext,tmp);
    else
    {
        strcat(tmp,ext+i); ext[i]=0;
        strcat(ext,tmp);
    }
}

int Organizer::exec(const char *name)
{
    INFO(eh<<name);
    
    if(SCParser::exec(name)==0) return 0;
    bindcommand(name,"setoverwrite",setoverwrite());
    //bindcommand(name,"setnolog",setnolog());
    bindcommand(name,"quit",quit());
    bindcommand(name,"sleep",getsleep());    
    return -1;
}

int Organizer::parse(FILE *file)
{
    int ret;
    if (file==NULL)
    {
        ERROR("parse: file==NULL!");
        return -1;
    }
    SysParse(file); //acquire dirname
    //opendir();
    //printSysInfo();
    ret=SCParser::parse(file);
    printEndInfo();
    closefiles();
    return ret;
}

int Organizer::parse_line(FILE *file)
{
    int ret;
    if (file==NULL)
    {
        ERROR("parse: file==NULL!");
        return -1;
    }
    ret=SCParser::parse(file);
//    if((!diropened)&&(strlen(dirname)>0))
//    {
//        INFO_Printf("I am in parse_line again\n");
//        opendir();
//        printSysInfo();
//        diropened=true;
//    }
    return ret;
}

#if 0
int Organizer::parse_(FILE *file)
{
    LFile *f;
    f = LFile::Open("A.log",LFile::O_Read);
    LFile::Close(f);
    return 0;
}
#endif

int Organizer::assignvar(int offset)
{
    int s;
    switch(vartype[curn])
    {
    case(INT): s=read_buffer("%d",int); break;
    case(LONG): s=read_buffer("%ld",long); break;
    case(DOUBLE): s=read_buffer("%lf",double); break;
    /* case(STRING): s=read_buffer("%s",char); break; */
    case(STRING): s=1;strcpy((char *)varptr[curn],buffer);break;
     default: FATAL("unknown vartype ("<<vartype[curn]<<")");
    }
    if(s==0)WARNING("Expression syntax error: ("
                    <<varname[curn]<<" = "<<buffer<<"), variable "
                    <<varname[curn]<<" unchanged");
    switch(vartype[curn])
    {
    case(INT): INFO_buffer(ah,*,int); break;
    case(LONG): INFO_buffer(ah,*,long); break;
    case(DOUBLE): INFO_buffer(ah,*,double); break;
    case(STRING): INFO_buffer(ah,,char); break;
    default: WARNING("unknown vartype ("<<vartype[curn]<<")");
    }

    if(vartype[curn]==STRING)
        if(strcmp(varname[curn],"dirname")==0) /* bug fix (added ==0) Wei 12/19/2006 */
         {
              welcome();
         }

//            if(!diropened)
//            {
//                opendir();
//                printSysInfo();
//                diropened=true;
//            }
            
//    if((!diropened)&&(strlen(dirname)>0))
//    {
//        opendir();
//        printSysInfo();
//        diropened=true;
//    }
    
    return (!s);
}

void Organizer::welcome()
{
   if(!diropened)
   {
       opendir();
       printSysInfo();
       diropened=true;
   }
}

void Organizer::quit()
{
    INFO(HIC"QUIT"NOR);
    printEndInfo();
#ifndef _USETK
#ifdef _PARALLEL
    MPI_Finalize();
#endif
    exit(0);
#endif
}

void Organizer::getsleep()
{
#ifndef _USETK
    if(sleepseconds>0)
        sleep((unsigned)sleepseconds);
#endif
}
void Organizer::init()//system use
{
    SCParser::init();
    bindvar("dirname",dirname,STRING);
    bindvar("sleepseconds",&sleepseconds,INT);
    strcpy(logfilename,"A.log");
}
int Organizer::SysParse(FILE *file)//system use
{
    while(readnextstring(file)==0)
    {
        if(bufferis("setoverwrite")||bufferis("setnolog")
           ||bufferis("dirname"))
        {
            parse_buffer(file);
            if(strlen(dirname)>0)return 0; //successfully get dirname
        }
        else
            FATAL("Organizer expect \"dirname = \" entry!");
    }
    return -1;
}
char * Organizer::currentime()
{
    static char s[100];
    time_t clock;
    clock=time(NULL);
    sprintf(s,"%s",ctime(&clock));
    s[strlen(s)-1]=0;
    return s;
}
void Organizer::printSysInfo()
{
    if(!renew)INFO(HIM"RUN "NOR<<"dirname = "<<dirname);
    else INFO(HIM"RE"HIM"RUN "NOR<<"dirname = "<<dirname);
    INFO(HIM"Begin"HIG" Time = "NOR<<currentime());
}
void Organizer::printEndInfo()
{
    struct rusage ru;
    int pid; char cmd[1000];
    INFO(HIM"End"HIG" Time = "NOR<<currentime());
    getrusage(RUSAGE_SELF,&ru);
        
    INFO_Printf(HIG"CPU time spent:"NOR"  %f s\n",
                ru.ru_utime.tv_sec+1e-6*ru.ru_utime.tv_usec);

#ifdef _USETCL
    /* remove tmp files */
    pid=getpid();
    sprintf(cmd,"rm -f /tmp/*.tmp%d.tcl\n",pid);
    //INFO_Printf(cmd);
    system(cmd);
#endif
    
}
void Organizer::closefiles(int zip, bool bg)
{
    LFile::Close(olog);
    _IO.Restore();
    if(zip==1)
        LFile::GZipFile(logfilename,bg);
    else
        LFile::BZipFile(logfilename,bg);
}

void SubHomeDir(char *fname, char *extname)
{
    char Home[400], *ph, *pt;
    ph=getenv("HOME"); 
    pt=strchr(fname, '~');
    if((pt==NULL)||(ph==NULL))
    {
       if(fname!=extname) strcpy(extname,fname);
//       INFO_Printf("~ not found in %s\n",fname);
       return;
    }
    *pt=0;
    strcpy(Home,fname);
    strcpy(Home+strlen(Home),ph);       
    strcpy(Home+strlen(Home),pt+1);
    strcpy(extname,Home);
    if(fname!=extname)  *pt='~';
    INFO_Printf("replace %s by %s\n",fname,extname);
}

int Organizer::opendir()
{
    char extname[200];
    SubHomeDir(dirname, extname);
        
    /* mkdir `dirname',cd `dirname',open `A.log',redirect `_pSys' */
    DUMP("mkdir "<<dirname<<"  overwrite="<<overwrite);
    if(mkdir(extname,S_IRWXU)!=0)
    {
        if(errno==EEXIST)
        {
            if(overwrite)
            {
                WARNING("directory "<<dirname<<" already exists");
                renew=true;
            }
            else
                FATAL("directory "<<dirname<<" already exists");
        }
        else
        {
            FATAL("directory open "<<dirname<<" failed");
            return -1;
        }
    }
    if(chdir(extname)!=0)
    {
        FATAL("cd "<<dirname<<" failed");
        return -1;
    }
    olog = LFile::Open(logfilename,LFile::O_WriteNew);
    if(!nolog)
        if (olog!=NULL)   LSys::Redirect(olog); /* redirect log stream */
    return 0;
}
