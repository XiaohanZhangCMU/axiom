/*
  scparser.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Wed Jan 21 23:38:40 2009

  FUNCTION  :  Input scripting parser module for MDFrame
*/

#include "general.h"
#include "scparser.h"

void SCParser::init()
{
//    INFO("SCParser::init()");
    bindvar("ncpu",&ncpu,INT);
    bindvar("shmsize",&shmsize,INT);
}

int SCParser::exec(const char *name)
{
    bindcommand_real(name,"abort",abortparser());
    return -1;
}
    
void SCParser::abortparser()
{
    willabort=true;
}

int SCParser::assignvar(int offset)
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
    case(INT): output_buffer(DUMP,*,int); break;
    case(LONG): output_buffer(DUMP,*,long); break;
    case(DOUBLE): output_buffer(DUMP,*,double); break;
    case(STRING): output_buffer(DUMP,,char); break;
    default: WARNING("unknown vartype ("<<vartype[curn]<<")");
    }
    return (!s);
}

void SCParser::dumpbuffer()
{
    DUMP("buffer="<<buffer);
}

bool SCParser::bufferis(const char *s)
{
    return ((strcmp(buffer,s)==0)?true:false);
}
bool SCParser::bufferbeginswith(char c)
{
    if(buffer[0]==c) return true;
    else return false;
}
bool SCParser::bufferendswith(char c)
{
    int L;
    L=strlen(buffer);
    if(L>=MAXTAG_BUF)
    {
        FATAL("SCParser::bufferbeginswith('"<<c<<"'), buffer too long!");
    }
    if(buffer[L-1]==c) return true;
    else return false;
}
bool SCParser::buffercontains(char c)
{
    int i,L;
    L=strlen(buffer);
    if(L>=MAXTAG_BUF)
    {
        FATAL("SCParser::bufferbeginswith('"<<c<<"'), buffer too long!");
    }
    for(i=0;i<L;i++)
        if(buffer[i]==c) return true;
    return false;
}

void SCParser::bindvar(const char *vn,void *p,int type)
{
    for(int i=0;i<varn;i++)
        if(strcmp(varname[i],vn)==0)
        {
            varptr[i]=p;
            vartype[i]=type;
            return;
        }
    sscanf(vn,"%s",varname[varn]);
    varptr[varn]=p;
    vartype[varn]=type;
    varn++;

    if(varn>MAXVAR)
        ERROR("Number of variables exceeds MAXVAR.");
}

/*
FILE * SCParser::getfilehandler(int argc, char *argv[])
{
    static FILE *fp;
    if(argc!=2) { INFO("Usage"<<argv[0]<<"script"); return NULL; }

    if(strcmp(argv[1],"-")==0) fp=stdin;
    else
    {
        fp=fopen(argv[1],"r");
        if(fp==NULL)
        {
            SYSERROR("File "<<argv[1]<<" open error!");
            return NULL;
        }
    }
    return fp;
}
*/
/* allow pre-processor */
FILE * SCParser::getfilehandler(int argc, char *argv[])
{
    static FILE *fp;
#ifndef NOPREPROCESS
    char fname[200], extname[200], com[200];
#endif
    //if(argc!=2) {
    if(argc<2) {
        INFO("Usage:\n    "<<argv[0]<<" (input script file)"
             <<    "\nor  "<<argv[0]<<"  - ");
        return NULL;
    }

    if(strcmp(argv[1],"-")==0) fp=stdin;
    else
    {
#ifndef NOPREPROCESS
        strcpy(fname,argv[1]);
        strcpy(extname,fname);
        sprintf(extname+strlen(extname),".ext");
#ifdef __rshlxcpp
        sprintf(com,"cp %s ~/tmp",fname); system(com);
        sprintf(com,"rsh lx01 \'cd ~/tmp ; cpp %s >! %s\'",fname,extname);
        system(com);
        sprintf(com,"cp ~/tmp/%s .",extname); system(com);
#else
        sprintf(com,"cpp %s > %s 2> /dev/null",fname,extname);
        printf("%s\n",com);
        system(com);
#endif
        fp=fopen(extname,"r");
        sprintf(com,"/bin/rm %s",extname);system(com);
#else
        fp=fopen(argv[1],"r");
#endif
        if(fp==NULL)
        {
            SYSERROR("File "<<argv[1]<<" open error!");
            return NULL;
        }
    }
    return fp;
}

int SCParser::readnextstring(FILE *file) //into buffer
{
    int i=0;
    int state=0;
    char c; int fs; 
    /* state
       0: string not begin, no comment
       1: string not begin, in comment
       2: string begin
       4: string end with comment
    */
    while (1)
    {
        c=cchar;
        
        fs=fscanf(file,"%c",&cchar);
        if(fs==EOF) break;

        if(c==-1) continue;
        
        switch (state)
        {
        case 0: /* 0: string not begin, not in comment */
            switch (c)
            {
            case ' ':
            case '\t':
            case '\n':
            case '\r': break;
            case '#': state=1;break;
            case '\"': state=5; break;
            default: buffer[i++]=c;state=2;break;
            }
        case 1: /* 1: string not begin, start comment */
            if (c=='\n') state=0; break;
        case 2: /* 2: string begin */
            switch (c)
            {
            case ' ':
            case '\t':
            case '\n':
            case '\r': buffer[i]=0; return 0;
            case '#': buffer[i]=0; state=4;break;
            default:
                if(i<MAXTAG_BUF-1) buffer[i++]=c;
                else  {
                    buffer[i++]=c;
                    buffer[i]=0;
                    return 0;
                }
                break;
            }
        case 4: /* 4: string end with comment */
            if (c=='\n') return 0;
            break;
        case 5: /* 5: begin quotation mark */
            switch (c)
            {
            case '\"': buffer[i]=0; return 0;
            default:
                if(i<MAXTAG_BUF-1) buffer[i++]=c;
                else {
                    buffer[i++]=c;
                    buffer[i]=0;
                    return 0;
                }
            }
        }
        /* new in version 2: characters like + [ ] always form
         * single character string */
        if(state==2)
        {
            switch(c)
            {
            case '=':
            case '[':
            case ']': buffer[i]=0; return 0;
            }
            switch(cchar)
            {
            case '=':
            case '[':
            case ']': buffer[i]=0; return 0;
            }
        }
    }
    return -1;
}

int SCParser::identify(const char *name)
{
    int i;
    char shortname[1000], *p, *q;

    shift=0;
    strcpy(shortname,name);
    p = strchr(shortname,'(');
    if(p!=NULL) 
    {
        *p=0; q = strchr(p+1,')');
        if(q!=NULL)
        {
            *q=0;
            sscanf(p+1,"%d",&shift);
            //INFO_Printf("name=%s shift=%d\n",name,shift);
        }
    }
    for(i=0;i<varn;i++)
    {
        if(strcmp(shortname,(const char *)varname[i])==0)
        {
            curn = i; return curn;
        }
    }
    curn = -1; return curn;
}

int SCParser::parse(FILE *file)
{
    //INFO_Printf("SCParser::parse\n");
    if (file==NULL)
    {
        ERROR("parse: file==NULL!");
        return -1;
    }
    while((readnextstring(file)==0)&&(!willabort))
    {
        parse_buffer(file);
    }
    return 0;
}

int SCParser::parse_buffer(FILE *file)
{
    int offset; bool complete; 
 id:
    /* variable name or command read in */
    identify(buffer);
    if(curn<0) /* it's a command to be executed */
    {
        int r=exec(buffer);
        if(r!=0) WARNING("unrecognized command " HIR <<buffer<< NOR " ignored");
        return r;
    }
    else
    { /* it is a variable to be assigned value */
        /* read it's value */
        if(readnextstring(file)!=0)
        {
            FATAL("unexpected end of file");
            return(-1);
        }
        /* expecting an `=' sign */
        if(!bufferis("="))
        {
            WARNING("expecting =, resume from identifier ["<<buffer<<"]");
            goto id;
        }
        if(readnextstring(file)!=0)
        { /* no value read in */
            FATAL("unexpected end of file");
            return(-1);
        }
        /* if single value read */
        if(!bufferis("["))
        {
            return assignvar();
        }
        else
        { /* if `[' read expecting an array */
            DUMP("Read Array...");
            offset=0; complete=false;
            while(!complete)
            {
                if(readnextstring(file)!=0)
                {
                    FATAL("unexpected end of file");
                    return(-1);
                }
                if(!bufferis("]"))
                {
                    assignvar(offset);
                    offset++;
                }
                else
                    complete=true;
            }
            return 0;
        }
    }
}


