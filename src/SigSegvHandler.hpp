#ifndef _PGRAPH_SIGSEGV_HANDLER_
#define _PGRAPH_SIGSEGV_HANDLER_

static void (*SigSegvOrig)(int) = NULL;

static void SigSegvHandler(int sig)
{
    fprintf(stderr,"(%d): Segmentation Violation ... pausing\n",
            getpid() );
    pause();
    if (SigSegvOrig == SIG_DFL) {
        signal(sig, SIG_DFL);
    }
    else if (SigSegvOrig == SIG_IGN) {
    }
    else {
        SigSegvOrig(sig);
    }
}

static void TrapSigSegv()
{
    if ((SigSegvOrig=signal(SIGSEGV, SigSegvHandler)) == SIG_ERR) {
        fprintf(stderr, "TrapSigSegv: error from signal setting SIGSEGV");
        exit(EXIT_FAILURE);
    }
}

#endif /* _PGRAPH_SIGSEGV_HANDLER_ */
