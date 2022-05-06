#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include  <math.h>

#include <nicklib.h>
#include <getpars.h>
#include <globals.h> 

#include "badpairs.h"
#include "admutils.h"
#include "mcio.h"  
#include "mcmcpars.h"  
#include "egsubs.h"  
#include "ldsubs.h"  
#include "regsubs.h" 
#include "qpsubs.h" 


#define WVERSION   "4010"   

// mincount stripped 
// basic idea: Model each admixed sample as a linear fit 
// of 2 parental pops.  Compute residual and rolloff (res*wt) 
// jackknife by chromosome
// weights added as option.  Latest I/O (mcio.c)
// timeoffset added 
// mesh + FFT convolutions.  slightly faster FFT
// zero popsize trapped.  runfit option
// admixlist introduced  
// snpname, seed passed into rolljack
// fft bug fixed (very low coverage) 
// PM: modified rolljack call to jackknife
// covariance output in col 1 of output file (dumpit) 
// specnames etc 
// zero snps working trapped 
// 3 sources with whichadm implemented.  1 => near adm , 2 => far adm
// numchrom added
// usehets added
// nrmsd added 
// warning given if admixlist not set.
// bugfix if numchrom > 22 
// better implementation of cwd() 

#define MAXSTR  512
#define MAXAPOPS  100

#define BINSIZE  .0005
#define MAXDIS   0.05  

extern int packmode ;
int zdipcorrmode = NO ; 
int regdim = 2 ; 
int selfsource = YES ; 
int samecoeffs = NO ; 
int singlesampleoutput = NO ;
int usehets = YES ; 

char *trashdir = "/var/tmp" ;
int qtmode = NO ;
int debug = NO ;

char **specnames = NULL ; 
int numspec=0 ; 

int *speccols ; 
char *specfile=NULL ;  
char *specstring=NULL ;  

Indiv **indivmarkers;
SNP **snpmarkers ;
SNP **snpm ;
int numsnps, numindivs ; 
int checkmap = YES ;
double *fc ; 

char *parname = NULL ;

char  *genotypename = NULL ;
char  *snpname = NULL ;
char *snppath = NULL ; 
char  *genooutfilename = NULL ;
char  *indoutfilename = NULL ;
char  *indivname = NULL ;
char *badsnpname = NULL ;
char *poplistname = NULL ;
char *admixpop = NULL ;
char *admixlist = NULL ;
char **admixpoplist = NULL  ;
int nadmixpop = 0; 
char *dumpname = NULL ;
char *loadname = NULL ;
char *oname = NULL ;
char *weightname = NULL ; 
char *timeoffsetname = NULL ; 
int crossmode = NO ;
int runmode = 0 ;
int ldmode = NO ;
int ransample = -1 ;
int flatmode = NO  ;
double chithresh = -1.0; // default 6.0 
int jackknife = YES ;
int qbin = 0 ;
int usevar = YES ;

FILE *ofile ;
int  *fchrom, *lchrom ;
int nums = 2,   ns2, ns4 ;
int numbigiter = 2 ;
int xchrom = -1 ;
int xnochrom = -1 ;
int seed = 0 ;

int noxdata = YES ;
int samelambda = NO ;
int fixlambda = NO ;
int flatinit = NO ;  // different default
int popestimate = NO ; 
	int nxlim = 100000 ;
int minparentcount = 10 ;
int norun = NO ;

double maxdis = MAXDIS   ;
double binsize = BINSIZE ; 

// if YES same across pops
int *firstsnp, *lastsnp ; 
int *firstqb, *lastqb ; 
int *qindex ;  

// expfit parameters 
int afffit = YES ; 
double lovalfit = .45 ;   
char *logfit = NULL ;
int runfit = YES ;

int whichadm = -1 ;  


char  unknowngender = 'U' ;

int readcommands(int argc, char **argv) ;
void setfc(SNP **snpm, int numsnps) ;
void sethmm(double *theta, double lambda, int ind) ;
int addr4(int a, int b, int c, int d) ;
double doalpha(int ind) ;
double dobeta(int ind) ;

void cleargcounts()   ;
void reest()   ;
void calcgcounts(int ind)  ;
void printpars()   ;
void printgams()   ;
void calcgtl(int ind) ;
int  setqbins(SNP **snpm, int numsnps, double binsize, int qbin)  ; 
double doscore(double **thet, double *lam) ;
void dumppars(char *dumpname,  double **wmix, double **theta, int nums, int numeg) ;
void dorc(double *ans, double *res, double binsize, double maxdis) ;
void dorc2(double *ans, double *res, double binsize, double maxdis) ;
void dorcx(double *ans, double *rr, double *dd, double binsize, double maxdis) ;
void dorcxx(double *ans, double *rr, double *dd, double binsize, double maxdis) ;
void dumpit(char *oname, double *ww, CORR **corrbins, int len, double bsize, char *hdr) ;
double cntit1(double *xc, SNP *cupt, Indiv **indm, int numindivs, int t) ; 
double cntit2(double *xc, SNP *cupt, SNP *cupt2, Indiv **indm, int numindivs, int t) ; 
double estyran(SNP **snpmarkers, int numsnps, int nxlim) ;
void cx1(double *xc, double *y1, double *y2) ;
double cminor(double *xx) ;
int calcldscore(SNP *cupt, SNP *cupt2, double *xscore) ;
void printzq(double *z0q, double *z1q, double *z2q, int numqbins) ;
void ddadd(double ***ddcbins, double *z0q, double *z1q, double *z2q, int numqbins, int diffmax)  ; 
void ddcorr(CORR ***corrjbins, double ***ddcbins, int numdbins, int qbin) ;   
void fixjcorr(CORR **corrbins, CORR ***corrjbins, int numbins)  ;

void slocv(double *cout, double *a1, double *b1, int n, int maxlag) ;
void fftcv(double *cout, double *a1, double *b1, int n, int maxlag) ;
void fftauto(double *cout, double *a1, int n, int maxlag) ;
void fftcv2(double *cout, double *cout2, double *a1, double *b1, int n, int maxlag) ;
void dorunfit() ;
int dowork(char **eglist, int numeg) ; 

int calcab(double *ya, double *yb, int *valids,  SNP **snpmarkers, int numsnps, 
  Indiv ** indivmarkers, int numindivs, int admindex, double *ccoeffs, int numeg) ;

int calcpred(double *target, double *pred, int *valids, SNP **snpmarkers, int numsnps, 
 Indiv ** indivmarkers, int numindivs, int admindex, double *zans) ;

double calccoeffs(double **pcoeffs, SNP **snpmarkers, int numsnps, Indiv **indivmarkers, 
 int numindivs, int *admindexes, int ncount, int numeg) ;
void getsn(char **specnames, int numspec, char *specstring)  ;

double dd0sum = 0 ; 
double dd0count = 0 ; 
int numdbins ; 

int main(int argc, char **argv)
{

  double ymem ;  
  int ret, i, j, k, m, s, copyret ; 
  int tt, t, chrom ; 
  int nignore, numrisks = 1 ;
  char **eglist, *sx ;
  int numeg ; 
  char ***xjobs ; 
  int  njobs ; 
  char *dir ;
  char *cwd  ;
  char sss[1024]  ;
  SNP *cupt ; 
  int fdescwd = -1 ;
  
  char **specsnps ; 
  double **xtemp, **xprobs=NULL ;  
  int nspecsnps = 0 ; 
  double *specprobs ; 
  double fcoeffs[3] ; 
  int numfc = 0 ; 




  ofile = stdout; 
  packmode = YES ;
  

  printnl() ;
  printf("## DATES.  Version %s\n", WVERSION) ;
  printnl() ;

  ZALLOC(speccols, 10, int) ;
  readcommands(argc, argv) ;

  if (parname == NULL) return 0 ;

  if (usevar) printf("usevar set!\n") ; 

  printnl() ;
  cputime(0) ;
  calcmem(0) ;
 
  if ((regdim==4) && (selfsource == NO)) regdim = 3 ;
  if (regdim<4) selfsource = NO ; 
  

  if (numspec>0) { 
    ZALLOC(specnames, numspec, char *) ; 
    getsn(specnames, numspec, specstring) ;
    for (k=0; k<numspec; ++k) { 
     printf("spec: %12s %d\n", specnames[k], speccols[k]) ;
    }
    t = numlines(specfile) ;
    ZALLOC(specsnps, t, char *) ; 
    s = numcols(specfile) - 1 ; 
    printf("number of data cols: %d\n", s) ;
    xtemp = initarray_2Ddouble(s, t, 0) ; 
    xprobs = initarray_2Ddouble(numspec, t, 0) ; 
    nspecsnps = getxxnames(&specsnps, xtemp, t, s, specfile) ;
    printf("nspecsnps: %d\n", nspecsnps) ;  
    for (j=0; j<numspec; j++) { 
     k = speccols[j] ;  
     copyarr(xtemp[k-1], xprobs[j], nspecsnps) ; // data format includes col 0 = snpname
    }
    free2D(&xtemp, s) ;
  }

  if (seed == 0) seed = seednum() ;
  SRAND(seed) ;
  printf("seed: %d\n", seed) ;

  if (snpname[0] != '/') { 
   snppath = realpath(snpname,  NULL) ;
   if (snppath == NULL) fatalx("realpath can't find snpname!\n") ;

   printf("snppath: %s\n", snppath) ; 
   fflush(stdout) ; 
  }

  numsnps = 
    getsnps(snpname, &snpmarkers, 0.0, badsnpname, &nignore, numrisks) ;

  numindivs = getindivs(indivname, &indivmarkers) ;


  printf("genotypename:  %s\n", genotypename) ;
  fflush(stdout) ;

   getgenos(genotypename, snpmarkers, indivmarkers, 
    numsnps, numindivs, nignore) ;


  for (i=0; i<numsnps; i++)  {

    cupt = snpmarkers[i] ; 
    cupt -> tagnumber = i ;
    chrom = cupt -> chrom ;
    if (chrom == 0) cupt -> ignore = YES ;
    if (chrom > numchrom)  cupt -> ignore = YES ;
    if ((xchrom>0) && (chrom != xchrom)) cupt -> ignore = YES ;
    if (chrom == xnochrom) cupt -> ignore = YES ;
    if (nspecsnps>0) cupt -> isfake = 1 ;  // temporary hack 
  }

  for (j=0; j<nspecsnps; j++) { 
   sx = specsnps[j] ; 
   k = snpindex(snpmarkers, numsnps, sx) ; 
   if (k<0) continue ; 
   cupt = snpmarkers[k] ; 
   cupt -> isfake = 0 ; 
   cupt -> isrfake = j ;  // overload free field
    ZALLOC(cupt -> modelscores, numspec, double) ;  
    specprobs = cupt -> modelscores ; 
    for (m=0; m<numspec; m++) { 
     specprobs[m] = xprobs[m][j] ;  
    }
  }

  if (xprobs!=NULL) free2D(&xprobs, numspec) ; 

  t = 0 ; 
  for (i=0; i<numsnps; i++)  {
    cupt = snpmarkers[i] ; 
    if (cupt -> isfake == 1) cupt -> ignore = YES ; 
    if (cupt -> ignore == NO) ++t ; 
  }
   printf("valid snps: %d\n", t) ;

   numsnps = rmsnps(snpmarkers, numsnps, NULL) ;   
/**
   numindivs = rmindivs(snpmarkers, numsnps, indivmarkers, numindivs) ;
*/  

   if (checkmap && (cmap(snpmarkers, numsnps) == NO)) { 
    printf("running DATES without a real map.  If you mean to do this set checkmap: NO\n") ;
    fatalx("no real map\n") ;
   }

  vclear(fcoeffs, -99, 3) ; 
  numfc = 0 ; 

  if (admixlist == NULL) { 
   printf("*** Stronhly advise use of admixlist:  Obsolescent parameters\n") ;  
   ZALLOC(eglist, numindivs+1, char *) ; 
   numeg = loadlist(eglist, poplistname) ;
   if (numeg != 2) fatalx("must have exactly 2 admixing populations! poplistname: %s\n", poplistname) ;
   eglist[numeg] = strdup(admixpop) ;
   fflush(stdout) ;
   ret = dowork(eglist, numeg) ; 
 }

  if (admixlist != NULL)  { 
   numfc = 0 ; 
   t = tt = numcols(admixlist) ;  // : stripped 
   if (tt==6) { 
    t = 4 ;  numfc = 2 ; 
   } 
   if (tt==8) { 
    t = 5 ;  numfc = 3 ; 
   } 
   
   numeg = t-2  ; 
   if ((numeg<2) || (numeg>3))  fatalx("bad number of columns in %s\n", admixlist) ;
   if (whichadm == -1) whichadm = 1; 
   ZALLOC(eglist, numeg+1, char *) ; 

// ZALLOC(cwd, 1025, char) ;
// getcwd(cwd, 1024) ; 
   fdescwd = open(".", O_RDONLY) ;

   njobs = numlines(admixlist) ; 
   ZALLOC(xjobs, tt, char **)  ; 
   for (k=0; k<tt; k++) { 
    ZALLOC(xjobs[k], njobs, char *) ; 
   } 

   njobs = getnames(&xjobs, njobs, tt, admixlist) ;

   for (i=0; i<njobs; i++) {
    for (k=0; k<numeg; ++k) { 
     eglist[k] = xjobs[k][i] ;  
    }

    eglist[numeg] = admixpop  = xjobs[numeg][i] ; 

    fc = NULL ; 
    if (numfc>0) { 
     for (j=0; j<numfc; ++j) { 
      sx = xjobs[numeg+j+2][i] ; 
      fcoeffs[j] = atof(sx) ; 
     }
     bal1(fcoeffs, numfc) ; 
     fc = fcoeffs ;
    } 

    freestring(&oname) ; 
    freestring(&logfit) ; 
    
    dir = xjobs[numeg+1][i] ; 

    sprintf(sss, "mkdir -p -m 0755 %s\n", dir) ;
    system(sss) ; 
    sprintf(sss, "chmod  0755 %s\n", dir) ;
    system(sss) ; 


    sprintf(sss, "cp %s %s", parname, dir) ; 
    copyret = system(sss) ; 
    // printf("copy param code: %d\n", copyret) ;

    chdir (dir) ; 
    printf("calling dowork: ") ; 
    for (k=0; k<=numeg; ++k) { 
     printf(" %s ", eglist[k]) ;
    }
    printnl() ; 
    fflush(stdout) ;
    ret = dowork(eglist, numeg) ;  
    fchdir(fdescwd) ; 
   }
  }

  ymem = calcmem(1)/1.0e6 ;
  printf("##end of dates: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;
  return ret ;

}
int dowork( char **eglist, int numeg)  
{

  char s1[MAXSTR], s2[MAXSTR] ;
  int x, x2, i, j, k, g, s, t, u, d, tt, iter ; 
  int a, b, c, j1, j2, k1, k2  ;
  SNP *cupt, *cupt2 ;
  Indiv *indx, **indadmix ;
  char dirname[MAXSTR] ;  

  double y, y1, y2, yy, ywt, z1, z2, znorm, yn, ytime, yres, yvar, ypred ;
  double *yy2 ;
  int ind, chrom, lastchrom, printit, num0=0 ;
  double xc[9], xww[3] ;
  char sss[256], *sx ;

  double *ttt ;
  double yscore, ybase ;
  double lambdabase ;
  double *scarray, *wacc, ywynn ;
  int ncount = 0 ;
  double lastscore = 0.0 ;
  double *tvecs, *pp ;
  double *w0, *w1, *w2, *res,  *ww, *ww1, *ww2, *dd, *wt, *wk  ;
  double *wa, *wb ;  
  double *yya, *yyb ; 
  int *sindex, *xindex ; 
  double *xnum, *xden, *xdenq ;
  double *znum, *zden ;
  double **xxnum, **xxden ;
  int numbins  ;  
  double co1, co2, dis, yran ;     
  double wtmean, wtnum ; 
  CORR **corrbins ;
  CORR ***corrsbins ;
  CORR *corrpt  ;
  CORR ***corrjbins ;  //jackknife
  CORR ****corrsjbins ;
  int nloop = 0 ;
  int nncount = 0 ;
  int nsnp1 = 0, nsnp2 = 0 ;
  double **jwt ;
  int numx ;  // snps in main double loop
  double timeoffset ;  
 
  int numqbins ; 
  double *z0q, *z1q, *z2q ; 
  double *zs0q, *zs1q, *zs2q ; 
  double **ddbins, ***ddcbins ; 
  double ***ddsbins, ****ddscbins ; 
  int diffmax ; 
  double ymem ;  
  int xcnt[10] ;  // admixing pop sizes
  int bcount[3] ;
  double bfreq ;  
  double *bf, *gf, ysig ; 
  int *valids ; 
  int *admindices ;
  int icount  ; 
  double **pcoeffs ;
  double zans[4] ; 
  int numadmgood ;
  double alpha, beta, *ccoeffs ; 
  int isok ; 

  char *admixpop ; 


  double **freqs ;  
  int ***counts, nrows, ncols ;
  int *zindex, *ztypes ;

  admixpop = eglist[numeg] ; 
  printnl() ;
  printstrings(eglist, numeg) ;  printf("admix: %s\n", admixpop) ;
  printnl() ;

  for (i=0; i<numsnps; i++) { 
   cupt = snpmarkers[i] ; 
   cupt -> ignore = NO ;
  }

  nrows = numindivs; ncols = numsnps ; 
  ZALLOC (zindex, numindivs, int);
  ZALLOC (ztypes, numindivs, int);
  for (i=0; i<numindivs; ++i) { 
   zindex[i] = i ; 
   indx = indivmarkers[i] ; 
   k = indxindex (eglist, numeg+1, indx->egroup);
   ztypes[i] = k;
  }

 ZALLOC(counts, numsnps, int **);
  for (k = 0; k < ncols; ++k) {
    counts[k] = initarray_2Dint (numeg+1, 2, 0);
  }
  freqs = initarray_2Ddouble(numsnps, numeg+1, -1) ; 

  countpops (counts, snpmarkers, zindex, ztypes, nrows, ncols);
  printf ("countpops called\n");

  for (k = 0; k < ncols; ++k) {
   cupt = snpmarkers[k] ; 
   isok = NO ; 
   for (j = 0; j <=  numeg; ++j) {
    y1 = counts[k][j][0] ; 
    y2 = counts[k][j][1] ; 
    yy = y1+y2 ; 
    if (yy<0.1) { 
     cupt -> ignore = YES ; 
     continue ; 
    }
    cupt -> cf[j] = freqs[k][j] = y2/yy ;  // flatten ?? .  Derived allele freq

    if (j<numeg) { 
     y = cupt -> cf[j] - cupt -> cf[0] ; 
     if (fabs(y) > .001) isok = YES ; 
    }   
  } 
  if (isok == NO) cupt -> ignore = YES ; 
  if (verbose) { 
   printf("zzvv %d %d %d ", k, cupt -> ignore, isok) ; 
   printmat(cupt -> cf, 1, numeg+1) ;  
  }
 } 

  free(ztypes) ; 
  free(zindex) ; 
  for (k = 0; k < ncols; ++k) {
    free2Dint(&counts[k], numeg+1) ;
  }
  free(counts) ;

  ZALLOC(admixpoplist, MAXAPOPS, char *) ;  
   if (admixpop == NULL) fatalx("no admixpop\n") ; 
    nadmixpop = 1 ;
    admixpoplist[0] = strdup(admixpop) ;

  if (oname == NULL) { 
    sprintf(sss, "%s.out", admixpop) ;  
    oname = strdup(sss) ;  
    printf("output tables: %s\n", oname) ;
  }


//   numsnps = rmsnps(snpmarkers, numsnps, NULL) ;   

   setfc(snpmarkers, numsnps) ; // does allocation
   numqbins = setqbins(snpmarkers, numsnps, binsize, qbin) ;       
   if (numqbins == 0) printf("no mesh\n") ; 
   else { 
    printf("allocating qindex. numqbins: %d\n", numqbins) ; 
    fflush(stdout) ; 
    ZALLOC(qindex, numqbins, int) ;
    printf("mesh used :: qbin: %d  numqbins: %d\n", qbin, numqbins) ;
    fflush(stdout) ; 
   }

   ncount = 0 ;
   ivzero(xcnt, numeg+1) ;
   ZALLOC(indadmix, numindivs, Indiv *) ; 
   for (i=0; i<numindivs; ++i) { 
    indx = indivmarkers[i] ;
    indx -> affstatus = -1 ;
    indx -> ignore = NO ; 
    t = indxindex(eglist, numeg+1, indx -> egroup) ;
    if (t>=0) { 
     indx -> affstatus = t ;
     ++xcnt[t] ;  
    }
    tt = indxindex(admixpoplist, nadmixpop, indx -> egroup) ;
    if (tt>=0) { 
     indx -> affstatus = 99 ;
     indadmix[ncount] = indx ;
     ++ncount ;   
    }
    if (indx -> affstatus  <  0) indx -> ignore = YES ;
    if (indx -> ignore) continue ;
   }
   if (ncount==1) selfsource = NO ; 
   if (selfsource == NO) regdim = MIN(regdim, 3) ;
   regdim = MAX(regdim, numeg) ;
   
    printf("regdim: %d\n", regdim) ; 
    for (u=0; u<numeg; ++u) {  
     sx = eglist[u] ; 
     if (xcnt[u]>0) continue ;  
     s = indxindex(specnames, numspec, sx) ; 
     if (s<0) fatalx("admixing pop %s has no samples\n", sx) ;
    }

   if (ncount == 0) fatalx("no admixed samples found!\n") ;
   printf("number admixed: %d numeg: %d\n", ncount, numeg) ;
   ZALLOC(admindices, ncount, int) ; 

   for (i=0; i<ncount; ++i) { 
    k  = indadmix[i] -> idnum ; 
    admindices[i] = k ;
   } 
   
   for (k=0; k<numeg; ++k) { 
     printf("group %d %15s samplesize: %4d\n", k, eglist[k], xcnt[k]) ;
   }

//   numindivs = rmindivs(snpmarkers, numsnps, indivmarkers, numindivs) ;

   printf("numsnps: %d  numindivs: %d\n", numsnps, numindivs) ;

   if (numsnps == 0)   { 
     printf("zero snps -- misspelled population? -- exiting\n") ; 
     return -1 ; 
   }

// allele counts for admixed pop
   ZALLOC(gf, numsnps, double) ;
   ZALLOC(bf, numsnps, double) ;

   ZALLOC(yya, numsnps, double) ;
   ZALLOC(yyb, numsnps, double) ;

   for (k=0; k<numsnps; ++k) { 
    cupt = snpmarkers[k] ; 
    cntit1(xc, cupt, indivmarkers, numindivs, 99) ;  
    fixit( cupt -> bcount, xc, 3) ; 
   }

   numbins = nnint(maxdis/binsize) + 5 ;
   ZALLOC(w0, numsnps, double) ;
   ZALLOC(w1, numsnps, double) ;
   ZALLOC(w2, numsnps, double) ;
   ZALLOC(wk, numsnps, double) ;
   ZALLOC(ww1, numsnps  + numbins, double) ;
   ZALLOC(ww2, numsnps  + numbins, double) ;
   ZALLOC(wa, numsnps, double) ;
   ZALLOC(wb, numsnps, double) ;
   ZALLOC(res, numsnps, double) ;
   ZALLOC(dd, numsnps, double) ;
   ZALLOC(wt, numsnps, double) ;

   ZALLOC(sindex, numsnps, int) ;
   ZALLOC(xindex, numsnps, int) ;

   ZALLOC(ww, numsnps +numbins, double) ;

   ZALLOC(xnum, numbins, double) ;
   ZALLOC(xden, numbins, double) ;
   ZALLOC(xdenq, numbins, double) ;
   ZALLOC(znum, numbins, double) ;
   ZALLOC(zden, numbins, double) ;
   ZALLOC(yy2, numbins, double) ;

    xxnum = initarray_2Ddouble(numchrom+1, numbins, 0) ;
    xxden = initarray_2Ddouble(numchrom+1, numbins, 1.0e-10) ;

    if (numqbins>0) { 
     ZALLOC(z0q, numqbins, double) ;  
     ZALLOC(z1q, numqbins, double) ;  
     ZALLOC(z2q, numqbins, double) ;  

     ZALLOC(zs0q, numqbins, double) ;  
     ZALLOC(zs1q, numqbins, double) ;  
     ZALLOC(zs2q, numqbins, double) ;  

     numdbins = numbins * qbin ; 
     ddbins = initarray_2Ddouble(7, numdbins, 0) ;   
     ZALLOC(ddcbins, numchrom+1, double **) ; 
     for (k=0; k<=numchrom; ++k) { 
      ddcbins[k] = initarray_2Ddouble(7, numdbins, 0) ;   
     }
     ZALLOC(ddsbins, ncount, double **) ;  
     ZALLOC(ddscbins, ncount, double ***) ;
     for (a=0; a<ncount;++a) { 
      ddsbins[a] = initarray_2Ddouble(7, numdbins, 0) ;   
      ZALLOC(ddscbins[a], numchrom+1, double **) ;
      for (k=0; k<=numchrom; ++k) { 
       ddscbins[a][k] = initarray_2Ddouble(7, numdbins, 0) ;   
      }
     }
    }

   ZALLOC(firstsnp, numchrom+1, int) ;
   ZALLOC(lastsnp, numchrom+1, int) ;
   ZALLOC(firstqb, numchrom+1, int) ; // first qbin in chromosome
   ZALLOC(lastqb, numchrom+1, int) ;

   vclear(ww, 1.0, numsnps) ;

   ZALLOC(corrbins, numbins, CORR *) ;

   for (x=0; x<numbins; ++x) { 
     ZALLOC(corrbins[x], 1, CORR) ;
     corrpt = corrbins[x] ; 
     clearcorr(corrpt) ;
   }

   ZALLOC(corrjbins, numchrom+1, CORR **) ;
   jwt = initarray_2Ddouble(numchrom+1, numbins, 0.0) ; 
   for (k=0; k<=numchrom; ++k) {
    ZALLOC(corrjbins[k], numbins, CORR *) ;

    for (x=0; x<numbins; ++x) { 
     ZALLOC(corrjbins[k][x], 1, CORR) ;
     corrpt = corrjbins[k][x] ; 
     clearcorr(corrpt) ;
    }
   }

   ZALLOC(corrsbins, ncount, CORR **) ; 
   ZALLOC(corrsjbins, ncount, CORR ***) ; 

  for (a=0; a<ncount; ++a) { 
   ZALLOC(corrsjbins[a], numchrom+1, CORR **) ;
   ZALLOC(corrsbins[a], numbins, CORR *) ;
   for (x=0; x<numbins; ++x) { 
     ZALLOC(corrsbins[a][x], 1, CORR) ;
     corrpt = corrsbins[a][x] ; 
     clearcorr(corrpt) ;
   }
   for (k=1; k<=numchrom; ++k) {
    ZALLOC(corrsjbins[a][k], numbins, CORR *) ;
    for (x=0; x<numbins; ++x) { 
     ZALLOC(corrsjbins[a][k][x], 1, CORR) ;
     corrpt = corrsjbins[a][k][x] ; 
     clearcorr(corrpt) ;
    }
   }
  }
   

 yn = (double) numsnps ;

 regdim = MAX(regdim, numeg) ; 
 pcoeffs = initarray_2Ddouble(ncount, regdim, 0) ; 

 if (fc == NULL) { 
   calccoeffs(pcoeffs, snpmarkers, numsnps, indivmarkers, numindivs, admindices, ncount, numeg) ;
 } 
 else { 
  for (k=0; k<ncount; ++k) {  
   copyarr(fc, pcoeffs[k], numeg) ;
  }
 }

   ZALLOC(valids, numsnps, int) ;

   for (icount=0; icount < ncount; ++icount) { 

    alpha = pcoeffs[icount][0] ; 
    k = admindices[icount] ; 
    ccoeffs = pcoeffs[icount] ; 
    bal1(ccoeffs, numeg) ; 
    t = calcab(yya, yyb, valids, snpmarkers, numsnps, indivmarkers, numindivs, k, ccoeffs, numeg) ; 
    if (t<0)  { 
      printf("*** admix fails :: ") ;
      printstrings(eglist, numeg) ;  printf("admix: %s\n", admixpop) ;
// should free lots of memory here  
      return -1 ;  
    }
   
    if (singlesampleoutput) {
     vzero(zs0q, numqbins) ;
     vzero(zs1q, numqbins) ;
     vzero(zs2q, numqbins) ;
    }

    k = admindices[icount] ; 
    indx = indivmarkers[k] ;
    if (indx -> ignore) continue ;  

    x = 0 ; 

    ivclear(sindex, -1, numsnps) ;
    ivclear(firstsnp, 1000*1000*1000, numchrom+1) ;
    ivclear(lastsnp, -1, numchrom+1) ;
    ivclear(firstqb, 1000*1000*1000, numchrom+1) ;
    ivclear(lastqb, -1, numchrom+1) ;
     
    for (j=0; j<numsnps; ++j) {  
     if (valids[j] == NO) continue ; 
     cupt = snpmarkers[j] ; 

     y = yya[j] - yyb[j] ; 
     
     wt[x] =  1.0  ;  
     res[x] = y  ; 

     t = cupt -> chrom ; 
     firstsnp[t] = MIN(firstsnp[t], j) ; 
     lastsnp[t] =  MAX(lastsnp[t], j) ; 

     if (weightname != NULL) wt[x] = cupt -> weight ;  
     if (isnan(wt[x])) fatalx("bad weight!\n") ;

     sindex[j] = x ;
     xindex[x] = j ; 
     ++x ; 
    }

    numx = x ; 
    printf("numx: %d\n", numx) ; 

    if (numx==0) fatalx("no valid snps!\n") ;


    vzero(z0q, numqbins) ;
    vzero(z1q, numqbins) ;
    vzero(z2q, numqbins) ;

   if (numqbins>0)  {

     ivclear(qindex, -1, numqbins) ;

     for (k1=0; k1<numx; ++k1)  { 

      j1 = xindex[k1] ;
      cupt = snpmarkers[j1] ;

      if (cupt -> ignore) continue ;
      if (valids[j1] == NO) continue ; 
      s = cupt -> tagnumber  ; 
      y = res[k1] * wt[k1]  ; 
      qindex[s] = k1 ; // just 1 snp maps

      t = cupt -> chrom ; 

      firstqb[t] = MIN(firstqb[t], s) ;
      lastqb[t] =  MAX(lastqb[t], s) ;

      z0q[s] += 1 ; 
      z1q[s] += y ; 
      z2q[s] += y*y ; 
    } 

   y = (double) qbin * maxdis / binsize ;
   diffmax = nnint(y) ; 
   ddadd(ddcbins, z0q, z1q, z2q, numqbins, diffmax) ;
   if (singlesampleoutput) { 
    ddadd(ddscbins[icount], z0q, z1q, z2q, numqbins, diffmax) ;
   }
  }

  if (numqbins == 0) printf("starting main loop. ID %s numsnps: %d\n", indx -> ID, numsnps) ;
  lastchrom = -1 ; 
  for (k1=0; k1<numx; ++k1)  { 
   if (numqbins > 0) break ; 
   j1 = xindex[k1] ;
   cupt = snpmarkers[j1] ;
   if (cupt -> ignore) continue ;
   if (valids[j1] == NO) continue ;
   chrom = cupt -> chrom ;
   printit = -1 ;
   if (lastchrom == -1) lastchrom = chrom ; 
   if (lastchrom != chrom) { 
    printit = chrom ; 
    lastchrom = chrom ;
   }
    for (k2=k1+1; k2 < numx; ++k2) { 
     j2 = xindex[k2] ;
     cupt2 = snpmarkers[j2] ;
     if (cupt2 -> ignore) continue ;
     if (valids[j2] == NO) continue ;
      t = cupt2 -> chrom - cupt -> chrom ;
      if (t != 0) break ;
      dis = cupt2 -> genpos - cupt -> genpos ;
      if (dis<0.0) fatalx("badbug\n") ;
      if (dis >= maxdis) break ;
      s = (int) (dis/binsize) ; // as in simpsim2 
      ytime = 1.0 ; 

      yy  = res[k1]*res[k2] ;                    
      ywt = wt[k1]*wt[k2] ;                     

      znum[s] += ywt * yy ;  
      zden[s] += ywt * ywt ;

      xnum[s] += ywt * yy ;  
      xden[s] += ywt * ywt ;

      t = ranmod(1000*1000) ; 
      if (t==-1)  {  
       printf("zza %12.6f ", dis) ;  
       printf("%d %12.6f ", s, y) ;  
       y1 = (double) s * binsize ;  
       printf("%12.6f ", exp(-timeoffset*y1) ) ; 
       printf(" %12.6f %12.6f", xnum[s]/znum[s]) ;  
       printf(" %12.6f %12.6f", xden[s]/zden[s]) ;  
       printnl() ;
      }
  

      yy2[s] += yy * yy ;

      y1 = res[k1]*wt[k1] ;
      y2 = res[k2]*wt[k2] ;
      
      corrpt = corrbins[s] ; 
      addcorr(corrpt, y1, y2) ;
      ++jwt[chrom][s] ;

      corrpt = corrsbins[icount][s] ; 
      addcorr(corrpt, y1, y2) ;

      for (k=1; k<= numchrom; ++k) {
        if (k==chrom) continue ;
        if (jackknife == NO) break ;
        corrpt = corrjbins[k][s] ;
        addcorr(corrpt, y1, y2) ;

        corrpt = corrsjbins[icount][k][s] ;
        addcorr(corrpt, y1, y2) ;

        xxnum[k][s] += ywt*yy ; 
        xxden[k][s] += ywt*ywt ; 
      }
   }
   printit = YES ; 
   if (printit > 0) { 
    printf("sample: %20s chrom: %3d done\n", indx -> ID, printit-1) ;  fflush(stdout) ;
   }
  }

  
   
   if (numqbins==0) printf("sample: %20s chrom: %3d done\n", indx -> ID, lastchrom) ;  fflush(stdout) ;
   y = cputime(1) ;  
   fflush(stdout) ; 


   printf("sample: %20s completed.  Total CPU since start: %9.3f\n", indx -> ID, y) ; fflush(stdout) ; 

 }

 free(valids) ; 

   numadmgood = 0 ; 
   for (icount = 0 ; icount < ncount; ++icount) { 
    k = admindices[icount] ; 
    indx = indivmarkers[k] ;
    if (indx -> ignore) continue ;  
    ++numadmgood ; 
   }


 if (numqbins>0) {  

   ddcorr(corrjbins, ddcbins,  numdbins, qbin) ;  
   fixjcorr(corrbins, corrjbins, numbins) ; // set up jackknifed corrs 

   vzero(ww, numbins) ;
   sprintf(sss, " ##Z-score and correlation:: %s  binsize: %12.6f", admixpop, binsize) ;
   dumpit(oname, ww, corrbins, numbins-5, binsize, sss) ;

   for (k=1; k<=numchrom; ++k) { 
    if (jackknife == NO) break ;
    sprintf(s1, "%s:%d", oname, k) ;
    sprintf(s2, "## Jackknife output: chrom %d", k) ;
    dumpit(s1, ww, corrjbins[k], numbins-5, binsize, s2) ;
   }

   for (icount = 0 ; icount < ncount; ++icount) { 
    if (singlesampleoutput == NO) break ; 
    k = admindices[icount] ; 
    indx = indivmarkers[k] ;
    if (indx -> ignore) continue ;  
    vzero(ww, numbins) ;
    ddcorr(corrsjbins[icount], ddscbins[icount],  numdbins, qbin) ;  
    fixjcorr(corrsbins[icount], corrsjbins[icount], numbins) ; // set up jackknifed corrs 
    sprintf(dirname, "%s", indx -> ID) ;  
    sprintf(sss, "mkdir -p -m 0755 %s\n", dirname) ;
    system(sss) ; 
    sprintf(sss, "chmod  0755 %s\n", dirname) ;
    system(sss) ; 
    chdir (dirname) ; 
    sprintf(sss, "ln -f ../%s ." , parname) ; 
    system(sss) ; 
    sprintf(sss, " ##Z-score and correlation:: %s  binsize: %12.6f :: sample %s", admixpop, binsize, indx -> ID) ;
    dumpit(oname, ww, corrsbins[icount], numbins-5, binsize, sss) ;
    for (k=1; k<=numchrom; ++k) { 
     if (jackknife == NO) break ;
     sprintf(s1, "%s:%d", oname, k) ;
     sprintf(s2, "## Jackknife output: chrom %d", k) ;
     dumpit(s1, ww, corrsjbins[icount][k], numbins-5, binsize, s2) ;
    }
    dorunfit() ;

    chdir("..") ;  
  }

 } 


  if (numqbins == 0) {

  vsp(xden, xden, 1.0e-20, numbins) ; 
  vsp(zden, zden, 1.0e-20, numbins) ; 
  vsp(yy2, yy2, 1.0e-10, numbins) ; 

   vsqrt(xden, xden, numbins) ;

   vvd(ww, xnum, xden, numbins) ;
   vvd(ww1, znum, zden, numbins) ;
   
   sprintf(sss, " ##Z-score and correlation:: %s  binsize: %12.6f", poplistname, binsize) ;
   dumpit(oname, ww, corrbins, numbins-5, binsize, sss) ;
  
   for (k=1; k<=numchrom; ++k) { 
    if (jackknife == NO) break ;
    sprintf(s1, "%s:%d", oname, k) ;
    sprintf(s2, "## Jackknife output: chrom %d", k) ;
    vvd(ww, xxnum[k], xxden[k], numbins) ;
    dumpit(s1, ww, corrjbins[k], numbins-5, binsize, s2) ;
   }

  }

  if (numadmgood == 0) fatalx("after analysis -- no valid admixed samples\n") ; 

  dorunfit() ;

   free(w0) ;
   free(w1) ;
   free(w2) ;
   free(wk) ;
   free(ww1) ;
   free(ww2) ;
   free(wa) ;
   free(wb) ;
   free(res) ;
   free(dd) ;
   free(wt) ;

   free(sindex) ;
   free(xindex) ;

   free(ww) ;

   free(xnum) ;
   free(xden) ;
   free(xdenq) ;
   free(znum) ;
   free(zden) ;
   free(yy2) ;

   free2D(&freqs, ncols) ;

    free2D(&xxnum, numchrom+1) ;
    free2D(&xxden, numchrom+1) ;

   free(firstsnp) ; 
   free(lastsnp) ; 
   free(firstqb) ; 
   free(lastqb) ; 

    if (numqbins>0) { 
     free(z0q) ; 
     free(z1q) ; 
     free(z2q) ; 
     

     for (k=1; k<=numchrom; ++k) { 
      free2D(&ddcbins[k], 7) ;
     }
     free(ddcbins) ; 
     free2D(&ddbins, 7) ; 
     for (a=0; a<ncount; ++a) {  
      for (k=1; k<=numchrom; ++k) { 
       free2D(&ddscbins[a][k], 7) ;
      }
      free(ddscbins[a]) ; 
     }
     free(ddscbins) ;
    }



   for (x=0; x<numbins; ++x) { 
     free(corrbins[x]) ;
   }

   free(corrbins) ;


   free2D(&jwt, numchrom+1) ;

   for (k=1; k<=numchrom; ++k) {

    for (x=0; x<numbins; ++x) { 
     free(corrjbins[k][x]) ; 
    }
    free(corrjbins[k]) ; 
   }

   free(corrjbins) ;


  for (a=0; a<ncount; ++a) { 
   for (x=0; x<numbins; ++x) { 
     free(corrsbins[a][x]) ;
   }

   free(corrsbins[a]) ;



   for (k=1; k<=numchrom; ++k) {

    for (x=0; x<numbins; ++x) { 
     free(corrsjbins[a][k][x]) ; 
    }
    free(corrsjbins[a][k]) ; 
   }

   free(corrsjbins[a]) ;
 }
 free(corrsbins) ;  
 free(corrsjbins) ;  

   free(bf) ; 
   free(gf) ;


  return 0;

}

void dorunfit() 
{  
  char s1[256], s2[256], s3[256], sss[256] ;

  if (runfit == NO) return ; 
   if (logfit == NULL) { 
    sprintf(s2, "%s:log", admixpop) ; 
    logfit = strdup(s2) ; 
  }

   if (jackknife) {
     sprintf(s1, "dates_jackknife -p %s -l %9.3f -z %s -r %d ", parname, lovalfit, admixpop, seed) ;
     if (snppath != NULL) { 
      sprintf(s3, "-m %s", snppath) ; 
      strcat(s1, s3) ;
     }
   }
   else  {
     sprintf(s1, "run_dates_expfit -p %s -l %9.3f -z %s -r %d ", parname, lovalfit, admixpop, seed) ;
   }
   
   if (afffit) { 
    strcat(s1, " -a") ; 
   }

   sprintf(s2, " >%s", logfit) ;
   strcpy(sss, s1) ; 
   strcat(sss, s2) ; 
   printf("calling fit!...\n") ; 
   printf("%s\n", sss) ; 
   fflush(stdout) ;
   system(sss) ;

}

void dumpit(char *oname, double *ww, CORR **corrbins, int len, double bsize, char *hdr) 
{
   double y, yreg, ycovar ; 
   int k, ret ;
   FILE *fff ;
   CORR *corrpt ;

   if (oname == NULL) return ;
   openit(oname, &fff, "w") ;

   if (hdr != NULL) fprintf(fff, "%s\n", hdr) ;
   for (k=0; k<len; ++k) { 

    y = (k+1) * bsize ;
    corrpt = corrbins[k] ;
    fprintf(fff, "%9.3f ", 100.0*y) ;  // CM
//  fprintf(fff, "%12.6f ", ww[k]) ;   
//  fprintf(fff, "%12.6f ", corrpt -> Z) ;
    if (runmode != 2) {
     calccorr(corrpt, 1, NO) ;
     ycovar = corrpt -> v12 ; 
     yreg = corrpt -> S12 / (corrpt -> S11 + 1.0e-20) ;
     if (isnan(yreg)) { 
      printf("bad corrpt!\n") ; 
      printcorr(corrpt) ; 
      fatalx("dumpit fails\n") ; 
     }
     fprintf(fff, "%15.9f ", ycovar) ; 
     fprintf(fff, "%12.6f ", yreg) ; 
     fprintf(fff, "%12.6f ", corrpt -> corr) ;   
     fprintf(fff, "%12.0f ", corrpt -> S0) ;
     if (corrpt -> S0 > 1.0e20) { 
      printf("bad corr!\n") ; 
      printcorr(corrpt) ; 
      fatalx("bad corr!\n") ; 
     }
     if (fabs(yreg) > 10.0)  { 
      printf("bad corr!\n") ; 
      printcorr(corrpt) ; 
      fatalx("bad yreg!\n") ; 
     }
    }
    fprintf(fff, "\n") ;

   }


   fclose(fff) ;

}
void
dorcx(double *ans, double *rr, double *dd, double binsize, double maxdis) 
// accumulate dot product into bins
{
   int x, x2, z, t, s, k ; 
   SNP *cupt, *cupt2 ;
   double dis ;
   double *xnum, *xd1, *xd2, y1, y2 ;
   double *rrr ;

   z =nnint(maxdis/binsize) ;  
   ZALLOC(xnum, z+1, double) ;
   ZALLOC(xd1, z+1, double) ;
   ZALLOC(xd2, z+1, double) ;

   vclear(xd1, 1.0e-10, z+1) ;
   vclear(xd2, 1.0e-10, z+1) ;

   for (x=0; x<numsnps; ++x)  {
    cupt = snpmarkers[x] ; 
    for (x2=x+1; x2<numsnps; ++x2)  {
      cupt2 = snpmarkers[x2] ; 
      t = cupt2 -> chrom - cupt -> chrom ;
      if (t != 0) break ;
      dis = cupt2 -> genpos - cupt -> genpos ;
      if (dis >= maxdis) break ;
      s = nnint(dis/binsize) ;
      y2 = dd[x]*dd[x2] ;
      for (k=0; k<numindivs; ++k) { 
       rrr = rr + k*numsnps ;

       switch (runmode) { 
        case 0:  
         y2 = dd[x]*dd[x2] ;
         y1 = rrr[x]*rrr[x2] ;
         break ;
      
        case 1:  
         y2 = rrr[x2]*dd[x2] ;
         y1 = rrr[x]*dd[x]  ;
         break ;

/**
        case 2:  
         y2  = 1.0 ;
         y1 = rrr[x]*rrr[x2]  ;
         break ;

        case 3:  
         y2  = 1.0 ;
         y1 = dd[x]*dd[x2] ;
         break ;
*/

        default: 
         fatalx("bad runmode\n") ; 
       }

       if (isnan(y1)) fatalx("bad y1\n") ;
       if (isnan(y2)) fatalx("bad y2\n") ;
       xnum[s] += y1*y2 ;
       xd1[s] += y1*y1 ;
       xd2[s] += y2*y2 ;
      }
    }
   }
  
   
   for (s=0; s <= z; ++s) {
    ans[s] = xnum[s]/sqrt(xd1[s]*xd2[s]) ;
   }

   free(xnum) ;
   free(xd1) ;
   free(xd2) ;
}

void
dorc2(double *ans, double *res, double binsize, double maxdis) 
// accumulate square of dot product into bins
{
   int x, x2, z, t ; 
   SNP *cupt, *cupt2 ;
   double dis, y ;

   z =nnint(maxdis/binsize) ;  
   vzero(ans, z+1) ;
   for (x=0; x<numsnps; ++x)  {
    cupt = snpmarkers[x] ; 
    for (x2=x+1; x2<numsnps; ++x2)  {
      cupt2 = snpmarkers[x2] ; 
      t = cupt2 -> chrom - cupt -> chrom ;
      if (t != 0) break ;
      dis = cupt2 -> genpos - cupt -> genpos ;
      if (dis >= maxdis) break ;
      z = nnint(dis/binsize) ;
      y  = res[x]*res[x2] ;
      ans[z] += y*y ;             
    }
   }
}

int readcommands(int argc, char **argv) 

{
  int i,haploid=0;
  phandle *ph ;
  char str[5000], *sx  ;
  char *tempname ;
  char *spt[MAXFF] ; 
  int n, k ;

  while ((i = getopt (argc, argv, "p:vV")) != -1) {

    switch (i)
      {

      case 'p':
	parname = strdup(optarg) ;
	break;

      case 'v':
	printf("version: %s\n", WVERSION) ; 
	break; 

      case 'V':
	verbose = YES ;
	break; 

      case '?':
	printf ("Usage: bad params.... \n") ;
	fatalx("bad params\n") ;
      }
  }

         
   if (parname == NULL) return -1 ;
   printf("parameter file: %s\n", parname) ;
   ph = openpars(parname) ;
   dostrsub(ph) ;

   getstring(ph, "genotypename:", &genotypename) ;
   getstring(ph, "snpname:", &snpname) ;
   getstring(ph, "indivname:", &indivname) ;
   getstring(ph, "poplistname:", &poplistname) ;
   getstring(ph, "weightname:", &weightname) ;
   getstring(ph, "timeoffset:", &timeoffsetname) ;
   getstring(ph, "admixpop:", &admixpop) ;
   getstring(ph, "admixlist:", &admixlist) ;
   getstring(ph, "badsnpname:", &badsnpname) ;
   getstring(ph, "dumpname:", &dumpname) ;
   getstring(ph, "loadname:", &loadname) ;
   getint(ph, "popestimate:", &popestimate) ;
   getint(ph, "ransample:", &ransample) ;
   getint(ph, "chrom:", &xchrom) ;
   getint(ph, "nochrom:", &xnochrom) ;
   getstring(ph, "output:",   &oname) ;
   getint(ph, "runmode:", &runmode) ;
   getint(ph, "ldmode:", &ldmode) ;
   getint(ph, "seed:", &seed) ;
   getint(ph, "nxlim:", &nxlim) ;
   getint(ph, "minparentcount:", &minparentcount) ;
   getint(ph, "norun:", &norun) ;
   getint(ph, "flatmode:", &flatmode) ;
   getint(ph, "zdipcorrmode:", &zdipcorrmode) ;
   getint(ph, "jackknife:", &jackknife) ;
   getint(ph, "usevar:", &usevar) ;
   getint(ph, "regdim:", &regdim) ;
   getint(ph, "selfsource:", &selfsource) ;
   getint(ph, "samecoeffs:", &samecoeffs) ;
   getint(ph, "singlesampleoutput:", &singlesampleoutput) ;
   getint(ph, "whichadm:", &whichadm) ;
   getint(ph, "numchrom:", &numchrom) ;
   getint(ph, "usehets:", &usehets) ;
   

   getdbl(ph, "maxdis:", &maxdis) ;
   getdbl(ph, "binsize:", &binsize) ;
   getdbl(ph, "chithresh:", &chithresh) ;
   getint(ph, "checkmap:", &checkmap) ;
   getint(ph, "qbin:", &qbin) ;
// parameters for runfit etc. 
   getint(ph, "runfit:", &runfit) ;
   getint(ph, "afffit:", &afffit) ;
   getdbl(ph, "lovalfit:", &lovalfit) ;

   getstring(ph, "specfile:", &specfile) ;
   sx = NULL ;
   getstring(ph, "speccols:", &sx) ;
   if (sx != NULL) { 
    numspec = splitupx(sx, spt, MAXFF, ':') ; 
    ZALLOC(speccols, numspec, int) ;
     
    for (k=0; k<numspec; ++k) {  
     speccols[k] = atoi(spt[k]) ; 
    }
    freestring(&sx) ; 
    freeup(spt, numspec) ;
   }
   getstring(ph, "specnames:", &specstring) ;


   writepars(ph) ;
   closepars(ph) ;

}

int      
setqbins(SNP **snpm, int numsnps, double binsize, int qbin) 
// return number of qbins ;  set tagnumber to qbin 
{
 double qb, ydis, lastgenpos ; 
 int chrom, s = 0, k ; 
 SNP *cupt ;

 if (qbin==0) return 0 ; 

 qb = binsize/(double) qbin ; 

 for (k=0; k<numsnps; ++k) { 
  cupt = snpm[k] ; 
  if (k==0) { 
   ydis = 0 ; 
   lastgenpos = cupt -> genpos ;
   s = floor(ydis/qb) ;  
   cupt -> tagnumber = s ; 
   chrom = cupt -> chrom ; 
   continue ; 
  }
  if (cupt -> chrom != chrom) { 
   ydis += 5 ; 
   chrom = cupt -> chrom ; 
  }
  else { 
   ydis += (cupt -> genpos - lastgenpos) ;
  }
  lastgenpos = cupt -> genpos ;
  s = floor(ydis/qb) ;  
  cupt -> tagnumber = s ;
 }

 return s+1 ;



}
void setfc(SNP **snpm, int numsnps) 
// also allocates fchrom , lchrom
{
  int i, j, chrom ;
  SNP *cupt ;
  double dis, totdis, pos1, pos2 ;

  ZALLOC(fchrom, numchrom+3, int) ;
  ZALLOC(lchrom, numchrom+3, int) ;

  ivclear(fchrom, 999999, numchrom+3) ;
  ivclear(lchrom, -999999, numchrom+3) ;
  
/* initialize real marker array */
  for (i=0; i<numsnps; i++) {
    cupt = snpm[i] ;
    if (cupt -> ignore) continue ;
    chrom = cupt -> chrom ;
//  if ((noxdata) && (chrom==23)) continue  ;
    fchrom[chrom] = MIN(fchrom[chrom],i) ;
    lchrom[chrom] = MAX(lchrom[chrom],i) ;
  }
}
double cntit1(double *xc, SNP *cupt, Indiv **indm, int numindivs, int t) 
{
  Indiv *indx ;
  int k, g ;

  vzero(xc, 3) ;
  if (cupt -> ignore) return -1 ;
  for (k=0; k<numindivs; ++k) { 
   indx = indm [k] ;
   if (indx -> ignore) continue ;
   if (indx -> affstatus  != t) continue ;
   g = getgtypes(cupt, k) ;
   if (g<0) continue ;
   ++xc[g] ;
  }
  return asum(xc,3) ;
}


double cntit2(double *xc, SNP *cupt, SNP *cupt2, Indiv **indm, int numindivs, int t) 
{
  Indiv *indx ;
  int k, e, f ;

  vzero(xc, 9) ;
  if (cupt -> ignore) return -1 ;
  if (cupt2 -> ignore) return -1 ;
  for (k=0; k<numindivs; ++k) { 
   indx = indm [k] ;
   if (indx -> ignore) continue ;
   if (indx -> affstatus  != t) continue ;
   e = getgtypes(cupt, k) ;
   if (e<0) continue ;
   f = getgtypes(cupt2, k) ;
   if (f<0) continue ;
   ++xc[3*e+f] ;
  }
  return asum(xc,9) ;

}

void cx1(double *xc, double *y1, double *y2) 
{
   double x1[3], x2[3] ;
   int e, f ;

   vzero(x1,3) ;
   vzero(x2,3) ;
   for (e=0; e<3; ++e) {  
    for (f=0; f<3; ++f) {  
     x1[e] += xc[3*e+f] ;
     x2[f] += xc[3*e+f] ;
    }
   }
  
   *y1 = cminor(x1) ;
   *y2 = cminor(x2) ;
    
}

double cminor(double *xx) 
{
  double y1, y2 ;
  y1 = xx[1]+2*xx[2] ; 
  y2 = 2*asum(xx, 3) - y1 ;
  return MIN(y1, y2) ;
}
void
printzq(double *z0q, double *z1q, double *z2q, int numqbins)
{
  int s, t ;

  printf("##printzq %9.0f\n", asum(z0q, numqbins) ) ;

  for (s=0; s<numqbins; ++s)  { 
   t = nnint(z0q[s]) ;  
   if (t==0) continue ;  
   printf("zq %9d %6.0f %9.3f %9.3f\n", s, z0q[s], z1q[s], z2q[s]) ;
 }

}


void
ddadd(double ***ddcbins, double *z0q, double *z1q, double *z2q, int numqbins, int diffmax) 

// this is the only potentially quadratic part of the code.   
{
  double *dd00, *dd01, *dd10, *dd11, *dd02, *dd20 ; 
  double **ddbins ;
  double a0, a1, a2, b0, b1, b2, y ; 
  int s, t, tmax, d, js, jt, k, slo, shi, len, j ; 
  SNP *cupts, *cuptt ;

  double *dans ;  
  double *dans2 ;  

  if (numqbins == 0) return ; 
  ZALLOC(dans, diffmax +10, double) ;
  ZALLOC(dans2, diffmax +10, double) ;
 

 for (k=1; k<=numchrom; ++k) { 

  slo = firstqb[k] ;
  shi = lastqb[k] ;

  ddbins = ddcbins[k] ;

  dd00 = ddbins[0] ; 
  dd01 = ddbins[1] ; 
  dd10 = ddbins[3] ; 
  dd11 = ddbins[4] ; 
  dd02 = ddbins[2] ; 
  dd20 = ddbins[6] ; 

  len = shi-slo +1 ;
  tmax = MIN(diffmax+1, len+1) ; 

  if (len <= 1) continue ;

  s = slo ; 
  
  fftauto(dans, z0q +s, len, diffmax) ; 
  vvp(dd00, dd00, dans, tmax) ;  

  for (j=0; j<=diffmax; ++j) { 
    y = fabs(dans[j]) ;
    if (isnan(y)) fatalx("bad ddadd\n") ;
    if (y> 1.0e20) { 
     printmat(z0q+s, 1, len) ; 
     fatalx("bad ddadd\n") ;
    }
  }

  fftcv2(dans, dans2, z0q +s, z1q +s, len, diffmax) ; 

  vvp(dd01, dd01, dans, tmax) ; 
  vvp(dd10, dd10, dans2, tmax) ; 

  for (j=0; j<=diffmax; ++j) { 
    if (isnan(dans[j])) fatalx("bad ddadd\n") ;
    if (isnan(dans2[j])) fatalx("bad ddadd\n") ;
  }

  fftauto(dans, z1q +s, len, diffmax) ; 
  vvp(dd11, dd11, dans, tmax) ; 

  fftcv2(dans, dans2, z0q +s, z2q +s, len, diffmax) ; 
  vvp(dd02, dd02, dans, tmax) ; 
  vvp(dd20, dd20, dans2, tmax) ; 
 }
 
 free(dans) ;
 free(dans2) ;

}

void ddcorr(CORR ***corrjbins, double ***ddcbins, int numdbins, int qbin)  
{
  double *dd00, *dd01, *dd10, *dd11, *dd02, *dd20 ; 
  int s, d, k ; 
  CORR *corrpt, **corrtemp ; 
  double ys, dbinsize ; 
  double **ddbins, **ddtemp ;   

  if (numdbins==0) return ; 

  ddbins = initarray_2Ddouble(7, numdbins, 0) ;

  for (k=1; k<=numchrom; ++k) { 
   corrtemp = corrjbins[k] ; 
   ddtemp = ddcbins[k] ; 

  dd00 = ddtemp[0] ; 
  dd01 = ddtemp[1] ; 
  dd10 = ddtemp[3] ; 
  dd11 = ddtemp[4] ; 
  dd02 = ddtemp[2] ; 
  dd20 = ddtemp[6] ; 

  dbinsize = binsize / (double) qbin ; 
  for (d=1; d<numdbins; ++d) {   // dangerous bend,  lag 0 should be 0 and is not computed correctly!
   if (dd00[d] < 0.5) continue ; 
   ys = (double) d  * dbinsize ;     
   s = (int) (ys/binsize) ; 
   corrpt = corrtemp[s] ; 
   addcorr2(corrpt, dd00[d], dd01[d], dd10[d], dd11[d], dd02[d], dd20[d]) ; 
  }
 
 }

}
void fixjcorr(CORR **corrbins, CORR ***corrjbins, int numbins) 
{
  int k, x ; 
  CORR *corrpt, *corrptj ; 
 
  for (x=0; x<numbins; ++x) { 
// step 1 calculate global corr
    corrpt = corrbins[x] ;  
    clearcorr(corrpt) ; // not really necessary 
    for (k=1; k<=numchrom; ++k) { 
     corrptj = corrjbins[k][x] ; 
     pluscorr(corrpt, corrpt, corrptj) ;         
    }
// step 2 calculate jackknifed correlations by subtraction
    for (k=1; k<=numchrom; ++k)  {  
     minuscorr(corrjbins[k][x], corrpt, corrjbins[k][x]) ;
    }
  }

}
void getsn(char **specnames, int numspec, char *specstring) 
{
  char *sx  ;
  int nsplit = 0 ;
  char *spt[MAXFF] ;

  if (specstring == NULL) fatalx("no specnames:\n") ; 
  nsplit = splitupx(specstring, spt, MAXFF, ':') ;
  if (nsplit != numspec) fatalx("bad specnames: %s\n", specstring) ;
  copystrings(spt, specnames, nsplit) ; 
}

int  solvitc(double *co, double *rr, int n, double *ans, double *ccc, double ccr)   
// taken from regressitc 
{
 int ret ;
 double  lambda, y1, y2 ;
 double *a1, *a2 ; 
 
 vzero(ans, n) ;
 ZALLOC(a1, n, double) ;
 ZALLOC(a2, n, double) ;

 ret = solvit(co, rr, n, a1) ;
 if (ret < 0) { 
  free(a1) ; free(a2) ; return -1 ; 
 }
 ret = solvit(co, ccc, n, a2) ;
 if (ret < 0) { 
  free(a1) ; free(a2) ; return -2 ; 
 }

// ans = a1 + lambda * a2
 y1 = vdot(a1, ccc, n) ;
 y2 = vdot(a2, ccc, n) ; // y2 > 0 if co is pos. def.
// ccr = y1 + lambda * y2

   lambda = (ccr-y1)/y2 ;
   vst(a2, a2, lambda, n) ;
   vvp(ans, a1, a2, n) ;

   free(a1) ; free(a2) ; 
   return 1 ; 

}

double calccoeffs3(double **pcoeffs, SNP **snpmarkers, int numsnps, Indiv **indivmarkers, 
 int numindivs, int *admindexes, int ncount) 
{

   double prod[4*4], prhs[4], ans[4] ;
   double tprod[4*4], trhs[4] ;
   double *coeffs, *ctrans, *rhs, *cx ; 
   int g, i, j, k, bcount[3], t, a  ;  
   double y1, y2, yscore, y ; 
   double *b1, *b2, *b3, *bkon ; 
   int *valids ; 
   SNP *cupt ;
   double ccc[4] ; 
   static int ncall = 0, ngood, ret ; 
   int admindex ; 
   Indiv *indx ; 
   double *w0, *w1, *w2, *ww1, *ww2 ;  
   double *xj, *xk ; 

   ++ncall ; 


   ZALLOC(coeffs, numsnps*4, double) ;
   ZALLOC(rhs, numsnps, double) ;
   ZALLOC(valids, numsnps, int) ; 
   
   cx = coeffs ; 
   b1 = cx ;     cx += numsnps ; 
   b2 = cx ;     cx += numsnps ; 
   b3 = cx ;       cx += numsnps ; 
   bkon = cx ;      
   
  vzero(tprod, 4*4) ;
  vzero(trhs, 4) ;

  for (a=0; a<ncount; ++a) { 
   
   vzero(pcoeffs[a], regdim) ; 
   admindex = admindexes[a] ;  
   indx = indivmarkers[admindex] ;
   ngood = 0 ;  
   ivclear(valids, YES, numsnps) ; 



   vzero(rhs, numsnps) ; 
   vzero(coeffs, numsnps*4) ; 
   
   for (j=0; j<numsnps; j++)  { 
     cupt = snpmarkers[j] ; 
     if (cupt -> ignore) { 
      valids[j] = NO ; 
      continue ;  
     }

     g = getgtypes(cupt, admindex) ; 
     if (g<0) { 
      valids[j] = NO ; 
      continue ;  
     }

     b1[j] = cupt -> cf[0]  ;  
     b2[j] = cupt -> cf[1] ; 
     b3[j] = cupt -> cf[2] ; 
  
     rhs[j] = (double) g / 2.0 ; 

     if (regdim==3) continue ; 

     copyiarr(cupt -> bcount, bcount, 3) ; 

   }

   t = intsum(valids, numsnps) ;
   if (t<100) { 
    printf("too little data for: %s number of valid snps %d\n", indx -> ID, t) ;
    indx -> ignore = YES ; 
    continue ;
   }

 vzero(prod, 4*4) ; 
 vzero(prhs, 4) ; 

 for (i=0; i<numsnps ; i++) {
  if (valids[i] == NO) {         
   b1[i] = b2[i] = b3[i] =  bkon[i] =  rhs[i] = 0 ; 
  }
 }
  for (j=0; j<regdim ; j++) {
   xj = coeffs + j*numsnps ; 
   prhs[j] = vdot(xj, rhs, numsnps) ;      
   for (k=j; k<regdim ; k++) {
    xk = coeffs + k*numsnps; 
    prod[j*regdim+k] = prod[k*regdim+j] =  vdot(xj, xk, numsnps) ;
   }
  }
   vvp(tprod, tprod, prod, 4*4) ; 
   vvp(trhs, trhs, prhs, 4) ; 

   vclear(ccc, 1, regdim) ; 
   if (samecoeffs) continue ; 
   ret = solvitc(prod, prhs, regdim, ans, ccc, 1) ;  
// printf("zzans: %d %s ", ret, indx->ID) ; printmat(ans, 1, regdim) ;
   printf("coeffs: %s", indx->ID) ; printmat(ans, 1, regdim) ;
   copyarr(ans, pcoeffs[a], regdim) ;

   if (regdim==2) { 
    fatalx("very bad bug\n") ;
   }

  }

   if (samecoeffs) { 
     vclear(ccc, 1, regdim) ; 
     ret = solvitc(tprod, trhs, regdim, ans, ccc, 1) ;  
     printf("coeffs(samecoeffs set): global") ; printmat(ans, 1, regdim) ;
     for (a=0; a<ncount; ++a) { 
      copyarr(ans, pcoeffs[a], regdim) ;
     }
   }
   free(coeffs) ; 
   free(rhs) ; 
   free(valids) ;
   
}


double calccoeffs(double **pcoeffs, SNP **snpmarkers, int numsnps, Indiv **indivmarkers, 
 int numindivs, int *admindexes, int ncount, int numeg) 
{

   double prod[4*4], prhs[4], ans[4] ;
   double tprod[4*4], trhs[4] ;
   double *coeffs, *ctrans, *rhs, *cx ; 
   int g, i, j, k, bcount[3], t, a  ;  
   double y1, y2, yscore, y ; 
   double *bleft, *bright, *bkon, *bpred ;
   int *valids ; 
   SNP *cupt ;
   double ccc[4] ; 
   static int ncall = 0, ngood, ret ; 
   int admindex ; 
   Indiv *indx ; 
   double *w0, *w1, *w2, *ww1, *ww2 ;  
   double *xj, *xk ; 

   ++ncall ; 
   if (numeg==3) return calccoeffs3(pcoeffs, snpmarkers, numsnps, indivmarkers, numindivs, admindexes, ncount) ;

   if (regdim==2) { 
    ZALLOC(ww1, numsnps, double) ;
    ZALLOC(ww2, numsnps, double) ;
   }

   ZALLOC(coeffs, numsnps*4, double) ;
   ZALLOC(rhs, numsnps, double) ;
   ZALLOC(valids, numsnps, int) ; 
   
   cx = coeffs ; 
   bleft = cx ;     cx += numsnps ; 
   bright = cx ;     cx += numsnps ; 
   bkon = cx ;       cx += numsnps ; 
   bpred = cx ;      
   
  vzero(tprod, 4*4) ;
  vzero(trhs, 4) ;

  for (a=0; a<ncount; ++a) { 
   
   vzero(pcoeffs[a], regdim) ; 
   admindex = admindexes[a] ;  
   indx = indivmarkers[admindex] ;
   ngood = 0 ;  
   ivclear(valids, YES, numsnps) ; 



   vzero(rhs, numsnps) ; 
   vzero(coeffs, numsnps*4) ; 
   
   for (j=0; j<numsnps; j++)  { 
     cupt = snpmarkers[j] ; 
     if (cupt -> ignore) { 
      valids[j] = NO ; 
      continue ;  
     }

     g = getgtypes(cupt, admindex) ; 
     if (g<0) { 
      valids[j] = NO ; 
      continue ;  
     }

     bleft[j]  = cupt -> cf[0]  ;  
     bright[j] = cupt -> cf[1] ; 
     rhs[j] = (double) g / 2.0 ; 

     if (regdim==2) continue ; 
     bkon[j] = 0.5 ; 
     if (regdim==3) continue ; 

     copyiarr(cupt -> bcount, bcount, 3) ; 

     --bcount[g] ; 
     if (bcount[g] < 0) fatalx("(calccoeffs) badbug!\n") ; 

     y1 = bcount[1] + 2*bcount[2] ;
     y2 = 2*intsum(bcount, 3) ; 
     if (y2<0.1) { 
      valids[j] = NO ; 
      continue ;  
     }
     bpred[j] = y1/y2 ; 
     ++ngood ; 
   }  
   t = intsum(valids, numsnps) ;
   if (t<100) { 
    printf("too little data for: %s number of valid snps %d\n", indx -> ID, t) ;
    indx -> ignore = YES ; 
    continue ;
   }

 vzero(prod, 4*4) ; 
 vzero(prhs, 4) ; 

 for (i=0; i<numsnps ; i++) {
  if (valids[i] == NO) {         
   bleft[i] = bright[i] = bkon[i] = bpred[i] = rhs[i] = 0 ; 
  }
 }
  for (j=0; j<regdim ; j++) {
   xj = coeffs + j*numsnps ; 
   prhs[j] = vdot(xj, rhs, numsnps) ;      
   for (k=j; k<regdim ; k++) {
    xk = coeffs + k*numsnps; 
    prod[j*regdim+k] = prod[k*regdim+j] =  vdot(xj, xk, numsnps) ;
   }
  }
   vvp(tprod, tprod, prod, 4*4) ; 
   vvp(trhs, trhs, prhs, 4) ; 

   vclear(ccc, 1, regdim) ; 
   if (samecoeffs) continue ; 
   ret = solvitc(prod, prhs, regdim, ans, ccc, 1) ;  
   printf("coeffs:  %s ", indx->ID) ; printmat(ans, 1, regdim) ;

   if (regdim==2) { 
    vvm(ww1, rhs, bright, numsnps) ; // target - right
    vvm(ww2, bleft, bright, numsnps) ; // left - right 
    y = vdot(ww1, ww2, numsnps) / vdot(ww2, ww2, numsnps) ;
    printf("coeff: %3d %20s %9.3f\n", i, indx -> ID, y) ;
   }
   copyarr(ans, pcoeffs[a], regdim) ;
  }

   if (samecoeffs) { 
     vclear(ccc, 1, regdim) ; 
     ret = solvitc(tprod, trhs, regdim, ans, ccc, 1) ;  
     printf("(samecoeffs): %s ", "global") ; printmat(ans, 1, regdim) ;
     for (a=0; a<ncount; ++a) { 
      copyarr(ans, pcoeffs[a], regdim) ;
     }
   }
   free(coeffs) ; 
   free(rhs) ; 
   free(valids) ;
   
}

double xcval(double yy[2], double *coef, double *zzz, int g, int which) 
{

  double zz[3] ;
  double y, xa, xb, xc, zden, yl, aa, bb, cc, bbb, ccc, xbb, q ; 
  int numeg = 2 ; 
  double yy0[2], yy1[2] ; 

  vclear(yy, -999, 2) ;

  if (g<0) return -9999 ;   

  if (which>0) numeg = 3 ; 

  copyarr(zzz, zz, 3) ; 

  if (g==1) { 
   yl = xcval(yy0,  coef, zzz, 0, which) ;
   yl += xcval(yy1,  coef, zzz, 2, which) ;
   yl /= 2.0 ; 
   vvp(yy, yy0, yy1, 2) ; 
   vst(yy, yy, 0.5, 2) ;
   return yl ; 
  }

  if (g==0) {
    vcompl(zz, zz, 3) ; 
  }

  zden = vdot(coef, zz, numeg) ;
  yl = log(2.0*zden) ; 
  if (!isfinite(yl)) fatalx("badbug (xcval)\n") ; 

 xa = zz[0] ; 
 xb = zz[1] ; 
 xc = zz[2] ; 

  aa = coef[0] ; 
  bb = coef[1] ; 
  cc = coef[2] ; 

     if (numeg == 2) {

      q = sqrt(aa * bb) ; 
      yy[0] = q*xa/zden ; 
      yy[1] = q*xb/zden ; 

      return yl ; 

   }
 

     if (which==1) {    // near adm event
       bbb = bb + cc ;
       q = sqrt(aa*bbb) ;
       xbb = (bb*xb + cc*xc ) / bbb ;
       yy[0] = q*xa/zden ;
       yy[1] = q*xbb/zden ;

       return yl ; 
     }

     if (which==2) {    // far adm event
       y = bb*cc/(bb+cc) ;
       q = sqrt(y) ;
       yy[0] = q*xb/zden ;
       yy[1] = q*xc/zden ;
       return yl ; 
     }





  return yl ; 
}

int calcab(double *ya, double *yb, int *valids,  SNP **snpmarkers, int numsnps, 
  Indiv ** indivmarkers, int numindivs, int admindex, double *coeffs, int numeg) 
{
  double alpha, beta ;         
  double xa, xb, xc, zden, pa, pb, a1, b1, c1, aa, bb, cc, bbb, ccc, xbb ; 
  double ynuma, ynumb, yden ; 
  double zz[3], coef[3] ; 
  double ylike = 0.0, y, ymin, q ;
  Indiv *indx ; 
  SNP *cupt ; 
  int j, g, t ; 
  long ncall = 0, which  ; 
  double yy[2] ;

  vzero(ya, numsnps) ; 
  vzero(yb, numsnps) ; 
  ivclear(valids, YES, numsnps) ; 
  alpha = coeffs[0] ; 
  beta = 1 - alpha ; 
  vzero(zz, 3) ; 

  ++ncall ; 
  indx = indivmarkers[admindex] ;

  vmaxmin(coeffs, numeg, NULL, &ymin) ;
  if (ymin<=0) {
     printf("*** warning: bad coeffs:: %9.3f ***\n", ymin) ;
     return -1 ; 
  } 

  vzero(coef, 3) ;
  copyarr(coeffs, coef, numeg) ; 
  bal1(coef, numeg) ;
  
  aa = coef[0] ; 
  bb = coef[1] ; 
  cc = coef[2] ; 

  which = whichadm ; 

  if (numeg == 2) which = 0 ;

   for (j=0; j<numsnps; j++)  { 

     cupt = snpmarkers[j] ; 
     if (cupt -> ignore) { 
      valids[j] = NO ; 
      continue ;  
     }

     copyarr(cupt -> cf, zz, numeg) ; 
     g = getgtypes(cupt, admindex) ; 
     if ((usehets == NO) && (g == 1)) g = -1 ;  // 
     if (g<0) { 
      valids[j] = NO ; 
      continue ;  
     } 
     ylike += xcval(yy, coef, zz, g, which) ;
     ya[j] = yy[0] ; 
     yb[j] = yy[1] ; 
  }


  printf(" calcab: %20s :: %9.3f\n",  indx -> ID, ylike) ; 
  return intsum(valids, numsnps) ;

}


int calcpred(double *target, double *pred, int *valids, SNP **snpmarkers, int numsnps, 
 Indiv ** indivmarkers, int numindivs, int admindex, double *zans)
{
   double *coeffs, *ctrans, *rhs, *cx, ans[4] ; 
   int g, j, bcount[3], t, k  ;  
   double y1, y2, yscore, y ; 
   double *bleft, *bright, *bkon, *bpred ;
   SNP *cupt ;
   double ccc[4] ; 
   static int ncall = 0, ngood ; 

   ++ncall ; 

   vzero(pred, numsnps) ; 
   vzero(target, numsnps) ; 
   ivclear(valids, YES, numsnps) ; 

   ZALLOC(coeffs, numsnps*4, double) ;
   ZALLOC(ctrans, numsnps*4, double) ;
   ZALLOC(rhs, numsnps, double) ;
   
   cx = coeffs ; 
   bleft = cx ;     cx += numsnps ; 
   bright = cx ;     cx += numsnps ; 
   bkon = cx ;       cx += numsnps ; 
   bpred = cx ;      
   
   ngood = 0 ;  
   for (j=0; j<numsnps; j++)  { 
     cupt = snpmarkers[j] ; 
     if (cupt -> ignore) { 
      valids[j] = NO ; 
      continue ;  
     }
     bleft[j]  = cupt -> cf[0]  ;  
     bright[j] = cupt -> cf[1] ; 
            
     bkon[j] = 0.5 ; 
     copyiarr(cupt -> bcount, bcount, 3) ; 
     g = getgtypes(cupt, admindex) ; 
     if (g<0) { 
      valids[j] = NO ; 
      continue ;  
     }
    if (regdim==4) { 
     --bcount[g] ; 
     if (bcount[g] < 0) fatalx("(calcpred) badbug!\n") ; 
     y1 = bcount[1] + 2*bcount[2] ;
     y2 = 2*intsum(bcount, 3) ; 
     if (y2<0.1) { 
      valids[j] = NO ; 
      continue ;  
     }
     bpred[j] = y1/y2 ; 
    }
     rhs[j] = (double) g / 2.0 ; 
     ++ngood ; 
   }  

   transpose(ctrans, coeffs, 4, numsnps) ;    
   for (j=0; j<numsnps; ++j) {  
     if (valids[j] == NO) { 
      vzero(ctrans + j*4, 4) ; 
      rhs[j] = 0 ; 
     } 
   } 
   if (zans == NULL) { 
    vclear(ccc, 1, 4) ; 
    yscore = regressitc(ans, ctrans, rhs, numsnps, 4, ccc, 1) ; 
    printf("regress: %20s ", indivmarkers[admindex] -> ID) ;  
    y = asum(ans, 4) ; 
    printmatx(ans, 1, 4) ;  printf(" ::  %9.3f", y) ;  printnl() ; 
    if (y<0.1) fatalx("bad regress\n") ; 
   }
   else { 
    copyarr(zans, ans, 4) ; 
   }

// impose constraint 
   for (j=0; j<numsnps; j++) { 
    if (valids[j] == 0) continue ; 
    y = vdot(ctrans+j*4, ans, regdim) ; 
    pred[j] = clip(y, 0, 1) ;
    target[j] = rhs[j] ;
   }
 
   free(coeffs) ; 
   free(ctrans) ; 
   free(rhs) ; 
   
   return intsum(valids, numsnps) ;

}
