/**************************************************************************************
    
    The biggest thing to remember here is that the matrices are all COLUMN-WISE ordered.
    In other words, they fill the first column first then the second etc.  This is
    the fortran way of doing it.

    Throughout this file I use the folowing indexers:
    I -- indexes the observations/records.
    J -- indexes the questions/variables.
    K -- indexes the states/pure types.
    L -- indexes the answers/responses.
**************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <pthread.h>
#include <sched.h>

#define BIG_NUM 1e25F
#define REAL_ZERO 1e-25F

int q_intcmp(const void *a, const void *b) {
  return(*((int*)a) - *((int*)b));
}

/* Warning!  v is replaced with the coded version of v! */
int int_coder(int *v, int N) {
  int I=1,J=0,n,*sv;

  sv = Calloc(N, int);

  memcpy(sv,v,N*sizeof(int));

  qsort(sv,(size_t)N,sizeof(int),q_intcmp);

  for(n=1;n<N;n++) {
    if(sv[n-1] != sv[n]) {
      if(sv[I-1] < 1) J++; /* do not code if < 1 */
      sv[I++] = sv[n];
    }
  }
  if(J) {
    memmove(sv,sv + J, I*sizeof(int));
    I -= J;
  }

  for(n=0;n<N;n++) /* Factors are indexed from 1 */
    if(v[n] > 0) 
      v[n] = 1 + ((int *)bsearch(v + n,sv,(size_t)I,sizeof(int),q_intcmp)) - sv;

  Free(sv);

  return(I);
}
/* Warning!  v is replaced with the coded version of v! */
int int_coder2(int *v, int N, int *sv) {
  int I=1,J=0,n;

  memcpy(sv,v,N*sizeof(int));

  qsort(sv,(size_t)N,sizeof(int),q_intcmp);

  for(n=1;n<N;n++) {
    if(sv[n-1] != sv[n]) {
      if(sv[I-1] < 1) J++; /* do not code if < 1 */
      sv[I++] = sv[n];
    }
  }
  if(J) {
    memmove(sv,sv + J, I*sizeof(int));
    I -= J;
  }

  for(n=0;n<N;n++) /* Factors are indexed from 1 */
    if(v[n] > 0) 
      v[n] = 1 + ((int *)bsearch(v + n,sv,(size_t)I,sizeof(int),q_intcmp)) - sv;

  return(I);
}

/* Put this where it needs to go */
//        R_CheckUserInterrupt();

/* get the list element named str, or return NULL */
SEXP getListElement(SEXP list, char *str) {
   SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
   int i;

   for (i = 0; i < length(list); i++)
       if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
           elmt = VECTOR_ELT(list, i);
           break;
       }
   return elmt;
}
SEXP mkSequence(double start, double end, double by) {
    unsigned long n,N = (unsigned long) fabs((end - start)/by) + 1;
    SEXP R_v;

    PROTECT(R_v = allocVector(REALSXP,N));
    for(n=0;n<N;n++) {
        REAL(R_v)[n] = start + n * by;
    }
    UNPROTECT(1);
    return(R_v);
}

double total_entropy(int *s, unsigned int I, unsigned int sK) {
  unsigned int i, k, *s_cnts;
  double total_entropy=0;
  
  s_cnts = Calloc(sK, unsigned int);
  for(i=0;i<I;i++) {
    if(s[i] > 0) s_cnts[s[i]-1]++;
  }
  for(k=0;k<sK;k++) {
    total_entropy += s_cnts[k]/(double)I * log(I/(double)s_cnts[k]);
  }
  Free(s_cnts);
  return(total_entropy);
}

/* Calculates the marginal entropy for each level of one factor.
 *  v3 -- uses less memory than v2 and is faster than v1
 */
int minfo_v3(int *s, int *v, int *w, unsigned int I, unsigned int sK, unsigned int vK, unsigned int wK, double *entropy) {
  int *vw, *svw;
  int i,sk,vk,wk,vwk;
  unsigned int *vw_cnt, *svw_cnt, *w_cnt, *sw_cnt,vwK;

  /* optimized vwk */
   vw = Calloc(I, int);
  svw = Calloc(I, int);

  for(i=0;i<I;i++) vw[i] = (w[i] > 0 && v[i] > 0) ? (w[i]-1) * vK + v[i] : -1;
  vwK = int_coder2(vw,I,svw);
  if(vwK > I) error("Levels greater than records! %d > %d\n",vwK,I);

//  printf("Saving % d bytes\n",((int)(2*wK*vK) - (int)(2*(vwK+1) + 2*I)) * (int)sizeof(int));
  /**********************/

    w_cnt = Calloc(wK,           unsigned int);
   vw_cnt = Calloc(vwK+1,        unsigned int);
   sw_cnt = Calloc(wK * sK,      unsigned int);
  svw_cnt = Calloc((vwK+1) * sK, unsigned int);
  
  for(i=0;i<I;i++) {
    if(w[i] > 0) { 
      w_cnt[w[i]-1]++;
      if(s[i] > 0) sw_cnt[(w[i]-1) * sK + s[i]-1]++;
      if(v[i] > 0) {
        vw_cnt[vw[i]-1]++;
        if(s[i] > 0) svw_cnt[(vw[i]-1) * sK + s[i]-1]++;
      }
    }
  }
  Free(vw);

  for(vk=0;vk<vK;vk++) {
    entropy[vk] = 0;
    for(wk=0;wk<wK;wk++) {
      vwk = wk * vK + vk+1;
      vwk = ((int *)bsearch(&vwk,svw,(size_t)vwK,sizeof(int),q_intcmp)) - svw;
      vwk = (0 > vwk || vwk > vwK) ? vwK : vwk; /* for those missing */
      for(sk=0;sk<sK;sk++) {
        if(svw_cnt[vwk*sK + sk] > 0)
          entropy[vk] += svw_cnt[vwk*sK + sk]/(double)I * log(vw_cnt[vwk]/(double)svw_cnt[vwk*sK + sk]);

        if(sw_cnt[wk * sK + sk] - svw_cnt[vwk*sK + sk] > 0)
          entropy[vk] += (sw_cnt[wk * sK + sk] - svw_cnt[vwk*sK + sk])/(double)I * log((w_cnt[wk] - vw_cnt[vwk])/(double)(sw_cnt[wk * sK + sk] - svw_cnt[vwk*sK + sk]));
      }
    }
  }

  Free(svw_cnt);
  Free( sw_cnt);
  Free( vw_cnt);
  Free(  w_cnt);
  Free(svw);
  return(vwK); 
}
   
/* Calculates the marginal entropy for each level of one factor.
 *  v2 -- Uses more memory but much faster
 */
int minfo_v2(int *s, int *v, int *w, unsigned int I, unsigned int sK, unsigned int vK, unsigned int wK, double *entropy) {
  int i,sk,vk,wk,vwk;
  unsigned int *vw_cnt, *svw_cnt, *w_cnt, *sw_cnt;

    w_cnt = Calloc(wK,           unsigned int);
   vw_cnt = Calloc(wK * vK,      unsigned int);
   sw_cnt = Calloc(wK * sK,      unsigned int);
  svw_cnt = Calloc(wK * vK * sK, unsigned int);
  
  for(i=0;i<I;i++) {
    if(w[i] > 0) { 
      w_cnt[w[i]-1]++;
      if(s[i] > 0) sw_cnt[(w[i]-1) * sK + s[i]-1]++;
      if(v[i] > 0) {
        vwk = (w[i]-1) * vK + v[i]-1;
        vw_cnt[vwk]++;
        if(s[i] > 0) svw_cnt[vwk * sK + s[i]-1]++;
      }
    }
  }

  /* Thread This */
  for(vk=0;vk<vK;vk++) {
    entropy[vk] = 0;
    for(wk=0;wk<wK;wk++) {
      vwk = wk * vK + vk;
      for(sk=0;sk<sK;sk++) {
        if(svw_cnt[vwk*sK + sk] > 0)
          entropy[vk] += svw_cnt[vwk*sK + sk]/(double)I * log(vw_cnt[vwk]/(double)svw_cnt[vwk*sK + sk]);

        if(sw_cnt[wk * sK + sk] - svw_cnt[vwk*sK + sk] > 0)
          entropy[vk] += (sw_cnt[wk * sK + sk] - svw_cnt[vwk*sK + sk])/(double)I * log((w_cnt[wk] - vw_cnt[vwk])/(double)(sw_cnt[wk * sK + sk] - svw_cnt[vwk*sK + sk]));
      }
    }
  }
  /***************/

  Free(svw_cnt);
  Free( sw_cnt);
  Free( vw_cnt);
  Free(  w_cnt);
  return(1); 
}
/* Calculates the marginal entropy for each level of one factor.
 *
 */
int minfo_v(int *s, int *v, int *w, unsigned int I, unsigned int sK, unsigned int vK, unsigned int wK, double *entropy) {
  int i,sk,vk,vwk,vwK=2*wK;
  unsigned int *vw_cnt, *svw_cnt;

  for(vk=0;vk<vK;vk++) {
    R_CheckUserInterrupt();
    entropy[vk] = 0;

     vw_cnt = Calloc(vwK,      unsigned int);
    svw_cnt = Calloc(vwK * sK, unsigned int);

    for(i=0;i<I;i++) {
      if(w[i] > 0) {
        vwk = (vk == (v[i]-1)) * wK + w[i]-1;
        vw_cnt[vwk]++;
        if(s[i] > 0) svw_cnt[vwk * sK + s[i]-1]++;
      }
    }

    for(vwk=0;vwk<vwK;vwk++) {
      for(sk=0;sk<sK;sk++) {
        if(svw_cnt[vwk*sK + sk] > 0)
          entropy[vk] += svw_cnt[vwk*sK + sk]/(double)I * log(vw_cnt[vwk]/(double)svw_cnt[vwk*sK + sk]);
      }
    }

    Free(svw_cnt);
    Free(vw_cnt);
  }
  return(1); 
}

SEXP R_minfo(SEXP R_svar, SEXP R_ivars, SEXP R_count, SEXP R_val, SEXP R_progress_flag, SEXP R_entropy_flag) {
  int *w;
  unsigned int vk,wK=1;
  unsigned int i,I=length(R_svar), j,J=nrows(R_ivars),k,sK,vK;
  SEXP R_rv, R_total_entropy;

  w   = Calloc(I, int); for(i=0;i<I;i++) w[i] = 1; /* We do 1 offset */

  sK = nlevels(R_svar);
  PROTECT(R_total_entropy = allocVector(REALSXP,1));
  REAL(R_total_entropy)[0] = total_entropy(INTEGER(R_svar),I,sK);

  if(INTEGER(R_count)[0] == 0) {
    SEXP R_entropy, R_tmp1;
    PROTECT(R_rv = allocVector(VECSXP,J));
    setAttrib(R_rv, R_NamesSymbol, getAttrib(R_ivars, R_NamesSymbol));

    for(j=0;j<J;j++) {
      PROTECT(R_tmp1 = VECTOR_ELT(R_ivars,j));
      vK = nlevels(R_tmp1);
      PROTECT(R_entropy = allocVector(REALSXP,vK));
      setAttrib(R_entropy, R_NamesSymbol, getAttrib(R_tmp1, R_LevelsSymbol));

      minfo_v2(INTEGER(R_svar), INTEGER(R_tmp1), w, I, sK, vK, wK, REAL(R_entropy));
      if(INTEGER(R_entropy_flag)[0] == 0) /* then convert to minfo */
        for(vk=0;vk<vK;vk++) REAL(R_entropy)[vk] = 1 - REAL(R_entropy)[vk]/REAL(R_total_entropy)[0];

      SET_VECTOR_ELT(R_rv,j,R_entropy);
      UNPROTECT(2);
    }
  } else {
    unsigned int k=0,l=0;
    SEXP R_factor_names, R_level_names, R_tmp1, R_tmp2;
    int min_j, min_k, *min_v, L=0, *rv_levels, *v_cnt, *s_cnt, l_count=INTEGER(R_count)[0];
    double *entropy, minfo=0, min_minfo=0, *rv_minfo, l_val =REAL(R_val)[0];
    const char **rv_factor_names, **rv_level_names;
    const char *min_level_name, *min_factor_name;

    PROTECT(R_factor_names = getAttrib(R_ivars, R_NamesSymbol));

    for(j=0;j<J;j++) L += nlevels(VECTOR_ELT(R_ivars,j)); 
    L = L < l_count ? L : l_count;

    rv_factor_names = Calloc(L+1, const char *);
    rv_level_names  = Calloc(L+1, const char *);
    rv_minfo        = Calloc(L+1, double);
    rv_levels       = Calloc(L+1, int);
    v_cnt           = Calloc(L+1, int);
    s_cnt           = Calloc(L+1, int);

    if(INTEGER(R_progress_flag)[0]) {
      printf("%5s %7s %6s %9s %9s %9s %s %s\n",
              "#","minfo","cells","svar","ivar","conv","factor","level");
      fflush(stdout);
    }

    k=0;for(i=0;i<I;i++) s_cnt[l] += INTEGER(R_svar)[i] == 2;
    v_cnt[l] = I;
    rv_factor_names[l] = rv_level_names[l] = "-";

    if(INTEGER(R_progress_flag)[0]) {
      printf("%5d %7.5f %6d %9d %9d %8.4f%% \"%s\" \"%s\"\n",
        l,min_minfo,wK,
        s_cnt[l], v_cnt[l], v_cnt[l] ? 100 * s_cnt[l]/(double)v_cnt[l] : 0.0,
        rv_factor_names[l],rv_level_names[l]);
      fflush(stdout);
    }

    while(l++ < L && minfo < l_val) {
      R_CheckUserInterrupt();
      min_minfo = BIG_NUM;
      for(j=0;j<J;j++) {
        PROTECT(R_tmp1 = VECTOR_ELT(R_ivars,j));
        PROTECT(R_level_names = getAttrib(R_tmp1, R_LevelsSymbol));
        vK = nlevels(R_tmp1);
        entropy = Calloc(vK, double);

        minfo_v2(INTEGER(R_svar), INTEGER(R_tmp1), w, I, sK, vK, wK, entropy);

        for(k=0;k<vK;k++) {
          if(entropy[k] > 0 && min_minfo > entropy[k]) {
            min_j = j;
            min_k = k;
            min_v = INTEGER(R_tmp1);
            min_factor_name = CHAR(STRING_ELT(R_factor_names, min_j));
            min_level_name  = CHAR(STRING_ELT(R_level_names, min_k));
            min_minfo = entropy[min_k];
          }
        }
        Free(entropy);
      }
      if(l > 0 && minfo == 1 - min_minfo / REAL(R_total_entropy)[0]) {
        UNPROTECT(2*J);
        break;
      } else {
        /* Now update w */
        for(i=0;i<I;i++) {
          w[i] = (w[i]-1) * 2 + (min_v[i]-1 == min_k) + 1;

          v_cnt[l] += (min_v[i]-1 == min_k);  /* and gather some statistics */
          s_cnt[l] += (min_v[i]-1 == min_k) * (INTEGER(R_svar)[i] > 1);
//          w_cnt[l] += (w[i] > 1);
//          sw_cnt[l] += (w[i] > 1) * (INTEGER(R_svar)[i] > 1);
        }
//        if(l > 1 && w_cnt[l] == w_cnt[l-1] && sw_cnt[l] == sw_cnt[l-1]) break;

        wK = int_coder(w,I);
        if(wK > I) error("Levels greater than records! %d > %d\n",wK,I);
        UNPROTECT(2*J);

        minfo = 1 - min_minfo / REAL(R_total_entropy)[0];
        if(INTEGER(R_entropy_flag)[0] == 0) min_minfo = minfo;

        rv_factor_names[l] = min_factor_name;
        rv_level_names[l]  = min_level_name;
        rv_minfo[l]        = min_minfo;
        rv_levels[l]       = wK;

        if(INTEGER(R_progress_flag)[0]) {
          printf("%5d %7.5f %6d %9d %9d %8.4f%% \"%s\" \"%s\"\n",
            l,min_minfo,wK,
            s_cnt[l], v_cnt[l], v_cnt[l] ? 100 * s_cnt[l]/(double)v_cnt[l] : 0.0,
            rv_factor_names[l],rv_level_names[l]);
          fflush(stdout);
        }
      }
    }
    UNPROTECT(1);
    k = 7;
    PROTECT(R_rv = allocVector(VECSXP,k));
    PROTECT(R_tmp2 = allocVector(STRSXP,k));

    k = 0;
    SET_STRING_ELT(R_tmp2,k,mkChar("pos"));
    SET_VECTOR_ELT(R_rv,k,allocVector(INTSXP,l));
    PROTECT(R_tmp1 = VECTOR_ELT(R_rv,k));
    for(j=0;j<l;j++) INTEGER(R_tmp1)[j] = j;
    setAttrib(R_rv, install("row.names"), R_tmp1);
    UNPROTECT(1);

    k = 1;
    SET_STRING_ELT(R_tmp2,k,mkChar("factor"));
    SET_VECTOR_ELT(R_rv,k,allocVector(STRSXP,l));
    PROTECT(R_tmp1 = VECTOR_ELT(R_rv,k));
    for(j=0;j<l;j++) SET_STRING_ELT(R_tmp1,j,mkChar(rv_factor_names[j]));
    UNPROTECT(1);

    k = 2;
    SET_STRING_ELT(R_tmp2,k,mkChar("level"));
    SET_VECTOR_ELT(R_rv,k,allocVector(STRSXP,l));
    PROTECT(R_tmp1 = VECTOR_ELT(R_rv,k));
    for(j=0;j<l;j++) SET_STRING_ELT(R_tmp1,j,mkChar(rv_level_names[j]));
    UNPROTECT(1);

    k = 3;
    SET_STRING_ELT(R_tmp2,k,mkChar("minfo"));
    SET_VECTOR_ELT(R_rv,k,allocVector(REALSXP,l));
    memcpy(REAL(VECTOR_ELT(R_rv,k)),rv_minfo,l * sizeof(double));

    k = 4;
    SET_STRING_ELT(R_tmp2,k,mkChar("cells"));
    SET_VECTOR_ELT(R_rv,k,allocVector(INTSXP,l));
    memcpy(INTEGER(VECTOR_ELT(R_rv,k)),rv_levels,l * sizeof(int));

    k = 5;
    SET_STRING_ELT(R_tmp2,k,mkChar("svar"));
    SET_VECTOR_ELT(R_rv,k,allocVector(INTSXP,l));
    memcpy(INTEGER(VECTOR_ELT(R_rv,k)),s_cnt,l * sizeof(int));

    k = 6;
    SET_STRING_ELT(R_tmp2,k,mkChar("ivar"));
    SET_VECTOR_ELT(R_rv,k,allocVector(INTSXP,l));
    memcpy(INTEGER(VECTOR_ELT(R_rv,k)),v_cnt,l * sizeof(int));

    setAttrib(R_rv, R_NamesSymbol, R_tmp2);
    UNPROTECT(1);

    PROTECT(R_tmp1 = allocVector(STRSXP, 1));
    SET_STRING_ELT(R_tmp1,0,mkChar("data.frame"));
    classgets(R_rv,R_tmp1);
    UNPROTECT(1);

    Free(rv_factor_names);
    Free(rv_level_names);
    Free(rv_minfo);
    Free(rv_levels);
    Free(v_cnt);
    Free(s_cnt);
  }
  setAttrib(R_rv,install("total_entropy"),R_total_entropy);

  UNPROTECT(2);
  return(R_rv);
}

/**************************************************************************************
    minfo -- calculates the mutual information coefficient (between 0 and 1) based
        upon shannon information for each variable.
**************************************************************************************/
int minfo_old(unsigned int *data, unsigned int *state, unsigned int I, unsigned int J, double *mi) {
  unsigned int i,j,k,K,l,L;
  unsigned int *state_counts, **var_counts, **tot_counts, *Lv;
  double tot_entropy=0;

  /* Calc maxes and make storage */
  K = state[0];
  for(i=1;i<I;i++) if(K < state[i]) K = state[i];
  K++;

  state_counts = Calloc(K,unsigned int);
  Lv           = Calloc(J,unsigned int);
  var_counts   = Calloc(J,unsigned int *);  /* list of vectors */
  tot_counts   = Calloc(J,unsigned int *);  /* list of matricies */

  for(j=0;j<J;j++) {
    Lv[j] = data[j*I + 0];
    for(i=1;i<I;i++) if(Lv[j] < data[j*I+i]) Lv[j] = data[j*I+i];
    Lv[j]++;
    var_counts[j] = Calloc(Lv[j],unsigned int);
    tot_counts[j] = Calloc(Lv[j]*K,unsigned int);
  }
  /* Tabulate the counts */
  for(i=0;i<I;i++) {
    R_CheckUserInterrupt();
    k = state[i];
    state_counts[k]++;
    for(j=0;j<J;j++) {
      L = Lv[j];
      l = data[j*I+i];
      var_counts[j][l]++;
      tot_counts[j][k*L+l]++;
    }
  }
  /* Calculate the entropy of each Variable */
  for(j=0;j<J;j++) {
    R_CheckUserInterrupt();
    mi[j] = 0;
    L = Lv[j]; 
    for(l=0;l<L;l++) { 
//      if(var_counts[j][l]) {
        for(k=0;k<K;k++) { if(tot_counts[j][k*L+l]) {
//            mi[j] += var_counts[j][l]/(double)I * 
//                     tot_counts[j][k*L+l]/(double)var_counts[j][l] * 
//                     log(var_counts[j][l]/(double)tot_counts[j][k*L+l]);
            mi[j] += tot_counts[j][k*L+l]/(double)I * 
                     log(var_counts[j][l]/(double)tot_counts[j][k*L+l]);
          }
        }
//      }
    }
  }
  /* Calculate the total entropy */
  for(k=0;k<K;k++) {
    if(state_counts[k]) {
     tot_entropy += state_counts[k]/(double)I * log(I/(double)state_counts[k]); 
    }
  }
  /* Calc the information */
  for(j=0;j<J;j++) {
    mi[j] = 1 - mi[j]/tot_entropy;
    if(mi[j] < 0) mi[j] = 0;
    if(mi[j] > 1) mi[j] = 1;
  }
  /* Free the memory */
  Free(state_counts);
  Free(Lv);
  for(j=0;j<J;j++) {
    Free(var_counts[j]);
    Free(tot_counts[j]);
  }
  Free(var_counts);
  Free(tot_counts);
  return(0);
}

int q_strcmp(const void *a, const void *b) {
    return(strcmp(*((char **)a),*((char **)b)));
}

unsigned int str_coder(char **v, unsigned int N, unsigned int *rv) {
    unsigned int I=1,n;
    char **sv;

    sv = Calloc(N,char *);
    memcpy(sv,v,N*sizeof(char *));

    qsort(v,(size_t)N,sizeof(char*),q_strcmp);

    for(n=1;n<N;n++)
        if(strcmp(v[n],v[n-1])) v[I++] = v[n];

    for(n=0;n<N;n++)
        rv[n] = ((char**)bsearch(sv + n,v,(size_t)I,sizeof(char*),q_strcmp)) - v;

    Free(sv);
    return(I);
}

int q_ushortcmp(const void *a, const void *b) {
    return(*((unsigned short *)a) - *((unsigned short *)b));
}

unsigned int uint_coder(unsigned int *v, unsigned int N, unsigned int *rv) {
    unsigned int I=1,n;
    unsigned int *sv;

    sv = Calloc(N,unsigned int);
    memcpy(sv,v,N*sizeof(unsigned int));

    qsort(sv,(size_t)N,sizeof(unsigned int),q_intcmp);

    for(n=1;n<N;n++)
        if(sv[n] != sv[n-1]) sv[I++] = sv[n];

    for(n=0;n<N;n++)
        rv[n] = ((unsigned int *)bsearch(v + n,sv,(size_t)I,sizeof(unsigned int),q_intcmp)) - sv;

    Free(sv);
    return(I);
}
unsigned short ushort_coder(unsigned short *v, unsigned int N, unsigned short *rv) {
    unsigned int I=1,n;
    unsigned short *sv;

    sv = Calloc(N,unsigned short);
    memcpy(sv,v,N*sizeof(unsigned short));

    qsort(v,(size_t)N,sizeof(unsigned short),q_ushortcmp);

    for(n=1;n<N;n++)
        if(v[n] != v[n-1]) v[I++] = v[n];

    for(n=0;n<N;n++)
        rv[n] = ((unsigned short *)bsearch(sv + n,v,(size_t)I,sizeof(unsigned short),q_ushortcmp)) - v;

    Free(sv);
    return(I);
}


SEXP R_minfo_old(SEXP R_data, SEXP R_state, SEXP R_iterate, SEXP R_progress, SEXP R_count, SEXP R_val) {
    unsigned int i,I = nrows(R_data), j,J = ncols(R_data);
    unsigned int *data, *state;
    char **s;

    SEXP R_v, R_names, R_mi;

    if(length(R_state) != I) 
        error("State vector is not the same length as rows in data! nrow(data) %d len(state) %d\n",I,length(R_state));

    s     = Calloc(I,char *);
    state = Calloc(I,unsigned int);

    PROTECT(R_names = VECTOR_ELT(getAttrib(R_data,R_DimNamesSymbol),1));
    if(isNull(R_names)) {
        UNPROTECT(1);
        PROTECT(R_names = coerceVector(mkSequence(1.0,(double)J,1.0),STRSXP));
    }

    PROTECT(R_v = coerceVector(R_state,STRSXP));
    for(i=0;i<I;i++) s[i] = (char *)CHAR(STRING_ELT(R_v,i));
    j = str_coder(s,I,state);
    UNPROTECT(1);

    PROTECT(R_v = coerceVector(R_data,STRSXP));
    data = Calloc(I*J,unsigned int);
    for(j=0;j<J;j++) {
        for(i=0;i<I;i++) s[i] = (char *)CHAR(STRING_ELT(R_v,j*I+i));
        i = str_coder(s,I,data+j*I);
    }
    UNPROTECT(1);

    if(INTEGER(R_iterate)[0]) {
        double *mi,max_val;
        unsigned int *itmp,stop=0,K,r;
        int *idx,max_idx;

        SEXP R_nnames;

        mi = Calloc(J,double);
        itmp = Calloc(I,unsigned int);
        idx  = Calloc(J,int);
        for(j=0;j<J;j++) idx[j] = j;

        r=0;
        while(!stop) {
            R_CheckUserInterrupt();
            minfo_old(data+r*I,state,I,J-r,mi+r);
            // now find the one with the highest minfo
            max_val = 0.0F;max_idx = 0;
            for(j=r;j<J;j++) if(max_val < mi[j]) {max_val = mi[j]; max_idx = j;}
            // swap that column to column "r"
            mi[r] = max_val;
            memcpy(itmp, data+max_idx*I, I*sizeof(unsigned int));
            memcpy(data+max_idx*I, data+r*I, I*sizeof(unsigned int));
            memcpy(data+r*I, itmp, I*sizeof(unsigned int));
            j = idx[max_idx]; idx[max_idx] = idx[r]; idx[r] = j;

            s[r] = (char *)CHAR(STRING_ELT(R_names,idx[r]));

            // find max data value
            K = data[r*I+0]; for(i=1;i<I;i++) if(K < data[r*I+i]) K = data[r*I+i]; K++;
            if(INTEGER(R_progress)[0]) Rprintf("%5u %.5f %s %d levels\n",r+1,mi[r],s[r],K);
//            Rprintf("     ");for(i=0;i<I;i++) Rprintf("%2d ",data[r*I+i]); Rprintf("\n");
            if(mi[r] < ((REAL(R_val)[0] < 1.0F) ? REAL(R_val)[0] : 1.0F) 
                && r+1 < ((INTEGER(R_count)[0] < J) ? INTEGER(R_count)[0] : J)) {
                for(j=r+1;j<J;j++)  {
//                    Rprintf("j %d idx[j] %d %s:\n",j,idx[j],CHAR(STRING_ELT(R_names,idx[j])));
//                    Rprintf("     ");for(i=0;i<I;i++) Rprintf("%2d ",data[j*I+i]); Rprintf("\n");
                    for(i=0;i<I;i++) itmp[i] = data[j*I+i] * K + data[r*I+i];
//                    Rprintf("     ");for(i=0;i<I;i++) Rprintf("%2d ",itmp[i]); Rprintf("\n");
                    uint_coder(itmp,I,data+j*I);
//                    Rprintf("     ");for(i=0;i<I;i++) Rprintf("%2d ",data[j*I+i]); Rprintf("\n");
                }
            } else stop = 1;
            r++;
        }
        PROTECT(R_mi = allocVector(REALSXP,r));
        memcpy(REAL(R_mi),mi,r*sizeof(double));

        PROTECT(R_nnames = allocVector(STRSXP,r));
        for(i=0;i<r;i++) SET_STRING_ELT(R_nnames,i,mkChar(s[i]));
        setAttrib(R_mi,R_NamesSymbol,R_nnames);
        UNPROTECT(1);

        Free(mi);
        Free(itmp);
        Free(idx);
    } else {
        PROTECT(R_mi = allocVector(REALSXP,J));
        minfo_old(data,state,I,J,REAL(R_mi));
        setAttrib(R_mi,R_NamesSymbol,R_names);
    }

    Free(s);
    Free(state);
    UNPROTECT(2);
    return(R_mi);
}
/**************************************************************************************
    Calculates the maximum INTEGER in each column and returns a vector
**************************************************************************************/
SEXP max_cols(SEXP x) {
    int i,j,tmp;
    int rows;
    int cols;
    SEXP max;

    if(isMatrix(x)) {
        rows = nrows(x);
        cols = ncols(x);
    } else {
        rows = length(x);
        cols = 1;
    }

    PROTECT(max = allocVector(INTSXP, cols));

    /* This takes care of all the NANs and NAs etc */
    for(j=0;j<cols;j++) {
        INTEGER(max)[j] = (INTEGER(x)[j*rows] > 0) ? INTEGER(x)[j*rows] : 0;
    }
    for(i=1;i<rows;i++) {
        for(j=0;j<cols;j++) {
            if((tmp = INTEGER(x)[j*rows + i]) > INTEGER(max)[j]) {
               INTEGER(max)[j] = tmp;
            }
        }
    }
    UNPROTECT(1);
    return(max);
}
            
    
SEXP minfo_old_old(SEXP data, SEXP state) {
    int i,I = nrows(data);
    int j,J = ncols(data);
    int k,K = INTEGER(max_cols(state))[0] + 1;
    int l,L;
    double tot_entropy = 0.0;
    SEXP minfo,state_counts,rec_counts,var_counts,tot_counts,Lv,tmp1,tmp2;

    if(!I || !J || !K || I != nrows(state)) error("Invalid data set.  I=%d J=%d K=%d\n",I,J,K);

    /* Allocate and Zero all the space */
    PROTECT(minfo        = allocVector(REALSXP, J)); for(j=0;j<J;j++) REAL(minfo)[j]           = 0.0;
    PROTECT(state_counts = allocVector(INTSXP, K));  for(k=0;k<K;k++) INTEGER(state_counts)[k] = 0;
    PROTECT(rec_counts   = allocVector(INTSXP, J));
    PROTECT(var_counts   = allocVector(VECSXP, J));
    PROTECT(tot_counts   = allocVector(VECSXP, J));
    PROTECT(Lv = max_cols(data)); 
    for(j=0;j<J;j++) {
        INTEGER(rec_counts)[j] = 0;

        L = INTEGER(Lv)[j] + 1;
        PROTECT(tmp1 = allocVector(INTSXP, L));
        PROTECT(tmp2 = allocMatrix(INTSXP, L, K));
        for(l=0;l<L;l++) {
            INTEGER(tmp1)[l] = 0;
            for(k=0;k<K;k++) {
                INTEGER(tmp2)[k*L + l] = 0;
            }
        }
        SET_VECTOR_ELT(var_counts,j,tmp1);
        SET_VECTOR_ELT(tot_counts,j,tmp2);
        UNPROTECT(2);
    }

    /* Tabulate the Counts */
    for(i=0;i<I;i++) {
        k = INTEGER(state)[i];
        if(k >= 0) {                                  /* This takes care of NA values */
            INTEGER(state_counts)[k]++;
            for(j=0;j<J;j++) {
                l = INTEGER(data)[j*I + i];
                if(l >= 0) {
                    tmp1 = VECTOR_ELT(var_counts,j);
                    tmp2 = VECTOR_ELT(tot_counts,j);
                    L = nrows(tmp2);

                    INTEGER(rec_counts)[j]++;
                    INTEGER(tmp1)[l]++;
                    INTEGER(tmp2)[k*L + l]++;
                }
            }
        }
    }

    /* Calculate the Entropy of each Variable */
    for(j=0;j<J;j++) {                                     /* i */
        i = INTEGER(rec_counts)[j];                        /* This is just I if there are no NA values */
        if(i) {
            tmp1 = VECTOR_ELT(var_counts,j);
            tmp2 = VECTOR_ELT(tot_counts,j);
            L = nrows(tmp2);
            for(l=0;l<L;l++) {                             /* j */
                if(INTEGER(tmp1)[l]) {
                    for(k=0;k<K;k++) {                     /* k */
                        if(INTEGER(tmp2)[k*L + l]) {
                            REAL(minfo)[j] += INTEGER(tmp1)[l]/(double)i *
                                INTEGER(tmp2)[k*L + l]/(double)INTEGER(tmp1)[l] *
                                log(INTEGER(tmp1)[l]/(double)INTEGER(tmp2)[k*L + l]);
                        }
                    }
                }
            }
        }
    }
    /* Calculate the total Entropy */
    for(k=0;k<K;k++) {
        if(INTEGER(state_counts)[k])
            tot_entropy += INTEGER(state_counts)[k]/(double)I * log(I/(double)INTEGER(state_counts)[k]);
    }
    if(tot_entropy < REAL_ZERO) error("Error calculating total entropy! e=%E\n",tot_entropy);

    /* Now calc the information */        
    for(j=0;j<J;j++) {
        REAL(minfo)[j] = 1 - REAL(minfo)[j]/tot_entropy;
        if(REAL(minfo)[j] < 0) REAL(minfo)[j] = 0;   /* This can happen with NA values! */
        if(REAL(minfo)[j] > 1) REAL(minfo)[j] = 1;
    }
        
    UNPROTECT(6);
    return(minfo);
}

SEXP sum_cols(SEXP R_X) {
    int m,M = nrows(R_X);
    int n,N = ncols(R_X);
    
    double *X = REAL(R_X);
    double *rv; 
    
    SEXP R_rv;
    
    PROTECT(R_rv = allocVector(REALSXP,M));
    rv = REAL(R_rv);

    
    for(m=0;m<M;m++) {
        rv[m] = 0.0;
        for(n=0;n<N;n++) {
            rv[m] += X[n*M+m];
        }   
    }       
    UNPROTECT(1);
    return(R_rv);
}   


/*******************************************************************
# lams.init initializes the lambdas using conditional probabilties given the state variable.
#    It then uses "auditor pesimisum" to re-adjust the probabilities in order to remove
#    any zero probabilities.  i.e.  if something happens 10/10 times then the probabilities
#    become 10/11 of happening and 1/11 of not happening.
#
#    data  -- matrix of factor levels.  Must be all integers.
#    state -- vector of factor levels.  Must be all integers.
#*******************************************************************/
SEXP gomlams_init(SEXP data, SEXP state) {
    int i,I=nrows(data);
    int j,J=ncols(data);
    int k,K=INTEGER(max_cols(state))[0] + 1;
    long l,L;
    double sum;
    SEXP lams,Lv,counts,lamM;

    /* Allocate all the space needed */
    PROTECT(counts = allocVector(INTSXP,K));
    for(k=0;k<K;k++) {
        INTEGER(counts)[k] = 0;
    }

    PROTECT(lams   = allocVector(VECSXP,J));
    PROTECT(Lv = max_cols(data)); 
    for(j=0;j<J;j++) {
        L = INTEGER(Lv)[j] + 1;
        PROTECT(lamM = allocMatrix(REALSXP,L,K));
        for(l=0;l<L;l++) {
            for(k=0;k<K;k++) {
                REAL(lamM)[k*L + l] = 0.0;
            }
        }
        SET_VECTOR_ELT(lams,j,lamM);
        UNPROTECT(1);
    }

    /* Tabulate the counts */
    for(i=0;i<I;i++) {
        k = INTEGER(state)[i];
        if(k >= 0) {
            INTEGER(counts)[k]++;
            for(j=0;j<J;j++) {
                l = INTEGER(data)[j*I + i];
                if(l >= 0) {
                    lamM = VECTOR_ELT(lams,j);
                    L   = nrows(lamM);
                    REAL(lamM)[k*L + l]++;
                }
            }
        }
    }
    /* Now make the conditional Probabilities */
    for(j=0;j<J;j++) {
        lamM = VECTOR_ELT(lams,j);
        L   = nrows(lamM);
        for(k=0;k<K;k++) {
            sum = 0;
            for(l=0;l<L;l++) {
                if(INTEGER(counts)[k]) {
                    if(REAL(lamM)[k*L + l] > REAL_ZERO) {
                       REAL(lamM)[k*L + l] /= INTEGER(counts)[k];
                    } else {
                       REAL(lamM)[k*L + l] = 1.0/INTEGER(counts)[k];
                    }
                    sum += REAL(lamM)[k*L + l];
                } else {
                    if(REAL(lamM)[k*L + l] > REAL_ZERO) {
                        error("State factor level %d never observed and something is wrong.\n",k);
                    } else {
                       REAL(lamM)[k*L + l] = 1.0/L;
                    }
                    sum += 1.0/L;
                }
            }
            if(sum < REAL_ZERO) error("InitLams:sum should NOT be zero: sum=%E threshold=%E j=%d k=%d K=%d\n",sum,REAL_ZERO,j,k,K);

            /* This does the "pessimistic auditor" adjustment */
            for(l=0;l<L;l++) {
                REAL(lamM)[k*L + l] /= sum;
            }
        }
    }
    UNPROTECT(3);
    return(lams);
}

void dmemset(double *m, double s, int N) {
    int n;
    
    for(n=0;n<N;n++) {
        m[n] = s;
    }
}   
void imemset(int *m, int s, int N) {
    int n;  
                
    for(n=0;n<N;n++) {
        m[n] = s;
    }   
}   

int isum(int *x, int N) {
    int n,sum=0;

    if(N<0) return(0);

    for(n=0;n<N;n++) {
        sum += x[n];
    }
    return(sum);
}

/* Returns a K*Lt mattrix of doubles */
SEXP ipflams_init(SEXP R_data, SEXP R_state, SEXP R_Lv) {
    int i,I=nrows(R_data);
    int j,J=ncols(R_data);
    int k,K=INTEGER(max_cols(R_state))[0] + 1;
    long l,L;

    int *data  = INTEGER(R_data);
    int *state = INTEGER(R_state);
    int *Lv    = INTEGER(R_Lv);
    int *counts,jsum;
    double *lams,*sums;

    SEXP R_lams;

    L = isum(Lv,J);
    /* Allocate all the space needed */
    PROTECT(R_lams = allocMatrix(REALSXP,L,K));
    lams = REAL(R_lams);
    dmemset(lams,0,L*K);

    counts = Calloc(K,int);
    imemset(counts,0,K); /* Zero them all */

    sums   = Calloc(J,double);
    dmemset(sums,0,J);

    /* Tabulate the counts */
    for(i=0;i<I;i++) {
        k = state[i];
        if(k >= 0) {
            counts[k]++;
            for(j=0;j<J;j++) {
                l = data[j*I+i];
                if(l >= 0) { /* avoid the NAs etc */
                    lams[k*L+l]++;
                }
            }
        }
    }
    /* Now make the conditional Probabilities */
    for(k=0;k<K;k++) {
        j=0;
        jsum=Lv[0];
        for(l=0;l<L;l++) {
            if(l >= jsum) {j++;jsum+=Lv[j];}

            if(counts[k]) {
                if(lams[k*L+l] > 0) {
                   lams[k*L+l] /= counts[k];
                } else {
                   lams[k*L+l] = 1/counts[k];
                }
            } else {
                if(lams[k*L+l] > 0) {
                    error("State factor level %d never observed and something is wrong.\n",k);
                } else {
                   lams[k*L+l] = 1/Lv[j];
                }
            }
            sums[j] += lams[k*L+l];
        }
    }
        /* This does the "pessimistic auditor" adjustment */
    for(k=0;k<K;k++) {
        j=0;
        jsum=Lv[0];
        for(l=0;l<L;l++) {
            if(l >= jsum) {j++;jsum+=Lv[j];}
            lams[k*L+l] = (Lv[j]>1) ? (1-lams[k*L+l]/sums[j])/(Lv[j] - 1) : 1;
        }
    }
    Free(counts);
    Free(sums);

    UNPROTECT(1);
    return(R_lams);
}
/*******************************************************************
# giks_init initializes the giks by summing each observation's lams
#     accross J,K then normalizing over K.
*******************************************************************/
SEXP giks_init(SEXP R_data, SEXP R_lams) {
    int i,I = nrows(R_data);
    int j,J = ncols(R_data);
    int k,K = ncols(R_lams);
    int l,L = nrows(R_lams);
    double sum;

    int    *data;
    double *lams;
    double *giks;
    SEXP R_giks, R_X, R_Y;

    PROTECT(R_giks = allocMatrix(REALSXP,I,K));
    giks = REAL(R_giks);
    PROTECT(R_X = coerceVector(R_data,INTSXP));
    data = INTEGER(R_X);
    PROTECT(R_Y = coerceVector(R_lams,REALSXP));
    lams = REAL(R_Y);

    for(i=0;i<I;i++) {
        for(k=0;k<K;k++) {     /* Intialize to zero... */ 
            giks[k*I+i] = 0.0; /* You'd think the allocMatrix would do this! */
        }

        sum = 0;
        for(j=0;j<J;j++) {
            l = data[j*I+i];
            if(l >= 0) {
                for(k=0;k<K;k++) {
                    giks[k*I+i] += lams[k*L+l];
                    sum         += lams[k*L+l];
                }
            }
        }
        if(sum < 0) 
            error("InitGiks:sum should NOT be zero: sum=%E i=%d\n",sum,i);
        for(k=0;k<K;k++) {
            giks[k*I+i] /= sum;
        }
    }
    UNPROTECT(3);
    return(R_giks);
}
SEXP ipfgiks_init(SEXP R_data, SEXP R_lams) {
    int i,I = nrows(R_data);
    int j,J = ncols(R_data);
    int k,K = ncols(R_lams);
    int l,L = nrows(R_lams);
    double sum;

    int    *data = INTEGER(R_data);
    double *lams = REAL(R_lams);
    double *giks;
    SEXP R_giks;

    PROTECT(R_giks = allocMatrix(REALSXP,I,K));
    giks = REAL(R_giks);

    for(i=0;i<I;i++) {
        for(k=0;k<K;k++) {     /* Intialize to zero... */ 
            giks[k*I+i] = 0.0; /* You'd think the allocMatrix would do this! */
        }

        sum = 0;
        for(j=0;j<J;j++) {
            l = data[j*I+i];
            if(l >= 0) {
                for(k=0;k<K;k++) {
                    giks[k*I+i] += exp(-lams[k*L+l]);
                    sum         += exp(-lams[k*L+l]);
                }
            }
        }
        if(sum < 0) 
            error("InitGiks:sum should NOT be zero: sum=%E i=%d\n",sum,i);
        for(k=0;k<K;k++) {
            giks[k*I+i] /= sum;
        }
    }
    UNPROTECT(1);
    return(R_giks);
}


SEXP combine(SEXP by, SEXP offset, SEXP list, SEXP delim) {
    int i,j,k=0;
    char *s1,*s2,*s3,*sd;
    SEXP tmplist,tmpby,tmpoffset;
    SEXP rv;

    PROTECT(rv = allocVector(STRSXP,choose(length(list)-INTEGER(offset)[0],INTEGER(by)[0])));

    if(INTEGER(by)[0] > 1) {
        sd = (char *)CHAR(STRING_ELT(delim,0));
        for(i=INTEGER(offset)[0];i<length(list)-INTEGER(by)[0]+1;i++) {
            s1 = (char *)CHAR(STRING_ELT(list,i));
            s2 = Calloc(strlen(s1)+strlen(sd)+1,char);
            strcpy(s2,s1);
            strcat(s2,sd);

            PROTECT(tmpby     = allocVector(INTSXP,1));
            PROTECT(tmpoffset = allocVector(INTSXP,1));

            INTEGER(tmpby)[0]     = INTEGER(by)[0]-1;
            INTEGER(tmpoffset)[0] = i+1;

            PROTECT(tmplist = combine(tmpby,tmpoffset,list,delim));

            for(j=0;j<length(tmplist);j++) {
                s1 = (char *)CHAR(STRING_ELT(tmplist,j));
                s3 = R_alloc(strlen(s2)+strlen(s1)+1,sizeof(char));
                strcpy(s3,s2);
                strcat(s3,s1);
                SET_STRING_ELT(rv,k++,mkChar(s3));
            }
            UNPROTECT(3);
            Free(s2);
        }
    } else {
        for(i=INTEGER(offset)[0];i<length(list);i++) {
            s1 = (char *)CHAR(STRING_ELT(list,i));
            s2 = R_alloc(strlen(s1)+1,sizeof(char));
            strcpy(s2,s1);
            SET_STRING_ELT(rv,k++,mkChar(s2));
        }
    }

    UNPROTECT(1);
    return(rv);
}

/* Combine -- non-recursive formula.
   Works on string matrix X
   Returns string matrix Y
*/
#define COMBINE_BUFF_SIZE 1024
SEXP R_combine(SEXP R_X, SEXP R_by, SEXP R_sep) {
    int i,I = nrows(R_X);
    int j,J = ncols(R_X);
    int m,*idx;
    int b,k,l,by = INTEGER(R_by)[0];
    char stop=0,*sbuff;

    SEXP R_Y;

    if(J <= by || by < 2) return(R_X);

//    Rprintf("I %d J %d choose %d\n",I,J,(int)choose((double)J,(double)by));

    PROTECT(R_Y = allocMatrix(STRSXP,I,(int)choose((double)J,(double)by)));
    PROTECT(R_X = coerceVector(R_X,STRSXP));

    sbuff = Calloc(COMBINE_BUFF_SIZE,char);
    idx   = Calloc(by,int);
    for(b=0;b<by;b++) idx[b] = b;
    
    m=0;
    while(!stop) {
//        Rprintf("m %d idx %d %d %d\n",m,idx[0],idx[1],idx[2]);
        for(i=0;i<I;i++) {
            R_CheckUserInterrupt();
            k=0; 
            while(k < COMBINE_BUFF_SIZE && (sbuff[k] = CHAR(STRING_ELT(R_X,idx[0]*I+i))[k])) 
                k++;
            for(b=1;b<by;b++) {
                l=0;
                while(k < COMBINE_BUFF_SIZE && (sbuff[k] = CHAR(STRING_ELT(R_sep,0))[l++])) 
                    k++;
                l=0;
                while(k < COMBINE_BUFF_SIZE && (sbuff[k] = CHAR(STRING_ELT(R_X,idx[b]*I+i))[l++])) 
                    k++;
            }  // Copies '\0' from last string into sbuff
            SET_STRING_ELT(R_Y, m*I+i, mkChar(sbuff));
        }
        m++;
        // This is where the magic takes place!
        idx[by-1]++;
        for(b=by-1;b>0;b--) if(idx[b] > J-(by-b)) idx[b-1]++;
        for(b=1;b<by;b++)   if(idx[b] > J-(by-b)) idx[b] = idx[b-1]+1;
        if(idx[0] > J-by) stop = 1;
    }

    Free(sbuff);
    Free(idx);
    UNPROTECT(2);
    return(R_Y);
}
    

void check_val(double *x) {
    if(*x < -BIG_NUM) *x = -BIG_NUM;
    else
    if(*x >  BIG_NUM) *x =  BIG_NUM;
    else
    if(isnan(*x))     *x =  0.0F;
}

/****************************************************************************
    pthread implementation!
****************************************************************************/
typedef struct {
    int *data, *counts;
    double *lams, *giks;
    int I,J,K,L;
} ipf_L_struct;
typedef struct {
    ipf_L_struct *iL;
    double *lik;
    int i1,i2;
} pipf_L_struct;
    
void *pipf_L(void *arg) {
    int i,j,k,l;
    double tmp1;
    pipf_L_struct *Ld = arg;

    for(i=Ld->i1;i<Ld->i2;i++) {
        for(j=0;j<Ld->iL->J;j++) {
            l = Ld->iL->data[j*Ld->iL->I+i];
            if(l != NA_INTEGER) {
                tmp1=0;
                for(k=0;k<Ld->iL->K;k++) {
                    tmp1 += Ld->iL->giks[k*Ld->iL->I+i] * 
                        Ld->iL->lams[k*Ld->iL->L+l] * 
                        Ld->iL->counts[i];
                }
                *Ld->lik += tmp1;
            }
        }
        sched_yield();
    }
    pthread_exit((void *) 0);
}
SEXP R_pipf_L(SEXP R_data, SEXP R_counts, SEXP R_lams, SEXP R_giks, SEXP R_threads) {
    int I = nrows(R_data);
    int J = ncols(R_data);
    int K = ncols(R_giks);
    int L = nrows(R_lams);
    int t,T = (*INTEGER(R_threads) < I) ? *INTEGER(R_threads) : I;
    int status,len,idx=0;

    double *lik;
    ipf_L_struct  iL;
    pipf_L_struct *Ld;
    pthread_attr_t attr;
    pthread_t *threads;

    SEXP R_Lik;

    lik = Calloc(T,double); // And sets to zero!
    Ld  = Calloc(T,pipf_L_struct);
    threads = Calloc(T,pthread_t);

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    pthread_attr_setschedpolicy(&attr, SCHED_RR);
//    pthread_attr_setscope(&attr, PTHREAD_SCOPE_PROCESS);
//    pthread_setconcurrency(T);
    
    iL.data   = INTEGER(R_data);
    iL.counts = INTEGER(R_counts);
    iL.lams   = REAL(R_lams);
    iL.giks   = REAL(R_giks);
    iL.I      = I;
    iL.J      = J;
    iL.K      = K;
    iL.L      = L;

    len = (int) ceil((double)I/(double)T);
//    Rprintf("R_pipf_L len %d\n",len);
    for(t=0;t<T;t++) {
        Ld[t].iL  = &iL;
        Ld[t].lik = &lik[t];
        Ld[t].i1  = idx; idx += len;
        Ld[t].i2  = (idx < I) ? idx : I;
        pthread_create(&threads[t], &attr, pipf_L, (void *) &Ld[t]);
    }
    for(t=0;t<T;t++)
        pthread_join(threads[t], (void **) &status);
    PROTECT(R_Lik = allocVector(REALSXP,1));
    *REAL(R_Lik) = 0.0F;
    for(t=0;t<T;t++)
        *REAL(R_Lik) += lik[t];
        
    Free(lik);
    Free(Ld);
    Free(threads);
    UNPROTECT(1);
    return(R_Lik);
}
    

/*****************************************************************************************
 * Log likelihood of the exponetial gom
*****************************************************************************************/
double ipf_L_raw(
    int *data,
    int *counts,
    double *lams,
    double *giks,
    int I,
    int J,
    int K,
    int L) {

    int i,j,k,l;
    double tmp1,Lik=0;

    for(i=0;i<I;i++) {
        for(j=0;j<J;j++) {
            l = data[j*I+i];
            if(l >= 0) {
                tmp1=0;
                for(k=0;k<K;k++) {
                    tmp1 += giks[k*I+i] * lams[k*L+l] * counts[i];
                }
                Lik += tmp1;
            }
        }
    }
    check_val(&Lik);
    return(Lik);
}


SEXP ipf_L(SEXP R_data, SEXP R_counts, SEXP R_lams, SEXP R_giks) {
    int I = nrows(R_data);
    int J = ncols(R_data);
    int K = ncols(R_giks);
    int L = nrows(R_lams);

    SEXP R_Lik;

    PROTECT(R_Lik = allocVector(REALSXP,1));
    REAL(R_Lik)[0] = ipf_L_raw(INTEGER(R_data),INTEGER(R_counts),REAL(R_lams),REAL(R_giks),I,J,K,L);

    UNPROTECT(1);
    return(R_Lik);
}

typedef struct {
    int *data, *counts;
    double *lams, *giks, *alphas, *betas;
    double w;
    int I,J,K,L,*Lv;
} ipf_Lc_struct;
typedef struct {
    ipf_Lc_struct *iLc;
    double *Lc,*rms;
    int i1,i2;
} pipf_Lc_struct;

void *pipf_Lc(void *arg) {
    int i,j,k,l,lsum;
    int I,J,K,L,*Lv;
    double tmp1,tmp2;
    pipf_Lc_struct *Lcd = arg;

    I = Lcd->iLc->I;
    J = Lcd->iLc->J;
    K = Lcd->iLc->K;
    L = Lcd->iLc->L;
    Lv = Lcd->iLc->Lv;

    for(i=Lcd->i1;i<Lcd->i2;i++) {
        lsum=0;
        for(j=0;j<J;j++) {
            l = Lcd->iLc->data[j*I+i];
            if(l != NA_INTEGER) {
                tmp1=0;
                for(k=0;k<K;k++) {
                    tmp1 += Lcd->iLc->giks[k*I+i] * Lcd->iLc->lams[k*L+l];
                }
                *Lcd->Lc += Lcd->iLc->counts[i] * tmp1;
            }
            tmp2=0;
            for(l=lsum;l<lsum+Lv[j];l++) {
                tmp1=0;
                for(k=0;k<K;k++) {
                    tmp1 += Lcd->iLc->giks[k*I+i] * Lcd->iLc->lams[k*L+l];
                }
                tmp2 += exp(-tmp1);
            }
            lsum += Lv[j];
            tmp2 -= 1;
            *Lcd->Lc += Lcd->iLc->betas[j*I+i] * tmp2;
            tmp2 *= tmp2;
            *Lcd->Lc += Lcd->iLc->w/2 * tmp2;
            *Lcd->rms += tmp2;
        }
        tmp1=0;
        for(k=0;k<K;k++) {
            tmp1 += Lcd->iLc->giks[k*I+i];
        }
        tmp1 -= 1;
        *Lcd->Lc += Lcd->iLc->alphas[i] * tmp1;
        tmp1 *= tmp1;
        *Lcd->Lc += Lcd->iLc->w/2 * tmp1;
        *Lcd->rms += tmp1;
        sched_yield();
    }
    pthread_exit((void *) 0);
}
/* Done this way to accomidate the optim command */
SEXP R_pipf_Lc(
            SEXP R_x, /* Contains Lams, Giks and Alphas, Betas */
            SEXP R_data, 
            SEXP R_counts, 
            SEXP R_w,
            SEXP R_K,
            SEXP R_Lv, SEXP R_threads) {
    int I     = nrows(R_data);
    int J     = ncols(R_data);
    int K     = INTEGER(R_K)[0];
    int L     = isum(INTEGER(R_Lv),J);
    int status,idx,len,t,T = (*INTEGER(R_threads) < I) ? *INTEGER(R_threads) : I;
    double *Lc,*rms;

    ipf_Lc_struct iLc;
    pipf_Lc_struct *Lcd;
    pthread_attr_t  attr;
    pthread_t *threads;

    SEXP R_Lc,R_rms;

    Lc  = Calloc(T,double); // zeroes too!
    rms = Calloc(T,double); // zeroes too!
    Lcd = Calloc(T,pipf_Lc_struct);
    threads = Calloc(T,pthread_t);

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    pthread_attr_setschedpolicy(&attr, SCHED_RR);
//    pthread_attr_setscope(&attr, PTHREAD_SCOPE_PROCESS);
//    pthread_setconcurrency(T);

    iLc.data   = INTEGER(R_data);
    iLc.counts = INTEGER(R_counts);
    iLc.lams   = REAL(R_x);
    iLc.giks   = iLc.lams + K*L;
    iLc.alphas = iLc.giks + K*I;
    iLc.betas  = iLc.alphas + I;
    iLc.w      = *REAL(R_w);
    iLc.I = I;
    iLc.J = J;
    iLc.K = K;
    iLc.L = L;
    iLc.Lv = INTEGER(R_Lv);

    idx=0;len = (int) ceil((double)I/(double)T);
//    Rprintf("R_pipf_Lc len %d\n",len);
    for(t=0;t<T;t++) {
        Lcd[t].iLc = &iLc;
        Lcd[t].Lc  = &Lc[t];
        Lcd[t].rms = &rms[t];
        Lcd[t].i1 = idx; idx += len;
        Lcd[t].i2 = (idx < I) ? idx : I;
        pthread_create(&threads[t], &attr, pipf_Lc, (void *) &Lcd[t]);
    }
    for(t=0;t<T;t++)
        pthread_join(threads[t], (void **) &status);
    PROTECT(R_Lc  = allocVector(REALSXP,1));
    PROTECT(R_rms = allocVector(REALSXP,1));
    *REAL(R_Lc)  = 0.0F;
    *REAL(R_rms) = 0.0F;
    for(t=0;t<T;t++) {
        *REAL(R_Lc) += Lc[t];
        *REAL(R_rms) += rms[t];
    }
    *REAL(R_rms) = sqrt(*REAL(R_rms)/(I+I*J));

    setAttrib(R_Lc, install("rms"), R_rms);

    Free(Lc);
    Free(rms);
    Free(Lcd);
    Free(threads);
    UNPROTECT(2);
    return(R_Lc);
}

/* Lagrangian of the Likelihood function -- includes constraints etc. */
double ipf_Lc_raw(
    int *data, 
    int *counts, 
    double *lams, 
    double *giks, 
    double *alphas,
    double *betas, 
    double w,
    int I, 
    int J, 
    int K, 
    int L,
    int *Lv,
    double *rms) {

    int i,j,k,l,lsum;
    double tmp1,tmp2,Lc=0;

    *rms=0;
    for(i=0;i<I;i++) {
        R_CheckUserInterrupt();
        lsum=0;
        for(j=0;j<J;j++) {
            l = data[j*I+i];
            if(l >= 0) {
                tmp1=0;
                for(k=0;k<K;k++) {
                    tmp1 += giks[k*I+i] * lams[k*L+l];
                }
                Lc += counts[i] * tmp1;
            }
            tmp2=0;
            for(l=lsum;l<lsum+Lv[j];l++) {
                tmp1=0;
                for(k=0;k<K;k++) {
                    tmp1 += giks[k*I+i] * lams[k*L+l];
                }
                tmp2 += exp(-tmp1);
            }
            lsum += Lv[j];
            tmp2 -= 1;
            Lc += betas[j*I+i] * tmp2;
            tmp2 = pow(tmp2,2);
            Lc += w/2 * tmp2;
            *rms += tmp2;
        }
        tmp1=0;
        for(k=0;k<K;k++) {
            tmp1 += giks[k*I+i];
        }
        tmp1 -= 1;
        Lc += alphas[i] * tmp1;
        tmp1 = pow(tmp1,2);
        Lc += w/2 * tmp1;
        *rms += tmp1;
    }

    *rms = sqrt(*rms/(I+I*J));
    check_val(rms);
    check_val(&Lc);

    return(Lc);
}
/* Done this way to accomidate the optim command */
SEXP ipf_Lc(
            SEXP R_x, /* Contains Lams, Giks and Alphas, Betas */
            SEXP R_data, 
            SEXP R_counts, 
            SEXP R_w,
            SEXP R_K,
            SEXP R_Lv) {
    int I     = nrows(R_data);
    int J     = ncols(R_data);
    int K     = INTEGER(R_K)[0];
    int L     = isum(INTEGER(R_Lv),J);

    double *lams,*giks,*alphas,*betas;

    SEXP R_Lc,R_rms;

    PROTECT(R_Lc  = allocVector(REALSXP,1));
    PROTECT(R_rms = allocVector(REALSXP,1));

    lams   = REAL(R_x);
    giks   = lams + K*L;
    alphas = giks + K*I;
    betas  = alphas + I;

    REAL(R_Lc)[0] = ipf_Lc_raw(
        INTEGER(R_data),INTEGER(R_counts),
        lams,giks,alphas,betas,REAL(R_w)[0],
        I,J,K,L,INTEGER(R_Lv),REAL(R_rms));


    setAttrib(R_Lc, install("rms"), R_rms);

    UNPROTECT(2);
    return(R_Lc);
}

/* Done this way to accomidate the optim command */
/* The LLc functions do the same as the Lc functions but assume the lamdas are known and don't change */
SEXP ipf_LLc(
            SEXP R_x, /* Contains Giks, Alphas and Betas */
            SEXP R_data, 
            SEXP R_counts, 
            SEXP R_lams,
            SEXP R_w,
            SEXP R_Lv) {
    int I     = nrows(R_data);
    int J     = ncols(R_data);
    int K     = ncols(R_lams);
    int L     = isum(INTEGER(R_Lv),J);

    double *giks,*alphas,*betas;

    SEXP R_Lc,R_rms;

    PROTECT(R_Lc  = allocVector(REALSXP,1));
    PROTECT(R_rms = allocVector(REALSXP,1));

    giks   = REAL(R_x);
    alphas = giks + K*I;
    betas  = alphas + I;

    REAL(R_Lc)[0] = ipf_Lc_raw(
        INTEGER(R_data),INTEGER(R_counts),
        REAL(R_lams),giks,alphas,betas,REAL(R_w)[0],
        I,J,K,L,INTEGER(R_Lv),REAL(R_rms));

//    Rprintf("Lc %f rms %f\n",REAL(R_Lc)[0],REAL(R_rms)[0]);

    setAttrib(R_Lc, install("rms"), R_rms);

    UNPROTECT(2);
    return(R_Lc);
}

typedef struct {
    int *data, *counts;
    double *lams, *giks, *alphas, *betas;
    double w;
    int I,J,K,L,*Lv;
    double *dlams, *dgiks, *dalphas, *dbetas;
} ipf_dLc_struct;
typedef struct {
    ipf_dLc_struct *idLc;
    int i1,i2;
} pipf_dLc_struct;

void *pipf_dLc1(void *arg) {
    int lsum;
    int i,i0,j,j0,k,k0,l,l0;
    int I,J,K,L,*Lv;
    double tmp1,tmp2,tmp3,tmp4;

    pipf_dLc_struct *dLcd = arg;

    I = dLcd->idLc->I;
    J = dLcd->idLc->J;
    K = dLcd->idLc->K;
    L = dLcd->idLc->L;
    Lv = dLcd->idLc->Lv;

    for(i0=dLcd->i1;i0<dLcd->i2;i0++) {
        lsum=0;
        for(j0=0;j0<J;j0++) {
            tmp2=0;
            for(l=lsum;l<lsum+Lv[j0];l++) {
                tmp1=0;
                for(k=0;k<K;k++) {
                    tmp1 += dLcd->idLc->giks[k*I+i0] * dLcd->idLc->lams[k*L+l];
                }
                tmp2 += exp(-tmp1);
            }
            lsum += Lv[j0];
            dLcd->idLc->dbetas[j0*I+i0] = (tmp2 - 1);
        }
        tmp4=0;
        for(k0=0;k0<K;k0++) {
            dLcd->idLc->dgiks[k0*I+i0] = dLcd->idLc->alphas[i0];
            lsum=0;
            for(j=0;j<J;j++) {
                l = dLcd->idLc->data[j*I+i0];
                if(l != NA_INTEGER)
                    dLcd->idLc->dgiks[k0*I+i0] += dLcd->idLc->lams[k0*L+l] * dLcd->idLc->counts[i0];
                tmp2=tmp3=0;
                for(l=lsum;l<lsum+Lv[j];l++) {
                    tmp1=0;
                    for(k=0;k<K;k++) {
                        tmp1 += dLcd->idLc->giks[k*I+i0] * dLcd->idLc->lams[k*L+l];
                    }
                    tmp2 += dLcd->idLc->lams[k0*L+l] * exp(-tmp1);
                    tmp3 +=                exp(-tmp1);
                }
                lsum += Lv[j];
                dLcd->idLc->dgiks[k0*I+i0] -= dLcd->idLc->betas[j*I+i0] * tmp2;
                dLcd->idLc->dgiks[k0*I+i0] -= dLcd->idLc->w * (tmp3-1) * tmp2;
            } 
            tmp1=0;
            for(k=0;k<K;k++) {
                tmp1 += dLcd->idLc->giks[k*I+i0];
            }
            dLcd->idLc->dgiks[k0*I+i0] += dLcd->idLc->w*(tmp1-1);
            tmp4 += dLcd->idLc->giks[k0*I+i0];
        }
        dLcd->idLc->dalphas[i0] = (tmp4 - 1);
        sched_yield();
    }
    pthread_exit((void *) 0);
}
void *pipf_dLc2(void *arg) {
    int lsum;
    int i,j,j0,k,k0,l,l0;
    int I,J,K,L,*Lv;
    double tmp1,tmp2,tmp3;

    pipf_dLc_struct *dLcd = arg;

    I = dLcd->idLc->I;
    J = dLcd->idLc->J;
    K = dLcd->idLc->K;
    L = dLcd->idLc->L;
    Lv = dLcd->idLc->Lv;

    for(k0=dLcd->i1;k0<dLcd->i2;k0++) {
        lsum=0;
        for(j0=0;j0<J;j0++) {
            for(l0=lsum;l0<lsum+Lv[j0];l0++) {
                dLcd->idLc->dlams[k0*L+l0] = 0;
                for(i=0;i<I;i++) {
                    l = dLcd->idLc->data[j0*I+i];
                    if(l != NA_INTEGER && l0 == l)
                        dLcd->idLc->dlams[k0*L+l0] += dLcd->idLc->giks[k0*I+i] * 
                            dLcd->idLc->counts[i];
                    tmp2=0;
                    for(k=0;k<K;k++) {
                        tmp2 += dLcd->idLc->giks[k*I+i] * dLcd->idLc->lams[k*L+l0];
                    }
                    tmp2 = exp(-tmp2);
                    dLcd->idLc->dlams[k0*L+l0] -= dLcd->idLc->betas[j0*I+i] * 
                        dLcd->idLc->giks[k0*I+i] * tmp2;
                    tmp3=0;
                    for(l=lsum;l<lsum+Lv[j0];l++) {
                        tmp1=0;
                        for(k=0;k<K;k++) {
                            tmp1 += dLcd->idLc->giks[k*I+i] * dLcd->idLc->lams[k*L+l];
                        }
                        tmp3 += exp(-tmp1);
                    }
                    dLcd->idLc->dlams[k0*L+l0] -= dLcd->idLc->w * 
                        dLcd->idLc->giks[k0*I+i] * tmp2 * (tmp3-1);
                    sched_yield();
                }
            }
            lsum += Lv[j0];
        }
    }
    pthread_exit((void *) 0);
}

SEXP R_pipf_dLc(SEXP R_x, SEXP R_data, SEXP R_counts, SEXP R_w, SEXP R_K, SEXP R_Lv, SEXP R_threads) {
    int I     = nrows(R_data);
    int J     = ncols(R_data);
    int K     = INTEGER(R_K)[0];
    int L     = isum(INTEGER(R_Lv),J);
    int status,idx,len;
    int t,T1 = (*INTEGER(R_threads) < I) ? *INTEGER(R_threads) : I;
    int   T2 = (*INTEGER(R_threads) < K) ? *INTEGER(R_threads) : K;

    ipf_dLc_struct idLc;
    pipf_dLc_struct *dLcd;
    pthread_attr_t attr;
    pthread_t *threads;

    SEXP R_dlams,R_dgiks,R_dalphas,R_dbetas,R_list,R_names;
    PROTECT(R_dlams   = allocMatrix(REALSXP,L,K));
    PROTECT(R_dgiks   = allocMatrix(REALSXP,I,K));
    PROTECT(R_dalphas = allocMatrix(REALSXP,I,1));
    PROTECT(R_dbetas  = allocMatrix(REALSXP,I,J));
    dLcd    = Calloc(T1 + T2,pipf_dLc_struct);
    threads = Calloc(T1 + T2,pthread_t);

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    pthread_attr_setschedpolicy(&attr, SCHED_RR);
//    pthread_attr_setscope(&attr, PTHREAD_SCOPE_PROCESS);
//    pthread_setconcurrency(T1+T2);

    idLc.data   = INTEGER(R_data);
    idLc.counts = INTEGER(R_counts);
    idLc.lams   = REAL(R_x);
    idLc.giks   = idLc.lams + L*K;
    idLc.alphas = idLc.giks + I*K;
    idLc.betas  = idLc.alphas + I;
    idLc.w      = *REAL(R_w);
    idLc.I = I;
    idLc.J = J;
    idLc.K = K;
    idLc.L = L;
    idLc.Lv = INTEGER(R_Lv);
    idLc.dlams   = REAL(R_dlams);
    idLc.dgiks   = REAL(R_dgiks);
    idLc.dalphas = REAL(R_dalphas);
    idLc.dbetas  = REAL(R_dbetas);

    idx=0;len = (int) ceil((double)I/(double)T1);
//    Rprintf("R_pipf_dLc Ilen %d\n",len);
    for(t=0;t<T1;t++) {
        dLcd[t].idLc = &idLc;
        dLcd[t].i1 = idx;
        idx += len;
        dLcd[t].i2 = (idx < I) ? idx : I;
//        Rprintf("Launching %d\n",t);
        pthread_create(&threads[t], &attr, pipf_dLc1, (void *) &dLcd[t]);
    }
    idx=0;len = (int) ceil((double)K/(double)T2);
//    Rprintf("R_pipf_dLc Klen %d\n",len);
    for(t=T1;t<T1+T2;t++) {
        dLcd[t].idLc = &idLc;
        dLcd[t].i1 = idx;
        idx += len;
        dLcd[t].i2 = (idx < K) ? idx : K;
//        Rprintf("Launching %d\n",t);
        pthread_create(&threads[t], &attr, pipf_dLc2, (void *) &dLcd[t]);
    }
    for(t=0;t<T1+T2;t++) {
//        Rprintf("Collecting %d\n",t);
        pthread_join(threads[t], (void **) &status);
    }

//    Rprintf("Collected threads and returning.\n");
    PROTECT(R_list = allocVector(VECSXP,4));
    SET_VECTOR_ELT(R_list,0,R_dlams);
    SET_VECTOR_ELT(R_list,1,R_dgiks);
    SET_VECTOR_ELT(R_list,2,R_dalphas);
    SET_VECTOR_ELT(R_list,3,R_dbetas);

    PROTECT(R_names = allocVector(STRSXP,4));
    SET_STRING_ELT(R_names,0,mkChar("dlams"));
    SET_STRING_ELT(R_names,1,mkChar("dgiks"));
    SET_STRING_ELT(R_names,2,mkChar("dalphas"));
    SET_STRING_ELT(R_names,3,mkChar("dbetas"));

    setAttrib(R_list, R_NamesSymbol, R_names);

    Free(dLcd);
    Free(threads);
    UNPROTECT(6);
    return(R_list);
}

void ipf_dLc_raw(
        int *data, int *counts,
        double *lams, double *giks, double *alphas, double *betas, double w,
        int I, int J, int K, int L, int *Lv,
        double *dlams, double *dgiks, double *dalphas, double *dbetas) {

    int lsum;
    int i,i0,j,j0,k,k0,l,l0;
    double tmp1,tmp2,tmp3,tmp4;

    for(i0=0;i0<I;i0++) {
        R_CheckUserInterrupt();
        lsum=0;
        for(j0=0;j0<J;j0++) {
            tmp2=0;
            for(l=lsum;l<lsum+Lv[j0];l++) {
                tmp1=0;
                for(k=0;k<K;k++) {
                    tmp1 += giks[k*I+i0] * lams[k*L+l];
                }
                tmp2 += exp(-tmp1);
            }
            lsum += Lv[j0];

            dbetas[j0*I+i0] = (tmp2 - 1);
            check_val(dbetas + j0*I+i0);
        }
        tmp4=0;
        for(k0=0;k0<K;k0++) {
            dgiks[k0*I+i0] = alphas[i0];
            lsum=0;
            for(j=0;j<J;j++) {
                l = data[j*I+i0];
                if(l >= 0)
                    dgiks[k0*I+i0] += lams[k0*L+l] * counts[i0];
                tmp2=tmp3=0;
                for(l=lsum;l<lsum+Lv[j];l++) {
                    tmp1=0;
                    for(k=0;k<K;k++) {
                        tmp1 += giks[k*I+i0] * lams[k*L+l];
                    }
                    tmp2 += lams[k0*L+l] * exp(-tmp1);
                    tmp3 +=                exp(-tmp1);
                }
                lsum += Lv[j];

                dgiks[k0*I+i0] -= betas[j*I+i0] * tmp2;
                dgiks[k0*I+i0] -= w * (tmp3-1) * tmp2;
            } 
            tmp1=0;
            for(k=0;k<K;k++) {
                tmp1 += giks[k*I+i0];
            }
            dgiks[k0*I+i0] += w*(tmp1-1);

            check_val(dgiks + k0*I+i0);

            tmp4 += giks[k0*I+i0];
        }
        dalphas[i0] = (tmp4 - 1);
        check_val(dalphas + i0);
    }
    for(k0=0;k0<K;k0++) {
        R_CheckUserInterrupt();
        lsum=0;
        for(j0=0;j0<J;j0++) {
            for(l0=lsum;l0<lsum+Lv[j0];l0++) {
                dlams[k0*L+l0] = 0;
                for(i=0;i<I;i++) {
                    l = data[j0*I+i];
                    if(l >= 0 && l0 == l)
                        dlams[k0*L+l0] += giks[k0*I+i] * counts[i];

                    tmp2=0;
                    for(k=0;k<K;k++) {
                        tmp2 += giks[k*I+i] * lams[k*L+l0];
                    }
                    tmp2 = exp(-tmp2);
                    dlams[k0*L+l0] -= betas[j0*I+i] * giks[k0*I+i] * tmp2;

                    tmp3=0;
                    for(l=lsum;l<lsum+Lv[j0];l++) {
                        tmp1=0;
                        for(k=0;k<K;k++) {
                            tmp1 += giks[k*I+i] * lams[k*L+l];
                        }
                        tmp3 += exp(-tmp1);
                    }

                    dlams[k0*L+l0] -= w * giks[k0*I+i] * tmp2 * (tmp3-1);
                }
                check_val(dlams + k0*L+l0);
            }
            lsum += Lv[j0];
        }
    }
}
SEXP ipf_dLc(SEXP R_x, SEXP R_data, SEXP R_counts, SEXP R_w, SEXP R_K, SEXP R_Lv) {
    int I     = nrows(R_data);
    int J     = ncols(R_data);
    int K     = INTEGER(R_K)[0];
    int L     = isum(INTEGER(R_Lv),J);

    double *lams,*giks,*alphas,*betas;

    SEXP R_dlams,R_dgiks,R_dalphas,R_dbetas,R_list,R_names;

    lams   = REAL(R_x);
    giks   = lams + L*K;
    alphas = giks + I*K;
    betas  = alphas + I;

    PROTECT(R_dlams   = allocMatrix(REALSXP,L,K));
    PROTECT(R_dgiks   = allocMatrix(REALSXP,I,K));
    PROTECT(R_dalphas = allocMatrix(REALSXP,I,1));
    PROTECT(R_dbetas  = allocMatrix(REALSXP,I,J));

    ipf_dLc_raw(INTEGER(R_data),INTEGER(R_counts),
        lams,giks,alphas,betas,REAL(R_w)[0],
        I,J,K,L,INTEGER(R_Lv),
        REAL(R_dlams),REAL(R_dgiks),REAL(R_dalphas),REAL(R_dbetas));

    PROTECT(R_list = allocVector(VECSXP,4));
    SET_VECTOR_ELT(R_list,0,R_dlams);
    SET_VECTOR_ELT(R_list,1,R_dgiks);
    SET_VECTOR_ELT(R_list,2,R_dalphas);
    SET_VECTOR_ELT(R_list,3,R_dbetas);

    PROTECT(R_names = allocVector(STRSXP,4));
    SET_STRING_ELT(R_names,0,mkChar("dlams"));
    SET_STRING_ELT(R_names,1,mkChar("dgiks"));
    SET_STRING_ELT(R_names,2,mkChar("dalphas"));
    SET_STRING_ELT(R_names,3,mkChar("dbetas"));

    setAttrib(R_list, R_NamesSymbol, R_names);

    UNPROTECT(6);
    return(R_list);
}

void ipf_dLLc_raw(  /* Just like above, but lams don't change -- have lams, trying to estimate giks  */
        int *data, int *counts,
        double *lams, double *giks, double *alphas, double *betas, double w,
        int I, int J, int K, int L, int *Lv,
        double *dgiks, double *dalphas, double *dbetas) {

    int lsum;
    int i0,j,j0,k,k0,l;
    double tmp1,tmp2,tmp3,tmp4;

    for(i0=0;i0<I;i0++) {
        lsum=0;
        for(j0=0;j0<J;j0++) {
            tmp2=0;
            for(l=lsum;l<lsum+Lv[j0];l++) {
                tmp1=0;
                for(k=0;k<K;k++) {
                    tmp1 += giks[k*I+i0] * lams[k*L+l];
                }
                tmp2 += exp(-tmp1);
            }
            lsum += Lv[j0];

            dbetas[j0*I+i0] = (tmp2 - 1);
            check_val(dbetas + j0*I+i0);
        }
        tmp4=0;
        for(k0=0;k0<K;k0++) {
            dgiks[k0*I+i0] = alphas[i0];
            lsum=0;
            for(j=0;j<J;j++) {
                l = data[j*I+i0];
                if(l >= 0)
                    dgiks[k0*I+i0] += lams[k0*L+l] * counts[i0];
                tmp2=tmp3=0;
                for(l=lsum;l<lsum+Lv[j];l++) {
                    tmp1=0;
                    for(k=0;k<K;k++) {
                        tmp1 += giks[k*I+i0] * lams[k*L+l];
                    }
                    tmp2 += lams[k0*L+l] * exp(-tmp1);
                    tmp3 +=                exp(-tmp1);
                }
                lsum += Lv[j];

                dgiks[k0*I+i0] -= betas[j*I+i0] * tmp2;
                dgiks[k0*I+i0] -= w * (tmp3-1) * tmp2;
            } 
            tmp1=0;
            for(k=0;k<K;k++) {
                tmp1 += giks[k*I+i0];
            }
            dgiks[k0*I+i0] += w*(tmp1-1);

            check_val(dgiks + k0*I+i0);

            tmp4 += giks[k0*I+i0];
        }
        dalphas[i0] = (tmp4 - 1);
        check_val(dalphas + i0);
    }
}
SEXP ipf_dLLc(SEXP R_x, SEXP R_data, SEXP R_counts, SEXP R_lams, SEXP R_w, SEXP R_Lv) {
    int I     = nrows(R_data);
    int J     = ncols(R_data);
    int K     = ncols(R_lams);
    int L     = isum(INTEGER(R_Lv),J);

    double *giks,*alphas,*betas;

    SEXP R_dgiks,R_dalphas,R_dbetas,R_list,R_names;

    giks   = REAL(R_x);
    alphas = giks + I*K;
    betas  = alphas + I;

    PROTECT(R_dgiks   = allocMatrix(REALSXP,I,K));
    PROTECT(R_dalphas = allocMatrix(REALSXP,I,1));
    PROTECT(R_dbetas  = allocMatrix(REALSXP,I,J));

    ipf_dLLc_raw(INTEGER(R_data),INTEGER(R_counts),
        REAL(R_lams),giks,alphas,betas,REAL(R_w)[0],
        I,J,K,L,INTEGER(R_Lv),
        REAL(R_dgiks),REAL(R_dalphas),REAL(R_dbetas));

    PROTECT(R_list = allocVector(VECSXP,3));
    SET_VECTOR_ELT(R_list,0,R_dgiks);
    SET_VECTOR_ELT(R_list,1,R_dalphas);
    SET_VECTOR_ELT(R_list,2,R_dbetas);

    PROTECT(R_names = allocVector(STRSXP,3));
    SET_STRING_ELT(R_names,0,mkChar("dgiks"));
    SET_STRING_ELT(R_names,1,mkChar("dalphas"));
    SET_STRING_ELT(R_names,2,mkChar("dbetas"));

    setAttrib(R_list, R_NamesSymbol, R_names);

    UNPROTECT(5);
    return(R_list);
}

SEXP gik_state(SEXP giks, SEXP n, SEXP sep, SEXP sym) {
    int i,I = nrows(giks);
    int k,K = ncols(giks);
    int depth = INTEGER(n)[0];
    int maxk;
    double tmpd;
    char *mysep = (char *)CHAR(STRING_ELT(sep,0));
    char *s1,*s2;

    SEXP state,g;

    PROTECT(state = allocVector(STRSXP,I));
    PROTECT(g     = allocMatrix(REALSXP,I,K));

    for(i=0;i<I;i++) {
        for(k=0;k<K;k++) {
            REAL(g)[k*I + i] = REAL(giks)[k*I + i];
        }
    }

    for(i=0;i<I;i++) {
        tmpd = REAL(g)[0*I + i];
        maxk = 0;
        for(k=1;k<K;k++) {
            if(REAL(g)[k*I + i] > tmpd) {
                tmpd = REAL(g)[k*I + i];
                maxk = k;
            }
        }
        SET_STRING_ELT(state,i,mkChar(CHAR(STRING_ELT(sym,maxk))));
        REAL(g)[maxk*I + i] /= K;
    }

    while(--depth >= 1) {
        for(i=0;i<I;i++) {
            tmpd = REAL(g)[0*I + i];
            maxk = 0;
            for(k=1;k<K;k++) {
                if(REAL(g)[k*I + i] > tmpd) {
                    tmpd = REAL(g)[k*I + i];
                    maxk = k;
                }
            }
            s1 = (char *)CHAR(STRING_ELT(state,i));
            s2 = R_alloc(strlen(s1)+strlen(mysep)+strlen(CHAR(STRING_ELT(sym,maxk)))+1,sizeof(char));
            sprintf(s2,"%s%s%s",s1,mysep,CHAR(STRING_ELT(sym,maxk)));
            SET_STRING_ELT(state,i,mkChar(s2));
            REAL(g)[maxk*I + i] /= K;
        }
    }
    UNPROTECT(2);
    return(state);
}

SEXP lams_gom(SEXP data, SEXP lams, SEXP giks, SEXP tol, SEXP max_iter) {
    int i,I     = nrows(data);
    int j,J     = ncols(data);
    int k1,k2,K = ncols(giks);
    int l1,l2,L;
    int iter=0;
    SEXP lamM,rv;
    double tmp1,tmp2,tmp3,tmp4;
    double err,newlam;

    do {
        iter++;
        err = 0.0;
        for(j=0;j<J;j++) {
            R_FlushConsole();
//            R_ProcessEvents();
            lamM = VECTOR_ELT(lams,j);
            L    = nrows(lamM);
            for(l1=0;l1<L;l1++) {
                for(k1=0;k1<K;k1++) {
                    tmp1=tmp3=0;
                    for(i=0;i<I;i++) {
                        l2 = INTEGER(data)[j*I + i];
                        if(l2 >= 0) {
                            tmp2=tmp4=0;
                            for(k2=0;k2<K;k2++) {
                                tmp2 += REAL(giks)[k2*I + i] * REAL(lamM)[k2*L + l1];
                                tmp4 += REAL(giks)[k2*I + i] * REAL(lamM)[k2*L + l2];
                            }
                                         tmp1 +=  REAL(giks)[k1*I + i] * REAL(lamM)[k1*L + l2] / tmp4;
                            if(l1 == l2) tmp3 +=  REAL(giks)[k1*I + i] * REAL(lamM)[k1*L + l1] / tmp2;
                        }
                    }
                    newlam = tmp3/tmp1;
                    err += pow(newlam - REAL(lamM)[k1*L + l1],2);
                    REAL(lamM)[k1*L + l1] = newlam;
                }
            }
        }
    } while((REAL(tol)[0] > 0 && err > REAL(tol)[0]) ||
            (INTEGER(max_iter)[0] > 0 && iter < INTEGER(max_iter)[0]));

    if(REAL(tol)[0] > 0 ) {
        PROTECT(rv = allocVector(INTSXP, 1));
        INTEGER(rv)[0] = iter;
    } else {
        PROTECT(rv = allocVector(REALSXP, 1));
        REAL(rv)[0] = err;
    }

    UNPROTECT(1);
    return(rv);
}

SEXP lams_igom(SEXP data, SEXP lams, SEXP giks, SEXP tol, SEXP max_iter, SEXP minfo) {
    int i,I     = nrows(data);
    int j,J     = ncols(data);
    int k1,k2,K = ncols(giks);
    int l1,l2,L;
    int iter=0;
    SEXP lamM,rv;
    double tmp1,tmp2,tmp3,tmp4;
    double err,newlam;

    do {
        iter++;
        err = 0.0;
        for(j=0;j<J;j++) {
            R_FlushConsole();
//            R_ProcessEvents();
            if(REAL(minfo)[j]) {   /* skip it if there's no information */
                lamM = VECTOR_ELT(lams,j);
                L    = nrows(lamM);
                for(l1=0;l1<L;l1++) {
                    for(k1=0;k1<K;k1++) {
                        tmp1=tmp3=0;
                        for(i=0;i<I;i++) {
                            l2 = INTEGER(data)[j*I + i];
                            if(l2 >= 0) {
                                tmp2=tmp4=0;
                                for(k2=0;k2<K;k2++) {
                                    tmp2 += REAL(giks)[k2*I + i] * REAL(lamM)[k2*L + l1];
                                    tmp4 += REAL(giks)[k2*I + i] * REAL(lamM)[k2*L + l2];
                                }
                                             tmp1 +=  REAL(minfo)[j] * REAL(giks)[k1*I + i] * REAL(lamM)[k1*L + l2] / tmp4;
                                if(l1 == l2) tmp3 +=  REAL(minfo)[j] * REAL(giks)[k1*I + i] * REAL(lamM)[k1*L + l1] / tmp2;
                            }
                        }
                        newlam = tmp3/tmp1;
                        err += pow(newlam - REAL(lamM)[k1*L + l1],2);
                        REAL(lamM)[k1*L + l1] = newlam;
                    }
                }
            }
        }
    } while((REAL(tol)[0] > 0 && err > REAL(tol)[0]) ||
            (INTEGER(max_iter)[0] > 0 && iter < INTEGER(max_iter)[0]));

    if(REAL(tol)[0] > 0 ) {
        PROTECT(rv = allocVector(INTSXP, 1));
        INTEGER(rv)[0] = iter;
    } else {
        PROTECT(rv = allocVector(REALSXP, 1));
        REAL(rv)[0] = err;
    }
    UNPROTECT(1);
    return(rv);
}

SEXP giks_gom(SEXP data, SEXP lams, SEXP giks, SEXP tol, SEXP max_iter) {
    int i,I     = nrows(data);
    int j,J     = ncols(data);
    int k1,k2,K = ncols(giks);
    int l,L;
    int iter=0;
    SEXP lamM,rv;
    double tmp1,tmp2,tot_J;
    double err,newgik;

    do {
        iter++;
        err = 0.0;
        for(i=0;i<I;i++) {
            R_FlushConsole();
//            R_ProcessEvents();
            for(k1=0;k1<K;k1++) {
                tmp1=tot_J=0;
                for(j=0;j<J;j++) {
                    l    = INTEGER(data)[j*I + i];
                    if(l >= 0) {
                        tot_J++;
                        lamM = VECTOR_ELT(lams,j);
                        L    = nrows(lamM);
                        tmp2=0;
                        for(k2=0;k2<K;k2++) {
                            tmp2 += REAL(giks)[k2*I + i] * REAL(lamM)[k2*L + l];
                        }
                        tmp1 += REAL(giks)[k1*I + i] * REAL(lamM)[k1*L + l] / tmp2;
                    }
                }
                newgik = tmp1/tot_J;
                err += pow(newgik - REAL(giks)[k1*I + i],2);
                REAL(giks)[k1*I + i] = newgik;
            }
        }
    } while((REAL(tol)[0] > 0 && err > REAL(tol)[0]) ||
            (INTEGER(max_iter)[0] > 0 && iter < INTEGER(max_iter)[0]));

    if(REAL(tol)[0] > 0 ) {
        PROTECT(rv = allocVector(INTSXP, 1));
        INTEGER(rv)[0] = iter;
    } else {
        PROTECT(rv = allocVector(REALSXP, 1));
        REAL(rv)[0] = err;
    }
    UNPROTECT(1);
    return(rv);
}
SEXP giks_igom(SEXP data, SEXP lams, SEXP giks, SEXP tol, SEXP max_iter, SEXP minfo) {
    int i,I     = nrows(data);
    int j,J     = ncols(data);
    int k1,k2,K = ncols(giks);
    int l,L;
    int iter=0;
    SEXP lamM,rv;
    double tmp1,tmp2,tot_J;
    double err,newgik;

    do {
        iter++;
        err = 0.0;
        for(i=0;i<I;i++) {
            R_FlushConsole();
//            R_ProcessEvents();
            for(k1=0;k1<K;k1++) {
                tmp1=tot_J=0;
                for(j=0;j<J;j++) {
                    l    = INTEGER(data)[j*I + i];
                    if(l >= 0) {
                        tot_J += REAL(minfo)[j];
                        lamM = VECTOR_ELT(lams,j);
                        L    = nrows(lamM);
                        tmp2=0;
                        for(k2=0;k2<K;k2++) {
                            tmp2 += REAL(giks)[k2*I + i] * REAL(lamM)[k2*L + l];
                        }
                        tmp1 += REAL(minfo)[j] * REAL(giks)[k1*I + i] * REAL(lamM)[k1*L + l] / tmp2;
                    }
                }
                newgik = tmp1/tot_J;
                err += pow(newgik - REAL(giks)[k1*I + i],2);
                REAL(giks)[k1*I + i] = newgik;
            }
        }
    } while((REAL(tol)[0] > 0 && err > REAL(tol)[0]) ||
            (INTEGER(max_iter)[0] > 0 && iter < INTEGER(max_iter)[0]));

    if(REAL(tol)[0] > 0 ) {
        PROTECT(rv = allocVector(INTSXP, 1));
        INTEGER(rv)[0] = iter;
    } else {
        PROTECT(rv = allocVector(REALSXP, 1));
        REAL(rv)[0] = err;
    }
    UNPROTECT(1);
    return(rv);
}

SEXP get_giks(int node_idx, int rec_idx, SEXP rho, int debug) {
    int i,I;
    char *var_name,*var_value,*p=NULL;
    SEXP R_node,R_tmp;

    PROTECT(R_node  = VECTOR_ELT(findVar(install("gtree"),rho),node_idx-1));
    if(R_node == R_NilValue) error("Bad gtree node");
//    var_name  = CHAR(STRING_ELT(coerceVector(getListElement(R_node,"var.name"),STRSXP),0));
    var_name  = (char *)CHAR(STRING_ELT(coerceVector(VECTOR_ELT(R_node,0),STRSXP),0));

    PROTECT(R_tmp = getListElement(findVar(install("rec"),rho), var_name));
    if(!strcmp(CHAR(STRING_ELT(getAttrib(R_tmp, install("class")),0)),"factor")) {
        i = INTEGER(R_tmp)[rec_idx];
        var_value = (char *)CHAR(STRING_ELT(getAttrib(R_tmp,install("levels")),i-1));
    } else
        var_value = (char *)CHAR(STRING_ELT(coerceVector(R_tmp,STRSXP),rec_idx));
    UNPROTECT(1);

//    PROTECT(R_tmp = getListElement(R_node,"values"));
    PROTECT(R_tmp = VECTOR_ELT(R_node,3));
    I = length(R_tmp);
    for(i=0;i<I;i++)
        if(strspn(
                p = (char *)CHAR(STRING_ELT(coerceVector(R_tmp,STRSXP),i)),
                var_value
            ) == strlen(var_value)) {
            p += strlen(var_value) + 1;
            break;
        }
    UNPROTECT(1);

    if(debug) Rprintf("node %d var %s val %s I %d i %d p %s\n",
        node_idx,var_name,var_value,I,i,p);

    UNPROTECT(1);
    if(i >= I || !strcmp(p,"**")) /* Nothing found or the end of the road */
//        return(getListElement(R_node,"giks"));
        return(VECTOR_ELT(R_node,5));
    else
        if(strspn(p,"*"))
//            return(getListElement(VECTOR_ELT(findVar(install("gleaf"),rho),atoi(p+1)-1),"giks"));
            return(VECTOR_ELT(VECTOR_ELT(findVar(install("gleaf"),rho),atoi(p+1)-1),1));
        else
            return(get_giks(atoi(p),rec_idx,rho,debug));
}

SEXP R_get_giks(SEXP R_node_idx, SEXP rho, SEXP R_debug) {
    return(get_giks(
        INTEGER(coerceVector(R_node_idx,INTSXP))[0], 
        0,
        rho, 
        INTEGER(coerceVector(R_debug   ,INTSXP))[0]
   ));
}

SEXP R_gikker(SEXP RI, SEXP RK, SEXP Rtrunk, SEXP Renv, SEXP Rdebug) {
    int i,I = INTEGER(RI)[0];
    int k,K = INTEGER(RK)[0];

    SEXP Rgiks,Rtmp;

    if(INTEGER(Rdebug)[0]) Rprintf("I %d K %d\n",I,K);
    PROTECT(Rgiks = allocMatrix(REALSXP,I,K));
    for(i=0;i<I;i++) {
        if(INTEGER(Rdebug)[0]) Rprintf("gikker for %d\n",i);
        PROTECT(Rtmp = get_giks(INTEGER(Rtrunk)[0],i,Renv,INTEGER(Rdebug)[0]));
        for(k=0;k<K;k++) {
            REAL(Rgiks)[k*I+i] = REAL(Rtmp)[k];
        }
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return(Rgiks);
}
SEXP R_gikker2(SEXP Ripf, SEXP Rrecs, SEXP Rdebug) {
    int i,I = length(VECTOR_ELT(Rrecs,0));
    int k,K = ncols(getListElement(Ripf,"giks"));
    int j,J,node_idx;
    char *p=NULL,*var_name,*var_value;
    char ntype,go,debug = INTEGER(Rdebug)[0];

    SEXP Rgiks,Rgik=R_NilValue,R_trunk_node,R_node,R_var,R_gtree,R_gleaf,R_val,R_attr,
        Rnode,Rntype,Rtypes,Rnidx,Rclass,Rnames;

    PROTECT(Rgiks = allocMatrix(REALSXP,I,K));

    PROTECT(Rnidx = allocVector(INTSXP,I));
    PROTECT(Rntype = allocVector(INTSXP,I));

    PROTECT(R_gtree = getListElement(getListElement(Ripf,"m"),"gtree"));
    PROTECT(R_gleaf = getListElement(getListElement(Ripf,"m"),"gleaf"));
    node_idx = length(R_gtree)-1;
    PROTECT(R_trunk_node = VECTOR_ELT(R_gtree,node_idx));
    if(debug) Rprintf("I %d K %d root_idx %d\n",I,K,node_idx);

    for(i=0;i<I;i++) {
        R_node = R_trunk_node;
        node_idx = length(R_gtree)-1;
        go = 1;
        ntype = 3;
        while(go) {
            var_name = (char *)CHAR(STRING_ELT(VECTOR_ELT(R_node,0),0)); // var.name
            R_var = getListElement(Rrecs, var_name);
            R_attr = getAttrib(R_var, install("class"));
            if(R_attr != R_NilValue && strcmp(CHAR(STRING_ELT(R_attr,0)),"factor") == 0) {
                j = INTEGER(R_var)[i];
                var_value = (char *)CHAR(STRING_ELT(getAttrib(R_var,install("levels")),j-1));
            } else {
                var_value = (char *)CHAR(STRING_ELT(coerceVector(R_var,STRSXP),i));
            }
            R_val = VECTOR_ELT(R_node,3); // values
            J = length(R_val);
            for(j=0;j<J;j++) {
                if(strspn(p = (char *)CHAR(STRING_ELT(R_val,j)), var_value) == strlen(var_value)) {
                    p += strlen(var_value) + 1;
                    break;
                }
            }    
            if(debug) Rprintf("i %d node %d var %s val %s J %d j %d p %s\n",
                i, node_idx, var_name, var_value, J, j, p);

            if(j >= J || !strcmp(p,"**")) { // nothing found or end of road
                go = 0; ntype = 2; Rgik = VECTOR_ELT(R_node,5);
            } else if(strspn(p,"*")) { // leaf node -- good!
                go = 0; ntype = 1; Rgik = VECTOR_ELT(VECTOR_ELT(R_gleaf,atoi(p+1)-1),1);
            } else { // else keep going!
                node_idx = atoi(p)-1;
                R_node = VECTOR_ELT(R_gtree,node_idx);
            }
        }
        INTEGER(Rntype)[i] = ntype;
        INTEGER(Rnidx)[i] = node_idx;
        for(k=0;k<K;k++) {
            REAL(Rgiks)[k*I+i] = REAL(Rgik)[k];
        }
    }
    PROTECT(Rtypes = allocVector(STRSXP,2));
    SET_STRING_ELT(Rtypes,0,mkChar("gleaf"));
    SET_STRING_ELT(Rtypes,1,mkChar("gtree"));
    setAttrib(Rntype, install("levels"), Rtypes);

    PROTECT(Rclass = allocVector(STRSXP,1));
    SET_STRING_ELT(Rclass,0,mkChar("factor"));
    classgets(Rntype, Rclass);
    PROTECT(Rnode = allocVector(VECSXP,2));
    SET_VECTOR_ELT(Rnode,0,Rnidx);
    SET_VECTOR_ELT(Rnode,1,Rntype);

    PROTECT(Rnames = allocVector(STRSXP,2));
    SET_STRING_ELT(Rnames,0,mkChar("index"));
    SET_STRING_ELT(Rnames,1,mkChar("type"));
    setAttrib(Rnode, R_NamesSymbol, Rnames);

    setAttrib(Rgiks, install("nodes"), Rnode);

    UNPROTECT(10);
    return(Rgiks);
}

/* Calculates distance from each row of X to each row of Y.
    returns vector of distance and indexes of Y's closest.
    skips columns of X that conatain NULL or NA.
*/
SEXP R_closest(SEXP X, SEXP Y) {
    int xi,xI = nrows(X), j, J = ncols(X);
    int yi,yI = nrows(Y);
    int *idx;
    double *dist;
    double dtmp1,dtmp2;

    SEXP R_list, R_tmp1, R_tmp2, RX, RY;

    if(J != ncols(Y)) error("mydist: %d != %d Number of columns not equal in X and Y\n",J,ncols(Y));

    dist = Calloc(yI,double);
    idx  = Calloc(yI,int);
    PROTECT(RX = coerceVector(X,REALSXP));
    PROTECT(RY = coerceVector(Y,REALSXP));

    PROTECT(R_list = allocVector(VECSXP,xI));

    for(xi=0;xi<xI;xi++) {
        for(yi=0;yi<yI;yi++) {
            idx[yi] = yi+1;
            dist[yi] = 0.0F;
            for(j=0;j<J;j++) {
                dtmp1 = REAL(RX)[j*xI+xi];
                if(R_FINITE(dtmp1)) {
                    dtmp2 = dtmp1 - REAL(RY)[j*yI+yi];
                    dist[yi] += dtmp2 * dtmp2;
                } 
            }
        }
        rsort_with_index(dist,idx,yI);
        j=1; while(dist[0] == dist[j]) j++;
        PROTECT(R_tmp1 = allocVector(REALSXP,j));
        for(yi=0;yi<j;yi++)
            REAL(R_tmp1)[yi] = sqrt(dist[yi]);

        PROTECT(R_tmp2 = allocVector(INTSXP,j));
        memcpy(INTEGER(R_tmp2),idx,j * sizeof(int));

        setAttrib(R_tmp1, R_NamesSymbol, R_tmp2);
        
        SET_VECTOR_ELT(R_list,xi,R_tmp1);
        UNPROTECT(2);
    }
    UNPROTECT(3);
    Free(dist);
    Free(idx);
    return(R_list);
}

SEXP R_postP(SEXP R_P, SEXP R_X) {
    int   M = nrows(R_X), n,N = ncols(R_X);
    int i,I = nrows(R_P), j,J = ncols(R_P);
    int *X;
    double *postP,*P;

    SEXP R_postP, R_P2, R_X2;

    if(I != M) error("Bad data I %d M %d\n",I,M);

    PROTECT(R_postP = allocVector(REALSXP,I));
    PROTECT(R_X2 = coerceVector(R_X,INTSXP));
    PROTECT(R_P2 = coerceVector(R_P, REALSXP));
    postP = REAL(R_postP);
    X = INTEGER(R_X2);
    P = REAL(R_P2);

    for(i=0;i<I;i++) {
        R_CheckUserInterrupt();
        postP[i] = 1.0F;
        for(n=0;n<N;n++) {
            j = X[n*I+i];
            if(j < 0 || j >= J) error("Bad value at %d %d in data j %d J %d\n",i,n,j,J);
            postP[i] *= P[j*I+i];
        }
    }

    
    UNPROTECT(3);
    return(R_postP);
}

SEXP R_postP2(SEXP R_P, SEXP R_X, SEXP R_States) {
    int m,M = nrows(R_X), n,N = ncols(R_X);
    int i,I = nrows(R_P), j,J = ncols(R_P);
    int s,S = INTEGER(R_States)[0];
    int *X;
    double tmpP,*P;

    SEXP R_I, R_S, R_postP, R_list, R_names, R_P2, R_X2;

    PROTECT(R_I = allocVector(INTSXP,M));
    PROTECT(R_S = allocVector(INTSXP,M));
    PROTECT(R_postP = allocVector(REALSXP,M));

    PROTECT(R_X2 = coerceVector(R_X,INTSXP));  X = INTEGER(R_X2);
    PROTECT(R_P2 = coerceVector(R_P,REALSXP)); P = REAL(R_P2);

    for(m=0;m<M;m++) {
        R_CheckUserInterrupt();
        REAL(R_postP)[m]=0.0F;
        for(s=0;s<S;s++) {
            for(i=0;i<I;i++) {
                tmpP = 1.0F;
                for(n=0;n<N;n++) {
                    j = X[n*M+m]+s;
                    if(j < 0 || j >= J) 
                        error("Bad value in data m %d s %d n %d j %d J %d\n",m,s,n,j,J);
                    tmpP *= P[j*I+i];
                }
                if(REAL(R_postP)[m] < tmpP) {
                    INTEGER(R_I)[m] = i;
                    INTEGER(R_S)[m] = s;
                    REAL(R_postP)[m] = tmpP;
                }
            }
        }
        INTEGER(R_I)[m] += 1; // Make it R compatible
        INTEGER(R_S)[m] += 1;
    }

    PROTECT(R_list  = allocVector(VECSXP,3));
    SET_VECTOR_ELT(R_list,0,R_I);
    SET_VECTOR_ELT(R_list,1,R_S);
    SET_VECTOR_ELT(R_list,2,R_postP);

    PROTECT(R_names = allocVector(STRSXP,3));
    SET_STRING_ELT(R_names,0,mkChar("I"));
    SET_STRING_ELT(R_names,1,mkChar("S"));
    SET_STRING_ELT(R_names,2,mkChar("postP"));

    setAttrib(R_list, R_NamesSymbol, R_names);
    
    UNPROTECT(7);
    return(R_list);
}

typedef struct {
    double *eta, *lams;
    int *svar,*Y;
    int I,J,K,L1,L2;
    double *P;
} wazuni_struct;
typedef struct {
    wazuni_struct *w;
    int i1,i2;
} pwazuni_struct;

void *pwazuni(void *arg) {
    int i,I,j,J,k1,k2,K,l,L1,L2;
    double *dtmp1,dtmp2,dtmp3;

    pwazuni_struct *pw = arg;

    I = pw->w->I;
    J = pw->w->J;
    K = pw->w->K;
    L1 = pw->w->L1;
    L2 = pw->w->L2;

    dtmp1 = Calloc(K,double);

    for(i=pw->i1;i<pw->i2;i++) {
        dtmp3 = 0.0F;
        for(j=0;j<J;j++) {
            pw->w->P[j*I+i] = 0.0F;
            for(k1=0;k1<K;k1++) {
                dtmp1[k1] = pw->w->eta[k1] + pw->w->lams[k1*L1 + pw->w->svar[j]];
                for(l=0;l<L2;l++) if(pw->w->Y[l*I + i] != NA_INTEGER)
                        dtmp1[k1] += pw->w->lams[k1*L1 + pw->w->Y[l*I + i]];
            }
            for(k1=0;k1<K;k1++) {
                dtmp2 = 1.0F;
                for(k2=0;k2<K;k2++) if(k2 != k1) dtmp2 *= (dtmp1[k2] - dtmp1[k1]);
                pw->w->P[j*I+i] += exp(-dtmp1[k1])/dtmp2;
            }
            dtmp3 += pw->w->P[j*I+i];
        }
        for(j=0;j<J;j++) pw->w->P[j*I+i] /= dtmp3;
        sched_yield();
    }
    Free(dtmp1);
    pthread_exit((void *) 0);
}

SEXP R_pwazuni(SEXP R_svar, SEXP R_Y, SEXP R_eta, SEXP R_lams, SEXP R_threads) {
    int I = nrows(R_Y);
    int J = length(R_svar);
    int K = ncols(R_lams), L1 = nrows(R_lams);
    int L2 = ncols(R_Y);
    int status,idx,len;
    int t,T = (*INTEGER(R_threads) < I) ? *INTEGER(R_threads) : I;

    wazuni_struct w;
    pwazuni_struct *pw;
    pthread_attr_t attr;
    pthread_t *threads;

    SEXP R_P, R_dimnames;

    PROTECT(R_P = allocMatrix(REALSXP, I, J));
    PROTECT(R_dimnames = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(R_dimnames, 0, VECTOR_ELT(getAttrib(R_Y, R_DimNamesSymbol),0));
    SET_VECTOR_ELT(R_dimnames, 1, getAttrib(R_svar, R_NamesSymbol));
    setAttrib(R_P, R_DimNamesSymbol, R_dimnames);
    pw      = Calloc(T,pwazuni_struct);
    threads = Calloc(T,pthread_t);

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    pthread_attr_setschedpolicy(&attr, SCHED_RR);
//    pthread_attr_setscope(&attr, PTHREAD_SCOPE_PROCESS);
//    pthread_setconcurrency(T);

    w.eta  = REAL(R_eta);
    w.lams = REAL(R_lams);
    w.svar = INTEGER(R_svar);
    w.Y    = INTEGER(R_Y);
    w.I    = I;
    w.J    = J;
    w.K    = K;
    w.L1   = L1;
    w.L2   = L2;
    w.P    = REAL(R_P);

    idx=0;len = (int) ceil((double)I/(double)T); 
    for(t=0;t<T;t++) {
        pw[t].w = &w;
        pw[t].i1 = idx;
        idx += len;
        pw[t].i2 = (idx < I) ? idx : I;
        pthread_create(&threads[t], &attr, pwazuni, (void *) &pw[t]);
    }
    for(t=0;t<T;t++) {
        pthread_join(threads[t],(void **) &status);
    }

    UNPROTECT(2);
    return(R_P);
}

SEXP R_wazuni(SEXP R_svar, SEXP R_Y, SEXP R_eta, SEXP R_lams) {
    int i,I = nrows(R_Y);
    int j,J = length(R_svar);
    int k1,k2,K = ncols(R_lams), L1 = nrows(R_lams);
    int l,L2 = ncols(R_Y);
    double *dtmp1,dtmp2,dtmp3;

    SEXP R_P, R_dimnames;

    dtmp1 = Calloc(K,double);

    PROTECT(R_P = allocMatrix(REALSXP, I, J));
    PROTECT(R_dimnames = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(R_dimnames, 0, VECTOR_ELT(getAttrib(R_Y, R_DimNamesSymbol),0));
    SET_VECTOR_ELT(R_dimnames, 1, getAttrib(R_svar, R_NamesSymbol));
    setAttrib(R_P, R_DimNamesSymbol, R_dimnames);

    for(i=0;i<I;i++) {
        dtmp3 = 0.0F;
        for(j=0;j<J;j++) {
            REAL(R_P)[j*I+i] = 0.0F;
            for(k1=0;k1<K;k1++) {
                dtmp1[k1] = REAL(R_eta)[k1] + REAL(R_lams)[k1*L1 + INTEGER(R_svar)[j]];
                for(l=0;l<L2;l++) if(INTEGER(R_Y)[l*I + i] != NA_INTEGER)
                        dtmp1[k1] += REAL(R_lams)[k1*L1 + INTEGER(R_Y)[l*I + i]];
            }
            for(k1=0;k1<K;k1++) {
                dtmp2 = 1.0F;
                for(k2=0;k2<K;k2++) if(k2 != k1) dtmp2 *= (dtmp1[k2] - dtmp1[k1]);
                REAL(R_P)[j*I+i] += exp(-dtmp1[k1])/dtmp2;
            }
            dtmp3 += REAL(R_P)[j*I+i];
        }
        for(j=0;j<J;j++) REAL(R_P)[j*I+i] /= dtmp3;
    }
    Free(dtmp1);
    UNPROTECT(2);
    return(R_P);
}

/* The following functions are for the Cutler archtypal model

*/
SEXP R_cutler_L(SEXP R_X, SEXP R_lams, SEXP R_giks) {
    int i,i2,I = nrows(R_X);
    int j,J    = ncols(R_X);
    int k,K    = ncols(R_giks);
    double dtmp1,dtmp2;

    SEXP R_L;

    PROTECT(R_L = allocVector(REALSXP,1));

//    Rprintf("I %d J %d K %d\n",I,J,K);

    *REAL(R_L) = 0.0F;
    for(i=0;i<I;i++) {
        for(j=0;j<J;j++) {
            dtmp1 = 0.0F;
            for(k=0;k<K;k++) {
                dtmp2 = 0.0F;
                for(i2=0;i2<I;i2++) {
                    dtmp2 += REAL(R_lams)[k*I+i2] * REAL(R_X)[j*I+i2];
                }
                dtmp1 += REAL(R_giks)[k*I+i] * dtmp1;
            }
            dtmp1 = REAL(R_X)[j*I+i] - dtmp1;
            *REAL(R_L) += dtmp1 * dtmp1;
        }
    }
    UNPROTECT(1);
    return(R_L);
}
SEXP R_cutler_Lc(SEXP R_x, SEXP R_X, SEXP R_w, SEXP R_K) {
    int i,i2,I = nrows(R_X);
    int j,J    = ncols(R_X);
    int k,K    = *INTEGER(R_K);
    double dtmp1,dtmp2,*lams,*giks,*alphas,*betas;

    SEXP R_Lc, R_rms;

    PROTECT(R_Lc  = allocVector(REALSXP,1));
    PROTECT(R_rms = allocVector(REALSXP,1));

    lams   = REAL(R_x);
    giks   = lams   + I*K;
    alphas = giks   + I*K;
    betas  = alphas + I;

    *REAL(R_Lc)  = 0.0F;
    *REAL(R_rms) = 0.0F;
    for(i=0;i<I;i++) {
        for(j=0;j<J;j++) {
            dtmp1 = 0.0F;
            for(k=0;k<K;k++) {
                dtmp2 = 0.0F;
                for(i2=0;i2<I;i2++) {
                    dtmp2 += lams[k*I+i2] * REAL(R_X)[j*I+i2];
                }
                dtmp1 += giks[k*I+i] * dtmp1;
            }
            dtmp1 = REAL(R_X)[j*I+i] - dtmp1;
            *REAL(R_Lc) += dtmp1 * dtmp1;
        }
        dtmp1 = 0.0F;
        for(k=0;k<K;k++) {
            dtmp1 += giks[k*I+i];
        }
        dtmp1 -= 1.0F;
        *REAL(R_Lc) += alphas[i] * dtmp1;
        *REAL(R_Lc) += *REAL(R_w) / 2.0F * dtmp1 * dtmp1;
        *REAL(R_rms) += dtmp1 * dtmp1;
    }
    for(k=0;k<K;k++) {
        dtmp1 = 0.0F;
        for(i=0;i<I;i++) {
            dtmp1 += lams[k*I+i];
        }
        dtmp1 -= 1.0F;
        *REAL(R_Lc) += betas[k] * dtmp1;
        *REAL(R_Lc) += *REAL(R_w) / 2.0F * dtmp1 * dtmp1;
        *REAL(R_rms) += dtmp1 * dtmp1;
    }
    *REAL(R_rms) = sqrt(*REAL(R_rms)/(I+K));
    setAttrib(R_Lc, install("rms"), R_rms);
    UNPROTECT(2);
    return(R_Lc);
}
/* See page 66 of my technical log */
SEXP R_cutler_dLc(SEXP R_x, SEXP R_X, SEXP R_w, SEXP R_K) {
    int i0,i,I = nrows(R_X),i2;
    int j0,j,J = ncols(R_X);
    int k0,k,K = *INTEGER(R_K);
    double dtmp1,dtmp2,dtmp3,dtmp4;
    double *lams,*giks,*alphas,*betas;

    SEXP R_dgiks, R_dlams, R_dalphas, R_dbetas, R_list, R_names;

    PROTECT(R_dgiks   = allocMatrix(REALSXP,I,K));
    PROTECT(R_dlams   = allocMatrix(REALSXP,I,K));
    PROTECT(R_dalphas = allocMatrix(REALSXP,I,1));
    PROTECT(R_dbetas  = allocMatrix(REALSXP,K,1));

    lams   = REAL(R_x);
    giks   = lams   + I*K;
    alphas = giks   + I*K;
    betas  = alphas + I;

    Rprintf("R_cutler_dLc Starting I %d J %d K %d\n",I,J,K);

    for(k0=0;k0<K;k0++) {
        REAL(R_dbetas)[k0] = 0.0F;
        for(i=0;i<I;i++) {
            REAL(R_dbetas)[k0] += lams[k0*I+i];
        }
        REAL(R_dbetas)[k0] -= 1.0F;
    }
    for(i0=0;i0<I;i0++) {
        REAL(R_dalphas)[i0] = 0.0F;
        for(k0=0;k0<K;k0++) {
            REAL(R_dalphas)[i0] += giks[k0*I+i0];
            REAL(R_dgiks)[k0*I+i0] = 0.0F;
            REAL(R_dlams)[k0*I+i0] = 0.0F;
            for(j=0;j<J;j++) {
                dtmp1 = 0.0F;
                for(k=0;k<K;k++) {
                    dtmp2 = 0.0F;
                    for(i=0;i<I;i++) {
                        dtmp2 += lams[k*I+i] * REAL(R_X)[j*I+i];
                    }
                    dtmp1 += giks[k*I+i0] * dtmp2;
                } 
                dtmp2 = 0.0F;
                for(i=0;i<I;i++) {
                    dtmp2 += lams[k*I+i] * REAL(R_X)[j*I+i];

                    dtmp3 = 0.0F;
                    for(k=0;k<K;k++) {
                        dtmp4 = 0.0F;
                        for(i2=0;i2<I;i2++) {
                            dtmp4 += REAL(R_X)[j*I+i2] * lams[k*I+i2];
                        }
                        dtmp3 += giks[k*I+i] * dtmp4;
                    }
/*dlams*/           REAL(R_dlams)[k0*I+i0] += 2.0 * (REAL(R_X)[j*I+i] - dtmp3) * 
                        -giks[k0*I+i] * REAL(R_X)[j*I+i0];

                }
/*dgiks*/       REAL(R_dgiks)[k0*I+i0] += 2.0F * (REAL(R_X)[j*I+i0] - dtmp1) * (-dtmp2);
            }

            dtmp1 = 0.0F;
            for(i=0;i<I;i++) {
                dtmp1 += lams[k0*I+i];
            }
            REAL(R_dlams)[k0*I+i0] += betas[k0] + *REAL(R_w) * (dtmp1 - 1.0F);

            dtmp1 = 0.0F;
            for(k=0;k<K;k++) {
                dtmp1 += giks[k*I+i0];
            }
            REAL(R_dgiks)[k0*I+i0] += alphas[i0] + *REAL(R_w) * (dtmp1 - 1.0F);
        }
        REAL(R_dalphas)[i0] -= 1.0F;
    }

    PROTECT(R_list = allocVector(VECSXP,4));
    SET_VECTOR_ELT(R_list,0,R_dlams);
    SET_VECTOR_ELT(R_list,1,R_dgiks);
    SET_VECTOR_ELT(R_list,2,R_dalphas);
    SET_VECTOR_ELT(R_list,3,R_dbetas);

    PROTECT(R_names = allocVector(STRSXP,4));
    SET_STRING_ELT(R_names,0,mkChar("dlams"));
    SET_STRING_ELT(R_names,1,mkChar("dgiks"));
    SET_STRING_ELT(R_names,2,mkChar("dalphas"));
    SET_STRING_ELT(R_names,3,mkChar("dbetas"));

    setAttrib(R_list, R_NamesSymbol, R_names);
    
    UNPROTECT(6);
    return(R_list);
}

SEXP R_read_ts(SEXP R_file) {
    int i,I=100480507;
    char *file = (char *)CHAR(STRING_ELT(R_file,0));
    FILE *fp;
    int mid,rid,score;

    SEXP R_m, R_names, R_dimnames, R_mid, R_rid, R_score;

    PROTECT(R_m = allocMatrix(INTSXP,I,3));

    fp = fopen(file, "r");
    if(fp == NULL) error("File %s didn't open %d\n",file,errno);
    i=0;
    while(i < I && fscanf(fp,"%d,%d,%d\n",&mid,&rid,&score)) {
        INTEGER(R_m)[0*I+i] = mid;
        INTEGER(R_m)[1*I+i] = rid;
        INTEGER(R_m)[2*I+i] = score;
        i++;
    }
    fclose(fp);

    PROTECT(R_names = allocVector(STRSXP,3));
    SET_STRING_ELT(R_names, 0, mkChar("mid"));
    SET_STRING_ELT(R_names, 1, mkChar("rid"));
    SET_STRING_ELT(R_names, 2, mkChar("score"));

    PROTECT(R_dimnames = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(R_dimnames, 0, R_NilValue);
    SET_VECTOR_ELT(R_dimnames, 1, R_names);
    setAttrib(R_m, R_DimNamesSymbol, R_dimnames);

    UNPROTECT(3);
    return(R_m);
}

SEXP R_read_tsmsort(SEXP R_file) {
    int i,j,I=17770,J=233000;
    FILE *fp;
    int *vbuff,*sbuff;
    char *file = (char *)CHAR(STRING_ELT(R_file,0));
    int rv,mid,rid,score,last_mid;

    SEXP R_list, R_names, R_l1, R_s, R_v;

    Rprintf("file %s I %d J %d\n",file,I,J);
    PROTECT(R_list  = allocVector(VECSXP,I));
    PROTECT(R_names = allocVector(INTSXP,I));
    vbuff = Calloc(J,int);
    sbuff = Calloc(J,int);

    fp = fopen(file, "r");
    if(fp == NULL) error("File %s didn't open %d\n",file,errno);

    last_mid = -1;
    i=j=0;
    while(i < I && (rv = fscanf(fp,"%d,%d,%d\n",&mid,&rid,&score))) {
        if(mid != last_mid || rv == EOF) {
            if(j) {
                PROTECT(R_l1 = allocVector(VECSXP,2));
                PROTECT(R_s  = allocVector(INTSXP,j));
                memcpy(INTEGER(R_s),sbuff,j * sizeof(int));
                SET_VECTOR_ELT(R_l1,0,R_s);

                PROTECT(R_v = allocVector(INTSXP,j));
                memcpy(INTEGER(R_v),vbuff,j * sizeof(int));
                SET_VECTOR_ELT(R_l1,1,R_v);
                
                SET_VECTOR_ELT(R_list,i,R_l1);
                INTEGER(R_names)[i] = last_mid;
                UNPROTECT(3);
                i++;
            }
            j = 0;
        }
        vbuff[j] = rid;
        sbuff[j] = score;
        j++;
        if(j >= J) error("Buffer overrun!\n");
        last_mid = mid;
    }
    fclose(fp);

    setAttrib(R_list, R_NamesSymbol, R_names);
    Free(vbuff);
    Free(sbuff);
    UNPROTECT(2);
    return(R_list);
}
SEXP R_read_tsrsort(SEXP R_file) {
    int i,j,I=480189,J=18000;
    FILE *fp;
    int *vbuff,*sbuff;
    char *file = (char *)CHAR(STRING_ELT(R_file,0));
    int rv,mid,rid,score,last_rid;

    SEXP R_list, R_names, R_l1, R_s, R_v;

    Rprintf("file %s I %d J %d\n",file,I,J);
    PROTECT(R_list  = allocVector(VECSXP,I));
    PROTECT(R_names = allocVector(INTSXP,I));
    vbuff = Calloc(J,int);
    sbuff = Calloc(J,int);

    fp = fopen(file, "r");
    if(fp == NULL) error("File %s didn't open %d\n",file,errno);

    last_rid = -1;
    i=j=0;
    while(i < I && (rv = fscanf(fp,"%d,%d,%d\n",&mid,&rid,&score))) {
        if(rid != last_rid || rv == EOF) {
            if(j) {
                PROTECT(R_l1 = allocVector(VECSXP,2));
                PROTECT(R_s  = allocVector(INTSXP,j));
                memcpy(INTEGER(R_s),sbuff,j * sizeof(int));
                SET_VECTOR_ELT(R_l1,0,R_s);

                PROTECT(R_v = allocVector(INTSXP,j));
                memcpy(INTEGER(R_v),vbuff,j * sizeof(int));
                SET_VECTOR_ELT(R_l1,1,R_v);
                
                SET_VECTOR_ELT(R_list,i,R_l1);
                INTEGER(R_names)[i] = last_rid;
                UNPROTECT(3);
                i++;
            }
            j = 0;
        }
        vbuff[j] = mid;
        sbuff[j] = score;
        j++;
        if(j >= J) error("Buffer overrun!\n");
        last_rid = rid;
    }
    fclose(fp);

    setAttrib(R_list, R_NamesSymbol, R_names);
    Free(vbuff);
    Free(sbuff);
    UNPROTECT(2);
    return(R_list);
}

#define STAT_COUNT 7
#define OVERALL_MEAN 3.60429

int coipf_stats(SEXP R_list, int idx1, int idx2, int *ibuff, double *dbuff, double *giks, int gI, int K, int count, int skip, double *stats) {
    int j,J,k;
    int *v1,*v2;
    double dtmp,totd;
    double stat1,stat2,dist1,dist2,dmin,dmax,wmean;

    J  =  length(VECTOR_ELT(VECTOR_ELT(R_list,idx1),0));
    v1 = INTEGER(VECTOR_ELT(VECTOR_ELT(R_list,idx1),0));
    memcpy(ibuff,v1,J * sizeof(int));
    v2 = INTEGER(VECTOR_ELT(VECTOR_ELT(R_list,idx1),1)); // These must already be "indexed"

//    Rprintf("v2 %d %d %d J %d\n",v2[0],v2[1],v2[2],J);
//    Rprintf("idx2 giks %f %f %f\n",giks[0*gI+idx2],    giks[1*gI+idx2],    giks[2*gI+idx2]);
//    Rprintf("v2  giks %f %f %f\n",giks[0*gI+v2[0]-1],giks[1*gI+v2[0]-1],giks[2*gI+v2[0]-1]);

    for(j=0;j<J;j++) {
        dbuff[j] = 0.0F;
        for(k=0;k<K;k++) {
            dtmp = giks[k*gI+idx2] - giks[k*gI+v2[j]-1]; // idexes are 1 offset
            dbuff[j] += dtmp*dtmp;
        }
        dbuff[j] = sqrt(dbuff[j]) + 1e-3F;  // + 1e-3F to keep it real
    }
    rsort_with_index(dbuff,ibuff,J);

    stat1=stat2=dist1=dist2=0.0F;
    dmin=dmax=dbuff[skip];
    J = (J < count+skip) ? J : count+skip;
    for(j=skip;j<J;j++) {
        stat1 += (double)ibuff[j];
        stat2 += (double)ibuff[j]*ibuff[j];
        dist1 += dbuff[j];
        dist2 += dbuff[j]*dbuff[j];
        if(dmin > dbuff[j]) dmin = dbuff[j]; // min
        if(dmax < dbuff[j]) dmax = dbuff[j]; // max
    }
    // Weighted Mean Calculations
    wmean = 0.0F;
    totd = (J-skip)*(dmin+dmax)-dist1;
    for(j=skip;j<J;j++) { 
        wmean += (double)ibuff[j] * (dmin+dmax-dbuff[j])/totd;
    }

//  For model building....
//    for(j=0;j<count;j++) {
//        if(j<J) Rprintf("%d ",ibuff[j]);
//        else    Rprintf("NA ");
//    }
//    for(j=0;j<count;j++) {
//        if(j<J) Rprintf("%f ",(dmin+dmax-dbuff[j])/totd);
//        else    Rprintf("NA ");
//    }
//    Rprintf("\n");
//        
/*
    if(1.0F > wmean ||  wmean > 5.0F || isnan(wmean)) {
        Rprintf("wmean %f j %d dmin %f dmax %f dist1 %f totd %f\n",wmean,j,dmin,dmax,dist1,totd);
        for(k=skip;k<J;k++)
            Rprintf("  %d %d %f %f\n",k,ibuff[k],(dmin+dmax-dbuff[k])/totd,dbuff[k]);
    }
*/
    J -= skip;
    if(J) {
        stats[0] = wmean;
        stats[1] = stat1/J;
        stats[2] = sqrt((stat2 - stat1*stat1/J)/J);
        stats[3] = dist1/J;
        stats[4] = sqrt((dist2 - dist1*dist1/J)/J);
        stats[5] = dmin;
        stats[6] = dmax;
    } else {
        stats[0] = stats[1] = stats[5] = stats[6] = OVERALL_MEAN;
        stats[2] = 0.0F;
        stats[3] = 0.0F;
        stats[4] = 0.0F;
    }

    return(J);
}

SEXP R_predict_coipf(SEXP R_X, SEXP R_model, SEXP R_mnum, SEXP R_rnum, SEXP R_skip) {
    int i, j, k, I = nrows(R_X), J, mJ,rJ;
    int ridx,midx;
    int mid,rid,score, mnum = *INTEGER(R_mnum), rnum = *INTEGER(R_rnum);
    int *ibuff;
    double *dbuff,*stats;

    SEXP R_list,R_names,R_pscore;
    SEXP R_mwmean,R_mstat1,R_mstat2,R_mdist1,R_mdist2,R_mdmin,R_mdmax;
    SEXP R_rwmean,R_rstat1,R_rstat2,R_rdist1,R_rdist2,R_rdmin,R_rdmax;
    SEXP R_ml, R_rl, R_mgiks, R_rgiks, R_rids;

    PROTECT(R_pscore = allocVector(REALSXP,I));
    PROTECT(R_mwmean = allocVector(REALSXP,I));
    PROTECT(R_mstat1 = allocVector(REALSXP,I));
    PROTECT(R_mstat2 = allocVector(REALSXP,I));
    PROTECT(R_mdist1 = allocVector(REALSXP,I));
    PROTECT(R_mdist2 = allocVector(REALSXP,I));
    PROTECT(R_mdmin  = allocVector(REALSXP,I));
    PROTECT(R_mdmax  = allocVector(REALSXP,I));
    PROTECT(R_rwmean = allocVector(REALSXP,I));
    PROTECT(R_rstat1 = allocVector(REALSXP,I));
    PROTECT(R_rstat2 = allocVector(REALSXP,I));
    PROTECT(R_rdist1 = allocVector(REALSXP,I));
    PROTECT(R_rdist2 = allocVector(REALSXP,I));
    PROTECT(R_rdmin  = allocVector(REALSXP,I));
    PROTECT(R_rdmax  = allocVector(REALSXP,I));

    PROTECT(R_ml    = getListElement(R_model,"ml"));
    PROTECT(R_rl    = getListElement(R_model,"rl"));
    PROTECT(R_mgiks = coerceVector(getListElement(R_model,"mgiks"),REALSXP));
    PROTECT(R_rgiks = coerceVector(getListElement(R_model,"rgiks"),REALSXP));
    PROTECT(R_rids  = coerceVector(duplicate(getAttrib(R_rl,R_NamesSymbol)),INTSXP));
    rJ = length(R_rl);
    mJ = length(R_ml);

    J = (mJ<rJ) ? rJ : mJ;
    dbuff = Calloc(J, double);
    ibuff = Calloc(J, int);
    stats = Calloc(STAT_COUNT, double);

//    Rprintf("I %d J %d rJ %d mJ %d num %d\n",I,J,rJ,mJ,num);
    for(i=0;i<I;i++) {
        mid   = INTEGER(R_X)[0*I+i];
        rid   = INTEGER(R_X)[1*I+i];
//        Rprintf("%d \n",INTEGER(R_X)[2*I+i]);
//        score = INTEGER(R_X)[2*I+i];

        midx = mid-1; // For zero based indexing!
        ridx = (int*)bsearch(&rid,INTEGER(R_rids),rJ,sizeof(int),q_intcmp) - INTEGER(R_rids);

//        Rprintf("i %.4d mid %d midx %d rid %d ridx %d score %d\n",i,mid,midx,rid,ridx,score);

        coipf_stats(R_ml,midx,ridx,ibuff,dbuff,REAL(R_rgiks),nrows(R_rgiks),ncols(R_rgiks),mnum,*INTEGER(R_skip),stats);
        REAL(R_mwmean)[i] = stats[0];
        REAL(R_mstat1)[i] = stats[1];
        REAL(R_mstat2)[i] = stats[2];
        REAL(R_mdist1)[i] = stats[3];
        REAL(R_mdist2)[i] = stats[4];
        REAL(R_mdmin)[i]  = stats[5];
        REAL(R_mdmax)[i]  = stats[6];

        coipf_stats(R_rl,ridx,midx,ibuff,dbuff,REAL(R_mgiks),nrows(R_mgiks),ncols(R_mgiks),rnum,*INTEGER(R_skip),stats);
        REAL(R_rwmean)[i] = stats[0];
        REAL(R_rstat1)[i] = stats[1];
        REAL(R_rstat2)[i] = stats[2];
        REAL(R_rdist1)[i] = stats[3];
        REAL(R_rdist2)[i] = stats[4];
        REAL(R_rdmin)[i]  = stats[5];
        REAL(R_rdmax)[i]  = stats[6];

        REAL(R_pscore)[i] = (REAL(R_mwmean)[i] + REAL(R_rwmean)[i])/2.0F;
        R_CheckUserInterrupt();
    }

    PROTECT(R_list = allocVector(VECSXP,15));
    SET_VECTOR_ELT(R_list,0,R_pscore);
    SET_VECTOR_ELT(R_list,1,R_mwmean);
    SET_VECTOR_ELT(R_list,2,R_mstat1);
    SET_VECTOR_ELT(R_list,3,R_mstat2);
    SET_VECTOR_ELT(R_list,4,R_mdist1);
    SET_VECTOR_ELT(R_list,5,R_mdist2);
    SET_VECTOR_ELT(R_list,6,R_mdmin);
    SET_VECTOR_ELT(R_list,7,R_mdmax);
    SET_VECTOR_ELT(R_list,8,R_rwmean);
    SET_VECTOR_ELT(R_list,9,R_rstat1);
    SET_VECTOR_ELT(R_list,10,R_rstat2);
    SET_VECTOR_ELT(R_list,11,R_rdist1);
    SET_VECTOR_ELT(R_list,12,R_rdist2);
    SET_VECTOR_ELT(R_list,13,R_rdmin);
    SET_VECTOR_ELT(R_list,14,R_rdmax);

    PROTECT(R_names = allocVector(STRSXP,15));
    SET_STRING_ELT(R_names,0,mkChar("pscore"));
    SET_STRING_ELT(R_names,1,mkChar("mwmean"));
    SET_STRING_ELT(R_names,2,mkChar("msmean"));
    SET_STRING_ELT(R_names,3,mkChar("msstd"));
    SET_STRING_ELT(R_names,4,mkChar("mdmean"));
    SET_STRING_ELT(R_names,5,mkChar("mdstd"));
    SET_STRING_ELT(R_names,6,mkChar("mdmin"));
    SET_STRING_ELT(R_names,7,mkChar("mdmax"));
    SET_STRING_ELT(R_names,8,mkChar("rwmean"));
    SET_STRING_ELT(R_names,9,mkChar("rsmean"));
    SET_STRING_ELT(R_names,10,mkChar("rsstd"));
    SET_STRING_ELT(R_names,11,mkChar("rdmean"));
    SET_STRING_ELT(R_names,12,mkChar("rdstd"));
    SET_STRING_ELT(R_names,13,mkChar("rdmin"));
    SET_STRING_ELT(R_names,14,mkChar("rdmax"));

    setAttrib(R_list, R_NamesSymbol, R_names);

    Free(dbuff);
    Free(stats);
    Free(ibuff);
    UNPROTECT(22);
    return(R_list);
}

typedef struct {
    int *X,mnum,rnum,skip,I;
    SEXP R_pscore;
    SEXP R_mwmean,R_mstat1,R_mstat2,R_mdist1,R_mdist2,R_mdmin,R_mdmax;
    SEXP R_rwmean,R_rstat1,R_rstat2,R_rdist1,R_rdist2,R_rdmin,R_rdmax;
    SEXP R_ml, R_rl, R_mgiks, R_rgiks, R_rids;
} coipf_struct;
typedef struct {
    coipf_struct *s;
    int i1,i2;
} pcoipf_struct;

void *pcoipf(void *arg) {
    int i,I,J,rJ,mJ,mid,rid,midx,ridx,*ibuff;
    double *dbuff,*stats;
    pcoipf_struct *a = arg;
    coipf_struct  *s = a->s;

    rJ = length(s->R_rl);
    mJ = length(s->R_ml);

    J = (mJ<rJ) ? rJ : mJ;
    dbuff = Calloc(J, double);
    ibuff = Calloc(J, int);
    stats = Calloc(STAT_COUNT, double);

    for(i=a->i1;i<a->i2;i++) {
        mid   = s->X[0*s->I+i];
        rid   = s->X[1*s->I+i];

        midx = mid-1; // For zero based indexing!
        ridx = (int*)bsearch(&rid,INTEGER(s->R_rids),rJ,sizeof(int),q_intcmp) - INTEGER(s->R_rids);

        coipf_stats(s->R_ml,midx,ridx,ibuff,dbuff,REAL(s->R_rgiks),nrows(s->R_rgiks),ncols(s->R_rgiks),s->mnum,s->skip,stats);
        REAL(s->R_mwmean)[i] = stats[0];
        REAL(s->R_mstat1)[i] = stats[1];
        REAL(s->R_mstat2)[i] = stats[2];
        REAL(s->R_mdist1)[i] = stats[3];
        REAL(s->R_mdist2)[i] = stats[4];
        REAL(s->R_mdmin)[i]  = stats[5];
        REAL(s->R_mdmax)[i]  = stats[6];

        coipf_stats(s->R_rl,ridx,midx,ibuff,dbuff,REAL(s->R_mgiks),nrows(s->R_mgiks),ncols(s->R_mgiks),s->rnum,s->skip,stats);
        REAL(s->R_rwmean)[i] = stats[0];
        REAL(s->R_rstat1)[i] = stats[1];
        REAL(s->R_rstat2)[i] = stats[2];
        REAL(s->R_rdist1)[i] = stats[3];
        REAL(s->R_rdist2)[i] = stats[4];
        REAL(s->R_rdmin)[i]  = stats[5];
        REAL(s->R_rdmax)[i]  = stats[6];

        REAL(s->R_pscore)[i] = (REAL(s->R_mwmean)[i] + REAL(s->R_rwmean)[i])/2.0F;
//        R_CheckUserInterrupt();
        sched_yield();
    }
    Free(dbuff);
    Free(ibuff);
    Free(stats);

    pthread_exit((void *) 0);
}

/* pthreaded version! */
SEXP R_predict_pcoipf(SEXP R_X, SEXP R_model, SEXP R_mnum, SEXP R_rnum, SEXP R_skip, SEXP R_threads) {
    int t,T;
    int status,len,idx=0;

    coipf_struct s;
    pcoipf_struct *a;
    pthread_attr_t attr;
    pthread_t *threads;

    SEXP R_list,R_names;

    s.X    = INTEGER(R_X);
    s.mnum  = *INTEGER(R_mnum);
    s.rnum  = *INTEGER(R_rnum);
    s.skip = *INTEGER(R_skip);
    s.I    = nrows(R_X);
    T = (*INTEGER(R_threads)<s.I) ? *INTEGER(R_threads) : s.I;

    PROTECT(s.R_pscore = allocVector(REALSXP,s.I));
    PROTECT(s.R_mwmean = allocVector(REALSXP,s.I));
    PROTECT(s.R_mstat1 = allocVector(REALSXP,s.I));
    PROTECT(s.R_mstat2 = allocVector(REALSXP,s.I));
    PROTECT(s.R_mdist1 = allocVector(REALSXP,s.I));
    PROTECT(s.R_mdist2 = allocVector(REALSXP,s.I));
    PROTECT(s.R_mdmin  = allocVector(REALSXP,s.I));
    PROTECT(s.R_mdmax  = allocVector(REALSXP,s.I));
    PROTECT(s.R_rwmean = allocVector(REALSXP,s.I));
    PROTECT(s.R_rstat1 = allocVector(REALSXP,s.I));
    PROTECT(s.R_rstat2 = allocVector(REALSXP,s.I));
    PROTECT(s.R_rdist1 = allocVector(REALSXP,s.I));
    PROTECT(s.R_rdist2 = allocVector(REALSXP,s.I));
    PROTECT(s.R_rdmin  = allocVector(REALSXP,s.I));
    PROTECT(s.R_rdmax  = allocVector(REALSXP,s.I));

    PROTECT(s.R_ml    = getListElement(R_model,"ml"));
    PROTECT(s.R_rl    = getListElement(R_model,"rl"));
    PROTECT(s.R_mgiks = coerceVector(getListElement(R_model,"mgiks"),REALSXP));
    PROTECT(s.R_rgiks = coerceVector(getListElement(R_model,"rgiks"),REALSXP));
    PROTECT(s.R_rids  = coerceVector(duplicate(getAttrib(s.R_rl,R_NamesSymbol)),INTSXP));
    
    threads = Calloc(T,pthread_t);
    a       = Calloc(T,pcoipf_struct);

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    pthread_attr_setschedpolicy(&attr, SCHED_RR);
//    pthread_attr_setscope(&attr, PTHREAD_SCOPE_PROCESS);
//    pthread_setconcurrency(T);

    len = (int) ceil((double)s.I/(double)T);
    for(t=0;t<T;t++) {
        a[t].s  = &s;
        a[t].i1 = idx; idx += len;
        a[t].i2 = (idx<s.I) ? idx : s.I;
        pthread_create(&threads[t], &attr, pcoipf, (void *) &a[t]);
    }
    for(t=0;t<T;t++) {
        pthread_join(threads[t], (void **) &status);
    }

    PROTECT(R_list = allocVector(VECSXP,15));
    SET_VECTOR_ELT(R_list,0,s.R_pscore);
    SET_VECTOR_ELT(R_list,1,s.R_mwmean);
    SET_VECTOR_ELT(R_list,2,s.R_mstat1);
    SET_VECTOR_ELT(R_list,3,s.R_mstat2);
    SET_VECTOR_ELT(R_list,4,s.R_mdist1);
    SET_VECTOR_ELT(R_list,5,s.R_mdist2);
    SET_VECTOR_ELT(R_list,6,s.R_mdmin);
    SET_VECTOR_ELT(R_list,7,s.R_mdmax);
    SET_VECTOR_ELT(R_list,8,s.R_rwmean);
    SET_VECTOR_ELT(R_list,9,s.R_rstat1);
    SET_VECTOR_ELT(R_list,10,s.R_rstat2);
    SET_VECTOR_ELT(R_list,11,s.R_rdist1);
    SET_VECTOR_ELT(R_list,12,s.R_rdist2);
    SET_VECTOR_ELT(R_list,13,s.R_rdmin);
    SET_VECTOR_ELT(R_list,14,s.R_rdmax);

    PROTECT(R_names = allocVector(STRSXP,15));
    SET_STRING_ELT(R_names,0,mkChar("pscore"));
    SET_STRING_ELT(R_names,1,mkChar("mwmean"));
    SET_STRING_ELT(R_names,2,mkChar("msmean"));
    SET_STRING_ELT(R_names,3,mkChar("msstd"));
    SET_STRING_ELT(R_names,4,mkChar("mdmean"));
    SET_STRING_ELT(R_names,5,mkChar("mdstd"));
    SET_STRING_ELT(R_names,6,mkChar("mdmin"));
    SET_STRING_ELT(R_names,7,mkChar("mdmax"));
    SET_STRING_ELT(R_names,8,mkChar("rwmean"));
    SET_STRING_ELT(R_names,9,mkChar("rsmean"));
    SET_STRING_ELT(R_names,10,mkChar("rsstd"));
    SET_STRING_ELT(R_names,11,mkChar("rdmean"));
    SET_STRING_ELT(R_names,12,mkChar("rdstd"));
    SET_STRING_ELT(R_names,13,mkChar("rdmin"));
    SET_STRING_ELT(R_names,14,mkChar("rdmax"));

    setAttrib(R_list, R_NamesSymbol, R_names);

    Free(threads);
    Free(a);
    UNPROTECT(22);
    return(R_list);
}

SEXP export_netflix(SEXP R_list, SEXP R_I, SEXP R_file) {
    int i, I = *INTEGER(R_I);
    char *file = (char *)CHAR(STRING_ELT(R_file,0));
    FILE *fp;
    int *mid, *rid, last_mid;
    double *pscore;

    SEXP R_mid, R_rid, R_pscore;

    fp = fopen(file, "w");
    if(fp == NULL) error("File %s didn't open %d\n",file,errno);

    R_mid = coerceVector(getListElement(R_list,"mid"),INTSXP);
    if(R_mid == R_NilValue) error("Can't find mid column\n");
    mid = INTEGER(R_mid);
    R_rid = coerceVector(getListElement(R_list,"rid"),INTSXP);
    if(R_rid == R_NilValue) error("Can't find rid column\n");
    rid = INTEGER(R_rid);
    R_pscore = coerceVector(getListElement(R_list,"pscore"),REALSXP);
    if(R_pscore == R_NilValue) error("Can't find pscore column\n");
    pscore = REAL(R_pscore);

    Rprintf("Writing %d records to file %s ... ",I,file);

    last_mid=-1;
    for(i=0;i<I;i++) {
        if(last_mid != mid[i]) fprintf(fp,"%d:\n",mid[i]);
        fprintf(fp,"%f\n",pscore[i]);
        last_mid = mid[i];
    }
    Rprintf("done.\n");
    fclose(fp);

    return(R_NilValue);
}
