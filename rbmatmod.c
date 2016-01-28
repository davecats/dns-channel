/*-O3 */
#define _GNU_SOURCE

#define _FILE_OFFSET_BITS 64

#define _LARGE_FILES


#include <unistd.h>

#include <stdlib.h>

#include <stdio.h>
/*** typedef __builtin_va_list __gnuc_va_list; ununderstood ***/
/*** typedef __gnuc_va_list va_list; ununderstood ***/

#include <fcntl.h>

#include <math.h>

#include <limits.h>

#include <float.h>

#include <string.h>


#include <sys/time.h>

#include <sys/types.h>

#include <setjmp.h>

#include <errno.h>

#include <signal.h>
struct arrptr{int l,h; size_t i; char *a;};
struct dynptr{int t; void* p;};
jmp_buf errorhook;
#if defined NO_fenv
  #define feenableexcept(fpe)
  #define feclearexcept(fpe)
#else
  #include <fenv.h>
  #define fpe FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO
#endif

/*!(
#if defined __linux__ && defined __i386__
  #include <fpu_control.h>
  fpu_control_t fpumask;
  #define fpudef fpumask=_FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM)
#else
  #define _FPU_SETCW(fpumask)
  #define fpudef
#endif
#if defined __linux__ && defined __SSE__
  #include <xmmintrin.h>
  unsigned int ssemask;
  #define ssedef ssemask=_MM_MASK_INVALID | _MM_MASK_DIV_ZERO | _MM_MASK_OVERFLOW
#else
  #define _MM_SET_EXCEPTION_MASK(ssemask)
  #define _MM_SET_EXCEPTION_STATE(ssemask)
  #define fpudef
#endif
!)*/
char ERRORMESSAGE_[(80+1)];

#define mmovefrom(var,buf,type) *(type *)(buf)=*var

#define mmoveto(var,buf,type) *var=*(type *)(buf)

#define printERR fprintf(stderr,"\r%.*s: PROGRAM HALTED  \n",(int)sizeof(ERRORMESSAGE_),ERRORMESSAGE_);fflush(stderr)

#define mainstart feenableexcept(fpe); {void (*sig)(int); if((sig=signal(SIGSEGV,trapsignal))!=SIG_DFL)signal(SIGSEGV,sig); if((sig=signal(SIGFPE,trapsignal))!=SIG_DFL)signal(SIGFPE,sig); if((sig=signal(SIGILL,trapsignal))!=SIG_DFL)signal(SIGILL,sig); if((sig=signal(SIGINT,trapsignal))!=SIG_DFL)signal(SIGINT,sig); else if(setjmp(errorhook)){   if(ERRORMESSAGE_[0]){     printERR;     return -1;   }else return 0; }}
int dynptrerr(int type);void *errmalloc();void ioerr(FILE *fil);void errfclose(FILE **f);void errfopen(FILE **f, const char *name, int mode);int scanrec(FILE *f, const char *format, void *var) ;int myfgets(char *name, char *var, char *varend, FILE *f) ;int mygetline(char *name, char **var, FILE *f) ;void trapsignal(int signum);

/*!*/
/*! This program computes the power spectral density */
/*! and cross spectral density of:*/
/*!*/
/*! (uu, vv, ww, uv, vw, uw)*/
/*!*/
/*! This program IS parallel*/
/*!*/

/*gamma=0*/
/*outinterv=0*/
struct COMPLEX_{double REAL_;double IMAG_;};void complex_1INV(struct COMPLEX_ *RESULT_,double x_REAL,double x_IMAG);
void complex_2EXP(struct COMPLEX_ *RESULT_,double x_REAL,double x_IMAG);
void complex_3SINH(struct COMPLEX_ *RESULT_,double x_REAL,double x_IMAG);
void complex_4COSH(struct COMPLEX_ *RESULT_,double x_REAL,double x_IMAG);
void complex_5TANH(struct COMPLEX_ *RESULT_,double x_REAL,double x_IMAG);
void complex_6COTH(struct COMPLEX_ *RESULT_,double x_REAL,double x_IMAG);
void complex_7LOG(struct COMPLEX_ *RESULT_,double x_REAL,double x_IMAG);
void complex_8POWER(struct COMPLEX_ *RESULT_,double x_REAL,double x_IMAG,double y_REAL,double y_IMAG);
void complex_9SQRT(struct COMPLEX_ *RESULT_,double x_REAL,double x_IMAG);
void complex_10CRAND(struct COMPLEX_ *RESULT_);
void complex_11CGAUSS(struct COMPLEX_ *RESULT_);
void fft_1IFT(int Rin_l,int Rin_h,size_t Rin_i,char *Rin__,int Rout_l,int Rout_h,size_t Rout_i,char *Rout__);
void fft_2FFT(int Rin_l,int Rin_h,size_t Rin_i,char *Rin__,int Rout_l,int Rout_h,size_t Rout_i,char *Rout__);
void fft_3RFT(int Rin_l,int Rin_h,size_t Rin_i,char *Rin__,int Rout_l,int Rout_h,size_t Rout_i,char *Rout__);
void fft_4HFT(int Rin_l,int Rin_h,size_t Rin_i,char *Rin__,int Rout_l,int Rout_h,size_t Rout_i,char *Rout__);
void fft_5IFTU(int Rin_l,int Rin_h,size_t Rin_i,char *Rin__,int Rout_l,int Rout_h,size_t Rout_i,char *Rout__);
void fft_6FFTU(int Rin_l,int Rin_h,size_t Rin_i,char *Rin__,int Rout_l,int Rout_h,size_t Rout_i,char *Rout__);
void fft_7RFTU(int Rin_l,int Rin_h,size_t Rin_i,char *Rin__,int Rout_l,int Rout_h,size_t Rout_i,char *Rout__);
void fft_8HFTU(int Rin_l,int Rin_h,size_t Rin_i,char *Rin__,int Rout_l,int Rout_h,size_t Rout_i,char *Rout__);
int fft_9FFTfit(int N_);
void fft_fft_10ReverseInc(int *K_,int N_);
int fft_fft_LASTN;
struct arrptr fft_fft_CEXP;
struct arrptr fft_fft_tempVEC;
struct arrptr fft_fft_RI;
char *fft_fft_RI2;
char *fft_fft_RI3;
void fft_fft_11SETUP(int N_);
double fft_fft_12C3;
/*C3=fft_fft_12C3*/
void fft_fft_13BTFLY(int N_,int M_);
void fft_fft_14BTFLYI(int N_);
void fft_1IFT(int Rin_l,int Rin_h,size_t Rin_i,char *Rin__,int Rout_l,int Rout_h,size_t Rout_i,char *Rout__);
void fft_2FFT(int Rin_l,int Rin_h,size_t Rin_i,char *Rin__,int Rout_l,int Rout_h,size_t Rout_i,char *Rout__);
void fft_3RFT(int Rin_l,int Rin_h,size_t Rin_i,char *Rin__,int Rout_l,int Rout_h,size_t Rout_i,char *Rout__);
void fft_4HFT(int Rin_l,int Rin_h,size_t Rin_i,char *Rin__,int Rout_l,int Rout_h,size_t Rout_i,char *Rout__);
void fft_5IFTU(int Rin_l,int Rin_h,size_t Rin_i,char *Rin__,int Rout_l,int Rout_h,size_t Rout_i,char *Rout__);
void fft_6FFTU(int Rin_l,int Rin_h,size_t Rin_i,char *Rin__,int Rout_l,int Rout_h,size_t Rout_i,char *Rout__);
void fft_7RFTU(int Rin_l,int Rin_h,size_t Rin_i,char *Rin__,int Rout_l,int Rout_h,size_t Rout_i,char *Rout__);
void fft_8HFTU(int Rin_l,int Rin_h,size_t Rin_i,char *Rin__,int Rout_l,int Rout_h,size_t Rout_i,char *Rout__);
double *fft_15REALIFIED(int x_l,int x_h,size_t x_i,char *x__,int y_);
double DotProduct(int a_l,int a_h,size_t a_i,char *a__,int b_l,int b_h,size_t b_i,char *b__);
void rbmat_1LeftMult(int c_l,int c_h,size_t c_i,char *c__,int A_l,int A_h,size_t A_i,int A__l,int A__h,size_t A__i,char *A___,int b_l,int b_h,size_t b_i,char *b__);
void RightMult(int c_l,int c_h,size_t c_i,char *c__,int a_l,int a_h,size_t a_i,char *a__,int B_l,int B_h,size_t B_i,int B__l,int B__h,size_t B__i,char *B___);
struct REALINVMAT_{int lo_;int hi_;};void rbmat_2MatEqu(int c_l,int c_h,size_t c_i,int c__l,int c__h,size_t c__i,char *c___,double a_);
void rbmat_3MatEqu(FILE *c_f,int c_l,int c_h,size_t c_i,int c__l,int c__h,size_t c__i,off_t c___,double a_);
void rbmat_4LeftMult(int c_l,int c_h,size_t c_i,char *c__,FILE *A_f,int A_l,int A_h,size_t A_i,int A__l,int A__h,size_t A__i,off_t A___,int b_l,int b_h,size_t b_i,char *b__);
void rbmat_5RightMult(int c_l,int c_h,size_t c_i,char *c__,int a_l,int a_h,size_t a_i,char *a__,FILE *B_f,int B_l,int B_h,size_t B_i,int B__l,int B__h,size_t B__i,off_t B___);
void rbmat_6LUdecomp(int AA_l,int AA_h,size_t AA_i,int AA__l,int AA__h,size_t AA__i,char *AA___);
void rbmat_7LUdecomp(FILE *A_f,int A_l,int A_h,size_t A_i,int A__l,int A__h,size_t A__i,off_t A___);
void rbmat_8LeftLDiv(int x_l,int x_h,size_t x_i,char *x__,int A_l,int A_h,size_t A_i,int A__l,int A__h,size_t A__i,char *A___,int t_l,int t_h,size_t t_i,char *t__);
void rbmat_9LeftUDiv(int x_l,int x_h,size_t x_i,char *x__,int A_l,int A_h,size_t A_i,int A__l,int A__h,size_t A__i,char *A___,int t_l,int t_h,size_t t_i,char *t__);
void rbmat_10LeftLUDiv(int x_l,int x_h,size_t x_i,char *x__,FILE *A_f,int A_l,int A_h,size_t A_i,int A__l,int A__h,size_t A__i,off_t A___,int t_l,int t_h,size_t t_i,char *t__);
void rbmat_11RightLUDiv(int x_l,int x_h,size_t x_i,char *x__,int t_l,int t_h,size_t t_i,char *t__,int A_l,int A_h,size_t A_i,int A__l,int A__h,size_t A__i,char *A___);
void rbmat_12RightLUDiv(int x_l,int x_h,size_t x_i,char *x__,int t_l,int t_h,size_t t_i,char *t__,FILE *A_f,int A_l,int A_h,size_t A_i,int A__l,int A__h,size_t A__i,off_t A___);
double Lanczos_norm_;
void rbmat_13Lanczos(int x_l,int x_h,size_t x_i,char *x__,void (*A_)(int,int,size_t,char *,int,int,size_t,char *),void (*AT_)(int,int,size_t,char *,int,int,size_t,char *),int y1_l,int y1_h,size_t y1_i,char *y1__,double eps_);
void rbmat_14PLU(int m_l,int m_h,size_t m_i,int m__l,int m__h,size_t m__i,char *m___,struct REALINVMAT_ *RESULT_);
void rbmat_15LeftLUDiv(int x_l,int x_h,size_t x_i,char *x__,struct REALINVMAT_ *m_,int t_l,int t_h,size_t t_i,char *t__);
void rbmat_16RightLUDiv(int x_l,int x_h,size_t x_i,char *x__,int t_l,int t_h,size_t t_i,char *t__,struct REALINVMAT_ *m_);
void rbmat_17INV(int mat_l,int mat_h,size_t mat_i,int mat__l,int mat__h,size_t mat__i,char *mat___,int RESULT_l,int RESULT_h,size_t RESULT_i,int RESULT__l,int RESULT__h,size_t RESULT__i,char *RESULT___);
double rbmat_18DET(int mat_l,int mat_h,size_t mat_i,int mat__l,int mat__h,size_t mat__i,char *mat___);
void rbmat_19Lanczos(int x_l,int x_h,size_t x_i,char *x__,int mat_l,int mat_h,size_t mat_i,int mat__l,int mat__h,size_t mat__i,char *mat___,int y_l,int y_h,size_t y_i,char *y__,double eps_);

struct rbmat_Lanczos_R_s20 {int l,h; size_t i;struct arrptr a;};
struct rbmat_Lanczos_R_s20 rbmat_Lanczos_R_Lanczos_mat;
void rbmat_Lanczos_R_21A(int v2_l,int v2_h,size_t v2_i,char *v2__,int v1_l,int v1_h,size_t v1_i,char *v1__);
void rbmat_Lanczos_R_22AT(int v2_l,int v2_h,size_t v2_i,char *v2__,int v1_l,int v1_h,size_t v1_i,char *v1__);
void rbmat_19Lanczos(int x_l,int x_h,size_t x_i,char *x__,int mat_l,int mat_h,size_t mat_i,int mat__l,int mat__h,size_t mat__i,char *mat___,int y_l,int y_h,size_t y_i,char *y__,double eps_);

#include <sys/mman.h>

#include <sys/wait.h>
/*** typedef struct sigevent ununderstood ***/
/*** {
 sigval_t sigev_value;
 int sigev_signo;
 int sigev_notify;

 union
 {
	int _pad[((64 / sizeof (int)) - 4)];

	
	__pid_t _tid;

	struct
	 {
	 void (*_function) (sigval_t);	
	 pthread_attr_t *_attribute;		
	 } _sigev_thread;
 } _sigev_un;
 } sigevent_t; ununderstood ***/

#include <sys/shm.h>


#define SHMPAGE 4194304
size_t shmavail;char *shmaddr;void *shmalloc(size_t size);sigset_t oldmask;void donothing(int signum);void setup_signal_USR1();


#include <sys/socket.h>

#include <netinet/in.h>
/*** typedef uint32_t in_addr_t; ununderstood ***/
/*** in_addr_t s_addr ununderstood ***/
/*** typedef uint16_t in_port_t; ununderstood ***/
/*** uint8_t	__u6_addr8[16] ununderstood ***/
/*** uint16_t __u6_addr16[8] ununderstood ***/
/*** uint32_t __u6_addr32[4] ununderstood ***/
/*** in_port_t sin_port ununderstood ***/
/*** in_port_t sin6_port ununderstood ***/
/*** uint32_t sin6_flowinfo ununderstood ***/
/*** uint32_t sin6_scope_id ununderstood ***/
/*** uint32_t gr_interface ununderstood ***/
/*** uint32_t gsr_interface ununderstood ***/
/*** uint32_t imsf_fmode ununderstood ***/
/*** uint32_t imsf_numsrc ununderstood ***/
/*** uint32_t gf_interface ununderstood ***/
/*** uint32_t gf_fmode ununderstood ***/
/*** uint32_t gf_numsrc ununderstood ***/
/*** extern uint32_t ntohl (uint32_t __netlong)  __attribute__ ((__const__)); ununderstood ***/
/*** extern uint16_t ntohs (uint16_t __netshort) ununderstood ***/
/*** __attribute__ ((__const__)); ununderstood ***/
/*** extern uint32_t htonl (uint32_t __hostlong) ununderstood ***/
/*** __attribute__ ((__const__)); ununderstood ***/
/*** extern uint16_t htons (uint16_t __hostshort) ununderstood ***/
/*** __attribute__ ((__const__)); ununderstood ***/
/*** uint32_t ip6m_mtu ununderstood ***/
/*** extern uint8_t *inet6_option_alloc (struct cmsghdr *__cmsg, int __datalen,
				 int __multx, int __plusy) ununderstood ***/

#include <netinet/tcp.h>
/*** # error "Adjust your <bits/endian.h> defines"
u_int16_t window ununderstood ***/
/*** enum ununderstood ***/
/*** {
 TCP_NO_QUEUE,
 TCP_RECV_QUEUE,
 TCP_SEND_QUEUE,
 TCP_QUEUES_NR,
}; ununderstood ***/

#include <netdb.h>
/*** typedef struct sigevent ununderstood ***/
/*** {
 sigval_t sigev_value;
 int sigev_signo;
 int sigev_notify;

 union
 {
	int _pad[((64 / sizeof (int)) - 4)];

	
	__pid_t _tid;

	struct
	 {
	 void (*_function) (sigval_t);	
	 pthread_attr_t *_attribute;		
	 } _sigev_thread;
 } _sigev_un;
 } sigevent_t; ununderstood ***/
/*** uint32_t n_net ununderstood ***/
int tcpserver(uint16_t port)
;int tcpclient(const char *hostname, uint16_t port) 
;int udpsocket(uint16_t myport, const char *hostname, uint16_t hostport) 
;

volatile int *barrier_;

/*!USE rtchecks*/
/*nsmp=4*/
int iproc_;
int nproc_;
/*bufsize=800*/
/*baseport=(IPPORT_USERRESERVED+111)*/
FILE *prev_;
FILE *next_;
/**/
void rbmatmod_1LUdecompStep(int A_l,int A_h,size_t A_i,char *A__){{
  double piv_;
if( (next_==NULL) ){ memset((ssize_t)sizeof(double)+(A_h-2)*A_i+A__,0,(ssize_t)sizeof(double)*2); (*(double *)(2*(ssize_t)sizeof(double)+(A_h-3)*A_i+A__))=0.;}else{
    ssize_t _2st;
_2st=(-2)*(ssize_t)sizeof(double);
  if(!(fread(_2st+(A_h-1)*A_i+A__,(ssize_t)sizeof(double)*(2-(-2)+1),1, next_ )==1&&fread(_2st+A_h*A_i+A__,(ssize_t)sizeof(double)*(2-(-2)+1),1, next_ )==1))ioerr( next_ );
  };
  
   {int i_=A_h-2  ;while(i_>=A_l){
     {int k_=2  ;while(k_>=1){
      char *j_3;char *j_4;
j_3=k_*(ssize_t)sizeof(double)+i_*A_i+A__;j_4=(i_+k_)*A_i+A__;
      piv_=(*(double *)(j_3));
       {int _t5= - (-2) ;do{{ j_3=-(ssize_t)sizeof(double)+j_3;j_4=-(ssize_t)sizeof(double)+j_4;  (*(double *)(j_3))-=piv_*(*(double *)(j_4) );}_t5--;}while(_t5>0);}
    k_-=1;};}
    piv_=1./(*(double *)(i_*A_i+A__));  (*(double *)(i_*A_i+A__))=piv_;
     {char *j_=(-2)*(ssize_t)sizeof(double)+i_*A_i+A__;int j_1=(-1)-(-2);  while(j_1>=0){ (*(double *)(j_))*=piv_;j_+=(ssize_t)sizeof(double);j_1--;};}
  i_-=1;};}
  if( (prev_==NULL) ){ ssize_t _7st;
_7st=(-2)*(ssize_t)sizeof(double);
memset(_7st+A_l*A_i+A__,0,(ssize_t)sizeof(double)*((-1)-(-2)+1)); (*(double *)((-2)*(ssize_t)sizeof(double)+(A_l+1)*A_i+A__))=0.;}else{
    ssize_t _8st;
_8st=(-2)*(ssize_t)sizeof(double);
fwrite(_8st+A_l*A_i+A__,(ssize_t)sizeof(double)*(2-(-2)+1),1,prev_);fwrite(_8st+(A_l+1)*A_i+A__,(ssize_t)sizeof(double)*(2-(-2)+1),1,prev_);
  };
}}

void rbmatmod_2LeftLUDivStep1(int x_l,int x_h,size_t x_i,char *x__,int A_l,int A_h,size_t A_i,char *A__,int t_l,int t_h,size_t t_i,char *t__){{
  if( !(next_==NULL) ){   if(!(fread(0+(double *)((x_h-1)*x_i+x__),(ssize_t)sizeof(double),1, next_ )==1&&fread(0+(double *)(x_h*x_i+x__),(ssize_t)sizeof(double),1, next_ )==1))ioerr( next_ );};
   {int i_=A_h-2  ;while(i_>=A_l){
    char *j_3;char *j_4;
double sum_;
j_3=2*(ssize_t)sizeof(double)+i_*A_i+A__;j_4=2*x_i+i_*x_i+x__;
    sum_=(*(double *)(i_*t_i+t__));
     {int _t5=2 ;while(_t5>0){ sum_-=(*(double *)(j_3))*(*(double *)(j_4)); j_3=-(ssize_t)sizeof(double)+j_3;j_4=-x_i+j_4;_t5--;};}
    (*(double *)(j_4))=sum_*(*(double *)(j_3));
  i_-=1;};}
  if( !(prev_==NULL) ){ fwrite((double*)((x_l+2)*x_i+x__),(ssize_t)sizeof(double),1,prev_);fwrite((double*)((x_l+3)*x_i+x__),(ssize_t)sizeof(double),1,prev_);};
}}

void rbmatmod_3LeftLUDivStep2(int x_l,int x_h,size_t x_i,char *x__,int A_l,int A_h,size_t A_i,char *A__){{
  if( !(prev_==NULL) ){   if(!(fread(0+(double *)(x_l*x_i+x__),(ssize_t)sizeof(double),1, prev_ )==1&&fread(0+(double *)((x_l+1)*x_i+x__),(ssize_t)sizeof(double),1, prev_ )==1))ioerr( prev_ );};
   {int i_=A_l  ;while(i_<=A_h){
    char *j_4;char *j_5;
double sum_;
j_4=(-2)*(ssize_t)sizeof(double)+i_*A_i+A__;j_5=(-2)*x_i+i_*x_i+x__;
    sum_=0.;
     {int _t6= - (-2) ;while(_t6>0){ sum_+=(*(double *)(j_4))*(*(double *)(j_5)); j_4=(ssize_t)sizeof(double)+j_4;j_5=x_i+j_5;_t6--;};}
    (*(double *)(j_5))-=sum_;
  i_+=1;};}
  if( !(next_==NULL) ){ fwrite((double*)((x_h-3)*x_i+x__),(ssize_t)sizeof(double),1,next_);fwrite((double*)((x_h-2)*x_i+x__),(ssize_t)sizeof(double),1,next_);};
}}

void rbmatmod_4RightLUDivStep1(int x_l,int x_h,size_t x_i,char *x__,int t_l,int t_h,size_t t_i,char *t__,int A_l,int A_h,size_t A_i,char *A__){{
  if( !(next_==NULL) ){   if(!(fread(0+(double *)((x_h-1)*x_i+x__),(ssize_t)sizeof(double),1, next_ )==1&&fread(0+(double *)(x_h*x_i+x__),(ssize_t)sizeof(double),1, next_ )==1))ioerr( next_ );};
   {int i_=A_h-2  ;while(i_>=A_l){
    double sum_;
sum_=(*(double *)(i_*t_i+t__));
     {int j_=(-2)  ;do{{ sum_-=(*(double *)(j_*(ssize_t)sizeof(double)+(i_-j_)*A_i+A__))*(*(double *)((i_-j_)*x_i+x__) );}j_+=1;}while(j_<= - 1);}
    (*(double *)(i_*x_i+x__))=sum_;
  i_-=1;};}
  if( !(prev_==NULL) ){ fwrite((double*)((x_l+2)*x_i+x__),(ssize_t)sizeof(double),1,prev_);fwrite((double*)((x_l+3)*x_i+x__),(ssize_t)sizeof(double),1,prev_);};
}}

void rbmatmod_5RightLUDivStep2(int x_l,int x_h,size_t x_i,char *x__,int A_l,int A_h,size_t A_i,char *A__){{
  if( !(prev_==NULL) ){   if(!(fread(0+(double *)(x_l*x_i+x__),(ssize_t)sizeof(double),1, prev_ )==1&&fread(0+(double *)((x_l+1)*x_i+x__),(ssize_t)sizeof(double),1, prev_ )==1&&fread(0+(double *)((x_l+2)*x_i+x__),(ssize_t)sizeof(double),1, prev_ )==1&&fread(0+(double *)((x_l+3)*x_i+x__),(ssize_t)sizeof(double),1, prev_ )==1))ioerr( prev_ );};
   {int i_=A_l  ;while(i_<=A_h-2){
    char *j_6;char *j_7;
double _8p;
/*p=_8p*/
j_6=i_*A_i+A__;j_7=i_*x_i+x__;
    _8p=(*(double *)(j_7))*(*(double *)(i_*A_i+A__));

    (*(double *)(j_7))=_8p;
     {int _t9=2 ;do{{ j_6=(ssize_t)sizeof(double)+j_6;j_7=x_i+j_7;  (*(double *)(j_7))-=(*(double *)(j_6))*_8p ;}_t9--;}while(_t9>0);}
  i_+=1;};}
  if( !(next_==NULL) ){
    fwrite((double*)((x_h-3)*x_i+x__),(ssize_t)sizeof(double),1,next_);fwrite((double*)((x_h-2)*x_i+x__),(ssize_t)sizeof(double),1,next_);fwrite((double*)((x_h-1)*x_i+x__),(ssize_t)sizeof(double),1,next_);fwrite((double*)(x_h*x_i+x__),(ssize_t)sizeof(double),1,next_);
    (*(double *)((x_h-1)*x_i+x__))*=(*(double *)((A_h-1)*A_i+A__));
    (*(double *)(x_h*x_i+x__))=((*(double *)(x_h*x_i+x__))-(*(double *)((ssize_t)sizeof(double)+(A_h-1)*A_i+A__))*(*(double *)((x_h-1)*x_i+x__)))*(*(double *)(A_h*A_i+A__));
  };  
}}
