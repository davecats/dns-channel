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
/**/
/*nsmp=4*/
int iproc_;
int nproc_;
/*bufsize=800*/
/*baseport=(IPPORT_USERRESERVED+111)*/
FILE *prev_;
FILE *next_;
void rbmatmod_1LUdecompStep(int A_l,int A_h,size_t A_i,char *A__);
void rbmatmod_2LeftLUDivStep1(int x_l,int x_h,size_t x_i,char *x__,int A_l,int A_h,size_t A_i,char *A__,int t_l,int t_h,size_t t_i,char *t__);
void rbmatmod_3LeftLUDivStep2(int x_l,int x_h,size_t x_i,char *x__,int A_l,int A_h,size_t A_i,char *A__);
void rbmatmod_4RightLUDivStep1(int x_l,int x_h,size_t x_i,char *x__,int t_l,int t_h,size_t t_i,char *t__,int A_l,int A_h,size_t A_i,char *A__);
void rbmatmod_5RightLUDivStep2(int x_l,int x_h,size_t x_i,char *x__,int A_l,int A_h,size_t A_i,char *A__);

int ny_;
int nx_;
int nz_;
double alfa0_;
double beta0_;
double a_;
double ymin_;
double ymax_;
double t_max_;
double dt_field_;
double dt_save_;
double u_conv_;
double u0_;
double un_;
double w_conv_;
double w0_;
double wn_;
double ni_;
double meanpx_;
double meanpz_;
double meanflowx_;
double meanflowz_;
double px_;
double corrpx_;
double pz_;
double corrpz_;
double flowx_;
double flowz_;
double deltat_;
double cflmax_;
double time_;
int time_from_restart_;
char *restart_file_;
void dnsdata_1read_initial_data(){{
 FILE *in_data_;
char c;
in_data_=NULL;errfopen(&in_data_,"dns.in",O_RDWR|O_CREAT);
    if(!(scanrec( in_data_ ," ny = %d",&ny_)>0&&scanrec( in_data_ ," nx = %d",&nx_)>0&&scanrec( in_data_ ," nz = %d",&nz_)>0&&scanrec( in_data_ ," alfa0 = %lg",&alfa0_)>0&&scanrec( in_data_ ," beta0 = %lg",&beta0_)>0&&scanrec( in_data_ ," ymin = %lg",&ymin_)>0&&scanrec( in_data_ ," ymax = %lg",&ymax_)>0&&scanrec( in_data_ ," a = %lg",&a_)>0&&scanrec( in_data_ ," ni = %lg",&ni_)>0))ioerr( in_data_ );  ni_=1./ni_;
 {do{  }while(   (scanrec( in_data_ ," meanpx = %lg",&meanpx_ )>0 ||scanrec( in_data_ ," meanflowx = %lg",&meanflowx_)>0));}
 {do{  }while(   (scanrec( in_data_ ," meanpz = %lg",&meanpz_ )>0 ||scanrec( in_data_ ," meanflowz = %lg",&meanflowz_)>0));}
    if(!(scanrec( in_data_ ," u_conv = %lg",&u_conv_)>0 &&scanrec( in_data_ ," w_conv = %lg",&w_conv_)>0))ioerr( in_data_ );
    if(!(scanrec( in_data_ ," u0 = %lg",&u0_)>0 &&scanrec( in_data_ ," un = %lg",&un_)>0 &&scanrec( in_data_ ," w0 = %lg",&w0_)>0 &&scanrec( in_data_ ," wn = %lg",&wn_)>0))ioerr( in_data_ );
 {do{  }while(   (scanrec( in_data_ ," deltat = %lg",&deltat_ )>0 ||scanrec( in_data_ ," cflmax = %lg",&cflmax_)>0));}
    if(!(scanrec( in_data_ ," t_max = %lg",&t_max_)>0 &&(scanrec( in_data_ ," time_from_restart = %c%*4[A-Za-z] ",&c)&&((time_from_restart_=(c=='T')||(c=='Y')||(c=='t')||(c=='y'))||(c=='F')||(c=='N')||(c=='f')||(c=='n'))) &&scanrec( in_data_ ," dt_field = %lg",&dt_field_)>0 &&scanrec( in_data_ ," dt_save = %lg",&dt_save_ )>0))ioerr( in_data_ );
 if( !   (mygetline(" restart_file =",&restart_file_ , in_data_ ))){ if(restart_file_)free(restart_file_);restart_file_=strndup("",(int)(strlen("")-1)+1);};
 if(in_data_)errfclose(&in_data_);
 if( (next_==NULL) ){
      fprintf(stdout,"""nproc=""");fprintf(stdout,"%d",nproc_);fprintf(stdout,"\t");fprintf(stdout,"""nsmp=""");fprintf(stdout,"%d",4);putc('\n',stdout);
      fprintf(stdout,"""nx=""");fprintf(stdout,"%d",nx_);fprintf(stdout,"\t" );fprintf(stdout,"""nz=""");fprintf(stdout,"%d",nz_);fprintf(stdout,"\t" );fprintf(stdout,"""ny=""");fprintf(stdout,"%d",ny_);fprintf(stdout,"\t" );fprintf(stdout,"""time=""");fprintf(stdout,"%g",time_);putc('\n',stdout);
      fprintf(stdout,"""meanflowx=""");fprintf(stdout,"%g",meanflowx_);fprintf(stdout,"\t" );fprintf(stdout,"""meanpx=""");fprintf(stdout,"%g",meanpx_);fprintf(stdout,"\t" );fprintf(stdout,"""meanflowz=""");fprintf(stdout,"%g",meanflowz_);fprintf(stdout,"\t" );fprintf(stdout,"""meanpz=""");fprintf(stdout,"%g",meanpz_);putc('\n',stdout);
      fprintf(stdout,"""ymin=""");fprintf(stdout,"%g",ymin_);fprintf(stdout,"\t" );fprintf(stdout,"""ymax=""");fprintf(stdout,"%g",ymax_);fprintf(stdout,"\t" );fprintf(stdout,"""a=""");fprintf(stdout,"%g",a_);fprintf(stdout,"\t" );fprintf(stdout,"""alfa0=""");fprintf(stdout,"%g",alfa0_);fprintf(stdout,"\t" );fprintf(stdout,"""beta0=""");fprintf(stdout,"%g",beta0_);fprintf(stdout,"\t" );fprintf(stdout,"""1/ni=""");fprintf(stdout,"%g",1./ni_);putc('\n',stdout);
      fprintf(stdout,"""u_conv=""");fprintf(stdout,"%g",u_conv_);fprintf(stdout,"\t" );fprintf(stdout,"""u0=""");fprintf(stdout,"%g",u0_);fprintf(stdout,"\t" );fprintf(stdout,"""un=""");fprintf(stdout,"%g",un_);fprintf(stdout,"\t" );fprintf(stdout,"""w_conv=""");fprintf(stdout,"%g",w_conv_);fprintf(stdout,"\t" );fprintf(stdout,"""w0=""");fprintf(stdout,"%g",w0_);fprintf(stdout,"\t" );fprintf(stdout,"""wn=""");fprintf(stdout,"%g",wn_);putc('\n',stdout);
      fprintf(stdout,"""deltat=""");fprintf(stdout,"%g",deltat_);fprintf(stdout,"\t" );fprintf(stdout,"""cflmax=""");fprintf(stdout,"%g",cflmax_);fprintf(stdout,"\t" );fprintf(stdout,"""t_max=""");fprintf(stdout,"%g",t_max_);fprintf(stdout,"\t" );fprintf(stdout,"""dt_save=""");fprintf(stdout,"%g",dt_save_);fprintf(stdout,"\t" );fprintf(stdout,"""dt_field=""");fprintf(stdout,"%g",dt_field_);putc('\n',stdout);
 };
if(in_data_)errfclose(&in_data_);}}

int dnsdata_2nyl;
/*nyl=dnsdata_2nyl*/
int dnsdata_3nyh;
/*nyh=dnsdata_3nyh*/
int dnsdata_4M;
int dnsdata_5l;
int dnsdata_6m;
int dnsdata_7h;
size_t dnsdata_8i;
ssize_t dnsdata_9st;
char *y_;
ssize_t dnsdata_11st;
struct dnsdata_s10{char d0_[(ssize_t)sizeof(double)*(2-(-2)+1)];char d1_[(ssize_t)sizeof(double)*(2-(-2)+1)];char d2_[(ssize_t)sizeof(double)*(2-(-2)+1)];char d4_[(ssize_t)sizeof(double)*(2-(-2)+1)];};int dnsdata_12M;
int dnsdata_13l;
int dnsdata_14m;
int dnsdata_15h;
size_t dnsdata_16i;
ssize_t dnsdata_17st;
char *derivatives_;
char d040_[(ssize_t)sizeof(double)*(2-(-2)+1)];
char d140_[(ssize_t)sizeof(double)*(2-(-2)+1)];
char d240_[(ssize_t)sizeof(double)*(2-(-2)+1)];
char d340_[(ssize_t)sizeof(double)*(2-(-2)+1)];
char d14m1_[(ssize_t)sizeof(double)*(2-(-2)+1)];
char d24m1_[(ssize_t)sizeof(double)*(2-(-2)+1)];
char d04n_[(ssize_t)sizeof(double)*(2-(-2)+1)];
char d14n_[(ssize_t)sizeof(double)*(2-(-2)+1)];
char d24n_[(ssize_t)sizeof(double)*(2-(-2)+1)];
char d14np1_[(ssize_t)sizeof(double)*(2-(-2)+1)];
char d24np1_[(ssize_t)sizeof(double)*(2-(-2)+1)];
char dnsdata_setup_derivatives_M[(ssize_t)sizeof(double)*(4+1)*(4+1)];
char dnsdata_setup_derivatives_t[(ssize_t)sizeof(double)*(4+1)];
int dnsdata_setup_derivatives_18M;
int dnsdata_setup_derivatives_19m;
int dnsdata_setup_derivatives_20;
int dnsdata_setup_derivatives_21;
int dnsdata_23l;
int dnsdata_24h;
size_t dnsdata_25i;
size_t dnsdata_26i;
ssize_t dnsdata_27st;
ssize_t dnsdata_28st;
char *D2vmat_;
char *etamat_;
char *D0mat_;
ssize_t dnsdata_31st;
char v0bc_[(ssize_t)sizeof(double)*(2-(-2)+1)*(-(-1)+1)];
char vnbc_[(ssize_t)sizeof(double)*(2-(-2)+1)*(1+1)];
char eta0bc_[(ssize_t)sizeof(double)*(2-(-2)+1)*(-(-1)+1)];
char etanbc_[(ssize_t)sizeof(double)*(2-(-2)+1)*(1+1)];
void dnsdata_32applybc_0(int eq_l,int eq_h,size_t eq_i,int eq__l,int eq__h,size_t eq__i,char *eq___,char *bc_) {{
   {int i_= - 1  ;do{{ (*(double *)(i_*eq__i+eq_i+eq___))-=(*(double *)((-2)*eq__i+eq_i+eq___))*(*(double *)(i_*(ssize_t)sizeof(double)-(ssize_t)sizeof(double)*(2-(-2)+1)+bc_))/(*(double *)((-2)*(ssize_t)sizeof(double)-(ssize_t)sizeof(double)*(2-(-2)+1)+bc_)) ;}i_+=1;}while(i_<=2);}
   {int i_=0  ;do{{ (*(double *)(i_*eq__i+eq_i+eq___))-=(*(double *)(-eq__i+eq_i+eq___))*(*(double *)(i_*(ssize_t)sizeof(double)+bc_))/(*(double *)(-(ssize_t)sizeof(double)+bc_)) ;}i_+=1;}while(i_<=2);}
   {int i_=0  ;do{{ (*(double *)((i_-1)*eq__i+2*eq_i+eq___))-=(*(double *)((-2)*eq__i+2*eq_i+eq___))*(*(double *)(i_*(ssize_t)sizeof(double)+bc_))/(*(double *)(-(ssize_t)sizeof(double)+bc_)) ;}i_+=1;}while(i_<=2);}
}}
void dnsdata_33applybc_n(int eq_l,int eq_h,size_t eq_i,int eq__l,int eq__h,size_t eq__i,char *eq___,char *bc_) {{
   {int i_= - 2  ;do{{ (*(double *)(i_*eq__i+(ny_-1)*eq_i+eq___))-=(*(double *)(2*eq__i+(ny_-1)*eq_i+eq___))*(*(double *)(i_*(ssize_t)sizeof(double)+(ssize_t)sizeof(double)*(2-(-2)+1)+bc_))/(*(double *)(2*(ssize_t)sizeof(double)+(ssize_t)sizeof(double)*(2-(-2)+1)+bc_)) ;}i_+=1;}while(i_<=1);}
   {int i_= - 2  ;do{{ (*(double *)(i_*eq__i+(ny_-1)*eq_i+eq___))-=(*(double *)(eq__i+(ny_-1)*eq_i+eq___))*(*(double *)(i_*(ssize_t)sizeof(double)+bc_))/(*(double *)((ssize_t)sizeof(double)+bc_)) ;}i_+=1;}while(i_<=0);}
   {int i_= - 2  ;do{{ (*(double *)((i_+1)*eq__i+(ny_-2)*eq_i+eq___))-=(*(double *)(2*eq__i+(ny_-2)*eq_i+eq___))*(*(double *)(i_*(ssize_t)sizeof(double)+bc_))/(*(double *)((ssize_t)sizeof(double)+bc_)) ;}i_+=1;}while(i_<=0);}
}}

/*! Integral in y direction*/
void dnsdata_34yintegr(double *RESULT_,int f_l,int f_h,size_t f_i,char *f__){{
   {int iy_=((dnsdata_2nyl / 2))*2+1  ;while(iy_<=dnsdata_3nyh ){
   double _35yp1;
/*yp1=_35yp1*/
double _36ym1;
/*ym1=_36ym1*/
double _37a1;
/*a1=_37a1*/
double _38a3;
/*a3=_38a3*/
double _39a2;
/*a2=_39a2*/
_35yp1=(*(double *)((iy_+1)*(ssize_t)sizeof(double)+y_))-(*(double *)(iy_*(ssize_t)sizeof(double)+y_));
  _36ym1=(*(double *)((iy_-1)*(ssize_t)sizeof(double)+y_))-(*(double *)(iy_*(ssize_t)sizeof(double)+y_) );

   _37a1=-1./3.*_36ym1+1./6.*_35yp1+1./6.*_35yp1*_35yp1/_36ym1;

   _38a3=1./3.*_35yp1-1./6.*_36ym1-1./6.*_36ym1*_36ym1/_35yp1;

   _39a2=_35yp1-_36ym1-_37a1-_38a3;

   (*RESULT_)+=_37a1*(*(double *)((iy_-1)*f_i+f__)) + _39a2*(*(double *)(iy_*f_i+f__)) + _38a3*(*(double *)((iy_+1)*f_i+f__))    ;
  iy_+= 2 ;};}
}}

struct VELOCITY_{struct COMPLEX_ u_;struct COMPLEX_ v_;struct COMPLEX_ w_;};struct MOMFLUX_{struct COMPLEX_ uu_;struct COMPLEX_ uv_;struct COMPLEX_ vv_;struct COMPLEX_ vw_;struct COMPLEX_ ww_;struct COMPLEX_ uw_;};/*!INLINE FUNCTION OS(INTEGER iy,i)=0.5*[d4(i)-2*k2*d2(i)+k2*k2*d0(i)] !Vittori*/
/*!INLINE FUNCTION SQ(INTEGER iy,i)=0.5*[d2(i)-k2*d0(i)]               !Vittori*/
int nxd_;
int nzd_;
int dnsdata_35h;
int dnsdata_36h;
size_t dnsdata_37i;
size_t dnsdata_38i;
char *Vd_;
int dnsdata_39h;
int dnsdata_40h;
size_t dnsdata_41i;
size_t dnsdata_42i;
char *VVd_;
/*maxtimelevels=1*/
struct rhstype_{struct COMPLEX_ eta_;struct COMPLEX_ D2v_;};struct VETA_{struct COMPLEX_ v_;struct COMPLEX_ eta_;};int dnsdata_43h;
int dnsdata_44l;
int dnsdata_45h;
size_t dnsdata_46i;
size_t dnsdata_47i;
size_t dnsdata_48i;
ssize_t dnsdata_49st;
ssize_t dnsdata_50st;
char *V_;
int dnsdata_51h;
int dnsdata_52l;
int dnsdata_53h;
int dnsdata_54M;
int dnsdata_55l;
int dnsdata_56m;
int dnsdata_57h;
size_t dnsdata_58i;
size_t dnsdata_59i;
size_t dnsdata_60i;
ssize_t dnsdata_61st;
ssize_t dnsdata_62st;
char *oldrhs_;
int dnsdata_63h;
int dnsdata_64l;
int dnsdata_65h;
size_t dnsdata_66i;
size_t dnsdata_67i;
size_t dnsdata_68i;
ssize_t dnsdata_69st;
char *memrhs_;
double u1zero_;
double w1zero_;
int ismp_;
int dnsdata_70h;
int dnsdata_71h;
int dnsdata_72l;
int dnsdata_73h;
off_t dnsdata_74i;
off_t dnsdata_75i;
off_t dnsdata_77i;
off_t dnsdata_78i;
off_t dnsdata_79i;
ssize_t dnsdata_80st;
ssize_t dnsdata_81st;
struct dnsdata_s76{char header_[(1023+1)];};FILE *diskimage_;
int dnsdata_82h;
int dnsdata_83h;
int dnsdata_84h;
int dnsdata_85l;
int dnsdata_86h;
off_t dnsdata_87i;
off_t dnsdata_88i;
off_t dnsdata_90i;
ssize_t dnsdata_91st;
off_t dnsdata_92i;
off_t dnsdata_93i;
off_t dnsdata_94i;
ssize_t dnsdata_95st;
ssize_t dnsdata_96st;
struct dnsdata_s89{int nyimage_;int nximage_;int nzimage_;double timage_;double yminimage_;double ymaximage_;double aimage_;double alfa0image_;double beta0image_;double niimage_;};FILE *diskfield_;
int cont_;
int outcont_;
FILE *time_file_;
int miny_;
int maxy_;
double cfl_;
double cflm_;
void dnsdata_97getcfl(){{
/*!nx OPPURE nxd?? Check tesi Ferro*/
double _98dx;
/*dx=_98dx*/
double _99dz;
/*dz=_99dz*/
_98dx=3.14159265358979323846/(alfa0_*(double)(nxd_));
  _99dz=2.*3.14159265358979323846/(beta0_*(double)(nzd_));

cfl_=0.;
 {int iy_=dnsdata_2nyl  ;while(iy_<=dnsdata_3nyh){
  double dy_;
double _127;
double _128;
double _129;
double _130M;
dy_=0.5*((*(double *)((iy_+1)*(ssize_t)sizeof(double)+y_))-(*(double *)((iy_-1)*(ssize_t)sizeof(double)+y_)));
  fflush(NULL); ismp_=0  ;while(ismp_<4-1&&fork()){ismp_+=1;;};{

     {int ix_=ismp_  ;while(ix_<=nx_ ){
      int _100h;
int _101h;
int _104l;
int _105h;
ssize_t _106st;
size_t _107i;
int _108l;
int _109l;
_100h=nz_;
_101h=nz_;
{char *_103_;
_103_=iy_*(ssize_t)sizeof(struct VELOCITY_)+ix_*dnsdata_47i+V_;
 {char *_102_=ix_*dnsdata_37i+Vd_;int _102_1=_100h; do{{ (*(struct VELOCITY_ *)(_102_))=(*(struct VELOCITY_ *)(_103_)); _103_ =dnsdata_46i+_103_;}_102_+=(ssize_t)sizeof(struct VELOCITY_);_102_1--;}while(_102_1>=0);}}
      _104l=nz_+1;
_105h=nzd_-nz_-1;
_106st=_104l*(ssize_t)sizeof(struct VELOCITY_);
_107i=(ssize_t)sizeof(struct VELOCITY_)*(_105h-_104l+1);
memset(_106st+ix_*dnsdata_37i+Vd_,0,_107i);
      _108l= - nz_;
_109l= - nz_;
{char *_111_;
_111_=_108l*dnsdata_46i+iy_*(ssize_t)sizeof(struct VELOCITY_)+ix_*dnsdata_47i+V_;
 {char *_110_=_108l*(ssize_t)sizeof(struct VELOCITY_)+nzd_*(ssize_t)sizeof(struct VELOCITY_)+ix_*dnsdata_37i+Vd_;int _110_1=(-1)-_108l; do{{ (*(struct VELOCITY_ *)(_110_))=(*(struct VELOCITY_ *)(_111_)); _111_ =dnsdata_46i+_111_;}_110_+=(ssize_t)sizeof(struct VELOCITY_);_110_1--;}while(_110_1>=0);}}
      { char *_112w;
_112w=ix_*dnsdata_37i+Vd_;
{fft_1IFT(0,dnsdata_36h,(ssize_t)sizeof(struct VELOCITY_),((char*)&(*(struct VELOCITY_ *)(_112w)).u_),0,dnsdata_36h,(ssize_t)sizeof(struct VELOCITY_),((char*)&(*(struct VELOCITY_ *)(_112w)).u_));}; {fft_1IFT(0,dnsdata_36h,(ssize_t)sizeof(struct VELOCITY_),((char*)&(*(struct VELOCITY_ *)(_112w)).v_),0,dnsdata_36h,(ssize_t)sizeof(struct VELOCITY_),((char*)&(*(struct VELOCITY_ *)(_112w)).v_));}; {fft_1IFT(0,dnsdata_36h,(ssize_t)sizeof(struct VELOCITY_),((char*)&(*(struct VELOCITY_ *)(_112w)).w_),0,dnsdata_36h,(ssize_t)sizeof(struct VELOCITY_),((char*)&(*(struct VELOCITY_ *)(_112w)).w_));};}
    ix_+= 4;};}
    if( ismp_==0 ){ int _113l;
int _114h;
ssize_t _115st;
size_t _116i;
_113l=nx_+1;
_114h=nxd_-1;
_115st=_113l*dnsdata_37i;
_116i=dnsdata_37i*(_114h-_113l+1);
memset(_115st+Vd_,0,_116i);};
    {
  int md;
while(!((*barrier_)==(ismp_)))sigsuspend(&oldmask);
  (*barrier_)=(md=((ismp_)+1) % (4),md>=0?md:md+(4));  kill(0,SIGUSR1);
  while(!((*barrier_)<=(ismp_)))sigsuspend(&oldmask);
};
     {int iz_=ismp_  ;do{{ { char *_117w;
_117w=iz_*(ssize_t)sizeof(struct VELOCITY_)+Vd_;
{fft_3RFT(0,dnsdata_35h,dnsdata_37i,((char*)&(*(struct VELOCITY_ *)(_117w)).u_),0,dnsdata_35h,dnsdata_37i,((char*)&(*(struct VELOCITY_ *)(_117w)).u_));}; {fft_3RFT(0,dnsdata_35h,dnsdata_37i,((char*)&(*(struct VELOCITY_ *)(_117w)).v_),0,dnsdata_35h,dnsdata_37i,((char*)&(*(struct VELOCITY_ *)(_117w)).v_));}; {fft_3RFT(0,dnsdata_35h,dnsdata_37i,((char*)&(*(struct VELOCITY_ *)(_117w)).w_),0,dnsdata_35h,dnsdata_37i,((char*)&(*(struct VELOCITY_ *)(_117w)).w_));} ;}}iz_+= 4;}while(iz_<=dnsdata_36h );}
    {
  int md;
while(!((*barrier_)==(ismp_)))sigsuspend(&oldmask);
  (*barrier_)=(md=((ismp_)+1) % (4),md>=0?md:md+(4));  kill(0,SIGUSR1);
  while(!((*barrier_)<=(ismp_)))sigsuspend(&oldmask);
};
  } if(ismp_<4-1)exit(0);;
 ismp_=0  ;while(ismp_<4-1){if(wait(NULL)<0){strcpy(ERRORMESSAGE_,"wait"); longjmp(errorhook,1);};ismp_+=1;;};
  /*! Un campione ogni due*/
   _127=-DBL_MAX;  {int _118i_=dnsdata_35h;do{{int _119i_=dnsdata_36h;do{{double _120M;
_120M=fabs((*(double *)(_119i_*(ssize_t)sizeof(struct VELOCITY_)+_118i_*dnsdata_37i+(char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(Vd_)).u_)).IMAG_))) ;
if(_127<_120M)_127=_120M;}_119i_--;}while(_119i_>=0);}_118i_--;}while(_118i_>=0);} _128=-DBL_MAX;  {int _121i_=dnsdata_35h;do{{int _122i_=dnsdata_36h;do{{double _123M;
_123M=fabs((*(double *)(_122i_*(ssize_t)sizeof(struct VELOCITY_)+_121i_*dnsdata_37i+(char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(Vd_)).v_)).IMAG_))) ;
if(_128<_123M)_128=_123M;}_122i_--;}while(_122i_>=0);}_121i_--;}while(_121i_>=0);} _129=-DBL_MAX;  {int _124i_=dnsdata_35h;do{{int _125i_=dnsdata_36h;do{{double _126M;
_126M=fabs((*(double *)(_125i_*(ssize_t)sizeof(struct VELOCITY_)+_124i_*dnsdata_37i+(char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(Vd_)).w_)).IMAG_))) ;
if(_129<_126M)_129=_126M;}_125i_--;}while(_125i_>=0);}_124i_--;}while(_124i_>=0);}_130M=_127/_98dx+_128/dy_+_129/_99dz;
cfl_=cfl_;if(cfl_<_130M)cfl_=_130M;
iy_+=1;};}
}}

double energy_;
double slice_energy_;
double diss_;
double slice_diss_; 
struct COMPLEX_ dudy_;
struct COMPLEX_ dvdy_;
struct COMPLEX_ dwdy_;
int dnsdata_98h;
int dnsdata_99l;
int dnsdata_100h;
size_t dnsdata_101i;
size_t dnsdata_102i;
ssize_t dnsdata_103st;
char *fieldbuf_;
int dnsdata_104h;
int dnsdata_105l;
int dnsdata_106h;
size_t dnsdata_107i;
size_t dnsdata_108i;
ssize_t dnsdata_109st;
char *velbuf_;
void dnsdata_110outstats(){{
outcont_+=1;  
if( outcont_>0 ){
	outcont_=0;  dnsdata_97getcfl();  cflm_=0.;  energy_=0.;  slice_energy_=0.;  diss_=0.;  slice_diss_=0.;

	 {int iy_=dnsdata_2nyl  ;do{  {int ix_=dnsdata_43h;do{{int iz_=dnsdata_44l;do{ if(!( (ix_==0 )&&( iz_==0))){ { 
	struct VELOCITY_ *_111w;
_111w=(struct VELOCITY_ *)(iy_*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_);
 slice_energy_ +=  1./2.* (((*_111w).u_.REAL_*(*_111w).u_.REAL_+(*_111w).u_.IMAG_*(*_111w).u_.IMAG_)+((*_111w).v_.REAL_*(*_111w).v_.REAL_+(*_111w).v_.IMAG_*(*_111w).v_.IMAG_)+((*_111w).w_.REAL_*(*_111w).w_.REAL_+(*_111w).w_.IMAG_*(*_111w).w_.IMAG_))*(double)((ix_==0?1:2))*((*(double *)((iy_+1)*(ssize_t)sizeof(double)+y_))-(*(double *)((iy_-1)*(ssize_t)sizeof(double)+y_)))*0.5 ;}}iz_++;}while(iz_<=dnsdata_45h);}ix_--;}while(ix_>=0);}iy_+=1;}while(iy_<=dnsdata_3nyh );}

	 {int iy_=dnsdata_2nyl  ;while(iy_<=dnsdata_3nyh){
		 {int ix_=1  ;while(ix_<=nx_ ){  {int iz_=dnsdata_44l;do{ {
			char *_112w;
double _113alfa;
/*alfa=_113alfa*/
double _114beta;
/*beta=_114beta*/
double _115r;
double _116i;
double _117r;
double _118r;
double _119i;
double _120r;
double _121r;
double _122i;
double _123r;
double _124i;
double _125r;
double _126r;
double _127i;
double _128r;
double _129r;
double _130i;
double _131r;
double _132i;
double _133r;
double _134r;
double _135i;
double _136r;
double _137r;
double _138i;
double _139r;
double _140i;
double _141r;
double _142i;
double _143r;
double _144i;
double _145r;
double _146i;
double _147r;
double _148i;
double _149r;
double _150i;
_112w=iz_*dnsdata_46i+ix_*dnsdata_47i+V_;
_113alfa=(double)(ix_)*alfa0_;
  _114beta=(double)(iz_)*beta0_;

			_115r=((*(struct COMPLEX_*)((iy_+1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).u_)).REAL_-((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).u_)).REAL_));
_116i=((*(struct COMPLEX_*)((iy_+1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).u_)).IMAG_-((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).u_)).IMAG_));
_117r=(1./(((*(double *)((iy_+1)*(ssize_t)sizeof(double)+y_))-(*(double *)(iy_*(ssize_t)sizeof(double)+y_))) ));
_118r=((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).u_)).REAL_-((*(struct COMPLEX_*)((iy_-1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).u_)).REAL_));
_119i=((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).u_)).IMAG_-((*(struct COMPLEX_*)((iy_-1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).u_)).IMAG_));
_120r=(1./(((*(double *)(iy_*(ssize_t)sizeof(double)+y_))-(*(double *)((iy_-1)*(ssize_t)sizeof(double)+y_)))));
_121r=((_115r*_117r)+(_118r*_120r));
_122i=((_116i*_117r)+(_119i*_120r));
 {register double temp=(0.5*_121r);dudy_ .IMAG_=(0.5*_122i);dudy_ .REAL_=temp;};
			_123r=((*(struct COMPLEX_*)((iy_+1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).v_)).REAL_-((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).v_)).REAL_));
_124i=((*(struct COMPLEX_*)((iy_+1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).v_)).IMAG_-((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).v_)).IMAG_));
_125r=(1./(((*(double *)((iy_+1)*(ssize_t)sizeof(double)+y_))-(*(double *)(iy_*(ssize_t)sizeof(double)+y_))) ));
_126r=((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).v_)).REAL_-((*(struct COMPLEX_*)((iy_-1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).v_)).REAL_));
_127i=((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).v_)).IMAG_-((*(struct COMPLEX_*)((iy_-1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).v_)).IMAG_));
_128r=(1./(((*(double *)(iy_*(ssize_t)sizeof(double)+y_))-(*(double *)((iy_-1)*(ssize_t)sizeof(double)+y_)))));
_129r=((_123r*_125r)+(_126r*_128r));
_130i=((_124i*_125r)+(_127i*_128r));
 {register double temp=(0.5*_129r);dvdy_ .IMAG_=(0.5*_130i);dvdy_ .REAL_=temp;};
			_131r=((*(struct COMPLEX_*)((iy_+1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).w_)).REAL_-((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).w_)).REAL_));
_132i=((*(struct COMPLEX_*)((iy_+1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).w_)).IMAG_-((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).w_)).IMAG_));
_133r=(1./(((*(double *)((iy_+1)*(ssize_t)sizeof(double)+y_))-(*(double *)(iy_*(ssize_t)sizeof(double)+y_))) ));
_134r=((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).w_)).REAL_-((*(struct COMPLEX_*)((iy_-1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).w_)).REAL_));
_135i=((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).w_)).IMAG_-((*(struct COMPLEX_*)((iy_-1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).w_)).IMAG_));
_136r=(1./(((*(double *)(iy_*(ssize_t)sizeof(double)+y_))-(*(double *)((iy_-1)*(ssize_t)sizeof(double)+y_)))));
_137r=((_131r*_133r)+(_134r*_136r));
_138i=((_132i*_133r)+(_135i*_136r));
 {register double temp=(0.5*_137r);dwdy_ .IMAG_=(0.5*_138i);dwdy_ .REAL_=temp;};
			_139r=(_113alfa*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).u_)).REAL_);
_140i=(_113alfa*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).u_)).IMAG_);
_141r=(_113alfa*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).v_)).REAL_);
_142i=(_113alfa*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).v_)).IMAG_);
_143r=(_113alfa*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).w_)).REAL_);
_144i=(_113alfa*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).w_)).IMAG_);
_145r=(_114beta*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).u_)).REAL_);
_146i=(_114beta*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).u_)).IMAG_);
_147r=(_114beta*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).v_)).REAL_);
_148i=(_114beta*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).v_)).IMAG_);
_149r=(_114beta*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).w_)).REAL_);
_150i=(_114beta*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_112w)).w_)).IMAG_);
 slice_diss_ +=  2.* 0.5*( (_139r*_139r+_140i*_140i)+(_141r*_141r+_142i*_142i)+(_143r*_143r+_144i*_144i) + 
                                                   (_145r*_145r+_146i*_146i)+(_147r*_147r+_148i*_148i)+(_149r*_149r+_150i*_150i) + 
                                                   (dudy_.REAL_*dudy_.REAL_+dudy_.IMAG_*dudy_.IMAG_)+(dvdy_.REAL_*dvdy_.REAL_+dvdy_.IMAG_*dvdy_.IMAG_)+(dwdy_.REAL_*dwdy_.REAL_+dwdy_.IMAG_*dwdy_.IMAG_) ) * ((*(double *)((iy_+1)*(ssize_t)sizeof(double)+y_))-(*(double *)((iy_-1)*(ssize_t)sizeof(double)+y_)))*0.5;
		}iz_++;}while(iz_<=dnsdata_45h);}ix_+=1;};}
		  {int iz_=dnsdata_44l;do{ if(!(iz_==0 )){{
                        char *_151w;
double _152beta;
/*beta=_152beta*/
double _153r;
double _154i;
double _155r;
double _156r;
double _157i;
double _158r;
double _159r;
double _160i;
double _161r;
double _162i;
double _163r;
double _164r;
double _165i;
double _166r;
double _167r;
double _168i;
double _169r;
double _170i;
double _171r;
double _172r;
double _173i;
double _174r;
double _175r;
double _176i;
double _177r;
double _178i;
double _179r;
double _180i;
double _181r;
double _182i;
_151w=iz_*dnsdata_46i+V_;
_152beta=(double)(iz_)*beta0_;

			_153r=((*(struct COMPLEX_*)((iy_+1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).u_)).REAL_-((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).u_)).REAL_));
_154i=((*(struct COMPLEX_*)((iy_+1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).u_)).IMAG_-((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).u_)).IMAG_));
_155r=(1./(((*(double *)((iy_+1)*(ssize_t)sizeof(double)+y_))-(*(double *)(iy_*(ssize_t)sizeof(double)+y_))) ));
_156r=((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).u_)).REAL_-((*(struct COMPLEX_*)((iy_-1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).u_)).REAL_));
_157i=((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).u_)).IMAG_-((*(struct COMPLEX_*)((iy_-1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).u_)).IMAG_));
_158r=(1./(((*(double *)(iy_*(ssize_t)sizeof(double)+y_))-(*(double *)((iy_-1)*(ssize_t)sizeof(double)+y_)))));
_159r=((_153r*_155r)+(_156r*_158r));
_160i=((_154i*_155r)+(_157i*_158r));
 {register double temp=(0.5*_159r);dudy_ .IMAG_=(0.5*_160i);dudy_ .REAL_=temp;};
			_161r=((*(struct COMPLEX_*)((iy_+1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).v_)).REAL_-((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).v_)).REAL_));
_162i=((*(struct COMPLEX_*)((iy_+1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).v_)).IMAG_-((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).v_)).IMAG_));
_163r=(1./(((*(double *)((iy_+1)*(ssize_t)sizeof(double)+y_))-(*(double *)(iy_*(ssize_t)sizeof(double)+y_))) ));
_164r=((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).v_)).REAL_-((*(struct COMPLEX_*)((iy_-1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).v_)).REAL_));
_165i=((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).v_)).IMAG_-((*(struct COMPLEX_*)((iy_-1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).v_)).IMAG_));
_166r=(1./(((*(double *)(iy_*(ssize_t)sizeof(double)+y_))-(*(double *)((iy_-1)*(ssize_t)sizeof(double)+y_)))));
_167r=((_161r*_163r)+(_164r*_166r));
_168i=((_162i*_163r)+(_165i*_166r));
 {register double temp=(0.5*_167r);dvdy_ .IMAG_=(0.5*_168i);dvdy_ .REAL_=temp;};
			_169r=((*(struct COMPLEX_*)((iy_+1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).w_)).REAL_-((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).w_)).REAL_));
_170i=((*(struct COMPLEX_*)((iy_+1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).w_)).IMAG_-((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).w_)).IMAG_));
_171r=(1./(((*(double *)((iy_+1)*(ssize_t)sizeof(double)+y_))-(*(double *)(iy_*(ssize_t)sizeof(double)+y_))) ));
_172r=((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).w_)).REAL_-((*(struct COMPLEX_*)((iy_-1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).w_)).REAL_));
_173i=((*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).w_)).IMAG_-((*(struct COMPLEX_*)((iy_-1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).w_)).IMAG_));
_174r=(1./(((*(double *)(iy_*(ssize_t)sizeof(double)+y_))-(*(double *)((iy_-1)*(ssize_t)sizeof(double)+y_)))));
_175r=((_169r*_171r)+(_172r*_174r));
_176i=((_170i*_171r)+(_173i*_174r));
 {register double temp=(0.5*_175r);dwdy_ .IMAG_=(0.5*_176i);dwdy_ .REAL_=temp;};
			_177r=(_152beta*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).u_)).REAL_);
_178i=(_152beta*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).u_)).IMAG_);
_179r=(_152beta*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).v_)).REAL_);
_180i=(_152beta*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).v_)).IMAG_);
_181r=(_152beta*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).w_)).REAL_);
_182i=(_152beta*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_151w)).w_)).IMAG_);
 slice_diss_ +=  0.5*( (_177r*_177r+_178i*_178i)+(_179r*_179r+_180i*_180i)+(_181r*_181r+_182i*_182i) + 
                                               (dudy_.REAL_*dudy_.REAL_+dudy_.IMAG_*dudy_.IMAG_)+(dvdy_.REAL_*dvdy_.REAL_+dvdy_.IMAG_*dvdy_.IMAG_)+(dwdy_.REAL_*dwdy_.REAL_+dwdy_.IMAG_*dwdy_.IMAG_) ) * ((*(double *)((iy_+1)*(ssize_t)sizeof(double)+y_))-(*(double *)((iy_-1)*(ssize_t)sizeof(double)+y_)))*0.5;
		}}iz_++;}while(iz_<=dnsdata_45h);}
	iy_+=1;};}
	if( !(prev_==NULL) ){ 
		  if(!(fread(0+&energy_,(ssize_t)sizeof(double),1, prev_ )==1&& fread(0+&cflm_,(ssize_t)sizeof(double),1, prev_ )==1&& fread(0+&diss_,(ssize_t)sizeof(double),1, prev_ )==1))ioerr( prev_ );  fflush(prev_);
		energy_+=slice_energy_;  diss_+=slice_diss_;  /*!cflm=MAX(cfl,cflm) */
		if( cfl_ >cflm_ ){  cflm_ = cfl_; };
  	}else{ 
		energy_=slice_energy_;  diss_=slice_diss_;  cflm_=cfl_;
	};
	if( !(next_==NULL) ){ fwrite(&energy_,(ssize_t)sizeof(double),1,next_); fwrite(&cflm_,(ssize_t)sizeof(double),1,next_); fwrite(&diss_,(ssize_t)sizeof(double),1,next_); fflush(next_);};
	if( cflmax_>0.){ 
		deltat_=cflmax_/cflm_;
		if( !(next_==NULL) ){   if(!(fread(0+&deltat_,(ssize_t)sizeof(double),1, next_ )==1))ioerr( next_ ); fflush(next_);};
		if( !(prev_==NULL) ){ fwrite(&deltat_,(ssize_t)sizeof(double),1,prev_); fflush(prev_);};
	};

	if( (next_==NULL) ){
	  double _183;
double _184;
double _185;
double _186;
 _183=0.; {int i_= - 2  ;do{{(*&_183)+= - (*(double *)(i_*(ssize_t)sizeof(double)+d14n_-((-2)*(ssize_t)sizeof(double))))*(*(struct VELOCITY_ *)((ny_-1+i_)*(ssize_t)sizeof(struct VELOCITY_)+V_)).u_.REAL_ ;}i_+=1;}while(i_<=2 );} _184=0.; {int i_= - 2  ;do{{(*&_184)+= - (*(double *)(i_*(ssize_t)sizeof(double)+d14n_-((-2)*(ssize_t)sizeof(double))))*(*(struct VELOCITY_ *)((ny_-1+i_)*(ssize_t)sizeof(struct VELOCITY_)+V_)).w_.REAL_ ;}i_+=1;}while(i_<=2 );} fprintf(stdout,"%1.9g\t%1.9g\t%1.9g\t%1.9g\t%1.9g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g",time_ ,u1zero_ ,_183
        	,w1zero_  ,_184
		,flowx_,px_+corrpx_,flowz_,pz_+corrpz_,cflm_*deltat_,deltat_,energy_,diss_);putc('\n',stdout);
	   _185=0.; {int i_= - 2  ;do{{(*&_185)+= - (*(double *)(i_*(ssize_t)sizeof(double)+d14n_-((-2)*(ssize_t)sizeof(double))))*(*(struct VELOCITY_ *)((ny_-1+i_)*(ssize_t)sizeof(struct VELOCITY_)+V_)).u_.REAL_ ;}i_+=1;}while(i_<=2 );} _186=0.; {int i_= - 2  ;do{{(*&_186)+= - (*(double *)(i_*(ssize_t)sizeof(double)+d14n_-((-2)*(ssize_t)sizeof(double))))*(*(struct VELOCITY_ *)((ny_-1+i_)*(ssize_t)sizeof(struct VELOCITY_)+V_)).w_.REAL_ ;}i_+=1;}while(i_<=2 );} fprintf( time_file_ ,"%1.9g\t%1.9g\t%1.9g\t%1.9g\t%1.9g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g",time_ 
        	,u1zero_ ,_185
	        ,w1zero_  ,_186
	        ,flowx_,px_+corrpx_,flowz_,pz_+corrpz_,cflm_*deltat_,deltat_,energy_,diss_);putc('\n', time_file_ );  fflush(time_file_);
	};
};
 
   	
if( (time_>0.)&&( (int)floor((time_+0.5*deltat_)/dt_save_) >(int)floor((time_-0.5*deltat_)/dt_save_) )){
    {int iy_=(dnsdata_2nyl-2);do{{ { struct VELOCITY_ *_187w;
_187w=(struct VELOCITY_ *)(iy_*(ssize_t)sizeof(struct VELOCITY_)+V_);
{register double temp=(*_187w).u_.REAL_-u_conv_;(*_187w).u_.IMAG_=(*_187w).u_.IMAG_;(*_187w).u_.REAL_=temp;}; {register double temp=(*_187w).w_.REAL_-w_conv_;(*_187w).w_.IMAG_=(*_187w).w_.IMAG_;(*_187w).w_.REAL_=temp;};}}iy_++;}while(iy_<=(dnsdata_3nyh+2));}
  if( !(prev_==NULL) ){  while(!feof( prev_)&&getc( prev_)!='\n'){};};
  errfopen(&diskimage_ ,"Dati.cart.out",O_RDWR|O_CREAT); 
  {
    if( (next_==NULL) ){
      char *_188; fprintf(stdout,"""Writing Dati.cart.out at time""\t%g" ,time_);putc('\n',stdout);
      fseeko(diskimage_,(-(0))+(size_t)(*(struct dnsdata_s76 *)(0)).header_ ,SEEK_SET);_188=malloc(snprintf(NULL,0,"         ny=""%d""       nx=""%d""       nz=""%d""\n"
"	 alfa0=""%g""     beta0=""%g""\n"
"         ymin=""%g""        ymax=""%g""          a=""%g""\n"
"         ni=""%g""       time=""%g""\n"
"         """,ny_,nx_,nz_,alfa0_,beta0_,ymin_,ymax_,a_,1./ni_,time_)+1);sprintf(_188,"         ny=""%d""       nx=""%d""       nz=""%d""\n"
"	 alfa0=""%g""     beta0=""%g""\n"
"         ymin=""%g""        ymax=""%g""          a=""%g""\n"
"         ni=""%g""       time=""%g""\n"
"         """,ny_,nx_,nz_,alfa0_,beta0_,ymin_,ymax_,a_,1./ni_,time_); fprintf( diskimage_,"%s",_188);putc('\n', diskimage_);
    free(_188);};

    {int iy_=  miny_  ;while(iy_<=maxy_){
       {int iz_=dnsdata_105l;do{{ int _189h;
int _190h;
_189h=nx_;
_190h=nx_;
{char *_192_;
_192_=iy_*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+V_;
 {char *_191_=iz_*(ssize_t)sizeof(struct VELOCITY_)+velbuf_;int _191_1=_189h; do{{ (*(struct VELOCITY_ *)(_191_))=(*(struct VELOCITY_ *)(_192_)); _192_ =dnsdata_47i+_192_;}_191_+=dnsdata_107i;_191_1--;}while(_191_1>=0);}} }iz_++;}while(iz_<=dnsdata_106h);}
        {int _193i_=dnsdata_71h;do{{if((ssize_t)sizeof(struct VELOCITY_)==(ssize_t)sizeof(struct VELOCITY_)&&(ssize_t)sizeof(struct VELOCITY_)==(ssize_t)sizeof(struct VELOCITY_)){if(!(fseeko(diskimage_,_193i_*(off_t)dnsdata_74i+iy_*(off_t)dnsdata_75i+(ssize_t)sizeof(struct dnsdata_s76)-dnsdata_81st+(size_t)(char*)(struct dnsdata_s76*)(0)+dnsdata_80st,SEEK_SET)==0&&fwrite(_193i_*dnsdata_107i+velbuf_+dnsdata_80st,(ssize_t)sizeof(struct VELOCITY_),dnsdata_73h-dnsdata_72l+1,diskimage_)==dnsdata_73h-dnsdata_72l+1))ioerr(diskimage_);}else{  {int _194i_=dnsdata_72l;do{{if(!(fseeko(diskimage_,_194i_*(off_t)(ssize_t)sizeof(struct VELOCITY_)+_193i_*(off_t)dnsdata_74i+iy_*(off_t)dnsdata_75i+(ssize_t)sizeof(struct dnsdata_s76)-dnsdata_81st+(size_t)(char*)(struct dnsdata_s76*)(0),SEEK_SET)==0&&fwrite((struct VELOCITY_*)(_194i_*(ssize_t)sizeof(struct VELOCITY_)+_193i_*dnsdata_107i+velbuf_),(ssize_t)sizeof(struct VELOCITY_),1,diskimage_)==1))ioerr(diskimage_);}_194i_++;}while(_194i_<=dnsdata_73h);}};}_193i_--;}while(_193i_>=0);}
   iy_+=1;};}

  if(diskimage_)errfclose(&diskimage_);
  if( !(next_==NULL) ){  putc('\n', next_);};
    {int iy_=(dnsdata_2nyl-2);do{{ { struct VELOCITY_ *_195w;
_195w=(struct VELOCITY_ *)(iy_*(ssize_t)sizeof(struct VELOCITY_)+V_);
{register double temp=(*_195w).u_.REAL_+u_conv_;(*_195w).u_.IMAG_=(*_195w).u_.IMAG_;(*_195w).u_.REAL_=temp;}; {register double temp=(*_195w).w_.REAL_+w_conv_;(*_195w).w_.IMAG_=(*_195w).w_.IMAG_;(*_195w).w_.REAL_=temp;};}}iy_++;}while(iy_<=(dnsdata_3nyh+2));}
}};
    
if( (time_>0.)&&( (int)floor((time_+0.5*deltat_)/dt_field_) >(int)floor((time_-0.5*deltat_)/dt_field_) )){
  char *field_name_;
cont_+=1;  field_name_=0;
   {char *tmp; tmp=malloc(snprintf(NULL,0,"""Field""%d"".fld""",cont_)+1); sprintf(tmp,"""Field""%d"".fld""",cont_); if(field_name_ )free(field_name_ );field_name_ =tmp;}
    {int iy_=(dnsdata_2nyl-2);do{{ { struct VELOCITY_ *_196w;
_196w=(struct VELOCITY_ *)(iy_*(ssize_t)sizeof(struct VELOCITY_)+V_);
{register double temp=(*_196w).u_.REAL_-u_conv_;(*_196w).u_.IMAG_=(*_196w).u_.IMAG_;(*_196w).u_.REAL_=temp;}; {register double temp=(*_196w).w_.REAL_-w_conv_;(*_196w).w_.IMAG_=(*_196w).w_.IMAG_;(*_196w).w_.REAL_=temp;};}}iy_++;}while(iy_<=(dnsdata_3nyh+2));}
  if( !(prev_==NULL) ){  while(!feof( prev_)&&getc( prev_)!='\n'){};};
  errfopen(&diskfield_ ,field_name_,O_RDWR|O_CREAT); 
  {
    int _198l;
int _199h;
int _200l;
int _201h;
ssize_t _202st;
int _204l;
int _205h;
int _206l;
int _207h;
ssize_t _208st;
if( (next_==NULL) ){
      double fwrite197_;
 fprintf(stdout,"""Writing field_file at time""\t%g" ,time_);putc('\n',stdout);
      if(!(fseeko(diskfield_,(size_t)&(*(struct dnsdata_s89*)(0)).nyimage_,SEEK_SET)==0&&fwrite(&ny_,(ssize_t)sizeof(int),1,diskfield_)==1))ioerr(diskfield_); if(!(fseeko(diskfield_,(size_t)&(*(struct dnsdata_s89*)(0)).nximage_,SEEK_SET)==0&&fwrite(&nx_,(ssize_t)sizeof(int),1,diskfield_)==1))ioerr(diskfield_); if(!(fseeko(diskfield_,(size_t)&(*(struct dnsdata_s89*)(0)).nzimage_,SEEK_SET)==0&&fwrite(&nz_,(ssize_t)sizeof(int),1,diskfield_)==1))ioerr(diskfield_);
      if(!(fseeko(diskfield_,(size_t)&(*(struct dnsdata_s89*)(0)).timage_,SEEK_SET)==0&&fwrite(&time_,(ssize_t)sizeof(double),1,diskfield_)==1))ioerr(diskfield_); if(!(fseeko(diskfield_,(size_t)&(*(struct dnsdata_s89*)(0)).yminimage_,SEEK_SET)==0&&fwrite(&ymin_,(ssize_t)sizeof(double),1,diskfield_)==1))ioerr(diskfield_); if(!(fseeko(diskfield_,(size_t)&(*(struct dnsdata_s89*)(0)).ymaximage_,SEEK_SET)==0&&fwrite(&ymax_,(ssize_t)sizeof(double),1,diskfield_)==1))ioerr(diskfield_);
      if(!(fseeko(diskfield_,(size_t)&(*(struct dnsdata_s89*)(0)).aimage_,SEEK_SET)==0&&fwrite(&a_,(ssize_t)sizeof(double),1,diskfield_)==1))ioerr(diskfield_); if(!(fseeko(diskfield_,(size_t)&(*(struct dnsdata_s89*)(0)).alfa0image_,SEEK_SET)==0&&fwrite(&alfa0_,(ssize_t)sizeof(double),1,diskfield_)==1))ioerr(diskfield_); if(!(fseeko(diskfield_,(size_t)&(*(struct dnsdata_s89*)(0)).beta0image_,SEEK_SET)==0&&fwrite(&beta0_,(ssize_t)sizeof(double),1,diskfield_)==1))ioerr(diskfield_); fwrite197_=1./ni_;if(!(fseeko(diskfield_,(size_t)&(*(struct dnsdata_s89*)(0)).niimage_,SEEK_SET)==0&&fwrite(&fwrite197_,(ssize_t)sizeof(double),1,diskfield_)==1))ioerr(diskfield_);
    };
    _198l=miny_;
_199h=maxy_;
_200l=miny_;
_201h=maxy_;
_202st=_198l*(ssize_t)sizeof(double);
if((ssize_t)sizeof(double)==(ssize_t)sizeof(double)&&(ssize_t)sizeof(struct VELOCITY_)==(ssize_t)sizeof(double)){if(!(fseeko(diskfield_,(ssize_t)sizeof(struct dnsdata_s89)-dnsdata_91st+(size_t)(char*)(struct dnsdata_s89*)(0)+_202st,SEEK_SET)==0&&fwrite((char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(V_)).u_)).REAL_+_202st,(ssize_t)sizeof(double),_199h-_198l+1,diskfield_)==_199h-_198l+1))ioerr(diskfield_);}else{  {int _203i_=_198l;do{{if(!(fseeko(diskfield_,_203i_*(off_t)(ssize_t)sizeof(double)+(ssize_t)sizeof(struct dnsdata_s89)-dnsdata_91st+(size_t)(char*)(struct dnsdata_s89*)(0),SEEK_SET)==0&&fwrite((double*)(_203i_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(V_)).u_)).REAL_),(ssize_t)sizeof(double),1,diskfield_)==1))ioerr(diskfield_);}_203i_++;}while(_203i_<=_199h);}};
    _204l=miny_;
_205h=maxy_;
_206l=miny_;
_207h=maxy_;
_208st=_204l*(ssize_t)sizeof(double);
if((ssize_t)sizeof(double)==(ssize_t)sizeof(double)&&(ssize_t)sizeof(struct VELOCITY_)==(ssize_t)sizeof(double)){if(!(fseeko(diskfield_,((ssize_t)sizeof(struct dnsdata_s89)+dnsdata_90i)-dnsdata_91st+(size_t)(char*)(struct dnsdata_s89*)(0)+_208st,SEEK_SET)==0&&fwrite((char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(V_)).w_)).REAL_+_208st,(ssize_t)sizeof(double),_205h-_204l+1,diskfield_)==_205h-_204l+1))ioerr(diskfield_);}else{  {int _209i_=_204l;do{{if(!(fseeko(diskfield_,_209i_*(off_t)(ssize_t)sizeof(double)+((ssize_t)sizeof(struct dnsdata_s89)+dnsdata_90i)-dnsdata_91st+(size_t)(char*)(struct dnsdata_s89*)(0),SEEK_SET)==0&&fwrite((double*)(_209i_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(V_)).w_)).REAL_),(ssize_t)sizeof(double),1,diskfield_)==1))ioerr(diskfield_);}_209i_++;}while(_209i_<=_205h);}};

     {int iy_=  miny_  ;while(iy_<=maxy_){
       {int ix_=dnsdata_98h;do{{int iz_=dnsdata_99l;do{ {
      struct VETA_ *_210w;
double _211r;
/*ialfa=0.,_211r*/
double _213r;
/*ibeta=0.,_213r*/
_210w=(struct VETA_ *)(iz_*(ssize_t)sizeof(struct VETA_)+ix_*dnsdata_101i+fieldbuf_);
 {register double temp=(*(struct VELOCITY_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).v_.REAL_;(*_210w).v_ .IMAG_=(*(struct VELOCITY_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).v_.IMAG_;(*_210w).v_ .REAL_=temp;};
      _211r=alfa0_*(double)(ix_);
  _213r=beta0_*(double)(iz_);

       {register double temp=(-_213r*(*(struct VELOCITY_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).u_.IMAG_)-(-_211r*(*(struct VELOCITY_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).w_.IMAG_);(*_210w).eta_ .IMAG_=(_213r*(*(struct VELOCITY_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).u_.REAL_)-(_211r*(*(struct VELOCITY_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).w_.REAL_);(*_210w).eta_ .REAL_=temp;};
     }iz_++;}while(iz_<=dnsdata_100h);}ix_--;}while(ix_>=0);}
      {int _215i_=dnsdata_84h;do{{if((ssize_t)sizeof(struct VETA_)==(ssize_t)sizeof(struct VETA_)&&(ssize_t)sizeof(struct VETA_)==(ssize_t)sizeof(struct VETA_)){if(!(fseeko(diskfield_,_215i_*(off_t)dnsdata_87i+iy_*(off_t)dnsdata_88i+(((ssize_t)sizeof(struct dnsdata_s89)+dnsdata_90i)+dnsdata_90i)-dnsdata_96st+(size_t)(char*)(struct dnsdata_s89*)(0)+dnsdata_95st,SEEK_SET)==0&&fwrite(_215i_*dnsdata_101i+fieldbuf_+dnsdata_95st,(ssize_t)sizeof(struct VETA_),dnsdata_86h-dnsdata_85l+1,diskfield_)==dnsdata_86h-dnsdata_85l+1))ioerr(diskfield_);}else{  {int _216i_=dnsdata_85l;do{{if(!(fseeko(diskfield_,_216i_*(off_t)(ssize_t)sizeof(struct VETA_)+_215i_*(off_t)dnsdata_87i+iy_*(off_t)dnsdata_88i+(((ssize_t)sizeof(struct dnsdata_s89)+dnsdata_90i)+dnsdata_90i)-dnsdata_96st+(size_t)(char*)(struct dnsdata_s89*)(0),SEEK_SET)==0&&fwrite((struct VETA_*)(_216i_*(ssize_t)sizeof(struct VETA_)+_215i_*dnsdata_101i+fieldbuf_),(ssize_t)sizeof(struct VETA_),1,diskfield_)==1))ioerr(diskfield_);}_216i_++;}while(_216i_<=dnsdata_86h);}};}_215i_--;}while(_215i_>=0);}
    iy_+=1;};}

  if(diskfield_)errfclose(&diskfield_);
  if( !(next_==NULL) ){  putc('\n', next_);};
    {int iy_=(dnsdata_2nyl-2);do{{ { struct VELOCITY_ *_217w;
_217w=(struct VELOCITY_ *)(iy_*(ssize_t)sizeof(struct VELOCITY_)+V_);
{register double temp=(*_217w).u_.REAL_+u_conv_;(*_217w).u_.IMAG_=(*_217w).u_.IMAG_;(*_217w).u_.REAL_=temp;}; {register double temp=(*_217w).w_.REAL_+w_conv_;(*_217w).w_.IMAG_=(*_217w).w_.IMAG_;(*_217w).w_.REAL_=temp;};}}iy_++;}while(iy_<=(dnsdata_3nyh+2));}
}if(field_name_)free(field_name_);};
}}
    
void dnsdata_111read_restart_file(){{
if( restart_file_[0]==0 ){
  time_=0.;  memset(dnsdata_50st+V_,0,dnsdata_48i);
  if( (next_==NULL) ){  {int iy_=(dnsdata_3nyh+2)-10  ;while(iy_<=(dnsdata_3nyh+2)){
      {int ix_=dnsdata_43h;do{{int iz_=dnsdata_44l;do{{ { struct VELOCITY_ *_112w;
double _113r;
struct COMPLEX_ _114;
double _115r;
struct COMPLEX_ _116;
double _117r;
struct COMPLEX_ _118;
_112w=(struct VELOCITY_ *)(iy_*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_);
_113r=((double)(rand())/2147483647.)*2.*3.14159265358979323846;
complex_2EXP(&_114,0.,_113r);{register double temp=(0.0001*_114.REAL_);(*_112w).u_.IMAG_=(0.0001*_114.IMAG_);(*_112w).u_.REAL_=temp;}; _115r=((double)(rand())/2147483647.)*2.*3.14159265358979323846;
complex_2EXP(&_116,0.,_115r);{register double temp=(0.0001*_116.REAL_);(*_112w).v_.IMAG_=(0.0001*_116.IMAG_);(*_112w).v_.REAL_=temp;}; _117r=((double)(rand())/2147483647.)*2.*3.14159265358979323846;
complex_2EXP(&_118,0.,_117r);{register double temp=(0.0001*_118.REAL_);(*_112w).w_.IMAG_=(0.0001*_118.IMAG_);(*_112w).w_.REAL_=temp;};}}iz_++;}while(iz_<=dnsdata_45h);}ix_--;}while(ix_>=0);}
     {int iz_=1  ;do{{ {register double temp=(*(struct VELOCITY_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+V_)).u_.REAL_;(*(struct VELOCITY_ *)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(-iz_)*dnsdata_46i+V_)).u_.IMAG_=-((*(struct VELOCITY_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+V_)).u_.IMAG_);(*(struct VELOCITY_ *)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(-iz_)*dnsdata_46i+V_)).u_.REAL_=temp;}; {register double temp=(*(struct VELOCITY_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+V_)).v_.REAL_;(*(struct VELOCITY_ *)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(-iz_)*dnsdata_46i+V_)).v_.IMAG_=-((*(struct VELOCITY_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+V_)).v_.IMAG_);(*(struct VELOCITY_ *)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(-iz_)*dnsdata_46i+V_)).v_.REAL_=temp;}; {register double temp=(*(struct VELOCITY_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+V_)).w_.REAL_;(*(struct VELOCITY_ *)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(-iz_)*dnsdata_46i+V_)).w_.IMAG_=-((*(struct VELOCITY_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+V_)).w_.IMAG_) ;(*(struct VELOCITY_ *)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(-iz_)*dnsdata_46i+V_)).w_.REAL_=temp;};}iz_+=1;}while(iz_<=nz_);}
  iy_+=1;};}};
   {int iy_=(dnsdata_2nyl-2)  ;do{{ { struct VELOCITY_ *_119w;
_119w=(struct VELOCITY_ *)(iy_*(ssize_t)sizeof(struct VELOCITY_)+V_);
(*_119w).u_.REAL_=3./2.*(1.-pow((1.-(*(double *)(iy_*(ssize_t)sizeof(double)+y_))),2)); (*_119w).u_.IMAG_=0.; memset(&(*_119w).v_,0,(ssize_t)sizeof(struct COMPLEX_)); (*_119w).w_.IMAG_=0.;}}iy_+=1;}while(iy_<=(dnsdata_3nyh+2));}
}else{
  if( (next_==NULL) ){  fprintf(stdout,"""Reading from restart_file: ""\t%s" ,restart_file_);putc('\n',stdout);};
  errfopen(&diskimage_ ,restart_file_,O_RDWR|O_CREAT);   {
    if( time_from_restart_ ){
      fseeko(diskimage_,(-(0))+(size_t)(*(struct dnsdata_s76 *)(0)).header_ ,SEEK_SET);   if(!(scanrec( diskimage_," ny = %d",&ny_)>0&&scanrec( diskimage_," nx = %d",&nx_)>0&&scanrec( diskimage_," nz = %d",&nz_)>0&&scanrec( diskimage_," alfa0 = %lg",&alfa0_)>0&&scanrec( diskimage_," beta0 = %lg",&beta0_)>0&&scanrec( diskimage_," ymin = %lg",&ymin_)>0&&scanrec( diskimage_," ymax = %lg",&ymax_)>0&&scanrec( diskimage_," a = %lg",&a_)>0&&scanrec( diskimage_," ni = %lg",&ni_)>0&&scanrec( diskimage_," time = %lg",&time_)>0))ioerr( diskimage_);  ni_=1./ni_;
      if( (next_==NULL) ){  fprintf(stdout,"""Starting at non-zero time=""%g",time_);putc('\n',stdout);};
    };
     {int iy_=  dnsdata_2nyl-2  ;while(iy_<=dnsdata_3nyh+2){
         {int _120i_=dnsdata_43h;do{{ssize_t _121st;
_121st=dnsdata_44l*(ssize_t)sizeof(struct VELOCITY_);
if(dnsdata_46i==(ssize_t)sizeof(struct VELOCITY_)&&(ssize_t)sizeof(struct VELOCITY_)==(ssize_t)sizeof(struct VELOCITY_)){if(!(fseeko(diskimage_,_120i_*(off_t)dnsdata_74i+iy_*(off_t)dnsdata_75i+(ssize_t)sizeof(struct dnsdata_s76)-dnsdata_81st+(size_t)(char*)(struct dnsdata_s76*)(0)+_121st,SEEK_SET)==0&&fread(_120i_*dnsdata_47i+iy_*(ssize_t)sizeof(struct VELOCITY_)+V_+_121st,(ssize_t)sizeof(struct VELOCITY_),dnsdata_45h-dnsdata_44l+1,diskimage_)==dnsdata_45h-dnsdata_44l+1))ioerr(diskimage_);}else{  {int _122i_=dnsdata_44l;do{{if(!(fseeko(diskimage_,_122i_*(off_t)(ssize_t)sizeof(struct VELOCITY_)+_120i_*(off_t)dnsdata_74i+iy_*(off_t)dnsdata_75i+(ssize_t)sizeof(struct dnsdata_s76)-dnsdata_81st+(size_t)(char*)(struct dnsdata_s76*)(0) ,SEEK_SET)==0&&fread((struct VELOCITY_ *)(_122i_*dnsdata_46i+_120i_*dnsdata_47i+iy_*(ssize_t)sizeof(struct VELOCITY_)+V_),(ssize_t)sizeof(struct VELOCITY_),1,diskimage_)==1))ioerr(diskimage_);}_122i_++;}while(_122i_<=dnsdata_45h);}};}_120i_--;}while(_120i_>=0);}
    iy_+=1;};}

  if(diskimage_)errfclose(&diskimage_);
}};
  {int iy_=(dnsdata_2nyl-2);do{{ { struct VELOCITY_ *_123w;
_123w=(struct VELOCITY_ *)(iy_*(ssize_t)sizeof(struct VELOCITY_)+V_);
{register double temp=(*_123w).u_.REAL_+u_conv_;(*_123w).u_.IMAG_=(*_123w).u_.IMAG_;(*_123w).u_.REAL_=temp;}; {register double temp=(*_123w).w_.REAL_+w_conv_;(*_123w).w_.IMAG_=(*_123w).w_.IMAG_;(*_123w).w_.REAL_=temp;};}}iy_++;}while(iy_<=(dnsdata_3nyh+2));}
}}

void dnsdata_112simple(struct COMPLEX_ *rhs_,int old_l,int old_h,size_t old_i,char *old__,double unkn_REAL,double unkn_IMAG,double impl_REAL,double impl_IMAG,double expl_REAL,double expl_IMAG){{
  double _113r;
_113r=(1./(deltat_));
{register double temp=(unkn_REAL*_113r)+expl_REAL;(*rhs_).IMAG_=(unkn_IMAG*_113r)+expl_IMAG;(*rhs_).REAL_=temp;};
}}

void dnsdata_113CN_AB(struct COMPLEX_ *rhs_,int old_l,int old_h,size_t old_i,char *old__,double unkn_REAL,double unkn_IMAG,double impl_REAL,double impl_IMAG,double expl_REAL,double expl_IMAG){{
  double _114r;
_114r=2./deltat_;
{register double temp=(_114r*unkn_REAL)+impl_REAL+(3.*expl_REAL)-((*(struct COMPLEX_*)(old_i+old__)).REAL_);(*rhs_).IMAG_=(_114r*unkn_IMAG)+impl_IMAG+(3.*expl_IMAG)-((*(struct COMPLEX_*)(old_i+old__)).IMAG_);(*rhs_).REAL_=temp;};
  {register double temp=expl_REAL;(*(struct COMPLEX_ *)(old_i+old__)).IMAG_=expl_IMAG;(*(struct COMPLEX_ *)(old_i+old__)).REAL_=temp;};
}}

void dnsdata_114RK1_rai(struct COMPLEX_ *rhs_,int old_l,int old_h,size_t old_i,char *old__,double unkn_REAL,double unkn_IMAG,double impl_REAL,double impl_IMAG,double expl_REAL,double expl_IMAG){{
  double _115r;
_115r=120./32./deltat_;
{register double temp=(_115r*unkn_REAL)+impl_REAL+(2.*expl_REAL);(*rhs_).IMAG_=(_115r*unkn_IMAG)+impl_IMAG+(2.*expl_IMAG);(*rhs_).REAL_=temp;};
  {register double temp=expl_REAL;(*(struct COMPLEX_ *)(old_i+old__)).IMAG_=expl_IMAG;(*(struct COMPLEX_ *)(old_i+old__)).REAL_=temp;};
}}
double dnsdata_115RK1_rai_coeff;

void dnsdata_116RK2_rai(struct COMPLEX_ *rhs_,int old_l,int old_h,size_t old_i,char *old__,double unkn_REAL,double unkn_IMAG,double impl_REAL,double impl_IMAG,double expl_REAL,double expl_IMAG){{
  double _117r;
double _118r;
double _119r;
_117r=120./(8.*deltat_);
_118r=50./8.;
_119r=34./8.;
{register double temp=(_117r*unkn_REAL)+impl_REAL+(_118r*expl_REAL)-(_119r*(*(struct COMPLEX_*)(old_i+old__)).REAL_);(*rhs_).IMAG_=(_117r*unkn_IMAG)+impl_IMAG+(_118r*expl_IMAG)-(_119r*(*(struct COMPLEX_*)(old_i+old__)).IMAG_);(*rhs_).REAL_=temp;};
  {register double temp=expl_REAL;(*(struct COMPLEX_ *)(old_i+old__)).IMAG_=expl_IMAG;(*(struct COMPLEX_ *)(old_i+old__)).REAL_=temp;};
}}
double dnsdata_117RK2_rai_coeff;

void dnsdata_118RK3_rai(struct COMPLEX_ *rhs_,int old_l,int old_h,size_t old_i,char *old__,double unkn_REAL,double unkn_IMAG,double impl_REAL,double impl_IMAG,double expl_REAL,double expl_IMAG){{
  double _119r;
double _120r;
double _121r;
_119r=120./(20.*deltat_);
_120r=90./20.;
_121r=50./20.;
{register double temp=(_119r*unkn_REAL)+impl_REAL+(_120r*expl_REAL)-(_121r*(*(struct COMPLEX_*)(old_i+old__)).REAL_);(*rhs_).IMAG_=(_119r*unkn_IMAG)+impl_IMAG+(_120r*expl_IMAG)-(_121r*(*(struct COMPLEX_*)(old_i+old__)).IMAG_);(*rhs_).REAL_=temp;};
}}
double dnsdata_119RK3_rai_coeff;

void dnsdata_120RK1_kom(struct COMPLEX_ *rhs_,int old_l,int old_h,size_t old_i,char *old__,double unkn_REAL,double unkn_IMAG,double impl_REAL,double impl_IMAG,double expl_REAL,double expl_IMAG){{
  double _121r;
_121r=1020./240./deltat_;
{register double temp=(_121r*unkn_REAL)+impl_REAL+(2.*expl_REAL);(*rhs_).IMAG_=(_121r*unkn_IMAG)+impl_IMAG+(2.*expl_IMAG);(*rhs_).REAL_=temp;};
  {register double temp=expl_REAL;(*(struct COMPLEX_ *)(old_i+old__)).IMAG_=expl_IMAG;(*(struct COMPLEX_ *)(old_i+old__)).REAL_=temp;};
}}
double dnsdata_121RK1_kom_coeff;

void dnsdata_122RK2_kom(struct COMPLEX_ *rhs_,int old_l,int old_h,size_t old_i,char *old__,double unkn_REAL,double unkn_IMAG,double impl_REAL,double impl_IMAG,double expl_REAL,double expl_IMAG){{
  double _123r;
double _124r;
double _125r;
_123r=1020./(32.*deltat_);
_124r=289./32.;
_125r=225./32.;
{register double temp=(_123r*unkn_REAL)+impl_REAL+(_124r*expl_REAL)-(_125r*(*(struct COMPLEX_*)(old_i+old__)).REAL_);(*rhs_).IMAG_=(_123r*unkn_IMAG)+impl_IMAG+(_124r*expl_IMAG)-(_125r*(*(struct COMPLEX_*)(old_i+old__)).IMAG_);(*rhs_).REAL_=temp;};
  {register double temp=expl_REAL;(*(struct COMPLEX_ *)(old_i+old__)).IMAG_=expl_IMAG;(*(struct COMPLEX_ *)(old_i+old__)).REAL_=temp;};
}}
double dnsdata_123RK2_kom_coeff;

void dnsdata_124RK3_kom(struct COMPLEX_ *rhs_,int old_l,int old_h,size_t old_i,char *old__,double unkn_REAL,double unkn_IMAG,double impl_REAL,double impl_IMAG,double expl_REAL,double expl_IMAG){{
  double _125r;
double _126r;
double _127r;
_125r=1020./(68.*deltat_);
_126r=25./4.;
_127r=17./4.;
{register double temp=(_125r*unkn_REAL)+impl_REAL+(_126r*expl_REAL)-(_127r*(*(struct COMPLEX_*)(old_i+old__)).REAL_);(*rhs_).IMAG_=(_125r*unkn_IMAG)+impl_IMAG+(_126r*expl_IMAG)-(_127r*(*(struct COMPLEX_*)(old_i+old__)).IMAG_);(*rhs_).REAL_=temp;};
  {register double temp=expl_REAL;(*(struct COMPLEX_ *)(old_i+old__)).IMAG_=expl_IMAG;(*(struct COMPLEX_ *)(old_i+old__)).REAL_=temp;};
}}
double dnsdata_125RK3_kom_coeff;

void dnsdata_126RK4_kom(struct COMPLEX_ *rhs_,int old_l,int old_h,size_t old_i,char *old__,double unkn_REAL,double unkn_IMAG,double impl_REAL,double impl_IMAG,double expl_REAL,double expl_IMAG){{
  double _127r;
double _128r;
double _129r;
_127r=1020./(170.*deltat_);
_128r=9./2.;
_129r=5./2.;
{register double temp=(_127r*unkn_REAL)+impl_REAL+(_128r*expl_REAL)-(_129r*(*(struct COMPLEX_*)(old_i+old__)).REAL_);(*rhs_).IMAG_=(_127r*unkn_IMAG)+impl_IMAG+(_128r*expl_IMAG)-(_129r*(*(struct COMPLEX_*)(old_i+old__)).IMAG_);(*rhs_).REAL_=temp;};
}}
double dnsdata_127RK4_kom_coeff;

