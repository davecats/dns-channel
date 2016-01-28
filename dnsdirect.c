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
void dnsdata_1read_initial_data();
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
void dnsdata_32applybc_0(int eq_l,int eq_h,size_t eq_i,int eq__l,int eq__h,size_t eq__i,char *eq___,char *bc_);
void dnsdata_33applybc_n(int eq_l,int eq_h,size_t eq_i,int eq__l,int eq__h,size_t eq__i,char *eq___,char *bc_);
void dnsdata_34yintegr(double *RESULT_,int f_l,int f_h,size_t f_i,char *f__);
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
void dnsdata_97getcfl();
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
void dnsdata_110outstats();
void dnsdata_111read_restart_file();
void dnsdata_112simple(struct COMPLEX_ *rhs_,int old_l,int old_h,size_t old_i,char *old__,double unkn_REAL,double unkn_IMAG,double impl_REAL,double impl_IMAG,double expl_REAL,double expl_IMAG);

void dnsdata_113CN_AB(struct COMPLEX_ *rhs_,int old_l,int old_h,size_t old_i,char *old__,double unkn_REAL,double unkn_IMAG,double impl_REAL,double impl_IMAG,double expl_REAL,double expl_IMAG);

void dnsdata_114RK1_rai(struct COMPLEX_ *rhs_,int old_l,int old_h,size_t old_i,char *old__,double unkn_REAL,double unkn_IMAG,double impl_REAL,double impl_IMAG,double expl_REAL,double expl_IMAG);
double dnsdata_115RK1_rai_coeff;

void dnsdata_116RK2_rai(struct COMPLEX_ *rhs_,int old_l,int old_h,size_t old_i,char *old__,double unkn_REAL,double unkn_IMAG,double impl_REAL,double impl_IMAG,double expl_REAL,double expl_IMAG);
double dnsdata_117RK2_rai_coeff;

void dnsdata_118RK3_rai(struct COMPLEX_ *rhs_,int old_l,int old_h,size_t old_i,char *old__,double unkn_REAL,double unkn_IMAG,double impl_REAL,double impl_IMAG,double expl_REAL,double expl_IMAG);
double dnsdata_119RK3_rai_coeff;

void dnsdata_120RK1_kom(struct COMPLEX_ *rhs_,int old_l,int old_h,size_t old_i,char *old__,double unkn_REAL,double unkn_IMAG,double impl_REAL,double impl_IMAG,double expl_REAL,double expl_IMAG);
double dnsdata_121RK1_kom_coeff;

void dnsdata_122RK2_kom(struct COMPLEX_ *rhs_,int old_l,int old_h,size_t old_i,char *old__,double unkn_REAL,double unkn_IMAG,double impl_REAL,double impl_IMAG,double expl_REAL,double expl_IMAG);
double dnsdata_123RK2_kom_coeff;

void dnsdata_124RK3_kom(struct COMPLEX_ *rhs_,int old_l,int old_h,size_t old_i,char *old__,double unkn_REAL,double unkn_IMAG,double impl_REAL,double impl_IMAG,double expl_REAL,double expl_IMAG);
double dnsdata_125RK3_kom_coeff;

void dnsdata_126RK4_kom(struct COMPLEX_ *rhs_,int old_l,int old_h,size_t old_i,char *old__,double unkn_REAL,double unkn_IMAG,double impl_REAL,double impl_IMAG,double expl_REAL,double expl_IMAG);
double dnsdata_127RK4_kom_coeff;

/*!USE rtchecks*/
/**/
int dnsdirect_1h;
int dnsdirect_2l;
int dnsdirect_3h;
struct dnsdirect_s4{struct COMPLEX_ u_;struct COMPLEX_ v_;struct COMPLEX_ w_;struct COMPLEX_ vy_;struct COMPLEX_ eta_;};size_t dnsdirect_5i;
size_t dnsdirect_6i;
ssize_t dnsdirect_7st;
char *bc0_;
char *bcn_;
void dnsdirect_8buildrhs(void (*timescheme_)(struct COMPLEX_ *,int,int,size_t,char *,double ,double ,double ,double ,double ,double ));
void dnsdirect_9linsolve(double lambda_);
void dnsdirect_10vetaTOuvw();
void dnsdirect_11calcp0(struct COMPLEX_ *RESULT_,int ix_,int iz_);
void dnsdirect_12calcpn(struct COMPLEX_ *RESULT_,int ix_,int iz_);
void dnsdirect_13computeflowrate(double lambda_);
void dnsdirect_14deriv(int f0_l,int f0_h,size_t f0_i,char *f0__,int f1_l,int f1_h,size_t f1_i,char *f1__);
void dnsdirect_15convolutions(int Vplane_l,int Vplane_h,size_t Vplane_i,int Vplane__l,int Vplane__h,size_t Vplane__i,char *Vplane___,int VVplane_l,int VVplane_h,size_t VVplane_i,int VVplane__l,int VVplane__h,size_t VVplane__i,char *VVplane___);

void dnsdirect_15convolutions(int Vplane_l,int Vplane_h,size_t Vplane_i,int Vplane__l,int Vplane__h,size_t Vplane__i,char *Vplane___,int VVplane_l,int VVplane_h,size_t VVplane_i,int VVplane__l,int VVplane__h,size_t VVplane__i,char *VVplane___){{
  /*! LOOP FOR ix=ismp TO nx BY nsmp*/
   {int ix_=(ismp_*(nx_+1) / 4 ) ;while(ix_<=((ismp_+1)*(nx_+1) / 4 )-1){
    int _16h;
int _17h;
int _20l;
int _21h;
ssize_t _22st;
size_t _23i;
int _24l;
int _25l;
_16h=nz_;
_17h=nz_;
{char *_19_;
_19_=ix_*Vplane_i+Vplane___;
 {char *_18_=ix_*dnsdata_37i+Vd_;int _18_1=_16h; do{{ (*(struct VELOCITY_ *)(_18_))=(*(struct VELOCITY_ *)(_19_)); _19_ =Vplane__i+_19_;}_18_+=(ssize_t)sizeof(struct VELOCITY_);_18_1--;}while(_18_1>=0);}}
    _20l=nz_+1;
_21h=nzd_-nz_-1;
_22st=_20l*(ssize_t)sizeof(struct VELOCITY_);
_23i=(ssize_t)sizeof(struct VELOCITY_)*(_21h-_20l+1);
memset(_22st+ix_*dnsdata_37i+Vd_,0,_23i);
    _24l= - nz_;
_25l= - nz_;
{char *_27_;
_27_=_24l*Vplane__i+ix_*Vplane_i+Vplane___;
 {char *_26_=_24l*(ssize_t)sizeof(struct VELOCITY_)+nzd_*(ssize_t)sizeof(struct VELOCITY_)+ix_*dnsdata_37i+Vd_;int _26_1=(-1)-_24l; do{{ (*(struct VELOCITY_ *)(_26_))=(*(struct VELOCITY_ *)(_27_)); _27_ =Vplane__i+_27_;}_26_+=(ssize_t)sizeof(struct VELOCITY_);_26_1--;}while(_26_1>=0);}}
    { char *_28w;
_28w=ix_*dnsdata_37i+Vd_;
{fft_5IFTU(0,dnsdata_36h,(ssize_t)sizeof(struct VELOCITY_),((char*)&(*(struct VELOCITY_ *)(_28w)).u_),0,dnsdata_36h,(ssize_t)sizeof(struct VELOCITY_),((char*)&(*(struct VELOCITY_ *)(_28w)).u_));}; {fft_5IFTU(0,dnsdata_36h,(ssize_t)sizeof(struct VELOCITY_),((char*)&(*(struct VELOCITY_ *)(_28w)).v_),0,dnsdata_36h,(ssize_t)sizeof(struct VELOCITY_),((char*)&(*(struct VELOCITY_ *)(_28w)).v_));}; {fft_5IFTU(0,dnsdata_36h,(ssize_t)sizeof(struct VELOCITY_),((char*)&(*(struct VELOCITY_ *)(_28w)).w_),0,dnsdata_36h,(ssize_t)sizeof(struct VELOCITY_),((char*)&(*(struct VELOCITY_ *)(_28w)).w_));};}
  ix_+=1;};}
  if( ismp_==0 ){ int _29l;
int _30h;
ssize_t _31st;
size_t _32i;
_29l=nx_+1;
_30h=nxd_-1;
_31st=_29l*dnsdata_37i;
_32i=dnsdata_37i*(_30h-_29l+1);
memset(_31st+Vd_,0,_32i);};
  {
  int md;
while(!((*barrier_)==(ismp_)))sigsuspend(&oldmask);
  (*barrier_)=(md=((ismp_)+1) % (4),md>=0?md:md+(4));  kill(0,SIGUSR1);
  while(!((*barrier_)<=(ismp_)))sigsuspend(&oldmask);
};
   {int iz_=(ismp_*(dnsdata_36h+1) / 4 ) ;do{{
    { char *_33w;
_33w=iz_*(ssize_t)sizeof(struct VELOCITY_)+Vd_;
{fft_7RFTU(0,dnsdata_35h,dnsdata_37i,((char*)&(*(struct VELOCITY_ *)(_33w)).u_),0,dnsdata_35h,dnsdata_37i,((char*)&(*(struct VELOCITY_ *)(_33w)).u_));}; {fft_7RFTU(0,dnsdata_35h,dnsdata_37i,((char*)&(*(struct VELOCITY_ *)(_33w)).v_),0,dnsdata_35h,dnsdata_37i,((char*)&(*(struct VELOCITY_ *)(_33w)).v_));}; {fft_7RFTU(0,dnsdata_35h,dnsdata_37i,((char*)&(*(struct VELOCITY_ *)(_33w)).w_),0,dnsdata_35h,dnsdata_37i,((char*)&(*(struct VELOCITY_ *)(_33w)).w_));};}
      {int ix_=dnsdata_35h;do{{ {
      struct VELOCITY_ *_34w;
struct MOMFLUX_ *_35w;
_34w=(struct VELOCITY_ *)(iz_*(ssize_t)sizeof(struct VELOCITY_)+ix_*dnsdata_37i+Vd_);
_35w=(struct MOMFLUX_ *)(iz_*VVplane__i+ix_*VVplane_i+VVplane___);
(*_35w).uu_.REAL_=(*_34w).u_.REAL_*(*_34w).u_.REAL_;  (*_35w).uu_.IMAG_=(*_34w).u_.IMAG_*(*_34w).u_.IMAG_;
      (*_35w).uv_.REAL_=(*_34w).u_.REAL_*(*_34w).v_.REAL_;  (*_35w).uv_.IMAG_=(*_34w).u_.IMAG_*(*_34w).v_.IMAG_;
      (*_35w).vv_.REAL_=(*_34w).v_.REAL_*(*_34w).v_.REAL_;  (*_35w).vv_.IMAG_=(*_34w).v_.IMAG_*(*_34w).v_.IMAG_;
      (*_35w).vw_.REAL_=(*_34w).v_.REAL_*(*_34w).w_.REAL_;  (*_35w).vw_.IMAG_=(*_34w).v_.IMAG_*(*_34w).w_.IMAG_;
      (*_35w).ww_.REAL_=(*_34w).w_.REAL_*(*_34w).w_.REAL_;  (*_35w).ww_.IMAG_=(*_34w).w_.IMAG_*(*_34w).w_.IMAG_;
      (*_35w).uw_.REAL_=(*_34w).u_.REAL_*(*_34w).w_.REAL_;  (*_35w).uw_.IMAG_=(*_34w).u_.IMAG_*(*_34w).w_.IMAG_;
    }}ix_--;}while(ix_>=0);}
    { char *_36w;
_36w=iz_*VVplane__i+VVplane___;
{fft_8HFTU(VVplane_l,VVplane_h,VVplane_i,((char*)&(*(struct MOMFLUX_ *)(_36w)).uu_),VVplane_l,VVplane_h,VVplane_i,((char*)&(*(struct MOMFLUX_ *)(_36w)).uu_));}; {fft_8HFTU(VVplane_l,VVplane_h,VVplane_i,((char*)&(*(struct MOMFLUX_ *)(_36w)).uv_),VVplane_l,VVplane_h,VVplane_i,((char*)&(*(struct MOMFLUX_ *)(_36w)).uv_));}; {fft_8HFTU(VVplane_l,VVplane_h,VVplane_i,((char*)&(*(struct MOMFLUX_ *)(_36w)).vv_),VVplane_l,VVplane_h,VVplane_i,((char*)&(*(struct MOMFLUX_ *)(_36w)).vv_));}; {fft_8HFTU(VVplane_l,VVplane_h,VVplane_i,((char*)&(*(struct MOMFLUX_ *)(_36w)).vw_),VVplane_l,VVplane_h,VVplane_i,((char*)&(*(struct MOMFLUX_ *)(_36w)).vw_));}; {fft_8HFTU(VVplane_l,VVplane_h,VVplane_i,((char*)&(*(struct MOMFLUX_ *)(_36w)).ww_),VVplane_l,VVplane_h,VVplane_i,((char*)&(*(struct MOMFLUX_ *)(_36w)).ww_));}; {fft_8HFTU(VVplane_l,VVplane_h,VVplane_i,((char*)&(*(struct MOMFLUX_ *)(_36w)).uw_),VVplane_l,VVplane_h,VVplane_i,((char*)&(*(struct MOMFLUX_ *)(_36w)).uw_));};}
  /*! FOR iz=ismp TO HI BY nsmp*/
  }iz_+=1;}while(iz_<=((ismp_+1)*(dnsdata_36h+1) / 4 )-1);}
  {
  int md;
while(!((*barrier_)==(ismp_)))sigsuspend(&oldmask);
  (*barrier_)=(md=((ismp_)+1) % (4),md>=0?md:md+(4));  kill(0,SIGUSR1);
  while(!((*barrier_)<=(ismp_)))sigsuspend(&oldmask);
};
   {int ix_=(ismp_*(nx_+1) / 4 ) ;do{{ { char *_37w;
_37w=ix_*VVplane_i+VVplane___;
{fft_6FFTU(VVplane__l,VVplane__h,VVplane__i,((char*)&(*(struct MOMFLUX_ *)(_37w)).uu_),VVplane__l,VVplane__h,VVplane__i,((char*)&(*(struct MOMFLUX_ *)(_37w)).uu_));}; {fft_6FFTU(VVplane__l,VVplane__h,VVplane__i,((char*)&(*(struct MOMFLUX_ *)(_37w)).uv_),VVplane__l,VVplane__h,VVplane__i,((char*)&(*(struct MOMFLUX_ *)(_37w)).uv_));}; {fft_6FFTU(VVplane__l,VVplane__h,VVplane__i,((char*)&(*(struct MOMFLUX_ *)(_37w)).vv_),VVplane__l,VVplane__h,VVplane__i,((char*)&(*(struct MOMFLUX_ *)(_37w)).vv_));}; {fft_6FFTU(VVplane__l,VVplane__h,VVplane__i,((char*)&(*(struct MOMFLUX_ *)(_37w)).vw_),VVplane__l,VVplane__h,VVplane__i,((char*)&(*(struct MOMFLUX_ *)(_37w)).vw_));}; {fft_6FFTU(VVplane__l,VVplane__h,VVplane__i,((char*)&(*(struct MOMFLUX_ *)(_37w)).ww_),VVplane__l,VVplane__h,VVplane__i,((char*)&(*(struct MOMFLUX_ *)(_37w)).ww_));}; {fft_6FFTU(VVplane__l,VVplane__h,VVplane__i,((char*)&(*(struct MOMFLUX_ *)(_37w)).uw_),VVplane__l,VVplane__h,VVplane__i,((char*)&(*(struct MOMFLUX_ *)(_37w)).uw_));};}
  /*! FOR ix=ismp TO nx BY nsmp*/
  }ix_+=1;}while(ix_<=((ismp_+1)*(nx_+1) / 4 )-1);}
}}
  
void dnsdirect_8buildrhs(void (*timescheme_)(struct COMPLEX_ *,int,int,size_t,char *,double ,double ,double ,double ,double ,double )){{
fflush(NULL); ismp_=0  ;while(ismp_<4-1&&fork()){ismp_+=1;;};{

  {int iy_=dnsdata_2nyl-4  ;while(iy_<=dnsdata_3nyh+2){
  if( iy_<=dnsdata_3nyh ){
   int md;
dnsdirect_15convolutions(0,dnsdata_43h,dnsdata_47i,dnsdata_44l,dnsdata_45h,dnsdata_46i,(iy_+2)*(ssize_t)sizeof(struct VELOCITY_)+V_,0,dnsdata_39h,dnsdata_41i,0,dnsdata_40h,(ssize_t)sizeof(struct MOMFLUX_)*(4+1),(md=((iy_+2)+1000)%5,md>=0?md:md+5)*(ssize_t)sizeof(struct MOMFLUX_)+VVd_);
   if( iy_>=dnsdata_2nyl ){
    /*! WITH derivatives(iy) LOOP FOR ix=ismp TO nx BY nsmp AND iz=-nz TO nz*/
    {struct dnsdata_s10 *_16w;
_16w=(struct dnsdata_s10 *)(iy_*(ssize_t)sizeof(struct dnsdata_s10)+derivatives_) ;
 {int ix_=(ismp_*(nx_+1) / 4 ) ;while(ix_<=((ismp_+1)*(nx_+1) / 4 )-1 ){ {int iz_= - nz_  ;while(iz_<=nz_){
      double _17alfa;
/*alfa=_17alfa*/
double _18beta;
/*beta=_18beta*/
double _19k2;
/*k2=_19k2*/
struct MOMFLUX_ *VVm2_;
struct MOMFLUX_ *VVm1_;
struct MOMFLUX_ *VV0_;
struct MOMFLUX_ *VV1_;
struct MOMFLUX_ *VV2_;           
_17alfa=alfa0_*(double)(ix_);
  _18beta=beta0_*(double)(iz_);

        
      _19k2=(_17alfa*_17alfa)+(_18beta*_18beta);

      
      if( iz_>=0 ){     
        VVm2_=(struct MOMFLUX_ *)((md=((iy_-2)+1000)%5,md>=0?md:md+5)*(ssize_t)sizeof(struct MOMFLUX_)+iz_*(ssize_t)sizeof(struct MOMFLUX_)*(4+1)+ix_*dnsdata_41i+VVd_); VVm1_=(struct MOMFLUX_ *)((md=((iy_-1)+1000)%5,md>=0?md:md+5)*(ssize_t)sizeof(struct MOMFLUX_)+iz_*(ssize_t)sizeof(struct MOMFLUX_)*(4+1)+ix_*dnsdata_41i+VVd_); VV0_=(struct MOMFLUX_ *)((md=(iy_+1000)%5,md>=0?md:md+5)*(ssize_t)sizeof(struct MOMFLUX_)+iz_*(ssize_t)sizeof(struct MOMFLUX_)*(4+1)+ix_*dnsdata_41i+VVd_); VV1_=(struct MOMFLUX_ *)((md=((iy_+1)+1000)%5,md>=0?md:md+5)*(ssize_t)sizeof(struct MOMFLUX_)+iz_*(ssize_t)sizeof(struct MOMFLUX_)*(4+1)+ix_*dnsdata_41i+VVd_); VV2_=(struct MOMFLUX_ *)((md=((iy_+2)+1000)%5,md>=0?md:md+5)*(ssize_t)sizeof(struct MOMFLUX_)+iz_*(ssize_t)sizeof(struct MOMFLUX_)*(4+1)+ix_*dnsdata_41i+VVd_);
      }else{        
        VVm2_=(struct MOMFLUX_ *)((md=((iy_-2)+1000)%5,md>=0?md:md+5)*(ssize_t)sizeof(struct MOMFLUX_)+(nzd_+iz_)*(ssize_t)sizeof(struct MOMFLUX_)*(4+1)+ix_*dnsdata_41i+VVd_); VVm1_=(struct MOMFLUX_ *)((md=((iy_-1)+1000)%5,md>=0?md:md+5)*(ssize_t)sizeof(struct MOMFLUX_)+(nzd_+iz_)*(ssize_t)sizeof(struct MOMFLUX_)*(4+1)+ix_*dnsdata_41i+VVd_); VV0_=(struct MOMFLUX_ *)((md=(iy_+1000)%5,md>=0?md:md+5)*(ssize_t)sizeof(struct MOMFLUX_)+(nzd_+iz_)*(ssize_t)sizeof(struct MOMFLUX_)*(4+1)+ix_*dnsdata_41i+VVd_); VV1_=(struct MOMFLUX_ *)((md=((iy_+1)+1000)%5,md>=0?md:md+5)*(ssize_t)sizeof(struct MOMFLUX_)+(nzd_+iz_)*(ssize_t)sizeof(struct MOMFLUX_)*(4+1)+ix_*dnsdata_41i+VVd_); VV2_=(struct MOMFLUX_ *)((md=((iy_+2)+1000)%5,md>=0?md:md+5)*(ssize_t)sizeof(struct MOMFLUX_)+(nzd_+iz_)*(ssize_t)sizeof(struct MOMFLUX_)*(4+1)+ix_*dnsdata_41i+VVd_);
      };	

      {
      char *_20w;
double _21r;
double _22r;
double _23r;
double _24r;
double _25r;
double _27r_26D0_uw_;
double _28i_26D0_uw_;
/*D0_uw_=_27r_26D0_uw_,_28i_26D0_uw_*/
double _29r;
double _30r;
double _31r;
double _32r;
double _33r;
double _35r_34D1_uw_;
double _36i_34D1_uw_;
/*D1_uw_=_35r_34D1_uw_,_36i_34D1_uw_*/
double _37r;
double _38r;
double _39r;
double _40r;
double _41r;
double _42r;
double _43i;
double _44r;
double _45r;
double _46r;
double _47r;
double _48r;
double _50r_49rhsu;
double _51i_49rhsu;
/*rhsu=_50r_49rhsu,_51i_49rhsu*/
double _52r;
double _53r;
double _54r;
double _55r;
double _56r;
double _57r;
double _58i;
double _59r;
double _60r;
double _61r;
double _62r;
double _63r;
double _64r;
double _65r;
double _66r;
double _67r;
double _68r;
double _69r;
double _70i;
double _72r_71rhsv;
double _73i_71rhsv;
/*rhsv=_72r_71rhsv,_73i_71rhsv*/
double _74r;
double _75r;
double _76r;
double _77r;
double _78r;
double _79r;
double _80r;
double _81r;
double _82r;
double _83r;
double _84r;
double _85i;
double _87r_86rhsw;
double _88i_86rhsw;
/*rhsw=_87r_86rhsw,_88i_86rhsw*/
double _89r;
double _90r;
double _91r;
double _93r;
double _94r;
double _95r;
double _96r;
double _97r;
double _98r;
double _99i;
double _100r;
double _101r;
double _102r;
double _103r;
double _104r;
double _105r;
double _106i;
double _107r;
double _108r;
double _109r;
double _110r;
double _111r;
double _112r;
double _113r;
double _114r;
double _115r;
double _116r;
double _117r;
double _118i;
double _119r;
double _120i;
struct COMPLEX_ _121;
/*!       	 0.5/ni*[ialfa*[ialfa*D(d1,uu)+D(d2,uv)+ibeta*D1_uw_]+ !Vittori*/
/*!       	 ibeta*[ialfa*D1_uw_+D(d2,vw)+ibeta*D(d1,ww)]]-k2*rhsv}!Vittori*/
      _20w=iy_*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_;
_21r=(*(double *)((-2)*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st));
_22r=(*(double *)(-(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st));
_23r=(*(double *)((*_16w).d0_-dnsdata_11st));
_24r=(*(double *)((ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st));
_25r=(*(double *)(2*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st));
_27r_26D0_uw_=((_21r*(*VVm2_).uw_.REAL_)+(_22r*(*VVm1_).uw_.REAL_)+(_23r*(*VV0_).uw_.REAL_)+(_24r*(*VV1_).uw_.REAL_)+(_25r*(*VV2_).uw_.REAL_));
_28i_26D0_uw_=((_21r*(*VVm2_).uw_.IMAG_)+(_22r*(*VVm1_).uw_.IMAG_)+(_23r*(*VV0_).uw_.IMAG_)+(_24r*(*VV1_).uw_.IMAG_)+(_25r*(*VV2_).uw_.IMAG_));

      _29r=(*(double *)((-2)*(ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_30r=(*(double *)(-(ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_31r=(*(double *)((*_16w).d1_-dnsdata_11st));
_32r=(*(double *)((ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_33r=(*(double *)(2*(ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_35r_34D1_uw_=((_29r*(*VVm2_).uw_.REAL_)+(_30r*(*VVm1_).uw_.REAL_)+(_31r*(*VV0_).uw_.REAL_)+(_32r*(*VV1_).uw_.REAL_)+(_33r*(*VV2_).uw_.REAL_));
_36i_34D1_uw_=((_29r*(*VVm2_).uw_.IMAG_)+(_30r*(*VVm1_).uw_.IMAG_)+(_31r*(*VV0_).uw_.IMAG_)+(_32r*(*VV1_).uw_.IMAG_)+(_33r*(*VV2_).uw_.IMAG_));

      _37r=(*(double *)((-2)*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st));
_38r=(*(double *)(-(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st));
_39r=(*(double *)((*_16w).d0_-dnsdata_11st));
_40r=(*(double *)((ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st));
_41r=(*(double *)(2*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st));
_42r=((_37r*(*VVm2_).uu_.REAL_)+(_38r*(*VVm1_).uu_.REAL_)+(_39r*(*VV0_).uu_.REAL_)+(_40r*(*VV1_).uu_.REAL_)+(_41r*(*VV2_).uu_.REAL_));
_43i=((_37r*(*VVm2_).uu_.IMAG_)+(_38r*(*VVm1_).uu_.IMAG_)+(_39r*(*VV0_).uu_.IMAG_)+(_40r*(*VV1_).uu_.IMAG_)+(_41r*(*VV2_).uu_.IMAG_));
_44r=(*(double *)((-2)*(ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_45r=(*(double *)(-(ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_46r=(*(double *)((*_16w).d1_-dnsdata_11st));
_47r=(*(double *)((ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_48r=(*(double *)(2*(ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_50r_49rhsu=(_17alfa*_43i)-((_44r*(*VVm2_).uv_.REAL_)+(_45r*(*VVm1_).uv_.REAL_)+(_46r*(*VV0_).uv_.REAL_)+(_47r*(*VV1_).uv_.REAL_)+(_48r*(*VV2_).uv_.REAL_))-(-_18beta*_28i_26D0_uw_);
_51i_49rhsu=(-_17alfa*_42r)-((_44r*(*VVm2_).uv_.IMAG_)+(_45r*(*VVm1_).uv_.IMAG_)+(_46r*(*VV0_).uv_.IMAG_)+(_47r*(*VV1_).uv_.IMAG_)+(_48r*(*VV2_).uv_.IMAG_))-(_18beta*_27r_26D0_uw_);

      _52r=(*(double *)((-2)*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st));
_53r=(*(double *)(-(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st));
_54r=(*(double *)((*_16w).d0_-dnsdata_11st));
_55r=(*(double *)((ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st));
_56r=(*(double *)(2*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st));
_57r=((_52r*(*VVm2_).uv_.REAL_)+(_53r*(*VVm1_).uv_.REAL_)+(_54r*(*VV0_).uv_.REAL_)+(_55r*(*VV1_).uv_.REAL_)+(_56r*(*VV2_).uv_.REAL_));
_58i=((_52r*(*VVm2_).uv_.IMAG_)+(_53r*(*VVm1_).uv_.IMAG_)+(_54r*(*VV0_).uv_.IMAG_)+(_55r*(*VV1_).uv_.IMAG_)+(_56r*(*VV2_).uv_.IMAG_));
_59r=(*(double *)((-2)*(ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_60r=(*(double *)(-(ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_61r=(*(double *)((*_16w).d1_-dnsdata_11st));
_62r=(*(double *)((ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_63r=(*(double *)(2*(ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_64r=(*(double *)((-2)*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st));
_65r=(*(double *)(-(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st));
_66r=(*(double *)((*_16w).d0_-dnsdata_11st));
_67r=(*(double *)((ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st));
_68r=(*(double *)(2*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st));
_69r=((_64r*(*VVm2_).vw_.REAL_)+(_65r*(*VVm1_).vw_.REAL_)+(_66r*(*VV0_).vw_.REAL_)+(_67r*(*VV1_).vw_.REAL_)+(_68r*(*VV2_).vw_.REAL_));
_70i=((_64r*(*VVm2_).vw_.IMAG_)+(_65r*(*VVm1_).vw_.IMAG_)+(_66r*(*VV0_).vw_.IMAG_)+(_67r*(*VV1_).vw_.IMAG_)+(_68r*(*VV2_).vw_.IMAG_));
_72r_71rhsv=(_17alfa*_58i)-((_59r*(*VVm2_).vv_.REAL_)+(_60r*(*VVm1_).vv_.REAL_)+(_61r*(*VV0_).vv_.REAL_)+(_62r*(*VV1_).vv_.REAL_)+(_63r*(*VV2_).vv_.REAL_))-(-_18beta*_70i);
_73i_71rhsv=(-_17alfa*_57r)-((_59r*(*VVm2_).vv_.IMAG_)+(_60r*(*VVm1_).vv_.IMAG_)+(_61r*(*VV0_).vv_.IMAG_)+(_62r*(*VV1_).vv_.IMAG_)+(_63r*(*VV2_).vv_.IMAG_))-(_18beta*_69r);

      _74r=(*(double *)((-2)*(ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_75r=(*(double *)(-(ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_76r=(*(double *)((*_16w).d1_-dnsdata_11st));
_77r=(*(double *)((ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_78r=(*(double *)(2*(ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_79r=(*(double *)((-2)*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st));
_80r=(*(double *)(-(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st));
_81r=(*(double *)((*_16w).d0_-dnsdata_11st));
_82r=(*(double *)((ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st));
_83r=(*(double *)(2*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st));
_84r=((_79r*(*VVm2_).ww_.REAL_)+(_80r*(*VVm1_).ww_.REAL_)+(_81r*(*VV0_).ww_.REAL_)+(_82r*(*VV1_).ww_.REAL_)+(_83r*(*VV2_).ww_.REAL_));
_85i=((_79r*(*VVm2_).ww_.IMAG_)+(_80r*(*VVm1_).ww_.IMAG_)+(_81r*(*VV0_).ww_.IMAG_)+(_82r*(*VV1_).ww_.IMAG_)+(_83r*(*VV2_).ww_.IMAG_));
_87r_86rhsw=(_17alfa*_28i_26D0_uw_)-((_74r*(*VVm2_).vw_.REAL_)+(_75r*(*VVm1_).vw_.REAL_)+(_76r*(*VV0_).vw_.REAL_)+(_77r*(*VV1_).vw_.REAL_)+(_78r*(*VV2_).vw_.REAL_))-(-_18beta*_85i);
_88i_86rhsw=(-_17alfa*_27r_26D0_uw_)-((_74r*(*VVm2_).vw_.IMAG_)+(_75r*(*VVm1_).vw_.IMAG_)+(_76r*(*VV0_).vw_.IMAG_)+(_77r*(*VV1_).vw_.IMAG_)+(_78r*(*VV2_).vw_.IMAG_))-(_18beta*_84r);

/*!     rhsu=0.5/ni*[-ialfa*D(d0,uu)-D(d1,uv)-ibeta*D0_uw_]   !Vittori*/
/*!     rhsv=0.5/ni*[-ialfa*D(d0,uv)-D(d1,vv)-ibeta*D(d0,vw)] !Vittori*/
/*!     rhsw=0.5/ni*[-ialfa*D0_uw_-D(d1,vw)-ibeta*D(d0,ww)]   !Vittori*/
_89r=( (*(double *)((-2)*(ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st))*(*(double *)((-2)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).v_))).IMAG_)))+(*(double *)(-(ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st))*(*(double *)(-(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).v_))).IMAG_)))+(*(double *)((*_16w).d2_-dnsdata_11st))*(*(double *)(((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).v_))).IMAG_)))+(*(double *)((ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st))*(*(double *)((ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).v_))).IMAG_)))+(*(double *)(2*(ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st))*(*(double *)(2*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).v_))).IMAG_))));
_90r=( (*(double *)((-2)*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double *)((-2)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).v_))).IMAG_)))+(*(double *)(-(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double *)(-(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).v_))).IMAG_)))+(*(double *)((*_16w).d0_-dnsdata_11st))*(*(double *)(((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).v_))).IMAG_)))+(*(double *)((ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double *)((ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).v_))).IMAG_)))+(*(double *)(2*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double *)(2*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).v_))).IMAG_))));
_91r=( ((*(double*)((-2)*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double*)((-2)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).v_))).REAL_)))+(*(double*)(-(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double*)(-(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).v_))).REAL_)))+(*(double*)((*_16w).d0_-dnsdata_11st))*(*(double*)(((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).v_))).REAL_)))+(*(double*)((ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double*)((ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).v_))).REAL_)))+(*(double*)(2*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double*)(2*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).v_))).REAL_)))));
_93r=(*(double *)((-2)*(ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_94r=(*(double *)(-(ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_95r=(*(double *)((*_16w).d1_-dnsdata_11st));
_96r=(*(double *)((ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_97r=(*(double *)(2*(ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_98r=((_93r*(*VVm2_).uu_.REAL_)+(_94r*(*VVm1_).uu_.REAL_)+(_95r*(*VV0_).uu_.REAL_)+(_96r*(*VV1_).uu_.REAL_)+(_97r*(*VV2_).uu_.REAL_));
_99i=((_93r*(*VVm2_).uu_.IMAG_)+(_94r*(*VVm1_).uu_.IMAG_)+(_95r*(*VV0_).uu_.IMAG_)+(_96r*(*VV1_).uu_.IMAG_)+(_97r*(*VV2_).uu_.IMAG_));
_100r=(*(double *)((-2)*(ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st));
_101r=(*(double *)(-(ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st));
_102r=(*(double *)((*_16w).d2_-dnsdata_11st));
_103r=(*(double *)((ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st));
_104r=(*(double *)(2*(ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st));
_105r=((-_17alfa*_99i)+((_100r*(*VVm2_).uv_.REAL_)+(_101r*(*VVm1_).uv_.REAL_)+(_102r*(*VV0_).uv_.REAL_)+(_103r*(*VV1_).uv_.REAL_)+(_104r*(*VV2_).uv_.REAL_))+(-_18beta*_36i_34D1_uw_));
_106i=((_17alfa*_98r)+((_100r*(*VVm2_).uv_.IMAG_)+(_101r*(*VVm1_).uv_.IMAG_)+(_102r*(*VV0_).uv_.IMAG_)+(_103r*(*VV1_).uv_.IMAG_)+(_104r*(*VV2_).uv_.IMAG_))+(_18beta*_35r_34D1_uw_));
_107r=(*(double *)((-2)*(ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st));
_108r=(*(double *)(-(ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st));
_109r=(*(double *)((*_16w).d2_-dnsdata_11st));
_110r=(*(double *)((ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st));
_111r=(*(double *)(2*(ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st));
_112r=(*(double *)((-2)*(ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_113r=(*(double *)(-(ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_114r=(*(double *)((*_16w).d1_-dnsdata_11st));
_115r=(*(double *)((ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_116r=(*(double *)(2*(ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st));
_117r=((_112r*(*VVm2_).ww_.REAL_)+(_113r*(*VVm1_).ww_.REAL_)+(_114r*(*VV0_).ww_.REAL_)+(_115r*(*VV1_).ww_.REAL_)+(_116r*(*VV2_).ww_.REAL_));
_118i=((_112r*(*VVm2_).ww_.IMAG_)+(_113r*(*VVm1_).ww_.IMAG_)+(_114r*(*VV0_).ww_.IMAG_)+(_115r*(*VV1_).ww_.IMAG_)+(_116r*(*VV2_).ww_.IMAG_));
_119r=((-_17alfa*_36i_34D1_uw_)+((_107r*(*VVm2_).vw_.REAL_)+(_108r*(*VVm1_).vw_.REAL_)+(_109r*(*VV0_).vw_.REAL_)+(_110r*(*VV1_).vw_.REAL_)+(_111r*(*VV2_).vw_.REAL_))+(-_18beta*_118i));
_120i=((_17alfa*_35r_34D1_uw_)+((_107r*(*VVm2_).vw_.IMAG_)+(_108r*(*VVm1_).vw_.IMAG_)+(_109r*(*VV0_).vw_.IMAG_)+(_110r*(*VV1_).vw_.IMAG_)+(_111r*(*VV2_).vw_.IMAG_))+(_18beta*_117r));
 memset(&_121,0,(ssize_t)sizeof(struct COMPLEX_));  {int i_=(-2);do{{double _92r;
_92r=(ni_*((*(double *)(i_*(ssize_t)sizeof(double)+(*_16w).d4_-dnsdata_11st))-2.*_19k2*(*(double *)(i_*(ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st))+_19k2*_19k2*(*(double *)(i_*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))));
{register double temp=(*&_121).REAL_+(_92r*(*(struct COMPLEX_*)(i_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_20w)).v_)).REAL_);(*&_121).IMAG_=(*&_121).IMAG_+(_92r*(*(struct COMPLEX_*)(i_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_20w)).v_)).IMAG_);(*&_121).REAL_=temp;};}i_++;}while(i_<=2);}(*timescheme_)(&(*((struct rhstype_ *)((md=(iy_+1000)%3,md>=0?md:md+3)*(ssize_t)sizeof(struct rhstype_)+iz_*dnsdata_66i+ix_*dnsdata_67i+memrhs_))).D2v_,1,1,(ssize_t)sizeof(struct rhstype_),(char*)&(*(struct rhstype_ *)(iy_*(ssize_t)sizeof(struct rhstype_)+iz_*dnsdata_58i+ix_*dnsdata_59i+oldrhs_)).D2v_,((*(double*)((-2)*(ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st))*(*(double*)((-2)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).v_))).REAL_)))+(*(double*)(-(ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st))*(*(double*)(-(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).v_))).REAL_)))+(*(double*)((*_16w).d2_-dnsdata_11st))*(*(double*)(((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).v_))).REAL_)))+(*(double*)((ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st))*(*(double*)((ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).v_))).REAL_)))+(*(double*)(2*(ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st))*(*(double*)(2*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).v_))).REAL_))))-(_19k2*_91r),_89r-(_19k2*_90r)
      		 ,_121.REAL_,_121.IMAG_
        	 ,(-_17alfa*_106i)+(-_18beta*_120i)-(_19k2*_72r_71rhsv),(_17alfa*_105r)+(_18beta*_119r)-(_19k2*_73i_71rhsv));
if( (ix_==0 )&&( iz_==0 )){
       /*! u media conservata in eta.REAL e w media in eta.IMAG*/
double _122r;
double _123r;
double _124r;
_122r=( (*(double *)((-2)*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double *)((-2)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(_20w)).w_)).REAL_)))+(*(double *)(-(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double *)(-(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(_20w)).w_)).REAL_)))+(*(double *)((*_16w).d0_-dnsdata_11st))*(*(double *)(((char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(_20w)).w_)).REAL_)))+(*(double *)((ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double *)((ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(_20w)).w_)).REAL_)))+(*(double *)(2*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double *)(2*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(_20w)).w_)).REAL_))));
_123r=ni_*( (*(double *)((-2)*(ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st))*(*(double *)((-2)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(_20w)).w_)).REAL_)))+(*(double *)(-(ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st))*(*(double *)(-(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(_20w)).w_)).REAL_)))+(*(double *)((*_16w).d2_-dnsdata_11st))*(*(double *)(((char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(_20w)).w_)).REAL_)))+(*(double *)((ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st))*(*(double *)((ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(_20w)).w_)).REAL_)))+(*(double *)(2*(ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st))*(*(double *)(2*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(_20w)).w_)).REAL_))));
_124r=((_87r_86rhsw)+pz_);
(*timescheme_)(&(*((struct rhstype_ *)((md=(iy_+1000)%3,md>=0?md:md+3)*(ssize_t)sizeof(struct rhstype_)+memrhs_))).eta_,1,1,(ssize_t)sizeof(struct rhstype_),(char*)&(*(struct rhstype_ *)(iy_*(ssize_t)sizeof(struct rhstype_)+oldrhs_)).eta_,((*(double*)((-2)*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double*)((-2)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(_20w)).u_)).REAL_)))+(*(double*)(-(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double*)(-(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(_20w)).u_)).REAL_)))+(*(double*)((*_16w).d0_-dnsdata_11st))*(*(double*)(((char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(_20w)).u_)).REAL_)))+(*(double*)((ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double*)((ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(_20w)).u_)).REAL_)))+(*(double*)(2*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double*)(2*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(_20w)).u_)).REAL_)))),_122r
                  ,ni_*((*(double*)((-2)*(ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st))*(*(double*)((-2)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(_20w)).u_)).REAL_)))+(*(double*)(-(ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st))*(*(double*)(-(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(_20w)).u_)).REAL_)))+(*(double*)((*_16w).d2_-dnsdata_11st))*(*(double*)(((char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(_20w)).u_)).REAL_)))+(*(double*)((ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st))*(*(double*)((ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(_20w)).u_)).REAL_)))+(*(double*)(2*(ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st))*(*(double*)(2*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(_20w)).u_)).REAL_)))),_123r
/*!                 0.5*[D2(u.REAL)]+0.5*D2(w.REAL)*I, !Vittori*/
                  ,_50r_49rhsu+px_,_124r);
}else{
double _125r;
double _126r;
double _127r;
double _128r;
struct COMPLEX_ _132;
_125r=( (*(double *)((-2)*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double *)((-2)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).u_))).IMAG_)))+(*(double *)(-(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double *)(-(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).u_))).IMAG_)))+(*(double *)((*_16w).d0_-dnsdata_11st))*(*(double *)(((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).u_))).IMAG_)))+(*(double *)((ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double *)((ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).u_))).IMAG_)))+(*(double *)(2*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double *)(2*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).u_))).IMAG_))));
_126r=( ((*(double*)((-2)*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double*)((-2)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).u_))).REAL_)))+(*(double*)(-(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double*)(-(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).u_))).REAL_)))+(*(double*)((*_16w).d0_-dnsdata_11st))*(*(double*)(((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).u_))).REAL_)))+(*(double*)((ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double*)((ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).u_))).REAL_)))+(*(double*)(2*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double*)(2*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).u_))).REAL_)))));
_127r=( (*(double *)((-2)*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double *)((-2)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).w_))).IMAG_)))+(*(double *)(-(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double *)(-(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).w_))).IMAG_)))+(*(double *)((*_16w).d0_-dnsdata_11st))*(*(double *)(((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).w_))).IMAG_)))+(*(double *)((ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double *)((ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).w_))).IMAG_)))+(*(double *)(2*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double *)(2*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).w_))).IMAG_))));
_128r=( ((*(double*)((-2)*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double*)((-2)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).w_))).REAL_)))+(*(double*)(-(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double*)(-(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).w_))).REAL_)))+(*(double*)((*_16w).d0_-dnsdata_11st))*(*(double*)(((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).w_))).REAL_)))+(*(double*)((ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double*)((ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).w_))).REAL_)))+(*(double*)(2*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))*(*(double*)(2*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_20w)).w_))).REAL_)))));
 memset(&_132,0,(ssize_t)sizeof(struct COMPLEX_));  {int i_=(-2);do{{double _129r;
double _130r;
double _131i;
_129r=(ni_*((*(double *)(i_*(ssize_t)sizeof(double)+(*_16w).d2_-dnsdata_11st))-_19k2*(*(double *)(i_*(ssize_t)sizeof(double)+(*_16w).d0_-dnsdata_11st))));
_130r=((-_18beta*(*(struct COMPLEX_*)(i_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_20w)).u_)).IMAG_)-(-_17alfa*(*(struct COMPLEX_*)(i_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_20w)).w_)).IMAG_));
_131i=((_18beta*(*(struct COMPLEX_*)(i_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_20w)).u_)).REAL_)-(_17alfa*(*(struct COMPLEX_*)(i_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_20w)).w_)).REAL_)) ;
{register double temp=(*&_132).REAL_+(_129r*_130r);(*&_132).IMAG_=(*&_132).IMAG_+(_129r*_131i);(*&_132).REAL_=temp;};}i_++;}while(i_<=2);}(*timescheme_)(&(*((struct rhstype_ *)((md=(iy_+1000)%3,md>=0?md:md+3)*(ssize_t)sizeof(struct rhstype_)+iz_*dnsdata_66i+ix_*dnsdata_67i+memrhs_))).eta_,1,1,(ssize_t)sizeof(struct rhstype_),(char*)&(*(struct rhstype_ *)(iy_*(ssize_t)sizeof(struct rhstype_)+iz_*dnsdata_58i+ix_*dnsdata_59i+oldrhs_)).eta_,(-_18beta*_125r)-(-_17alfa*_127r),(_18beta*_126r)-(_17alfa*_128r)
        	   ,_132.REAL_,_132.IMAG_
                   ,(-_18beta*_51i_49rhsu)-(-_17alfa*_88i_86rhsw),(_18beta*_50r_49rhsu)-(_17alfa*_87r_86rhsw));
};
    }iz_+=1;};}ix_+=1;};}}
   };
  };
  if( iy_-2>=dnsdata_2nyl ){
     {int ix_=(ismp_*(nx_+1) / 4 ) ;do{  {int iz_=dnsdata_44l;do{{ int md;
{ struct VELOCITY_ *_133w;
struct rhstype_ *_134w;
_133w=(struct VELOCITY_ *)((iy_-2)*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_);
_134w=((struct rhstype_ *)((md=((iy_-2)+1000)%3,md>=0?md:md+3)*(ssize_t)sizeof(struct rhstype_)+iz_*dnsdata_66i+ix_*dnsdata_67i+memrhs_));
{register double temp=(*_134w).eta_.REAL_;(*_133w).u_.IMAG_=(*_134w).eta_.IMAG_;(*_133w).u_.REAL_=temp;}; {register double temp=(*_134w).D2v_.REAL_;(*_133w).v_.IMAG_=(*_134w).D2v_.IMAG_;(*_133w).v_.REAL_=temp;};}
    /*! FOR ix=ismp TO nx BY nsmp AND ALL iz*/
    }iz_++;}while(iz_<=dnsdata_45h);}ix_+=1;}while(ix_<=((ismp_+1)*(nx_+1) / 4 )-1 );}
  };
 iy_+=1;};}
} if(ismp_<4-1)exit(0);;
 ismp_=0  ;while(ismp_<4-1){if(wait(NULL)<0){strcpy(ERRORMESSAGE_,"wait"); longjmp(errorhook,1);};ismp_+=1;;};
}}

void dnsdirect_11calcp0(struct COMPLEX_ *RESULT_,int ix_,int iz_){{
  double _16alfa;
/*alfa=_16alfa*/
double _17beta;
/*beta=_17beta*/
double _18k2;
/*k2=_18k2*/
_16alfa=alfa0_*(double)(ix_);
  _17beta=beta0_*(double)(iz_);

    
  _18k2=(_16alfa*_16alfa)+(_17beta*_17beta);

  {
  char *_19w;
_19w=iz_*dnsdata_46i+ix_*dnsdata_47i+V_;
if( _18k2>0.){
    double _20r;
double _21r;
double _22i;
double _23r;
double _24r;
double _25i;
double _26r;
double _27r;
double _28i;
double _29r;
double _30r;
double _31i;
double _32r;
double _33r;
double _34i;
double _35r;
double _36i;
double _37r;
double _38r;
double _39i;
double _40r;
double _41r;
double _42i;
double _43r;
double _44r;
double _45i;
double _46r;
double _47r;
double _48i;
double _49r;
double _50r;
double _51i;
double _52r;
double _53i;
double _54r;
double _55r;
double _56i;
_20r=(*(double *)((-2)*(ssize_t)sizeof(double)+d240_-((-2)*(ssize_t)sizeof(double))));
_21r=(*(struct COMPLEX_*)(-(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).u_))).REAL_;
_22i=(*(struct COMPLEX_*)(-(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).u_))).IMAG_;
_23r=(*(double *)(-(ssize_t)sizeof(double)+d240_-((-2)*(ssize_t)sizeof(double))));
_24r=(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_19w)).u_))).REAL_;
_25i=(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_19w)).u_))).IMAG_;
_26r=(*(double *)(d240_-((-2)*(ssize_t)sizeof(double))));
_27r=(*(struct COMPLEX_*)((ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).u_))).REAL_;
_28i=(*(struct COMPLEX_*)((ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).u_))).IMAG_;
_29r=(*(double *)((ssize_t)sizeof(double)+d240_-((-2)*(ssize_t)sizeof(double))));
_30r=(*(struct COMPLEX_*)(2*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).u_))).REAL_;
_31i=(*(struct COMPLEX_*)(2*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).u_))).IMAG_;
_32r=(*(double *)(2*(ssize_t)sizeof(double)+d240_-((-2)*(ssize_t)sizeof(double))));
_33r=(*(struct COMPLEX_*)(3*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).u_))).REAL_;
_34i=(*(struct COMPLEX_*)(3*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).u_))).IMAG_;
_35r=((_20r*_21r)+(_23r*_24r)+(_26r*_27r)+(_29r*_30r)+(_32r*_33r));
_36i=((_20r*_22i)+(_23r*_25i)+(_26r*_28i)+(_29r*_31i)+(_32r*_34i));
_37r=(*(double *)((-2)*(ssize_t)sizeof(double)+d240_-((-2)*(ssize_t)sizeof(double))));
_38r=(*(struct COMPLEX_*)(-(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).w_))).REAL_;
_39i=(*(struct COMPLEX_*)(-(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).w_))).IMAG_;
_40r=(*(double *)(-(ssize_t)sizeof(double)+d240_-((-2)*(ssize_t)sizeof(double))));
_41r=(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_19w)).w_))).REAL_;
_42i=(*(struct COMPLEX_*)(((char*)&(*(struct VELOCITY_*)(_19w)).w_))).IMAG_;
_43r=(*(double *)(d240_-((-2)*(ssize_t)sizeof(double))));
_44r=(*(struct COMPLEX_*)((ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).w_))).REAL_;
_45i=(*(struct COMPLEX_*)((ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).w_))).IMAG_;
_46r=(*(double *)((ssize_t)sizeof(double)+d240_-((-2)*(ssize_t)sizeof(double))));
_47r=(*(struct COMPLEX_*)(2*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).w_))).REAL_;
_48i=(*(struct COMPLEX_*)(2*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).w_))).IMAG_;
_49r=(*(double *)(2*(ssize_t)sizeof(double)+d240_-((-2)*(ssize_t)sizeof(double))));
_50r=(*(struct COMPLEX_*)(3*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).w_))).REAL_;
_51i=(*(struct COMPLEX_*)(3*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).w_))).IMAG_;
_52r=((_37r*_38r)+(_40r*_41r)+(_43r*_44r)+(_46r*_47r)+(_49r*_50r));
_53i=((_37r*_39i)+(_40r*_42i)+(_43r*_45i)+(_46r*_48i)+(_49r*_51i));
_54r= - ni_/_18k2;
_55r=((-_16alfa*_36i)+(-_17beta*_53i));
_56i=((_16alfa*_35r)+(_17beta*_52r));
{register double temp=(_54r*_55r);(*RESULT_).IMAG_=(_54r*_56i);(*RESULT_).REAL_=temp;};
  }else{
    memset(RESULT_,0,(ssize_t)sizeof(struct COMPLEX_));
  };
}}}
void dnsdirect_12calcpn(struct COMPLEX_ *RESULT_,int ix_,int iz_){{
  double _16alfa;
/*alfa=_16alfa*/
double _17beta;
/*beta=_17beta*/
double _18k2;
/*k2=_18k2*/
_16alfa=alfa0_*(double)(ix_);
  _17beta=beta0_*(double)(iz_);

    
  _18k2=(_16alfa*_16alfa)+(_17beta*_17beta);

  {
  char *_19w;
_19w=iz_*dnsdata_46i+ix_*dnsdata_47i+V_;
if( _18k2>0.){
    double _20r;
double _21r;
double _22i;
double _23r;
double _24r;
double _25i;
double _26r;
double _27r;
double _28i;
double _29r;
double _30r;
double _31i;
double _32r;
double _33r;
double _34i;
double _35r;
double _36i;
double _37r;
double _38r;
double _39i;
double _40r;
double _41r;
double _42i;
double _43r;
double _44r;
double _45i;
double _46r;
double _47r;
double _48i;
double _49r;
double _50r;
double _51i;
double _52r;
double _53i;
double _54r;
double _55r;
double _56i;
_20r=(*(double *)((-2)*(ssize_t)sizeof(double)+d24n_-((-2)*(ssize_t)sizeof(double))));
_21r=(*(struct COMPLEX_*)((ny_-3)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).u_))).REAL_;
_22i=(*(struct COMPLEX_*)((ny_-3)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).u_))).IMAG_;
_23r=(*(double *)(-(ssize_t)sizeof(double)+d24n_-((-2)*(ssize_t)sizeof(double))));
_24r=(*(struct COMPLEX_*)((ny_-2)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).u_))).REAL_;
_25i=(*(struct COMPLEX_*)((ny_-2)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).u_))).IMAG_;
_26r=(*(double *)(d24n_-((-2)*(ssize_t)sizeof(double))));
_27r=(*(struct COMPLEX_*)((ny_-1)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).u_))).REAL_;
_28i=(*(struct COMPLEX_*)((ny_-1)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).u_))).IMAG_;
_29r=(*(double *)((ssize_t)sizeof(double)+d24n_-((-2)*(ssize_t)sizeof(double))));
_30r=(*(struct COMPLEX_*)(ny_*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).u_))).REAL_;
_31i=(*(struct COMPLEX_*)(ny_*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).u_))).IMAG_;
_32r=(*(double *)(2*(ssize_t)sizeof(double)+d24n_-((-2)*(ssize_t)sizeof(double))));
_33r=(*(struct COMPLEX_*)((ny_+1)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).u_))).REAL_;
_34i=(*(struct COMPLEX_*)((ny_+1)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).u_))).IMAG_;
_35r=((_20r*_21r)+(_23r*_24r)+(_26r*_27r)+(_29r*_30r)+(_32r*_33r));
_36i=((_20r*_22i)+(_23r*_25i)+(_26r*_28i)+(_29r*_31i)+(_32r*_34i));
_37r=(*(double *)((-2)*(ssize_t)sizeof(double)+d24n_-((-2)*(ssize_t)sizeof(double))));
_38r=(*(struct COMPLEX_*)((ny_-3)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).w_))).REAL_;
_39i=(*(struct COMPLEX_*)((ny_-3)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).w_))).IMAG_;
_40r=(*(double *)(-(ssize_t)sizeof(double)+d24n_-((-2)*(ssize_t)sizeof(double))));
_41r=(*(struct COMPLEX_*)((ny_-2)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).w_))).REAL_;
_42i=(*(struct COMPLEX_*)((ny_-2)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).w_))).IMAG_;
_43r=(*(double *)(d24n_-((-2)*(ssize_t)sizeof(double))));
_44r=(*(struct COMPLEX_*)((ny_-1)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).w_))).REAL_;
_45i=(*(struct COMPLEX_*)((ny_-1)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).w_))).IMAG_;
_46r=(*(double *)((ssize_t)sizeof(double)+d24n_-((-2)*(ssize_t)sizeof(double))));
_47r=(*(struct COMPLEX_*)(ny_*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).w_))).REAL_;
_48i=(*(struct COMPLEX_*)(ny_*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).w_))).IMAG_;
_49r=(*(double *)(2*(ssize_t)sizeof(double)+d24n_-((-2)*(ssize_t)sizeof(double))));
_50r=(*(struct COMPLEX_*)((ny_+1)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).w_))).REAL_;
_51i=(*(struct COMPLEX_*)((ny_+1)*(ssize_t)sizeof(struct VELOCITY_)+((char*)&(*(struct VELOCITY_*)(_19w)).w_))).IMAG_;
_52r=((_37r*_38r)+(_40r*_41r)+(_43r*_44r)+(_46r*_47r)+(_49r*_50r));
_53i=((_37r*_39i)+(_40r*_42i)+(_43r*_45i)+(_46r*_48i)+(_49r*_51i));
_54r= - ni_/_18k2;
_55r=((-_16alfa*_36i)+(-_17beta*_53i));
_56i=((_16alfa*_35r)+(_17beta*_52r));
{register double temp=(_54r*_55r);(*RESULT_).IMAG_=(_54r*_56i);(*RESULT_).REAL_=temp;};
  }else{
    memset(RESULT_,0,(ssize_t)sizeof(struct COMPLEX_));
  };
}}}

void dnsdirect_9linsolve(double lambda_){{
 {int ix_=0  ;while(ix_<=nx_){
  double _16alfa;
/*alfa=_16alfa*/
_16alfa=alfa0_*(double)(ix_);
  
    {int iz_=dnsdata_23l;do{
    double _17beta;
/*beta=_17beta*/
double _18k2;
/*k2=_18k2*/
_17beta=beta0_*(double)(iz_);
  
    _18k2=(_16alfa*_16alfa)+(_17beta*_17beta);

     {int iy_=dnsdata_2nyl  ;while(iy_<=dnsdata_3nyh ){  {int i_=(-2);do{ {
      struct dnsdata_s10 *_19w;
_19w=(struct dnsdata_s10 *)(iy_*(ssize_t)sizeof(struct dnsdata_s10)+derivatives_);
(*(double *)(i_*(ssize_t)sizeof(double)+iy_*(ssize_t)sizeof(double)*(2-(-2)+1)+iz_*dnsdata_25i+D2vmat_))=lambda_*((*(double *)(i_*(ssize_t)sizeof(double)+(*_19w).d2_-dnsdata_11st))-_18k2*(*(double *)(i_*(ssize_t)sizeof(double)+(*_19w).d0_-dnsdata_11st)))-(ni_*((*(double *)(i_*(ssize_t)sizeof(double)+(*_19w).d4_-dnsdata_11st))-2.*_18k2*(*(double *)(i_*(ssize_t)sizeof(double)+(*_19w).d2_-dnsdata_11st))+_18k2*_18k2*(*(double *)(i_*(ssize_t)sizeof(double)+(*_19w).d0_-dnsdata_11st))));
      (*(double *)(i_*(ssize_t)sizeof(double)+iy_*(ssize_t)sizeof(double)*(2-(-2)+1)+iz_*dnsdata_25i+etamat_))=lambda_*(*(double *)(i_*(ssize_t)sizeof(double)+(*_19w).d0_-dnsdata_11st))-(ni_*((*(double *)(i_*(ssize_t)sizeof(double)+(*_19w).d2_-dnsdata_11st))-_18k2*(*(double *)(i_*(ssize_t)sizeof(double)+(*_19w).d0_-dnsdata_11st)))) ;
    }i_++;}while(i_<=2);}iy_+=1;};}
    /*! condizioni al contorno*/
    if( (prev_==NULL) ){ {
      struct dnsdirect_s4 *_20w;
double _21r;
double _22r;
double _23r;
double _24i;
double _25r;
double _26r;
double _27r;
double _28i;
double _29r;
double _30r;
double _31r;
double _32i;
double _33r;
double _34r;
double _35r;
double _36i;
double _37r;
double _38r;
double _39r;
double _40i;
double _41r;
double _42r;
double _43r;
double _44i;
_20w=(struct dnsdirect_s4 *)(iz_*(ssize_t)sizeof(struct dnsdirect_s4)+ix_*dnsdirect_5i+bc0_);
if( (ix_==0 )&&( iz_==0 )){
        (*_20w).u_.REAL_=u_conv_+u0_;(*_20w).u_.IMAG_=0.;
        {register double temp=(*_20w).w_.REAL_+w_conv_+w0_;(*_20w).w_.IMAG_=(*_20w).w_.IMAG_;(*_20w).w_.REAL_=temp;};
 	memset(&(*_20w).v_,0,(ssize_t)sizeof(struct COMPLEX_));  memset(&(*_20w).vy_,0,(ssize_t)sizeof(struct COMPLEX_));  {register double temp=(*_20w).u_.REAL_+(-(*_20w).w_.IMAG_);(*_20w).eta_.IMAG_=(*_20w).u_.IMAG_+((*_20w).w_.REAL_);(*_20w).eta_.REAL_=temp;};
      }else{
        {register double temp=(_16alfa*(*_20w).u_.IMAG_)-(-_17beta*(*_20w).w_.IMAG_);(*_20w).vy_.IMAG_=(-_16alfa*(*_20w).u_.REAL_)-(_17beta*(*_20w).w_.REAL_);(*_20w).vy_.REAL_=temp;};
        {register double temp=(-_17beta*(*_20w).u_.IMAG_)-(-_16alfa*(*_20w).w_.IMAG_);(*_20w).eta_.IMAG_=(_17beta*(*_20w).u_.REAL_)-(_16alfa*(*_20w).w_.REAL_);(*_20w).eta_.REAL_=temp;};
      };
      _21r=(*(double *)((-2)*(ssize_t)sizeof(double)+v0bc_-(-(ssize_t)sizeof(double)*(2-(-2)+1)+dnsdata_11st)));
_22r=(1./((*(double *)((-2)*(ssize_t)sizeof(double)-(ssize_t)sizeof(double)*(2-(-2)+1)+v0bc_-(-(ssize_t)sizeof(double)*(2-(-2)+1)+dnsdata_11st)) )));
_23r=((*_20w).vy_.REAL_*_22r);
_24i=((*_20w).vy_.IMAG_*_22r);
{register double temp=(*_20w).v_.REAL_-(_21r*_23r);(*_20w).v_.IMAG_=(*_20w).v_.IMAG_-(_21r*_24i);(*_20w).v_.REAL_=temp;};
      dnsdata_32applybc_0(dnsdata_2nyl,(dnsdata_3nyh+2),(ssize_t)sizeof(double)*(2-(-2)+1),(-2),2,(ssize_t)sizeof(double),iz_*dnsdata_25i+D2vmat_,v0bc_-((-1)*(ssize_t)sizeof(double)*(2-(-2)+1)+dnsdata_11st));
      dnsdata_32applybc_0(dnsdata_2nyl,(dnsdata_3nyh+2),(ssize_t)sizeof(double)*(2-(-2)+1),(-2),2,(ssize_t)sizeof(double),iz_*dnsdata_25i+etamat_,eta0bc_-((-1)*(ssize_t)sizeof(double)*(2-(-2)+1)+dnsdata_11st))   ;
      _25r=(*(double *)((-2)*(ssize_t)sizeof(double)+(ssize_t)sizeof(double)*(2-(-2)+1)+iz_*dnsdata_25i+D2vmat_));
_26r=(1./((*(double *)((-2)*(ssize_t)sizeof(double)-(ssize_t)sizeof(double)*(2-(-2)+1)+v0bc_-(-(ssize_t)sizeof(double)*(2-(-2)+1)+dnsdata_11st)))));
_27r=((*_20w).vy_.REAL_*_26r);
_28i=((*_20w).vy_.IMAG_*_26r);
_29r=(*(double *)(-(ssize_t)sizeof(double)+(ssize_t)sizeof(double)*(2-(-2)+1)+iz_*dnsdata_25i+D2vmat_));
_30r=(1./((*(double *)(-(ssize_t)sizeof(double)+v0bc_-(-(ssize_t)sizeof(double)*(2-(-2)+1)+dnsdata_11st)))));
_31r=((*_20w).v_.REAL_*_30r);
_32i=((*_20w).v_.IMAG_*_30r);
{register double temp=(*(struct VELOCITY_*)((ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).v_.REAL_-(_25r*_27r)-(_29r*_31r);(*(struct VELOCITY_ *)((ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).v_.IMAG_=(*(struct VELOCITY_*)((ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).v_.IMAG_-(_25r*_28i)-(_29r*_32i);(*(struct VELOCITY_ *)((ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).v_.REAL_=temp;};
      _33r=(*(double *)((-2)*(ssize_t)sizeof(double)+2*(ssize_t)sizeof(double)*(2-(-2)+1)+iz_*dnsdata_25i+D2vmat_));
_34r=(1./((*(double *)(-(ssize_t)sizeof(double)+v0bc_-(-(ssize_t)sizeof(double)*(2-(-2)+1)+dnsdata_11st)))));
_35r=((*_20w).v_.REAL_*_34r);
_36i=((*_20w).v_.IMAG_*_34r);
{register double temp=(*(struct VELOCITY_*)(2*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).v_.REAL_-(_33r*_35r);(*(struct VELOCITY_ *)(2*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).v_.IMAG_=(*(struct VELOCITY_*)(2*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).v_.IMAG_-(_33r*_36i);(*(struct VELOCITY_ *)(2*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).v_.REAL_=temp;};
      _37r=(*(double *)(-(ssize_t)sizeof(double)+(ssize_t)sizeof(double)*(2-(-2)+1)+iz_*dnsdata_25i+etamat_));
_38r=(1./((*(double *)(-(ssize_t)sizeof(double)+eta0bc_-(-(ssize_t)sizeof(double)*(2-(-2)+1)+dnsdata_11st))          )));
_39r=((*_20w).eta_.REAL_*_38r);
_40i=((*_20w).eta_.IMAG_*_38r);
{register double temp=(*(struct VELOCITY_*)((ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).u_.REAL_-(_37r*_39r);(*(struct VELOCITY_ *)((ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).u_.IMAG_=(*(struct VELOCITY_*)((ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).u_.IMAG_-(_37r*_40i);(*(struct VELOCITY_ *)((ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).u_.REAL_=temp;};
      _41r=(*(double *)((-2)*(ssize_t)sizeof(double)+2*(ssize_t)sizeof(double)*(2-(-2)+1)+iz_*dnsdata_25i+etamat_));
_42r=(1./((*(double *)(-(ssize_t)sizeof(double)+eta0bc_-(-(ssize_t)sizeof(double)*(2-(-2)+1)+dnsdata_11st)))));
_43r=((*_20w).eta_.REAL_*_42r);
_44i=((*_20w).eta_.IMAG_*_42r);
{register double temp=(*(struct VELOCITY_*)(2*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).u_.REAL_-(_41r*_43r);(*(struct VELOCITY_ *)(2*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).u_.IMAG_=(*(struct VELOCITY_*)(2*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).u_.IMAG_-(_41r*_44i);(*(struct VELOCITY_ *)(2*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).u_.REAL_=temp;};
    }};
    if( (next_==NULL) ){ {
      struct dnsdirect_s4 *_45w;
double _46r;
double _47r;
double _48r;
double _49i;
double _50r;
double _51r;
double _52r;
double _53i;
double _54r;
double _55r;
double _56r;
double _57i;
double _58r;
double _59r;
double _60r;
double _61i;
double _62r;
double _63r;
double _64r;
double _65i;
double _66r;
double _67r;
double _68r;
double _69i;
_45w=(struct dnsdirect_s4 *)(iz_*(ssize_t)sizeof(struct dnsdirect_s4)+ix_*dnsdirect_5i+bcn_);
if( (ix_==0 )&&( iz_==0 )){
        (*_45w).u_.REAL_=u_conv_+un_;(*_45w).u_.IMAG_=0.;
        {register double temp=(*_45w).w_.REAL_+w_conv_+wn_;(*_45w).w_.IMAG_=(*_45w).w_.IMAG_;(*_45w).w_.REAL_=temp;};
	memset(&(*_45w).v_,0,(ssize_t)sizeof(struct COMPLEX_));  memset(&(*_45w).vy_,0,(ssize_t)sizeof(struct COMPLEX_));  {register double temp=(*_45w).u_.REAL_+(-(*_45w).w_.IMAG_);(*_45w).eta_.IMAG_=(*_45w).u_.IMAG_+((*_45w).w_.REAL_);(*_45w).eta_.REAL_=temp;};
      }else{
        {register double temp=(_16alfa*(*_45w).u_.IMAG_)-(-_17beta*(*_45w).w_.IMAG_);(*_45w).vy_.IMAG_=(-_16alfa*(*_45w).u_.REAL_)-(_17beta*(*_45w).w_.REAL_);(*_45w).vy_.REAL_=temp;};
        {register double temp=(-_17beta*(*_45w).u_.IMAG_)-(-_16alfa*(*_45w).w_.IMAG_);(*_45w).eta_.IMAG_=(_17beta*(*_45w).u_.REAL_)-(_16alfa*(*_45w).w_.REAL_);(*_45w).eta_.REAL_=temp;};
      };
      _46r=(*(double *)(2*(ssize_t)sizeof(double)+vnbc_-dnsdata_11st));
_47r=(1./((*(double *)(2*(ssize_t)sizeof(double)+(ssize_t)sizeof(double)*(2-(-2)+1)+vnbc_-dnsdata_11st))));
_48r=((*_45w).vy_.REAL_*_47r);
_49i=((*_45w).vy_.IMAG_*_47r);
{register double temp=(*_45w).v_.REAL_-(_46r*_48r);(*_45w).v_.IMAG_=(*_45w).v_.IMAG_-(_46r*_49i);(*_45w).v_.REAL_=temp;};
      dnsdata_33applybc_n(dnsdata_2nyl,(dnsdata_3nyh+2),(ssize_t)sizeof(double)*(2-(-2)+1),(-2),2,(ssize_t)sizeof(double),iz_*dnsdata_25i+D2vmat_,vnbc_-dnsdata_11st);
      dnsdata_33applybc_n(dnsdata_2nyl,(dnsdata_3nyh+2),(ssize_t)sizeof(double)*(2-(-2)+1),(-2),2,(ssize_t)sizeof(double),iz_*dnsdata_25i+etamat_,etanbc_-dnsdata_11st);
      _50r=(*(double *)(2*(ssize_t)sizeof(double)+(ny_-1)*(ssize_t)sizeof(double)*(2-(-2)+1)+iz_*dnsdata_25i+D2vmat_));
_51r=(1./((*(double *)(2*(ssize_t)sizeof(double)+(ssize_t)sizeof(double)*(2-(-2)+1)+vnbc_-dnsdata_11st))));
_52r=((*_45w).vy_.REAL_*_51r);
_53i=((*_45w).vy_.IMAG_*_51r);
_54r=(*(double *)((ssize_t)sizeof(double)+(ny_-1)*(ssize_t)sizeof(double)*(2-(-2)+1)+iz_*dnsdata_25i+D2vmat_));
_55r=(1./((*(double *)((ssize_t)sizeof(double)+vnbc_-dnsdata_11st))));
_56r=((*_45w).v_.REAL_*_55r);
_57i=((*_45w).v_.IMAG_*_55r);
{register double temp=(*(struct VELOCITY_*)((ny_-1)*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).v_.REAL_-(_50r*_52r)-(_54r*_56r);(*(struct VELOCITY_ *)((ny_-1)*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).v_.IMAG_=(*(struct VELOCITY_*)((ny_-1)*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).v_.IMAG_-(_50r*_53i)-(_54r*_57i);(*(struct VELOCITY_ *)((ny_-1)*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).v_.REAL_=temp;};
      _58r=(*(double *)(2*(ssize_t)sizeof(double)+(ny_-2)*(ssize_t)sizeof(double)*(2-(-2)+1)+iz_*dnsdata_25i+D2vmat_));
_59r=(1./((*(double *)((ssize_t)sizeof(double)+vnbc_-dnsdata_11st))));
_60r=((*_45w).v_.REAL_*_59r);
_61i=((*_45w).v_.IMAG_*_59r);
{register double temp=(*(struct VELOCITY_*)((ny_-2)*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).v_.REAL_-(_58r*_60r);(*(struct VELOCITY_ *)((ny_-2)*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).v_.IMAG_=(*(struct VELOCITY_*)((ny_-2)*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).v_.IMAG_-(_58r*_61i);(*(struct VELOCITY_ *)((ny_-2)*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).v_.REAL_=temp;};
      _62r=(*(double *)((ssize_t)sizeof(double)+(ny_-1)*(ssize_t)sizeof(double)*(2-(-2)+1)+iz_*dnsdata_25i+etamat_));
_63r=(1./((*(double *)((ssize_t)sizeof(double)+etanbc_-dnsdata_11st))));
_64r=((*_45w).eta_.REAL_*_63r);
_65i=((*_45w).eta_.IMAG_*_63r);
{register double temp=(*(struct VELOCITY_*)((ny_-1)*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).u_.REAL_-(_62r*_64r);(*(struct VELOCITY_ *)((ny_-1)*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).u_.IMAG_=(*(struct VELOCITY_*)((ny_-1)*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).u_.IMAG_-(_62r*_65i);(*(struct VELOCITY_ *)((ny_-1)*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).u_.REAL_=temp;};
      _66r=(*(double *)(2*(ssize_t)sizeof(double)+(ny_-2)*(ssize_t)sizeof(double)*(2-(-2)+1)+iz_*dnsdata_25i+etamat_));
_67r=(1./((*(double *)((ssize_t)sizeof(double)+etanbc_-dnsdata_11st)          )));
_68r=((*_45w).eta_.REAL_*_67r);
_69i=((*_45w).eta_.IMAG_*_67r);
{register double temp=(*(struct VELOCITY_*)((ny_-2)*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).u_.REAL_-(_66r*_68r);(*(struct VELOCITY_ *)((ny_-2)*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).u_.IMAG_=(*(struct VELOCITY_*)((ny_-2)*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).u_.IMAG_-(_66r*_69i);(*(struct VELOCITY_ *)((ny_-2)*(ssize_t)sizeof(struct VELOCITY_)+iz_*dnsdata_46i+ix_*dnsdata_47i+V_)).u_.REAL_=temp;};
    }};
    rbmatmod_1LUdecompStep(dnsdata_2nyl,(dnsdata_3nyh+2),(ssize_t)sizeof(double)*(2-(-2)+1),iz_*dnsdata_25i+D2vmat_);  rbmatmod_1LUdecompStep(dnsdata_2nyl,(dnsdata_3nyh+2),(ssize_t)sizeof(double)*(2-(-2)+1),iz_*dnsdata_25i+etamat_  );
    {
    char *_70w;
_70w=iz_*dnsdata_46i+ix_*dnsdata_47i+V_;
rbmatmod_2LeftLUDivStep1((dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(struct VELOCITY_),(char*)&(*(struct COMPLEX_ *)(&(*(struct VELOCITY_ *)(_70w)).v_)).REAL_,dnsdata_2nyl,(dnsdata_3nyh+2),(ssize_t)sizeof(double)*(2-(-2)+1),iz_*dnsdata_25i+D2vmat_,(dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(struct VELOCITY_),(char*)&(*(struct COMPLEX_ *)(&(*(struct VELOCITY_ *)(_70w)).v_)).REAL_);
    rbmatmod_2LeftLUDivStep1((dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(struct VELOCITY_),(char*)&(*(struct COMPLEX_ *)(&(*(struct VELOCITY_ *)(_70w)).v_)).IMAG_,dnsdata_2nyl,(dnsdata_3nyh+2),(ssize_t)sizeof(double)*(2-(-2)+1),iz_*dnsdata_25i+D2vmat_,(dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(struct VELOCITY_),(char*)&(*(struct COMPLEX_ *)(&(*(struct VELOCITY_ *)(_70w)).v_)).IMAG_)    ;
    rbmatmod_2LeftLUDivStep1((dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(struct VELOCITY_),(char*)&(*(struct COMPLEX_ *)(&(*(struct VELOCITY_ *)(_70w)).u_)).REAL_,dnsdata_2nyl,(dnsdata_3nyh+2),(ssize_t)sizeof(double)*(2-(-2)+1),iz_*dnsdata_25i+etamat_,(dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(struct VELOCITY_),(char*)&(*(struct COMPLEX_ *)(&(*(struct VELOCITY_ *)(_70w)).u_)).REAL_)       ;
    rbmatmod_2LeftLUDivStep1((dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(struct VELOCITY_),(char*)&(*(struct COMPLEX_ *)(&(*(struct VELOCITY_ *)(_70w)).u_)).IMAG_,dnsdata_2nyl,(dnsdata_3nyh+2),(ssize_t)sizeof(double)*(2-(-2)+1),iz_*dnsdata_25i+etamat_,(dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(struct VELOCITY_),(char*)&(*(struct COMPLEX_ *)(&(*(struct VELOCITY_ *)(_70w)).u_)).IMAG_)    ;
  }iz_++;}while(iz_<=dnsdata_24h);}
  {if( !(prev_==NULL) ){ fflush(prev_);};};
    {int iz_=dnsdata_44l;do{
    {
    char *_71w;
_71w=iz_*dnsdata_46i+ix_*dnsdata_47i+V_;
rbmatmod_3LeftLUDivStep2((dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(struct VELOCITY_),(char*)&(*(struct COMPLEX_ *)(&(*(struct VELOCITY_ *)(_71w)).v_)).REAL_,dnsdata_2nyl,(dnsdata_3nyh+2),(ssize_t)sizeof(double)*(2-(-2)+1),iz_*dnsdata_25i+D2vmat_);
    rbmatmod_3LeftLUDivStep2((dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(struct VELOCITY_),(char*)&(*(struct COMPLEX_ *)(&(*(struct VELOCITY_ *)(_71w)).v_)).IMAG_,dnsdata_2nyl,(dnsdata_3nyh+2),(ssize_t)sizeof(double)*(2-(-2)+1),iz_*dnsdata_25i+D2vmat_);
    rbmatmod_3LeftLUDivStep2((dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(struct VELOCITY_),(char*)&(*(struct COMPLEX_ *)(&(*(struct VELOCITY_ *)(_71w)).u_)).REAL_,dnsdata_2nyl,(dnsdata_3nyh+2),(ssize_t)sizeof(double)*(2-(-2)+1),iz_*dnsdata_25i+etamat_);
    rbmatmod_3LeftLUDivStep2((dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(struct VELOCITY_),(char*)&(*(struct COMPLEX_ *)(&(*(struct VELOCITY_ *)(_71w)).u_)).IMAG_,dnsdata_2nyl,(dnsdata_3nyh+2),(ssize_t)sizeof(double)*(2-(-2)+1),iz_*dnsdata_25i+etamat_);
    if( (next_==NULL) ){
      struct COMPLEX_ _73;
double _74r;
double _75i;
double _76r;
struct COMPLEX_ _78;
double _79r;
double _80i;
double _81r;
struct COMPLEX_ _83;
double _84r;
double _85i;
double _86r;
struct COMPLEX_ _88;
double _89r;
 memset(&_73,0,(ssize_t)sizeof(struct COMPLEX_)); {int i_= - 2  ;do{{double _72r;
_72r=(*(double *)(i_*(ssize_t)sizeof(double)+vnbc_-dnsdata_11st) );
{register double temp=(*&_73).REAL_+((*(struct COMPLEX_*)((ny_-1+i_)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).v_)).REAL_*_72r);(*&_73).IMAG_=(*&_73).IMAG_+((*(struct COMPLEX_*)((ny_-1+i_)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).v_)).IMAG_*_72r);(*&_73).REAL_=temp;};}i_+=1;}while(i_<=0);}_74r=((*(struct dnsdirect_s4*)(iz_*(ssize_t)sizeof(struct dnsdirect_s4)+ix_*dnsdirect_5i+bcn_)).v_.REAL_-_73.REAL_);
_75i=((*(struct dnsdirect_s4*)(iz_*(ssize_t)sizeof(struct dnsdirect_s4)+ix_*dnsdirect_5i+bcn_)).v_.IMAG_-_73.IMAG_);
_76r=(1./((*(double *)((ssize_t)sizeof(double)+vnbc_-dnsdata_11st))));
{register double temp=(_74r*_76r);(*(struct COMPLEX_ *)(ny_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).v_)  ).IMAG_=(_75i*_76r);(*(struct COMPLEX_ *)(ny_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).v_)  ).REAL_=temp;};
       memset(&_78,0,(ssize_t)sizeof(struct COMPLEX_)); {int i_= - 2  ;do{{double _77r;
_77r=(*(double *)(i_*(ssize_t)sizeof(double)+(ssize_t)sizeof(double)*(2-(-2)+1)+vnbc_-dnsdata_11st) );
{register double temp=(*&_78).REAL_+((*(struct COMPLEX_*)((ny_-1+i_)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).v_)).REAL_*_77r);(*&_78).IMAG_=(*&_78).IMAG_+((*(struct COMPLEX_*)((ny_-1+i_)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).v_)).IMAG_*_77r);(*&_78).REAL_=temp;};}i_+=1;}while(i_<=1);}_79r=((*(struct dnsdirect_s4*)(iz_*(ssize_t)sizeof(struct dnsdirect_s4)+ix_*dnsdirect_5i+bcn_)).vy_.REAL_-_78.REAL_);
_80i=((*(struct dnsdirect_s4*)(iz_*(ssize_t)sizeof(struct dnsdirect_s4)+ix_*dnsdirect_5i+bcn_)).vy_.IMAG_-_78.IMAG_);
_81r=(1./((*(double *)(2*(ssize_t)sizeof(double)+(ssize_t)sizeof(double)*(2-(-2)+1)+vnbc_-dnsdata_11st))));
{register double temp=(_79r*_81r);(*(struct COMPLEX_ *)((ny_+1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).v_)).IMAG_=(_80i*_81r);(*(struct COMPLEX_ *)((ny_+1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).v_)).REAL_=temp;};
       memset(&_83,0,(ssize_t)sizeof(struct COMPLEX_)); {int i_= - 2  ;do{{double _82r;
_82r=(*(double *)(i_*(ssize_t)sizeof(double)+etanbc_-dnsdata_11st) );
{register double temp=(*&_83).REAL_+((*(struct COMPLEX_*)((ny_-1+i_)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).u_)).REAL_*_82r);(*&_83).IMAG_=(*&_83).IMAG_+((*(struct COMPLEX_*)((ny_-1+i_)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).u_)).IMAG_*_82r);(*&_83).REAL_=temp;};}i_+=1;}while(i_<=0);}_84r=((*(struct dnsdirect_s4*)(iz_*(ssize_t)sizeof(struct dnsdirect_s4)+ix_*dnsdirect_5i+bcn_)).eta_.REAL_-_83.REAL_);
_85i=((*(struct dnsdirect_s4*)(iz_*(ssize_t)sizeof(struct dnsdirect_s4)+ix_*dnsdirect_5i+bcn_)).eta_.IMAG_-_83.IMAG_);
_86r=(1./((*(double *)((ssize_t)sizeof(double)+etanbc_-dnsdata_11st))));
{register double temp=(_84r*_86r);(*(struct COMPLEX_ *)(ny_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).u_)).IMAG_=(_85i*_86r);(*(struct COMPLEX_ *)(ny_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).u_)).REAL_=temp;};
       memset(&_88,0,(ssize_t)sizeof(struct COMPLEX_)); {int i_= - 2  ;do{{double _87r;
_87r=(*(double *)(i_*(ssize_t)sizeof(double)+(ssize_t)sizeof(double)*(2-(-2)+1)+etanbc_-dnsdata_11st) );
{register double temp=(*&_88).REAL_+((*(struct COMPLEX_*)((ny_-1+i_)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).u_)).REAL_*_87r);(*&_88).IMAG_=(*&_88).IMAG_+((*(struct COMPLEX_*)((ny_-1+i_)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).u_)).IMAG_*_87r);(*&_88).REAL_=temp;};}i_+=1;}while(i_<=1);}_89r=(1./((*(double *)(2*(ssize_t)sizeof(double)+(ssize_t)sizeof(double)*(2-(-2)+1)+etanbc_-dnsdata_11st))));
{register double temp=((-_88.REAL_)*_89r);(*(struct COMPLEX_ *)((ny_+1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).u_)).IMAG_=((-_88.IMAG_)*_89r);(*(struct COMPLEX_ *)((ny_+1)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).u_)).REAL_=temp;};
    };
    if( (prev_==NULL) ){
      struct COMPLEX_ _91;
double _92r;
double _93i;
double _94r;
struct COMPLEX_ _96;
double _97r;
double _98i;
double _99r;
struct COMPLEX_ _101;
double _102r;
double _103i;
double _104r;
struct COMPLEX_ _106;
double _107r;
 memset(&_91,0,(ssize_t)sizeof(struct COMPLEX_)); {int i_=0  ;do{{double _90r;
_90r=(*(double *)(i_*(ssize_t)sizeof(double)+v0bc_-(-(ssize_t)sizeof(double)*(2-(-2)+1)+dnsdata_11st)) );
{register double temp=(*&_91).REAL_+((*(struct COMPLEX_*)((1+i_)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).v_)).REAL_*_90r);(*&_91).IMAG_=(*&_91).IMAG_+((*(struct COMPLEX_*)((1+i_)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).v_)).IMAG_*_90r);(*&_91).REAL_=temp;};}i_+=1;}while(i_<=2);}_92r=((*(struct dnsdirect_s4*)(iz_*(ssize_t)sizeof(struct dnsdirect_s4)+ix_*dnsdirect_5i+bc0_)).v_.REAL_-_91.REAL_);
_93i=((*(struct dnsdirect_s4*)(iz_*(ssize_t)sizeof(struct dnsdirect_s4)+ix_*dnsdirect_5i+bc0_)).v_.IMAG_-_91.IMAG_);
_94r=(1./((*(double *)(-(ssize_t)sizeof(double)+v0bc_-(-(ssize_t)sizeof(double)*(2-(-2)+1)+dnsdata_11st)))));
{register double temp=(_92r*_94r);(*(struct COMPLEX_ *)(&(*(struct VELOCITY_*)(_71w)).v_) ).IMAG_=(_93i*_94r);(*(struct COMPLEX_ *)(&(*(struct VELOCITY_*)(_71w)).v_) ).REAL_=temp;};
       memset(&_96,0,(ssize_t)sizeof(struct COMPLEX_)); {int i_= - 1  ;do{{double _95r;
_95r=(*(double *)(i_*(ssize_t)sizeof(double)-(ssize_t)sizeof(double)*(2-(-2)+1)+v0bc_-(-(ssize_t)sizeof(double)*(2-(-2)+1)+dnsdata_11st)) );
{register double temp=(*&_96).REAL_+((*(struct COMPLEX_*)((1+i_)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).v_)).REAL_*_95r);(*&_96).IMAG_=(*&_96).IMAG_+((*(struct COMPLEX_*)((1+i_)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).v_)).IMAG_*_95r);(*&_96).REAL_=temp;};}i_+=1;}while(i_<=2);}_97r=((*(struct dnsdirect_s4*)(iz_*(ssize_t)sizeof(struct dnsdirect_s4)+ix_*dnsdirect_5i+bc0_)).vy_.REAL_-_96.REAL_);
_98i=((*(struct dnsdirect_s4*)(iz_*(ssize_t)sizeof(struct dnsdirect_s4)+ix_*dnsdirect_5i+bc0_)).vy_.IMAG_-_96.IMAG_);
_99r=(1./((*(double *)((-2)*(ssize_t)sizeof(double)-(ssize_t)sizeof(double)*(2-(-2)+1)+v0bc_-(-(ssize_t)sizeof(double)*(2-(-2)+1)+dnsdata_11st)))));
{register double temp=(_97r*_99r);(*(struct COMPLEX_ *)(-(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).v_)).IMAG_=(_98i*_99r);(*(struct COMPLEX_ *)(-(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).v_)).REAL_=temp;};
       memset(&_101,0,(ssize_t)sizeof(struct COMPLEX_)); {int i_=0  ;do{{double _100r;
_100r=(*(double *)(i_*(ssize_t)sizeof(double)+eta0bc_-(-(ssize_t)sizeof(double)*(2-(-2)+1)+dnsdata_11st)) );
{register double temp=(*&_101).REAL_+((*(struct COMPLEX_*)((1+i_)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).u_)).REAL_*_100r);(*&_101).IMAG_=(*&_101).IMAG_+((*(struct COMPLEX_*)((1+i_)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).u_)).IMAG_*_100r);(*&_101).REAL_=temp;};}i_+=1;}while(i_<=2);}_102r=((*(struct dnsdirect_s4*)(iz_*(ssize_t)sizeof(struct dnsdirect_s4)+ix_*dnsdirect_5i+bc0_)).eta_.REAL_-_101.REAL_);
_103i=((*(struct dnsdirect_s4*)(iz_*(ssize_t)sizeof(struct dnsdirect_s4)+ix_*dnsdirect_5i+bc0_)).eta_.IMAG_-_101.IMAG_);
_104r=(1./((*(double *)(-(ssize_t)sizeof(double)+eta0bc_-(-(ssize_t)sizeof(double)*(2-(-2)+1)+dnsdata_11st)))));
{register double temp=(_102r*_104r);(*(struct COMPLEX_ *)(&(*(struct VELOCITY_*)(_71w)).u_) ).IMAG_=(_103i*_104r);(*(struct COMPLEX_ *)(&(*(struct VELOCITY_*)(_71w)).u_) ).REAL_=temp;};
       memset(&_106,0,(ssize_t)sizeof(struct COMPLEX_)); {int i_= - 1  ;do{{double _105r;
_105r=(*(double *)(i_*(ssize_t)sizeof(double)-(ssize_t)sizeof(double)*(2-(-2)+1)+eta0bc_-(-(ssize_t)sizeof(double)*(2-(-2)+1)+dnsdata_11st)) );
{register double temp=(*&_106).REAL_+((*(struct COMPLEX_*)((1+i_)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).u_)).REAL_*_105r);(*&_106).IMAG_=(*&_106).IMAG_+((*(struct COMPLEX_*)((1+i_)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).u_)).IMAG_*_105r);(*&_106).REAL_=temp;};}i_+=1;}while(i_<=2);}_107r=(1./((*(double *)((-2)*(ssize_t)sizeof(double)-(ssize_t)sizeof(double)*(2-(-2)+1)+eta0bc_-(-(ssize_t)sizeof(double)*(2-(-2)+1)+dnsdata_11st)))));
{register double temp=((-_106.REAL_)*_107r);(*(struct COMPLEX_ *)(-(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).u_)).IMAG_=((-_106.IMAG_)*_107r);(*(struct COMPLEX_ *)(-(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_71w)).u_)).REAL_=temp;};
    };
  }iz_++;}while(iz_<=dnsdata_45h);}
  {if( !(next_==NULL) ){ fflush(next_);};};
ix_+=1;};} 
}}

void dnsdirect_14deriv(int f0_l,int f0_h,size_t f0_i,char *f0__,int f1_l,int f1_h,size_t f1_i,char *f1__){{
if( (prev_==NULL) ){
   (*(double *)(f1__))=0.; {int i_= - 2  ;do{{(*(double *)(f1__))+=(*(double *)(i_*(ssize_t)sizeof(double)+d140_-((-2)*(ssize_t)sizeof(double))))*(*(double *)((1+i_)*f0_i+f0__)) ;}i_+=1;}while(i_<=2);}
   (*(double *)(-f1_i+f1__))=0.; {int i_= - 2  ;do{{(*(double *)(-f1_i+f1__))+=(*(double *)(i_*(ssize_t)sizeof(double)+d14m1_-((-2)*(ssize_t)sizeof(double))))*(*(double *)((1+i_)*f0_i+f0__)) ;}i_+=1;}while(i_<=2);}
};
if( (next_==NULL) ){  
   (*(double *)(ny_*f1_i+f1__))=0.; {int i_= - 2  ;do{{(*(double *)(ny_*f1_i+f1__))+=(*(double *)(i_*(ssize_t)sizeof(double)+d14n_-((-2)*(ssize_t)sizeof(double))))*(*(double *)((ny_-1+i_)*f0_i+f0__)) ;}i_+=1;}while(i_<=2);}
   (*(double *)((ny_+1)*f1_i+f1__))=0.; {int i_= - 2  ;do{{(*(double *)((ny_+1)*f1_i+f1__))+=(*(double *)(i_*(ssize_t)sizeof(double)+d14np1_-((-2)*(ssize_t)sizeof(double))))*(*(double *)((ny_-1+i_)*f0_i+f0__)) ;}i_+=1;}while(i_<=2);}
};
 {int i_=dnsdata_2nyl  ;do{{ {struct dnsdata_s10 *_16w;
_16w=(struct dnsdata_s10 *)(i_*(ssize_t)sizeof(struct dnsdata_s10)+derivatives_) ;
(*(double *)(i_*f1_i+f1__))=( (*(double *)((-2)*(ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st))*(*(double *)((-2)*f0_i+(i_*f0_i+f0__)))+(*(double *)(-(ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st))*(*(double *)(-f0_i+(i_*f0_i+f0__)))+(*(double *)((*_16w).d1_-dnsdata_11st))*(*(double *)((i_*f0_i+f0__)))+(*(double *)((ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st))*(*(double *)(f0_i+(i_*f0_i+f0__)))+(*(double *)(2*(ssize_t)sizeof(double)+(*_16w).d1_-dnsdata_11st))*(*(double *)(2*f0_i+(i_*f0_i+f0__)))) ;}}i_+=1;}while(i_<=dnsdata_3nyh);}
if( (prev_==NULL) ){
  { struct dnsdata_s10 *_17w;
_17w=(struct dnsdata_s10 *)((ssize_t)sizeof(struct dnsdata_s10)+derivatives_);
(*(double *)(f1_i+f1__))-=((*(double *)(-(ssize_t)sizeof(double)+(*_17w).d0_-dnsdata_11st))*(*(double *)(f1__))+(*(double *)((-2)*(ssize_t)sizeof(double)+(*_17w).d0_-dnsdata_11st))*(*(double *)(-f1_i+f1__)));}
  { struct dnsdata_s10 *_18w;
_18w=(struct dnsdata_s10 *)(2*(ssize_t)sizeof(struct dnsdata_s10)+derivatives_);
(*(double *)(2*f1_i+f1__))-=(*(double *)((-2)*(ssize_t)sizeof(double)+(*_18w).d0_-dnsdata_11st))*(*(double *)(f1__));}
};
if( (next_==NULL) ){  
  { struct dnsdata_s10 *_19w;
_19w=(struct dnsdata_s10 *)((ny_-1)*(ssize_t)sizeof(struct dnsdata_s10)+derivatives_);
(*(double *)((ny_-1)*f1_i+f1__))-=((*(double *)((ssize_t)sizeof(double)+(*_19w).d0_-dnsdata_11st))*(*(double *)(ny_*f1_i+f1__))+(*(double *)(2*(ssize_t)sizeof(double)+(*_19w).d0_-dnsdata_11st))*(*(double *)((ny_+1)*f1_i+f1__)));}
  { struct dnsdata_s10 *_20w;
_20w=(struct dnsdata_s10 *)((ny_-2)*(ssize_t)sizeof(struct dnsdata_s10)+derivatives_);
(*(double *)((ny_-2)*f1_i+f1__))-=(*(double *)(2*(ssize_t)sizeof(double)+(*_20w).d0_-dnsdata_11st))*(*(double *)(ny_*f1_i+f1__));}
};
rbmatmod_2LeftLUDivStep1(f1_l,f1_h,f1_i,f1__,dnsdata_2nyl,(dnsdata_3nyh+2),(ssize_t)sizeof(double)*(2-(-2)+1),D0mat_,f1_l,f1_h,f1_i,f1__);
}}

void dnsdirect_10vetaTOuvw(){{
/*! Remember: eta=+I*beta*u-I*alfa*w*/
{ {char *_17_;
_17_=(dnsdata_2nyl-2)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(V_)).u_)).IMAG_;
 {char *_16_=(dnsdata_2nyl-2)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(V_)).w_)).REAL_;int _16_1=(dnsdata_3nyh+2)-(dnsdata_2nyl-2); do{{ (*(double *)(_16_))=(*(double *)(_17_)); _17_ =(ssize_t)sizeof(struct VELOCITY_)+_17_;}_16_+=(ssize_t)sizeof(struct VELOCITY_);_16_1--;}while(_16_1>=0);}} { {char *_18_=(dnsdata_2nyl-2)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(V_)).u_)).IMAG_;int _18_1=(dnsdata_3nyh+2)-(dnsdata_2nyl-2); do{{ (*(double *)(_18_))=0.;}_18_+=(ssize_t)sizeof(struct VELOCITY_);_18_1--;}while(_18_1>=0);}} { {char *_19_=(dnsdata_2nyl-2)*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(V_)).w_)).IMAG_;int _19_1=(dnsdata_3nyh+2)-(dnsdata_2nyl-2); do{{ (*(double *)(_19_))=0.;}_19_+=(ssize_t)sizeof(struct VELOCITY_);_19_1--;}while(_19_1>=0);}}}
  {int ix_=dnsdata_43h;do{{int iz_=dnsdata_44l;do{ if(!((ix_==0 )&&( iz_==0))){
  {
  char *_20w;
_20w=iz_*dnsdata_46i+ix_*dnsdata_47i+V_;
dnsdirect_14deriv((dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(struct VELOCITY_),(char*)&(*(struct COMPLEX_ *)(&(*(struct VELOCITY_ *)(_20w)).v_)).REAL_,(dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(struct VELOCITY_),(char*)&(*(struct COMPLEX_ *)(&(*(struct VELOCITY_ *)(_20w)).w_)).REAL_);
  dnsdirect_14deriv((dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(struct VELOCITY_),(char*)&(*(struct COMPLEX_ *)(&(*(struct VELOCITY_ *)(_20w)).v_)).IMAG_,(dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(struct VELOCITY_),(char*)&(*(struct COMPLEX_ *)(&(*(struct VELOCITY_ *)(_20w)).w_)).IMAG_);
}}iz_++;}while(iz_<=dnsdata_45h);}ix_--;}while(ix_>=0);}
{if( !(prev_==NULL) ){ fflush(prev_);};};
  {int ix_=dnsdata_43h;do{
  double _21alfa;
/*alfa=_21alfa*/
_21alfa=alfa0_*(double)(ix_);

    {int iz_=dnsdata_44l;do{ if(!((ix_==0 )&&( iz_==0))){
    {
    char *_22w;
double _23beta;
/*beta=_23beta*/
double _24k2;
/*k2=_24k2*/
_22w=iz_*dnsdata_46i+ix_*dnsdata_47i+V_;
rbmatmod_3LeftLUDivStep2((dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(struct VELOCITY_),(char*)&(*(struct COMPLEX_ *)(&(*(struct VELOCITY_ *)(_22w)).w_)).REAL_,dnsdata_2nyl,(dnsdata_3nyh+2),(ssize_t)sizeof(double)*(2-(-2)+1),D0mat_);
    rbmatmod_3LeftLUDivStep2((dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(struct VELOCITY_),(char*)&(*(struct COMPLEX_ *)(&(*(struct VELOCITY_ *)(_22w)).w_)).IMAG_,dnsdata_2nyl,(dnsdata_3nyh+2),(ssize_t)sizeof(double)*(2-(-2)+1),D0mat_);
    {if( !(next_==NULL) ){ fflush(next_);};};
    _23beta=beta0_*(double)(iz_);
  _24k2=(_21alfa*_21alfa)+(_23beta*_23beta);

      {int iy_=(dnsdata_2nyl-2);do{{ double _25r;
double _26i;
double _27r;
double _28r;
double _29i;
/*temp=(-_29i),_28r*/
double _31r;
double _32i;
double _33r;
double _34r;
double _35i;
_25r=((_21alfa*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_22w)).w_)).REAL_)-(_23beta*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_22w)).u_)).REAL_));
_26i=((_21alfa*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_22w)).w_)).IMAG_)-(_23beta*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_22w)).u_)).IMAG_));
_27r=(1./(_24k2));
_28r=(_25r*_27r);
_29i=(_26i*_27r);

      _31r=((_23beta*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_22w)).w_)).REAL_)+(_21alfa*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_22w)).u_)).REAL_));
_32i=((_23beta*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_22w)).w_)).IMAG_)+(_21alfa*(*(struct COMPLEX_*)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_22w)).u_)).IMAG_));
_33r=(1./(_24k2 ));
_34r=(_31r*_33r);
_35i=(_32i*_33r);
{register double temp=(-_35i);(*(struct COMPLEX_ *)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_22w)).w_)).IMAG_=_34r;(*(struct COMPLEX_ *)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_22w)).w_)).REAL_=temp;};
      {register double temp=(-_29i);(*(struct COMPLEX_ *)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_22w)).u_)).IMAG_=_28r;(*(struct COMPLEX_ *)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(_22w)).u_)).REAL_=temp;};
    }iy_++;}while(iy_<=(dnsdata_3nyh+2));}
  }}iz_++;}while(iz_<=dnsdata_45h);}
ix_--;}while(ix_>=0);}
{if( !(next_==NULL) ){ fflush(next_);};};
}}

void dnsdirect_13computeflowrate(double lambda_){{
  size_t _16i;
ssize_t _17st;
char *ucorr_;
struct _s25{double u_;double w_;double uc_;double wc_;};struct _s25 fr_;   
_16i=(ssize_t)sizeof(double)*((dnsdata_3nyh+2)-(dnsdata_2nyl-2)+1);
_17st=(dnsdata_2nyl-2)*(ssize_t)sizeof(double);

ucorr_=malloc(_16i);if(ucorr_==NULL)errmalloc();ucorr_-=_17st;memset(_17st+ucorr_,0,_16i);	
   {int iy_=dnsdata_2nyl  ;do{{ (*(double *)(iy_*(ssize_t)sizeof(double)+ucorr_))=1.;}iy_+=1;}while(iy_<=dnsdata_3nyh   );}
   {int iy_=dnsdata_2nyl  ;while(iy_<=dnsdata_3nyh ){{
     struct dnsdata_s10 *_18w;
_18w=(struct dnsdata_s10 *)(iy_*(ssize_t)sizeof(struct dnsdata_s10)+derivatives_);
{  {int _19i_=(-2);do{{(*(double *)(_19i_*(ssize_t)sizeof(double)+iy_*(ssize_t)sizeof(double)*(2-(-2)+1)+etamat_))=lambda_*(*(double *)(_19i_*(ssize_t)sizeof(double)+(*_18w).d0_-dnsdata_11st))-ni_*(*(double *)(_19i_*(ssize_t)sizeof(double)+(*_18w).d2_-dnsdata_11st)) ;}_19i_++;}while(_19i_<=2);}}
  }iy_+=1;};}   
  if( (prev_==NULL) ){ dnsdata_32applybc_0(dnsdata_2nyl,(dnsdata_3nyh+2),(ssize_t)sizeof(double)*(2-(-2)+1),(-2),2,(ssize_t)sizeof(double),etamat_,eta0bc_-((-1)*(ssize_t)sizeof(double)*(2-(-2)+1)+dnsdata_11st));};
  if( (next_==NULL)  ){ dnsdata_33applybc_n(dnsdata_2nyl,(dnsdata_3nyh+2),(ssize_t)sizeof(double)*(2-(-2)+1),(-2),2,(ssize_t)sizeof(double),etamat_,etanbc_-dnsdata_11st);};
  rbmatmod_1LUdecompStep(dnsdata_2nyl,(dnsdata_3nyh+2),(ssize_t)sizeof(double)*(2-(-2)+1),etamat_);
  rbmatmod_2LeftLUDivStep1((dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(double),ucorr_,dnsdata_2nyl,(dnsdata_3nyh+2),(ssize_t)sizeof(double)*(2-(-2)+1),etamat_,(dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(double),ucorr_);
  {if( !(prev_==NULL) ){ fflush(prev_);};};
  rbmatmod_3LeftLUDivStep2((dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(double),ucorr_,dnsdata_2nyl,(dnsdata_3nyh+2),(ssize_t)sizeof(double)*(2-(-2)+1),etamat_);
  {if( !(next_==NULL) ){ fflush(next_);};};
  if( (next_==NULL) ){
    double _21;
double _22;
 _21=0.; {int i_= - 2  ;do{{(*&_21)+=(*(double *)((ny_-1+i_)*(ssize_t)sizeof(double)+ucorr_))*(*(double *)(i_*(ssize_t)sizeof(double)+etanbc_-dnsdata_11st) );}i_+=1;}while(i_<=0);}(*(double *)(ny_*(ssize_t)sizeof(double)+ucorr_)  )=( - _21)/(*(double *)((ssize_t)sizeof(double)+etanbc_-dnsdata_11st));
     _22=0.; {int i_= - 2  ;do{{(*&_22)+=(*(double *)((ny_-1+i_)*(ssize_t)sizeof(double)+ucorr_))*(*(double *)(i_*(ssize_t)sizeof(double)+(ssize_t)sizeof(double)*(2-(-2)+1)+etanbc_-dnsdata_11st) );}i_+=1;}while(i_<=1);}(*(double *)((ny_+1)*(ssize_t)sizeof(double)+ucorr_))=( - _22)/(*(double *)(2*(ssize_t)sizeof(double)+(ssize_t)sizeof(double)*(2-(-2)+1)+etanbc_-dnsdata_11st));
  };
  if( (prev_==NULL) ){
    double _23;
double _24;
 _23=0.; {int i_=0  ;do{{(*&_23)+=(*(double *)((1+i_)*(ssize_t)sizeof(double)+ucorr_))*(*(double *)(i_*(ssize_t)sizeof(double)+eta0bc_-(-(ssize_t)sizeof(double)*(2-(-2)+1)+dnsdata_11st)) );}i_+=1;}while(i_<=2);}(*(double *)(ucorr_) )=( - _23)/(*(double *)(-(ssize_t)sizeof(double)+eta0bc_-(-(ssize_t)sizeof(double)*(2-(-2)+1)+dnsdata_11st)));
     _24=0.; {int i_= - 1  ;do{{(*&_24)+=(*(double *)((1+i_)*(ssize_t)sizeof(double)+ucorr_))*(*(double *)(i_*(ssize_t)sizeof(double)-(ssize_t)sizeof(double)*(2-(-2)+1)+eta0bc_-(-(ssize_t)sizeof(double)*(2-(-2)+1)+dnsdata_11st)) );}i_+=1;}while(i_<=2);}(*(double *)(-(ssize_t)sizeof(double)+ucorr_))=( - _24)/(*(double *)((-2)*(ssize_t)sizeof(double)-(ssize_t)sizeof(double)*(2-(-2)+1)+eta0bc_-(-(ssize_t)sizeof(double)*(2-(-2)+1)+dnsdata_11st)));
  };
  
  {
    {int iy_=(dnsdata_2nyl-2);do{{ memset((struct COMPLEX_ *)(iy_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct VELOCITY_*)(V_)).v_),0,(ssize_t)sizeof(struct COMPLEX_)); }iy_++;}while(iy_<=(dnsdata_3nyh+2));}
  if( !(prev_==NULL) ){   if(!(fread(0+&fr_ ,(ssize_t)sizeof(struct _s25),1, prev_ )==1))ioerr( prev_ );}else{ memset(&fr_,0,(ssize_t)sizeof(struct _s25));};
  dnsdata_34yintegr(&fr_.u_ ,(dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(struct VELOCITY_),(char*)&(*(struct COMPLEX_ *)(&(*(struct VELOCITY_ *)(V_)).u_)).REAL_);
  dnsdata_34yintegr(&fr_.w_ ,(dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(struct VELOCITY_),(char*)&(*(struct COMPLEX_ *)(&(*(struct VELOCITY_ *)(V_)).w_)).REAL_);
  dnsdata_34yintegr(&fr_.uc_ ,(dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(double),ucorr_)      ;
  if( !(next_==NULL) ){ fwrite(&fr_,(ssize_t)sizeof(struct _s25),1,next_); fflush(next_);   if(!(fread(0+&fr_,(ssize_t)sizeof(struct _s25),1, next_ )==1))ioerr( next_ );};
  if( !(prev_==NULL) ){ fwrite(&fr_,(ssize_t)sizeof(struct _s25),1,prev_); fflush(prev_);};
   px_ = 6.*ni_/fr_.u_   /*!Constant Power Input*/;
/*!  pz = (1-gamma)*6*ni/fr.w   !Constant Power Input*/
  if( meanflowx_ != -99.){ 
    px_=0.;  corrpx_=(meanflowx_+u_conv_*(ymax_-ymin_)-fr_.u_)/fr_.uc_  /*! Constant Q*/;
    {  {int _27i_=(dnsdata_2nyl-2);do{{(*(double *)(_27i_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(V_)).u_)).REAL_))+=corrpx_*(*(double *)(_27i_*(ssize_t)sizeof(double)+ucorr_)) ;}_27i_++;}while(_27i_<=(dnsdata_3nyh+2));}}
  };
  if( meanflowz_ != -99.){ 
    pz_=0.;  corrpz_=(meanflowz_+w_conv_*(ymax_-ymin_)-fr_.w_)/fr_.uc_  /*! Constant Q*/;
    {  {int _29i_=(dnsdata_2nyl-2);do{{(*(double *)(_29i_*(ssize_t)sizeof(struct VELOCITY_)+(char*)&(*(struct COMPLEX_*)(&(*(struct VELOCITY_*)(V_)).w_)).REAL_))+=corrpz_*(*(double *)(_29i_*(ssize_t)sizeof(double)+ucorr_)) ;}_29i_++;}while(_29i_<=(dnsdata_3nyh+2));}}
  };
  if( meanpx_ != -99.){  px_ = meanpx_ /*! Constant Px*/;};
  if( meanpz_ != -99.){  pz_ = meanpz_ /*! Constant Pz*/;};
  if( !(prev_==NULL) ){   if(!(fread(0+&fr_ ,(ssize_t)sizeof(struct _s25),1, prev_ )==1))ioerr( prev_ );}else{ memset(&fr_,0,(ssize_t)sizeof(struct _s25));};
  dnsdata_34yintegr(&fr_.u_ ,(dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(struct VELOCITY_),(char*)&(*(struct COMPLEX_ *)(&(*(struct VELOCITY_ *)(V_)).u_)).REAL_);
  dnsdata_34yintegr(&fr_.w_ ,(dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(struct VELOCITY_),(char*)&(*(struct COMPLEX_ *)(&(*(struct VELOCITY_ *)(V_)).w_)).REAL_);
  dnsdata_34yintegr(&fr_.uc_ ,(dnsdata_2nyl-2),(dnsdata_3nyh+2),(ssize_t)sizeof(double),ucorr_)      ;
  if( !(next_==NULL) ){ fwrite(&fr_,(ssize_t)sizeof(struct _s25),1,next_); fflush(next_);   if(!(fread(0+&fr_,(ssize_t)sizeof(struct _s25),1, next_ )==1))ioerr( next_ );};
  if( !(prev_==NULL) ){ fwrite(&fr_,(ssize_t)sizeof(struct _s25),1,prev_); fflush(prev_);};
  flowx_=fr_.u_-u_conv_*(ymax_-ymin_);  flowz_=fr_.w_-w_conv_*(ymax_-ymin_);
  if( (prev_==NULL) ){ 
      u1zero_=0.; {int i_= - 2  ;do{{(*&u1zero_)+=(*(double *)(i_*(ssize_t)sizeof(double)+d140_-((-2)*(ssize_t)sizeof(double))))*(*(struct VELOCITY_ *)((i_+1)*(ssize_t)sizeof(struct VELOCITY_)+V_)).u_.REAL_ ;}i_+=1;}while(i_<=2);}
      w1zero_=0.; {int i_= - 2  ;do{{(*&w1zero_)+=(*(double *)(i_*(ssize_t)sizeof(double)+d140_-((-2)*(ssize_t)sizeof(double))))*(*(struct VELOCITY_ *)((i_+1)*(ssize_t)sizeof(struct VELOCITY_)+V_)).w_.REAL_ ;}i_+=1;}while(i_<=2);}
  };
  if( !(prev_==NULL) ){   if(!(fread(0+&u1zero_,(ssize_t)sizeof(double),1, prev_ )==1&& fread(0+&w1zero_,(ssize_t)sizeof(double),1, prev_ )==1))ioerr( prev_ ); fflush(prev_);};
  if( !(next_==NULL)  ){ fwrite(&u1zero_,(ssize_t)sizeof(double),1,next_); fwrite(&w1zero_,(ssize_t)sizeof(double),1,next_); fflush(next_);};
}free(ucorr_+_17st);}}

