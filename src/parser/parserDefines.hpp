

// #ifdef __AVX
// #define BYTE_ALIGNMENT 32
// #else
// #define BYTE_ALIGNMENT 16
// #endif





#define MAX_TIP_EV     0.999999999 /* max tip vector value, sum of EVs needs to be smaller than 1.0, otherwise the numerics break down */
#define smoothings     32          /* maximum smoothing passes through tree */
#define iterations     10          /* maximum iterations of iterations per insert */
#define newzpercycle   1           /* iterations of makenewz per tree traversal */
#define nmlngth        256         /* number of characters in species name */
#define deltaz         0.00001     /* test of net branch length change in update */
#define defaultz       0.9         /* value of z assigned as starting point */
#define unlikely       -1.0E300    /* low likelihood for initialization */


#define SUMMARIZE_LENGTH -3
#define SUMMARIZE_LH     -2
#define NO_BRANCHES      -1

#define MASK_LENGTH 32
#define GET_BITVECTOR_LENGTH(x) ((x % MASK_LENGTH) ? (x / MASK_LENGTH + 1) : (x / MASK_LENGTH))

#define zmin       1.0E-15  /* max branch prop. to -log(zmin) (= 34) */
#define zmax (1.0 - 1.0E-6) /* min branch prop. to 1.0-zmax (= 1.0E-6) */

#define twotothe256  \
  115792089237316195423570985008687907853269984665640564039457584007913129639936.0  
                                                     /*  2**256 (exactly)  */

#define minlikelihood  (1.0/twotothe256)
#define minusminlikelihood -minlikelihood




/* 18446744073709551616.0 */

/*4294967296.0*/

/* 18446744073709551616.0 */

/*  2**64 (exactly)  */
/* 4294967296 2**32 */

#define badRear         -1

// #define NUM_BRANCHES     16

#define TRUE             1
#define FALSE            0



#define LIKELIHOOD_EPSILON 0.0000001

#define AA_SCALE 10.0
#define AA_SCALE_PLUS_EPSILON 10.001

/* ALPHA_MIN is critical -> numerical instability, eg for 4 discrete rate cats                    */
/* and alpha = 0.01 the lowest rate r_0 is                                                        */
/* 0.00000000000000000000000000000000000000000000000000000000000034878079110511010487             */
/* which leads to numerical problems Table for alpha settings below:                              */
/*                                                                                                */
/* 0.010000 0.00000000000000000000000000000000000000000000000000000000000034878079110511010487    */
/* 0.010000 yielded nasty numerical bugs in at least one case !                                   */
/* 0.020000 0.00000000000000000000000000000044136090435925743185910935350715027016962154188875    */
/* 0.030000 0.00000000000000000000476844846859006690412039180149775802624789852441798419292220    */
/* 0.040000 0.00000000000000049522423236954066431210260930029681736928018820007024736185030633    */
/* 0.050000 0.00000000000050625351310359203371872643495343928538368616365517027588794007897377    */
/* 0.060000 0.00000000005134625283884191118711474021861409372524676086868566926568746566772461    */
/* 0.070000 0.00000000139080650074206434685544624965062437960128249869740102440118789672851562    */
/* 0.080000 0.00000001650681201563587066858709818343436959153791576682124286890029907226562500    */
/* 0.090000 0.00000011301977332931251259273962858978301859735893231118097901344299316406250000    */
/* 0.100000 0.00000052651925834844387815526344648331402709118265192955732345581054687500000000    */


#define ALPHA_MIN    0.02
#define ALPHA_MAX    1000.0

#define RATE_MIN     0.0000001
#define RATE_MAX     1000000.0

#define INVAR_MIN    0.0001
#define INVAR_MAX    0.9999

#define TT_MIN       0.0000001
#define TT_MAX       1000000.0

#define FREQ_MIN     0.001

/* 
   previous values between 0.001 and 0.000001

   TO AVOID NUMERICAL PROBLEMS WHEN FREQ == 0 IN PARTITIONED MODELS, ESPECIALLY WITH AA 
   previous value of FREQ_MIN was: 0.000001, but this seemed to cause problems with some 
   of the 7-state secondary structure models with some rather exotic small toy test datasets,
   on the other hand 0.001 caused problems with some of the 16-state secondary structure models

   For some reason the frequency settings seem to be repeatedly causing numerical problems
   
*/

#define ITMAX 100



#define SHFT(a,b,c,d)                (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b)                    ((b) > 0.0 ? fabs(a) : -fabs(a))

#define ABS(x)    (((x)<0)   ?  (-(x)) : (x))
#define MIN(x,y)  (((x)<(y)) ?    (x)  : (y))
#define MAX(x,y)  (((x)>(y)) ?    (x)  : (y))
#define NINT(x)   ((int) ((x)>0 ? ((x)+0.5) : ((x)-0.5)))

#ifdef _USE_FPGA_LOG
extern double log_approx (double input);
#define LOG(x)  log_approx(x)
#else
#define LOG(x)  log(x)
#endif


#ifdef _USE_FPGA_EXP
extern double exp_approx (double x);
#define EXP(x)  exp_approx(x)
#else
#define EXP(x)  exp(x)
#endif


#define LOGF(x) logf(x)


#define PointGamma(prob,alpha,beta)  PointChi2(prob,2.0*(alpha))/(2.0*(beta))

// #define programName        "the phylip file parser for ExaML"
// #define programVersion     "1.0.1"
// #define programDate        "March 29 2013"


#define  TREE_EVALUATION            0
#define  BIG_RAPID_MODE             1
#define  CALC_BIPARTITIONS          3
#define  SPLIT_MULTI_GENE           4
#define  CHECK_ALIGNMENT            5
#define  PER_SITE_LL                6
#define  PARSIMONY_ADDITION         7
#define  CLASSIFY_ML                9
#define  DISTANCE_MODE              11
#define  GENERATE_BS                12
#define  COMPUTE_ELW                13
#define  BOOTSTOP_ONLY              14
#define  COMPUTE_LHS                17
#define  COMPUTE_BIPARTITION_CORRELATION 18
#define  THOROUGH_PARSIMONY         19
#define  COMPUTE_RF_DISTANCE        20
#define  MORPH_CALIBRATOR           21
#define  CONSENSUS_ONLY             22
#define  MESH_TREE_SEARCH           23
#define  FAST_SEARCH                24
#define  MORPH_CALIBRATOR_PARSIMONY 25
#define  SH_LIKE_SUPPORTS           28

#define M_GTRCAT         1
#define M_GTRGAMMA       2
#define M_BINCAT         3
#define M_BINGAMMA       4
#define M_PROTCAT        5
#define M_PROTGAMMA      6
#define M_32CAT          7
#define M_32GAMMA        8
#define M_64CAT          9
#define M_64GAMMA        10


#define DAYHOFF    0
#define DCMUT      1
#define JTT        2
#define MTREV      3
#define WAG        4
#define RTREV      5
#define CPREV      6
#define VT         7
#define BLOSUM62   8
#define MTMAM      9
#define LG         10
#define MTART      11
#define MTZOA      12
#define PMB        13
#define HIVB       14
#define HIVW       15
#define JTTDCMUT   16
#define FLU        17 
#define AUTO       18
#define LG4        19
#define GTR        20  /* GTR always needs to be the last one */

#define NUM_PROT_MODELS 21

/* bipartition stuff */

#define BIPARTITIONS_ALL       0
#define GET_BIPARTITIONS_BEST  1
#define DRAW_BIPARTITIONS_BEST 2
#define BIPARTITIONS_BOOTSTOP  3
#define BIPARTITIONS_RF  4



/* bootstopping stuff */

#define BOOTSTOP_PERMUTATIONS 100
#define START_BSTOP_TEST      10

#define FC_THRESHOLD          99
#define FC_SPACING            50
#define FC_LOWER              0.99
#define FC_INIT               20

#define FREQUENCY_STOP 0
#define MR_STOP        1
#define MRE_STOP       2
#define MRE_IGN_STOP   3

#define MR_CONSENSUS 0
#define MRE_CONSENSUS 1
#define STRICT_CONSENSUS 2



/* bootstopping stuff end */


#define TIP_TIP     0
#define TIP_INNER   1
#define INNER_INNER 2

#define MIN_MODEL        -1
#define BINARY_DATA      0
#define DNA_DATA         1
#define AA_DATA          2
#define SECONDARY_DATA   3
#define SECONDARY_DATA_6 4
#define SECONDARY_DATA_7 5
#define GENERIC_32       6
#define GENERIC_64       7
#define MAX_MODEL        8

#define SEC_6_A 0
#define SEC_6_B 1
#define SEC_6_C 2
#define SEC_6_D 3
#define SEC_6_E 4

#define SEC_7_A 5
#define SEC_7_B 6
#define SEC_7_C 7
#define SEC_7_D 8
#define SEC_7_E 9
#define SEC_7_F 10

#define SEC_16   11
#define SEC_16_A 12
#define SEC_16_B 13
#define SEC_16_C 14
#define SEC_16_D 15
#define SEC_16_E 16
#define SEC_16_F 17
#define SEC_16_I 18
#define SEC_16_J 19
#define SEC_16_K 20

#define ORDERED_MULTI_STATE 0
#define MK_MULTI_STATE      1
#define GTR_MULTI_STATE     2





#define CAT         0
#define GAMMA       1
#define GAMMA_I     2


#define PCF 32


#define REARR_SETTING 1
#define FAST_SPRS     2
#define SLOW_SPRS     3
