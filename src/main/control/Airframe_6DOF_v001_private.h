/*
 * Airframe_6DOF_v001_private.h
 *
 * Trial License - for use to evaluate programs for possible purchase as
 * an end-user only.
 *
 * Code generation for model "Airframe_6DOF_v001".
 *
 * Model version              : 1.0
 * Simulink Coder version : 9.3 (R2020a) 18-Nov-2019
 * C source code generated on : Wed Apr 29 23:00:24 2020
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: ARM Compatible->ARM Cortex
 * Code generation objective: Execution efficiency
 * Validation result: Not run
 */

#ifndef RTW_HEADER_Airframe_6DOF_v001_private_h_
#define RTW_HEADER_Airframe_6DOF_v001_private_h_
#include "rtwtypes.h"
#include "multiword_types.h"

/* Private macros used by the generated code to access rtModel */
#ifndef rtmIsMajorTimeStep
# define rtmIsMajorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
# define rtmIsMinorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmSetTPtr
# define rtmSetTPtr(rtm, val)          ((rtm)->Timing.t = (val))
#endif

#ifndef ATMOS_TYPEDEF

typedef enum { COESA = 1, MILHDBK310, MILSTD210C } AtmosTypeIdx;

typedef enum { PROFILE = 1, ENVELOPE } ModelIdx;

typedef enum { HIGHTEMP = 1, LOWTEMP, HIGHDENSITY,
  LOWDENSITY, HIGHPRESSURE, LOWPRESSURE } VarIdx;

typedef enum { PP1 = 1, PP10 } PPercentIdx;

typedef enum { K5 = 1, K10, K20, K30, K40 } PAltIdx;

typedef enum { EXTREME = 1, P1, P5, P10, P20 } EPercentIdx;

#define ATMOS_TYPEDEF
#endif                                 /* ATMOS_TYPEDEF */

#ifndef ATMOS_DEFINE
#define PRESSURE0                      101325.0                  /*  N/m^2                  */
#define TEMPERATURE0                   288.15                    /*  K                      */
#define GRAV_CONST                     9.80665                   /*  m/s^2                  */
#define MOL_WT                         28.9644                   /*  kg/kgmol (air)         */
#define R_HAT                          8314.32                   /*  J/kgmol.K (gas const.) */
#define GAMMA                          1.4                       /*  (specific heat ratio) */
#define GMR                            ( GRAV_CONST * MOL_WT / R_HAT )
#define ATMOS_DEFINE
#endif                                 /* ATMOS_DEFINE */

#ifndef COESA76_DEFINE_DATA

/* 1976 COESA atmosphere model */
#define NUM1976PTS                     8

static real_T altitude76[NUM1976PTS] = {/* in meters (m) */
  0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0, 84852.0 };

static real_T tempGradient76[NUM1976PTS] = {/* in K/m  */
  (-0.0065), 0.0, 0.0010, 0.0028, 0.0, -0.0028, -0.0020, -0.0020 };

#define COESA76_DEFINE_DATA
#endif                                 /* COESA76_DEFINE_DATA */

extern real_T rt_urand_Upu32_Yd_f_pw_snf(uint32_T *u);
extern real_T rt_nrand_Upu32_Yd_f_pw_snf(uint32_T *u);
extern real_T rt_powd_snf(real_T u0, real_T u1);
extern real_T rt_atan2d_snf(real_T u0, real_T u1);
extern void rt_mrdivide_U1d1x3_U2d_9vOrDY9Z(const real_T u0[3], const real_T u1
  [9], real_T y[3]);
extern real_T rt_modd_snf(real_T u0, real_T u1);
void InitCalcAtmosCOESA(real_T *temperature76, real_T *pressureRatio76);
void CalcAtmosCOESA(const real_T *altitude, real_T *temp, real_T *pressure,
                    real_T *density, real_T *speedofsound, real_T *temperature76,
                    real_T *pressureRatio76, int_T numPoints);
extern uint32_T plook_bincpa(real_T u, const real_T bp[], uint32_T maxIndex,
  real_T *fraction, uint32_T *prevIndex);
extern real_T intrp2d_la(const uint32_T bpIndex[], const real_T frac[], const
  real_T table[], const uint32_T stride, const uint32_T maxIndex[]);
extern uint32_T binsearch_u32d_prevIdx(real_T u, const real_T bp[], uint32_T
  startIndex, uint32_T maxIndex);

/* private model entry point functions */
extern void Airframe_6DOF_v001_derivatives(void);

#endif                            /* RTW_HEADER_Airframe_6DOF_v001_private_h_ */
