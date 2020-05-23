/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: PartialImplemAWES3DOF_P1_B.h
 *
 * Code generated for Simulink model 'PartialImplemAWES3DOF_P1_B'.
 *
 * Model version                  : 1.21
 * Simulink Coder version         : 9.3 (R2020a) 18-Nov-2019
 * C/C++ source code generated on : Fri May 22 19:31:49 2020
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: ARM Compatible->ARM Cortex
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_PartialImplemAWES3DOF_P1_B_h_
#define RTW_HEADER_PartialImplemAWES3DOF_P1_B_h_
#include <string.h>
#ifndef PartialImplemAWES3DOF_P1_B_COMMON_INCLUDES_
# define PartialImplemAWES3DOF_P1_B_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#endif                         /* PartialImplemAWES3DOF_P1_B_COMMON_INCLUDES_ */

#include "PartialImplemAWES3DOF_P1_B_types.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContStateDisabled
# define rtmGetContStateDisabled(rtm)  ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
# define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
# define rtmGetContStates(rtm)         ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
# define rtmSetContStates(rtm, val)    ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
# define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
# define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
# define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
# define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetIntgData
# define rtmGetIntgData(rtm)           ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
# define rtmSetIntgData(rtm, val)      ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
# define rtmGetOdeF(rtm)               ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
# define rtmSetOdeF(rtm, val)          ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
# define rtmGetOdeY(rtm)               ((rtm)->odeY)
#endif

#ifndef rtmSetOdeY
# define rtmSetOdeY(rtm, val)          ((rtm)->odeY = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
# define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
# define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
# define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
# define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetZCCacheNeedsReset
# define rtmGetZCCacheNeedsReset(rtm)  ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
# define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
# define rtmGetdX(rtm)                 ((rtm)->derivs)
#endif

#ifndef rtmSetdX
# define rtmSetdX(rtm, val)            ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
# define rtmGetStopRequested(rtm)      ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
# define rtmSetStopRequested(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
# define rtmGetStopRequestedPtr(rtm)   (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
# define rtmGetT(rtm)                  (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTPtr
# define rtmGetTPtr(rtm)               ((rtm)->Timing.t)
#endif

/* Block signals (default storage) */
typedef struct {
  real_T Subtract;                     /* '<Root>/Subtract' */
  real_T Sum;                          /* '<Root>/Sum' */
  real_T FilterCoefficient;            /* '<S39>/Filter Coefficient' */
  real_T Sum_i;                        /* '<S45>/Sum' */
  real_T IntegralGain;                 /* '<S33>/Integral Gain' */
} B_PartialImplemAWES3DOF_P1_B_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  real_T UD_DSTATE;                    /* '<S3>/UD' */
} DW_PartialImplemAWES3DOF_P1_B_T;

/* Continuous states (default storage) */
typedef struct {
  real_T Plant_CSTATE[2];              /* '<Root>/Plant' */
  real_T Controller_CSTATE;            /* '<S1>/Controller' */
  real_T Integrator_CSTATE;            /* '<S36>/Integrator' */
  real_T Filter_CSTATE;                /* '<S31>/Filter' */
} X_PartialImplemAWES3DOF_P1_B_T;

/* State derivatives (default storage) */
typedef struct {
  real_T Plant_CSTATE[2];              /* '<Root>/Plant' */
  real_T Controller_CSTATE;            /* '<S1>/Controller' */
  real_T Integrator_CSTATE;            /* '<S36>/Integrator' */
  real_T Filter_CSTATE;                /* '<S31>/Filter' */
} XDot_PartialImplemAWES3DOF_P1_T;

/* State disabled  */
typedef struct {
  boolean_T Plant_CSTATE[2];           /* '<Root>/Plant' */
  boolean_T Controller_CSTATE;         /* '<S1>/Controller' */
  boolean_T Integrator_CSTATE;         /* '<S36>/Integrator' */
  boolean_T Filter_CSTATE;             /* '<S31>/Filter' */
} XDis_PartialImplemAWES3DOF_P1_T;

/* Invariant block signals (default storage) */
typedef struct {
  const real_T Altitude_Selector2;     /* '<Root>/Altitude_Selector2' */
  const real_T Theta_Angle_Selector;   /* '<Root>/Theta_Angle_Selector' */
  const real_T Gain1;                  /* '<S2>/Gain1' */
} ConstB_PartialImplemAWES3DOF__T;

#ifndef ODE3_INTG
#define ODE3_INTG

/* ODE3 Integration Data */
typedef struct {
  real_T *y;                           /* output */
  real_T *f[3];                        /* derivatives */
} ODE3_IntgData;

#endif

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  real_T AltitudeCommand;              /* '<Root>/Altitude Command' */
  real_T AltError;                     /* '<Root>/AltError' */
  real_T Altitude;                     /* '<Root>/Altitude' */
  real_T AltitudeError;                /* '<Root>/Altitude Error' */
} ExtY_PartialImplemAWES3DOF_P1_T;

/* Real-time Model Data Structure */
struct tag_RTM_PartialImplemAWES3DOF_T {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;
  X_PartialImplemAWES3DOF_P1_B_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[5];
  real_T odeF[3][5];
  ODE3_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    time_T stepSize0;
    uint32_T clockTick1;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

/* Block signals (default storage) */
extern B_PartialImplemAWES3DOF_P1_B_T PartialImplemAWES3DOF_P1_B_B;

/* Continuous states (default storage) */
extern X_PartialImplemAWES3DOF_P1_B_T PartialImplemAWES3DOF_P1_B_X;

/* Block states (default storage) */
extern DW_PartialImplemAWES3DOF_P1_B_T PartialImplemAWES3DOF_P1_B_DW;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY_PartialImplemAWES3DOF_P1_T PartialImplemAWES3DOF_P1_B_Y;
extern const ConstB_PartialImplemAWES3DOF__T PartialImplemAWES3DOF_P1_ConstB;/* constant block i/o */

/* Model entry point functions */
extern void PartialImplemAWES3DOF_P1_B_initialize(void);
extern void PartialImplemAWES3DOF_P1_B_step(void);
extern void PartialImplemAWES3DOF_P1_B_terminate(void);

/* Real-time Model object */
extern RT_MODEL_PartialImplemAWES3DO_T *const PartialImplemAWES3DOF_P1_B_M;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'PartialImplemAWES3DOF_P1_B'
 * '<S1>'   : 'PartialImplemAWES3DOF_P1_B/Continous Controller'
 * '<S2>'   : 'PartialImplemAWES3DOF_P1_B/Degrees to Radians'
 * '<S3>'   : 'PartialImplemAWES3DOF_P1_B/Continous Controller/Difference Filter'
 * '<S4>'   : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller'
 * '<S5>'   : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Anti-windup'
 * '<S6>'   : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/D Gain'
 * '<S7>'   : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Filter'
 * '<S8>'   : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Filter ICs'
 * '<S9>'   : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/I Gain'
 * '<S10>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Ideal P Gain'
 * '<S11>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Ideal P Gain Fdbk'
 * '<S12>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Integrator'
 * '<S13>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Integrator ICs'
 * '<S14>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/N Copy'
 * '<S15>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/N Gain'
 * '<S16>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/P Copy'
 * '<S17>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Parallel P Gain'
 * '<S18>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Reset Signal'
 * '<S19>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Saturation'
 * '<S20>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Saturation Fdbk'
 * '<S21>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Sum'
 * '<S22>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Sum Fdbk'
 * '<S23>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Tracking Mode'
 * '<S24>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Tracking Mode Sum'
 * '<S25>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Tsamp - Integral'
 * '<S26>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Tsamp - Ngain'
 * '<S27>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/postSat Signal'
 * '<S28>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/preSat Signal'
 * '<S29>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Anti-windup/Passthrough'
 * '<S30>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/D Gain/Internal Parameters'
 * '<S31>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Filter/Cont. Filter'
 * '<S32>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Filter ICs/Internal IC - Filter'
 * '<S33>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/I Gain/Internal Parameters'
 * '<S34>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Ideal P Gain/Passthrough'
 * '<S35>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Ideal P Gain Fdbk/Disabled'
 * '<S36>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Integrator/Continuous'
 * '<S37>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Integrator ICs/Internal IC'
 * '<S38>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/N Copy/Disabled'
 * '<S39>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/N Gain/Internal Parameters'
 * '<S40>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/P Copy/Disabled'
 * '<S41>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Parallel P Gain/Internal Parameters'
 * '<S42>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Reset Signal/Disabled'
 * '<S43>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Saturation/Passthrough'
 * '<S44>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Saturation Fdbk/Disabled'
 * '<S45>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Sum/Sum_PID'
 * '<S46>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Sum Fdbk/Disabled'
 * '<S47>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Tracking Mode/Disabled'
 * '<S48>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Tracking Mode Sum/Passthrough'
 * '<S49>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Tsamp - Integral/Passthrough'
 * '<S50>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/Tsamp - Ngain/Passthrough'
 * '<S51>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/postSat Signal/Forward_Path'
 * '<S52>'  : 'PartialImplemAWES3DOF_P1_B/Continous Controller/PID Controller/preSat Signal/Forward_Path'
 */
#endif                            /* RTW_HEADER_PartialImplemAWES3DOF_P1_B_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
