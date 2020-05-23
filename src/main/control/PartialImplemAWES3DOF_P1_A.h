/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: PartialImplemAWES3DOF_P1_A.h
 *
 * Code generated for Simulink model 'PartialImplemAWES3DOF_P1_A'.
 *
 * Model version                  : 1.156
 * Simulink Coder version         : 9.3 (R2020a) 18-Nov-2019
 * C/C++ source code generated on : Fri May 22 18:47:53 2020
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: ARM Compatible->ARM Cortex
 * Code generation objective: Execution efficiency
 * Validation result: Not run
 */

#ifndef RTW_HEADER_PartialImplemAWES3DOF_P1_A_h_
#define RTW_HEADER_PartialImplemAWES3DOF_P1_A_h_
#include <string.h>
#include <math.h>
#ifndef PartialImplemAWES3DOF_P1_A_COMMON_INCLUDES_
# define PartialImplemAWES3DOF_P1_A_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#endif                         /* PartialImplemAWES3DOF_P1_A_COMMON_INCLUDES_ */

#include "PartialImplemAWES3DOF_P1_A_types.h"
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtGetInf.h"

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
  real_T Transpose[9];                 /* '<S12>/Transpose' */
  real_T rtb_Transpose_m[9];
  real_T sincos_o1_j[3];               /* '<S8>/sincos' */
  real_T sincos_o2[3];                 /* '<S15>/sincos' */
  real_T Product_b[3];                 /* '<S1>/Product' */
  real_T Sum1[2];                      /* '<S7>/Sum1' */
  real_T TmpSignalConversionAtPositionIn[2];/* '<S7>/Transform  to Earth Axes' */
  real_T frac[2];
  real_T q;                            /* '<S7>/q' */
  real_T Sum4;                         /* '<Root>/Sum4' */
  real_T Product2;                     /* '<S9>/Product2' */
  real_T Va;                           /* '<S2>/Va' */
  real_T TrigonometricFunction;        /* '<S1>/Trigonometric Function' */
  real_T u2rhoV2;                      /* '<S19>/1//2rhoV^2' */
  real_T CM;                           /* '<S22>/CM' */
  real_T Gain_j;                       /* '<S6>/Gain' */
  real_T Product1_l;                   /* '<S6>/Product1' */
  real_T UnaryMinus;                   /* '<S11>/Unary Minus' */
  real_T CM_el;                        /* '<S23>/CM_el' */
  real_T gain2_p;                      /* '<S20>/gain2 ' */
} B_PartialImplemAWES3DOF_P1_A_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  uint32_T m_bpIndex;                  /* '<S5>/Thrust' */
  int_T Uw_IWORK;                      /* '<S7>/U,w' */
} DW_PartialImplemAWES3DOF_P1_A_T;

/* Continuous states (default storage) */
typedef struct {
  real_T theta;                        /* '<S7>/Theta' */
  real_T u[2];                         /* '<S7>/U,w' */
  real_T Xe[2];                        /* '<S7>/Position' */
  real_T q;                            /* '<S7>/q' */
} X_PartialImplemAWES3DOF_P1_A_T;

/* Periodic continuous state vector (global) */
typedef int_T PeriodicIndX_PartialImplemAWE_T[1];
typedef real_T PeriodicRngX_PartialImplemAWE_T[2];

/* State derivatives (default storage) */
typedef struct {
  real_T theta;                        /* '<S7>/Theta' */
  real_T u[2];                         /* '<S7>/U,w' */
  real_T Xe[2];                        /* '<S7>/Position' */
  real_T q;                            /* '<S7>/q' */
} XDot_PartialImplemAWES3DOF_P1_T;

/* State disabled  */
typedef struct {
  boolean_T theta;                     /* '<S7>/Theta' */
  boolean_T u[2];                      /* '<S7>/U,w' */
  boolean_T Xe[2];                     /* '<S7>/Position' */
  boolean_T q;                         /* '<S7>/q' */
} XDis_PartialImplemAWES3DOF_P1_T;

/* Invariant block signals (default storage) */
typedef struct {
  const real_T Gain[3];                /* '<S1>/Gain' */
  const real_T Gain1;                  /* '<S2>/Gain1' */
  const real_T Thrust;                 /* '<S5>/Thrust' */
} ConstB_PartialImplemAWES3DOF__T;

#ifndef ODE3_INTG
#define ODE3_INTG

/* ODE3 Integration Data */
typedef struct {
  real_T *y;                           /* output */
  real_T *f[3];                        /* derivatives */
} ODE3_IntgData;

#endif

/* Constant parameters (default storage) */
typedef struct {
  /* Pooled Parameter (Expression: aeroCoeff.alpha_vec)
   * Referenced by:
   *   '<S22>/Prelookup'
   *   '<S23>/Prelookup1'
   */
  real_T pooled3[10];

  /* Expression: aeroCoeff.CD
   * Referenced by: '<S22>/CD'
   */
  real_T CD_Table[10];

  /* Expression: aeroCoeff.CL
   * Referenced by: '<S22>/CL'
   */
  real_T CL_Table[10];

  /* Expression: aeroCoeff.def_vec
   * Referenced by: '<S23>/Prelookup'
   */
  real_T Prelookup_BreakpointsData[5];

  /* Expression: aeroCoeff.CD_el
   * Referenced by: '<S23>/CD_el'
   */
  real_T CD_el_Table[50];

  /* Expression: aeroCoeff.CL_el
   * Referenced by: '<S23>/CL_el'
   */
  real_T CL_el_Table[5];

  /* Expression: aeroCoeff.CM
   * Referenced by: '<S22>/CM'
   */
  real_T CM_Table[10];

  /* Expression: aeroCoeff.CM_el
   * Referenced by: '<S23>/CM_el'
   */
  real_T CM_el_Table[5];
} ConstP_PartialImplemAWES3DOF__T;

/* External inputs (root inport signals with default storage) */
typedef struct {
  real_T AltCmd;                       /* '<Root>/Elevator Command' */
} ExtU_PartialImplemAWES3DOF_P1_T;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  real_T Altitude;                     /* '<Root>/Altitude' */
} ExtY_PartialImplemAWES3DOF_P1_T;

/* Real-time Model Data Structure */
struct tag_RTM_PartialImplemAWES3DOF_T {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;
  X_PartialImplemAWES3DOF_P1_A_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[6];
  real_T odeF[3][6];
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
    boolean_T firstInitCondFlag;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

/* Block signals (default storage) */
extern B_PartialImplemAWES3DOF_P1_A_T PartialImplemAWES3DOF_P1_A_B;

/* Continuous states (default storage) */
extern X_PartialImplemAWES3DOF_P1_A_T PartialImplemAWES3DOF_P1_A_X;

/* Block states (default storage) */
extern DW_PartialImplemAWES3DOF_P1_A_T PartialImplemAWES3DOF_P1_A_DW;

/* External inputs (root inport signals with default storage) */
extern ExtU_PartialImplemAWES3DOF_P1_T PartialImplemAWES3DOF_P1_A_U;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY_PartialImplemAWES3DOF_P1_T PartialImplemAWES3DOF_P1_A_Y;
extern const ConstB_PartialImplemAWES3DOF__T PartialImplemAWES3DOF_P1_ConstB;/* constant block i/o */

/* Constant parameters (default storage) */
extern const ConstP_PartialImplemAWES3DOF__T PartialImplemAWES3DOF_P1_ConstP;

/* Model entry point functions */
extern void PartialImplemAWES3DOF_P1_A_initialize(void);
extern void PartialImplemAWES3DOF_P1_A_step(void);
extern void PartialImplemAWES3DOF_P1_A_terminate(void);

/* Real-time Model object */
extern RT_MODEL_PartialImplemAWES3DO_T *const PartialImplemAWES3DOF_P1_A_M;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<S14>/Unit Conversion' : Unused code path elimination
 * Block '<S21>/Airspeed' : Unused code path elimination
 * Block '<S25>/Product' : Unused code path elimination
 * Block '<S25>/Product1' : Unused code path elimination
 * Block '<S25>/Sum' : Unused code path elimination
 * Block '<Root>/Altitude_' : Unused code path elimination
 * Block '<S29>/Constant' : Unused code path elimination
 * Block '<S29>/Constant1' : Unused code path elimination
 * Block '<S29>/Constant2' : Unused code path elimination
 * Block '<S31>/Reshape1' : Unused code path elimination
 * Block '<S31>/pg,qg,rg' : Unused code path elimination
 * Block '<S6>/Product2' : Unused code path elimination
 * Block '<S16>/Reshape (9) to [3x3] column-major' : Reshape block reduction
 * Block '<S13>/Unit Conversion' : Eliminated nontunable gain of 1
 * Block '<S17>/Reshape (9) to [3x3] column-major' : Reshape block reduction
 * Block '<S30>/Reshape1' : Reshape block reduction
 * Block '<S31>/Reshape' : Reshape block reduction
 */

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
 * '<Root>' : 'PartialImplemAWES3DOF_P1_A'
 * '<S1>'   : 'PartialImplemAWES3DOF_P1_A/3DOF EOM'
 * '<S2>'   : 'PartialImplemAWES3DOF_P1_A/Aerodynamic Model'
 * '<S3>'   : 'PartialImplemAWES3DOF_P1_A/Degrees to Radians'
 * '<S4>'   : 'PartialImplemAWES3DOF_P1_A/Environment Model'
 * '<S5>'   : 'PartialImplemAWES3DOF_P1_A/Propulsion'
 * '<S6>'   : 'PartialImplemAWES3DOF_P1_A/TetherForce'
 * '<S7>'   : 'PartialImplemAWES3DOF_P1_A/3DOF EOM/3DOF (Body Axes)'
 * '<S8>'   : 'PartialImplemAWES3DOF_P1_A/3DOF EOM/Rotation Angles to Direction Cosine Matrix'
 * '<S9>'   : 'PartialImplemAWES3DOF_P1_A/3DOF EOM/3DOF (Body Axes)/Calculate qdot'
 * '<S10>'  : 'PartialImplemAWES3DOF_P1_A/3DOF EOM/3DOF (Body Axes)/Determine Force,  Mass & Inertia'
 * '<S11>'  : 'PartialImplemAWES3DOF_P1_A/3DOF EOM/3DOF (Body Axes)/Gravity'
 * '<S12>'  : 'PartialImplemAWES3DOF_P1_A/3DOF EOM/3DOF (Body Axes)/Transform  to Earth Axes'
 * '<S13>'  : 'PartialImplemAWES3DOF_P1_A/3DOF EOM/3DOF (Body Axes)/Velocity Conversion'
 * '<S14>'  : 'PartialImplemAWES3DOF_P1_A/3DOF EOM/3DOF (Body Axes)/Velocity Conversion1'
 * '<S15>'  : 'PartialImplemAWES3DOF_P1_A/3DOF EOM/3DOF (Body Axes)/Transform  to Earth Axes/Rotation Angles to Direction Cosine Matrix'
 * '<S16>'  : 'PartialImplemAWES3DOF_P1_A/3DOF EOM/3DOF (Body Axes)/Transform  to Earth Axes/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S17>'  : 'PartialImplemAWES3DOF_P1_A/3DOF EOM/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S18>'  : 'PartialImplemAWES3DOF_P1_A/Aerodynamic Model/Aero Forces and Moments'
 * '<S19>'  : 'PartialImplemAWES3DOF_P1_A/Aerodynamic Model/Dynamic Pressure'
 * '<S20>'  : 'PartialImplemAWES3DOF_P1_A/Aerodynamic Model/Elevator Forces and Moments'
 * '<S21>'  : 'PartialImplemAWES3DOF_P1_A/Aerodynamic Model/Incidence  & Airspeed'
 * '<S22>'  : 'PartialImplemAWES3DOF_P1_A/Aerodynamic Model/Lookup Aero Coeffs'
 * '<S23>'  : 'PartialImplemAWES3DOF_P1_A/Aerodynamic Model/Lookup Elevator Aero Coeffs'
 * '<S24>'  : 'PartialImplemAWES3DOF_P1_A/Aerodynamic Model/Dynamic Pressure/dot'
 * '<S25>'  : 'PartialImplemAWES3DOF_P1_A/Aerodynamic Model/Incidence  & Airspeed/dot'
 * '<S26>'  : 'PartialImplemAWES3DOF_P1_A/Aerodynamic Model/Lookup Aero Coeffs/Subsystem'
 * '<S27>'  : 'PartialImplemAWES3DOF_P1_A/Aerodynamic Model/Lookup Elevator Aero Coeffs/Subsystem'
 * '<S28>'  : 'PartialImplemAWES3DOF_P1_A/Aerodynamic Model/Lookup Elevator Aero Coeffs/Subsystem1'
 * '<S29>'  : 'PartialImplemAWES3DOF_P1_A/Environment Model/Constant Atmosphere'
 * '<S30>'  : 'PartialImplemAWES3DOF_P1_A/Environment Model/Constant Gravity'
 * '<S31>'  : 'PartialImplemAWES3DOF_P1_A/Environment Model/Constant Wind'
 */
#endif                            /* RTW_HEADER_PartialImplemAWES3DOF_P1_A_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
