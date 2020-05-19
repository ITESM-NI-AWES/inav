/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep.h
 *
 * Code generated for Simulink model 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep'.
 *
 * Model version                  : 1.155
 * Simulink Coder version         : 9.3 (R2020a) 18-Nov-2019
 * C/C++ source code generated on : Tue May  5 13:04:57 2020
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_h_
#define RTW_HEADER_FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_h_
#include <string.h>
#include <math.h>
#ifndef FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_COMMON_INCLUDES_
# define FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#endif
     /* FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_COMMON_INCLUDES_ */

#include "FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_types.h"
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtGetInf.h"
#include "common/time.h"
/* Macros for accessing real-time model data structure */
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
extern uint32_t stepNumber;
/* Block signals (default storage) */
typedef struct {
  real_T q;                            /* '<S13>/q' */
  real_T Sum1[2];                      /* '<S13>/Sum1' */
  real_T TmpSignalConversionAtPositionIn[2];/* '<S13>/Transform  to Earth Axes' */
  real_T Sum4;                         /* '<Root>/Sum4' */
  real_T Product2;                     /* '<S15>/Product2' */
  real_T startTime;                    /* '<S9>/startTime' */
  real_T startTime_a;                  /* '<S10>/startTime' */
} B_FinWing_55_AWES_Drone_Reel__T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  uint32_T m_bpIndex;                  /* '<S8>/ThrustX' */
  int_T Uw_IWORK;                      /* '<S13>/U,w' */
} DW_FinWing_55_AWES_Drone_Reel_T;

/* Continuous states (default storage) */
typedef struct {
  real_T theta;                        /* '<S13>/Theta' */
  real_T u[2];                         /* '<S13>/U,w' */
  real_T Xe[2];                        /* '<S13>/Position' */
  real_T q;                            /* '<S13>/q' */
} X_FinWing_55_AWES_Drone_Reel__T;

/* Periodic continuous state vector (global) */
typedef int_T PeriodicIndX_FinWing_55_AWES__T[1];
typedef real_T PeriodicRngX_FinWing_55_AWES__T[2];

/* State derivatives (default storage) */
typedef struct {
  real_T theta;                        /* '<S13>/Theta' */
  real_T u[2];                         /* '<S13>/U,w' */
  real_T Xe[2];                        /* '<S13>/Position' */
  real_T q;                            /* '<S13>/q' */
} XDot_FinWing_55_AWES_Drone_Re_T;

/* State disabled  */
typedef struct {
  boolean_T theta;                     /* '<S13>/Theta' */
  boolean_T u[2];                      /* '<S13>/U,w' */
  boolean_T Xe[2];                     /* '<S13>/Position' */
  boolean_T q;                         /* '<S13>/q' */
} XDis_FinWing_55_AWES_Drone_Re_T;

/* Invariant block signals (default storage) */
typedef struct {
  const real_T Gain[3];                /* '<S1>/Gain' */
  const real_T Gain1;                  /* '<S2>/Gain1' */
  const real_T ThrustX;                /* '<S8>/ThrustX' */
} ConstB_FinWing_55_AWES_Drone__T;

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
   *   '<S29>/Prelookup'
   *   '<S30>/Prelookup1'
   */
  real_T pooled3[10];

  /* Expression: aeroCoeff.CD
   * Referenced by: '<S29>/CD'
   */
  real_T CD_Table[10];

  /* Expression: aeroCoeff.CL
   * Referenced by: '<S29>/CL'
   */
  real_T CL_Table[10];

  /* Expression: aeroCoeff.def_vec
   * Referenced by: '<S30>/Prelookup'
   */
  real_T Prelookup_BreakpointsData[5];

  /* Expression: aeroCoeff.CD_el
   * Referenced by: '<S30>/CD_el'
   */
  real_T CD_el_Table[50];

  /* Expression: aeroCoeff.CL_el
   * Referenced by: '<S30>/CL_el'
   */
  real_T CL_el_Table[5];

  /* Expression: aeroCoeff.CM
   * Referenced by: '<S29>/CM'
   */
  real_T CM_Table[10];

  /* Expression: aeroCoeff.CM_el
   * Referenced by: '<S30>/CM_el'
   */
  real_T CM_el_Table[5];
} ConstP_FinWing_55_AWES_Drone__T;

/* External inputs (root inport signals with default storage) */
typedef struct {
  real_T ElevonPitch;                  /* '<Root>/ElevonPitch' */
} ExtU_FinWing_55_AWES_Drone_Re_T;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  real_T Altitude;                     /* '<Root>/Altitude' */
} ExtY_FinWing_55_AWES_Drone_Re_T;

/* Real-time Model Data Structure */
struct tag_RTM_FinWing_55_AWES_Drone_T {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;
  X_FinWing_55_AWES_Drone_Reel__T *contStates;
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
extern B_FinWing_55_AWES_Drone_Reel__T FinWing_55_AWES_Drone_Reel_Ou_B;

/* Continuous states (default storage) */
extern X_FinWing_55_AWES_Drone_Reel__T FinWing_55_AWES_Drone_Reel_Ou_X;

/* Block states (default storage) */
extern DW_FinWing_55_AWES_Drone_Reel_T FinWing_55_AWES_Drone_Reel_O_DW;

/* External inputs (root inport signals with default storage) */
extern ExtU_FinWing_55_AWES_Drone_Re_T FinWing_55_AWES_Drone_Reel_Ou_U;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY_FinWing_55_AWES_Drone_Re_T FinWing_55_AWES_Drone_Reel_Ou_Y;
extern const ConstB_FinWing_55_AWES_Drone__T FinWing_55_AWES_Drone_Re_ConstB;/* constant block i/o */

/* Constant parameters (default storage) */
extern const ConstP_FinWing_55_AWES_Drone__T FinWing_55_AWES_Drone_Re_ConstP;

/* Model entry point functions */
extern void FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_initialize(void);
extern void FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_step(void);
extern void FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_terminate(void);

/* Real-time Model object */
extern RT_MODEL_FinWing_55_AWES_Dron_T *const FinWing_55_AWES_Drone_Reel_O_M;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<S20>/Unit Conversion' : Unused code path elimination
 * Block '<S28>/Airspeed' : Unused code path elimination
 * Block '<S32>/Product' : Unused code path elimination
 * Block '<S32>/Product1' : Unused code path elimination
 * Block '<S32>/Sum' : Unused code path elimination
 * Block '<S3>/Clock1' : Unused code path elimination
 * Block '<S3>/Gain' : Unused code path elimination
 * Block '<S3>/Output' : Unused code path elimination
 * Block '<S3>/Product' : Unused code path elimination
 * Block '<S3>/Product1' : Unused code path elimination
 * Block '<S3>/Product2' : Unused code path elimination
 * Block '<S3>/Sum' : Unused code path elimination
 * Block '<S3>/deltaFreq' : Unused code path elimination
 * Block '<S3>/initialFreq' : Unused code path elimination
 * Block '<S3>/targetTime' : Unused code path elimination
 * Block '<Root>/Constant Steady State' : Unused code path elimination
 * Block '<Root>/Constant Steady State1' : Unused code path elimination
 * Block '<S4>/Gain1' : Unused code path elimination
 * Block '<S6>/Diff' : Unused code path elimination
 * Block '<S6>/UD' : Unused code path elimination
 * Block '<S36>/Constant' : Unused code path elimination
 * Block '<S36>/Constant1' : Unused code path elimination
 * Block '<S36>/Constant2' : Unused code path elimination
 * Block '<S38>/Reshape1' : Unused code path elimination
 * Block '<S38>/pg,qg,rg' : Unused code path elimination
 * Block '<Root>/Multiport Switch' : Unused code path elimination
 * Block '<Root>/Positive Step' : Unused code path elimination
 * Block '<Root>/Power' : Unused code path elimination
 * Block '<Root>/Pulse Generator' : Unused code path elimination
 * Block '<S9>/Clock' : Unused code path elimination
 * Block '<S9>/Constant' : Unused code path elimination
 * Block '<S9>/Look-Up Table1' : Unused code path elimination
 * Block '<S9>/Math Function' : Unused code path elimination
 * Block '<S9>/Output' : Unused code path elimination
 * Block '<S9>/Sum' : Unused code path elimination
 * Block '<S10>/Clock' : Unused code path elimination
 * Block '<S10>/Constant' : Unused code path elimination
 * Block '<S10>/Look-Up Table1' : Unused code path elimination
 * Block '<S10>/Math Function' : Unused code path elimination
 * Block '<S10>/Output' : Unused code path elimination
 * Block '<S10>/Sum' : Unused code path elimination
 * Block '<Root>/Saturation' : Unused code path elimination
 * Block '<Root>/Scope' : Unused code path elimination
 * Block '<Root>/Scope1' : Unused code path elimination
 * Block '<Root>/Scope2' : Unused code path elimination
 * Block '<Root>/SignalCmnd' : Unused code path elimination
 * Block '<Root>/Sine Trajectory' : Unused code path elimination
 * Block '<Root>/Subtract' : Unused code path elimination
 * Block '<S11>/Product2' : Unused code path elimination
 * Block '<S12>/Gain' : Unused code path elimination
 * Block '<S39>/Airspeed' : Unused code path elimination
 * Block '<S39>/Incidence' : Unused code path elimination
 * Block '<S40>/Product' : Unused code path elimination
 * Block '<S40>/Product1' : Unused code path elimination
 * Block '<S40>/Sum' : Unused code path elimination
 * Block '<S12>/Scope' : Unused code path elimination
 * Block '<S22>/Reshape (9) to [3x3] column-major' : Reshape block reduction
 * Block '<S19>/Unit Conversion' : Eliminated nontunable gain of 1
 * Block '<S23>/Reshape (9) to [3x3] column-major' : Reshape block reduction
 * Block '<S37>/Reshape1' : Reshape block reduction
 * Block '<S38>/Reshape' : Reshape block reduction
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
 * '<Root>' : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep'
 * '<S1>'   : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/3DOF EOM'
 * '<S2>'   : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Aerodynamic Model'
 * '<S3>'   : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Chirp Signal'
 * '<S4>'   : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Degrees to Radians1'
 * '<S5>'   : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Degrees to Radians2'
 * '<S6>'   : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Difference Filter'
 * '<S7>'   : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Environment Model'
 * '<S8>'   : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Propulsion'
 * '<S9>'   : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Repeating Sequence'
 * '<S10>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Repeating Sequence1'
 * '<S11>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/TetherForce'
 * '<S12>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Visualization'
 * '<S13>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/3DOF EOM/3DOF (Body Axes)'
 * '<S14>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/3DOF EOM/Rotation Angles to Direction Cosine Matrix'
 * '<S15>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/3DOF EOM/3DOF (Body Axes)/Calculate qdot'
 * '<S16>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/3DOF EOM/3DOF (Body Axes)/Determine Force,  Mass & Inertia'
 * '<S17>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/3DOF EOM/3DOF (Body Axes)/Gravity'
 * '<S18>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/3DOF EOM/3DOF (Body Axes)/Transform  to Earth Axes'
 * '<S19>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/3DOF EOM/3DOF (Body Axes)/Velocity Conversion'
 * '<S20>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/3DOF EOM/3DOF (Body Axes)/Velocity Conversion1'
 * '<S21>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/3DOF EOM/3DOF (Body Axes)/Transform  to Earth Axes/Rotation Angles to Direction Cosine Matrix'
 * '<S22>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/3DOF EOM/3DOF (Body Axes)/Transform  to Earth Axes/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S23>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/3DOF EOM/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S24>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Aerodynamic Model/Aero Forces and Moments'
 * '<S25>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Aerodynamic Model/Degrees to Radians'
 * '<S26>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Aerodynamic Model/Dynamic Pressure'
 * '<S27>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Aerodynamic Model/Elevator Forces and Moments'
 * '<S28>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Aerodynamic Model/Incidence  & Airspeed'
 * '<S29>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Aerodynamic Model/Lookup Aero Coeffs'
 * '<S30>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Aerodynamic Model/Lookup Elevator Aero Coeffs'
 * '<S31>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Aerodynamic Model/Dynamic Pressure/dot'
 * '<S32>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Aerodynamic Model/Incidence  & Airspeed/dot'
 * '<S33>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Aerodynamic Model/Lookup Aero Coeffs/Subsystem'
 * '<S34>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Aerodynamic Model/Lookup Elevator Aero Coeffs/Subsystem'
 * '<S35>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Aerodynamic Model/Lookup Elevator Aero Coeffs/Subsystem1'
 * '<S36>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Environment Model/Constant Atmosphere'
 * '<S37>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Environment Model/Constant Gravity'
 * '<S38>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Environment Model/Constant Wind'
 * '<S39>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Visualization/Incidence  & Airspeed'
 * '<S40>'  : 'FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep/Visualization/Incidence  & Airspeed/dot'
 */
#endif  /* RTW_HEADER_FinWing_55_AWES_Drone_Reel_Out_Test_V_03_1_fixedstep_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
