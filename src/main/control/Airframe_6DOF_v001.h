/*
 * Airframe_6DOF_v001.h
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

#ifndef RTW_HEADER_Airframe_6DOF_v001_h_
#define RTW_HEADER_Airframe_6DOF_v001_h_
#include <string.h>
#include <math.h>
#include <float.h>
#ifndef Airframe_6DOF_v001_COMMON_INCLUDES_
# define Airframe_6DOF_v001_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#endif                                 /* Airframe_6DOF_v001_COMMON_INCLUDES_ */

#include "Airframe_6DOF_v001_types.h"

/* Shared type includes */
#include "multiword_types.h"
#include "rtGetNaN.h"
#include "rt_nonfinite.h"
#include "rt_defines.h"
#include "rtGetInf.h"
#include "common/time.h"
#include "flight/imu.h"
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
  real_T Sum1;                         /* '<S52>/Sum1' */
  real_T Product[4];                   /* '<S71>/Product' */
  real_T SFunction_o1;                 /* '<S50>/S-Function' */
  real_T SFunction_o2;                 /* '<S50>/S-Function' */
  real_T SFunction_o3;                 /* '<S50>/S-Function' */
  real_T SFunction_o4;                 /* '<S50>/S-Function' */
  real_T Product2[3];                  /* '<S8>/Product2' */
  real_T TmpSignalConversionAtphithetaps[3];/* '<S7>/phidot thetadot psidot' */
  real_T Sum[3];                       /* '<S6>/Sum' */
  real_T Product_b[3];                 /* '<S14>/Product' */
  real_T w[2];                         /* '<S77>/w' */
  real_T LwgV1[2];                     /* '<S77>/Lwg//V 1' */
  real_T w_c[2];                       /* '<S77>/w ' */
  real_T w_n[2];                       /* '<S76>/w' */
  real_T w_d[2];                       /* '<S76>/w ' */
  real_T w1[2];                        /* '<S76>/w 1' */
  real_T w_a[2];                       /* '<S75>/w' */
  real_T w1_a[2];                      /* '<S75>/w1' */
  real_T w_l[2];                       /* '<S74>/w' */
  real_T UnaryMinus[2];                /* '<S74>/Unary Minus' */
  real_T w_cz[2];                      /* '<S73>/w' */
  real_T sigma_w[2];                   /* '<S72>/sigma_w' */
  real_T w_d5[2];                      /* '<S72>/w' */
} B_Airframe_6DOF_v001_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  real_T NextOutput[4];                /* '<S71>/White Noise' */
  real_T SFunction_temp_table[8];      /* '<S50>/S-Function' */
  real_T SFunction_pres_table[8];      /* '<S50>/S-Function' */
  real_T Product2_DWORK4[9];           /* '<S8>/Product2' */
  uint32_T PreLookUpIndexSearchprobofexcee;
                         /* '<S78>/PreLook-Up Index Search  (prob of exceed)' */
  uint32_T m_bpIndex;                  /* '<S4>/ThrustX' */
  uint32_T PreLookUpIndexSearchaltitude_DW;
                               /* '<S78>/PreLook-Up Index Search  (altitude)' */
  uint32_T RandSeed[4];                /* '<S71>/White Noise' */
  int8_T ifHeightMaxlowaltitudeelseifHei;
  /* '<S67>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */
  int8_T ifHeightMaxlowaltitudeelseifH_k;
  /* '<S66>/if Height < Max low altitude  elseif Height > Min isotropic altitude ' */
  boolean_T Hwgws_MODE;                /* '<S62>/Hwgw(s)' */
  boolean_T Hvgws_MODE;                /* '<S62>/Hvgw(s)' */
  boolean_T Hugws_MODE;                /* '<S62>/Hugw(s)' */
  boolean_T Hrgw_MODE;                 /* '<S61>/Hrgw' */
  boolean_T Hqgw_MODE;                 /* '<S61>/Hqgw' */
  boolean_T Hpgw_MODE;                 /* '<S61>/Hpgw' */
} DW_Airframe_6DOF_v001_T;

/* Continuous states (default storage) */
typedef struct {
  real_T Xe[3];                        /* '<S6>/xe,ye,ze' */
  real_T phi[3];                       /* '<S7>/phi theta psi' */
  real_T Ubody[3];                     /* '<S6>/ub,vb,wb' */
  real_T p[3];                         /* '<S6>/p,q,r ' */
  real_T wg_p1_CSTATE[2];              /* '<S77>/wg_p1' */
  real_T wg_p2_CSTATE[2];              /* '<S77>/wg_p2' */
  real_T vg_p1_CSTATE[2];              /* '<S76>/vg_p1' */
  real_T vgw_p2_CSTATE[2];             /* '<S76>/vgw_p2' */
  real_T ug_p_CSTATE[2];               /* '<S75>/ug_p' */
  real_T rgw_p_CSTATE[2];              /* '<S74>/rgw_p' */
  real_T qgw_p_CSTATE[2];              /* '<S73>/qgw_p' */
  real_T pgw_p_CSTATE[2];              /* '<S72>/pgw_p' */
} X_Airframe_6DOF_v001_T;

/* Periodic continuous state vector (global) */
typedef int_T PeriodicIndX_Airframe_6DOF_v0_T[3];
typedef real_T PeriodicRngX_Airframe_6DOF_v0_T[6];

/* State derivatives (default storage) */
typedef struct {
  real_T Xe[3];                        /* '<S6>/xe,ye,ze' */
  real_T phi[3];                       /* '<S7>/phi theta psi' */
  real_T Ubody[3];                     /* '<S6>/ub,vb,wb' */
  real_T p[3];                         /* '<S6>/p,q,r ' */
  real_T wg_p1_CSTATE[2];              /* '<S77>/wg_p1' */
  real_T wg_p2_CSTATE[2];              /* '<S77>/wg_p2' */
  real_T vg_p1_CSTATE[2];              /* '<S76>/vg_p1' */
  real_T vgw_p2_CSTATE[2];             /* '<S76>/vgw_p2' */
  real_T ug_p_CSTATE[2];               /* '<S75>/ug_p' */
  real_T rgw_p_CSTATE[2];              /* '<S74>/rgw_p' */
  real_T qgw_p_CSTATE[2];              /* '<S73>/qgw_p' */
  real_T pgw_p_CSTATE[2];              /* '<S72>/pgw_p' */
} XDot_Airframe_6DOF_v001_T;

/* State disabled  */
typedef struct {
  boolean_T Xe[3];                     /* '<S6>/xe,ye,ze' */
  boolean_T phi[3];                    /* '<S7>/phi theta psi' */
  boolean_T Ubody[3];                  /* '<S6>/ub,vb,wb' */
  boolean_T p[3];                      /* '<S6>/p,q,r ' */
  boolean_T wg_p1_CSTATE[2];           /* '<S77>/wg_p1' */
  boolean_T wg_p2_CSTATE[2];           /* '<S77>/wg_p2' */
  boolean_T vg_p1_CSTATE[2];           /* '<S76>/vg_p1' */
  boolean_T vgw_p2_CSTATE[2];          /* '<S76>/vgw_p2' */
  boolean_T ug_p_CSTATE[2];            /* '<S75>/ug_p' */
  boolean_T rgw_p_CSTATE[2];           /* '<S74>/rgw_p' */
  boolean_T qgw_p_CSTATE[2];           /* '<S73>/qgw_p' */
  boolean_T pgw_p_CSTATE[2];           /* '<S72>/pgw_p' */
} XDis_Airframe_6DOF_v001_T;

/* Invariant block signals (default storage) */
typedef struct {
  const real_T Selector[9];            /* '<S8>/Selector' */
  const real_T Selector2[9];           /* '<S8>/Selector2' */
  const real_T Product3;               /* '<S33>/Product3' */
  const real_T Sum[3];                 /* '<S27>/Sum' */
  const real_T UnitConversion_a;       /* '<S98>/Unit Conversion' */
  const real_T ThrustX;                /* '<S4>/ThrustX' */
} ConstB_Airframe_6DOF_v001_T;

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
  /* Expression: h_vec
   * Referenced by: '<S78>/PreLook-Up Index Search  (altitude)'
   */
  real_T PreLookUpIndexSearchaltitude_Br[12];

  /* Expression: sigma_vec'
   * Referenced by: '<S78>/Medium//High Altitude Intensity'
   */
  real_T MediumHighAltitudeIntensity_Tab[84];

  /* Computed Parameter: MediumHighAltitudeIntensity_max
   * Referenced by: '<S78>/Medium//High Altitude Intensity'
   */
  uint32_T MediumHighAltitudeIntensity_max[2];
} ConstP_Airframe_6DOF_v001_T;

/* External inputs (root inport signals with default storage) */
typedef struct {
  real_T Aileron;                      /* '<Root>/AileronCmd' */
  real_T Elevator;                     /* '<Root>/ElevatorCmd' */
  real_T Rudder;                       /* '<Root>/RudderCmd' */
} ExtU_Airframe_6DOF_v001_T;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  real_T StatesOut[18];                /* '<Root>/StatesOut' */
} ExtY_Airframe_6DOF_v001_T;

/* Real-time Model Data Structure */
struct tag_RTM_Airframe_6DOF_v001_T {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;
  X_Airframe_6DOF_v001_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[28];
  real_T odeF[3][28];
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
    struct {
      uint8_T TID[3];
    } TaskCounters;

    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[3];
  } Timing;
};

/* Block signals (default storage) */
extern B_Airframe_6DOF_v001_T Airframe_6DOF_v001_B;

/* Continuous states (default storage) */
extern X_Airframe_6DOF_v001_T Airframe_6DOF_v001_X;

/* Block states (default storage) */
extern DW_Airframe_6DOF_v001_T Airframe_6DOF_v001_DW;

/* External inputs (root inport signals with default storage) */
extern ExtU_Airframe_6DOF_v001_T Airframe_6DOF_v001_U;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY_Airframe_6DOF_v001_T Airframe_6DOF_v001_Y;
extern const ConstB_Airframe_6DOF_v001_T Airframe_6DOF_v001_ConstB;/* constant block i/o */

/* Constant parameters (default storage) */
extern const ConstP_Airframe_6DOF_v001_T Airframe_6DOF_v001_ConstP;

/* Model entry point functions */
extern void Airframe_6DOF_v001_initialize(void);
extern void Airframe_6DOF_v001_step(timeUs_t currentTimeUs);
extern void Airframe_6DOF_v001_terminate(void);

/* Real-time Model object */
extern RT_MODEL_Airframe_6DOF_v001_T *const Airframe_6DOF_v001_M;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<S12>/Unit Conversion' : Unused code path elimination
 * Block '<S13>/Unit Conversion' : Unused code path elimination
 * Block '<S57>/Unit Conversion' : Unused code path elimination
 * Block '<S58>/Unit Conversion' : Unused code path elimination
 * Block '<S59>/Unit Conversion' : Unused code path elimination
 * Block '<S53>/Incidence' : Unused code path elimination
 * Block '<S123>/Airspeed' : Unused code path elimination
 * Block '<S123>/Incidence' : Unused code path elimination
 * Block '<S123>/Product' : Unused code path elimination
 * Block '<S123>/Sideslip' : Unused code path elimination
 * Block '<S126>/Product' : Unused code path elimination
 * Block '<S126>/Product1' : Unused code path elimination
 * Block '<S126>/Product2' : Unused code path elimination
 * Block '<S126>/Sum' : Unused code path elimination
 * Block '<S5>/State Plots' : Unused code path elimination
 * Block '<S17>/Reshape (9) to [3x3] column-major' : Reshape block reduction
 * Block '<S19>/Reshape1' : Reshape block reduction
 * Block '<S19>/Reshape2' : Reshape block reduction
 * Block '<S20>/Reshape1' : Reshape block reduction
 * Block '<S20>/Reshape2' : Reshape block reduction
 * Block '<S8>/Reshape' : Reshape block reduction
 * Block '<S8>/Reshape1' : Reshape block reduction
 * Block '<S11>/Unit Conversion' : Eliminated nontunable gain of 1
 * Block '<S14>/Reshape1' : Reshape block reduction
 * Block '<S14>/Reshape2' : Reshape block reduction
 * Block '<S27>/coefAdjust' : Eliminated nontunable gain of 1
 * Block '<S55>/Unit Conversion' : Eliminated nontunable gain of 1
 * Block '<S56>/Unit Conversion' : Eliminated nontunable gain of 1
 * Block '<S84>/Reshape' : Reshape block reduction
 * Block '<S84>/Reshape1' : Reshape block reduction
 * Block '<S86>/Reshape' : Reshape block reduction
 * Block '<S92>/Reshape' : Reshape block reduction
 * Block '<S92>/Reshape1' : Reshape block reduction
 * Block '<S94>/Reshape' : Reshape block reduction
 * Block '<S122>/Reshape (9) to [3x3] column-major' : Reshape block reduction
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
 * '<Root>' : 'Airframe_6DOF_v001'
 * '<S1>'   : 'Airframe_6DOF_v001/6DOF EOM'
 * '<S2>'   : 'Airframe_6DOF_v001/Aerodynamic Model'
 * '<S3>'   : 'Airframe_6DOF_v001/Environment Model'
 * '<S4>'   : 'Airframe_6DOF_v001/Propulsion'
 * '<S5>'   : 'Airframe_6DOF_v001/Visualization'
 * '<S6>'   : 'Airframe_6DOF_v001/6DOF EOM/6DOF (Euler Angles)'
 * '<S7>'   : 'Airframe_6DOF_v001/6DOF EOM/6DOF (Euler Angles)/Calculate DCM & Euler Angles'
 * '<S8>'   : 'Airframe_6DOF_v001/6DOF EOM/6DOF (Euler Angles)/Calculate omega_dot'
 * '<S9>'   : 'Airframe_6DOF_v001/6DOF EOM/6DOF (Euler Angles)/Determine Force,  Mass & Inertia'
 * '<S10>'  : 'Airframe_6DOF_v001/6DOF EOM/6DOF (Euler Angles)/Vbxw'
 * '<S11>'  : 'Airframe_6DOF_v001/6DOF EOM/6DOF (Euler Angles)/Velocity Conversion'
 * '<S12>'  : 'Airframe_6DOF_v001/6DOF EOM/6DOF (Euler Angles)/Velocity Conversion1'
 * '<S13>'  : 'Airframe_6DOF_v001/6DOF EOM/6DOF (Euler Angles)/Velocity Conversion2'
 * '<S14>'  : 'Airframe_6DOF_v001/6DOF EOM/6DOF (Euler Angles)/transform to Inertial axes '
 * '<S15>'  : 'Airframe_6DOF_v001/6DOF EOM/6DOF (Euler Angles)/Calculate DCM & Euler Angles/Rotation Angles to Direction Cosine Matrix'
 * '<S16>'  : 'Airframe_6DOF_v001/6DOF EOM/6DOF (Euler Angles)/Calculate DCM & Euler Angles/phidot thetadot psidot'
 * '<S17>'  : 'Airframe_6DOF_v001/6DOF EOM/6DOF (Euler Angles)/Calculate DCM & Euler Angles/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S18>'  : 'Airframe_6DOF_v001/6DOF EOM/6DOF (Euler Angles)/Calculate omega_dot/3x3 Cross Product'
 * '<S19>'  : 'Airframe_6DOF_v001/6DOF EOM/6DOF (Euler Angles)/Calculate omega_dot/I x w'
 * '<S20>'  : 'Airframe_6DOF_v001/6DOF EOM/6DOF (Euler Angles)/Calculate omega_dot/I x w1'
 * '<S21>'  : 'Airframe_6DOF_v001/6DOF EOM/6DOF (Euler Angles)/Calculate omega_dot/3x3 Cross Product/Subsystem'
 * '<S22>'  : 'Airframe_6DOF_v001/6DOF EOM/6DOF (Euler Angles)/Calculate omega_dot/3x3 Cross Product/Subsystem1'
 * '<S23>'  : 'Airframe_6DOF_v001/6DOF EOM/6DOF (Euler Angles)/Vbxw/Subsystem'
 * '<S24>'  : 'Airframe_6DOF_v001/6DOF EOM/6DOF (Euler Angles)/Vbxw/Subsystem1'
 * '<S25>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Actuator Aerodynamic Coefficients'
 * '<S26>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Aerodynamic Coefficients'
 * '<S27>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Aerodynamic Forces and Moments'
 * '<S28>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Dynamic Pressure'
 * '<S29>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Incidence, Sideslip, & Airspeed'
 * '<S30>'  : 'Airframe_6DOF_v001/Aerodynamic Model/ToBus'
 * '<S31>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Actuator Aerodynamic Coefficients/Aileron'
 * '<S32>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Actuator Aerodynamic Coefficients/Elevator'
 * '<S33>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Actuator Aerodynamic Coefficients/Flap'
 * '<S34>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Actuator Aerodynamic Coefficients/Rudder'
 * '<S35>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Aerodynamic Coefficients/Body Rate Damping'
 * '<S36>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Aerodynamic Coefficients/Datum Coefficients'
 * '<S37>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Aerodynamic Coefficients/Body Rate Damping/p'
 * '<S38>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Aerodynamic Coefficients/Body Rate Damping/q'
 * '<S39>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Aerodynamic Coefficients/Body Rate Damping/r'
 * '<S40>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Aerodynamic Forces and Moments/3x3 Cross Product'
 * '<S41>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Aerodynamic Forces and Moments/CG-CP Transformation'
 * '<S42>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Aerodynamic Forces and Moments/Force Transformation'
 * '<S43>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Aerodynamic Forces and Moments/Moment Transformation'
 * '<S44>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Aerodynamic Forces and Moments/3x3 Cross Product/Subsystem'
 * '<S45>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Aerodynamic Forces and Moments/3x3 Cross Product/Subsystem1'
 * '<S46>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Dynamic Pressure/dot'
 * '<S47>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Incidence, Sideslip, & Airspeed/Subsystem'
 * '<S48>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Incidence, Sideslip, & Airspeed/Subsystem1'
 * '<S49>'  : 'Airframe_6DOF_v001/Aerodynamic Model/Incidence, Sideslip, & Airspeed/dot'
 * '<S50>'  : 'Airframe_6DOF_v001/Environment Model/COESA Atmosphere Model'
 * '<S51>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))'
 * '<S52>'  : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA'
 * '<S53>'  : 'Airframe_6DOF_v001/Environment Model/Incidence  & Airspeed'
 * '<S54>'  : 'Airframe_6DOF_v001/Environment Model/Rotation Angles to Direction Cosine Matrix'
 * '<S55>'  : 'Airframe_6DOF_v001/Environment Model/COESA Atmosphere Model/Density Conversion'
 * '<S56>'  : 'Airframe_6DOF_v001/Environment Model/COESA Atmosphere Model/Length Conversion'
 * '<S57>'  : 'Airframe_6DOF_v001/Environment Model/COESA Atmosphere Model/Pressure Conversion'
 * '<S58>'  : 'Airframe_6DOF_v001/Environment Model/COESA Atmosphere Model/Temperature Conversion'
 * '<S59>'  : 'Airframe_6DOF_v001/Environment Model/COESA Atmosphere Model/Velocity Conversion'
 * '<S60>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Angle Conversion'
 * '<S61>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Filters on angular rates'
 * '<S62>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Filters on velocities'
 * '<S63>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Length Conversion'
 * '<S64>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Length Conversion1'
 * '<S65>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/RMS turbulence  intensities'
 * '<S66>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Select angular rates'
 * '<S67>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Select velocities'
 * '<S68>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Turbulence scale lengths'
 * '<S69>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Velocity Conversion'
 * '<S70>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Velocity Conversion2'
 * '<S71>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/White Noise'
 * '<S72>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Filters on angular rates/Hpgw'
 * '<S73>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Filters on angular rates/Hqgw'
 * '<S74>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Filters on angular rates/Hrgw'
 * '<S75>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Filters on velocities/Hugw(s)'
 * '<S76>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Filters on velocities/Hvgw(s)'
 * '<S77>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Filters on velocities/Hwgw(s)'
 * '<S78>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/RMS turbulence  intensities/High Altitude Intensity'
 * '<S79>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/RMS turbulence  intensities/Low Altitude Intensity'
 * '<S80>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Select angular rates/Interpolate  rates'
 * '<S81>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Select angular rates/Low altitude  rates'
 * '<S82>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Select angular rates/Medium//High  altitude rates'
 * '<S83>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Select angular rates/Merge Subsystems'
 * '<S84>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Select angular rates/Interpolate  rates/wind to body transformation'
 * '<S85>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Select angular rates/Interpolate  rates/wind to body transformation/convert to earth coords'
 * '<S86>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Select angular rates/Low altitude  rates/wind to body transformation'
 * '<S87>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Select angular rates/Low altitude  rates/wind to body transformation/convert to earth coords'
 * '<S88>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Select velocities/Interpolate  velocities'
 * '<S89>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Select velocities/Low altitude  velocities'
 * '<S90>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Select velocities/Medium//High  altitude velocities'
 * '<S91>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Select velocities/Merge Subsystems'
 * '<S92>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Select velocities/Interpolate  velocities/wind to body transformation'
 * '<S93>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Select velocities/Interpolate  velocities/wind to body transformation/convert to earth coords'
 * '<S94>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Select velocities/Low altitude  velocities/wind to body transformation'
 * '<S95>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Select velocities/Low altitude  velocities/wind to body transformation/convert to earth coords'
 * '<S96>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Turbulence scale lengths/Low altitude scale length'
 * '<S97>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Turbulence scale lengths/Medium//High altitude scale length'
 * '<S98>'  : 'Airframe_6DOF_v001/Environment Model/Dryden Wind Turbulence Model  (Continuous (+q +r))/Turbulence scale lengths/Medium//High altitude scale length/Length Conversion'
 * '<S99>'  : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA/LatLong wrap'
 * '<S100>' : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA/LatLong wrap1'
 * '<S101>' : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA/LongLat_offset'
 * '<S102>' : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA/pos_deg'
 * '<S103>' : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA/LatLong wrap/Latitude Wrap 90'
 * '<S104>' : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA/LatLong wrap/Wrap Longitude'
 * '<S105>' : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA/LatLong wrap/Latitude Wrap 90/Compare To Constant'
 * '<S106>' : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA/LatLong wrap/Latitude Wrap 90/Wrap Angle 180'
 * '<S107>' : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA/LatLong wrap/Latitude Wrap 90/Wrap Angle 180/Compare To Constant'
 * '<S108>' : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA/LatLong wrap/Wrap Longitude/Compare To Constant'
 * '<S109>' : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA/LatLong wrap1/Latitude Wrap 90'
 * '<S110>' : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA/LatLong wrap1/Wrap Longitude'
 * '<S111>' : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA/LatLong wrap1/Latitude Wrap 90/Compare To Constant'
 * '<S112>' : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA/LatLong wrap1/Latitude Wrap 90/Wrap Angle 180'
 * '<S113>' : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA/LatLong wrap1/Latitude Wrap 90/Wrap Angle 180/Compare To Constant'
 * '<S114>' : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA/LatLong wrap1/Wrap Longitude/Compare To Constant'
 * '<S115>' : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA/LongLat_offset/Find Radian//Distance'
 * '<S116>' : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA/LongLat_offset/rotation_rad'
 * '<S117>' : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA/LongLat_offset/Find Radian//Distance/Angle Conversion2'
 * '<S118>' : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA/LongLat_offset/Find Radian//Distance/denom'
 * '<S119>' : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA/LongLat_offset/Find Radian//Distance/e'
 * '<S120>' : 'Airframe_6DOF_v001/Environment Model/Flat Earth to LLA/LongLat_offset/Find Radian//Distance/e^4'
 * '<S121>' : 'Airframe_6DOF_v001/Environment Model/Incidence  & Airspeed/dot'
 * '<S122>' : 'Airframe_6DOF_v001/Environment Model/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S123>' : 'Airframe_6DOF_v001/Visualization/Incidence, Sideslip, & Airspeed'
 * '<S124>' : 'Airframe_6DOF_v001/Visualization/Incidence, Sideslip, & Airspeed/Subsystem'
 * '<S125>' : 'Airframe_6DOF_v001/Visualization/Incidence, Sideslip, & Airspeed/Subsystem1'
 * '<S126>' : 'Airframe_6DOF_v001/Visualization/Incidence, Sideslip, & Airspeed/dot'
 */
#endif                                 /* RTW_HEADER_Airframe_6DOF_v001_h_ */
