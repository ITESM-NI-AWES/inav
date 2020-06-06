/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: controlfinalproject_controller_only.h
 *
 * Code generated for Simulink model 'controlfinalproject_controller_only'.
 *
 * Model version                  : 1.10
 * Simulink Coder version         : 9.3 (R2020a) 18-Nov-2019
 * C/C++ source code generated on : Fri Jun  5 22:51:29 2020
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_controlfinalproject_controller_only_h_
#define RTW_HEADER_controlfinalproject_controller_only_h_
#ifndef controlfinalproject_controller_only_COMMON_INCLUDES_
# define controlfinalproject_controller_only_COMMON_INCLUDES_
#include "rtwtypes.h"
#endif                /* controlfinalproject_controller_only_COMMON_INCLUDES_ */

#include "controlfinalproject_controller_only_types.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

/* Block states (default storage) for system '<Root>' */
typedef struct {
  real_T DiscreteAutopilotForwardDiffere;
                          /* '<Root>/Discrete Autopilot Forward Differencing' */
} DW_controlfinalproject_contro_T;

/* External inputs (root inport signals with default storage) */
typedef struct {
  real_T ErrorSignal;                  /* '<Root>/Error Signal' */
} ExtU_controlfinalproject_cont_T;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  real_T ManipulatedVariable;          /* '<Root>/Manipulated Variable' */
} ExtY_controlfinalproject_cont_T;

/* Real-time Model Data Structure */
struct tag_RTM_controlfinalproject_c_T {
  const char_T * volatile errorStatus;
};

/* Block states (default storage) */
extern DW_controlfinalproject_contro_T controlfinalproject_controll_DW;

/* External inputs (root inport signals with default storage) */
extern ExtU_controlfinalproject_cont_T controlfinalproject_controlle_U;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY_controlfinalproject_cont_T controlfinalproject_controlle_Y;

/* Model entry point functions */
extern void controlfinalproject_controller_only_initialize(void);
extern void controlfinalproject_controller_only_step(void);
extern void controlfinalproject_controller_only_terminate(void);

/* Real-time Model object */
extern RT_MODEL_controlfinalproject__T *const controlfinalproject_controll_M;

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
 * '<Root>' : 'controlfinalproject_controller_only'
 */
#endif                   /* RTW_HEADER_controlfinalproject_controller_only_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
