
/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: embededcodeexample1.h
 *
 * Code generated for Simulink model 'embededcodeexample1'.
 *
 * Model version                  : 1.10
 * Simulink Coder version         : 9.3 (R2020a) 18-Nov-2019
 * C/C++ source code generated on : Tue Apr 28 11:32:12 2020
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: ARM Compatible->ARM Cortex
 * Code generation objectives:
 *    1. Execution efficiency
 *    2. Traceability
 * Validation result: Not run
 */

#ifndef RTW_HEADER_embededcodeexample1_h_
#define RTW_HEADER_embededcodeexample1_h_
#include "control/rtwtypes.h"
#include <string.h>
#ifndef embededcodeexample1_COMMON_INCLUDES_
# define embededcodeexample1_COMMON_INCLUDES_
#include "control/rtwtypes.h"
#endif                                /* embededcodeexample1_COMMON_INCLUDES_ */
#include "common/time.h"
/* Model Code Variants */

/* Macros for accessing real-time model data structure */

/* Block signals and states (default storage) for system '<Root>' */
typedef struct {
  real32_T Product;                    /* '<S2>/Product' */
  real32_T Add;                        /* '<S2>/Add' */
} DW;

/* External inputs (root inport signals with default storage) */
typedef struct {
  real32_T coolant_temp;               /* '<Root>/coolant_temp' */
} ExtU;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  real32_T Temperature;                /* '<Root>/Temperature' */
  boolean_T Error;                     /* '<Root>/Error' */
} ExtY;

/* Parameters (default storage) */
struct P_ {
  real32_T DefaultTemp_Value;          /* Computed Parameter: DefaultTemp_Value
                                        * Referenced by: '<S1>/DefaultTemp'
                                        */
  real32_T Scaling_Value;              /* Computed Parameter: Scaling_Value
                                        * Referenced by: '<Root>/Scaling'
                                        */
  real32_T Offset_Value;               /* Computed Parameter: Offset_Value
                                        * Referenced by: '<Root>/Offset'
                                        */
  real32_T MaxTemp_Value;              /* Computed Parameter: MaxTemp_Value
                                        * Referenced by: '<S1>/MaxTemp'
                                        */
};

/* Parameters (default storage) */
typedef struct P_ P;

/* Block parameters (default storage) */
extern P rtP;

/* Block signals and states (default storage) */
extern DW rtDW;

/* External inputs (root inport signals with default storage) */
extern ExtU rtU;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY rtY;

/* Model entry point functions */
extern void embededcodeexample1_initialize(void);
extern void embededcodeexample1_step(timeUs_t currentTimeUs);

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<Root>/Constant' : Unused code path elimination
 * Block '<Root>/Display' : Unused code path elimination
 * Block '<Root>/Display1' : Unused code path elimination
 * Block '<Root>/Display2' : Unused code path elimination
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
 * '<Root>' : 'embededcodeexample1'
 * '<S1>'   : 'embededcodeexample1/OutOfRange'
 * '<S2>'   : 'embededcodeexample1/ScaleAndOffset'
 */

/*-
 * Requirements for '<Root>': embededcodeexample1
 */
#endif                                 /* RTW_HEADER_embededcodeexample1_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */

