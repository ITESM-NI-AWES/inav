
/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: embededcodeexample1.c
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

#include "control/embededcodeexample1.h"
#include "common/time.h"
/* Block signals and states (default storage) */
DW rtDW;

/* External inputs (root inport signals with default storage) */
ExtU rtU;

/* External outputs (root outports fed by signals with default storage) */
ExtY rtY;

/* Model step function */
void embededcodeexample1_step(timeUs_t currentTimeUs)
{
    UNUSED(currentTimeUs);
  /* Product: '<S2>/Product' incorporates:
   *  Constant: '<Root>/Scaling'
   *  Inport: '<Root>/coolant_temp'
   */
  static float increase = 0.0;
  increase ++;
  if (increase > 300) {
      increase = 0.0;
  }

  rtU.coolant_temp = increase;

  rtDW.Product = rtU.coolant_temp * rtP.Scaling_Value;

  /* Sum: '<S2>/Add' incorporates:
   *  Constant: '<Root>/Offset'
   */
  rtDW.Add = rtDW.Product - rtP.Offset_Value;

  /* Outport: '<Root>/Error' incorporates:
   *  Constant: '<S1>/MaxTemp'
   *  RelationalOperator: '<S1>/GreaterThan'
   */
  rtY.Error = (rtDW.Add > rtP.MaxTemp_Value);

  /* Switch: '<S1>/Switch1' incorporates:
   *  Outport: '<Root>/Error'
   */
  if (rtY.Error) {
    /* Outport: '<Root>/Temperature' incorporates:
     *  Constant: '<S1>/DefaultTemp'
     */
    rtY.Temperature = rtP.DefaultTemp_Value;
  } else {
    /* Outport: '<Root>/Temperature' */
    rtY.Temperature = rtDW.Add;
  }

  /* End of Switch: '<S1>/Switch1' */
}

/* Model initialize function */
void embededcodeexample1_initialize(void)
{
  /* Registration code */

  /* states (dwork) */
  (void) memset((void *) &rtDW, 0,
                sizeof(DW));

  /* external inputs */
  rtU.coolant_temp = 0.0F;

  /* external outputs */
  (void) memset((void *)&rtY, 0,
                sizeof(ExtY));
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */

