/*
 * Airframe_6DOF_v001_data.c
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

#include "Airframe_6DOF_v001.h"
#include "Airframe_6DOF_v001_private.h"

/* Invariant block signals (default storage) */
const ConstB_Airframe_6DOF_v001_T Airframe_6DOF_v001_ConstB = {
  { 5788.0, 0.0, -117.6, 0.0, 6928.9, 0.0, 117.6, 0.0, 11578.3 },/* '<S8>/Selector' */

  { 5788.0, 0.0, -117.6, 0.0, 6928.9, 0.0, 117.6, 0.0, 11578.3 },/* '<S8>/Selector2' */
  0.0,                                 /* '<S33>/Product3' */

  { 0.0, 0.0, 0.0 },                   /* '<S27>/Sum' */
  1749.9999999999998,                  /* '<S98>/Unit Conversion' */
  2300.0                               /* '<S4>/ThrustX' */
};

/* Constant parameters (default storage) */
const ConstP_Airframe_6DOF_v001_T Airframe_6DOF_v001_ConstP = {
  /* Expression: h_vec
   * Referenced by: '<S78>/PreLook-Up Index Search  (altitude)'
   */
  { 500.0, 1750.0, 3750.0, 7500.0, 15000.0, 25000.0, 35000.0, 45000.0, 55000.0,
    65000.0, 75000.0, 80000.0 },

  /* Expression: sigma_vec'
   * Referenced by: '<S78>/Medium//High Altitude Intensity'
   */
  { 3.2, 2.2, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.2, 3.6, 3.3,
    1.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.6, 6.9, 7.4, 6.7, 4.6, 2.7,
    0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 8.6, 9.6, 10.6, 10.1, 8.0, 6.6, 5.0, 4.2, 2.7,
    0.0, 0.0, 0.0, 11.8, 13.0, 16.0, 15.1, 11.6, 9.7, 8.1, 8.2, 7.9, 4.9, 3.2,
    2.1, 15.6, 17.6, 23.0, 23.6, 22.1, 20.0, 16.0, 15.1, 12.1, 7.9, 6.2, 5.1,
    18.7, 21.5, 28.4, 30.2, 30.7, 31.0, 25.2, 23.1, 17.5, 10.7, 8.4, 7.2 },

  /* Computed Parameter: MediumHighAltitudeIntensity_max
   * Referenced by: '<S78>/Medium//High Altitude Intensity'
   */
  { 11U, 6U }
};
