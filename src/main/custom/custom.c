#include "io/serial.h"
#include "drivers/serial.h"
#include "common/time.h"
#include "common/printf.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAXIMUM_STRING_SIZE 50
static serialPort_t * serialPort = NULL;

void customSerialTest (timeUs_t currentTimeUs){
    UNUSED(currentTimeUs);
    char inBuf[MAXIMUM_STRING_SIZE];
    char outBuf[MAXIMUM_STRING_SIZE];
    uint8_t i = 0;
    static double serialValue = (double)0.0;
    bool startReading = false;

    while (serialRxBytesWaiting(serialPort)) {
        char c = serialRead(serialPort);
        if (i >= MAXIMUM_STRING_SIZE) break;
        if (startReading){
            inBuf[i] = c;
            i++;
        }
        if(c == '<')startReading= true;
        if (c == '>')break;
    }
    if(strlen(inBuf)>9){
        serialValue = atof(inBuf);
        serialValue += (double)0.01;
        tfp_sprintf(outBuf,"%10f", serialValue);
        serialPrint(serialPort, outBuf);
    }
}
void customSerialTest_Init (void){
    serialPort = openCustomSerialPort(NULL, NULL, baudRates[BAUD_921600], MODE_RX | MODE_TX, SERIAL_NOT_INVERTED | SERIAL_STOPBITS_1);
}