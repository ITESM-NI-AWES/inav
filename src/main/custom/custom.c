#include "io/serial.h"
#include "drivers/serial.h"
#include "drivers/time.h"
#include "common/time.h"
#include "common/printf.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "controlfinalproject_controller_only.h"

#define MAXIMUM_STRING_SIZE 50

typedef enum {
    STATE_CS_PAYLOAD = 0,
    STATE_CS_SYNC,
    STATE_CS_FINISHED
} customSerialState_e;

static serialPort_t * serialPort = NULL;
bool scheduleIsMet = false;
uint8_t stepCounter = 1;
static volatile bool serialIsUsed = false;

void customSerialComm(void){
    char inBuf[MAXIMUM_STRING_SIZE];
    char outBuf[MAXIMUM_STRING_SIZE];
    uint8_t i = 0;
    customSerialState_e state = STATE_CS_SYNC;

    while (serialRxBytesWaiting(serialPort)) {
        char c = serialRead(serialPort);
        if (state != STATE_CS_FINISHED){
            if (i >= MAXIMUM_STRING_SIZE) break;
            switch(c){
                case '>':
                    if (state == STATE_CS_PAYLOAD){
                        state = STATE_CS_FINISHED;
                        controlfinalproject_controlle_U.ErrorSignal = atof(inBuf);
                        controlfinalproject_controller_only_step();
                        tfp_sprintf(outBuf,"%20f", controlfinalproject_controlle_Y.ManipulatedVariable);
                        serialPrint(serialPort, outBuf);
                    }
                    break;
                case '<':
                    state = STATE_CS_PAYLOAD;
                    break;
                default:
                    if (state == STATE_CS_PAYLOAD){
                        inBuf[i] = c;
                        i++;
                    }
                    break;
            }
        }
        if (state == STATE_CS_FINISHED)break;
    }
    
}


void customSerialTest (timeUs_t currentTimeUs){
    UNUSED(currentTimeUs);
    
    if (!serialIsUsed){
        serialIsUsed = true;
        customSerialComm();
        serialIsUsed = false;
    }
    
}
void customSerialTest_Init (void){
    //static customSerialFrameData_t customSerialFrameData;
    serialPort = openCustomSerialPort(NULL, NULL, baudRates[BAUD_921600], MODE_RX | MODE_TX, SERIAL_NOT_INVERTED | SERIAL_STOPBITS_1);
}