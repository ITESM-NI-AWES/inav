#include "io/serial.h"
#include "drivers/serial.h"
#include "common/time.h"
#include <stdint.h>

static serialPort_t * serialPort = NULL;

void customSerialTest (timeUs_t currentTimeUs){
    UNUSED(currentTimeUs);
    static uint8_t serialValue = 0;
    //char buf[10] = " Hola!!\n";
    if (serialRxBytesWaiting(serialPort)) {
        serialValue = serialRead(serialPort);
        if (serialValue<254){
            serialValue++;
        } else{
            serialValue = 0;
        }
        serialWrite(serialPort,serialValue);
        //serialPrint(serialPort, buf);
    }
}
void customSerialTest_Init (void){
    serialPort = openCustomSerialPort(NULL, NULL, baudRates[BAUD_115200], MODE_RX | MODE_TX, SERIAL_NOT_INVERTED | SERIAL_STOPBITS_1 | SERIAL_PARITY_NO );
    //serialPort = openCSerialPort(portConfig->identifier, FUNCTION_NONE, NULL, NULL, baudRates[BAUD_115200], MODE_TX, SERIAL_NOT_INVERTED | SERIAL_STOPBITS_1 | SERIAL_PARITY_NO );
    //serialPort = uartOpen(USART2, NULL, NULL, baudRates[BAUD_115200], MODE_TX, SERIAL_NOT_INVERTED | SERIAL_STOPBITS_1 | SERIAL_PARITY_NO);
}