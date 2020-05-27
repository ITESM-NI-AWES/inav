#include "io/serial.h"
#include "common/time.h"

static serialPort_t * serialPort = NULL;
static serialPortConfig_t * portConfig;

void customSerialTest (timeUs_t currentTimeUs){
    UNUSED(currentTimeUs);
    char buf[10] = " Hola!!\n";
    if (serialIsConnected(serialPort)) {
            //serialWrite(serialPort,2);
            serialPrint(serialPort, buf);
        }
}
void customSerialTest_Init (void){
    serialPort = openCustomSerialPort(NULL, NULL, baudRates[BAUD_115200], MODE_TX, SERIAL_NOT_INVERTED | SERIAL_STOPBITS_1 | SERIAL_PARITY_NO );
    //serialPort = openCSerialPort(portConfig->identifier, FUNCTION_NONE, NULL, NULL, baudRates[BAUD_115200], MODE_TX, SERIAL_NOT_INVERTED | SERIAL_STOPBITS_1 | SERIAL_PARITY_NO );
    //serialPort = uartOpen(USART2, NULL, NULL, baudRates[BAUD_115200], MODE_TX, SERIAL_NOT_INVERTED | SERIAL_STOPBITS_1 | SERIAL_PARITY_NO);
}