#ifndef ENDIANNESS_H
#define ENDIANNESS_H

int endian(void); // return 0 if little endian
void SwapBytesInt(void *pv);
void SwapBytesFloat(void *pv);
void SwapBytesDouble(void *pv);

#endif
