//Adds free function to call from other languages
#include<stdlib.h>

void FreeMem(void* Buffer)
{
    free(Buffer);
}