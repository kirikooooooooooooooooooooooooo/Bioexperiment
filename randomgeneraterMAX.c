#include <stdio.h>
#include <stdlib.h>

int GetRandom(int min, int max)
{
    return min + (int)(rand() * (max - min + 1.0) / (1.0 + RAND_MAX));
}

int main(void)
{
 
}