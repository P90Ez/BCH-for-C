# BCH for C

This implementation is a modified version of Robert Morelos-Zaragoza's encoder/decoder for binary BCH codes in C (Version 3.1). You can find the original implementation [here](https://www.eccpage.com/bch3.c).

For copyright and licensing information, please have a look at the link provided above.

## Disclaimer

The goal of my work was to modify the existing implementation to be usable on a microcontroler (STM32F401RE in my case). I optimized the memory usage / allocations as much as possible, but I did not dare to touch the underlying algorithms. 

The author is not responsible for any malfunctioning of this program, nor for any damage caused by it.

## Usage

> Note: You can find an example in `main.c`.

```c

// 1. Initialize the module (m: Galois field order, t: Number of correctable errors). See header for more information
BCH_Init(m, t);

// 2. encode your k-bit long word. (You can retrieve the length of k via BCH_GetWordLength())
uint8_t* CodeWord = BCH_Encode(Word);

// 3. transmit codeword / add artificial noise / etc.

// 4. decode you n-bit long codeword.
uint8_t* Word = BCH_Decode(CodeWord);

// 5. make sure to free all memory
BCH_Cleanup();
free(CodeWord);
free(Word);
```
