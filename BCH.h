// This implementation is a modified version of Robert Morelos-Zaragoza's encoder/decoder for binary BCH codes in C (Version 3.1)
// https://www.eccpage.com/bch3.c
// For copyright and licensing information, please have a look at the link provided above.
// The author is not responsible for any malfunctioning of this program, nor for any damage caused by it.

#ifndef BCH_INCLUDED
#define BCH_INCLUDED
#include <stdint.h>
#include <stdbool.h>

/// <summary>
/// Initializes the BCH module for a (n, k, t) BCH code, where n = (2**m) + 1 and k is calculated via t.
/// </summary>
/// <param name="m">Galois field order (m).</param>
/// <param name="t">Number of correctable errors (t).</param>
/// <returns>True on success, false otherwise.</returns>
bool BCH_Init(uint8_t const m, uint8_t const t);

/// <summary>
/// Cleanup. Frees all internally used memory.
/// </summary>
void BCH_Cleanup();

/// <summary>
/// Encodes a word with word length k to a code with length n.
/// </summary>
/// <param name="Word">Word to encode. Must be k-bits long.</param>
/// <returns>Encoded codeword with length n-bits on success. NullPtr otherwise.</returns>
uint8_t* BCH_Encode(uint8_t const* const Word);

/// <summary>
/// Decodes a code with length n to a word with length k.
/// </summary>
/// <param name="CodeWord">Code with length n-bits.</param>
/// <returns>Decoded word with length k-bits on success. NullPtr otherwise.</returns>
uint8_t* BCH_Decode(uint8_t const* const CodeWord);

/// <summary>
/// Getter for code length (n).
/// </summary>
/// <returns>Code length (n).</returns>
int BCH_GetCodeLength();

/// <summary>
/// Getter for word length (k).
/// </summary>
/// <returns>Word length (k).</returns>
int BCH_GetWordLength();

//Helper Functions
uint8_t GetBitAtIndex(uint8_t const* const Buffer, uint32_t const Index);
void SetBitAtIndex(uint8_t* const Buffer, uint32_t const Index, uint8_t const BitValue);

#endif