#include "BCH.h"
#include <stdio.h>
#include <stdlib.h>
#include "math.h"

int main()
{
	//m=9 and t=55 -> BCH(511, 130)
	if (!BCH_Init(9, 55)) printf("Init failed");

	int k = BCH_GetWordLength();
	int n = BCH_GetCodeLength();

	uint32_t const KBytes = ceil(k / 8.0);

	uint8_t* data = calloc(KBytes, sizeof(uint8_t));

	for (int i = 0; i < k; i++)
		SetBitAtIndex(data, i, ((i * k) % 37) % 2);

	//encode word
	uint8_t* codeword = BCH_Encode(data);
	if (codeword == 0) printf("encode failed");

	printf("Codeword: ");
	for (int i = 0; i < ceil(n / 8.0); i++) {
		printf("%2x ", codeword[i]);
	}
	printf("\n");

	//Add artificial noise/errors for testing
	int numerr = 6;
	int errpos[6];
	errpos[0] = 1;
	errpos[1] = 12;
	errpos[2] = 11;
	errpos[3] = 15;
	errpos[4] = 19;
	errpos[5] = 25;

	if (numerr > 0)
	{
		for (int i = 0; i < numerr; i++)
		{
			SetBitAtIndex(codeword, errpos[i], GetBitAtIndex(codeword, errpos[i]) ^ 1);
		}
	}
	printf("Codeword with errors: ");
	for (int i = 0; i < ceil(n/8.0); i++) {
		printf("%2x ", codeword[i]);
	}
	printf("\n");

	//Decode codeword
	uint8_t* decodeddata = BCH_Decode(codeword);
	if(decodeddata == 0) printf("decode failed!");

	printf("Recovered data: ");
	for (int i = 0; i < KBytes; i++) {
		printf("%2x ", decodeddata[i]);
	}
	printf("\n");

	//Check if decoding was successful
	int decerror = 0;
	for (int i = 0; i < k; i++)
	{
		if (GetBitAtIndex(data, i) != GetBitAtIndex(decodeddata, i)) decerror++;
	}

	if (decerror)
		printf("There were %d decoding errors in message positions\n", decerror);
	else
		printf("Succesful decoding\n");

	//Cleanup - make sure to free the values returned by BCH_Encode and BCH_Decode yourself!
	BCH_Cleanup();
	free(codeword);
	free(decodeddata);
	free(data);
}
