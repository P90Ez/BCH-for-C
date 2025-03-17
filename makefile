CC = gcc
CFLAGS = -shared
OUTPUT_DIR = build
OUTPUT = $(OUTPUT_DIR)/BCH.dll
SRC = BCH.c
HEADERS = BCH.h

all: $(OUTPUT)

$(OUTPUT): $(SRC) $(HEADERS) | $(OUTPUT_DIR)
	$(CC) $(CFLAGS) -o $(OUTPUT) $(SRC)

$(OUTPUT_DIR):
	mkdir $(OUTPUT_DIR)
