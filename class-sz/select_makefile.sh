#!/bin/bash

# Detect OS and architecture
UNAME_S=$(uname -s)
UNAME_M=$(uname -m)

echo "UNAME_S: $UNAME_S"
echo "UNAME_M: $UNAME_M"

# Select the appropriate Makefile
if [ "$UNAME_S" == "Darwin" ] && [ "$UNAME_M" == "arm64" ]; then
    echo "Using M1 Makefile"
    cp Makefile_m1 Makefile
elif [ "$UNAME_S" == "Linux" ]; then
    echo "Using Linux Makefile"
    cp Makefile_linux Makefile
    
    # Print environment variables
    echo "Environment Variables:"
    echo "C_INCLUDE_PATH: $C_INCLUDE_PATH"
    echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
    echo "LIBRARY_PATH: $LIBRARY_PATH"
    echo "CFLAGS: $CFLAGS"
    echo "LDFLAGS: $LDFLAGS"
elif [[ "$UNAME_S" == *"MINGW"* || "$UNAME_S" == *"CYGWIN"* || "$UNAME_S" == *"MSYS"* ]]; then
    echo "Using Windows Makefile"
    cp Makefile_windows Makefile
    
    # Optionally, print some environment variables if needed for Windows
    echo "Environment Variables:"
    echo "C_INCLUDE_PATH: $C_INCLUDE_PATH"
    echo "LIBRARY_PATH: $LIBRARY_PATH"
    echo "CFLAGS: $CFLAGS"
    echo "LDFLAGS: $LDFLAGS"
else
    echo "OS not supported yet: $UNAME_S. See README.md for pre-M1 mac. Otherwise, please get in touch or open an issue on the CLASS-SZ github."
    echo "UNAME_S: $UNAME_S"
    echo "UNAME_M: $UNAME_M"
    echo "Will try with Makefile_m1x"
    cp Makefile_m1 Makefile
fi
