#!/bin/bash

# Detect OS and architecture
UNAME_S=$(uname -s)
UNAME_M=$(uname -m)

# Select the appropriate Makefile
if [ "$UNAME_S" == "Darwin" ] && [ "$UNAME_M" == "arm64" ]; then
    cp Makefile_m1 Makefile
elif [ "$UNAME_S" == "Linux" ]; then
    cp Makefile_linux Makefile
else
    echo "OS not supported yet: $UNAME_S. See README.md for preM1 mac. Otherwise, please get in touch or open an issue on the CLASS-SZ github."
    exit 1
fi

