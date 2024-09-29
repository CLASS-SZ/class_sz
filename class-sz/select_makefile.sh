#!/bin/bash

# Detect OS and architecture
UNAME_S=$(uname -s)
UNAME_M=$(uname -m)

echo "UNAME_S: $UNAME_S"
echo "UNAME_M: $UNAME_M"

# Select the appropriate Makefile
if [ "$UNAME_S" == "Darwin" ] && [ "$UNAME_M" == "arm64" ]; then
    echo "using M1 Makefile"
    cp Makefile_m1 Makefile
elif [ "$UNAME_S" == "Linux" ]; then
    echo "using Linux Makefile"
    cp Makefile_linux Makefile
else
    echo "OS not supported yet: $UNAME_S. See README.md for preM1 mac. Otherwise, please get in touch or open an issue on the CLASS-SZ github."
    echo "UNAME_S: $UNAME_S"
    echo "UNAME_M: $UNAME_M"
    echo "will try with Makefile_m1x"
    cp Makefile_m1 Makefile
fi

