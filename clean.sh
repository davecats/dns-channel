#!/bin/bash
find . -name "*.c" -exec rm {} \;
find . -name "*.d" -exec rm {} \;
find . -name "*.o" -exec rm {} \;
find . -name "*.m~" -exec rm {} \;
