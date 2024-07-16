#!/bin/bash

set -ex

cmake . -B build -DSTELLA_ENABLE_TESTS=on
cmake --build build -j1 --target check
