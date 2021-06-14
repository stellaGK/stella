#!/bin/bash

set -ex

export GK_SYSTEM=gnu_ubuntu
make -I Makefiles -j2
