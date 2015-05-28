#!/bin/bash
set -e
$(boot2docker shellinit 2> /dev/null)
docker run -v /Users/rich/Documents/Projects/veg/plant:/src -it dockertest/plant-test $*
