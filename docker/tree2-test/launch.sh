#!/bin/bash
set -e
$(boot2docker shellinit 2> /dev/null)
docker run -v /Users/rich/Documents/Projects/veg/tree2:/src -it dockertest/tree2-test $*
