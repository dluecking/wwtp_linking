#!/bin/bash

cat intermediate/proteins/lc/*.faa \
    intermediate/proteins/plv/*.faa \
    intermediate/proteins/vph/*.faa \
    intermediate/proteins/gv/*.faa \
    > intermediate/proteins/all_proteins.faa