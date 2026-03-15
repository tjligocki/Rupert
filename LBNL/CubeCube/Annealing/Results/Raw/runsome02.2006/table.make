#!/bin/csh -f
awk -f table.make1.awk < runsome.out \
| sort +0n -1 +1n -2                 \
| awk -f table.make2.awk > table
