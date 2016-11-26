#!/bin/bash
cd output
rm result.gif

for a in [0-9]*.gif; do
    mv $a `printf %04d.%s ${a%.*} ${a##*.}`
done

convert -delay 5 -loop 0 -layers OptimizeFrame *.gif result.gif
#convert -delay 5 -loop 0 -coalesce -layers RemoveDups *.gif result.gif
