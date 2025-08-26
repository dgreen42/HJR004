#!/bin/bash -ue
minimap2 --version | sed 's/^/minimap2,/' >> tool_versions.txt
samtools --version | head -n1 >> tool_versions.txt
R --version | head -n1 >> tool_versions.txt
