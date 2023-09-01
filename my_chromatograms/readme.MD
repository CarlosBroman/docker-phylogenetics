# Sanger Viewer

This code allows for reading chromatograms from ABI files using Python.

The sequences are processed using Biopython
The plot are generated with matplotlib
Rendered in a webpage with JavaScript (work in progress)

## Initial idea:
https://biopython.org/wiki/ABI_traces

## Based on:
https://github.com/ponnhide/sangerseq_viewer/tree/main

## Using the plotly template
https://plotly.com/python/range-slider/


## Also took information from Ugene source code:
"The provided code snippet defines the implementation of the DNAChromatogram class within the UGENE bioinformatics software. This class represents DNA chromatogram data and provides methods for accessing chromatogram values, performing equality comparisons between chromatograms, and storing related information such as trace lengths and sequence lengths."
https://github.com/ugeneunipro/ugene/blob/master/src/corelibs/U2Core/src/datatype/DNAChromatogram.cpp


### For more information about the ABI file visit:
https://projects.nfstc.org/workshops/resources/articles/ABIF_File_Format.pdf

https://github.com/biopython/biopython/blob/master/Bio/SeqIO/AbiIO.py