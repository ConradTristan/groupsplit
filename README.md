# groupsplit
Repository for the results of the Algorithm Engineering course at Friedrich Schiller University Jena.

groupsplit is a program which takes PPM image files as an input and outputs grayscale images which have grouped together continuous color segments (e.g. text, background) into a single color.

Usage:

groupsplit -f inputfile -o outputfile \[optional args\]

Required arguments:

-f inputfile

Specifies the input file. Must be a file in PPM format.

-o outputfile

Specifies the file to write the output to.

Optional arguments:

-s blockSize

Defines the size of the blocks the program splits the input image into. Defaults to 16.

-g groupSize

Defines the amount of groups each block can assign pixels to. Defaults to 4.

-t threshold

Defines the sum threshold (R+G+B) in which pixels are assigned to the same group. Defaults to 50.

-c maxColors

Defines the amount of output colors the program maps the input to. Defaults to 2 (binarization).
