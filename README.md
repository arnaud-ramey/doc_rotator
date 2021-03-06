                  +----------------------+
                  |      doc_rotator     |
                  +----------------------+

[![Build Status](https://travis-ci.org/arnaud-ramey/doc_rotator.svg)](https://travis-ci.org/arnaud-ramey/doc_rotator)

A simple tool for rectifying tilted scanned images.

License :                  see the LICENSE file.
Authors :                  see the AUTHORS file.
How to build the program:  see the INSTALL file.

It first detects the lines segments in the document.
These lines can be made by text characters, drawings, etc.
The method used is the Hough transform for line segments, cf the Wikipedia page about it.
http://en.wikipedia.org/wiki/Hough_transform

However, not all lines correspond to the text lines, there also are some
outliers, such as a frame around the page, or simply false detections. To
determine the predominant orientation, a RANSAC filter is applied on the
bunch of line segments (cf the RANSAC Wikipedia page for details:
http://en.wikipedia.org/wiki/RANSAC ).

The most "popular" orientation is then used to rectify the document.

This application was coded in C++ using OpenCV.
The source code is online at https://sites.google.com/site/rameyarnaud/media/books/scanned-doc-rotater .
________________________________________________________________________________

How to use the program
________________________________________________________________________________
To display the help, just launch the program in a terminal.
It will display the help of the program.

Note: this can also be done with ImageMagick:
$ mkdir deskewed --parents ; for F in `ls *.jpeg`; do echo $F ; convert $F -deskew 40% "deskewed/$F" ; done
:)
