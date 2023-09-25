#******************************************************************
# Purpose: modules for some most widely used constants

# Author: George Kuczera
# Copyright, George Kuczera, 2005-2010. All rights reserved.
# Modified by Michael Leonard Dec 2022 - Without permission

ik4=4  # SELECTED_INT_KIND(9)   # Integer kind
ik=8   # SELECTED_INT_KIND(9)   # Integer kind
rk4=4  # Real*4 kind
rk=8   # Real kind

# Machine constants
hugeRe  = .Machine$double.xmax    # largest real on machine
tinyRe  = .Machine$double.xmin    # smallest real on machine
epsRe   = .Machine$double.eps     # smallest real epsilon on machine
minExp  = .Machine$double.min.exp # smallest exponent
maxExp  = .Machine$double.max.exp # biggest exponent
hugeInt = .Machine$integer.max    # largest integer on machine

# Other constants
zero=0.0; half=0.5; one=1.0; two=2.0; four=4.0
pi = 3.141592653589793238462643383279502884197
twoPi = 6.283185307179586476925286766559005768394

  # Script constants and structures
maxScriptLines=700; nCharScriptLine=60; scriptSize=maxScriptLines*nCharScriptLine; codeSize=3*scriptSize; maxnScriptProc=40

scriptProcType=list()
scriptProcType$nCode=0
scriptProcType$name=""
scriptProcType$nameUC=""
scriptProcType$description=""
scriptProcType$formula=""
scriptProcType$code=rep(0,codeSize)

