## Running bioezy from command line

1. cd to bzy.py location
2. (optional) workon venv...
3. py (python3) 
4. import bzy OR from bzy import *
5. call function like normal: eg - pattern_count('ata','ta')

## Capture input and output files

import sys
inFile = sys.argv[1] -- can take in any amount of arguments OR something like input_file.txt. NOTE sys.argv[0] = file name 
outFile = sys.argv[2]

Then you can read in all your data, do your manipulations, and write out the results:

with open(inFile,'r') as i:
    lines = i.readlines()

processedLines = manipulateData(lines)

with open(outFile,'w') as o:
    for line in processedLines:
        o.write(line)

You can call this program by running py script.py input_file.txt output_file.txt. It will run script on input_file.txt and send the stdout to output_file.txt

## Upload new version on PyPi

a. Delete all files in the dist folder.

b. Update the version number in the setup.py file.

c. Re-create the wheels: python3 setup.py sdist bdist_wheel (be in setup directory)

d. Re-upload the new files: twine upload dist/*

## Upload new version to Github format _v x.x.x u 0.0.x_

_version x.x.x update 0.0.x_

a. git status

b. git add -a (if all changes are ready to be committed)

c. git commit -m "version x.x.x update 0.0.x"

d. git pull

e. git push 

