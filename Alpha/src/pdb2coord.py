#!/usr/bin/python3
#pdb2coord.py
import re, sys

#VdWRadius in angstroms
radii=dict(
    {"C":1.5,
    "F":1.2,
    "H":0.4,
    "N":1.10,
    "O":1.05,
    "S":1.6})
def stripPDB(includeRadius=False):
    filename=str(sys.argv[1])
    fil=open(filename,'rt')
    endpos=fil.seek(0, 2)
    fil.seek(0)
    outbuff=[]
    count=0

    #Cycle through the file
    while (fil.tell()+2) < endpos:
        line=fil.readline()
        ret=cleanLine(line, includeRad=includeRadius)
        if ret is not None:
            outbuff.append(ret)
            count+=1
        
    fil.close()
    #If an output file is specified
    if "-o" in sys.argv:
        for i in range(len(outbuff)):
            outbuff[i]+='\n'

        fil=open(sys.argv[sys.argv.index("-o")+1],'wt')
        fil.write("3\n"+str(count)+"\n")
        fil.writelines(outbuff)
    else: #No output file
        print("3\n"+str(count)+"\n")
        for line in outbuff:
            print(line)

        

def cleanLine(line, includeRad=False):
    if "ATOM" in line:
        #print(line)
        newline=line[32:55:1]
        #print(newline)
        if includeRad:
            pt=line[77:78].strip(" ");
            if pt in radii:
                newline = newline+" {0}".format(radii[pt])
            else:
                raise IOError("Error: Unknown Element type found: {}".format(pt))
        #print(newline)
        elif "-s" in sys.argv:
            pt=line[77:78].strip(" ");
            print(newline+" {0}".format(radii[pt]))
        return newline
    else:
        return None

        
if __name__=="__main__":
    if ("-r" in sys.argv):
        #print("including radii")
        #print(radii)
        stripPDB(includeRadius=True);
    elif len(sys.argv)>=2:
        stripPDB()
    else:
        print(
        """Usage: python pdb2coord.py <inputfile> <opts>
      Options:
           -r              Output radii to the primary output
           -o  <filename>  File to direct primary output. Defaults to stdout
           -s              Output radii to stdout, but not to primary output.
        """)

