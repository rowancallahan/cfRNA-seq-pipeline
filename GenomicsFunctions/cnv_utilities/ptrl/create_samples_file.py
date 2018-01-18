#!/usr/bin/env python

### Create the samples file that is used by the PTRL CNV Rscript.
### File should look like this:
### DNA	Sample	SampleClass
### PG2-58-013	PG2.58.013.normal	Normal
### Usage: python create_samples_file.py <file_list> <outfile>

import sys

def main():
    
    handle_filelist = open(sys.argv[1], 'rU')
    handle_out = open(sys.argv[2], 'w')

    ### Write header.
    handle_out.write('\t'.join(["DNA","Sample", "SampleClass"]))
    handle_out.write('\n')

    with handle_filelist as filelist:
        for entry in filelist:
            dna = entry.split('.')[:-1][0] 
            sample = filter(str.isalnum, entry)
            new_sample = sample + ".normal"
            sampleclass = "Normal"
            handle_out.write('\t'.join([dna, new_sample, sampleclass]))
            handle_out.write('\n')

    handle_out.close()

if __name__ == "__main__":
    main()
