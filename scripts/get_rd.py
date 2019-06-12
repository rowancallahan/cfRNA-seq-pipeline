import pandas as pd

files = snakemake.input
out = snakemake.output[0]

def split_exon_fraction(files, out):
    sample_list = []
    for file in files:
        sample = file.split('.')[0]
        results = {}
        equal = False
        for line in open(file):
            if '==' in line:
                equal = True

            if not equal:
                line = line.split()
                key ='_'.join(line[:2])
                results[key] = float(line[-1])

            if equal:
                if not '==' in line and not'Group' in line:
                    line = line.split()
                    key = line[0]
                    results[key] = float(line[2])
                else:
                    continue

        exon_fraction = (results['CDS_Exons'] + results["5'UTR_Exons"] + results["3'UTR_Exons"]) / results["Total_Assigned"]
        intron_fraction = results['Introns'] / results['Total_Assigned']
        sample_list.append([sample, exon_fraction, intron_fraction])

    frame = pd.DataFrame(sample_list)
    frame.columns = ['Sample','Exon','Intron']
    frame.to_csv(out,sep='\t',index=False)

split_exon_fraction(files, out)
