# vim: ft=python
import os
import collections
import glob
import itertools
import pandas
import pysam
from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.Seq import Seq


configfile: "hpv_config.yaml"
workdir: os.environ['PWD']
shell.executable('bash')

PANEL = config['panel_type']
hpv_type = config['hpv_type']
hpv_ref = config['hpv_reference']['universal'] # currently, all bams are initially mapped to alpha reference

if PANEL == 'single':
    amp_bed = config['amplicon_bed']['single'] %hpv_type
    len_bed = config['len_bed']['single'] %hpv_type
else: # universal and ffpe don't have wildcards in path string
    amp_bed = config['amplicon_bed'][PANEL]
    len_bed = config['len_bed'][PANEL]

cov_dev = config["cov_dev"] # this is used for glu coverage script
GENES = ['E1', 'E2', 'E4', 'E5', 'E6', 'E7', 'E8', 'L1', 'L2']

bed = {}
# create a dictionary of type and UNPADDED genome length (only universal will have multiple keys)
# e.g. bed['HPV16_Ref'] = 7906
len_file = open(len_bed, 'r')
for line in len_file:
    (type, length) = (line.split()[0], int(line.split()[2]))
    bed[type] = length - 400
len_file.close()

## ---- The parser will have to be customized for each run ---- ##
def parse_sampleID(filename):
    base = filename.split('/')[-1] 
    # For filenames IonXpress_046_SC075985.tmap.bam
    if base.startswith('I'):
        return base.split('_')[2].split('.')[0]
    # otherwise PAP1111_2001_IonXpress_20.tmap.bam
    else:
        return base.split('_')[0]

bamfiles = sorted(glob.glob(config['tmap_path']), key=parse_sampleID)

d = {}
for key, value in itertools.groupby(bamfiles, parse_sampleID):
    if key.startswith(config['cohort']):   # comment this out if you want to keep blanks
        d[key] = list(value)

sampleIDs = d.keys()

TARGETS =   ['reports/%s_N-%d.fasta' %(config['deliver_proj'], config['fasta_n']),
            'reports/type_summary.tsv',
            ]

# eventually work these into the config yaml
include: 'qc_Snakefile' # creates the fastqc summaries and multiqc report
TARGETS += ['reports/filtered_read_count.tsv']
include: 'annotation_Snakefile' # annotates vcf and creates snpeff multiqc report
TARGETS += ['multiqc/snpeff_report.html']
include: 'peptide_Snakefile' # creates amino acid fasta files
TARGETS += ['aa_convert/%s_aa.fasta' %config['deliver_proj']]

# These rules run on the host node and are not submitted to the cluster.
localrules: all

#--------------------------------------------------------------------------
rule all:
    input: TARGETS

def link_rule_input_files(wildcards):
    return d[wildcards.sampleID]

#--------------------------------------------------------------------------
rule link:
    input: link_rule_input_files
    output: 'bams/{sampleID}.bam'
    params:
        bam = 'temp/{sampleID}/{sampleID}.temp.bam',
        sam = 'temp/{sampleID}/{sampleID}.temp.sam',
        head = 'temp/{sampleID}/{sampleID}.header'
    run:
        if (len(input) > 1):
            shell('mkdir -p %s' %os.path.dirname(params.bam))
            shell('samtools merge {params.bam} {input}')

            # all SM fields in @RG must be identical (this is generally only for pre-2016 runs)
            # create a samfile with the fixed header
            sam = pysam.Samfile(params.bam, 'rb')
            header = sam.header
            for i in header['RG']:
                i['SM'] = wildcards.sampleID
            outfile = pysam.Samfile(params.sam, 'wh', header=header)

            # add the reads from the original merged bam
            shell('samtools view {params.bam} >> {params.sam}')

            # convert back to bam
            shell('samtools view -h -b {params.sam} > {output}')
            #shell('rm {params.bam} {params.sam}') 
        else:
            shell('cd bams; ln -s ../{input} {wildcards.sampleID}.bam && touch -h {wildcards.sampleID}.bam')

#--------------------------------------------------------------------------
rule mapq_filter:
    input: rules.link.output
    output: 'mapq_filter/{sampleID}.filtered.bam'
    threads: 2
    params: # change these to temp snakemake files eventually
        mapq = int(config["aq_filter"]),
        temp = 'temp/{sampleID}/{sampleID}.aq.temp',
        pre =  'temp/{sampleID}/{sampleID}.sort.temp'
    run:
        shell('mkdir -p %s' %os.path.dirname(params.temp))
        shell('samtools view -h -q {params.mapq} {input} | samtools view -bT {hpv_ref} -o {params.temp}')
        shell('samtools sort -o {output} -@ {threads} -T {params.pre} {params.temp}')
        shell('samtools index {output}; rm {params.temp}')

#--------------------------------------------------------------------------
rule variant_call:
    input: rules.mapq_filter.output
    output:
        'tvc/{sampleID}/TSVC_variants.vcf',
        'tvc/{sampleID}/{sampleID}.ptrim.bam'
    threads: 2
    params:
        pipe = config["vc_pipe"],
        out = ' tvc/{sampleID}',
        param = config["vc_param"],
        vc_bin = config["vc_bin"],
    run:
        shell('python {params.pipe} \
        --input-bam {input} \
        --postprocessed-bam {output[1]} \
        --primer-trim-bed {amp_bed} \
        --reference-fasta {hpv_ref} \
        --num-threads {threads} \
        --output-dir {params.out} \
        --parameters-file {params.param} \
        --bin-dir {params.vc_bin} \
        --region-bed {len_bed}')

#--------------------------------------------------------------------------
rule adjust_padding:
    input: rules.variant_call.output[0]
    output: 'tvc_vcf/{sampleID}.tvc_no_pad.vcf'
    params: temp = '{sampleID}.temp.vcf'
    run:
        vcf = open(input[0], 'r')
        outfile = open(output[0], 'w')
        need_sort = False

        for line in vcf:
            if line.startswith('#'):
                outfile.write(line)
            else:
                type = line.split()[0]
                hpv_len = bed[type]
                loc = line.split()[1]
                if int(loc) > hpv_len:
                    new_loc = int(loc) - hpv_len
                    outfile.write(line.replace(loc, str(new_loc), 1))
                    need_sort = True
                else:
                    outfile.write(line)
        vcf.close()
        outfile.close()

        if need_sort == True:
            shell('vcf-sort -c {output} > {params.temp}')
            shell('mv {params.temp} {output}')

#--------------------------------------------------------------------------
rule hpv_bam:   # removes human reads
    input: 'tvc/{sampleID}/{sampleID}.ptrim.bam'
    output: 'ptrim_hpv/{sampleID}.hpv.bam'
    run:
        shell('samtools view -h -L {len_bed} {input} | samtools view -bS -o {output}')
        shell('samtools index {output}')

#--------------------------------------------------------------------------
rule pileup:
    input: rules.hpv_bam.output
    output: 'pileup/{sampleID}.pileup'
    run:
        shell('samtools mpileup -f {hpv_ref} -l {len_bed} {input} > {output}')

#--------------------------------------------------------------------------
rule fasta:
    input:
        rules.pileup.output,
        rules.adjust_padding.output
    output:
        'fasta/{sampleID}_HPV%s.fasta' %hpv_type
    run:
        # note vcf header is 0 because it comes after we've skipped 70 lines (as opposed to being the actual line 71)
        df = pandas.read_table(input[1], skiprows=70, header=0)    
        df2 = pandas.read_table(input[0], names=['chrom', 'pos', 'nt', 'cov', 'qual1', 'qual2'], sep='\t')
        types = list(set(df['#CHROM'].tolist() + df2['chrom'].tolist()))

        # create a fasta file for each HPV type found in the sample
        for hpv in types:
            print(hpv)
            num = hpv.replace('HPV', '').replace('_Ref', '')  # convert HPV16_Ref to 16
            # pull out sequence for each HPV type
            # someday make this more efficient by reading them all ahead of time
            fa = config['hpv_reference']['single'] %(num, num) # these are padded
            fa_handle = open(fa, 'r')

            # only high risk types need fastas
            # TODO - somday make a list of the high risk types (maybe in config?)
            if os.path.isfile(fa) == False:
                continue

            seq = ''
            for record in SeqIO.parse(fa_handle, 'fasta'):
                seq = str(record.seq)
                break # this also takes just the first record, it's just longer
            fa_handle.close()

            # now start looking for SNPs and deletions as per original pipeline
            dt = df[df['#CHROM'] == hpv].copy()
            newseq = seq[:len(seq)-400] # start with the ref sequence and add SNPs
            for idx, row in dt.iterrows():
                (pos, ref, alt) = (int(row['POS']), row['REF'], row['ALT'].split(',')[0])
                # Look for SNPs       
                if len(ref) == len(alt):
                    counter = 0
                    while counter < len(ref):
                        newseq = newseq[:pos-1+counter] + alt[counter] + newseq[pos+counter:]
                        counter += 1

                # Look for deletions
                elif (len(ref) > len(alt)) and (len(alt) > 1):
                    counter = 1 # the first nt is the same as the reference, so start at +1
                    while counter < len(ref):
                        newseq = newseq[:pos-1] + '-' + newseq[pos:]
                        counter += 1

                # Skip insertions and other types of variation
                else:
                    continue

            # now check pileup to make sure there was enough coverage at each location
            dp = df2[df2['chrom'] == hpv].copy()

            # Account for zero coverage in pileup (position is missing) by adding zeros
            dp.set_index('pos', inplace=True)
            allpos = list(range(len(seq)+1))
            allpos.pop(0)
            dp = dp.reindex(allpos).fillna(0) 
            dp['cov'] = dp['cov'].astype(int)

            # take padding into account
            dp['adj_pos'] = dp.index
            dp['adj_pos'] = dp['adj_pos'].apply(lambda x: (int(x) - (len(seq)-400)) if x > (len(seq)-400) else int(x))
            x = dp.groupby('adj_pos')['cov'].sum()

            # base calls with less than min_depth are called as N
            x = x[x < config['min_read']]
            for pos, depth in x.iteritems():
                newseq = newseq[:pos-1] + 'N' + newseq[pos:]

            # output a fasta file for each type in each sample
            outfile = open('fasta/%s_HPV%s.fasta' %(wildcards.sampleID, num), 'w')
            outfile.write('>%s_HPV%s\n' %(wildcards.sampleID, num))
            outfile.write(newseq + '\n')
            outfile.close()


#--------------------------------------------------------------------------
rule fasta_cat:
    input: expand('fasta/{sampleID}_HPV%s.fasta' %hpv_type, sampleID=sampleIDs)
    output: 'reports/%s.fasta' %config['deliver_proj']
    run:
        shell('cat {input} > {output}')

#--------------------------------------------------------------------------
rule fasta_n:
    input: rules.fasta_cat.output
    output: 'reports/%s_N-%d.fasta' %(config['deliver_proj'], config['fasta_n'])
    run:
        maxN = int(config['fasta_n'])
        keep = []
        for record in SeqIO.parse(input[0], 'fasta'):
            d = collections.Counter(record.seq) # creates a dictionary w/counts for each character found
            if 'N' in d.keys():
                if (float(d['N'])/len(record.seq) * 100) < maxN:
                    keep.append(record)
            else:
                keep.append(record)

        # SeqIO.write doesn't let you set wrap width, so use FastaIO directly
        outfile = open(output[0], 'w')
        fasta_out = FastaIO.FastaWriter(outfile, wrap=None)
        fasta_out.write_file(keep)
        outfile.close()

#--------------------------------------------------------------------------
rule type_summary:
    input: expand('pileup/{sampleID}.pileup', sampleID=d.keys())
    output: 'reports/type_summary.tsv'
    run:
        stacks = []
        # iterate through each sample pileup and pull out which types were found
        for sample in input:
            df = pandas.read_table(sample, names=['chrom', 'pos', 'nt', 'cov', 'qual1', 'qual2'], sep='\t')
            df['sampleID'] = sample.split('/')[-1].split('.')[0]
            # note this does not consider padding yet!
            df = df[df['cov'] >= int(config['min_read'])]
            x = df.groupby(['sampleID', 'chrom'])['pos'].count()
            y = x.unstack()
            stacks.append(y) 

        dfs = pandas.concat(stacks).fillna(0)
        dfs.to_csv(output[0], sep='\t')


