
from annotlib import CSN, Reference, Transcript, CSNException
import gzip
import subprocess
import os

scriptdir = os.path.dirname(os.path.realpath(__file__))


########################################################################################################################

# Process variants with CAVA
def processWithCAVA(submission_ID, transcriptdb_file, CSNs, ref_fn):
    ret = []

    # Create input file for CAVA
    inputvcf = open(submission_ID+'/cavain.vcf', 'w')
    cavainput_counter = 0
    for csn in CSNs:
        if csn.error is None:
            cavainput_counter += 1
            inputvcf.write(str(csn.vcf)+'\n')
    inputvcf.close()

    # Create config file for CAVA
    createCAVAconfigFile(submission_ID, transcriptdb_file, ref_fn)

    # Create transcript list file for CAVA
    createTranscriptListFile(submission_ID, CSNs)

    # Run CAVA
    if cavainput_counter > 0:
        stdout = subprocess.check_output([scriptdir + '/cava/cava.py', '-i', submission_ID + '/cavain.vcf', '-c', submission_ID + '/config.txt', '-o', submission_ID + '/cavaout'])

        # Read in CAVA output
        cavaoutput = readCAVAoutput(submission_ID)

    # Create return list by adding CAVA results
    for csn in CSNs:

        if csn.error is None:
            key = csn.vcf.chrom + ':' + str(csn.vcf.pos) + '_' + csn.vcf.ref + '_' + csn.vcf.alts[0] + '_' + csn.transcript
            record = cavaoutput[key]
            cavacsn = record['csn']
            if 'p.' in cavacsn: cava_p = cavacsn[cavacsn.find('p.')+2:]
            else: cava_p = ''

            if cava_p != str(csn.protein):
                record['error'] = 'X p. part (%s) incorrect, does not agree with CAVA annotation (%s)' % (csn.protein, cava_p)
            elif csn.CSNstring != cavacsn:
                if csn.CSNstring == record['altcsn']: record['error'] = 'X Indel incorrectly aligned to 5\' end (should be '+cavacsn+')'
                else: record['error'] = 'Uncharacterized error'
            else:
                record['chrom'] = csn.vcf.chrom
                record['pos'] = csn.vcf.pos
                record['ref'] = csn.vcf.ref
                record['alt'] = csn.vcf.alts[0]
                record['error'] = u'\u2713'+' CSN is correct'
        else:
            record = dict()
            record['error'] = 'X '+csn.error

        record['csn_input'] = csn.CSNstring

        record['transcript_input'] = csn.transcriptID
        ret.append(record.copy())

    return ret

# Create transcript list file for CAVA
def createTranscriptListFile(submission_ID, CSNs):
    transcripts = set()
    for csn in CSNs:
        if csn.error is None: transcripts.add(csn.transcript)
    transcripts = list(transcripts)
    out = open(submission_ID + '/transcripts.txt','w')
    for t in transcripts: out.write(t+'\n')
    out.close()

# Create CAVA config file
def createCAVAconfigFile(submission_ID, transcriptdb_file, ref_fn):
    config = open(submission_ID + '/config.txt', 'w')
    config.write('@reference = '+ref_fn+'\n')
    config.write('@ensembl = '+transcriptdb_file+'\n')
    config.write('@outputformat =  TSV\n')
    config.write('@transcriptlist = '+submission_ID + '/transcripts.txt'+'\n')
    config.close()

# Read CAVA output file
def readCAVAoutput(submission_ID):
    ret = dict()
    for line in open(submission_ID + '/cavaout.txt'):
        line = line.strip()
        if line == '' or line.startswith('#') or line.startswith('ID'): continue
        cols = line.split('\t')
        key = cols[1] + ':' + cols[2] + '_' + cols[3] + '_' + cols[4] + '_' + cols[8]
        ret[key] = {'type': cols[7], 'trans': cols[8], 'gene': cols[9], 'geneid': cols[10], 'loc': cols[12], 'class': cols[17], 'csn':cols[13], 'altcsn':cols[-3]}
    return ret

# Read transcript database file
def readTranscriptFile(transcriptdb_file):
    transcriptdb = gzip.open(transcriptdb_file, 'r')
    ret = dict()
    for line in transcriptdb:
        line = line.strip()
        if line == '': continue
        cols = line.split('\t')
        ret[cols[0]] = line
    return ret

# Write input to file
def writeInputToFile(submission_ID, transcriptIDs, CSNstrings):
    out = open(submission_ID+'/input.txt', 'w')
    for i in range(len(transcriptIDs)): out.write(CSNstrings[i]+'\t'+transcriptIDs[i]+'\n')
    out.close()

# Write output to file
def writeOutputToFile(submission_ID, results):
    out = open(submission_ID + '/output.txt', 'w')
    header = ['CSN','TRANSCRIPT','ERROR','CHROM','POS','REF','ALT','TYPE','GENE','LOC','CLASS']
    out.write('\t'.join(header) + '\n')
    for res in results:
        if res['error'].endswith('CSN is correct'): errfield = '.'
        else: errfield = '\"'+res['error'][2:]+'\"'

        output = [res['csn_input'], res['transcript_input'], errfield]

        if errfield != '.':
            output += ['.']*8
        else:
            output += [res['chrom'], str(res['pos']), res['ref'], res['alt'], res['type'], res['gene'], res['loc'], res['class']]
        out.write('\t'.join(output) + '\n')
    out.close()

# Run analysis
def run(submission_ID, CSNstrings, transcriptIDs, transdb, ref_fn):
    transcriptdb_file = scriptdir + '/transdbs/' + transdb[:-3]+'gz'

    reference = Reference(ref_fn)
    transcripts = readTranscriptFile(transcriptdb_file)

    writeInputToFile(submission_ID, transcriptIDs, CSNstrings)

    CSNs = []
    for i in range(len(transcriptIDs)):
        transcriptID, CSNstring = transcriptIDs[i], CSNstrings[i]
        csn = CSN(CSNstring, transcriptID)
        try:
            csn.parseString(CSNstring)
            if transcriptID not in transcripts.keys(): raise CSNException('Transcript ID (%s) not found in transcript database' % transcriptID)
            transcript = Transcript(transcripts[transcriptID])
            csn.convert(reference, transcript)
        except CSNException as e:
            csn.error = e.message
        CSNs.append(csn)

    ret = processWithCAVA(submission_ID, transcriptdb_file, CSNs, ref_fn)
    writeOutputToFile(submission_ID, ret)
    return ret