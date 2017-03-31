from flask import Flask, render_template, request, jsonify, send_file
import os
from os import listdir
from os.path import isfile, join
import process

application = Flask(__name__)

########################################################################################################################

# Pre-process input
def processInput(input_records):
    CSNstrings = []
    ENSTIDs = []
    for record in input_records:
        record = record.strip()
        if record == '': continue
        cols = record.split()
        if len(cols) != 2: continue
        CSNstrings.append(cols[0])
        ENSTIDs.append(cols[1])

    if len(CSNstrings) > 500:
        CSNstrings = CSNstrings[:500]
        ENSTIDs = ENSTIDs[:500]

    return CSNstrings, ENSTIDs

# Create output directory
def makeOutputDir():
    workingdir = os.getcwd()
    i = 0
    while True:
        i += 1
        if not os.path.isdir(workingdir+'/submissions/submission_'+str(i)):
            os.mkdir(workingdir+'/submissions/submission_'+str(i))
            return workingdir+'/submissions/submission_'+str(i)

# Get dictionary of transcript database files
def getTranscriptDBs():
    ret = dict()
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    filenames = [f for f in listdir(scriptdir + '/transdbs/') if f.endswith('.txt') and isfile(join(scriptdir + '/transdbs/', f))]
    for fn in filenames:
        name = ''
        for line in open(scriptdir + '/transdbs/' + fn):
            line = line.strip()
            if line.startswith('#'):
                name = line[line.find('Ensembl release'):]
                break
        ret[name] = fn
    return ret

# Search in transcript database
def getTranscriptsForSearch(transcriptdb_fn, searchtxt):
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    ret = []
    for line in open(scriptdir+'/transdbs/'+transcriptdb_fn):
        line = line.strip()
        if line == '' or line.startswith('#') or line.startswith('ENSG\t'): continue
        cols = line.split('\t')
        gene = cols[1]
        enst = cols[2]
        if searchtxt in gene or searchtxt in enst: ret.append((gene,enst))
    return ret

# Read CSN Lookup configuration file
def readConfigFile():
    ret = dict()
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    for line in open(scriptdir+'/config.txt'):
        line = line.strip()
        if line == '' or line.startswith('#'): continue
        [key, value] = line.split('=')
        key = key.strip()
        value = value.strip()
        ret[key] = value
    return ret

########################################################################################################################

dir = ''

@application.route('/', methods=['GET', 'POST'])
def index():
    global dir

    transcriptDBs = getTranscriptDBs()
    transcriptDB_names = sorted(transcriptDBs.keys())

    if request.method == 'POST':
        input_records = request.form['inputfield'].split('\n')
        input_records = [str(x) for x in input_records]
        selected_db = request.form['selected_transcriptdb']
        transcriptdb = transcriptDBs[selected_db]

        dir = makeOutputDir()

        CSNstrings, ENSTIDs = processInput(input_records)
        if len(CSNstrings) == 0: return render_template('index.html', transDBs = transcriptDB_names)

        refconfig = readConfigFile()
        gbuild = selected_db[selected_db.find('GRC'):-1]
        ref_fn = refconfig[gbuild]

        results = process.run(dir, CSNstrings, ENSTIDs, transcriptdb, ref_fn)
        return render_template('results.html', results = results)

    return render_template('index.html', transDBs = transcriptDB_names)

@application.route('/manual', methods=['GET'])
def manual():
    return render_template('manual.html')

@application.route('/transcripts', methods=['GET', 'POST'])
def transcripts():

    transcriptDBs = getTranscriptDBs()
    transcriptDB_names = sorted(transcriptDBs.keys())

    if request.method == 'POST':
        data = request.json
        searchtxt = data['searchtxt']
        selected_db = data['selected_db']

        transcriptdb = transcriptDBs[selected_db]
        return jsonify(getTranscriptsForSearch(transcriptdb, searchtxt))

    initdbfn = transcriptDBs[transcriptDB_names[0]]
    return render_template('transcripts.html', initlist=getTranscriptsForSearch(initdbfn, 'BRCA1'), transDBs = transcriptDB_names)

@application.route('/download', methods=['POST'])
def download():
    global dir
    return send_file(dir+'/output.txt', attachment_filename='output.txt', as_attachment=True, mimetype='text/plain')

@application.route('/ie', methods=['GET'])
def ie():
    return render_template('ie.html')


if __name__ == "__main__":
    application.debug = True
    application.run()
