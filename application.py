from flask import Flask, render_template, request, jsonify, send_file
import os
from os import listdir
from os.path import isfile, join
import logging
import process

application = Flask(__name__)


def create_output_directory():
    """
    Make sure we have somewhere to put the output of each
    CSN lookup
    """
    submissions_dir = os.getcwd() + "/submissions"

    if not os.path.exists(submissions_dir):
        os.mkdir(submissions_dir)


def log_function_enter_and_exit(the_function):
    def wrapper(*args, **kwargs):

        application.logger.debug(
            "Entering function {}".format(the_function.func_name)
        )

        ret_val = the_function(*args, **kwargs)

        application.logger.debug(
            "Exited function {}".format(the_function.func_name)
        )

        return ret_val

    return wrapper


@log_function_enter_and_exit
def processInput(input_records):
    application.logger.debug("Processing input records {}".format(input_records))

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


@log_function_enter_and_exit
def makeOutputDir():
    workingdir = os.getcwd()
    i = 0
    while True:
        i += 1
        if not os.path.isdir(workingdir+'/submissions/submission_'+str(i)):
            os.mkdir(workingdir+'/submissions/submission_'+str(i))
            return workingdir+'/submissions/submission_'+str(i)


@log_function_enter_and_exit
def getTranscriptDBs():
    ret = {}
    script_dir = os.path.dirname(os.path.realpath(__file__))
    file_names = []

    for f in listdir(script_dir + '/transdbs/'):
        if f.endswith('.txt'):
            application.logger.debug("Adding transcript db {}".format(f))
            file_names.append(f)

    for fn in file_names:
        for line in open(script_dir + '/transdbs/' + fn):
            line = line.strip()
            if line.startswith('#'):
                name = line[line.find('Ensembl release'):]
                ret[name] = fn
                application.logger.debug(
                    "Adding ENSEMBL release {} from file {}".format(name, fn)
                )
                break
    return ret


@log_function_enter_and_exit
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


@log_function_enter_and_exit
def load_config_from_file():
    script_dir = os.path.dirname(os.path.realpath(__file__))

    for line in open(script_dir + '/config.txt'):
        line = line.strip()

        if line == '' or line.startswith('#'):
            continue

        [key, value] = line.split('=')
        key = key.strip()
        value = value.strip()
        application.config[key] = value


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

        gbuild = selected_db[selected_db.find('GRC'):-1]
        ref_fn = application.config[gbuild]

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

    application.debug = False
    application.logger.info("Setting up loggers")

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(pathname)s - Line %(lineno)s - %(message)s")
    load_config_from_file()
    create_output_directory()

    file_handler = logging.FileHandler(application.config['log_file_name'], 'wa')
    file_handler.setFormatter(formatter)
    file_handler.setLevel(logging.DEBUG)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    stream_handler.setLevel(logging.DEBUG)

    application.logger.addHandler(file_handler)
    application.logger.addHandler(stream_handler)
    application.logger.setLevel(logging.DEBUG)
    application.logger.debug("Finished setting up logger")
    application.logger.debug("Running CSNLookup server")
    application.run()
