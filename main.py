from Bio import SeqIO
from Bio.Seq import Seq
from flask import Flask, flash, request, redirect, url_for, render_template,send_file,make_response
import urllib.request
from werkzeug.utils import secure_filename
from os.path import join, dirname, realpath
import pandas as pd
from io import StringIO,BytesIO
from time import perf_counter
from werkzeug.wrappers import Response

app = Flask(__name__)

# UPLOAD_FOLDER = 'static/uploads/'
UPLOAD_FOLDER = join(dirname(realpath(__file__)), 'static/uploads/')
app.secret_key = "secret key"
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024
ALLOWED_EXTENSIONS = set(['fasta'])

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@app.route('/')
def home():
    return render_template('index.html')
Wild_type = ['Wild_type',
                 {69: 'H', 70: 'V', 145: 'Y', 501: 'N', 570: 'A', 614: 'D', 681: 'P', 716: 'T', 982: 'S', 1118: 'D',
                  18: 'L', 80: 'D', 215: 'D', 242: 'L', 243: 'A', 244: 'L', 246: 'R', 417: 'K', 484: 'E', 701: 'A'}]
Alpha = ['B.1.1.7_Alpha',
         {69: 'S', 70: 'T', 78: 'H', 80: 'P', 145: 'N', 501: 'G', 611: 'G', 570: 'T', 614: 'C', 678: 'H', 681: 'A',
          713: 'I', 716: 'T', 979: 'A', 982: 'D', 1115: 'H', 1118: 'F'}]
Beta = ['B.1.351_Beta',
        {18: 'F', 80: 'A', 215: 'G', 242: 'H', 243: 'I', 244: 'S', 246: 'L', 414: 'N', 481: 'K', 498: 'Y', 501: 'G',
         611: 'G', 614: 'C', 698: 'V', 701: 'S'}]
Gamma = ['B.1.1.248_P1_Gamma', {417: 'T', 484: 'K', 501: 'Y'}]
Zeta = ['P.2_Zeta', {484: 'K', 614: 'G'}]
Kappa = ['B.1.617.1_Kappa',
         {19: 'R', 142: 'D', 157: 'V', 156: 'G', 450: 'R', 476: 'T', 482: 'Q', 612: 'G', 679: 'R'}]
Delta = ['B.1.617.2_Delta',
         {19: 'R', 156: 'G', 157: 'V', 415: 'K', 450: 'R', 476: 'K', 482: 'E', 612: 'G', 679: 'R'}]
Delta_Plus = ['B.1.617.2_Delta_Plus',
              {19: 'R', 156: 'G', 157: 'V', 415: 'N', 450: 'R', 476: 'K', 482: 'E', 612: 'G', 679: 'R'}]
Delta_Plus2 = ['B.1.617.3_2_Delta_Plus2',
               {19: 'R', 156: 'G', 157: 'V', 450: 'R', 476: 'T', 482: 'Q', 612: 'G', 679: 'R'}]
Omicron = ['B.1.1.529_Omicron',
           {67: 'V', 69: 'S', 93: 'I', 140: 'D', 141: 'H', 142: 'K', 143: 'N', 209: 'E', 210: 'P', 211: 'E',
            368: 'L', 370: 'P', 372: 'F'}]
Epsilon = ['B.1.427_Epsilon', {452: 'R', 614: 'G'}]
Epsilon2 = ['B.1.429_Epsilon', {13: 'I', 152: 'C', 452: 'R', 614: 'G'}]
Iota = ['B.1.526_Iota', {5: '"', 477: 'N', 484: 'K', 701: 'V', 95: 'I', 253: 'G', 614: 'G'}]
Eta = ['B.1.525_Eta', {69: 'S', 70: 'G', 477: 'N', 484: 'K', 701: 'V', 95: 'I', 253: 'G', 614: 'G', 888: 'L'}]
Theta = ['P.3_Theta',
         {141: 'Y', 142: 'H', 143: 'K', 481: 'K', 498: 'Y', 1173: 'F', 611: 'G', 678: 'H', 1089: 'K', 1098: 'Y'}]
Delta_lineage = ['Delta_lineage', {156: 'G', 157: 'V'}]
Omicron_lineage = ['Omicron_lineage',
                   {67: 'V', 69: 'S', 93: 'I', 140: 'D', 141: 'H', 142: 'K', 143: 'N', 209: 'E', 210: 'P'}]
listofVariants = [Wild_type, Alpha, Beta, Gamma, Zeta, Kappa, Delta, Delta_Plus, Delta_Plus2, Omicron, Epsilon,
                  Epsilon2, Iota, Eta, Theta, Delta_lineage, Omicron_lineage]

def variant_from_seq(entry):
    dictofresuls = {}
    try:
        my_seq = entry.seq
        pos1 = my_seq.rfind("ATGTTTGTT")
        pos2 = my_seq.find("TCAAATTACATTACAC")
        my_prot_seq = my_seq[pos1:pos2].translate()
        # print(my_prot_seq)
        if pos1 < 0 or pos2 < 0:
            return "Not found"
        for variant in listofVariants:
            matchscore = 0
            # print(variant[0])
            for key, value in variant[1].items():
                # print(key,value)
                if my_prot_seq[key - 1] == value:  # Minus one to adjust for pythonic counting
                    matchscore += 1
            dictofresuls[variant[0]] = (round(matchscore / len(variant[1]) * 100, 2))
        seq_result_series=pd.Series(data=dictofresuls,index=dictofresuls.keys())
        # print(seq_result_series)
        return seq_result_series
    except ValueError:
        return "invalid input"

@app.route('/', methods=['POST'])
def upload_image():
    global resultsdf
    resultsdf = pd.DataFrame()
    t1_start = perf_counter()
    if 'file' not in request.files:
        flash('No file part')
        return redirect(request.url)
    filels = request.files.getlist("file")
    for file in filels:
        # print(file.filename)
        a=(file.read().decode("utf-8") )
        b=StringIO(a)
        for record in SeqIO.parse(b, "fasta"):
            # print(record.id, record.seq[0:50])
            tempseries=variant_from_seq(record)
            resultsdf[str(record.id)]=tempseries
            # print('series',tempseries)
            # resultsdf=resultsdf.copy()
    resultsdf=resultsdf.T
    print(resultsdf)
    htmltable = resultsdf.to_html()
    render_template('index.html',table=htmltable)
    t1_stop = perf_counter()
    print(f"Elapsed time: {t1_stop-t1_start} seconds")
    flash(f"Elapsed time: {t1_stop-t1_start} seconds")
    # print(type(resultsdf.to_csv()))
    return render_template('index.html', table=htmltable)
@app.route('/download')
def download():
    global resultsdf
    try:
        resp = make_response(resultsdf.to_csv())
        resp.headers["Content-Disposition"] = f"attachment; filename=results.csv"
        resp.headers["Content-Type"] = "text/csv"
        resultsdf=pd.DataFrame()
        return resp
    except:
        flash(f"Submit files first")
        return render_template('index.html')


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=8080, debug=True)
