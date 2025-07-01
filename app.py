from flask import Flask, render_template, request, redirect, url_for, session
from werkzeug.utils import secure_filename
import os
from classes import FastaParser, MitochondrialDNA, SequenceMotif, SequenceAlignment
app = Flask(__name__)
app.secret_key = 'your_secret_key_here'  # Necessario per usare session

# Configurazione upload
UPLOAD_FOLDER = 'uploads'
ALLOWED_EXTENSIONS = {'fasta'}

if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        return redirect(request.url)
   
    file = request.files['file']
    if file.filename == '':
        return redirect(request.url)
   
    if file and allowed_file(file.filename):
        filename = secure_filename(file.filename)
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)
       
        # Parse del file FASTA
        parser = FastaParser()
        parser.parse_file(filepath)
       
        # Salva il DataFrame nella sessione
        session['filepath'] = filepath
       
        # Calcola le statistiche
        df = parser.get_DataFrame()
        stats = []
        for _, row in df.iterrows():
            mito_dna = MitochondrialDNA(row['Identifier'], row['Description'], row['Sequence'])
            stats.append({
                'ID': row['Identifier'],
                'Description': row['Description'],
                'Sequence': row['Sequence'],
                'Length': mito_dna.length(),
                'GC Content': f"{mito_dna.gc_content():.2f}%"
            })
       
        return render_template('stats.html', stats=stats)
   
    return redirect(request.url)

@app.route('/motif', methods=['GET', 'POST'])
def motif():
    if request.method == 'POST':
        parser = FastaParser()
        parser.parse_file(session['filepath'])
        df = parser.get_DataFrame()
        motif_analyzer = SequenceMotif(df)
       
        motif_input = request.form['motif']
       
        try:
            # Prova a interpretare l'input come "sequence, motif_length, minimum"
            seq_idx, motif_length, minimum = map(int, motif_input.split(','))
            if 0 <= seq_idx < len(df):
                sequence = df.iloc[seq_idx]['Sequence']
                results = motif_analyzer.extract_motifs(sequence, motif_length, minimum)
                return render_template('motif.html', results=results.to_string().split('\n'))
        except ValueError:
            # Se non è possibile, tratta l'input come un motif da cercare
            results = motif_analyzer.find_motif(motif_input)
            return render_template('motif.html', results=results.to_string().split('\n'))
   
    return render_template('motif.html')

@app.route('/align', methods=['GET', 'POST'])
def align():
    if request.method == 'POST':
        parser = FastaParser()
        parser.parse_file(session['filepath'])
        df = parser.get_DataFrame()
       
        try:
            seq1_idx = int(request.form['seq1'])
            seq2_idx = int(request.form['seq2'])
           
            if 0 <= seq1_idx < len(df) and 0 <= seq2_idx < len(df):
                seq1 = df.iloc[seq1_idx]['Sequence']
                seq2 = df.iloc[seq2_idx]['Sequence']
               
                alignment = SequenceAlignment(seq1, seq2)
                alignment.perform_alignment()
                result = alignment.format_alignment()
                score = alignment.alignments_score()
               
                message = f"Alignment Score: {score}"
                return render_template('align.html', message=message, alignment_result=result)
            else:
                return render_template('align.html', message="Invalid sequence indices")
        except ValueError:
            return render_template('align.html', message="Please enter valid numbers")
   
    return render_template('align.html', message="Enter sequence indices to align")

if __name__ == '__main__':
    app.run(debug=True)