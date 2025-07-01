'''flask in cui provo le route motif divise'''
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
def home():
    # Retrieve the variable from the session
    uploaded_file = session.get('filepath', None)
    return render_template('home.html', uploaded_file=uploaded_file)
    

@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        return redirect(url_for('home'))
   
    file = request.files['file']
    if file.filename == '':
        return redirect(url_for('home'))
   
    if file and allowed_file(file.filename):
        filename = secure_filename(file.filename)
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)
       
        # Parse del file FASTA
        parser = FastaParser()
        parser.parse_file(filepath)
        
        
       
        # Salva il DataFrame nella sessione
        session['filepath'] = filepath
        return redirect(url_for('home'))
#non si può render_template in una route per un'altra route , si può solo essere redirectati  
    
@app.route('/dataframe', methods=['GET', 'POST'])
def dataframe():
    if 'filepath' not in session:
        return redirect(url_for('home'))
    
    parser = FastaParser()
    parser.parse_file(session['filepath'])
    df = parser.get_DataFrame()
    data = []
    for i in range(len(df)):
        row=parser.get_row(i)
        mito_dna=MitochondrialDNA(*row)
        mito_attributes=mito_dna.get_attributes_value()
        identifier, description, sequence=mito_attributes
        data.append({ 'ID': identifier,
                'Description': description,
                'Sequence': sequence})

    return render_template('dataframe.html', data=data)


@app.route('/stats', methods=['GET', 'POST'])
def stats():
    parser = FastaParser()
    parser.parse_file(session['filepath'])
    # Calcola le statistiche
    df = parser.get_DataFrame()
    
    stats = []
    
    for i in range(len(df)):
        row=parser.get_row(i)
        mito_dna=MitochondrialDNA(*row)
        mito_attributes=mito_dna.get_attributes_value()
        identifier, description, sequence=mito_attributes
        stats.append({ 'ID': identifier,
                'Description': description,
                'Length': mito_dna.length(),
                'GC Content': f"{mito_dna.gc_content():.2f}%"
        })
        submit_indexes=request.form.get('submit_indexes')
        
        subseq=None
        
        if request.method == 'POST':
            if submit_indexes:
                 seq_idx, start, end=map(int, submit_indexes.split(','))
                 mito_DNA= MitochondrialDNA(*df.iloc[seq_idx])
                 subseq=mito_DNA.extract_subseq_by_indexing(start, end)
            
    return render_template('stats.html', stats=stats, subseq=subseq)
    #return render_template('stats.html')
#return redirect(request.url)

@app.route('/motif', methods=['GET', 'POST'])
def motif():
    search_motif=None
    find_motif=None
    
    if request.method == 'POST':
        parser = FastaParser()
        parser.parse_file(session['filepath'])
        df = parser.get_DataFrame()
        motif_analyzer = SequenceMotif(df)
           
        
        motif_input_search = request.form.get("submit_search")
        motif_input_find = request.form.get('submit_find')
        
        search_motif = []
        if motif_input_search:
            if motif_input_search:
                seq_idx, motif_length, minimum = map(int, motif_input_search.split(','))
                extraction= motif_analyzer.extract_motifs(seq_idx, motif_length, minimum)
            for _, row in extraction.iterrows():
                search_motif.append({
            'Motif': row['Motif'],
            'Indexes': row['Indexes'],
            'Repetitions': len(row['Indexes']),
        })
                
            
        find_motif = []
        if motif_input_find:
            results = motif_analyzer.find_motif(motif_input_find)
            find_motif = results.to_string().split('\n')
        
        return render_template('motif.html', search_motif=search_motif, find_motif=find_motif)

    return render_template('motif.html')
    


@app.route('/align', methods=['GET', 'POST'])
def align():
    error=None
    if request.method == 'POST':
        parser = FastaParser()
        parser.parse_file(session['filepath'])
        df = parser.get_DataFrame()
       
        try:
            seq1_idx = int(request.form['seq1'])
            seq2_idx = int(request.form['seq2'])
            
            if seq1_idx==seq2_idx:
                message='You are trying to align the same sequence, the alignment score is 100%:'
                return render_template('align.html', message=message)
           
            elif 0 <= seq1_idx < len(df) and 0 <= seq2_idx < len(df):
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