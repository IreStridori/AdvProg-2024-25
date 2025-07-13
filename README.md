# Advanced Programming project 2024/2025

## **1. Introduction**
The project focuses on the development of a software which is capable of reading and structuring a FASTA file and perform different analysis on mitochondrial DNA sequences from multiple species, while providing to the user an interacting, user-friendly, Web Interface. The project also relies on scalability, so that the possibility to expand it in the future is concrete.
 
### **Technologies used**
- Programming Language: Python (Pandas, Biopython and Flask libraries used)
- Back-end: Flask
- Front-end: HTML

---

## **2. System design**
In this repository you will find these files:
- `CRCcards.pdf`    <- Class-Responsibility-Collaborator cards
- `README.md` 	<- project documentation
- `UMLdiagram.pdf`	<- UML diagram
- `classes.py`  <- Python file with the algorithms and database handling
- `app.py`      <- Flask backend
- 5 html files       <- HTML templates for the web interface
- `file.fasta`    <- Example FASTA file for testing

### **Web Interface structure**
The things you should have in your folder when running the Web interface 

	--| app.py
	--| classes.py
	--| uploads/
 	--| _pycache_/
  	--| file.fasta
	--| templates/
		--| home.html
		--| dataframe.html
	 	--| stats.html
	  	--| motif.html
	   	--| align.html
     

---

## 3. Implementation: classes used and their responsibilities
### 3.1 Class FileParser (Abstract Superclass)
- **Purpose**: Provides an reusable structure for parsing different types of files. The only method it has is `parse_file(self, file)`, an abstract method for reading and processing data from a file path. The superclass ensures a that all subclasses implement the parse method.  

### 3.2 Class FastaParser (Concrete Subclass of DataParser)  
- **Purpose**: Implements specific parsing for FASTA files and organizes the files's genomic data in a structured table with an identifier, a description and an associated sequence. It implements the `parse_file(self, file)` method from FileParser.

- **Specific responsibilities**:
	- Reads a FASTA file line by line.
 	- Organizes the different file's objects into a list of lists with three keys:  
		- **ID**: The sequence identifier (e.g., _NC_012920.1_).  
		- **Description**: The description of the organism (e.g., _Homo sapiens mitochondrion, complete genome_).
		- **Sequence**: The nucleotide sequence itself (e.g., _GATCACAGGT..._).  
	- Converts the sequences into a Pandas DataFrame for better organization and management and returns it
 	- Gives out summary of the dataset via a specific method `describe`
  	- Returns a specific sequence of the DataFrame and its details 

- **Error management**: The class uses an `ensure_data_loaded` decorator which verifies if a file was uploaded by the user. In case it was not, it raises an error, informs the user with a specific message and halts the session. This decorator is used in two methods which require the uploaded data to properly function.

### 3.3 Class GenomicEntity (Concrete Superclass)
- **Purpose**: Models a generic genomic sequence by setting common attributes (identifier, description, sequence) and provides a common method for sequence length calculation.

### 3.4 Class MitochondrialDNA (Subclass of GenomicEntity) 
   - **Purpose**: Extends `GenomicEntity` to add specific functionalities for DNA sequences.  
   - **Specific responsibilities**:
     	- Extract a subsequence by taking as input parameters a start and an end index
     	- GC content analysis of a sequence

### 3.5 Class MotifAnalyser (Abstract Class)
- **Purpose**: Provides an reusable structure for working with different types of data and inserting different types of inputs (e.g., I could give as motif a tuple, a list, a string...). Defines the structure to work with (`data`) and an abstract method `find_motif(self, motif)` for serching a motif in a certain sequence. The superclass ensures a that all subclasses implement the find_motif method.

### 3.6 Class SequenceMotif (Concrete Subclass of MotifAnalyser)  
   - **Purpose**: extends `MotifAnalyser` and is responsible for identifying recurrent genetic patterns that may be biologically significant. It allows the use to insert as input a precise or not subsequence and it extracts it. It also provides insights into their frequency and occurrence across multiple sequences.  
   - **Methods**: 2 different ways of analysing a motif were implemented.
     - `extract_motifs(self, seq_idx, motif_length, minimum)`: if the user does not have a specific motif to search, this method extracts all possible subsequences of length `motif_length` from the given `seq_idx` which repeats at least `minimum` times. Returns a Pandas DataFrame with motif counts.
       
     - `find_motif(self, motif)`: if the user has a specific generic motif, this searches the subsequence across all sequences in the dataset. Returns a Pandas DataFrame listing:
		- **Motif** searched 
 		- Sequence **Identifier**   
		- **Number of occurrences** per sequence  

### 3.7 Sequence Alignment (Concrete Class)
- **Purpose**: using Biopython’s `PairwiseAligner` module, it manages to handle pairwise sequence alignment. It compares two sequences and calculates their similarity, producing a structured alignment output.  
- **Methods and responsibilities**:
	- Initializes the alignment object with **two input sequences (`seq1`, `seq2`)**.  
	- Uses **global alignment mode** by default.  
	- Executes **pairwise sequence alignment**.  
	- Stores the alignment results.  
	- `format_alignment(self, n=1)` formats and returns the **top `n` alignments** in a readable string format.  
	- Helps visualize how sequences match or differ.  
	- Returns the **alignment score**, a numerical measure of sequence similarity.

---

## 4. Expected inputs and outputs for classes

### FileParser 
`parse_file(self, file)`
- **Input**: file path.
- **Output**: None. Abstract method implemented by subclasses to parse file data.

### FastaParser
`ensure_data_loaded(func)` (decorator)
- **Input**: function/method that requires data to be loaded before execution.
- **Output**: decorated version of the input function that checks `self._data` before executing. If `self._data` is `None`, it raises a ValueError with a message. If data is loaded, it proceeds to execute the original function with all its arguments.

`parse_file(self, file)`
- **Input**: file path
- **Output**: None. It stores the data in the instance's `self._data` protected attribute as a DataFrame.

`get_DataFrame(self)`
- **Input**: no inputs.
- **Output**: Pandas DataFrame stored in `self._data` 

`get_summary(self)`
- **Input**: no inputs.
- **Output**: summary DataFrame (using pandas' `describe()` method). Raises an error if no data is loaded.

`get_row(self, index)`
- **Input**: integer index representing a row in the DataFrame.
- **Output**: list containing Identifier, Description, and Sequence for the specified row index. Raises an error if the index is out of bounds.

```python
parser = FastaParser()
parsed_file=parser.parse_file('file.fasta')
df=parser.get_DataFrame()
summary_file=parser.get_summary()
sequence=parser.get_row(1)
```

### GenomicEntity
`get_attributes_value(self)`
- **Input**: no inputs.
- **Output**: tuple of identifier, description, sequence.

`length(self)`
- **Input**: no inputs.
- **Output**: integer length of the sequence.

### MitochondrialDNA
`gc_content(self)`
- **Input**: no inputs.
- **Output**: float number representing GC content of the sequence.

`extract_subseq_by_indexing(self, start, end)`
- **Input**: two integers representing the start and end indices.
- **Output**: substring of the sequence from start to end included.

```python
dna=MitochondrialDNA(*sequence)
dna_attributes=dna.get_attributes_value()
dna_gc=dna.gc_content()
seq_idx=1
DNA=MitochondrialDNA(*df.iloc[seq_idx])
```

### MotifAnalyser
`find_motif(self, motif)`
- **Input**: string of a motif.
- **Output**: None. Abstract method implemented by subclasses.

### SequenceMotif
`extract_motifs(self, seq_idx, motif_length, minimum)`
- **Input**: integer zero-based index of the sequence to analyze, integer length of motifs to extract, integer treshold of occurrences for a motif to be included.
- **Output**: pandas DataFrame with motifs, indices and occurrences.

`find_motif(self, motif)`
- **Input**: string motif to search for.
- **Output**: pandas DataFrame with columns Identifier, Motif and Occurrences, showing how many times the motif appears in each sequence in the dataset.

```python
motif_analyzer = SequenceMotif(df)
seq_idx=1
motif_length=4
minimum=6
extraction=motif_analyzer.extract_motifs(seq_idx, motif_length, minimum)
motif_input_find='ATGGT'
found_motif = motif_analyzer.find_motif(motif_input_find)
```

### SequenceAlignment
`perform_alignment(self)`
- **Input**: no inputs.
- **Output**: returns the alignment results from the PairwiseAligner module.

`format_alignment(self, n=1)`
- **Input**: optional integer specifying how many alignments to format.
- **Output**: string representation of the first n alignment results.

`alignments_score(self)`
- **Input**: no inputs.
- **Output**: float score representing the quality of the alignment.

```python
seq1='ATGGT'
seq2='ATGT'
alignment = SequenceAlignment(seq1, seq2)
alignment_performed=alignment.perform_alignment()
result = alignment.format_alignment()
score = alignment.alignments_score()
```
---

## 5 Back-end: Flask framework functionalities
### Setup and imports
```python
from flask import Flask, render_template, request, redirect, url_for, session
from werkzeug.utils import secure_filename
import os
from classes import FastaParser, MitochondrialDNA, SequenceMotif, SequenceAlignment
```
- Imports necessary Flask components for web functionality
- Uses werkzeug.utils for secure file handling
- Imports custom Python classes that handle the biological aspects: (FastaParser, MitochondrialDNA, SequenceMotif, SequenceAlignment)

### App creation
```python
app = Flask(__name__)
app.secret_key = 'your_secret_key_here'  

# Upload configuration
UPLOAD_FOLDER = 'uploads'
ALLOWED_EXTENSIONS = {'fasta'}

if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
```

- The secret key is needed for session security
- Configures an upload folder for FASTA files
- Restricts uploads to only .fasta files
- Creates the upload directory if it doesn't exist

### Routes:
**Home route**:
```python
@app.route('/')
def home():
    # Retrieve the variable from the session
    uploaded_file = session.get('filepath', None)
    return render_template('home.html', uploaded_file=uploaded_file)
```
- Checks if a file has been previously uploaded by retrieving the filepath from the session
- Passes this information to the home html template

**File Upload Route**:
```python
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
       
        # Parse FASTA file
        parser = FastaParser()
        parser.parse_file(filepath)
       
        # Save filepath in session
        session['filepath'] = filepath
        return redirect(url_for('home'))
```
- Validates the file submission (presence and file type), if something is wrong it goes back to the homepage.
- Secures the filename to prevent security issues and saves it to the precreated upload directory
- Parses the file immediately using the FastaParser class
- Stores the filepath in the session for later access
- Redirects back to the home page

**DataFrame Route**
```python
@app.route('/dataframe', methods=['GET', 'POST'])
def dataframe():
    if 'filepath' not in session:
        return redirect(url_for('home'))
    
    parser = FastaParser()
    parser.parse_file(session['filepath'])
    df = parser.get_DataFrame()
    data = []
    for _, row in df.iterrows():
        mito_dna = MitochondrialDNA(row['Identifier'], row['Description'], row['Sequence'])
        data.append({
            'ID': row['Identifier'],
            'Description': row['Description'],
            'Sequence': row['Sequence']
        })
    return render_template('dataframe.html', data=data)
```
- Redirects to home if no file has been uploaded (no filepath in session)
- Parses the uploaded file to generate the DataFrame
- Iterates through each row, creating MitochondrialDNA objects
- Constructs a list of dictionaries with sequence information: ID, description, sequence
- Passes this data to the dataframe html template for displaying 

**Statistics Route**
```python
@app.route('/stats', methods=['GET', 'POST'])
def stats():
    parser = FastaParser()
    parser.parse_file(session['filepath'])
    # Calculate statistics
    df = parser.get_DataFrame()
    stats = []
    for _, row in df.iterrows():
        mito_dna = MitochondrialDNA(row['Identifier'], row['Description'], row['Sequence'])
        stats.append({
            'ID': row['Identifier'],
            'Description': row['Description'],
            'Length': mito_dna.length(),
            'GC Content': f"{mito_dna.gc_content():.2f}%"
        })
    return render_template('stats.html', stats=stats)
```
- Uses the dataframe and the MitochondrialDNA class to calculate: sequence length and GC content (with two decimal places)
- Passes the statistics to the stats template for display

**Motif Analysis Route**
```python
@app.route('/motif', methods=['GET', 'POST'])
def motif():
    search=None
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
            seq_idx, motif = motif_input_search.split(',')
            seq = df.iloc[int(seq_idx)]['Sequence']
            search = motif_analyzer.search_motif(seq, motif)
            
        find_motif = []
        if motif_input_find:
            results = motif_analyzer.find_motif(motif_input_find)
            find_motif = results.to_string().split('\n')
        
        return render_template('motif.html', search_motif=search, find_motif=find_motif)

    return render_template('motif.html')
```
Uses the SequenceMotif class to perform two different types of analysis depending on the type of form data:
- `Search` operation:
  - Parses the input to get sequence index and motif
  - Retrieves the sequence and performs the search
- `Find` operation:
  - Searches all sequences for the specified motif

Formats the results for display

**Sequence Alignment Route**
```python
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
```
- Once received the sequences from the form, implements error handling (if there is an error returns a clear messagge") for:
	- Invalid numeric inputs
	- Out-of-range sequence indices
- Uses the SequenceAlignment class to:
	- Perform the alignment algorithm
	- Format the alignment for visualization
	- Calculate the alignment score
- Returns the alignment result and score to the template

**Application Execution**
```python
if __name__ == '__main__':
    app.run(debug=True)
```
- Runs the Flask application when the script is executed directly and enables debug mode
  
---

## 6 Front-end: HTML strucutre and functionalities
We focused on a functional simple design rather than an aesthetic one given the scientific/research purposes. The navigation structure is:
**Upload → Parse → Visualize → Analyze (through statistics, motif searches, or alignments)**

### home.html 
It is the starting page for the application. Welcomes the user, explains what the tool does and enables the loading of a FASTA file.
- **Upload**:
  - Contains a form with a file input that accepts only .fasta files
  - Uses `POST` method to send the file to the "/upload" endpoint
  - Includes `multipart/form-data` encoding necessary for file uploads
  - Checks if a file has been uploaded and displays a link to visualize the parsed data as a DataFrame

### dataframe.html
It displays the parsed FASTA data in a tabular format and allows different operations.
- **Navigation**: Provides links to other possible features:
  - Statistics visualization ("/stats")
  -  Motif searching ("/motif")
  -  Pairwise sequence alignment ("/align")
  -  Return to home page
- **Data presentation**:
  - Renders a table with three columns: ID, Description and Sequence
  - Iterate through each sequence entry
  - Each row represents one sequence from the FASTA file with its metadata

### stats.html
It provides statistical analysis of the FASTA sequences:
- **Table structure**: Displays detailed metrics for each sequence:
  - ID and Description (identifiers from the FASTA file)
  - Length (number of nucleotides in each sequence)
  - GC Content (percentage of guanine and cytosine bases, an important genomic characteristic)
- **Data rendering**: loops through the "stats" array and populate the table
- **Navigation**: Includes links to the other possible operations 

### motifs.html
Enables users to search for conserved sequence patterns (motifs) within the genomic data and provides links to other possible operations.
- **Search options**: Provides two different search forms (both use POST method to submit to the "/motif" endpoint)
  - First form: Allows users to search by sequence index
  - Second form: Allows users to search by sequence motif, if known
- **Results display**:
  - For sequence index searches: Shows a Pandas Dataframe with 3 columns: motif, number of occurrences and offsets at which it was found
  - For motif pattern searches: Renders a list of sequences of the database in which this motif is present and how many times

### align.html
Provides sequence alignment capabilities and the link to other possible operations
- **User input**: Allows selection of two sequences by their indices in the DataFrame for pairwise global alignment
- **Form**: Uses numeric inputs and submits to "/align" endpoint using `POST` method
- **Results display**: 
  - Displays alignment results in a preformatted text block
  - The preformatted tag ensures that spacing and formatting of the alignment is preserved
	
---

## 7. Running procedures
Here are the commands to run the project from the terminal of your computer:
-  Install the dependencies. Use: `pip install dependency_name` to install the required libraries (flask, pandas, python3, biopython are needed).
-  Navigate to the project folder. Use the cd command to move into the directory containing your project files: `cd path_to_directory`
-  Run the Flask app with: `python3 app.py`

By default, the web app will be available at localhost:5000.
If you encounter issues (e.g. the app is slow to start), it’s possible that port 5000 is already in use. In that case, try running the app on a different port using: `flask run --port=8080`

---

## 8. Conclusions

### Current limitations and future possible improvements of the system:
- **File format**: The application only supports uploading files in FASTA format, limiting its usability for datasets stored in other common bioinformatics formats (e.g., GenBank or EMBL). This feature can be extended by using the superclass FileParser that we implemented.
- **Basic error handling**: For example, the system may not adequately handle malformed FASTA files or unexpected user inputs, which can lead to uninformative error messages or crashes.
- **Session-Based file storage**: The use of Flask sessions to store the file path is simple but not ideal for multi-user environments or larger datasets.
- **Limited analysis**: the software is focused on basic sequence statistics (length and GC content), motif search, and pairwise alignment. More advanced analyses, such as phylogenetic tree construction or structural predictions, are not supported. This can be improved by integrating additional bioinformatics tools and algorithms, such as for BLAST searches or multiple sequence alignments.
- **User Interface with static visualization**: The web interface is functional but basic. There is limited visualization for sequence alignments, motif occurrences, and other analytical results. Using JavaScript we could upgrade the web interface with interactive visualizations (e.g., graphical displays of alignments) to facilitate data interpretation.

### Lessons we learned during the project:
- **Modular Design is Crucial**: Organizing the application into distinct classes helped maintain the code easier to manage, understand and extend.
- **User-Centered Design Considerations**: the need for an intuitive interface and clear visualizations was clear from the beginning. 
- **Scalability role**: Handling large biological datasets is challenging. The project reinforced the importance of planning things with scalability and abstraction, both in terms of data processing and user queries.
