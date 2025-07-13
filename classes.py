from abc import ABC, abstractmethod
import pandas as pd
from Bio.Align import PairwiseAligner


def ensure_data_loaded(func):
    """Decoratore per validare la presenza del dataset caricato"""
    def wrapper(self, *args, **kwargs):
        if self._data is None:  
            raise ValueError("Nessun dataset caricato. Caricare un file FASTA prima di eseguire l'operazione.")
        return func(self, *args, **kwargs)
    return wrapper

class FileParser(ABC):
    """Classe astratta per il parsing dei files"""
    def __init__(self):
        self._data = None
        
    @abstractmethod
    def parse_file(self, file):
        pass


class FastaParser(FileParser): #eredita init da FileParser
    """ Classe per il parsing del file FASTA """
    
    def parse_file(self, file: str) -> None:
        sequences = []
        f = open(file, "r")
        identifier, description, sequence = None, None, '' 
        for line in f:
            line = line.strip()  # Rimuove eventuali spazi vuoti o caratteri di nuova linea
            if line.startswith(">"):  # Identifica una nuova sequenza
                if identifier is not None:  # Se c'è già una sequenza in corso, la aggiungiamo alla lista
                    sequences.append((identifier, description, sequence))
                sequence = ""
                header_parts = line[1:].split(" ", 1)  # divide la stringa in due parti usando lo spazio " " come separatore, ma solo al primo spazio che trova.
                identifier = header_parts[0]  
                description = header_parts[1] if len(header_parts) > 1 else ""  # Se presente, prende la descrizione
            else:
                sequence += line  # Aggiunge la riga
        #sequences = sequences[1:] #rimuove il primo elemento che è none,none, ""
        sequences.append((identifier, description, sequence)) 

        f.close()
        self._data = pd.DataFrame(sequences, columns=["Identifier", "Description", "Sequence"])

    def get_DataFrame(self):
        """Restituisce il DataFrame contenente i dati FASTA."""
        return self._data

    @ensure_data_loaded #se i dati non sono stati caricati non verrà eseguito.
    def get_summary(self):
        return self._data.describe(include="all") #restituisce un riepilogo del dataset: numero di valori unici per ogni colonna, frequenza degli identificatori, lunghezza media delle sequenze
        

    @ensure_data_loaded
    def get_row(self, index: int ):
        """Restituisce la sequenza e i dettagli della riga indicata."""
        if index >= len(self._data) or index < 0:
            raise IndexError(f"Index {index} out of range.")
        row = self._data.iloc[index].tolist()
        return row


class GenomicEntity():
    """ Superclasse per sequenze """
    def __init__(self, identifier, description, sequence):
        self._identifier = identifier #protected così che puoi accedere da mithochondrial dna
        self._description = description
        self._sequence = sequence

    def get_attributes_value(self):
        return self._identifier, self._description, self._sequence

    def length(self):
        return len(self._sequence)


class MitochondrialDNA(GenomicEntity):
    """ Classe per rappresentare il DNA mitocondriale """
    
    def gc_content(self):
        """Calcola il contenuto GC della sequenza."""
        g_count = self._sequence.count("G")
        c_count = self._sequence.count("C")
        return (g_count + c_count) / len(self._sequence) * 100

    def extract_subseq_by_indexing(self, start, end):
        return self._sequence[start:end + 1]
    

class MotifAnalyser(ABC):
    """ CLASSE ASTRATTA per la ricerca di mmotifi perché potrei voler lavorare con altre strutture dati oltre al DataFrame"""
    def __init__(self, data):
        self._data = data #nel nostro caso il dato sarà il dataframe

    @abstractmethod
    def find_motif(self, motif):
        pass


class SequenceMotif(MotifAnalyser): #inerita init da superclasse
    """ Classe per l'analisi dei motivi genetici """
    
    def extract_motifs(self, seq_idx, motif_length, minimum):
        """ Estrae tutte le sottosequenze di lunghezza motif_length. Filtra i motivi che compaiono più di minimum volte. Restituisce un DataFrame con i motivi e la loro frequenza """
        motifs=[]
        motifs_dic={}
        sequence= self._data.iloc[seq_idx]['Sequence']
        for i in range(len(sequence) - motif_length + 1):
            motif=sequence[i:i + motif_length]
            if motif in motifs_dic:
                motifs_dic[motif].append(i)
            else:
                motifs_dic[motif]=[i] #so that then the frequency should be just the length of the list of initial indexes
        md=motifs_dic.copy() #otherwise it changes size while itinerating if we itinerate over the original
        for k, v in md.items():
            if len(v)<= minimum:
                del motifs_dic[k]
            else:
                motifs.append([k, v])
        motifs_df=pd.DataFrame(motifs, columns=['Motif', 'Indexes'])
                
        return motifs_df
    
    
    def find_motif(self, motif):
        """ Cerca un motivo genetico specifico in tutte le sequenze del dataframe e ti dice quante volte l'ha trovato per ogni sequenza """
        motif = motif.upper()
        results = []
        for _, row in self._data.iterrows(): #!!!index
            count = row["Sequence"].count(motif)
            results.append((row["Identifier"], motif, count))
        return pd.DataFrame(results, columns=["Identifier", "Motif", "Occurrences"])


class SequenceAlignment:
    def __init__(self, seq1:str, seq2:str): #Inizializza l'oggetto di allineamento con due sequenze.
        self.__seq1 = seq1
        self.__seq2 = seq2
        self.__aligner = PairwiseAligner()
        self.__aligner.mode = 'global'  # Imposta l'allineamento globale di default
        self.__alignments = None # Variabile per memorizzare i risultati dell'allineamento

    def perform_alignment(self): 
        """ Esegue l'allineamento """
        self.__alignments = self.__aligner.align(self.__seq1, self.__seq2)
        return self.__alignments

    def format_alignment(self, n=1) -> str: 
        """ Formatta e restituisce i primi N risultati di allineamento """
        results = []
        for alignment in self.__alignments[n]:
            results.append(str(alignment))
        return "\n".join(results)

    def alignments_score(self) -> float: 
        """ Restituisce il punteggio dell'allineamento, misura di quanto siano simili le due sequenze. """
        return self.__alignments.score
