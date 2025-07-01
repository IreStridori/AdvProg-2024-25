'''versione 9 mattina da GITHUB'''
from abc import ABC, abstractmethod
import pandas as pd
from Bio.Align import PairwiseAligner
from collections import Counter

class FileParser(ABC):
    def __init__(self):
        self._data = None
        
    @abstractmethod
    def parse_file(self, file):
        pass


# Classe per il parsing del file FASTA
class FastaParser:

    def __init__(self):
        self._data = None

    @staticmethod  #forse così? così almeno come argomento gli va bene solo function
    # Decoratore per validare la presenza del dataset caricato
    def ensure_data_loaded(func):
        def wrapper(self, *args, **kwargs):
            if self._data is None:  # o si cambia come è scritto il decoratore oppure va messo dentro la classe perchè altrimenti non può accedere direttamente al self.__data; in ogni caso dove era stato messso a decorare un'altra funzione non andava bene
                raise ValueError("Nessun dataset caricato. Caricare un file FASTA prima di eseguire l'operazione.")
            return func(self, *args, **kwargs)
        return wrapper

    def parse_file(self, file: str) -> None:
        """Effettua il parsing del file FASTA e memorizza i dati in un DataFrame."""
        sequences = []
        f = open(file, "r")
        identifier, description, sequence = None, None, '' #crea queste tre variabili
        for line in f:
            line = line.strip()  # Rimuove eventuali spazi vuoti o caratteri di nuova linea
            if line.startswith(">"):  # Identifica una nuova sequenza
                if identifier is not None:  # Se c'è già una sequenza in corso, la aggiungiamo alla lista
                    sequences.append((identifier, description, sequence))
                sequence = ""
                header_parts = line[1:].split(" ", 1)  # divide la stringa in due parti usando lo spazio " " come separatore, ma solo al primo spazio che trova.
                identifier = header_parts[0]  # Primo elemento
                description = header_parts[1] if len(header_parts) > 1 else ""  # Se presente, prende la descrizione
            else:
                sequence += line  # Aggiunge la riga
        #sequences = sequences[1:] #rimuove il primo elemento che è none,none, ""
        #dobbiamo riaggiungere  if identifier != None?
        sequences.append((identifier, description, sequence)) #serve a riaggiornare le sequenze

        f.close()
        self._data = pd.DataFrame(sequences, columns=["Identifier", "Description", "Sequence"])

    """Restituisce il DataFrame contenente i dati FASTA."""
    def get_DataFrame(self):
        return self._data

    @ensure_data_loaded #se i dati non sono stati caricati non verrà eseguito.
    def get_summary(self):
        """describe(include="all") restituisce un riepilogo del dataset.: numero di valori unici per ogni colonna, frequenza degli identificatori, lunghezza media delle sequenze"""
        return self._data.describe(include="all")

    @ensure_data_loaded #aggiunto
    def get_seq(self, index):
        """Restituisce la sequenza e i dettagli della riga indicata."""
        if index >= len(self._data) or index < 0:
            raise IndexError(f"Indice {index} fuori dai limiti).")
        seq = self._data.iloc[index].tolist()
        return seq
# superclasse per sequenze
class GenomicEntity(ABC):
    def __init__(self, identifier, description, sequence):
        self._identifier = identifier #protected così che puoi accedere da mithochondrial dna
        self._description = description
        self._sequence = sequence

    def get_attributes_value(self):
        return self._identifier, self._description, self._sequence

    def length(self):
        return len(self._sequence)


# Classe per rappresentare il DNA mitocondriale
class MitochondrialDNA(GenomicEntity):

    def gc_content(self):
        """Calcola il contenuto GC della sequenza."""
        g_count = self._sequence.count("G")
        c_count = self._sequence.count("C")
        return (g_count + c_count) / len(self._sequence) * 100

    def extract_subseq_by_indexing(self, start, end):
        return self._sequence[start:end + 1]
    
#CLASSE ASTRATTA PERCHè POTREI NON VOLER LAVORARE CON UN DATA FRAME MA CON UNA STRINGA, O QUALSIASI COSA, QUINDI FACCIAMO IN MODO CHE IL NOSTRO SIA UN CASO PARTICOLARE
class MotifAnalyser(ABC):
    def __init__(self, data):
        self._data = data #nel nostro caso il dato sarà il dataframe

    @abstractmethod
    def find_motif(self, motif):
        pass

# Classe per l'analisi dei motivi genetic
class SequenceMotif(MotifAnalyser): #inerita init da superclasse
    #Estrae tutte le subsequenze di lunghezza motif_length dalla sequence e restituisce una lista di motivi.
    # Conta le occorrenze di ciascun motivo (usa Counter() di collections). Filtra i motivi che compaiono più di minimum volte. Restituisce un DataFrame con i motivi e la loro frequenza.
    def extract_motifs(self, sequence, motif_length, minimum):
        motifs = [sequence[i:i + motif_length] for i in range(len(sequence) - motif_length + 1)]
        motif_counts = dict(Counter(motifs))
        recurrent_motifs = {motif: count for motif, count in motif_counts.items() if count >= minimum}
        motifs_df = pd.DataFrame(list(recurrent_motifs.items()), columns=["motif", "count"]).set_index("motif")
        return motifs_df

    #Cerca un motivo genetico specifico in tutte le sequenze del dataframe e ti dice quante volte l'ha trovato per ogni sequenza
    def find_motif(self, motif):
        
        motif = motif.upper()
        results = []
        for _, row in self._data.iterrows(): #!!!index
            count = row["Sequence"].count(motif)
            results.append((row["Identifier"], motif, count))
        return pd.DataFrame(results, columns=["Identifier", "Motif", "Occurrences"])
    
# Classe per l'allineamento
class SequenceAlignment:
    def __init__(self, seq1, seq2): #Inizializza l'oggetto di allineamento con due sequenze.
        self.__seq1 = seq1
        self.__seq2 = seq2
        self.__aligner = PairwiseAligner()
        self.__aligner.mode = 'global'  # Imposta l'allineamento globale di default
        self.__alignments = None # Variabile per memorizzare i risultati dell'allineamento

    def perform_alignment(self): #esegue l'allineamento
        self.__alignments = self.__aligner.align(self.__seq1, self.__seq2)
        return self.__alignments

    def format_alignment(self, n=1): #Formatta e restituisce i primi N risultati di allineamento
        results = []
        for alignment in self.__alignments[n]:
            results.append(str(alignment))
        return "\n".join(results)

    def alignments_score(self): #Questa funzione restituisce il punteggio dell'allineamento (ovvero la qualità dell'allineamento stesso), che è una misura di quanto siano simili le due sequenze.
        return self.__alignments.score
