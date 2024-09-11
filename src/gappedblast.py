"""This module defines the `GappedBlast` class."""

from Bio.Align import substitution_matrices
from collections import defaultdict
import logging
import multiprocessing
from typing import List


from alignment import Alignment, worker_gapped_extension
from config import Config
from database import Database
from sequence import Sequence

# Configure logger
logger = logging.getLogger(__name__)


class GappedBlast:
    """A class to represent the custom BLAST 2 process.

    Attributes
    ----------
    query_fasta : str
        Path to the FASTA file of the query.
    db_fasta : str
        Path to the FASTA file of the database.
    db : Database
        An instance of the `Database` class containing sequence records.
    query : Sequence
        An instance of the `Sequence` class representing the query sequence.
    output : str
        Path to the output file.
    hits : Dict[str, List[Tuple[str, str]]]
        A dictionary where each key is a sequence ID and the value is a list of hits.
    """

    def __init__(self, params):
        """Initialize the GappedBlast instance.

        Parameters
        ----------
        params : Namespace
            Object containing the parameters of a blast run.
        """
        self.db_fasta = params.db
        self.query_fasta = params.query
        self.output = params.output
        if params.k != 3:
            Config.update_param("K", params.k)
        if params.matrix != "blosum62":
            Config.update_param(
                "MATRIX", substitution_matrices.load(params.matrix.upper())
            )
        if params.evalue != 0.001:
            Config.update_param("EVALUE", params.evalue)

    def load_data(self):
        """Load data from the provided FASTA files.

        Notes
        -----
        Attempts to load the database and query sequences.
        If loading fails, prompts the user for new file paths.
        """
        logger.info("Gapped-BLAST: Loading data...")
        ok = False
        while not ok:
            try:
                self.db = Database(self.db_fasta)
                self.query = Sequence.from_fasta(self.query_fasta)
                ok = True
            except FileNotFoundError as e:
                logger.error(f"\033[31m{e}\033[0m")
            except ValueError as e:
                logger.error(f"\033[31m{e}\033[0m")
            finally:
                if not ok:
                    self.db_fasta = input(
                        "\033[33mPlease re-enter the path to the database file: \033[0m"
                    )
                    self.query_fasta = input(
                        "\033[33mPlease re-enter the path to the query file: \033[0m"
                    )

    def hits_detection(self):
        """Detects hits between the query sequence and the database.

        Notes
        -----
        A hit is detected if the alignment score is greater than `T`.
        """
        logger.info("Gapped-BLAST: Indexation...")
        self.db.load_index()
        logger.info("Gapped-BLAST: Searching hits...")
        hits = defaultdict(list)
        q_words = self.query.words
        index = self.db.index
        size = len(q_words)
        for i, q_word in enumerate(q_words):
            if i % 10 == 0 and i != 0:
                logger.info(f"{i}/{size} query words processed...")
            for db_word, db_position in index.items():
                score = Alignment.compute_ungapped_score(q_word, db_word)
                if score <= Config.T:
                    continue
                for seq_id, _ in db_position:
                    hits[seq_id].append((q_word, db_word))
        self.hits = hits

    def ungapped_extension(self) -> List[Alignment]:
        """Extend the hits without allowing for gaps.

        Returns
        -------
        List[Alignment]
            A list of HSP to extend with gaps.

        Notes
        -----
        For each sequence in the database with detected hits, this method
        extends the hits without allowing for gaps. The goal is to
        determine the HSP that have the most potential for gapped extension.
        """
        logger.info("Gapped-BLAST: Extending hits without gaps...")
        alignments = []
        for sequence_index, hits in self.hits.items():
            logger.info(f"Extending {sequence_index}/{len(self.hits)}")
            q_record = self.query
            db_record = self.db.records[sequence_index]
            db_record.id = sequence_index
            alignments += Alignment.extend_to_hsp(
                q_record, db_record, hits
            )
        return alignments

    def gapped_extension(
        self, ungapped_alignments: List[Alignment]
    ) -> List[Alignment]:
        """Extend all the HSP previously found with gaps.

        Parameters
        ----------
        ungapped_alignments : List[Alignment]
            A list of HSP to extend with gaps.
        """
        logger.info("Gapped-BLAST: Extending hits with gaps...")
        i = 0
        gapped_alignments = []
        for hsp in ungapped_alignments:
            if i % 10 == 0 and i != 0:
                logger.info(f"{i} alignments have been extended...")
            seed = hsp.find_best_seed()
            forward = Alignment(
                hsp.seq_a[seed[0] :],
                hsp.seq_b[seed[1] :],
                0,
                0,
                1,
            ).needleman_wunsch_local_alignment()
            sub_a = hsp.seq_a[: seed[0] - 1]
            sub_b = hsp.seq_b[: seed[1] - 1]
            backward = Alignment(
                sub_a[::-1],
                sub_b[::-1],
                0,
                0,
                1,
            ).needleman_wunsch_local_alignment()
            backward.seq_a = backward.seq_a[::-1]
            backward.seq_b = backward.seq_b[::-1]
            gapped_alignment = forward.merge(backward)
            gapped_alignments.append(gapped_alignment)
            i += 1
        return gapped_alignments

    def parallel_gapped_extension(
        self, ungapped_alignments: List[Alignment]
    ) -> List[Alignment]:
        """Optimized version of gapped_extension() using multiprocessing.

        Notes
        -----
        Each process is a worker that extends several HSP with gaps.

        """
        logger.info("Gapped-BLAST: Extending hits with gaps...")
        workers = multiprocessing.Pool(None)
        returned_alignments = [
            workers.apply_async(func=worker_gapped_extension, args=[hsp])
            for hsp in ungapped_alignments
        ]
        i = 0
        all_alignments = []
        size = len(ungapped_alignments)
        for alignment in returned_alignments:
            if i % 10 == 0 and i != 0:
                logger.info(f"{i}/{size} alignments have been extended...")
            all_alignments.append(alignment.get())
            i += 1
        workers.close()
        workers.join()
        return all_alignments

    def run(self):
        """Execute the BLAST process."""
        self.load_data()
        self.hits_detection()
        ungapped_alignments = self.ungapped_extension()
        gapped_alignments = self.parallel_gapped_extension(
            ungapped_alignments
        )
