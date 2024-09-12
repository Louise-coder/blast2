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
    hits : Dict[int, List[Tuple[str, str]]]
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
        logger.info("Loading data...")
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
        logger.info("Indexation...")
        self.db.load_index()
        logger.info("Searching hits...")
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
        logger.info("Extending hits without gaps...")
        alignments = []
        for sequence_index, hits in self.hits.items():
            if sequence_index != 0 and sequence_index % 10 == 0:
                logger.info(
                    f"{sequence_index}/{len(self.hits)} alignments have been extended..."
                )
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

        Returns
        -------
        List[Alignment]
            A list of HSP extended with gaps.

        Notes
        -----
        For each HSP, this method extends the alignment with gaps in both:
        - The forward direction.
        - The backward direction.
        """
        logger.info("Extending hits with gaps...")
        i = 0
        gapped_alignments = []
        for hsp in ungapped_alignments:
            seed = hsp.find_best_seed()
            forward = Alignment(
                hsp.seq_a[seed[0] :],
                hsp.seq_b[seed[1] :],
                0,
                0,
                1,
            ).needleman_wunsch_local_alignment()
            sub_a = hsp.seq_a[: seed[0]]
            sub_b = hsp.seq_b[: seed[1]]
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
            gapped_alignment.seq_id = hsp.seq_id
            gapped_alignments.append(gapped_alignment)
            i += 1
            logger.info(f"{i} alignments have been extended...")
        return gapped_alignments

    def parallel_gapped_extension(
        self, ungapped_alignments: List[Alignment]
    ) -> List[Alignment]:
        """Optimized version of gapped_extension() using multiprocessing.

        Notes
        -----
        Each process is a worker that extends several HSP with gaps.

        """
        logger.info("Extending hits with gaps...")
        workers = multiprocessing.Pool(None)
        returned_alignments = [
            workers.apply_async(func=worker_gapped_extension, args=[hsp])
            for hsp in ungapped_alignments
        ]
        i = 0
        all_alignments = []
        size = len(ungapped_alignments)
        for alignment in returned_alignments:
            all_alignments.append(alignment.get())
            i += 1
            logger.info(f"{i}/{size} alignments have been extended...")
        workers.close()
        workers.join()
        return all_alignments

    def display_results(self, gapped_alignments: List[Alignment]):
        """Display the results of the BLAST process.

        Parameters
        ----------
        gapped_alignments : List[Alignment]
            The list of local alignments extended with gaps.

        Notes
        -----
        The results are displayed in the console.
        """
        logger.info("Collecting results...")
        q_len = len(self.query)
        db_len = sum([len(record.seq) for record in self.db.records])
        print("\nGAPPED BLASTP\n\n\n")
        print(
            f"Database:  {self.db_fasta}\n{len(self.db.records)} sequences; {db_len} total letters\n\n\n"
        )
        print(
            f"\033[32mQuery= {self.query.name} {self.query.description}\033[0m\n"
        )
        print(f"Length={q_len}\n")
        print("Sequences producing significant alignments:\n")
        for alignment in gapped_alignments:
            alignment.compute_statistics(q_len, db_len)
        i = 0
        best_alignments = Alignment.keep_best_alignments(gapped_alignments)
        best_alignments.sort(key=lambda x: x.evalue)
        for alignment in best_alignments:
            if alignment.evalue > Config.EVALUE or i > Config.NB_RESULTS:
                break
            print(
                f"{self.db.records[alignment.seq_id].name:<30}\tScore={alignment.normalized_score:5.1f} bits\tE-value={alignment.evalue:.1e}"
            )
            i += 1
        print("\n\n")
        i = 0
        for alignment in best_alignments:
            if alignment.evalue > Config.EVALUE or i > Config.NB_RESULTS:
                break
            alignment.display_results(
                self.query, self.db.records[alignment.seq_id]
            )
            i += 1
        if not best_alignments:
            print("***** No hits found *****\n")

    def get_results(self, gapped_alignments: List[Alignment]):
        """Display the results of the BLAST process.

        Parameters
        ----------
        gapped_alignments : List[Alignment]
            The list of local alignments extended with gaps.

        Notes
        -----
        The results are gathered in a text file.
        """
        logger.info("Collecting results...")
        q_len = len(self.query)
        db_len = sum([len(record.seq) for record in self.db.records])
        content = "GAPPED BLASTP\n\n\n\n"
        content += f"Database:  {self.db_fasta}\n{len(self.db.records)} sequences; {db_len} total letters\n\n\n\n"
        content += f"Query= {self.query.name} {self.query.description}\n\n"
        content += f"Length={q_len}\n\n"
        content += "Sequences producing significant alignments:\n\n"
        for alignment in gapped_alignments:
            alignment.compute_statistics(q_len, db_len)
        i = 0
        best_alignments = Alignment.keep_best_alignments(gapped_alignments)
        best_alignments.sort(key=lambda x: x.evalue)
        for alignment in best_alignments:
            if alignment.evalue > Config.EVALUE or i > Config.NB_RESULTS:
                break
            content += f"{self.db.records[alignment.seq_id].name:<30}\tScore={alignment.normalized_score:5.1f} bits\tE-value={alignment.evalue:.1e}\n"
            i += 1
        content += "\n\n"
        i = 0
        for alignment in best_alignments:
            if alignment.evalue > Config.EVALUE or i > Config.NB_RESULTS:
                break
            content += alignment.get_results(
                self.query, self.db.records[alignment.seq_id]
            )
            i += 1
        if not best_alignments:
            content += "***** No hits found *****\n"
        with open(self.output, "w") as out:
            out.write(content)

    def run(self):
        """Execute the BLAST process."""
        self.load_data()
        self.hits_detection()
        ungapped_alignments = self.ungapped_extension()
        gapped_alignments = self.parallel_gapped_extension(
            ungapped_alignments
        )
        if self.output:
            self.get_results(gapped_alignments)
        else:
            self.display_results(gapped_alignments)

    def run_with_time(self):
        """Same as run but with time measurement for each step."""
        import time

        start = time.time()
        self.load_data()
        end = time.time()
        print("Load_Data : ", end - start)
        start = time.time()
        self.hits_detection()
        end = time.time()
        print("Hits_detection : ", end - start)
        start = time.time()
        ungapped_alignments = self.ungapped_extension()
        end = time.time()
        print("Ungapped extension : ", end - start)
        start = time.time()
        gapped_alignments = self.parallel_gapped_extension(
            ungapped_alignments
        )
        end = time.time()
        print("Gapped extension : ", end - start)
        if self.output:
            self.get_results(gapped_alignments)
        else:
            self.display_results(gapped_alignments)
