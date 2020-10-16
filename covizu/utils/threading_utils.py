#threading
import threading


class AlignerThread(threading.Thread):
    def __init__(self, in_queue, out_queue):
        threading.Thread.__init__(self)
        self.in_queue = in_queue
        self.out_queue = out_queue

    def run(self):
        # unpack data payload
        header, seq, database, refseq = self.in_queue.get()
        # open db to check if sequences already aligned
        cursor, conn = open_connection(database)
        alignedseq = find_seq(conn, seq, refseq)
        if alignedseq == 0:
            self.in_queue.task_done()
            return 1
        self.out_queue.put((header, seq, alignedseq))
        conn.close()
        self.in_queue.task_done()
        return 0


def iterate_fasta_threaded(fasta, ref, database = 'data/gsaid.db'):
    """
    Function to iterate through fasta file
    :param cursor: sqlite database handler
    :param fasta: file containing sequences
    :param ref: file containing reference sequence, default-> NC_04552.fa
    """
    handle = gotoh2.iter_fasta(open(fasta))
    _, refseq = gotoh2.convert_fasta(open(ref))[0]

    in_queue = queue.Queue()
    out_queue = queue.Queue()

    # stream records into queue
    for h, s in handle:
        in_queue.put((h, s, database, refseq))

    # pairwise align sequences in queue
    for seq in range(0,in_queue.qsize()):
        AlignerThread(in_queue, out_queue).start()

    # block until queue is processed
    in_queue.join()

    # update database with aligned sequences
    cursor, conn = open_connection(database)
    while not out_queue.empty():
        header, seq, alignedseq = out_queue.get()
        insert_seq(cursor, seq, header, alignedseq)
        conn.commit()
        out_queue.task_done()
    conn.close()
    out_queue.join()


def iterate_handle_threaded(handle, ref, database = 'data/gsaid.db'):
    """
    Function to iterate through handle passed by Chunky bot
    :param cursor: sqlite database handler
    :param handle: list containing tuples (header, raw seq)
    :param ref: file containing reference sequence, default-> NC_04552.fa
    """

    _, refseq = gotoh2.convert_fasta(open(ref))[0]

    in_queue = queue.Queue()
    out_queue = queue.Queue()

    # stream records into queue
    for h, s in handle:
        in_queue.put((h, s, database, refseq))

    # pairwise align sequences in queue
    for seq in range(0,in_queue.qsize()):
        AlignerThread(in_queue, out_queue).start()

    # block until queue is processed
    in_queue.join()

    # update database with aligned sequences
    cursor, conn = open_connection(database)
    while not out_queue.empty():
        header, seq, alignedseq = out_queue.get()
        insert_seq(cursor, seq, header, alignedseq)
        conn.commit()
        out_queue.task_done()
    conn.close()
    out_queue.join()

