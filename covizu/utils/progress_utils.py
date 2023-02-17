import sys
from datetime import datetime

# FIXME: can this be implemented with Python's logging module?
callback_verbosity = {'DEBUG': 3, 'INFO': 2, 'WARN': 1, 'ERROR': 0}


class Callback:
    """
    Progress monitoring
    """
    def __init__(self, my_rank=None, nprocs=None, t0=None, verbosity="INFO",
                 debug_flag=u"\U0001F991",  # squid
                 info_flag=u"\U0001F3C4",  # surfer
                 warn_flag=u"\U0001F6A7",  # construction sign
                 error_flag=u"\U0001F9A0"  # virus
                 ):
        self.t0 = datetime.now() if t0 is None else datetime.fromtimestamp(t0)
        self.last_msg_length = 0
        self.verbosity = callback_verbosity.get(verbosity, 2)

        # progress monitoring in MPI context
        self.mpi_state = '' if nprocs is None else '[{}/{}]'.format(my_rank, nprocs)
        self.flags = {'INFO': info_flag, 'WARN': warn_flag, 'ERROR': error_flag, 'DEBUG': debug_flag}

    def callback(self, msg, level="INFO", replace=False):
        verbosity = callback_verbosity.get(level, 2)
        if verbosity <= self.verbosity:
            # TODO: colored text with "\x1b[1;31;40m {}\x1b[0m"
            if replace:
                sys.stdout.write('\b'*self.last_msg_length)
            self.last_msg_length = sys.stdout.write(
                '{} [{}]{} {}'.format(
                    ' ' if level not in self.flags else self.flags[level],
                    datetime.now() - self.t0,
                    self.mpi_state, msg
                )
            )
            if not replace:
                sys.stdout.write('\n')
            sys.stdout.flush()
