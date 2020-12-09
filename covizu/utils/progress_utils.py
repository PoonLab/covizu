import sys
from datetime import datetime


class Callback:
    """
    Progress monitoring
    """
    def __init__(self, my_rank=None, nprocs=None, t0=datetime.now(),
                 info=u"\U0001F9FF", warn=u"\U0001F6A7", error=u"\U0001F9A0"):
        self.t0 = t0
        self.last_msg_length = 0

        # progress monitoring in MPI context
        self.my_rank = my_rank
        self.nprocs = nprocs
        self.flags = {'INFO': info, 'WARN': warn, 'ERROR': error}

    def callback(self, msg, level="INFO", replace=False):
        # TODO: colored text with "\x1b[1;31;40m {}\x1b[0m"
        # TODO:
        if replace:
            sys.stdout.write('\b'*self.last_msg_length)
        self.last_msg_length = sys.stdout.write(
            '{} [{}] {}'.format(
                ' ' if level not in self.flags else self.flags[level],
                datetime.now() - self.t0,
                msg
            )
        )
        if not replace:
            sys.stdout.write('\n')
