import sys
from datetime import datetime


class Callback:
    """
    Progress monitoring
    TODO: implement message levels (WARN, ERROR)
    """
    def __init__(self):
        self.t0 = datetime.now()
        self.last_msg_length = 0

    def callback(self, msg, replace=False):
        if replace:
            sys.stdout.write('\b'*self.last_msg_length)
        self.last_msg_length = sys.stdout.write(
            u"\U0001F9A0" + ' [{}] {}'.format(datetime.now() - self.t0, msg)
        )
        if not replace:
            sys.stdout.write('\n')
