import sys
from subprocess import Popen, PIPE

class Blat:
    """Run blat"""

    def __init__(self, db, query, output = None, **kwargs):
        self.db = db
        self.query = query
        if output:
            self.output = output
        else:
            self.output = sys.stdout
        self.step = 5
        self.rep = 2253
        self.min_score = 0
        self.min_identity = 0
        for key in kwargs:
            if key == "step":
                self.step = kwargs[key]
            elif key == "rep":
                self.rep = kwargs[key]
            elif key == "min_score":
                self.min_score = kwargs[key]
            elif key == "min_identity":
                self.min_identity = kwargs[key]

    def run(self):
        args = "blat -stepSize={0} -repMatch={1} -minScore={2} " +\
                "-minIdentity={3} -noHead {4} {5} {6}"
        args = args.format(
                self.step,
                self.rep,
                self.min_score,
                self.min_identity,
                self.db,
                self.query,
                self.output
            )
        blat = Popen(args, shell=True, stdout=PIPE)
        self.stdout, self.stderr = blat.communicate()
        if self.stderr:
            raise IOError, "blat did not complete successfully"

