from monty.json import MSONable, MontyDecoder, MontyEncoder

import json

class BaseSchema(MSONable):

    @classmethod
    def from_file(cls, fname):
        with open(fname, "rb") as f:
            return json.load(f, cls=MontyDecoder)
    
    def to_file(self, fname):
        with open(fname, "w+") as f:
            json.dump(self.as_dict(), f, cls=MontyEncoder)
    