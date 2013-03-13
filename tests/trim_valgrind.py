#!/usr/bin/env python
import re

class ValgrindReader:
    def __init__(self, filename):
        con = file(filename, 'r')
        m = re.match("^==([0-9]+)== ", con.next())
        if not m:
            raise Exception("Does not look like a valgrind file")
        self.suppress_prefix = '==%s== \n' % m.group(1)
        self.errs = []
        while True:
            next_error = self._get_1_error(con)
            if next_error:
                next_error = ''.join(next_error)
                if next_error not in self.errs:
                    self.errs.append(next_error)
            else:
                break
        con.close()

    def _get_1_error(self, con):
        ret = []
        for s in con:
            if s == self.suppress_prefix:
                s = con.next()
                if s == "{\n":
                    ret.append(s)
                    while True:
                        s = con.next()
                        ret.append(s)
                        if s == "}\n":
                            break
                    return ret

    def write_out(self, filename):
        con = open(filename, 'w+')
        con.writelines(self.errs)
        con.close()

if __name__ == "__main__":
    s = ValgrindReader("valgrind.out")
    s.write_out("valgrind.supp")
