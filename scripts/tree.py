import rpy2.robjects as robjects
robjects.r("library(tree)")

class TreeModel():
    def __init__(self, time_disturbance, slope):
        self.time_disturbance = time_disturbance
        self.slope = slope
        self.sys = None
        self.landscape = None
        self.landscape_approximate = None
    def set_residents(self, lma, seed_rain, equilibrium=True):
        assert len(lma) == len(seed_rain)
        f = robjects.r['for_python']
        self.sys = f(self.time_disturbance,
                     self.slope,
                     robjects.FloatVector(lma),
                     robjects.FloatVector(seed_rain),
                     equilibrium)
        self.landscape = self.sys['make_landscape']()
        self.landscape_approximate = None
    def fitness(self, lma):
        res = self.landscape(robjects.FloatVector(lma))
        return r2list(res)
    def setup_fitness_approximate(self, method="gp"):
        f_approximate = robjects.r['fitness_landscape_approximate']
        method = robjects.StrVector([method])
        self.landscape_approximate = f_approximate(self.sys, method)
    def fitness_approximate(self, lma):
        if self.landscape_approximate is None:
            self.setup_fitness_approximate()
        res = self.landscape_approximate(robjects.FloatVector(lma))
        return r2list(res)
    def grid(self, n=50):
        res = robjects.r['seq_log_range'](self.bounds, n)
        return r2list(res)
    @property
    def bounds(self):
        return r2list(self.sys['bounds'])
    @property
    def gp(self):
        return GP(robjects.r['environment'](self.landscape_approximate)['res'])

# The underlying GP works on *log* X.  This is a bit of a pain, but
# it is what it is.
#
# The R interface to the GP code is very basic at the moment, but
# it'll be easy enough to pass things through here if they're useful.
class GP():
    def __init__(self, obj):
        if type(obj) is not robjects.environments.Environment or \
             robjects.r['class'](obj)[1] != 'gpreg':
            raise Exception("Does not look like a GP object")
        self.gp = obj
    @property
    def X(self):
        return r2list(self.gp['X'])
    @property
    def y(self):
        return r2list(self.gp['y'])
    def predict(self, x):
        res = self.gp['predict1'](x)
        return (res[0], robjects.r['attr'](res, 'variance')[0])

# Not really clear what the right way to do this is.  With numpy:
#   import rpy2.robjects.numpy2ri as rpyn
#   vector=rpyn.ri2numpy(vector_R)
def r2list(obj):
    return [i for i in obj]
