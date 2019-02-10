from . import estimation

class Frb:
    def __init__(self, name, dm, dm_err=None):
        self.name = name
        self.dm = dm
        if dm_err:
            self.dm_err = dm_err

    def __repr__(self):
        return 'Frb({0})'.format(vars(self))

    def calc_redshift(self, method='batten2019'):

        z, z_err = estimation.calc_redshift(self.dm, self.dm_err, method)

        self.z = z
        self.z_err = z_err

        return z, z_err
