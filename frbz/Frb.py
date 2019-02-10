from . import redshift_estimation

class Frb:
    def __init__(self, name, dm, dm_err=None):
        self.dm = dm
        self.name = name
        if dm_err:
            self.dm_err = dm_err

    def __repr__(self):
        return 'Frb({0}, {1}, {2})'.format(self.name, self.dm, self.dm_err)

    def calc_redshift(self, method='batten2019'):

        z, z_err = redshift_estimation.calc_redshift(self.dm, self.dm_err, method)
  
        return z, z_err

new = Frb("Utmost1", 300, dm_err=2)
print(new)
print(new.calc_redshift())
