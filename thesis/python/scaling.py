import numpy as np
import scipy.constants as const
from table import (
    make_table,
    make_full_table,
    make_composed_table,
    make_SI,
    write,
    search_replace_within_file,
)

m = const.physical_constants["electron mass"]
m = m[0]
e = const.physical_constants["elementary charge"]
e = e[0]
hbar = const.physical_constants["Planck constant over 2 pi"]
hbar = hbar[0]
V0 = e*0.1768

tau = hbar/V0
xi = np.sqrt(hbar**2 / m / V0)

print(tau)
print(xi)

write('tex_files/tau.tex', make_SI(tau * 1e15, r'\second', exp='e-15', figures=2))
write('tex_files/xi.tex', make_SI(xi * 1e10, r'\meter', exp='e-10', figures=2))
