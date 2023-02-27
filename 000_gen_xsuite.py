from cpymad.madx import Madx
from matplotlib import pyplot as plt

import xobjects as xo
import xpart as xp
import xtrack as xt

#################################
# Setting up MAD-X lattice      #
#################################

mad_thick = Madx()
mad_thick.call('hsr.madx')
mad_thick.call('ir12-ps.madx')
mad_thick.beam(energy=41.0,particle='antiproton')
mad_thick.input('''
dsw_ir10, k0 = -dwarm_ir10_41->k0;
dsw_ir12, k0 = +dwarm_ir12_41->k0;
''')
mad_thick.call('rhic-str-041-05.madx')
mad_thick.call('ir2-str-041-05.madx')
mad_thick.call('ir4-str-041-05.madx')
mad_thick.call('ir6-str-041-05.madx')
mad_thick.call('ir8c-str-041-05.madx')
mad_thick.call('ir10-str-041-05.madx')
mad_thick.call('ir12-str-041-05.madx')
seq_name = 'hsr41'
mad_thick.use(seq_name)
tw_thick = mad_thick.twiss()
seq_thick = mad_thick.sequence[seq_name]

mad_thick.input(f'''
save, sequence={seq_name}, file="{seq_name}_no_expr.seq", noexpr=true;
''')

mad_thin = Madx()
mad_thin.call(f'{seq_name}_no_expr.seq')
mad_thin.beam(energy=41.0,particle='antiproton')
mad_thin.use(seq_name)
tw_after_reload = mad_thin.twiss()

n_slices_sbend = 4
mad_thin.input(f'''
select, flag=makethin, clear;
select, flag=makethin, class=rbend, slice = 4, thick = false;
select, flag=makethin, class=sbend, slice = 4, thick = false;
select, flag=makethin, class=quadrupole, slice = 10, thick = false;
select, flag=makethin, class=sextupole, slice = 4, thick=false;
makethin, sequence={seq_name}, style=teapot, makedipedge=true;
''')
mad_thin.use(seq_name)
seq_thin = mad_thin.sequence[seq_name]
#tw_thin = mad_thin.twiss(betx=1, bety=1, alfx=0, alfy=0, dx=0, dy=0)

for ii, ee in list(enumerate(seq_thin.elements))[::100]:
    nn = ee.name
    print(nn)
    mad_thin.input(f'use,sequence=hsr41,range=#s/{nn};')
    mad_thin.twiss(betx=1, bety=1)

# This works
ii = 3490
nn = seq_thin.elements[ii].name
mad_thin.input(f'use,sequence=hsr41,range=#s/{nn};')
mad_thin.twiss(betx=1, bety=1)

# This does not work
ii = 3491
nn = seq_thin.elements[ii].name
mad_thin.input(f'use,sequence=hsr41,range=#s/{nn};')
mad_thin.twiss(betx=1, bety=1)

for ee in seq_thick.expanded_elements:
    print(ee.name, ee.base_type.name)
    if ee.base_type.name == 'sbend':
        break

prrrr
#################################
# Build line from MAD-X lattice #
#################################

line = xt.Line.from_madx_sequence(sequence=mad_thin.sequence[seq_name],
           deferred_expressions=False, install_apertures=False,
           apply_madx_errors=False)
line.particle_ref = xp.Particles(gamma=seq_thin.beam.gamma,
                                 mass0=xp.PROTON_MASS_EV,
                                 q0=seq_thin.beam.charge)

twiss_table_xt = xt.Tracker(line=line.copy()).twiss()
plt.figure(0)
plt.plot(twiss_table['s'],twiss_table_xt['betx'],'xb')
plt.plot(twiss_table['s'],twiss_table_xt['bety'],'xg')

#context = xo.ContextCpu()



plt.show()
