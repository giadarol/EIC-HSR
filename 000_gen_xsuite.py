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
seq_thin = mad_thin.sequence[seq_name]

# Makethin does not work with zero angle on the bends (putting a very small one)
bends_with_no_angle= []
for ee in seq_thin.elements:
     if hasattr(ee, 'k0') and ee.k0 != 0 and ee.angle == 0:
        ee.angle = 1e-20
        bends_with_no_angle.append(ee.name)

tw_after_reload = mad_thin.twiss()

n_slices_sbend = 4
mad_thin.input(f'''
select, flag=makethin, clear;
select, flag=makethin, class=quadrupole, slice = 20, thick=false;
select, flag=makethin, class=sextupole, slice = 1, thick=true;

select, flag=makethin, class=sbend, slice=0;
select, flag=makethin, pattern=yi_*, slice=10, thick=false;
select, flag=makethin, pattern=yo_*, slice=10, thick=false;
select, flag=makethin, pattern=bi_*, slice=10, thick=false;
select, flag=makethin, pattern=d5, slice=10, thick=false;

select, flag=makethin, pattern=dsw_ir10, slice=30; thick=false;
select, flag=makethin, pattern=dwarm_ir10_41, slice=30; thick=false;
select, flag=makethin, pattern=dwarm_ir12_41, slice=30; thick=false;
select, flag=makethin, pattern=dsw_ir12, slice=30; thick=false;

select, flag=makethin, pattern=bxds9m2, slice=30, thick=false;
select, flag=makethin, pattern=bxdsds04, slice=30, thick=false;
select, flag=makethin, pattern=bxds9m1, slice=30, thick=false;
select, flag=makethin, pattern=bxds02disp1, slice=30, thick=false;
select, flag=makethin, pattern=bxds01b, slice=30, thick=false;
select, flag=makethin, pattern=bxds01a, slice=30, thick=false;
select, flag=makethin, pattern=bxsp01, slice=30, thick=false;
select, flag=makethin, pattern=bxus9m1, slice=30, thick=false;
select, flag=makethin, pattern=bxus9m2, slice=30, thick=false;
select, flag=makethin, pattern=bxus9m3, slice=30, thick=false;


!select, flag=makethin, pattern=dwarm3
!select, flag=makethin, pattern=dwarm4, slice=10; thick=false;
!select, flag=makethin, pattern=b3pf, slice=30; thick=false;
!select, flag=makethin, pattern=b2pf, slice=30; thick=false;
!select, flag=makethin, pattern=b1apf, slice=30; thick=false;
!select, flag=makethin, pattern=b1pf, slice=30; thick=false;
!select, flag=makethin, pattern=b0apf, slice=30; thick=false;
!select, flag=makethin, pattern=b0pf, slice=30; thick=false;

makethin, sequence={seq_name}, style=teapot, makedipedge=true;
''')
mad_thin.use(seq_name)

tw_thin = mad_thin.twiss()

# tw_thin = mad_thin.twiss(
#      x=tw_thick.x[0], px=tw_thick.px[0], y=tw_thick.y[0], py=tw_thick.py[0],
#      betx=tw_thick.betx[0], bety=tw_thick.bety[0], alfx=tw_thick.alfx[0],
#      alfy=tw_thick.alfy[0]
# )

#################################
# Build line from MAD-X lattice #
#################################

# line = xt.Line.from_madx_sequence(sequence=mad_thin.sequence[seq_name],
#            deferred_expressions=False, install_apertures=False,
#            apply_madx_errors=False)
# line.particle_ref = xp.Particles(gamma=seq_thin.beam.gamma,
#                                  mass0=xp.PROTON_MASS_EV,
#                                  q0=seq_thin.beam.charge)
# line.build_tracker()
# tw_xs = line.twiss(method='4d')

plt.close('all')

plt.figure()
ax1 = plt.subplot(2,1,1)
plt.plot(tw_thick['s'], tw_thick['betx'], label='betx thick')
plt.plot(tw_thick['s'], tw_thick['bety'], label='bety thick')
plt.plot(tw_thin['s'], tw_thin['betx'], label='betx thin')
plt.plot(tw_thin['s'], tw_thin['bety'], label='bety thin')
# plt.plot(tw_xs['s'], tw_xs['betx'], label='betx xsuite')
# plt.plot(tw_xs['s'], tw_xs['bety'], label='bety xsuite')

ax2 = plt.subplot(2,1,2, sharex=ax1)
plt.plot(tw_thick['s'], tw_thick['x'], label='x thick')
plt.plot(tw_thick['s'], tw_thick['y'], label='y thick')
plt.plot(tw_thin['s'], tw_thin['x'], label='x thin')
plt.plot(tw_thin['s'], tw_thin['y'], label='y thin')

plt.legend()
plt.show()

#context = xo.ContextCpu()



plt.show()
