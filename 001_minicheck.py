from cpymad.madx import Madx
import numpy as np
import matplotlib.pyplot as plt

config = {
    'on_translation': 1,
    'on_rotation': 1,
    'on_k1': 1,
    'on_shifts': 1,
    'kill_orbit': False,
    'n_slices_bend': 1000
}

sequence_src = '''

option, thin_cf=false;

none = 0;
start_check: marker,l= 0,kmax= 0,kmin= 0,calib= 0,polarity= 0,type,apertype="circle",aperture={ 0},aper_offset={ 0},aper_tol={ 0},aper_vx={ -1},aper_vy={ -1},slot_id= 0,assembly_id= 0,mech_sep= 0,v_pos= 0,magnet= 0,model= -1,method= -1
,exact= -1,nst= -1,fringe= 0,bend_fringe=false,kill_ent_fringe=false,kill_exi_fringe=false,dx= 0,dy= 0,ds= 0,dtheta= 0,dphi= 0,dpsi= 0,aper_tilt= 0,comments;
ptb_b0pf: translation,dx= -0.03022642196;
prb_b0pf: yrotation,angle= 0.02670439316;
b0pf: sbend,l= 1.2,k0= -0.008660083518, k1:= -0.05940545468*on_k1 + corr_k1, angle=1e-30;
pte_b0pf: translation,dx= -0.0007490490899;
pre_b0pf: yrotation,angle= -0.025;
ptb_sol_f: translation,dx= -0.04999479183;
prb_sol_f: yrotation,angle= 0.025;
sol_half: drift,l= 2;
star_detect_f: sol_half;
pre_sol_f: yrotation,angle= -0.025;
ip6w: marker;
hsr41: sequence, l = 8;
start_check, at = 0;
ptb_b0pf, at = 0.8851317104, dx = -0.03022642196 ;
prb_b0pf, at = 0.8851317104;
b0pf, at = 1.48513171;
pte_b0pf, at = 2.08513171, dx = -0.0007490490899 ;
pre_b0pf, at = 2.08513171;
ptb_sol_f, at = 5.88756965, dx = -0.04999479183 ;
prb_sol_f, at = 5.88756965;
star_detect_f, at = 6.88756965;
pre_sol_f, at = 7.88756965;
ip6w, at = 7.88756965;
endsequence;
'''


new_lines = []
for ll in sequence_src.split('\n'):
    if 'translation' in ll:
        ll = ll.replace(';', "*on_translation;")
        ll = ll.replace('dx=', "dx:=")
    elif 'yrotation' in ll:
        ll = ll.replace(';', "*on_rotation;")
        ll = ll.replace('angle=', "angle:=")
    elif 'dx = ' in ll:
        ll = ll.replace(';', "*on_shifts;")
        ll = ll.replace('dx = ', "dx:=")
    new_lines.append(ll)
new_lines.append(f"on_translation={config['on_translation']};")
new_lines.append(f"on_rotation={config['on_rotation']};")
new_lines.append(f"on_k1={config['on_k1']};")
new_lines.append(f"on_shifts={config['on_shifts']};")
sequence_src = '\n'.join(new_lines)

tw_init = {'betx': 84.02557752667167,
            'alfx': 14.23087869857309,
            'bety': 737.8650722816444,
            'alfy': 56.03568671202535,
            'x': 0.014288215914637994,
            'px': -0.009838644402766588,
            'y': 0.0,
            'py': 0.0}

if config['kill_orbit']:
    tw_init['x'] = 0
    tw_init['px'] = 0
    tw_init['y'] = 0
    tw_init['py'] = 0

mad_thick = Madx()
mad_thick.input(sequence_src)
mad_thick.beam(energy=41.0,particle='antiproton')
mad_thick.use(sequence='hsr41')
tw_thick = mad_thick.twiss(**tw_init)


mad_thin = Madx()
mad_thin.input(sequence_src)
mad_thin.beam(energy=41.0,particle='antiproton')
mad_thin.use(sequence='hsr41')

mad_thin.input(f'''
select, flag=makethin, clear;
select, flag=makethin, class=quadrupole, slice = 20, thick=false;
select, flag=makethin, pattern=b0pf, slice={config['n_slices_bend']}, thick=false;

makethin, sequence=hsr41, style=teapot, makedipedge=true;
''')
mad_thin.use('hsr41')

tw_thin = mad_thin.twiss(**tw_init)

bety_thin_on_thick_check = np.interp(tw_thick.s, tw_thin['s'], np.sqrt(tw_thin['bety']))**2
alfy_thin_on_thick_check = np.interp(tw_thick.s, tw_thin['s'], tw_thin['alfy'])

plt.close('all')
plt.figure()
ax1 = plt.subplot(211)
plt.plot(tw_thick['s'][:-2], (tw_thick['bety']/bety_thin_on_thick_check -1)[:-2], '.-')
plt.ylabel(r'$\Delta \beta_y / \beta_y$')

ax2 = plt.subplot(212, sharex=ax1)
plt.plot(tw_thick['s'], tw_thick.x)
plt.plot(tw_thin['s'], tw_thin.x)

plt.ylabel(r'$x$')


s_b0pf = tw_thick.dframe().loc[tw_thick.name=='b0pf:1', 's'].values[0]

plt.axvline(s_b0pf, color='r')
plt.subplots_adjust(left=.15)
plt.show()

