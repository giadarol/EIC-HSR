from cpymad.madx import Madx
import numpy as np
import matplotlib.pyplot as plt

sequence_src = '''
    none = 0;
    start_check: marker,l= 0,kmax= 0,kmin= 0,calib= 0,polarity= 0,type,apertype="circle",aperture={ 0},aper_offset={ 0},aper_tol={ 0},aper_vx={ -1},aper_vy={ -1},slot_id= 0,assembly_id= 0,mech_sep= 0,v_pos= 0,magnet= 0,model= -1,method= -1
    ,exact= -1,nst= -1,fringe= 0,bend_fringe=false,kill_ent_fringe=false,kill_exi_fringe=false,dx= 0,dy= 0,ds= 0,dtheta= 0,dphi= 0,dpsi= 0,aper_tilt= 0,comments;
    ptb_b0apf: translation,dx= -0.0003241247841;
    prb_b0apf: yrotation,angle= 0.007208408396;
    b0apf: sbend,l= 0.6,k0= 0.01818264338;
    pte_b0apf: translation,dx= -0.003353390115;
    pre_b0apf: yrotation,angle= -0.005049849171;
    ptb_b0pf: translation,dx= -0.03022642196;
    prb_b0pf: yrotation,angle= 0.02670439316;
    b0pf: sbend,l= 1.2,k0= -0.008660083518,k1= -0.05940545468;
    pte_b0pf: translation,dx= -0.0007490490899;
    pre_b0pf: yrotation,angle= -0.025;
    ptb_sol_f: translation,dx= -0.04999479183;
    prb_sol_f: yrotation,angle= 0.025;
    sol_half: solenoid,l= 2;
    star_detect_f: sol_half;
    pre_sol_f: yrotation,angle= -0.025;
    ip6w: marker;
    hsr41: sequence, l = 9;
    start_check, at = 0;
    ptb_b0apf, at = 0.3838603422, dx = -0.0003241247841 ;
    prb_b0apf, at = 0.3838603422;
    b0apf, at = 0.6838603422;
    pte_b0apf, at = 0.9838603422, dx = -0.003353390115 ;
    pre_b0apf, at = 0.9838603422;
    ptb_b0pf, at = 1.88513171, dx = -0.03022642196 ;
    prb_b0pf, at = 1.88513171;
    b0pf, at = 2.48513171;
    pte_b0pf, at = 3.08513171, dx = -0.0007490490899 ;
    pre_b0pf, at = 3.08513171;
    ptb_sol_f, at = 6.88756965, dx = -0.04999479183 ;
    prb_sol_f, at = 6.88756965;
    star_detect_f, at = 7.88756965;
    pre_sol_f, at = 8.88756965;
    ip6w, at = 8.88756965;
    endsequence;
'''

tw_init = {'betx': 114.90943340961454,
            'alfx': 16.652973701933835,
            'bety': 853.8161081185903,
            'alfy': 60.12421586470439,
            'x': 0.015189573975841864,
            'px': 0.003229353952276041,
            'y': 0.0,
            'py': 0.0}



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
select, flag=makethin, class=sextupole, slice = 1, thick=true;
select, flag=makethin, class=sbend, slice=0;

!select, flag=makethin, pattern=b0apf, slice=1000, thick=false;
select, flag=makethin, pattern=b0pf, slice=1000, thick=false;

makethin, sequence=hsr41, style=teapot, makedipedge=false;
''')
mad_thin.use('hsr41')

tw_thin = mad_thin.twiss(**tw_init)

bety_thin_on_thick_check = np.interp(tw_thick.s, tw_thin['s'], tw_thin['bety'])
alfy_thin_on_thick_check = np.interp(tw_thick.s, tw_thin['s'], tw_thin['alfy'])

plt.figure()
plt.plot(tw_thick['s'], tw_thick['bety']/bety_thin_on_thick_check -1)

s_b0apf = tw_thick.dframe().loc[tw_thick.name=='b0apf:1', 's'].values[0]
s_b0pf = tw_thick.dframe().loc[tw_thick.name=='b0pf:1', 's'].values[0]

plt.axvline(s_b0apf, color='k')
plt.axvline(s_b0pf, color='r')

plt.show()

