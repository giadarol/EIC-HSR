mad_src = '''
start_mk: marker;
yrot: yrotation,angle= -0.025;
end_mk: marker;
test_sequence: sequence, l = 2;
start_mk, at = 0;
yrot, at = 2e-6;
end_mk, at = 4e-6;
endsequence;
'''

from cpymad.madx import Madx

mad = Madx()
mad.input(mad_src)
mad.beam()
mad.use('test_sequence')
mad.twiss(betx=1, bety=1)

px_test = -0.02490454681197178
tw_test = mad.twiss(betx=1, bety=1, x=0, px=px_test)
px_out_twiss = tw_test.dframe()[tw_test.name=='end_mk:1'].px.values[0]

mad.input(f'''
track, onepass=true;
start, x=0, px={px_test}, y=0, py=0, t=0, pt=0;
observe, place=end_mk;
run, turns=1;
endtrack;
''')

print(f'''
px from twiss: {px_out_twiss}
px from track: {mad.table['track.obs0001.p0001'].px[-1]}
''')
# Output:
#    px from twiss: 0.00010063136828029548
#    px from track: 9.287801779160956e-05