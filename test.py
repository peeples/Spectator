
from astropy.table import Table, Column
from astropy import units as u
import numpy as np


def test_table():
    g1 = Table()
    dummy = np.zeros(17)
    col_dummy = Column(name='age',data=dummy+11.,unit=u.yr)
    g1.add_column(col_dummy)
    col_dummy = Column(name='bar',data=dummy)
    g1.add_column(col_dummy)
    col_dummy = Column(name='foo',data=dummy)
    g1.add_column(col_dummy)
    g1['blah'] = "this is blah     "
    galaxies = {}
    for i in range(3):
        galaxies[i] = g1
    galaxies[2]['age'][12] = 4 #  this is really a key line for subscriptint 

    print "DONE BUILDING" 
   
    print galaxies[0]['blah'] 
    g1['blah'] = 'changed it' 
    print galaxies[0]['blah']
    g1['blah'] = 'changed it again' 
    print galaxies[0]['blah']

    galaxies[1] = g1 
    print galaxies[2]['age'][12] #  this is really a key line for subscripting 
    print galaxies[2]['age'][0:-1] 
    print galaxies[0]['foo']

    print galaxies.keys() 

    galaxies['jt'] = 'jason tumlinson' 
    print galaxies.keys() 

    print galaxies['jt']  
