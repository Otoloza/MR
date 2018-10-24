#/usr/bin/env python

import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from scipy import interpolate


# Some definitions
Rs = 6.96e10 #solar radius
Ms = 1.9891e33 #solar mass


#----------------------------
# Boris (wd.sm)
#----------------------------

def mass_radius_sm(temp, log):
   '''
		do i = 1,dimen(t)-1 {
		   if (t[$i-1] <= $(TEFF) && $(TEFF) <= t[$i]) {
		      set t1 = t[$i-1]
		      set t2 = t[$i]
		      define ti1  $($i-1)
		      define ti2  $i
		   }
		}
		#
		set logg = {7 7.5 8 8.5 9 9.5}
		set dimen(m1) = dimen(logg)
		set dimen(m2) = dimen(logg)
		set dimen(r1) = dimen(logg)
		set dimen(r2) = dimen(logg)
		#
		do i = 0,dimen(logg)-1 {
		   set m1[$i] = m$(logg[$i])[$ti1]
		   set m2[$i] = m$(logg[$i])[$ti2]
		   set r1[$i] = r$(logg[$i])[$ti1]
		   set r2[$i] = r$(logg[$i])[$ti2]
		}mr.py~
		#
		#
		spline logg m1 LOGG ml
		spline logg m2 LOGG mh
		#
		set _mwd = ml + (mh-ml)*(TEFF-t1)/(t2-t1)	
		#
		spline logg r1 LOGG rl
		spline logg r2 LOGG rh
		set _rwd = rl + (rh-rl)*(TEFF-t1)/(t2-t1)

   '''
   data = np.genfromtxt('/home/astro/phrnbk/python_modules/MR/mr.dat', unpack=True, comments='#', names=True)
   t    = data['Teff']
   logg = np.array([7.0,7.5,8.0,8.5,9.0,9.5])

   if (temp in t) and (log in logg):
      return data['r'+str(log).replace(".","_")][np.where(t==temp)]/Rs, data['m'+str(log).replace(".","_")][np.where(t==temp)]

   else:
      teff_l, teff_u = max(t[t < temp]), min(t[t > temp])
      m_row_l = [data['m'+str(logg[i]).replace(".","_")][np.where(t==teff_l)] for i in range(len(logg))]
      m_row_u = [data['m'+str(logg[i]).replace(".","_")][np.where(t==teff_u)] for i in range(len(logg))]
      r_row_l = [data['r'+str(logg[i]).replace(".","_")][np.where(t==teff_l)] for i in range(len(logg))]
      r_row_u = [data['r'+str(logg[i]).replace(".","_")][np.where(t==teff_u)] for i in range(len(logg))]

      m_int_l = interpolate.interp1d(logg, m_row_l, kind='cubic')
      ml      = m_int_l(log)
      m_int_u = interpolate.interp1d(logg, m_row_u, kind='cubic')
      mu      = m_int_u(log)

      #interpolating the mass bilinear 
      mwd     = ml + (mu-ml)*(temp-teff_l)/(teff_u-teff_l)

      r_int_l = interpolate.interp1d(logg, r_row_l, kind='cubic')
      rl      = r_int_l(log)
      r_int_u = interpolate.interp1d(logg, r_row_u, kind='cubic')
      ru      = r_int_u(log)

      # interpolating the radius (bilinear)
      rwd     = rl + (ru-rl)*(temp-teff_l)/(teff_u-teff_l)
      return rwd/Rs , mwd
