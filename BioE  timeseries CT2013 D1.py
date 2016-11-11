# -*- coding: utf-8 -*-
"""
Create emissions timeseries CTM 2013

@author: Jordan Capnerhurst 2016
"""

# import things
import matplotlib.pyplot as plt
import numpy as np
import netCDF4


# Set colour map
cool = cm = plt.get_cmap('jet')

# Import data
m1 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/JAN.nc', 'r')  # January
m2 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/FEB.nc', 'r')  # Febuary
m3 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/MAR.nc', 'r')  # March
m4 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/APR.nc', 'r')  # April
m5 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/MAY.nc', 'r')  # May
m6 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/JUN.nc', 'r')  # June
m7 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/JUL.nc', 'r')  # July
m8 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/AUG.nc', 'r')  # August
m9 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/SEP.nc', 'r')  # September
m10 = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/OCT.nc', 'r')  # October
mnov = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/NOV.nc', 'r')  # November
mdec = netCDF4.Dataset('O:/Honours_data/OEHYear/Complete year/DEC.nc', 'r')  # Decmber


# Import variables
# Store_bio

# [ time, source, lon, lat ] Jan
# plot Daily average
m11 = m1.variables['store_Bio'][:, 14:24, 0, :, :]
m12 = m1.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
jan = np.concatenate((m11, m12), axis=1)/81

# average
janav = np.mean(jan, axis=(0, 1))

janavres = np.reshape(janav, 3600)

# Reshape for month
janres = jan.reshape(744, 60, 60) 

###########################################

# [ time, source, lon, lat ] Feb Bigger due to comparison to 2011 ~~~ Temp
# plot Daily average
mt11 = m1.variables['temp_a'][:, 14:24, 0, :, :]
mt12 = m1.variables['temp_a'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
jant = np.concatenate((mt11, mt12), axis=1)

# average
jantav = np.mean(jant, axis=(0, 1))

jantavres = np.reshape(jantav, 3600)

# Reshape for month
jantres = jant.reshape(744, 60, 60) 

stdfebt = np.std(jantres, axis=(1, 2))

jantresav = np.mean(jantres, axis=(1,2))


##############################################################################

# [ time, source, lon, lat ] Feb Bigger due to comparison to 2011
# plot Daily average
m21 = m2.variables['store_Bio'][:, 14:24, 0, :, :]
m22 = m2.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
feb = np.concatenate((m21, m22), axis=1)/81

# average
febav = np.mean(feb, axis=(0, 1))

febavres = np.reshape(febav, 3600)

# Reshape for month
febres = feb.reshape(672, 60, 60) 

stdfeb = np.std(febres, axis=(1, 2))

febresav = np.mean(febres, axis=(1,2))

###########################################

# [ time, source, lon, lat ] Feb Bigger due to comparison to 2011 ~~~ Temp
# plot Daily average
mt21 = m2.variables['temp_a'][:, 14:24, 0, :, :]
mt22 = m2.variables['temp_a'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
febt = np.concatenate((mt21, mt22), axis=1)

# average
febtav = np.mean(febt, axis=(0, 1))

febtavres = np.reshape(febtav, 3600)

# Reshape for month
febtres = febt.reshape(672, 60, 60) 

stdfebt = np.std(febtres, axis=(1, 2))

febtresav = np.mean(febtres, axis=(1,2))

##############################################################################

# [ time, source, lon, lat ] Mar
# plot Daily average
m31 = m3.variables['store_Bio'][:, 14:24, 0, :, :]
m32 = m3.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
mar = np.concatenate((m31, m32), axis=1)/81

# average
marav = np.mean(mar, axis=(0, 1))

# Reshape for month
marres = mar.reshape(744, 60, 60) 

maravres = marav.reshape(3600)

###########################################

# [ time, source, lon, lat ] Feb Bigger due to comparison to 2011 ~~~ Temp
# plot Daily average
mt31 = m3.variables['temp_a'][:, 14:24, 0, :, :]
mt32 = m3.variables['temp_a'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
mart = np.concatenate((mt31, mt32), axis=1)

# average
martav = np.mean(mart, axis=(0, 1))

martavres = np.reshape(martav, 3600)

# Reshape for month
martres = mart.reshape(744, 60, 60) 

stdfebt = np.std(martres, axis=(1, 2))

martresav = np.mean(martres, axis=(1,2))

##############################################################################

# [ time, source, lon, lat ] Apr
# plot Daily average
m41 = m4.variables['store_Bio'][:, 14:24, 0, :, :]
m42 = m4.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
apr = np.concatenate((m41, m42), axis=1)/81

# average
aprav = np.mean(apr, axis=(0, 1))

# Reshape for month
aprres = apr.reshape(720, 60, 60) 

apravres = aprav.reshape(3600)

###########################################

# [ time, source, lon, lat ] Feb Bigger due to comparison to 2011 ~~~ Temp
# plot Daily average
mt41 = m4.variables['temp_a'][:, 14:24, 0, :, :]
mt42 = m4.variables['temp_a'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
aprt = np.concatenate((mt41, mt42), axis=1)

# average
aprtav = np.mean(aprt, axis=(0, 1))

aprtavres = np.reshape(aprtav, 3600)

# Reshape for month
aprtres = aprt.reshape(720, 60, 60) 

stdaprt = np.std(aprtres, axis=(1, 2))

aprtresav = np.mean(aprtres, axis=(1,2))

##############################################################################

# [ time, source, lon, lat ] May
# plot Daily average
m51 = m5.variables['store_Bio'][:, 14:24, 0, :, :]
m52 = m5.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
may = np.concatenate((m51, m52), axis=1)/81

# average
mayav = np.mean(may, axis=(0, 1))

# Reshape for month
mayres = may.reshape(744, 60, 60) 

mayavres = mayav.reshape(3600)

###########################################

# [ time, source, lon, lat ] Feb Bigger due to comparison to 2011 ~~~ Temp
# plot Daily average
mt51 = m5.variables['temp_a'][:, 14:24, 0, :, :]
mt52 = m5.variables['temp_a'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
mayt = np.concatenate((mt51, mt52), axis=1)

# average
maytav = np.mean(mayt, axis=(0, 1))

maytavres = np.reshape(maytav, 3600)

# Reshape for month
maytres = mayt.reshape(744, 60, 60) 

stdmayt = np.std(maytres, axis=(1, 2))

maytresav = np.mean(maytres, axis=(1,2))

##############################################################################

# [ time, source, lon, lat ] june
# plot Daily average
m61 = m6.variables['store_Bio'][:, 14:24, 0, :, :]
m62 = m6.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
jun = np.concatenate((m61, m62), axis=1)/81

# average
junav = np.mean(jun, axis=(0, 1))

# Reshape for month
junres = jun.reshape(720, 60, 60) 

junavres = junav.reshape(3600)

###########################################

# [ time, source, lon, lat ] Feb Bigger due to comparison to 2011 ~~~ Temp
# plot Daily average
mt61 = m6.variables['temp_a'][:, 14:24, 0, :, :]
mt62 = m6.variables['temp_a'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
junt = np.concatenate((mt61, mt62), axis=1)

# average
juntav = np.mean(junt, axis=(0, 1))

juntavres = np.reshape(juntav, 3600)

# Reshape for month
juntres = junt.reshape(720, 60, 60) 

stdjunt = np.std(juntres, axis=(1, 2))

juntresav = np.mean(juntres, axis=(1,2))

##############################################################################

# [ time, source, lon, lat ] July
# plot Daily average
m71 = m7.variables['store_Bio'][:, 14:24, 0, :, :]
m72 = m7.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
jul = np.concatenate((m71, m72), axis=1)/81

# average
julav = np.mean(jul, axis=(0, 1))

# Reshape for month
julres = jul.reshape(720, 60, 60) 

julavres = julav.reshape(3600)

###########################################

# [ time, source, lon, lat ] Feb Bigger due to comparison to 2011 ~~~ Temp
# plot Daily average
mt71 = m7.variables['temp_a'][:, 14:24, 0, :, :]
mt72 = m7.variables['temp_a'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
jult = np.concatenate((mt71, mt72), axis=1)

# average
jultav = np.mean(jult, axis=(0, 1))

jultavres = np.reshape(jultav, 3600)

# Reshape for month
jultres = jult.reshape(720, 60, 60) 

stdjult = np.std(jultres, axis=(1, 2))

jultresav = np.mean(jultres, axis=(1,2))

##############################################################################

# [ time, source, lon, lat ] aug
# plot Daily average
m81 = m8.variables['store_Bio'][:, 14:24, 0, :, :]
m82 = m8.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
aug = np.concatenate((m81, m82), axis=1)/81

# average
augav = np.mean(aug, axis=(0, 1))

# Reshape for month
augres = aug.reshape(744, 60, 60) 

augavres = augav.reshape(3600)

###########################################

# [ time, source, lon, lat ] Feb Bigger due to comparison to 2011 ~~~ Temp
# plot Daily average
mt81 = m8.variables['temp_a'][:, 14:24, 0, :, :]
mt82 = m8.variables['temp_a'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
augt = np.concatenate((mt81, mt82), axis=1)

# average
augtav = np.mean(augt, axis=(0, 1))

augtavres = np.reshape(augtav, 3600)

# Reshape for month
augtres = augt.reshape(744, 60, 60) 

stdaugt = np.std(augtres, axis=(1, 2))

augtresav = np.mean(augtres, axis=(1,2))

##############################################################################

# [ time, source, lon, lat ] sep
# plot Daily average
m91 = m9.variables['store_Bio'][:, 14:24, 0, :, :]
m92 = m9.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
sep = np.concatenate((m91, m92), axis=1)/81

# average
sepav = np.mean(sep, axis=(0, 1))

# Reshape for month
sepres = sep.reshape(720, 60, 60) 

sepavres = sepav.reshape(3600)

###########################################

# [ time, source, lon, lat ] Feb Bigger due to comparison to 2011 ~~~ Temp
# plot Daily average
mt91 = m9.variables['temp_a'][:, 14:24, 0, :, :]
mt92 = m9.variables['temp_a'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
sept = np.concatenate((mt91, mt92), axis=1)

# average
septav = np.mean(sept, axis=(0, 1))

septavres = np.reshape(septav, 3600)

# Reshape for month
septres = sept.reshape(720, 60, 60) 

stdsept = np.std(septres, axis=(1, 2))

septresav = np.mean(septres, axis=(1,2))

##############################################################################

# [ time, source, lon, lat ] oct
# plot Daily average
m101 = m10.variables['store_Bio'][:, 14:24, 0, :, :]
m102 = m10.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
oct = np.concatenate((m101, m102), axis=1)/81

# average
octav = np.mean(oct, axis=(0, 1))

# Reshape for month
octres = oct.reshape(696, 60, 60) 

octavres = octav.reshape(3600)

###########################################

# [ time, source, lon, lat ] Feb Bigger due to comparison to 2011 ~~~ Temp
# plot Daily average
mt101 = m10.variables['temp_a'][:, 14:24, 0, :, :]
mt102 = m10.variables['temp_a'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
octt = np.concatenate((mt101, mt102), axis=1)

# average
octtav = np.mean(octt, axis=(0, 1))

octtavres = np.reshape(octtav, 3600)

# Reshape for month
octtres = octt.reshape(696, 60, 60) 

stdoctt = np.std(octtres, axis=(1, 2))

octtresav = np.mean(octtres, axis=(1,2))

##############################################################################

# [ time, source, lon, lat ] nov
# plot Daily average
m111 = mnov.variables['store_Bio'][:, 14:24, 0, :, :]
m112 = mnov.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
nov = np.concatenate((m111, m112), axis=1)/81

# average
novav = np.mean(nov, axis=(0, 1))

# Reshape for month
novres = nov.reshape(720, 60, 60) 

novavres = novav.reshape(3600)

###########################################

# [ time, source, lon, lat ] Feb Bigger due to comparison to 2011 ~~~ Temp
# plot Daily average
mt111 = mnov.variables['temp_a'][:, 14:24, 0, :, :]
mt112 = mnov.variables['temp_a'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
novt = np.concatenate((mt111, mt112), axis=1)

# average
novtav = np.mean(novt, axis=(0, 1))

novtavres = np.reshape(novtav, 3600)

# Reshape for month
novtres = novt.reshape(720, 60, 60) 

stdnovt = np.std(novtres, axis=(1, 2))

novtresav = np.mean(novtres, axis=(1,2))


##############################################################################

# [ time, source, lon, lat ] dec
# plot Daily average
m121 = mdec.variables['store_Bio'][:, 14:24, 0, :, :]
m122 = mdec.variables['store_Bio'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
dec = np.concatenate((m121, m122), axis=1)/81

# average
decav = np.mean(dec, axis=(0, 1))

# Reshape for month
decres = dec.reshape(720, 60, 60)

decavres = decav.reshape(3600) 

###########################################

# [ time, source, lon, lat ] Feb Bigger due to comparison to 2011 ~~~ Temp
# plot Daily average
mt121 = mdec.variables['temp_a'][:, 14:24, 0, :, :]
mt122 = mdec.variables['temp_a'][:, 0:14, 0, :, :]

# concatenate back to 24 hour period
dect = np.concatenate((mt121, mt122), axis=1)

# average
dectav = np.mean(dect, axis=(0, 1))

dectavres = np.reshape(dectav, 3600)

# Reshape for month
dectres = dect.reshape(720, 60, 60) 

stddect = np.std(dectres, axis=(1, 2))

dectresav = np.mean(dectres, axis=(1,2))

###################################################################################################################
# concatanate all months for Bio emissions 

year = np.concatenate((janres,febres,marres,aprres,mayres,junres,julres,augres,sepres,octres,novres,decres))

# std devs 
stdevkoeh1 = np.std(year, axis=(1, 2))
stdeight = stdevkoeh1[0:8665:12]

yearmean = np.average(year, axis=(1, 2))

#8 hour avareage 
yeareight = yearmean[0:8665:12]

#size = np.arange(0, 7, 1)




##############################################################################
# Define map parameters

max1 = 9
min1 = 0
trans = 0.5

##############################################################################

# Other variables
at = (m2.variables['lndtype'][0, :, :])
st = (m2.variables['soiltype'][0, :, :])
lai = (m1.variables['lai'][0, :, :])
laires = st.reshape(3600)
t = (m2.variables['temp_a'][:, 0:24, 0, :, :])
g = (m2.variables['skin_temp'][:, 0:24, :, :])

##############################################################################

# X axis dates instead of times for monthly 
date = np.arange(yeareight.shape[0])  # assume  delta time between data is 1 hour
date1 = (date/54.)  # use days instead of hours
ar = np.array(date1)

plt.figure(figsize=(15, 5))

plt.plot(ar, yeareight, linestyle='-', linewidth=0.6, c='b',
         label='Biogenic Emissions')

# Create standard devation fill
plt.fill_between(ar,  yeareight-stdeight,  yeareight+stdeight, alpha=0.2, edgecolor='black', 
                 facecolor='blue', linewidth=0.2, label='1 Standard devation' )
         

plt.title("CTM Domain 3 3x3 Leaf Area Index\
\n Sydney Metropolitan Region Feburary 2011", fontsize=10)

plt.xlabel('Month')
plt.ylabel('Average Total Biogenic Emissions (kg/$km^2$/hr)')

plt.ylim(0, 10)
plt.xlim(1, 12)

# tick
plt.xticks(range(1, 12, 1 ), [str(i) for i in range(1, 12, 1)])
plt.yticks(np.arange(min(yearmean)-0.01, max(yearmean)+5, 1))  # use floats for y value

plt.show()
# plt.savefig('./bmap_syd.png')
