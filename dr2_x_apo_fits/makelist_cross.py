import pyfits 
filein = '/Users/ness/new_laptop/Apogee_DR12_play/allStar-v601.fits'
filein = '/home/ness/new_laptop/Apogee_DR12/allStar-v601.fits'

hdulist = pyfits.open(filein)
datain = hdulist[1].data
apstarid= datain['APSTAR_ID'] 
targetid= datain['TARGET_ID'] 
vscatter= datain['VSCATTER'] 
aspcap_id= datain['ASPCAP_ID'] 
fields = datain['FIELD'] 
loc = datain['LOCATION_ID']
params = datain['PARAM']
t = params[:,0]
g = params[:,1]
feh = datain['FE_H']
M = params[:,3]
SNR = datain["SNR"]
aspcapflag = datain["ASPCAPFLAG"]
starflag = datain["STARFLAGS"]
alpha = params[:,-1]
alphatest = alpha
fehtest = feh 
rot = []
for each in aspcapflag:
  rot.append(each & 10**10) 
apid = []
for each in apstarid:
  apid.append(each.split('.')[-1]) 
apid = array(apid) 
#take1 = logical_and(alphatest > -0.1,logical_and(abs(fehtest) < 3,logical_and(SNR > 285, SNR < 305)))
#take1 = logical_and(alphatest > -0.4,logical_and(abs(fehtest) < 3,logical_and(SNR > 285, SNR < 305)))
#take2 = logical_and(vscatter <1,logical_and(aspcapflag ==  0, vscatter  < 1))
#take3 = loc != 1
#ind = logical_and(take3, logical_and(take1,take2)) 
#ind = fehtest > -100000
Lamost_id = []
b= open("2MASSID_edit.txt", 'r')
b2= open("files.list", 'r')
bl = b.readlines() 
bl2 = b2.readlines() 
cl = []
cl2 =[]
for each in bl:
  cl.append(each.strip())
for each in bl2:
  cl2.append(each.strip())

desig = [] 
for each in cl2:
  a = pyfits.open(each) 
  head = a[0].header
  desig.append(head['DESIG'])

raval = []
decval =[]
for each in desig:
  raval.append(float(each.split('LAMOST')[1][2:8]) )
  decval.append(float(each.split('LAMOST')[1][11:18]) )

raval_2mass = []
decval_2mass =[]
for each in cl:
  raval_2mass.append(float(each[2:8])) 
  decval_2mass.append(float(each[10:17])) 

cltake = []
cl2take = [] 
cl = array(cl)
cl2 = array(cl2) 
for one,two,cl2val,clval in zip(raval,decval,cl2,cl):
  diff1=(one - array(raval_2mass)) 
  diff2=(two - array(decval_2mass))
  difftot=(diff1**2 + diff2**2)
  takediff = argsort(difftot)[0] 
  cl2take.append(cl2val)
  cltake.append(cl[takediff])

# cl2take = 2mass id , 
# cltake = LAMOST ID 
savetxt("2MASSID_ordered.txt", cl2take, fmt = "%s") 
savetxt("LAMOSTID_ordered.list", cltake, fmt = "%s") 
########################## BELOW IS WHERE READ IN AND GET THE VALUES OUT 

b= open("2MASSID_ordered.txt", 'r')
b2= open("LAMOSTID_ordered.list", 'r')
bl = b.readlines() 
bl2 = b2.readlines() 
cl = [] 
for each in bl2:
  cl.append(each.split()[0]) 
ind = [] 
for each in cl:
  if each in apid:
    ind.append(list(apid).index(each)) 

loctake = loc[ind]
idtake = apid[ind]  
idtake_LAMOST =  cl2 
t_take = t[ind]
g_take = g[ind]
M_take = M[ind]
feh_take = feh[ind]
alpha_take = alpha[ind]
snr_take = SNR[ind]
vscatter_take = vscatter[ind]
rot = array(rot) 
rot_take = rot[ind]
aspcap_take = aspcapflag[ind]
starflag_take = starflag[ind]

text2 ='/home/ness/new_laptop/Apogee_DR12/data.sdss3.org/sas/dr12/apogee/spectro/redux/r5/stars/l25_6d/v603/'
text3 = '/aspcapStar-r5-v603-2M' 
text4 = '.fits'
listin = [] 
listin_apogee = [] 
for a,b in zip(loctake, idtake):
  listin.append(text2+str(a)+text3+str(b)+text4) 
  listin_apogee.append(str(a)+"_"+str(b) ) 
data = zip(listin_apogee,cltake, cl2take, t_take, g_take, M_take, feh_take, alpha_take,snr_take,vscatter_take,rot_take,aspcap_take,starflag_take) 
data_more = zip(listin, t_take, g_take, M_take, feh_take, alpha_take,snr_take,vscatter_take,rot_take,aspcap_take,starflag_take) 
savetxt('training_LAMOST_ordered.list', data, fmt = "%s" ) 

take1 = logical_and(alpha_take > -0.4,logical_and(abs(feh_take) < 3,logical_and(snr_take > 400, snr_take < 12000305)))
take2 = logical_and(vscatter_take <1,logical_and(aspcap_take ==  0, vscatter_take  < 1))
ind2 = logical_and(take1,take2) 
data = array(data) 
data2 = data[ind2]
#data3 = data_more[ind2]
savetxt("training_LAMOST_highsnr_ordered.txt", data2, fmt = "%s")
#savetxt("training_LAMOST_highsnr_names.txt", data3, fmt = "%s")

#data3 = zip(listin, t_take, g_take, M_take, feh_take)
#savetxt("LAMOST_Apogeevals_tgmfeh.txt", data3, fmt = "%s") 

