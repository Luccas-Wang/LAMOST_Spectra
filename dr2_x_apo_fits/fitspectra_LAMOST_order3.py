import pyfits 
import glob
from scipy import optimize as opt
import pickle 
#from glob import glob
import scipy
from scipy import interpolate 
#listin = 'files.txt' 
normed_training_data = 'normed_data.pickle'
LARGE = 200.
#listin = 'test.txt' 
#a = open(listin, 'r')
#c = open("fileout.txt", 'r')
#savenames = c.readlines() 
#b = pyfits.open(a) 
#head = b[0].header 
#redshift = head['Z']
#start_wl_log10 =  head['COEFF0']
#diff_wl = head['CD1_1']
#start_pix = head['CRPIX1']
#SNR = a[0].header['SNR']
#ydata = a[0].data[0]
#xdata_in = a[0].data[1]
#redshift = head['Z']
#shift = redshift*xdata_in
#xdata = xdata_in + shift
#ivar = a[0].data[2]
#ysigma = 1/ivar**0.5

def weighted_median(values, weights, quantile):
    """weighted_median

    keywords
    --------

    values: ndarray
        input values

    weights: ndarray
        weights to apply to each value in values

    quantile: float
        quantile selection

    returns
    -------
    val: float
        median value
    """
    sindx = np.argsort(values)
    cvalues = 1. * np.cumsum(weights[sindx])
    cvalues = cvalues / cvalues[-1]
    foo = sindx[cvalues > quantile]
    if len(foo) == 0:
        return values[0]
    indx = foo[0]
    return values[indx]

def continuum_normalize(dataall, SNRall, delta_lambda=50):
    """continuum_normalize

    keywords
    --------

    dataall: ndarray, shape=(Nlambda, Nstar, 3)
        wavelengths, flux densities, errors

    delta_lambda:
        half-width of median region in angstroms


    returns
    -------
    continuum:     (Nlambda, Nstar)
        continuum level

    .. note::

        * does a lot of stuff *other* than continuum normalization

    .. todo::

        * bugs: for loops!
    """
    Nlambda, Nstar, foo = dataall.shape
    continuum = np.zeros((Nlambda, Nstar))
   
   # sanitize inputs
    for jj in range(Nstar):
        bad_a = np.logical_or(np.isnan(dataall[:, jj, 1]) ,np.isinf(dataall[:,jj, 1]))
        bad_b = np.logical_or(dataall[:, jj, 2] <= 0. , np.isnan(dataall[:, jj, 2]))
        bad = np.logical_or( np.logical_or(bad_a, bad_b) , np.isinf(dataall[:, jj, 2]))
        dataall[bad, jj, 1] = 0.
        dataall[bad, jj, 2] = np.Inf 
        continuum = np.zeros((Nlambda, Nstar))
    assert foo == 3
    for star in range(Nstar):
        print "get_continuum(): working on star" ,star
        for ll, lam in enumerate(dataall[:, 0, 0]):
            if dataall[ll, star, 0] != lam:
                print dataall[ll,star,0], lam , dataall[ll,0,0] 
                print ll, star 
                #print ll+1, star+1, dataall[ll+1, star+1, 0], dataall[ll+1,0,0] 
                assert False
            indx = (np.where(abs(dataall[:, star, 0] - lam) < delta_lambda))[0]
            ivar = 1. / (dataall[indx, star, 2] ** 2)
            ivar = np.array(ivar)
            q = 0.90
            continuum[ll, star] = weighted_median(dataall[indx, star, 1], ivar, q)
    for jj in range(Nstar):
        bad = np.where(continuum[:,jj] <= 0) 
        continuum[bad,jj] = 1.
        dataall[:, jj, 1] /= continuum[:,jj]
        dataall[:, jj, 2] /= continuum[:,jj]
        dataall[bad,jj, 1] = 1. 
        dataall[bad,jj, 2] = LARGE 
        bad = np.where(dataall[:, jj, 2] > LARGE) 
        dataall[bad,jj, 1] = 1. 
        dataall[bad,jj, 2] = LARGE 
    return dataall  


def get_normalized_training_data_reference():
  if glob.glob(normed_training_data): 
        file_in2 = open(normed_training_data, 'r') 
        dataall, metaall, labels, cluster_name, ids = pickle.load(file_in2)
        file_in2.close()
        return dataall, metaall, labels, Ametaall, cluster_name, ids
  fn = 'training_LAMOST.txt'  
  T_est,g_est,feh_est= np.loadtxt(fn, usecols = (1,2,3), unpack =1) 
  labels = ["teff", "logg", "feh"]
  a = open(fn, 'r') 
  al = a.readlines() 
  bl = []
  cluster_name = [] 
  ids = []
  for each in al:
    bl.append(each.split()[0]) 
    cluster_name.append(each.split()[1]) 
    ids.append(each.split()[0].split('-2M')[-1].split('.fits')[0])

  SNRall = np.zeros((len(bl))) 
  for jj,each in enumerate(bl):
    each = each.strip('\n')
    a = pyfits.open(each) 
    b = pyfits.getheader(each) 
    SNRall[jj] = a[0].header['SNR']
    start_wl =  a[1].header['CRVAL1']
    diff_wl = a[1].header['CDELT1']
    print np.atleast_2d(a[1].data).shape
    if jj == 0:
      nmeta = len(labels)
      nlam = len(a[1].data)
    val = diff_wl*(nlam) + start_wl 
    wl_full_log = np.arange(start_wl,val, diff_wl) 
    ydata = (np.atleast_2d(a[1].data))[0] 
    ydata_err = (np.atleast_2d(a[2].data))[0] 
    ydata_flag = (np.atleast_2d(a[3].data))[0] 
    assert len(ydata) == nlam
    ydata = np.array(ydata)
    ydata_err = np.array(ydata_err)
    starname2 = each.split('.fits')[0]+'.txt'
    sigma = (np.atleast_2d(a[2].data))[0]# /y1
    if jj == 0:
      npix = len(xdata) 
      dataall = np.zeros((npix, len(bl), 3))
      metaall = np.ones((len(bl), nmeta))
      Ametaall = np.ones((len(bl), nmeta))
    if jj > 0:
      assert xdata[0] == dataall[0, 0, 0]

    dataall[:, jj, 0] = xdata
    dataall[:, jj, 1] = ydata
    dataall[:, jj, 2] = sigma

    for k in range(0,len(bl)): 
        # must be synchronised with labels 
      metaall[k,0] = T_est[k] 
      metaall[k,1] = g_est[k] 
      metaall[k,2] = feh_est[k] 
  dataall = continuum_normalize(dataall,SNRall) #dataall

  file_in = open(normed_training_data, 'w')  
  pickle.dump((dataall, metaall, labels, cluster_name, ids),  file_in)
  file_in.close()
  return dataall, metaall, labels , cluster_name, ids


def get_normalized_test_data(testfile_in):
  """
    inputs
    ------
    testfile: str
        the file in with the list of fits files want to test - if normed, move on,
        if not normed, norm it
    if not noisify carry on as normal, otherwise do the noise tests

    returns
    -------
    testdata:
  """
  a = open(testfile_in, 'r')
  al = a.readlines()
  #d = open(testfile_out, 'r') 
  #dl = d.readlines() 
  namein = []
  nameout = []
  for each in al:
    namein.append(each.split()[2]) 
    nameout.append(each.split()[0]) 

  filename = testfile_in.split('.list')[0]
  if glob.glob(str(filename)+'.pickle'):
    file_in2 = open(str(filename)+'.pickle', 'r') 
    testdata,nameout_all = pickle.load(file_in2)
    file_in2.close()
    return testdata, nameout_all 
 # 
  SNRall = np.zeros(len(namein))
  nameout_all = [] 
  for jj,each in enumerate(namein):
    a = pyfits.open(each) 
    xdata_in = a[0].data[2] 
    redshift = a[0].header['Z']
    shift = redshift*xdata_in
    xdata= xdata_in - shift 
    ydata = a[0].data[0] 
    ysigma = 1/(a[0].data[1])**0.5
    newxgrid = arange(3900,8800,0.85) 
    ydata = resample_grid(xdata, ydata, newxgrid) 
    ysigma = resample_grid(xdata, ysigma, newxgrid) 
    xdata = newxgrid
    len_data = len(a[0].data[0]) 
    wl_low = 3700.
    wl_high = 9000.
    if jj == 0:
      #nlam = len(a[0].data[0])#[logical_and(xdata > wl_low, xdata < wl_high)])
      nlam = len(xdata) 
      testdata = np.zeros((nlam, len(namein), 3))
    SNR = median(ydata/ysigma)  
    SNRall[jj] = SNR
    testdata[:, jj, 0] = xdata#[logical_and(xdata > wl_low , xdata < wl_high)]
    testdata[:, jj, 1] = ydata#[logical_and(xdata > wl_low , xdata < wl_high)]
    testdata[:, jj, 2] = ysigma#[logical_and(xdata > wl_low , xdata < wl_high)]
    nameout_take = nameout[jj]
    nameout_all.append(nameout_take)
  testdata = continuum_normalize(testdata,SNRall) # testdata
  #file_in = open(str(nameout_take)+'.pickle', 'w')  
  #file_in2 = open(str(nameout_take)+'_SNR.pickle', 'w')
  #file_in = open("allstars"+'.pickle', 'w')  
  #file_in2 = open("allstars"+'_SNR.pickle', 'w')
  file_in = open(filename+'.pickle', 'w')  
  file_in2 = open(filename+'_SNR.pickle', 'w')
  pickle.dump((testdata,nameout_all),  file_in)
  pickle.dump(SNRall,  file_in2)
  file_in.close()
  file_in2.close()
  return testdata, nameout_all # not yet implemented but at some point should probably save ids into the normed pickle file 

def resample_grid(xdata, ydata,newxgrid):
  f = interpolate.interp1d(xdata, ydata)
  new_ydata= f(newxgrid)
  return new_ydata


def get_normalized_training_data(testfile_in):
  """
    inputs
    ------
    testfile: str
        the file in with the list of fits files want to test - if normed, move on,
        if not normed, norm it
    if not noisify carry on as normal, otherwise do the noise tests

    returns
    -------
    testdata:
  """
  if glob.glob(normed_training_data): 
      file_in2 = open(normed_training_data, 'r') 
      dataall, metaall, labels, cluster_name, ids = pickle.load(file_in2)
      file_in2.close()
  return dataall, metaall, labels, cluster_name, ids
  #fn = 'training_LAMOST_highsnr_ordered.txt'  
  fn = testfile_in 
  T_est,g_est,feh_est= np.loadtxt(fn, usecols = (3,4,5), unpack =1) 
  labels = ["teff", "logg", "feh"]
  
  a = open(testfile_in, 'r')
  al = a.readlines()
  namein = []
  nameout = []
  for each in al:
    namein.append(each.split()[2]) 
    nameout.append(each.split()[0]) 
  
 # 
  SNRall = np.zeros(len(namein))
  nameout_all = [] 
  namein_all = [] 
  for jj,each in enumerate(namein):
    a = pyfits.open(each) 
    xdata_in = a[0].data[2] 
    redshift = a[0].header['Z']
    shift = redshift*xdata_in
    xdata= xdata_in - shift 
    ydata = a[0].data[0] 
    ysigma = 1/(a[0].data[1])**0.5
    newxgrid = arange(3900,8800,0.85) 
    ydata = resample_grid(xdata, ydata, newxgrid) 
    ysigma = resample_grid(xdata, ysigma, newxgrid) 
    xdata = newxgrid
    len_data = len(a[0].data[0]) 
    wl_low = 3700.
    wl_high = 9000.
    if jj == 0:
      #nlam = len(a[0].data[0])#[logical_and(xdata > wl_low, xdata < wl_high)])
      nmeta = len(labels)
      nlam = len(xdata) 
      testdata = np.zeros((nlam, len(namein), 3))
      metaall = np.ones((len(namein), nmeta))
    SNR = median(ydata/ysigma)  
    SNRall[jj] = SNR
    testdata[:, jj, 0] = xdata#[logical_and(xdata > wl_low , xdata < wl_high)]
    testdata[:, jj, 1] = ydata#[logical_and(xdata > wl_low , xdata < wl_high)]
    testdata[:, jj, 2] = ysigma#[logical_and(xdata > wl_low , xdata < wl_high)]
    nameout_take = nameout[jj]
    namein_take = namein[jj]
    nameout_all.append(nameout_take)
    namein_all.append(namein_take)
  for k in range(0,len(namein)): 
    metaall[k,0] = T_est[k] 
    metaall[k,1] = g_est[k] 
    metaall[k,2] = feh_est[k] 

  testdata = continuum_normalize(testdata,SNRall) # testdata
  file_in = open("normed_data"+'.pickle', 'w')  
  #file_in = open(normed_training_data, 'w')  
  pickle.dump((testdata, metaall, labels, nameout_all, namein_all),  file_in)
  file_in.close()
  return testdata, nameout_all # not yet implemented but at some point should probably save ids into the normed pickle file 

def resample_grid(xdata, ydata,newxgrid):
  f = interpolate.interp1d(xdata, ydata)
  new_ydata= f(newxgrid)
  return new_ydata

def infer_labels_nonlinear(fn_pickle,testdata, ids, fout_pickle, weak_lower,weak_upper):
#def infer_labels(fn_pickle,testdata, fout_pickle, weak_lower=0.935,weak_upper=0.98):
    """
    """
    file_in = open(fn_pickle, 'r') 
    dataall, metaall, labels, offsets, coeffs, covs, scatters,chis,chisq = pickle.load(file_in)
    file_in.close()
    nstars = (testdata.shape)[1]
    nlabels = len(labels)
    Params_all = np.zeros((nstars, nlabels))
    MCM_rotate_all = np.zeros((nstars, np.shape(coeffs)[1]-1, np.shape(coeffs)[1]-1.))
    covs_all = np.zeros((nstars,nlabels, nlabels))
    for jj in range(0,nstars):
      #if np.any(testdata[:,jj,0] != dataall[:, 0, 0]):
      if np.any(abs(testdata[:,jj,0] - dataall[:, 0, 0]) > 0.0001): 
          print testdata[range(5),jj,0], dataall[range(5),0,0]
          assert False
      xdata = testdata[:,jj,0]
      ydata = testdata[:,jj,1]
      ysigma = testdata[:,jj,2]
      ydata_norm = ydata  - coeffs[:,0] # subtract the mean 
      f = ydata_norm 
      t,g,feh = metaall[:,0], metaall[:,1], metaall[:,2]
      x0,x1,x2,x3,x4,x5,x6,x7,x8,x9 = coeffs[:,0], coeffs[:,1], coeffs[:,2], coeffs[:,3], coeffs[:,4], coeffs[:,5], coeffs[:,6] ,coeffs[:,7], coeffs[:,8], coeffs[:,9] 
      x10,x11,x12 = coeffs[:,10], coeffs[:,11], coeffs[:,12] 
      Cinv = 1. / (ysigma ** 2 + scatters ** 2)
      Params,covs = nonlinear_invert(f, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,x11,x12,1/Cinv**0.5 ) 
      Params = Params+offsets 
      coeffs_slice = coeffs[:,-12:]
      MCM_rotate = np.dot(coeffs_slice.T, Cinv[:,None] * coeffs_slice)
      Params_all[jj,:] = Params 
      MCM_rotate_all[jj,:,:] = MCM_rotate 
      covs_all[jj,:,:] = covs
    filein = fout_pickle.split('_tags') [0] 
    if logical_or(filein == 'self_2nd_order', filein[0:4] == 'self'):
      file_in = open(fout_pickle, 'w')  
      file_normed = normed_training_data.split('.pickle')[0]
      chi2 = get_goodness_fit(fn_pickle, file_normed, Params_all )
      chi2_def = chi2/(len(xdata)*1. - 3.) #7209. #8575-363-3 # len(xdata)*1.
      pickle.dump((Params_all, covs_all,chi2_def,ids),  file_in)
      file_in.close()
    else: 
      chi2 = get_goodness_fit(fn_pickle, filein, Params_all )
      chi2_def = chi2/len(xdata)*1.
      file_in = open(fout_pickle, 'w')  
      pickle.dump((Params_all, covs_all, chi2_def, ids),  file_in)
      file_in.close()
    return Params_all , MCM_rotate_all

def func(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, a, b, c):
    f = (0 
         + x1*a 
         + x2*b 
         + x3*c 
         + x4* a**2# 
         + x5 * a * b
         + x6 * a * c 
         + x7*b**2
         + x8  * b * c 
         + x9*c**2
         + x10*a**3 
         + x11*b**3 
         + x12*c**3 ) 
    return f

# thankyou stack overflow for the example below on how to use the optimse function  
def nonlinear_invert(f, x1, x2, x3, x4, x5, x6, x7, x8, x9 ,x10, x11, x12, sigmavals):
    def wrapped_func(observation_points, a, b, c):
        x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12  = observation_points
        return func(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12,  a, b, c)

    xdata = np.vstack([x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12])
    model, cov = opt.curve_fit(wrapped_func, xdata, f, sigma = sigmavals, maxfev = 2000)#absolute_sigma = True)  is not an option in my version of scipy will upgrade scipy
    return model, cov 

def get_goodness_fit(fn_pickle, filein, Params_all ):
    fd = open(fn_pickle,'r')
    dataall, metaall, labels, offsets, coeffs, covs, scatters, chis, chisq = pickle.load(fd) 
    fd.close() 
    file_with_star_data = str(filein)+".pickle"
    file_normed = normed_training_data.split('.pickle')[0]
    if filein != file_normed: 
      f_flux= open(file_with_star_data, 'r') 
      flux,num = pickle.load(f_flux) 
    if filein == file_normed: 
      f_flux = open('normed_data.pickle', 'r') 
      flux, metaall, labelthings,  cluster_name, ids = pickle.load(f_flux)
    f_flux.close() 
    labels = Params_all 
    nlabels = shape(labels)[1]
    nstars = shape(labels)[0]
    features_data = np.ones((nstars, 1))
    offsets = np.mean(labels, axis = 0) 
    features_data = np.hstack((features_data, labels - offsets)) 
    newfeatures_data = np.array([np.outer(m, m)[np.triu_indices(nlabels)] for m in (labels - offsets)])
    newfeatures2_data = [a**3 for a in ( labels - offsets)] 
    features_data = np.hstack((features_data, newfeatures_data, newfeatures2_data))
    chi2_all = np.zeros(nstars) 
    for jj in range(nstars):
        model_gen = np.dot(coeffs,features_data.T[:,jj]) 
        data_star = flux[:,jj,1] 
        Cinv = 1. / (flux[:,jj, 2] ** 2 + scatters ** 2)  # invvar slice of data
        chi2 = sum( (Cinv) * (data_star - np.dot(coeffs, features_data.T[:,jj]))**2) 
        chi2_all[jj] = chi2
    return chi2_all 

def train(dataall, metaall, order, fn, cluster_name, logg_cut=100., teff_cut=0., leave_out=None):
    """
    - `leave out` must be in the correct form to be an input to `np.delete`
    """
    if leave_out is not None: #
        dataall = np.delete(dataall, [leave_out], axis = 1) 
        metaall = np.delete(metaall, [leave_out], axis = 0) 
   
    nstars, nmeta = metaall.shape
    
    offsets = np.mean(metaall, axis=0)
    features = np.ones((nstars, 1))
    if order >= 1:
        features = np.hstack((features, metaall - offsets)) 
    if order >= 2:
        newfeatures = np.array([np.outer(m, m)[np.triu_indices(nmeta)] for m in (metaall - offsets)])
        newfeatures2 = [a**3 for a in ( metaall - offsets)] 
        features = np.hstack((features, newfeatures, newfeatures2))

    blob = do_regressions(dataall, features)
    coeffs = np.array([b[0] for b in blob])
    covs = np.array([np.linalg.inv(b[1]) for b in blob])
    chis = np.array([b[2] for b in blob])
    chisqs = np.array([np.dot(b[2],b[2]) - b[3] for b in blob]) # holy crap be careful
    scatters = np.array([b[4] for b in blob])

    fd = open(fn, "w")
    pickle.dump((dataall, metaall, labels, offsets, coeffs, covs, scatters,chis,chisqs), fd)
    fd.close()
    return 

def do_regressions(dataall, features):
    """
    """
    nlam, nobj, ndata = dataall.shape
    nobj, npred = features.shape
    featuresall = np.zeros((nlam,nobj,npred))
    featuresall[:, :, :] = features[None, :, :]
    return map(do_one_regression, dataall, featuresall)

def do_one_regression_at_fixed_scatter(data, features, scatter):
    """
    Parameters
    ----------
    data: ndarray, [nobjs, 3]
        wavelengths, fluxes, invvars

    meta: ndarray, [nobjs, nmeta]
        Teff, Feh, etc, etc

    scatter:


    Returns
    -------
    coeff: ndarray
        coefficients of the fit

    MTCinvM: ndarray
        inverse covariance matrix for fit coefficients

    chi: float
        chi-squared at best fit

    logdet_Cinv: float
        inverse of the log determinant of the cov matrice
        :math:`\sum(\log(Cinv))`
    """
    # least square fit
    #pick = logical_and(data[:,1] < np.median(data[:,1]) + np.std(data[:,1])*3. , data[:,1] >  median(data[:,1]) - np.std(data[:,1])*3.)#5*std(data[:,1]) ) 
    Cinv = 1. / (data[:, 2] ** 2 + scatter ** 2)  # invvar slice of data
    M = features
    MTCinvM = np.dot(M.T, Cinv[:, None] * M) # craziness b/c Cinv isnt a matrix
    x = data[:, 1] # intensity slice of data
    MTCinvx = np.dot(M.T, Cinv * x)
    try:
        coeff = np.linalg.solve(MTCinvM, MTCinvx)
    except np.linalg.linalg.LinAlgError:
        print MTCinvM, MTCinvx, data[:,0], data[:,1], data[:,2]
        print features
    assert np.all(np.isfinite(coeff)) 
    chi = np.sqrt(Cinv) * (x - np.dot(M, coeff)) 
    logdet_Cinv = np.sum(np.log(Cinv)) 
    return (coeff, MTCinvM, chi, logdet_Cinv )

def do_one_regression(data, metadata):
    """
    does a regression at a single wavelength to fit calling the fixed scatter routine
    # inputs:
    """
    ln_s_values = np.arange(np.log(0.0001), 0., 0.5)
    chis_eval = np.zeros_like(ln_s_values)
    for ii, ln_s in enumerate(ln_s_values):
        foo, bar, chi, logdet_Cinv = do_one_regression_at_fixed_scatter(data, metadata, scatter = np.exp(ln_s))
        chis_eval[ii] = np.sum(chi * chi) - logdet_Cinv
    if np.any(np.isnan(chis_eval)):
        s_best = np.exp(ln_s_values[-1])
        return do_one_regression_at_fixed_scatter(data, metadata, scatter = s_best) + (s_best, )
    lowest = np.argmin(chis_eval)
    #if lowest == 0 or lowest == len(ln_s_values) + 1:
    if lowest == 0 or lowest == len(ln_s_values)-1:
        s_best = np.exp(ln_s_values[lowest])
        return do_one_regression_at_fixed_scatter(data, metadata, scatter = s_best) + (s_best, )
    #print data
    #print metadata
    #print "LOWEST" , lowest
    ln_s_values_short = ln_s_values[np.array([lowest-1, lowest, lowest+1])]
    chis_eval_short = chis_eval[np.array([lowest-1, lowest, lowest+1])]
    z = np.polyfit(ln_s_values_short, chis_eval_short, 2)
    f = np.poly1d(z)
    fit_pder = np.polyder(z)
    fit_pder2 = pylab.polyder(fit_pder)
    s_best = np.exp(np.roots(fit_pder)[0])
    return do_one_regression_at_fixed_scatter(data, metadata, scatter = s_best) + (s_best, )


#fpickle2 = "coeffs_2nd_order.pickle"
#dataall, metaall, labels, cluster_name, ids = get_normalized_training_data('training_LAMOST_highsnr_ordered.txt')
#if not glob.glob(fpickle2):
#  train(dataall, metaall, 2,  fpickle2,  cluster_name, logg_cut= 40.,teff_cut = 0.)
self_flag = 2. 
self_flag = 1. 
if self_flag == 2:
#  field = "self_2nd_order_"
#  file_in = open(normed_training_data, 'r') 
#  testdataall, metaall, labels,  cluster_name,ids = pickle.load(file_in)
#  file_in.close() 
  testmetaall, inv_covars = infer_labels_nonlinear("coeffs_2nd_order.pickle", testdataall, ids, field+"tags.pickle",-10.950,10.99) 

if self_flag == 1:
  #testdataall, ids = get_normalized_test_data('training_LAMOST_ordered_test.list')
  #testmetaall, inv_covars = infer_labels_nonlinear("coeffs_2nd_order.pickle", testdataall, ids, "training_LAMOST_ordered_test_tags.pickle",0.00,1.40) 
  #testdataall, ids = get_normalized_test_data('training_LAMOST_ordered.list')
  #testmetaall, inv_covars = infer_labels_nonlinear("coeffs_2nd_order.pickle", testdataall, ids, "training_LAMOST_ordered_tags.pickle",0.00,1.40) 
  testdataall, ids = get_normalized_test_data('training_LAMOST_ordered_notintraining.list')
  testmetaall, inv_covars = infer_labels_nonlinear("coeffs_2nd_order.pickle", testdataall, ids, "training_LAMOST_ordered_notintraining_tags.pickle",0.00,1.40) 
  #testdataall, ids = get_normalized_test_data('training_LAMOST_highsnr_ordered.txt')
