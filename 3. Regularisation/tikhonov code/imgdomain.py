# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from IPython import get_ipython

# plotting in an external window
get_ipython().run_line_magic('matplotlib', 'qt') 

def showimg(img):
    """ show the each coil img """
    if np.ndim(img) == 3:
        (ny,nx,nc) = np.shape(img)
        
        fig, ax = plt.subplots(1,nc,figsize = [18,3])
        n = 0
        slc = 0
        for _ in range(nc):
            ax[n].imshow(img[:,:,slc], 'gray')
            ax[n].set_xticks([])
            ax[n].set_yticks([])
            n += 1
            slc += 1
           
        fig.subplots_adjust(wspace=0, hspace=0)
        plt.show()
        
    else:
        nc = 1
        plt.imshow(img,'gray')
    
    
    
def getcoilsen(img):
    """ get the coil sensitivity by using root sum squared(RSS) method """
    (ny,nx,nc) = np.shape(img)
    
    coil_sensitivity = np.zeros_like(img, dtype='complex')
    
    rss_img = np.sqrt(np.sum(abs(img)**2,2))
    
    for slice_idx in range(nc):
        coil_sensitivity[:,:,slice_idx] = img[:,:,slice_idx] / rss_img
    
    return coil_sensitivity


def intercoil_noise_covariance(img, patch_size):
    (ny,nx,nc) = np.shape(img)
    
    vectorised_noise = np.zeros((patch_size**2,nc), dtype = 'complex')
    
    for coil_idx in range(nc):
        vectorised_noise[:,coil_idx] = np.ndarray.flatten(img[nx-patch_size:nx,nx-patch_size:nx,coil_idx])
        
    noise_cov = np.matmul(np.matrix.getH(vectorised_noise),vectorised_noise)
    
    return noise_cov
        

# =============================================================================
#def lcurve(A_tilde,y_tilde,x_zero,nc,r_fac):
#     U,S,Vh = np.linalg.svd(A_tilde)
#     
#     iteration = 0
#     lmbd = 0
#     tot_err = 0
#     while 1:
#         err_old = tot_err
#         
#         model_error = abs(np.sum(((lmbd**2/(S**2 + lmbd**2))[:,None] * (U[0:r_fac] @ y_tilde))**2, axis = 0))
#         prior_error = abs(np.sum(((S**2/(S**2+lmbd**2))[:,None] * (((U[0:r_fac] @ y_tilde) / S[:,None]) - x_zero))**2, axis = 0))
#         
#         tot_err = model_error + prior_error
#         lmbd_old = lmbd
#         
#         if ((tot_err - err_old > 0) and iteration > 2) or iteration > 10 or lmbd > 1 :
#             return lmbd_old
#             break
#         
#         else:
#             iteration += 1
#             
# #            diff_model_error = abs(np.sum( ((S**2/(lmbd**2 + S**2)**2)[:,None] * (U[0:r_fac] @ y_tilde) )**2 , axis = 0 ))
# #            diff_prior_error = abs(np.sum( ( (S**2/(S**2+lmbd**2))[:,None] * (((U[0:r_fac] @ y_tilde) / S[:,None]) - x_zero) )**2, axis = 0))
# #            diff_tot_err = diff_model_error + diff_prior_error
#              
# #            lmbd = lmbd + (10**7)*(tot_err / diff_tot_err)
#             lmbd = lmbd + 0.0001
# =============================================================================